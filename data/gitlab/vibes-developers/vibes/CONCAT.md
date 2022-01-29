FHI-vibes
===

Welcome to `FHI-vibes`, a `python` package for calculating, analyzing, and understanding the vibrational properties of solids from first principles. `FHI-vibes` is intended to seamlessly bridge between the harmonic approximation and fully anharmonic molecular dynamics simulations. To this end, `FHI-vibes` builds on several [existing packages](https://vibes-developers.gitlab.io/vibes/Credits/) and interfaces them in a consistent and user-friendly fashion. 

In the documentation and tutorials, knowledge of first-principles electronic-structure theory as well as proficiency with _ab initio_ codes such as [FHI-aims](https://aimsclub.fhi-berlin.mpg.de/) and high-performance computing are assumed. Additional experience with Python, the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/), or [Phonopy](https://atztogo.github.io/phonopy/) is helpful, but not needed.

`FHI-vibes` provides the following features:

- Geometry optimization via [ASE](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#module-ase.optimize),
- harmonic phonon calculations via [Phonopy](https://atztogo.github.io/phonopy/),
- molecular dynamics simulations in [NVE](https://wiki.fysik.dtu.dk/ase/ase/md.html#constant-nve-simulations-the-microcanonical-ensemble), [NVT](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin), and [NPT](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.nptberendsen) ensembles,
- [harmonic sampling](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.115504), and
- [anharmonicity quantification](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.4.083809).

Most of the functionality is high-throughput ready via [fireworks](https://materialsproject.github.io/fireworks/#).

## Overview

- [Installation](https://vibes-developers.gitlab.io/vibes/Installation)
- [Tutorial](https://vibes-developers.gitlab.io/vibes/Tutorial/0_intro)
- [Documentation](https://vibes-developers.gitlab.io/vibes/Documentation/0_intro)
- [Credits](https://vibes-developers.gitlab.io/vibes/Credits)
- [References](https://vibes-developers.gitlab.io/vibes/References)


## News

- `FHI-vibes` got [published in JOSS](https://joss.theoj.org/papers/10.21105/joss.02671)!
- [Our anharmonicity measure got published!](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.4.083809)
- [… the best is yet to come.](https://www.youtube.com/watch?v=B-Jq26BCwDs)

## Changelog

#### v1.0.3

- update dependencies to allow `phonopy` versions up to 2.8

#### v1.0.2

- First official release after passing the [JOSS review](https://github.com/openjournals/joss-reviews/issues/2671).
- Several additions to the documentation.

#### v1.0.0a10

- Enable conversion of trajectories to `ase.io.Trajectory` files for viewing with ASE [(!37)](https://gitlab.com/vibes-developers/vibes/-/merge_requests/37)
- Important fix for running NPT dynamics [(!36)](https://gitlab.com/vibes-developers/vibes/-/merge_requests/36)
- We have a changelog now!
# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the maintainers of this repository. This will make life easier for everyone.

## Report Issues

Please use the [issue tracker](https://gitlab.com/vibes-developers/vibes/-/issues) to report issues. Please try to answer these questions:

- Has this issue been discussed before? Please have a quick look at the existing issues. If not:
- What is the issue? What is the expected behavior?
- Is the problem reproducible? Please provide a _minimal_ example.


## Contribute Code via Merge Request

In order to contribute code to `FHI-vibes`, please follow the usual steps for [preparing and creating a merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html). A few remarks regarding our guidelines for code and code style:

- We use [black](https://black.readthedocs.io/en/stable/) with default settings and [isort](https://pycqa.github.io/isort/) for formatting the code. The settings for `isort` are included in `setup.cfg`.
- Please _document_ and _test_ your changes. Tests are found in `vibes/tests` and written with [pytest](https://docs.pytest.org/en/stable/). 
- Please use [google-type docstrings](https://google.github.io/styleguide/pyguide.html) for your functions. Optionally you can use type hints, but we currently don't enforce this.
- We loosely keep track of code coverage, please try not to decrease coverage when contributing new code.---
title: 'FHI-vibes: _Ab Initio_ Vibrational Simulations'
tags:
  - Python
  - Physics
  - Phonons
  - Transport
authors:
  - name: Florian Knoop
    orcid: 0000-0002-7132-039X
    affiliation: 1
  - name: Thomas A. R. Purcell
    orcid: 0000-0003-4564-7206
    affiliation: 1
  - name: Matthias Scheffler
    affiliation: 1
  - name: Christian Carbogno
    orcid: 0000-0003-0635-8364
    affiliation: 1
affiliations:
 - name: Fritz Haber Institute of the Max Planck Society, Berlin, Germany
   index: 1
date: July 2020
bibliography: paper.bib
---

# Introduction

The vibrational motion of nuclei determines many important properties of materials, including their thermodynamic equilibrium and non-equilibrium properties. Accurately assessing the nuclear dynamics and the associated material properties is therefore an important task for computational materials scientists in a broad range of sub-fields. Of particular importance are simulation techniques that build on first-principles electronic-structure simulations and thereby allow to systematically investigate the virtually infinite space of materials, including those systems for which little or no experimental data is hitherto available [@Curtarolo2013]. This allows one
to design novel and improved materials with optimal properties for many applications, e.g., high-performance thermal insulators for gas and airplane turbines [@Evans2008], organic semicondcutors with long-term phase stabilities [@Salzillo2016], thermoelectric generators [@Snyder2008], and improved thermal management systems [@Tian2019].

Essentially, there are two distinct routes towards assessing vibrational properties:
In perturbative _lattice dynamics_ techniques, the potential-energy surface on which the nuclei move is approximated with a Taylor expansion around the equilibrium structure. 
Typically, one starts from a second-order expansion, i.e., the _harmonic approximation_, which allows for an analytic solution of the equations of motion [@Dove1993] 
and thus for a straightforward evaluation of observables (thermodynamic expectation values). Higher-order terms in the Taylor expansion can be accounted for perturbatively.
Conversely, _molecular dynamics_ (MD) based approaches account for the full, non-perturbative potential-energy surface _without_ approximating the actual interactions. This requires
one to solve the equations of motion numerically by propagating the atoms in time; physical properties can then be extracted as time  averages of properly chosen observables [@Tuckerman2010].
Although both _lattice dynamics_ and _molecular dynamics_ techniques aim at computing the same physical observables, the involved methodologies, formalisms, and challenges are quite different.
Accordingly, both methodologies also have different strengths and weaknesses: For instance, performing and analyzing MD simulations is typically computationally and conceptually more challenging,
whereas perturbative lattice dynamics calculations inherently rely on approximations that are hard to validate.

To date, a variety of different software packages exists at different degrees of sophistication in both fields. Prominent examples are the _phonopy_ code [@Togo2015] for performing _lattice dynamics_
calculations using Parlinski's finite-difference formalism [@Parlinski1997] and the _i-PI_ code [@Kapil2019] for performing classical MD and quantum-mechanical path-integral MD simulations.
Both packages interface with a variety of first-principles codes like *VASP* [@Kresse1996], *QuantumEspresso* [@Giannozzi2009], *Abinit* [@Gonze2020], *FHI-aims* [@Blum2009], and several others.


# Statement of need
To date, there is no software solution that allows for the seamless bridging and interlinking of _lattice dynamics_ and  _molecular dynamics_ based approaches, despite the fact that actual material science studies can profit in accuracy and efficiency by exploiting both approaches. For instance,
potential use cases include 
accelerating _MD_ calculations by starting from harmonic equilibrium configurations [@West2006], 
analyzing _MD_ simulations in terms of harmonic phonons [@Turney2009], 
investigating the range of validity of the perturbative expansion used in _lattice dynamics_ [@Knoop2020], 
and overcoming finite-size and finite-time effects in _ab initio_ Green Kubo simulations of the thermal conductivity [@Carbogno2016]. 
Given the wide opportunities for application, the aspect of _integration_, i.e., the ability to utilize different methodologies from distinct codes in a flexible fashion using a consistent user interface, is paramount. In particular, this is a prerequisite for automatizing these workflows to enable hierarchical high-throughput screening of whole material classes in a systematic fashion. For example, such a workflow would start from geometry optimizations followed by a study of harmonic properties for many materials, so to single out candidate materials for more involved, fully anharmonic aiMD simulation techniques. Along these lines, let us mention that providing descriptive input and output files is a prerequisite for sharing raw data and results in a transparent and interpretable way in the spirit of open science and the FAIR Principles [@Draxl2018]. On top of that, tracking the provenance [@AiiDA] across different codes facilitates the repurposing and analysis of the obtained data.

# Summary

_FHI-vibes_ is a _python_ package that allows for such an integrated workflow. It uses the _Atomistic Simulation Environment (ASE)_ [@Larsen2017] as a backend in order to represent materials and to connect to various first-principles codes. Via _ASE_, _FHI-vibes_ provides a flexible framework for geometry optimization and MD, and connects to external codes like _spglib_ [@Togo2018], _phonopy_ [@Togo2015], _phono3py_ [@Togo2015b], and _hiphive_ [@Eriksson2019] that implement lattice dynamics techniques based on the harmonic approximation. For all these tasks, _FHI-vibes_ provides defined input files and a command line interface to set up and run calculations on local machines and clusters using the _slurm_ submission system. The output is organized in self-contained and descriptive output files that enable a straightforward exchange of the data obtained with different methodologies.
For advanced analysis, it provides an API fully compatible  with _ASE_ as well as _numpy_ [@Walt2011], _pandas_ [@McKinney2011], and _xarray_ [@Hoyer2017]; several user-friendly utilities allow to perform the most common postprocessing tasks within the command-line interface, such as providing comprehensive summaries of MD simulations or phonon calculations.

_FHI-vibes_ provides a connection to *FireWorks* [@Jain2015], a workflow management system for running simulation workflows on extensive sets of materials in high-throughput fashion. _FHI-vibes_ is tightly integrated with *FHI-aims* [@Blum2009] to perform energy and force calculations, but extending the functionality to any calculator available via *ASE* is straightforward.

_FHI-vibes_ was used to produce the results in [@Knoop2020].

## Features

To facilitate the scientific studies described in the [statement of need](#statement-of-need), _FHI-vibes 1.0_ offers the following main features:

- Free and symmetry-constrained geometry optimization, 

- harmonic phonon calculations, 

- molecular dynamics simulations, 

- harmonic sampling, and 

- anharmonicity quantification. 

An extensive user guide including tutorials and a reference documentation for these features is available at [`vibes.fhi-berlin.mpg.de`](http://vibes.fhi-berlin.mpg.de/). As demonstrated in a dedicated tutorial, the tasks can be easily combined and tailored to define workflows for high-throughput screening of material space.

The codebase and user interface of *FHI-vibes* are designed as a modular framework such that more advanced features and workflows are straightforward to add in the future.

# Acknowledgements
The authors would like to thank Roman Kempt and Marcel Hülsberg for testing and providing valuable feedback. F.K. would like to thank Marcel Langer and Zhenkun Yuan for feedback and Ask Hjorth Larsen for valuable discussions. T.P. would like to thank the Alexander von Humboldt Foundation for their support through the Alexander von Humboldt Postdoctoral Fellowship Program. This project was supported by TEC1p (the European Research Council (ERC) Horizon 2020 research and innovation programme, grant agreement No. 740233), BigMax (the Max Planck Society’s Research Network on Big-Data-Driven Materials-Science), and the NOMAD pillar of the FAIR-DI e.V. association.

# References
# Test

* Supercells with 100 atoms
* created with  `make_supercell geometry.in -n 100 --deviation 0`
lammps Green Kubo Ar
===

[reference heat flux](https://lammps.sandia.gov/doc/compute_heat_flux.html)

## Background

$\begin{aligned} \mathbf{J} &=\frac{1}{V}\left[\sum_{i} e_{i} \mathbf{v}_{i}-\sum_{i} \mathbf{S}_{i} \mathbf{v}_{i}\right] \\ &=\frac{1}{V}\left[\sum_{i} e_{i} \mathbf{v}_{i}+\sum_{i<j}\left(\mathbf{f}_{i j} \cdot \mathbf{v}_{j}\right) \mathbf{x}_{i j}\right] \\ &=\frac{1}{V}\left[\sum_{i} e_{i} \mathbf{v}_{i}+\frac{1}{2} \sum_{i<j}\left(\mathbf{f}_{i j} \cdot\left(\mathbf{v}_{i}+\mathbf{v}_{j}\right)\right) \mathbf{x}_{i j}\right] \end{aligned}$

$\kappa=\frac{V}{k_{B} T^{2}} \int_{0}^{\infty}\left\langle J_{x}(0) J_{x}(t)\right\rangle d t=\frac{V}{3 k_{B} T^{2}} \int_{0}^{\infty}\langle\mathbf{J}(0) \cdot \mathbf{J}(t)\rangle d t$

## Flags

`compute ID group-ID heat/flux ke-ID pe-ID stress-ID`

This compute calculates 

- 6 quantities and stores them in a 6-component vector. The 

- first 3 components are the x, y, z components of the full heat flux vector, i.e. (Jx, Jy, Jz). 

- The next 3 components are the x, y, z components of just the convective portion of the flux, i.e. the first term in the equation for J above.

The heat flux can be output every so many timesteps (e.g. via the *thermo_style* custom command). Then as a post-processing operation, an **auto-correlation** can be performed, its **integral estimated**, and the **Green-Kubo formula above evaluated**.

The _fix ave/correlate command_ can calculate the auto-correlation. The *trap() function* in the variable command can calculate the integral.

## Output info

This compute calculates a 

- global vector of length 6 
  
  - (total heat flux vector, followed by convective heat flux vector), which 

- can be accessed by indices 1-6. These values

- can be used by any command that uses global vector values from a compute as input. See the Howto output doc page for an overview of LAMMPS output options.

The vector values calculated by this compute are 

- “extensive”, meaning they scale with the number of atoms in the simulation. They 

- can be divided by the appropriate volume to get a flux, which would then be an “intensive” value, meaning independent of the number of atoms in the simulation. Note that 

- if the compute is “all”, then the *appropriate volume to divide by is the simulation box volume*. However, 

- if a sub-group is used, it should be the volume containing those atoms.

#### Thermo command

[reference thermo](https://lammps.sandia.gov/doc/thermo.html)

```bash
thermo N
```

Compute and print thermodynamic info every _N_ timesteps

```bash
thermo_style style args
```

### fix ave/correlate command

[reference ave/correlate](https://lammps.sandia.gov/doc/fix_ave_correlate.html)

```bash
fix ID group-ID ave/correlate Nevery Nrepeat Nfreq value1 ...
```

The Nevery, Nrepeat, and Nfreq arguments specify on what timesteps the input values will be used to calculate correlation data.

| argument  | meaning                                                  |
| --------- | -------------------------------------------------------- |
| `Nevery`  | use input values every this many timesteps               |
| `Nrepeat` | # of correlation time windows to accumulate              |
| `Nfreq`   | calculate time window averages every this many timesteps |

#### Example

For example, if **Nevery=10**, **Nrepeat=200**, and **Nfreq=2000**, then 

- values on timesteps 0,10,20,…,2000 will be used to compute the final averages on timestep 2000 
- 200 averages will be computed: Cij(0), Cij(10), Cij(20), ..., and Cij(2000). 
- Cij(30) on timestep 2000 will be the average of 199 samples, namely 
  - Vi(0)*Vj(30), Vi(10)*Vj(40), …, Vi(1980)*Vj(1990), Vi(1990)*Vj(2000)
- Cij(30) on timestep 6000 will be the average of 599 samples, namely 
  - Vi(0)*Vj(30), Vi(10)*Vj(40), …, Vi(5980)*Vj(5990), Vi(5990)*Vj(6000)
- and so on and so on

**ave**

If the *ave* setting is running, then the accumulation is never zeroed. Thus the output of correlation data at any timestep is the 

- average over samples accumulated every *Nevery* steps since the fix was defined.

it can only be restarted by deleting the fix via the unfix command, or by re-defining the fix by re-specifying it.

## Parameter

| Parameter/Variable | Value                                                       |
| ------------------ | ----------------------------------------------------------- |
| `vol`              | _a_ = 5.476, 4x4x4 Box -> 9943.92 Å**3                      |
| `scale`            | `28.7289125255705` [=`convert` * `dt`  / (`V * kB * T**2`)] |
| `convert`          | `4.83166430676946e-16`                                      |
| `kB`               | `1.3806504e-23`                                             |
| `dt`               | `4 * 10` (because step `s = 10`)                            |

**remarks**
volume: _a_ = 5.476, 4x4x4 Box -> 9943.92 Å**3

### Sample input script

```bash
# Sample LAMMPS input script for thermal conductivity of solid Ar

units       real
variable    T equal 70
variable    V equal vol
variable    dt equal 4.0
variable    s equal 10      # sample interval     = Nevery
variable    p equal 200     # correlation length  = Nrepeat
variable    d equal $p*$s   # dump interval       = Nfreq

# convert from LAMMPS real units to SI

variable    kB equal 1.3806504e-23    # [J/K] Boltzmann
variable    kCal2J equal 4186.0/6.02214e23
variable    A2m equal 1.0e-10
variable    fs2s equal 1.0e-15
variable    convert equal ${kCal2J}*${kCal2J}/${fs2s}/${A2m}

# setup problem

dimension    3
boundary     p p p
lattice      fcc 5.376 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region       box block 0 4 0 4 0 4
create_box   1 box
create_atoms 1 box
mass         1 39.948
pair_style   lj/cut 13.0
pair_coeff   * * 0.2381 3.405
timestep     ${dt}
thermo       $d

# equilibration and thermalization

velocity     all create $T 102486 mom yes rot yes dist gaussian
fix          NVT all nvt temp $T $T 10 drag 0.2
run          8000

# thermal conductivity calculation, switch to NVE if desired

#unfix       NVT
#fix         NVE all nve

reset_timestep 0
compute      myKE all ke/atom
compute      myPE all pe/atom
compute      myStress all stress/atom NULL virial
compute      flux all heat/flux myKE myPE myStress
variable     Jx equal c_flux[1]/vol
variable     Jy equal c_flux[2]/vol
variable     Jz equal c_flux[3]/vol
fix          JJ all ave/correlate $s $p $d &
             c_flux[1] c_flux[2] c_flux[3] type auto file J0Jt.dat ave running
variable     scale equal ${convert}/${kB}/$T/$T/$V*$s*${dt}
variable     k11 equal trap(f_JJ[3])*${scale}
variable     k22 equal trap(f_JJ[4])*${scale}
variable     k33 equal trap(f_JJ[5])*${scale}
thermo_style custom step temp v_Jx v_Jy v_Jz v_k11 v_k22 v_k33
run          100000
variable     k equal (v_k11+v_k22+v_k33)/3.0
variable     ndens equal count(all)/vol
print        "average conductivity: $k[W/mK] @ $T K, ${ndens}
```
Trajectory
===

## Indices and Labels

- atom labels: Capital letters starting from $`I`$: $`I, J, K, L, ...`$

- Coordinate labels: small letters starting form $`a`$: $`a, b, c, d, ...`$
  
  - instead of the greek letters $`\alpha, \beta, \gamma, \delta, ...`$

- Time label: just `time`

## Examples

| Observable                                           | dimensions                   |
| ---------------------------------------------------- | ---------------------------- |
| temperature, pressure, energy                        | [`time`]                     |
| stress                                               | [`time`, `a`, `b`]           |
| positions, velocities, forces, heat flux             | [`time`, `I`, `a`]           |
| stresses                                             | [`time`, `I`, `a`, `b`]      |
| heat flux autocorrelation function, cumulative kappa | [`time`, `I`, `J`, `a`, `b`] |
Green Kubo
===

## Thermal Conductivity

$`\kappa^{\alpha, \beta} \equiv \sum_{i, j} \kappa_{i, j}^{\alpha, \beta} = \frac{1}{k_\text{b} T^2 V} \sum_{i, j} \int \left\langle J_i (\tau) J_j \right\rangle ~ \text{d} \tau`$
vibes.slurm
===

Heavily inspired by [sirmarcels](https://gitlab.com/sirmarcel) [jobro](https://gitlab.com/sirmarcel/jobro) project.

Maybe we join forces one day :fist:

```
args = (
    "name",
    "logfile",
    "mail_type",
    "mail_address",
    "nodes",
    "cores",
    "timeout",
    "queue",
    "command",
    "tag",
)
```

vibes run
===

type `vibes run aims`, `vibes run phonopy`, or `vibes run md` to run one of the example input files

- `aims.in`

- `phonopy.in`

- `md.in`
# References

## How to cite `FHI-vibes`?

[Knoop et al., (2020). FHI-vibes: Ab Initio Vibrational Simulations. Journal of Open Source Software, 5(56), 2671](https://doi.org/10.21105/joss.02671)

```
@article{Knoop2020,
  doi = {10.21105/joss.02671},
  url = {https://doi.org/10.21105/joss.02671},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {56},
  pages = {2671},
  author = {Florian Knoop and Thomas A. R. Purcell and Matthias Scheffler and Christian Carbogno},
  title = {FHI-vibes: _Ab Initio_ Vibrational Simulations},
  journal = {Journal of Open Source Software}
}
```

## Work that was performed using `FHI-vibes`

### "Anharmonicity measure for materials"
[F. Knoop, T.A.R. Purcell, M. Scheffler, and C. Carbogno, Phys. Rev. Materials **4**, 083809 (2020)](https://doi.org/10.1103/PhysRevMaterials.4.083809), [arXiv:2006.14672](https://arxiv.org/abs/2006.14672)

```
@article{PhysRevMaterials.4.083809,
title = {Anharmonicity measure for materials},
author = {Knoop, Florian and Purcell, Thomas A. R. and Scheffler, Matthias and Carbogno, Christian},
journal = {Phys. Rev. Materials},
volume = {4},
issue = {8},
pages = {083809},
numpages = {12},
year = {2020},
month = {Aug},
publisher = {American Physical Society},
doi = {10.1103/PhysRevMaterials.4.083809},
url = {https://link.aps.org/doi/10.1103/PhysRevMaterials.4.083809}
}
```
FHI-vibes
===

Welcome to `FHI-vibes`, a `python` package for calculating, analyzing, and understanding the vibrational properties of solids from first principles. `FHI-vibes` is intended to seamlessly bridge between the harmonic approximation and fully anharmonic molecular dynamics simulations. To this end, `FHI-vibes` builds on several [existing packages](https://vibes-developers.gitlab.io/vibes/Credits/) and interfaces them in a consistent and user-friendly fashion. 

In the documentation and tutorials, knowledge of first-principles electronic-structure theory as well as proficiency with _ab initio_ codes such as [FHI-aims](https://aimsclub.fhi-berlin.mpg.de/) and high-performance computing are assumed. Additional experience with Python, the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/), or [Phonopy](https://atztogo.github.io/phonopy/) is helpful, but not needed.

`FHI-vibes` provides the following features:

- Geometry optimization via [ASE](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#module-ase.optimize),
- harmonic phonon calculations via [Phonopy](https://atztogo.github.io/phonopy/),
- molecular dynamics simulations in [NVE](https://wiki.fysik.dtu.dk/ase/ase/md.html#constant-nve-simulations-the-microcanonical-ensemble), [NVT](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin), and [NPT](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.nptberendsen) ensembles,
- [harmonic sampling](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.115504), and
- [anharmonicity quantification](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.4.083809).

Most of the functionality is high-throughput ready via [fireworks](https://materialsproject.github.io/fireworks/#).

## Overview

- [Installation](https://vibes-developers.gitlab.io/vibes/Installation)
- [Tutorial](https://vibes-developers.gitlab.io/vibes/Tutorial/0_intro)
- [Documentation](https://vibes-developers.gitlab.io/vibes/Documentation/0_intro)
- [Credits](https://vibes-developers.gitlab.io/vibes/Credits)
- [References](https://vibes-developers.gitlab.io/vibes/References)


## News

- `FHI-vibes` got [published in JOSS](https://joss.theoj.org/papers/10.21105/joss.02671)!
- [Our anharmonicity measure got published!](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.4.083809)
- [… the best is yet to come.](https://www.youtube.com/watch?v=B-Jq26BCwDs)

## Changelog

#### v1.0.3

- update dependencies to allow `phonopy` versions up to 2.8

#### v1.0.2

- First official release after passing the [JOSS review](https://github.com/openjournals/joss-reviews/issues/2671).
- Several additions to the documentation.

#### v1.0.0a10

- Enable conversion of trajectories to `ase.io.Trajectory` files for viewing with ASE [(!37)](https://gitlab.com/vibes-developers/vibes/-/merge_requests/37)
- Important fix for running NPT dynamics [(!36)](https://gitlab.com/vibes-developers/vibes/-/merge_requests/36)
- We have a changelog now!
# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the maintainers of this repository. This will make life easier for everyone.

## Report Issues

Please use the [issue tracker](https://gitlab.com/vibes-developers/vibes/-/issues) to report issues. Please try to answer these questions:

- Has this issue been discussed before? Please have a quick look at the existing issues. If not:
- What is the issue? What is the expected behavior?
- Is the problem reproducible? Please provide a _minimal_ example.


## Contribute Code via Merge Request

In order to contribute code to `FHI-vibes`, please follow the usual steps for [preparing and creating a merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html). A few remarks regarding our guidelines for code and code style:

- We use [black](https://black.readthedocs.io/en/stable/) with default settings and [isort](https://pycqa.github.io/isort/) for formatting the code. The settings for `isort` are included in `setup.cfg`.
- Please _document_ and _test_ your changes. Tests are found in `vibes/tests` and written with [pytest](https://docs.pytest.org/en/stable/). 
- Please use [google-type docstrings](https://google.github.io/styleguide/pyguide.html) for your functions. Optionally you can use type hints, but we currently don't enforce this.
- We loosely keep track of code coverage, please try not to decrease coverage when contributing new code.# Credits

`FHI-vibes` would not be possible without the following packages:

- The [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
- [FHI-aims: FHI _ab initio_ molecular simulations](https://aimsclub.fhi-berlin.mpg.de/)
- [Phonopy](https://atztogo.github.io/phonopy/) and [Phono3py](https://atztogo.github.io/phono3py/)
- [fireworks](https://materialsproject.github.io/fireworks/#)
- [hiPhive — High-order force constants for the masses](https://hiphive.materialsmodeling.org/index.html)
- [The SciPy  Stack](https://www.scipy.org/)
- [mkdocs](https://www.mkdocs.org/) and [mkdocs-material](https://squidfunk.github.io/mkdocs-material/)

### How to cite these packages:

Please make sure to give credit to the right people when using `FHI-vibes`:

- [How to cite ASE](https://wiki.fysik.dtu.dk/ase/faq.html#how-should-i-cite-ase)
- [How to cite FHI-aims](https://aimsclub.fhi-berlin.mpg.de/aims_publications.php)
- [How to cite phonopy](https://phonopy.github.io/phonopy/citation.html)
- [How to cite phono3py](https://phonopy.github.io/phono3py/citation.html)
- [How to cite fireworks](https://materialsproject.github.io/fireworks/#citing-fireworks)
- [How to cite hiphive](https://hiphive.materialsmodeling.org/credits.html)
- [How to cite FHI-vibes](References.md#how-to-cite-fhi-vibes)# Installation

## Prerequisites

- A working`python3.6+` environment, e.g., provided by [anaconda](https://docs.conda.io/en/latest/miniconda.html).

- A working `fortran` compiler, e.g., obtained by:
    - `apt-get install gfortran` in Debian-derived systems, or
    - `conda install -c conda-forge fortran-compiler` when `conda` is used.

- If you want to use `FHI-aims` for running _ab initio_ calculations, make sure you have a recent version that supports the iPi socket communication (this is the default for any version newer than the `200112_2` release when using the [CMake build system](https://aims-git.rz-berlin.mpg.de/aims/FHIaims/-/wikis/CMake-Tutorial)).


## Install `vibes`

`FHI-vibes` can be installed simply via pip:

```bash
pip install --user fhi-vibes
```

The `--user` option makes sure that the installation occurs in your homefolder under `~/.local/bin`, as typically
necessary on computing clusters. Please make sure that the `~/.local/bin` folder is listed in your `PATH`.

**If you run into problems, please have a look at our [troubleshooting section.](#Troubleshooting)**

## Configuration

Configure `vibes` by creating a `~/.vibesrc` configuration file template in the home directory. To this end, first run

```
vibes template configuration vibes > ~/.vibesrc
```

and edit the configuration file as described below:

### `basissetloc`

The `basissetloc` should point to the folder containing FHI-aims' species defaults, e.g., the `/path/to/FHIaims/species_defaults` folder.

### `aims_command`

The `aims_command` should be an executable script that takes care of setting up the environment and then running FHI-aims, for example a file called `/path/to/FHIaims/run_aims.sh`  that looks roughly like this (depends on you system!):

```
#!/bin/bash -l

ulimit -s unlimited
export OMP_NUM_THREADS=1

module purge
module load intel impi mkl
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MKLROOT}/lib/intel64/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${INTEL_HOME}/lib/intel64/

srun /path/to/FHIaims/build/aims.x
```

The script `run_aims.sh` has to be marked as executable, e.g., by running `chmod +x /path/to/FHIaims/run_aims.sh`.

**You're now good to go!**

## Autocompletion

To activate autocompletion of `vibes` subcommands, add this to your `.bashrc`:

```bash
eval "$(_VIBES_COMPLETE=source vibes)"
```

and source it.

If you use the `fishshell`, add a file `~/.config/fish/completions/vibes.fish` containing

```bash
eval (env _VIBES_COMPLETE=source-fish vibes)
```



## Troubleshooting

- `ModuleNotFoundError: No module named 'importlib_resources'`
    - Solution: `pip install importlib_resources dataclasses`
- `RuntimeError: Click will abort further execution because Python 3 was configured to use ASCII as encoding for the environment. Consult https://click.palletsprojects.com/python3/ for mitigation steps`
    - Solution:  `export LC_ALL=C.UTF-8 ; export LANG=C.UTF-8`
- `-bash: vibes: command not found`
    - Solution: `export PATH=$PATH:~/.local/bin`
- `ImportError: numpy.core.multiarray failed to import`
    - Solution: `pip install numpy -U` (or `conda update numpy` if you use conda)
- Various version conflicts
    - Consider using a [virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), e.g., via `conda create -n py38 -c anaconda python=3.8 numpy scipy`

- `/tmp/pip-build-env-xak_2vfy/overlay/lib/python3.7/site-packages/numpy/core/_multiarray_umath.cpython-37m-x86_64-linux-gnu.so: failed to map segment from shared object: Operation not permitted`
    - This might happen on HPC systems with limited access rights. The solution is to provide a writable `tmp` folder, e.g. via `mkdir ~/condatmp && export TMPDIR=~/condatmp/`

If your problem is not listed here, please [file an issue in our issue tracker](https://gitlab.com/vibes-developers/vibes/-/issues).
Output Files
===

## Trajectories

### `trajectory.son`

This file contains metadata and calculation results for a set of related calculations that, e.g., come from a geometry optimization, a phonopy  calculation, or a molecular dynamics run. The file format is [`son`](https://flokno.github.io/son/), a slight extension to `json` allowing to add data sequentially.

??? example "Example `trajectory.son`"
    ```
    {"MD": {
      "type": "molecular-dynamics",
      ...},
    "calculator": {
      "calculator": "Aims",
      "calculator_parameters": {
        "xc": "pbesol",
        "k_grid": [2, 2, 2],
        ...}},
    "atoms": {
      "cell": 
        [[ 8.33141234000000e+00, -8.33141234000000e+00,  0.00000000000000e+00],
         [ 8.33141234000000e+00,  8.33141234000000e+00,  0.00000000000000e+00],
         [ 0.00000000000000e+00,  0.00000000000000e+00,  1.24971185100000e+01]],
      "positions": 
        [[ 0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00],
         [ 4.16570617000000e+00,  0.00000000000000e+00,  0.00000000000000e+00],
       …}
    “primitive”: { …
    …}
    ===
    {"atoms": {
      "info": {
        "nsteps": 0,
        "dt":  4.91134739423203e-01,
        "aims_uuid": "D985353A-F8FD-4635-A939-E129A7E6E146"},
      "positions": 
        [[ 0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00],
         [ 4.16570617000000e+00,  0.00000000000000e+00,  0.00000000000000e+00], …}
    ---
    ...
    ```

### `trajectory.nc`

A [`NetCDF`](https://www.unidata.ucar.edu/software/netcdf/) file containing and [`xarray.Dataset`](http://xarray.pydata.org/en/stable/io.html?highlight=netcdf#netcdf) with post-processed data

??? example "Load `trajectory.nc`"
    ```python
    >>> import xarray as xr
    >>> 
    >>> ds = xr.open_dataset('trajectory.nc')
    >>> 
    >>> print(ds)
    <xarray.Dataset>
    Dimensions:              (I: 180, a: 3, b: 3, time: 10001)
    Coordinates:
      * time                 (time) float64 0.0 2.0 4.0 6.0 ... 2e+04 2e+04 2e+04
    Dimensions without coordinates: I, a, b
    Data variables:
        positions            (time, I, a) float64 ...
        displacements        (time, I, a) float64 ...
        velocities           (time, I, a) float64 ...
        momenta              (time, I, a) float64 ...
        forces               (time, I, a) float64 ...
        energy_kinetic       (time) float64 ...
        energy_potential     (time) float64 ...
        stress               (time, a, b) float64 ...
        stress_kinetic       (time, a, b) float64 ...
        stress_potential     (time, a, b) float64 ...
        temperature          (time) float64 ...
        cell                 (time, a, b) float64 ...
        positions_reference  (I, a) float64 ...
        lattice_reference    (a, b) float64 ...
        pressure             (time) float64 ...
        pressure_kinetic     (time) float64 ...
        pressure_potential   (time) float64 ...
    Attributes:
        name:             trajectory
        system_name:      O3Sb2
        natoms:           180
        time_unit:        fs
        timestep:         1.9999999999999978
        nsteps:           10000
        symbols:          ['O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',...
        masses:           [ 15.999  15.999  15.999  15.999  15.999  15.999  15.99...
        atoms_reference:  {"pbc": [true, true, true],\n"cell": \n[[ 1.50410061436...
        atoms_primitive:  {"pbc": [true, true, true],\n"cell": \n[[ 5.01366871456...
        atoms_supercell:  {"pbc": [true, true, true],\n"cell": \n[[ 1.50410061436...
        volume:           3051.2387320953862
        raw_metadata:     {"MD": {\n  "type": "molecular-dynamics",\n  "md-type":...
        hash:             097714462c68b9f8cbdf08e6a29c0bfda7922c01
    ```



## Numerical Data

### `.dat` files

Files with `.dat` suffix are 1D or 2D arrays that can be read with [`numpy.loadtxt`](https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html):

??? example "Example `frequencies.dat`"
    ```
    0.000000000000000000e+00
    0.000000000000000000e+00
    0.000000000000000000e+00
    9.328537621518436795e-01
    1.319254442145902706e+00
    1.615750112078766065e+00
    1.865707524303687359e+00
    2.085924425237942970e+00
    2.285015721907639907e+00
    2.468099064244776208e+00
    2.638508884291805412e+00
    2.798561286455531594e+00
    ```

??? example "Load `frequencies.dat`"
    ```python
    >>> import numpy as np
    >>> 
    >>> data = np.loadtxt('frequencies.dat')
    >>> 
    >>> print(data)
    [0.         0.         0.         0.93285376 1.31925444 1.61575011
     1.86570752 2.08592443 2.28501572 2.46809906 2.63850888 2.79856129]
    ```
    
### `.csv` files
Files `.csv` suffix are standard [`comma-separated values`](https://en.wikipedia.org/wiki/Comma-separated_values) that can be parsed, e.g., with `pandas`.

??? example "Example `trajectory.csv`"
    ```
    time,temperature,energy_kinetic,energy_potential,pressure_kinetic,pressure_potential,pressure
    0.0,286.2001210048815,6.658958660176626,-13092419.8403197,0.0014549191863471355,-0.0001421961672687167,0.0013127230190784188
    1.9999999999999978,289.5694460043815,6.7373520438187935,-13092419.9182694,0.0014720473956910023,-0.00012147676349475005,0.0013505706321962523
    3.9999999999999956,292.46310420107915,6.80467818694143,-13092419.9839958,0.0014867575181546967,,
    5.999999999999994,294.9143609337456,6.86171106725956,-13092420.0383063,0.0014992186605137458,,
    ```

??? example  "Load `trajectory.csv`"
    ```python
    >>> import pandas as pd
    >>> 
    >>> df = pd.read_csv('trajectory.csv')
    >>> 
    >>> print(df)
              time  temperature  energy_kinetic  energy_potential  pressure_kinetic  pressure_potential  pressure
    0          0.0   286.200121        6.658959     -1.309242e+07          0.001455           -0.000142  0.001313
    1          2.0   289.569446        6.737352     -1.309242e+07          0.001472           -0.000121  0.001351
    2          4.0   292.463104        6.804678     -1.309242e+07          0.001487                 NaN       NaN
    
    [3 rows x 7 columns]
    ```

### `.json` files
These are plain [`JSON`](https://www.json.org/) files that  can be parsed with the [python builtin `json` module](https://docs.python.org/3/library/json.html)

??? example "Example `md_describe.json`"
    ```
    {
     "time": {
      "count": 10001.0,
      "mean": 9999.999999999987,
      "std": 5774.368710084239,
      "min": 0.0,
      "25%": 4999.999999999994,
      "50%": 9999.999999999987,
      "75%": 14999.999999999984,
      "max": 19999.999999999975
     },
     "temperature": {
      "count": 10001.0,
      "mean": 309.14540993328325,
      "std": 14.54050396748833,
      "min": 257.7711666913883,
      "25%": 299.2299218510854,
      "50%": 309.2638007517673,
      "75%": 319.0663317618751,
      "max": 355.05955146050223
     },
     "energy_kinetic": {
      "count": 10001.0,
      "mean": 7.19282192299974,
      "std": 0.3383108800851539,
      "min": 5.997508095931853,
      "25%": 6.962120325100335,
      "50%": 7.1955764975375205,
      "75%": 7.4236499467457895,
      "max": 8.261096699662168
     },
     "energy_potential": {
      "count": 10001.0,
      "mean": -13092416.82839676,
      "std": 2.2319201670206765,
      "min": -13092420.901731301,
      "25%": -13092419.756083699,
      "50%": -13092415.5125696,
      "75%": -13092415.1902583,
      "max": -13092414.1603322
     },
     "pressure_kinetic": {
      "count": 10001.0,
      "mean": 0.0015715632359058732,
      "std": 7.391771228878895e-05,
      "min": 0.00131039852390555,
      "25%": 0.001521157129150454,
      "50%": 0.0015721650842653186,
      "75%": 0.001621996965507344,
      "max": 0.0018049711226602926
     },
     "pressure_potential": {
      "count": 1032.0,
      "mean": 0.00027321549417593297,
      "std": 0.002184381533766492,
      "min": -0.006481330360441526,
      "25%": -0.001147865102920749,
      "50%": 0.00023375660603530354,
      "75%": 0.001799169124854924,
      "max": 0.00778187950662052
     },
     "pressure": {
      "count": 1032.0,
      "mean": 0.001844596382543767,
      "std": 0.0021801900840728618,
      "min": -0.004976051568943327,
      "25%": 0.00041292059984573813,
      "50%": 0.0018072058785148773,
      "75%": 0.0033777412146081182,
      "max": 0.009323241222005936
     }
    }
    ```
    
??? example "Load `md_describe.json`"
    ```python
    >>> import json
    >>> 
    >>> data = json.load(open('md_describe.json'))
    >>> 
    >>> pprint(data)
    {'energy_kinetic': {'25%': 6.962120325100335,
                        '50%': 7.1955764975375205,
                        '75%': 7.4236499467457895,
                        'count': 10001.0,
                        'max': 8.261096699662168,
                        'mean': 7.19282192299974,
                        'min': 5.997508095931853,
                        'std': 0.3383108800851539},
     'energy_potential': {'25%': -13092419.756083699,
                          '50%': -13092415.5125696,
                          '75%': -13092415.1902583,
                          'count': 10001.0,
                          'max': -13092414.1603322,
                          'mean': -13092416.82839676,
                          'min': -13092420.901731301,
                          'std': 2.2319201670206765},
     'pressure': {'25%': 0.00041292059984573813,
                  '50%': 0.0018072058785148773,
                  '75%': 0.0033777412146081182,
                  'count': 1032.0,
                  'max': 0.009323241222005936,
                  'mean': 0.001844596382543767,
                  'min': -0.004976051568943327,
                  'std': 0.0021801900840728618},
     'pressure_kinetic': {'25%': 0.001521157129150454,
                          '50%': 0.0015721650842653186,
                          '75%': 0.001621996965507344,
                          'count': 10001.0,
                          'max': 0.0018049711226602926,
                          'mean': 0.0015715632359058732,
                          'min': 0.00131039852390555,
                          'std': 7.391771228878895e-05},
     'pressure_potential': {'25%': -0.001147865102920749,
                            '50%': 0.00023375660603530354,
                            '75%': 0.001799169124854924,
                            'count': 1032.0,
                            'max': 0.00778187950662052,
                            'mean': 0.00027321549417593297,
                            'min': -0.006481330360441526,
                            'std': 0.002184381533766492},
     'temperature': {'25%': 299.2299218510854,
                     '50%': 309.2638007517673,
                     '75%': 319.0663317618751,
                     'count': 10001.0,
                     'max': 355.05955146050223,
                     'mean': 309.14540993328325,
                     'min': 257.7711666913883,
                     'std': 14.54050396748833},
     'time': {'25%': 4999.999999999994,
              '50%': 9999.999999999987,
              '75%': 14999.999999999984,
              'count': 10001.0,
              'max': 19999.999999999975,
              'mean': 9999.999999999987,
              'min': 0.0,
              'std': 5774.368710084239}}
    ```



## Force Constants

### `FORCE_CONSTANTS`

These are force constants in the [`phonopy` format](https://phonopy.github.io/phonopy/input-files.html?highlight=force_const#force-constants-and-force-constants-hdf5) in the compact form `(n_primitive, n_supercell, 3, 3)`. They can be parsed with `phonopy.file_IO.parse_FORCE_CONSTANTS`.

### `FORCE_CONSTANTS_remapped`

These are force constants mapped to `(3 * n_supercell, 3 * n_supercell)` shape. They can be parsed with `numpy.loadtxt` similar to [`.dat.` files](#dat-files).

# Calculator Setup

FHI-vibes can set up any [ASE calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#module-ase.calculators) for performing a calculation by providing the  calculator class `name` and the respective `parameters` in the input file. If a `module` is specified, `vibes` will attempt to import the calculator from that module instead of `ase`. This can be used to work with custom calculators that are not (yet) included in `ase`.


## Example

```
...
[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4
...
```

This would set up a [Lennard Jones calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/others.html#lennard-jones) with a `sigma` value of 3.4 and default parameters otherwise.

For a non-`ase` calculator, this would be:

```
...
[calculator]
name:                          MyCalculator
module:                        mymodule

[calculator.parameters]
a:                             1.23
...
```

`vibes` will then attempt to import `MyCalculator` from `mymodule` and instantiate it with `a=1.23`. 

## Sections

### `[calculator]`

This section specifies which `ase.Calculator` should be set up and how.

#### `name`

The name of the [ASE calculator class name](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#supported-calculators).

Note that for non-`ase` calculators, `name` must be spelled identically to the class name in the module, i.e. typically `CamelCase`.

#### `module` (optional)

If specified, `vibes` will run `from module import name` to obtain the calculator class, instead of importing it from `ase`. 

### `[calculator.parameters]`

These are the keywords used 1:1 to set up the ASE calculator:

```python
cls = get_calculator_class(settings.calculator.get("name"))

calculator = cls(**settings.calculator.get("parameters"))
```

## Options for `FHI-aims`

FHI-vibes is most tightly integrated with the FHI-aims calculator and provides some extra features for performing _ab initio_ calculations with FHI-aims. A minimal input section to set up an FHI-aims calculator looks like this:

```
[calculator]
name:                          aims

[calculator.parameters]
xc:                            pw-lda

[calculator.kpoints]
density:                       3.5

[calculator.basissets]
default:                       intermediate
fallback:                      light
# system specific
# O:                           tight
# Ga:                          intermediate

[calculator.socketio]
port:                          12345
```

### `[calculator.parameters]`

These keywords correspond one-to-one to the FHI-aims keywords that  are written to `control.in`. Keyword-only arguments like `vdw_correction_hirshfeld` or `use_gpu` should be given with the value `true`:

```
[calculator.parameters]
xc:                            pw-lda
vdw_correction_hirshfeld:      true
use_gpu:                       true
...
```



### `[calculator.kpoints]` (optional)

#### `density`

Instead of giving a `k_grid` explicitly, FHI-vibes can compute `k_grid` such that the density of kpoints does not fall below this value in $\require{mediawiki-texvc} \AA^{-3}$ . This is optional, including `k_grid` in `[calculator.parameters]` is equally valid.

### `[calculator.basissets]`

Specify which basissets to use.

#### `default`

The default basis set to use, can be `light`, `intermediate`,`tight`, or `really_tight`.

#### `fallback`

The fallback option in case the specified basis set could not be found (`intermediate` basis sets are currently not compiled for each element)

#### Species dependent

The basis set can be given per chemical species by including the species and its desired basis set (uncomment, e.g., `O` in the example above.)

### `[calculator.socketio]` (optional)

Set up socket communication via [`SocketIOCalculator`](https://wiki.fysik.dtu.dk/ase/ase/calculators/socketio/socketio.html?highlight=socketio#ase.calculators.socketio.SocketIOCalculator). This has the potential to speed up calculations since a complete restart of FHI-aims after each completed SCF cycle is avoided. This feature is optional but recommended to use when performing calculations for related structures, e.g., during molecular dynamics simulations or phonon calculations.

#### `host`

The IP address to access the socket. Default is `localhost` and will only have to be modified for certain architectures (this will likely be made clear in the system documentation).

#### `port`

The socket port to use.

- `null`: don't use the socket.
- `auto`: Automatically select a port that is not currently in use by the `host` or registered in `/etc/services`
- `1024`-`65535`: If available use this port (if it's not already being used). We recommend using `auto` for all calculations.

#### `unixsocket`

Filename for the unix socket. If this is active TCP/IP socket will not be used (not recommended, but maybe necessary on some systems)!!! info
	An hands-on example for setting up and running a `phonopy` calculation can be found in the [Tutorial](../Tutorial/2_phonopy.md).

vibes supports [phonon calculations with the finite differences method](../Tutorial/2_phonopy_intro.md#phonons-harmonic-vibrations-in-solids) by setting up a `phonopy.in` file. A minimal `phonopy.in` would look like

```fo
[files]
geometry:                      geometry.in

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4

[phonopy]
supercell_matrix:              [2, 2, 2]
```

for performing a phonopy calculation for the structure in `geometry.in` with a Lennard-Jones calculator and a $2 \times 2 \times 2$  supercell.

??? info "Click: Default values for the complete list of supported keywords"
    ```
    [phonopy]
    supercell_matrix:              [1, 1, 1]
    displacement:                  0.01
    is_diagonal:                   False
    is_plusminus:                  auto
    symprec:                       1e-05
    q_mesh:                        [45, 45, 45]
    workdir:                       phonopy
    ```

## Sections

### `[phonopy]`
#### `supercell_matrix`
`list`: A $3 \times 1$ or $9 \times 1$ array specifying the [supercell matrix](../Tutorial/2_phonopy.md#supercell-matrix-supercell_matrix) used for setting up the supercell.

#### `displacement`

`float`: the distance in $\require{mediawiki-texvc} \AA$ used for the finite displacement.

#### `is_diagonal`

`True/False`: corresponds to the  `phonopy` settings tag [`DIAG`](https://phonopy.github.io/phonopy/setting-tags.html#diag)

#### `is_plusminus`

`True/False/auto`: corresponds to the  `phonopy` settings tag [`PM`](https://phonopy.github.io/phonopy/setting-tags.html#pm)

#### `symprec`

`float`: tolerance for symmetry detection

#### `q_mesh`

`list`: the q-points mesh used for postprocessing, e.g., density of states.

#### `workdir`

The working directory for running the calculations and storing output files.
# Slurm Submission

`vibes` can submit jobs on cluster using the `slurm` submission system. The configuration goes via the task input file by adding a `[slurm]` section. A template can be generated via

```
vibes template slurm >> your_input.in
```

??? info "Task input file with `[slurm]` section"

	```
	...
	[slurm]
	name:         test
	tag:          vibes
	mail_type:    all
	mail_address: your@mail.com
	nodes:        1
	cores:        32
	queue:        express
	timeout:      30
	```

## Submit a job

A job is submitted by using `vibes submit` instead of `vibes run`. For example, you would run `vibes submit md md.in` on a cluster to submit the calculation to the queue instead of running `vibes run md md.in` locally. `vibes submit` will submit the job to the cluster according to the specification it finds in the `[slurm]` section, see [below](#the-slurm-section).

The command will do the following:

- Write a `submit.sh` script according to the specifications found in the `[slurm]` section,
- submit the job to the queue via `sbatch submit.sh`, and
- log the time of submission and the job ID to a file called `.submit.log`.

## Restart a job

`vibes` supports autmatic restarts if the job will take longer than the available walltime, which is often the case during [_ab initio_ molecular dynamics simulations](../Tutorial/3_md_ab_initio.md). Restarts can be requested by adding

```
...
[restart]
command: the command to restart the calculation, e.g., `vibes submit md md.in`
```

to your task input file. This will run the command specified here shortly before the walltime is over.

## The `[slurm`] section

### `name`

`str`: The name of the job.

### `tag`

`str`: A common tag added to the job `name` to cluster sets of jobs.

### `mail_type`

`str`:  The [slurm mail type usually specified via `--mail-type`](https://slurm.schedmd.com/sbatch.html)

### `mail_address`

`str`: The mail address that slurm will send notifications to.

### `nodes`

`int`: The number of nodes to be used for the job (job dependent).

### `cores`

`int`: The number of cores per node to be used for the job (**machine dependent!**).

### `queue`

`str`: The name of the queue to be submitted to.

### `timeout`

`int`: Walltime for the job in minute. **Queue and cluster dependent!**Documentation
===

Here we give an overview of the features of `FHI-vibes` and document the basic usage principles.

In particular, you will find

- [Units](units.md)
- [Input files](input_files.md)
- [Output files](output_files.md)
- [Calculator setup](calculator_setup.md)
- Tasks and workflows, i.e.,
    - [Geometry optimization](relaxation.md)
    - [Phonon calculations](phonopy.md)
    - [Molecular dynamics simulations](md.md)
- and [High-Throughput workflows.](../High_Throughput/Documentation/0_overview.md)

Detailed introductions into specific tasks can be found in the [Tutorial](../Tutorial/0_intro.md).



## Design Philosophy

The general design principles are

- calculations should be defined with human-readable input files,
- a calculation should produce a self-contained output file that _includes_ metadata describing the calculation, as well as the calculated results.
- Output files should be easy to parse by a computer.## Units
The units used for input and output files are:

- Energy: ${\rm eV}$
- Length: $\require{mediawiki-texvc} {\rm \AA}$
- Time: ${\rm fs}$
- Frequency: ${\rm THz}$
- Temperature: ${\rm K}$
- Force: ${\rm eV} / {\rm \AA}$
- Stress: ${\rm eV} / {\rm \AA}^3$ # Input File Format

## Geometry input files
`FHI-vibes` uses the `FHI-aims` geometry description format `geometry.in`. A detailed documentation of the file format can be found [here](https://doi.org/10.6084/m9.figshare.12413477.v1).

### Example
An example `geometry.in` for fcc-silicon reads
```
lattice_vector 0.00 2.72 2.72
lattice_vector 2.72 0.00 2.72
lattice_vector 2.72 2.72 0.00

atom_frac 0.00 0.00 0.00 Si
atom_frac 0.25 0.25 0.25 Si
```

## Task input files
For performing a specific task, say, a geometry optimization, `FHI-vibes` employs single input files describing the task, e.g., a `relaxation.in` file describing the optimization according to the [documentation](relaxation.md).

### The `jconfigparser` syntax
The input files are parsed using [`jconfigparser`](https://pypi.org/project/jconfigparser/) . `jconfigparser` is an extension to the `python` [standard library `configparser`](https://docs.python.org/3/library/configparser.html) with the following additions:

- Nested section names separated with `.` ,
- values are parsed by `json`,
- repeated keywords possible,
- [`configparser.ExtendedInterpolation`](https://docs.python.org/3/library/configparser.html#configparser.ExtendedInterpolation) is used per default.

#### Example

An example for an input file for running a geometry optimization:

```
[files]
geometry:                      geometry.in

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4

[relaxation]
driver:                        BFGS
fmax:                          0.001
workdir:                       ${calculator:name}.relaxation

[relaxation.kwargs]
maxstep:                       0.2
```

This file will be parsed to a nested dictionary:

```
settings = {
    "calculator": {"name": "lj", "parameters": {"sigma": 3.4}},
    "files": {"geometry": "geometry.in"},
    "relaxation": {
        "driver": "BFGS",
        "fmax": 0.001,
        "kwargs": {"maxstep": 0.2},
        "workdir": "lj.relaxation",
    },
}
```

### `[files]` Section

This section contains filenames. 

#### `geometry`

`geometry` gives the name of the geometry input file to be used for a calculation:

```python
file = settings.files.get("geometry")

atoms = ase.io.read(file)
```

If there is just one geometry necessary for the given task and it is stored in `geometry.in`, this section can be omitted altogether.

#### `geometries`

Via `geometries`, a wildcard expression for finding geometry files for computation can be given, e.g. `geometries: samples/geometry.in.*` would specifiy to run a calculation for all geometry input files found in the folder `samples`.

```python
files = sorted(glob(settings.files.get("geometries")))

atoms_to_compute = [ase.io.read(file) for file in files]
```


#### `primitive`
Give a reference primitive cell in a file, e.g., `primitive: geometry.in.primitive`

#### `supercell`
Give a reference supercell in a file, e.g., `supercell: geometry.in.supercell`

#### Example
Example for specifying to run a job for the structure in `geometry.in`, while attaching a reference primitive and supercell to the output trajectory:

```
[files]
geometry:                      geometry.in
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell

...
```# Command Line Interface (CLI)

FHI-vibes comes with a command line interface (CLI) for

- creating input files (`vibes  template`),
- informing about input and output files (`vibes info`),
- running calculations (`vibes run` and `vibes submit`),
- processing calculations (`vibes output`), and
- performing several other tasks like converting output files (`vibes utils`).

Practical examples for working with the CLI are found in the [tutorials](../Tutorial/0_intro.md). Each of the the sub-commands has its own `--help` for additional information.


## `vibes template`

Create template files

```
$ vibes template --help

Usage: vibes template [OPTIONS] COMMAND [ARGS]...

  provide template input files for tasks and workflows

Options:
  -h, --help         Show this message and exit.

Commands:
  calculator     Calculator templates: aims, lj
  configuration  Configuration templates: .vibesrc, .fireworksrc
  md             provide template input for MD simulation (default: NVE)
  phonopy        provide template input for phonopy workflow.
  relaxation     provide template input for relaxation workflow.
  slurm          provide template slurm settings
```

The templates are printed to screen and can be piped to a file with `| tee`, `>` or `>>`.

## `vibes info`

```
$ vibes info --help

Usage: vibes info [OPTIONS] COMMAND [ARGS]...

  inform about content of a file

Options:
  -h, --help  Show this message and exit.

Commands:
  anharmonicity   Compute sigmaA for trajectory dataset in FILE
  csv             show contents of csv FILE
  geometry        inform about a structure in a geometry input file
  md              inform about MD simulation in FILE
  netcdf          show contents of netCDF FILE
  phonopy         inform about a phonopy calculation based on the input FILE
  relaxation      summarize geometry optimization in FILE
  settings        write the settings in FILE *including* the configuration
  trajectory      print metadata from trajectory in FILE

```

## `vibes run`

```
vibes run --help

Usage: vibes run [OPTIONS] COMMAND [ARGS]...

  run a vibes workflow

Options:
  -h, --help  Show this message and exit.

Commands:
  md           run an MD simulation from FILE (default: md.in)
  phonopy      run a phonopy calculation from FILE (default: phonopy.in)
  relaxation   run an relaxation from FILE (default: relaxation.in)
  singlepoint  run singlepoint calculations from FILE (default: aims.in)
```



## `vibes submit`

```
$ vibes submit --help

Usage: vibes submit [OPTIONS] COMMAND [ARGS]...

  submit a vibes workflow to slurm

Options:
  --dry
  -h, --help  Show this message and exit.

Commands:
  md           submit MD simulation from FILE (default: md.in)
  phonopy      submit a phonopy calculation from FILE (default: phonopy.in)
  relaxation   submit relaxation from FILE (default: relaxation.in)
  singlepoint  submit singlepoint calculations from FILE (default: aims.in)
```

## `vibes output`

```
$ vibes output --help

Usage: vibes output [OPTIONS] COMMAND [ARGS]...

  produce output of vibes workfow

Options:
  -h, --help  Show this message and exit.

Commands:
  phonopy          perform phonopy postprocess for trajectory in FILE
  trajectory (md)  write trajectory data in FILE to xarray.Dataset
```

## `vibes utils`

```
$ vibes utils --help

Usage: vibes utils [OPTIONS] COMMAND [ARGS]...

  tools and utilities

Options:
  -h, --help  Show this message and exit.

Commands:
  backup                backup FOLDER to TARGET
  create-samples        create samples from geometry in FILENAME
  force-constants (fc)  utils for working with force constants
  geometry              utils for manipulating structures (wrap, refine, etc.)
  hash                  create sha hash for FILE
  make-supercell        create a supercell of desired shape or size
  trajectory            trajectory utils
```

!!! info
	A hands-on example for setting up and running a molecular dynamics run [Tutorial](../Tutorial/3_md_intro.md).

vibes supports running molecular dynamics simulations via [ASE](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md) in NVE, NVT, and NPT ensembles. A minimal `md.in` for running 1000 steps of [Velocity Verlet dynamics](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.verlet) with a Lennard Jones calculator would be

```fo
[files]
geometry:                      geometry.in

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4

[md]
driver:                        VelocityVerlet
timestep:                      1
maxsteps:                      1000
```

??? info "Click: Default values for the complete list of supported keywords"
	NVE dynamics using Velocity Verlet propagation (`vibes template md`)
    ```
    [md]
    driver:                        VelocityVerlet
    timestep:                      1
    maxsteps:                      1000
    compute_stresses:              False
    workdir:                       md

    [md.kwargs]
    logfile:                       md.log
    ```
    NVT ensemble using a [Langevin](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin) thermostat (`vibes template md --nvt`):
    ```
    [md]
    driver:                        Langevin
    timestep:                      1
    maxsteps:                      1000
    compute_stresses:              False
    workdir:                       md
    
    [md.kwargs]
    temperature:                   300
    friction:                      0.02
    logfile:                       md.log
    ```
    NPT ensemble using a [Berendsen](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.nptberendsen) thermostat and barostat (`vibes template md --npt`):
    ```
    [md]
    driver:                        NPTBerendsen
    timestep:                      1
    maxsteps:                      1000
    compute_stresses:              False
    workdir:                       md
    
    [md.kwargs]
    temperature:                   300
    taut:                          500.0
    taup:                          1000.0
    pressure:                      1.01325
    compressibility:               4.57e-05
    logfile:                       md.log
    inhomogeneous:                 False
    ```


## Sections

### `[md]`
#### `driver`
`str`:  The MD algorithm. Either `VelocityVerlet` (NVE ensemble), `Langevin` (NVT ensemble), or `NPTBerendsen` (NPT ensemble).

#### `timestep`

`float`: the integration timestep in fs.

#### `maxsteps`

`int`:  the number of timesteps to perform.

#### `compute_stresses`

`bool/int`: specify whether to compute stress during the MD simulation (more costly in _ab initio_ MD).

#### `workdir`

`str`: the working directory for the MD simulation.

### `[md.kwargs]`
These are keyword arguments  that go straight to the ASE class implementing the MD algorithm, e.g.,

```python
cls = Langevin

md = cls(**settings.md.kwargs)
```

The keywords are documented in ASE:

- [`VelocityVerlet`](https://wiki.fysik.dtu.dk/ase/ase/md.html#ase.md.verlet.VelocityVerlet)
- [`Langevin`](https://wiki.fysik.dtu.dk/ase/ase/md.html#ase.md.langevin.Langevin)
- [`NPTBerendsen`](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.nptberendsen)

For `NPTBerendsen`, we add the `inhomogeneous` keyword, which decides if the `Inhomogeneous_NPTBerendsen` driver is used instead of `NPTBerendsen`. With `Inhomogeneous_NPTBerendsen`, the basis vectors are scaled independently, i.e. the size of the unit cell can change in three directions, but the angles remain constant. By default, `inhomogeneous=False`.

!!! warning
	In ASE, the `temperature` is usually given as an energy in ${\rm eV}$. In FHI-vibes, we use Kelvin consistently. `temperature: 300` thus corresponds to setting the thermostat to $300\,{\rm K}$.!!! info
	An hands-on example for setting up and running a relaxation can be found in the [Tutorial](../Tutorial/1_geometry_optimization.md).

vibes supports geometry optimization by setting up a `relaxation.in` file. A minimal `relaxation.in` would look like

```fo
[files]
geometry:                      geometry.in

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4

[relaxation]
driver:                        BFGS
fmax:                          0.001
```

for performing a BFGS optimization of the structure found in `geometry.in` until forces are converged below $\require{mediawiki-texvc} 1\,\text{meV}/\AA$.

??? info "Click: Default values for the complete list of supported keywords"
    ```
    [relaxation]
    driver:                        BFGS
    fmax:                          0.001
    unit_cell:                     True
    fix_symmetry:                  False
    hydrostatic_strain:            False
    constant_volume:               False
    scalar_pressure:               0.0
    decimals:                      12
    symprec:                       1e-05
    workdir:                       relaxation

    [relaxation.kwargs]
    maxstep:                       0.2
    logfile:                       relaxation.log
    restart:                       bfgs.restart
    ```

## Sections

### `[relaxation]`
Instructions to set up a geometry optimization workflow using an [ASE optimizer class](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#module-ase.optimize).

#### `driver`
Currently only [BFGS](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#bfgs) is supported, which is Quasi-Newton method using the [Broyden–Fletcher–Goldfarb–Shanno algoritm](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) to obtain an estimation of the Hessian.

#### `fmax`

`float`: Maximum residual force in $\require{mediawiki-texvc} \text{eV}/\AA$ (for the stress components: $\text{eV}/\AA^3$ ).

#### `unit_cell`

`True/False`: relax the unit cell using [`ase.constraints.ExpCellFilter`](https://wiki.fysik.dtu.dk/ase/ase/constraints.html?highlight=expcellfilter#ase.constraints.ExpCellFilter)

#### `fix_symmetry`

`True/False`: keep the spacegroup of the system fixed using [`ase.constraints.FixSymmetry`](https://wiki.fysik.dtu.dk/ase/dev/ase/constraints.html?highlight=fixsymmetry#ase.spacegroup.symmetrize.FixSymmetry)

#### `hydrostatic_strain`

`True/False`: apply isotropic pressure instead of stress for cell deformation, see [here](https://wiki.fysik.dtu.dk/ase/ase/constraints.html?highlight=expcellfilter#ase.constraints.ExpCellFilter)

#### `constant_volume`

`True/False`: keep volume constant, see [here](https://wiki.fysik.dtu.dk/ase/ase/constraints.html?highlight=expcellfilter#ase.constraints.ExpCellFilter)

#### `scalar_pressure`

`float`: apply external pressure given in $\text{eV}/\AA^3$

#### `mask`

`list`: `(6, 1)` shaped mask to enable/disable relaxation of strain components in Voigt notation, e.g., to allow for shape-preserving relation. See [ASE documentation of `mask` keyword in `ExpCellFilter` class for details](https://wiki.fysik.dtu.dk/ase/ase/constraints.html?highlight=expcellfilter#ase.constraints.ExpCellFilter). **Using the `mask` keyword will alter the effective stress used to relax the lattice and changes the behavior of `fmax`. Check!**.

#### `decimals`

`int`: number of digits used to round positions before writing `geometry.in.next_step`

#### `symprec`

`float`: symmetry precision used for detecting space group when `fix_symmetry=True`.

#### `workdir`

The working directory for running the calculations and storing output files.

### `[relaxation.kwargs]`

These keywords are used 1:1 by the ASE optimizer class, e.g.

```pytho
cls = ase.optimize.BFGS

optimzer = cls(**settings.relaxation.get("kwargs"))
```

#### `maxstep`

`float`: largest allowed move

#### `logfile`

`str`: logfile for the relaxation

#### `restart`

`str`: use this file to store restart information# Molecular Dynamics

!!! info
	We assume that you are familiar with the basics of molecular dynamics simulations and you are here to learn how to perform them with `FHI-vibes`.

## Recap

Molecular Dynamics (MD) simulations aim at exploring the dynamical properties of a system defined by the Hamiltonian $\mathcal H$

$$
\begin{align}
\mathcal{H}(\mathbf{R}, \mathbf{P})=\sum_{I} \frac{\mathbf{P}_{I}^{2}}{2 M_{I}}+\mathcal{V}(\mathbf{R})~,
\label{eq:H}
\end{align}
$$

where ${\bf R} = \{ {\bf R}_I \}$ denotes the atomic positions, ${\bf P} = \{ {\bf P}_I \}$ the atomic momenta, $M_I$ the atomic mass, and $\mathcal V ({\bf R})$ is a many-body potential, e.g., given by an empirical force-field, or an _ab initio_ energy functional. The dynamical evolution of the  system is described by the equations of motion,

$$
\begin{align}
M_{I} \ddot{\mathbf{R}}_{I}(t)=\mathbf{F}_{I}(t)=-\nabla_{I} \mathcal{V}(\mathbf{R}(t))~,
\label{eq:Newton}
\end{align}
$$

where the acceleration $\ddot{\mathbf{R}}_{I}(t)$ of an atom at time $t$ is given by the atomic force ${\bf F}_I (t)$, i.e., the inverse gradient of the many-body potential $\mathcal V ({\bf R})$.

Equation $\eqref{eq:Newton}$ is solved numerically from a given initial condition $\{{\bf R} (t_0), {\bf P} (t_0) \}$, for example with the [Velocity Verlet integrator](https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet).

The dynamical information encoded in the phase-space trajectory $\Gamma (t) = \{{\bf R} (t), {\bf P} (t) \}$ can then be used to evaluate expectation values of observables, which are typically functions of phase-space points:

$$
\begin{align}
	\langle O\rangle
	%&= \frac{1}{\mathcal{Z}}
	%\int \mathrm{d} \mathbf{R} \mathrm{d} \mathbf{P} ~
	%	\mathrm{e}^{-\beta \mathcal{H}(\mathbf{R}, {\bf P})}
	%	O(\mathbf{R}, {\bf P})
	%	\label{eq:O1} \\
	&=
	\lim_{T \rightarrow \infty} \frac{1}{T}
	\int_0^T {\rm d} t ~
	O(\mathbf{R} (t), {\bf P} (t))
	\label{eq:O2}~,
	\end{align}
$$

where Eq. $\eqref{eq:O2}$ holds for [ergodic systems](https://en.wikipedia.org/wiki/Ergodicity).

Compared to approximate treatments of the nuclear dynamics, e.g., the [harmonic approximation discussed earlier](2_phonopy_intro.md), MD has the advantage that it accounts for the full potential~$\mathcal{V}(\mathbf{R})$ and
thus also for _anharmonic effects_ (contributions not captured by a the harmonic approximation).

### Example: Pressure

For example, $O$ could be the instantaneous pressure given by

$$
\begin{align}
\def\d\{{\rm d}}
p({\bf R}, {\bf P})
	= - \left. \left( \frac{\d \mathcal H ({\bf R}, {\bf P})}{\d V} \right)\right|_T
	= \frac{1}{3V} \sum_I \frac{{\bf P}_I^2}{M_I}
	+ \frac{1}{V} \frac{\d \mathcal V ({\bf R})}{\d V}
	~,
	\label{eq:p}
\end{align}
$$

where $\d / \d V$ denotes a volume derivative. For evaluating the thermodynamic expectation value of Eq. $\eqref{eq:p}$, we first note that $p$ decouples into two contributions

$$
\begin{align}
	p({\bf R}, {\bf P}) = p_{\rm Kin} ({\bf P}) + p_{\rm Pot} ({\bf R})~
	\label{eq:p2}
\end{align}
$$

which can be evaluated independently of each other. We thus have

$$
\begin{align}
	\left\langle p \right\rangle
		= \left\langle p_{\rm Kin} \right\rangle
		+ \left\langle p_{\rm Pot} \right\rangle~,
	\label{eq:p3}
\end{align}
$$

where $p_{\rm Kin}$ denotes the kinetic ideal gas contribution, which yields the familiar [[HansenMcDonald]](references.md#HansenMcDonald)

$$
\left\langle p_{\rm Kin} \right\rangle = \frac{N k_{\rm B} T}{V}~,
$$

and $p_{\rm Pot}$ denotes the potential contribution which needs to be evaluated for the configurations ${\bf R}(t)$ generated by Eq. $\eqref{eq:Newton}$, so that

$$
\begin{align}
\left\langle p_{\rm Pot} \right\rangle
	= \lim_{N_{\rm t} \rightarrow \infty} \frac{1}{N_{\rm t}}
	\sum_n^{N_{\rm t}}
	\frac{1}{V} \frac{\d \mathcal V \left({\bf R} (t_n) \right)}{\d V}
\label{eq:<pPot>}
\end{align}
$$

is the expectation value of pressure in a simulation with discrete time steps $t_n$. Since consecutive time steps in an MD simulation are necessarily correlated, Eq. $\eqref{eq:<pPot>}$ can converge quite slowly with the number of time steps $N_{\rm t}$. We come back to this in the next tutorials.

## Canonical Ensemble: Thermostats
A simulation as described above models a [_microcanonical ensemble_](https://en.wikipedia.org/wiki/Microcanonical_ensemble), where particle number $N$, volume $V$, and energy $E$ are conserved.

In order to simulate a [_canonical ensemble_](https://en.wikipedia.org/wiki/Canonical_ensemble), where instead of the energy $E$ the temperature $T$ is the thermodynamic variable, one typically models the system including an interaction with a fictitious heat bath which allows to control the target temperature of the  system.

### Langevin thermostat
`FHI-vibes` uses a Langevin thermostat for canonical sampling. In Langevin dynamics, modified equations of motion are used with

$$
\begin{align}
	\dot{\bf P}_I(t)
		= {\bf F}_I
		- \gamma {\bf P}_I(t)
		+ \sqrt{2 M_I \gamma T} \xi(t)
	\label{eq:Langevin}~,
\end{align}
$$

where $\gamma$ is a friction parameter and $\xi (t)$ is a white-noise term obeying $\langle\xi(t) \xi(0)\rangle=2 k_{\rm B} T \gamma \delta(t)$.
# References

1. <a name="Knoop2020"></a> [F. Knoop, T.A.R. Purcell, M. Scheffler, and C. Carbogno, (2020), arXiv:2006.14672](https://arxiv.org/abs/2006.14672)
2. <a name="Baroni2001"></a> [S. Baroni, S. de Gironcoli, and A. Dal Corso, Rev. Mod. Phys. **73**, 515 (2001).](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.73.515)
3. <a name="AshcroftMermin"></a> N. W. Ashcroft and N. D. Mermin, _Solid State Physics_, Saunders College Publishing (1976).
4. <a name="BornHuang"></a>  M. Born, and K. Huang, _Dynamical Theory of Crystal Lattices_, Oxford University Press (1962).
5. <a name="Togo2015"></a> [Atsushi Togo and Isao Tanaka, Scr. Mater., **108**, 1-5 (2015)](https://phonopy.github.io/phonopy/)
6. <a name="Parlinski1997"></a> [K. Parlinski, Z. Q. Li, and Y. Kawazoe, Phys. Rev. Lett. **78**, 4063 (1997).](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.78.4063)
7. <a name="West2006"></a> [D. West, and S. K. Estreicher, Phys. Rev. Lett. **96**, 115504 (2006)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.115504)
8. <a name="Dove1993"></a> M. Dove, _Introduction to Lattice Dynamics_, Cambridge University Press (1993)
9. <a name="HansenMcDonald"></a> J.P. Hansen, and I.R. McDonald, _Theory of simple liquids_, Elsevier, 1990.
10. <a name="Nocedal"></a> J. Nocedal, and S. Wright, _Numerical optimization_, Springer Science & Business Media, 2006.
# Postprocessing

In this tutorial, you learn how to access the data you obtained from the MD calculation. It further showcases the interplay of `xarray`, `pandas`, `numpy` , and `matplotlib` for convenient exploration of data. We begin by showing how to convert the `xarray.Dataset` to plain `numpy` arrays if you prefer to work that way.

!!! warning
	We assume you successfully ran the simulation of LDA-Silicon at $300\,{\rm K}$ from the [previous chapter](3_md_ab_initio.md) and created a `trajectory.nc` dataset. You can can also take the reference data from [here](https://gitlab.com/vibes-developers/vibes-tutorial-files/-/tree/master/3_molecular_dynamics/ab_initio/si_8).

## Load and Inspect Dataset

The data in `trajectory.nc` is stored in an [`xarray.Dataset`](http://xarray.pydata.org/en/stable/data-structures.html#dataset), which is a collection of labeled `numpy` arrays. Let's open the file and inspect the content:

```python
import xarray as xr

# load the trajectory dataset
dataset = xr.load_dataset("trajectory.nc")

print(dataset)
```

??? info "`print(dataset)`"

    ```
    <xarray.Dataset>
    Dimensions:                    (I: 8, Ia: 24, Jb: 24, a: 3, b: 3, time: 2501)
    Coordinates:
      * time                       (time) float64 0.0 4.0 8.0 ... 9.996e+03 1e+04
    Dimensions without coordinates: I, Ia, Jb, a, b
    Data variables:
        positions                  (time, I, a) float64 0.0 0.0 0.0 ... 4.073 1.298
        displacements              (time, I, a) float64 -2.98e-11 ... -0.05708
        velocities                 (time, I, a) float64 0.03397 0.01561 ... 0.04307
        momenta                    (time, I, a) float64 0.9542 0.4384 ... 1.21
        forces                     (time, I, a) float64 0.08174 0.08174 ... 0.237
        energy_kinetic             (time) float64 0.3102 0.2771 ... 0.3916 0.3169
        energy_potential           (time) float64 -6.299e+04 ... -6.299e+04
        stress                     (time, a, b) float64 -0.001528 ... -0.005702
        stress_kinetic             (time, a, b) float64 -0.001839 ... -0.002492
        stress_potential           (time, a, b) float64 0.0003113 ... -0.003209
        temperature                (time) float64 300.0 268.0 215.7 ... 378.7 306.4
        cell                       (time, a, b) float64 5.419 2.4e-11 ... 5.419
        positions_reference        (I, a) float64 2.98e-11 1.084e-11 ... 4.064 1.355
        lattice_reference          (a, b) float64 5.419 -2.4e-11 ... 4e-12 5.419
        force_constants_remapped   (Ia, Jb) float64 13.7 0.0 0.0 ... 9.099e-11 13.7
        forces_harmonic            (time, I, a) float64 0.0818 0.08181 ... 0.217
        energy_potential_harmonic  (time) float64 0.002015 0.02751 ... 0.3128 0.3856
        sigma_per_sample           (time) float64 0.006452 0.03701 ... 0.1663 0.1836
        pressure                   (time) float64 0.0009887 0.001096 ... nan 0.00426
        pressure_kinetic           (time) float64 0.0013 0.001161 ... 0.001328
        pressure_potential         (time) float64 -0.0003113 -6.477e-05 ... 0.002932
        aims_uuid                  (time) object '57870CD4-D00E-43EB-A6F6-D2C87EEDBD0C' ... 'FF5DB1AC-F1C5-433F-AE30-913BF062C6B6'
    Attributes:
        name:             trajectory
        system_name:      Si
        natoms:           8
        time_unit:        fs
        timestep:         4.000000000000006
        nsteps:           2500
        symbols:          ['Si', 'Si', 'Si', 'Si', 'Si', 'Si', 'Si', 'Si']
        masses:           [28.085 28.085 28.085 28.085 28.085 28.085 28.085 28.085]
        atoms_reference:  {"pbc": [true, true, true],\n"cell": \n[[ 5.41850551468...
        atoms_primitive:  {"pbc": [true, true, true],\n"cell": \n[[-1.00000000000...
        atoms_supercell:  {"pbc": [true, true, true],\n"cell": \n[[ 5.41850551468...
        volume:           159.08841208433154
        raw_metadata:     {"MD": {\n  "type": "molecular-dynamics",\n  "md-type":...
        hash:             2d33a63b08cd6018441fa85ece2ba97d357eefc8
        sigma:            0.1561093848741265
    ```

Each of the `Data variables` can be accessed as an attribute and converted to a plain `numpy` array by calling `.data`. To access e.g. the positions, we do

```python
positions_as_ndarray = dataset.positions.data
```

The arrays can be extracted and written to file e.g. via [`numpy.savetxt`](https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html) or by other means.

## Example: Analyze Pressure

As an example on how to perform postprocess directly from the `xarray` dataset, we will now evaluate the potential pressure observed during the simulation, as [introduced earlier](3_md_intro.md#example-pressure). We will use [xarray](http://xarray.pydata.org/) and [pandas](https://pandas.pydata.org/) for the analysis. For interactive data exploration, we recommend to run the code in a [jupyter notebook](https://jupyter.org/) and play around with the suggested parameters like windows sizes etc.

### Load trajectory dataset

We first load the trajectory dataset and visualize the temperature:

```python
import xarray as xr

# load the trajectory dataset "trajectory.nc" from disk into the xarray.Dataset
dataset  = xr.load_dataset("trajectory.nc")

# extract temperature and potential pressure from all the data and convert to pandas.DataFrame
df_temperature_pressure = all_data[["temperature", "pressure_potential"]].to_dataframe()

# attach a moving average (width=200 timesteps) of the temperature
df_temperature_pressure["temperature_mean"] = df_temperature_pressure.temperature.rolling(window=200).mean()

# plot temperature and temperature_mean as function of time
ax = df_temperature_pressure[["temperature", "temperature_mean"]].plot()

ax.set_xlabel("Time (fs)")
ax.set_ylabel("Temperature (K)")
```

??? info "`df_temperature_pressure.plot`"
	![image](assets/md_temperature.png)

Since the calculation starts with all atoms located at their equilibrium positions, the initial potential energy is zero and the kinetic energy corresponds to ~300K, since we have setup the velocities
using the Maxwell-Boltzmann distribution. In the first 500 steps of the trajectory, the kinetic energy is partially converted to potential energy at. In turn, the temperature drops from $300\,{\rm K}$ to about $150\,{\rm K}$.
The missing thermal energy to obtain a temperature of 300K is then gradually provided by the thermostat, bringing the nuclear temperature back to $\sim 300\,{\rm K}$ after a few $\rm ps$.

### Discard thermalization period
We can remove the thermalization period from the simulation data, e.g., by [shifting the dataframe](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.shift.html):

```python
# discard 500 steps (2ps) of thermalization
shift = 500

df_temperature_pressure = df_temperature_pressure.shift(-shift).dropna()

df_temperature_pressure[['temperature', 'temperature_mean']].plot()
```

??? info "`df_temperature_pressure.plot` after removing 2500 simulation steps ($5\,{\rm ps}$)"
	![image](assets/md_temperature_thermalized.png)

### Inspect the pressure
We are now ready to inspect the pressure observed in the simulation and plot it including its cumulative average:

```python
from ase.units import GPa

p = df_temperature_pressure.pressure_potential / GPa

ax = p.plot(alpha=0.75)

p.expanding().mean().plot(ax=ax, color="k")
```

??? info "Plot pressure"
	![image](assets/md_pressure.png)


### Expectation value and convergence estimation

[As discussed earlier](3_md_intro.md), the expectation value of the pressure is given by the mean of the observed pressures,

$$
\begin{align}
\left\langle p_{\rm Pot} \right\rangle
	= \lim_{N_{\rm t} \rightarrow \infty} \frac{1}{N_{\rm t}}
	\sum_n^{N_{\rm t}}
	p_{\rm Pot}({\bf R} (t_n))~.
\label{eq:<pPot>}
\end{align}
$$

In our finite simulation, $N_{\rm t} = 2000 < \infty$, so that

$$
\begin{align}
\left\langle p_{\rm Pot} \right\rangle
= \left\langle p_{\rm Pot} \right\rangle_{N_t = 2000} + \Delta~,
\label{eq:p_final}
\end{align}
$$

where $\left\langle p_{\rm Pot} \right\rangle_{N_t = 2000} = -0.076\,{\rm GPa}$ is the mean pressure observed during the finite simulation, and $\Delta$ is the (unknown) difference to the fully converged expectation value.
Although, full converge would require an infinite trajetcory length and is thus formally never reachable, one can get arbitrarily close in practice and estimate the magnitude of the error $\Delta$.

We estimate this error by computing $\sigma_{\langle p \rangle}$, the [_standard error of the mean_](https://en.wikipedia.org/wiki/Standard_error):

$$
\begin{align}
\Delta \approx \sigma_{\langle p \rangle} = \frac{\sigma_p}{\sqrt{\tilde N_t}}~,
\label{eq:sigma_O}
\end{align}
$$

where $\sigma_p$ is the standard deviation of the pressure distribution observed during the simulation, and $\tilde N_t$ is an estimate of the number of _uncorrelated_ samples provided by the simulation. To this end, we estimate

$$
\begin{align}
\tilde N_t = N_t / \tau~,
\label{eq:N}
\end{align}
$$

where $\tau$ is the correlation time for the pressure.
The most straightforward way to compute $\tau$ is to evaluate the [autocorrelation function](https://en.wikipedia.org/wiki/Autocorrelation) and estimate its decay time:

```python
# estimate correlation time
import pandas as pd
from scipy import signal as si

# substract the mean pressure
pp = p - p.mean()

# get the autocorrelation function from
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.correlate.html
corr = si.correlate(pp, pp)[len(pp) - 1 :]

# normalize to C(0) = 1
corr /= corr[0]

# create as pandas.Series for plotting
s = pd.Series(corr).rolling(min_periods=0, window=10).mean()
ax = s.plot()

# estimate correlation time from the drop below 0.1
tau = s.index.where(s < 0.1).min()
ax.axvline(tau)

ax.set_xlim(0, 100)
ax.set_title(f"$\\tau$ is {int(tau)} steps")
```

??? info "Plot pressure autocorrelation function"
	![image](assets/md_autocorr.png)

In the  present example, the observable decorrelates after about 10 time steps ($\equiv 40\,{\rm fs}$). We therefore estimate the number of uncorrelated samples to be

$$
\begin{align*}
	\tilde N_t = N_t / 10 \approx 200
\end{align*}
$$

The standard deviation of the pressure distribution is

$$
\begin{align}
	\sigma_p = 0.239\,{\rm GPa}~,
\end{align}
$$

so that according to Eq. $\eqref{eq:sigma_O}$,

$$
\sigma_{\langle p \rangle} = \frac{0.239}{\sqrt{200}}\,{\rm GPa} \approx 0.053\,{\rm GPa}~.
$$


The final result for the pressure according to Eq. $\eqref{eq:p_final}$ is

$$
\begin{align*}
	\langle p_{\rm Pot} (300\,{\rm K}) \rangle = (-0.076 \pm 0.053)\,{\rm GPa}~,
\end{align*}
$$

which means that our result is converged within an estimated precision of $70\,\%$. **Remark:** This does _not_ mean that the true expectation lies within the given range. The estimated error is to be understood in the sense of a [confidence interval](https://en.wikipedia.org/wiki/Confidence_interval#Practical_example). The size of the error signals that the calculation is not fully converged and more sampling time would be necessary to report the observed pressure with confidence. You find reference for a total simulation time of $40\,{\rm ps}$ [here](https://gitlab.com/vibes-developers/vibes-tutorial-files/-/tree/master/3_molecular_dynamics/ab_initio/si_8_longer). How did the value and the error change?

Physics question: The observed potential pressure is negative. _Why?_ Do you expect a positive or negative lattice expansion at $300\,{\rm K}$?

??? info "Code snippet to compute the mean and the error estimator"

    ```python
    mean = p.mean()
    std = p.std()
    err = std / (len(p) / tau) ** 0.5

    print(f"Mean:  {mean:.5f} GPa")
    print(f"Std.:  {std:.5f} GPa")
    print(f"Error: {err:.5f} GPa ({abs(err / mean) * 100:.2f} %)")
    ```

### More examples

For more examples on how to directly work with the trajectory dataset in `trajectory.nc`, please have a look at  the [ASE Workshop Tutorial](https://gitlab.com/flokno/ase_workshop_tutorial_19) which analyzes _ab initio_ MD data for a perovskite.
# Tutorial

!!! warning "Warnings"

	- All tutorials assume you have a background in (_ab initio_) calculations and in the fundamental theory of vibrations in solids. The background theory is written down to establish a common notation;
	    and introduction to these topics can be found various textbooks, some of which you find in the [references](references.md).
	- The settings used throughout the tutorials are chosen in order to allow for smooth calculations. They are _not_ sufficient for producing publication-ready scientific results.
	- We assume you have [installed](../Installation.md) and [configured](../Installation.md#configuration) `FHI-vibes` successfully.
	- We assume you have some working experience with `python`, `numpy`, and `jupyter-notebook`.

In this tutorial, we introduce the functionality of `FHI-vibes` with hands-on examples.

### Outline

The following tutorials are available:

- [Single point calculations](0_singlepoint.md)
- [Geometry optimization](1_geometry_optimization.md)
- [Phonon calculations](2_phonopy_intro.md)
- [Molecular dynamics](3_md_intro.md)
- [Harmonic sampling](4_statistical_sampling.md)
- [Anharmonicity quantification](5_anharmonicity_quantification.md)
- [High-Throughput workflows](../High_Throughput/Tutorial/0_configuring_fw_for_vibes.md)

All tutorials discuss fcc-Silicon, which is treated at the _ab initio_ level using FHI-aims and the LDA exchange-correlation functional.
Although only very moderate computational ressources are needed (mostly few minutes of runtime on a modern multi-core node), the tutorials
are intended to be run on a computing cluster, since they aim at showcasing typical
FHI-vibes usage.

!!! info
	We assume that you are familiar with running *FHI-aims* for performing _ab initio_ calculations.



Typically, FHI-vibes requires two files: One describes the geometry of the system
following  syntax used for `geometry.in` files in FHI-aims. For fcc-Silicon in the
primitive unit cell, it thus reads

??? info "Geometry input file `geometry.in`"
    ```
    lattice_vector 0.0000000000000000 2.7149999999999999 2.7149999999999999
    lattice_vector 2.7149999999999999 0.0000000000000000 2.7149999999999999
    lattice_vector 2.7149999999999999 2.7149999999999999 0.0000000000000000
    atom_frac 0.0000000000000000 0.0000000000000000 0.0000000000000000 Si
    atom_frac 0.2500000000000000 0.2500000000000000 0.2500000000000000 Si
    ```

The second files describes the computational tasks and the numerical parameters used
in the calculation. Accordingly, it also contains a `calculator` section that specifies
that FHI-aims shall be used at the LDA level of theory and which numerical settings shall
be used for the Silicon example in the tutorial. This section thus reads:

??? info "`calculator` section for task in put file"
    ```
	[files]
	geometry:                      geometry.in

    [calculator]
    name:                          aims
    socketio:                      true

    [calculator.parameters]
    xc:                            pw-lda
    compute_forces:                true

    [calculator.kpoints]
    density:                       2

    [calculator.basissets]
    Si:                            light
    ```

Let's walk through the settings once:

- `[files]`
	- `geometry: geometry.in`: read the input geometry from `geometry.in`.
- `[calculator]`
	- `name: aims` means that `FHI-aims` will be used as explained [here](../Documentation/calculator_setup.md#calculator).
	- `socketio: true` means that the socketIO communication will be used. This will speed up the computation of related structures.
- `[calculator.parameters]`: these are settings that go directly to `control.in`
	- `xc: pw-lda` means that the pw-LDA exchange-correlation functional will be used.
	- `compute_forces: true` means that forces will be computed.
- `[calculator.kpoints]`: this is an optional way of setting k-point grids based on a target density (specifying `k_grid: X Y Z` in `[calculator.parameters]` is also possible!)
	- `density: 2` use a k-point density of at least 2 per $\require{mediawiki-texvc} \AA^{-3}$.
- `[calculator.basissets]`: Details on which basissets to use
	- `Si: light`: use _light default_ basis sets for silicon.

More details for each keyword can be found in the documentation.


!!! info
    For experimenting, testing, and rapid code-developing, it is often useful to use a force-field instead of an
    _ab initio_ calculator. For instance, one can use the Lennard-Jones (LJ) potential available within ASE to run all
    tutorials for LJ-Argon. Such calculations take only seconds, even on older laptops. The required geometry description
    and calculator settings are given below.


??? info "Geometry in put file `geometry.in`"
    ```
    lattice_vector 0.0000000000000000 2.6299999999999999 2.6299999999999999
    lattice_vector 2.6299999999999999 0.0000000000000000 2.6299999999999999
    lattice_vector 2.6299999999999999 2.6299999999999999 0.0000000000000000
    atom 0.0000000000000000 0.0000000000000000 0.0000000000000000 Ar
    ```
??? info "`calculator` section for task input file"
    ```
    	[files]
		geometry:                      geometry.in

        [calculator]
        name:                          lj

        [calculator.parameters]
        # parameters for LJ Argon
        sigma:    3.405
        epsilon:  0.010325
        rc:       8.0
    ```
<a name="1_GeometryOptimization"></a>

In this tutorial, you will learn how to perform a geometry optimization with `FHI-vibes`.

!!! info
	We give explicit references for LDA-Silicon. When using LJ-Argon, the only difference lies the definition of the calculator in the `[calculator]` section, and the respective structure defined in `geometry.in`.

## Define Inputs

For starting the relaxation, we use the `geometry.in` file for Silicon in the primitive unit cell  discussed in the [introduction](0_intro.md#test-systems) and copy it to the actual (empty) directory in which we are working. Generate a task input file for running a relaxation by copying the [calculator information for your test system](0_intro.md#test-systems) to a file called `relaxation.in`. Next, use the command line interface (CLI) of `FHI-vibes` to obtain default settings for performing the relaxation and appending them to the input file:

```
vibes template relaxation >> relaxation.in
```

In case of LDA-Silicon with `FHI-aims` calculator, the newly generated input file `relaxation.in` should look like this:

??? info "`relaxation.in`"
    ```
	[files]
	geometry:                      geometry.in

	[calculator]
    name:                          aims
    socketio:                      True

    [calculator.parameters]
    xc:                            pw-lda

    [calculator.kpoints]
    density:                       2

    [calculator.basissets]
    default:                       light

    [relaxation]
    driver:                        BFGS
    fmax:                          0.001
    unit_cell:                     True
    fix_symmetry:                  False
    hydrostatic_strain:            False
    constant_volume:               False
    scalar_pressure:               0.0
    decimals:                      12
    symprec:                       1e-05
    workdir:                       relaxation

    [relaxation.kwargs]
    maxstep:                       0.2
    logfile:                       relaxation.log
    restart:                       bfgs.restart
    ```

The settings file template you just generated contains all the necessary settings to set up and run a geometry optimization with `FHI-vibes` using `FHI-aims` as the force/stress calculator.
`FHI-vibes` will perform a [BFGS optimization of the structure as implemented in ASE](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#bfgs).
The keywords are explained in the [documentation](../Documentation/relaxation.md).

## Run calculation
You can start an interactive calculation with `vibes run relaxation` or by incorporating this command
in the respectives submission file (see Sec. Singlepoint).
We suggest pipe the output, e.g., like this:

```
vibes run relaxation | tee log.relaxation
```

`FHI-vibes` will create a working directory with the default name `relaxation` and will handle running the `FHI-aims` calculations necessary to perform the geometry optimization.
The log should read like that:

??? info "`relaxation.log`"
    ```
      [vibes.run]    run relaxation workflow with settings from relaxation.in

      [relaxation]   ** /draco/u/christia/Codes/vibes_v2/tutorials/GR/relaxation/trajectory.son does not exist, nothing to prepare
      [calculator]   Update aims k_grid with kpt density of 3 to [8, 8, 8]
      [calculator]   .. add `sc_accuracy_rho: 1e-06` to parameters (default)
      [calculator]   .. add `relativistic: atomic_zora scalar` to parameters (default)
      [calculator]   .. add `compensate_multipole_errors: False` to parameters (default)
      [calculator]   .. add `output_level: MD_light` to parameters (default)
      [calculator]   Add basisset `light` for atom `Si` to basissets folder.
      [calculator]   Calculator: aims
      [calculator]   settings:
      [calculator]     xc: pw-lda
      [calculator]     compute_forces: True
      [calculator]     k_grid: [8, 8, 8]
      [calculator]     sc_accuracy_rho: 1e-06
      [calculator]     relativistic: atomic_zora scalar
      [calculator]     compensate_multipole_errors: False
      [calculator]     output_level: MD_light
      [calculator]     compute_analytical_stress: True
      [calculator]     use_pimd_wrapper: ('localhost', 10011)
      [calculator]     aims_command: /u/christia/Codes/vibes_v2/run_aims.sh
      [calculator]     species_dir: /draco/u/christia/Codes/vibes_v2/tutorials/GR/relaxation/basissets
      [relaxation]   filter settings:
      [relaxation]     hydrostatic_strain: False
      [relaxation]     constant_volume: False
      [relaxation]     scalar_pressure: 0.0
      [relaxation]   driver: BFGS
      [relaxation]   settings:
      [relaxation]     type: optimization
      [relaxation]     optimizer: BFGS
      [relaxation]     maxstep: 0.2
      [socketio]     Use SocketIO with host localhost and port 10011
      [relaxation]   filter settings:
      [relaxation]     hydrostatic_strain: False
      [relaxation]     constant_volume: False
      [relaxation]     scalar_pressure: 0.0
      [relaxation]   Start step 0
      [relaxation]   Step 0 finished.
      [relaxation]   .. residual force:  0.000 meV/AA
      [relaxation]   .. residual stress: 289.641 meV/AA**3
      [vibes]        .. Space group:     Fd-3m (227)
      [relaxation]   clean atoms before logging
      [relaxation]   .. log
      [relaxation]   Step 1 finished.
      [relaxation]   .. residual force:  0.000 meV/AA
      [relaxation]   .. residual stress: 3.463 meV/AA**3
      [vibes]        .. Space group:     Fd-3m (227)
      [relaxation]   clean atoms before logging
      [relaxation]   .. log
      [relaxation]   Step 2 finished.
      [relaxation]   .. residual force:  0.000 meV/AA
      [relaxation]   .. residual stress: 0.039 meV/AA**3
      [vibes]        .. Space group:     Fd-3m (227)
      [relaxation]   clean atoms before logging
      [relaxation]   .. log
      [relaxation]   Relaxation converged.
      [relaxation]   done.
    ```

You will find the FHI-aims in- and output in `relaxation/calculation/`, the final converged structure in `relaxation/geometry.in.next_step`, and a summary of the relaxtion path in `relaxation/relaxation.log`.

For a detailed summary of the relaxation path, you may run

```
vibes info relaxation relaxation/trajectory.son
```

??? info "Output"
    ```
    Relaxation info for relaxation/trajectory.son:
    fmax:             1.000e+00 meV/AA
    # Step |   Free energy   |   F-F(1)   | max. force |  max. stress |  Volume  |  Spacegroup  |
    #      |       [eV]      |    [meV]   |  [meV/AA]  |  [meV/AA^3]  |  [AA^3]  |              |

        1    -15748.20070140     -0.605222       0.0000         0.3084     39.773   Fd-3m (227)
    --> converged.
    ```
# Harmonic Sampling

!!! warning "Warnings"

	- **Never ever rely on estimates obtained by the following method _without checking for the relevance of anharmonic effects_ (see [next tutorial](5_anharmonicity_quantification.md)).**
	- The tutorial assumes you are familiar with performing [phonon calculations](2_phonopy.md).

### Thermodynamic sampling via MD simulations

As outlined in [the MD tutorial](3_md_intro.md), MD simulations provide a way to sample phase space in order to obtain thermodynamic expectation values of observables,

$$
\begin{align}
	\langle O \rangle
	=
	\lim_{t_0 \rightarrow \infty} \frac{1}{t_0}
	\int_0^{t_0} {\rm d} t ~
	O \left( \mathbf{R} (t), {\bf P} (t)\right)
	\label{eq:O2}~,
\end{align}
$$

where  $O \left( \mathbf{R} (t), {\bf P} (t) \right)$ denotes the instantaneous value of the observable $O$ obtained for the phase-space configuration $\{{\bf R} (t),  {\bf P} (t) \}$ evaluated with respect to the the full many-body Hamiltonian $\mathcal H ({\bf P}, {\bf R})$. For example, $O$ could be the instantaneous pressure [as used in the MD tutorial](3_md_postprocess.md).

### Circumventing the dynamical propagation by generating samples from a harmonic distribution

In MD simulations, the many-body potential $\mathcal V$ needs to be evaluated _both_ for propagating the dynamical evolution of the system, and for evaluating the observable along the trajectory. On a conceptual level, only the evaluation of the observable is necessary to provide an estimate of the thermodynamic expectation value _if representative samples can be generated by other means_.

Following [[West2006]](references.md#West2006), we now present a way to create these samples based on the [harmonic approximation](2_phonopy_intro.md). In this approximation, the potential is given by

$$
\begin{align}
	\mathcal V^{(2)} ({\bf R}) 
	= \frac{1}{2} \sum_{I, J}
		\Delta {\bf R}_I \cdot \Phi_{IJ} \Delta {\bf R}_J~,
	\label{eq:V2}
\end{align}
$$

with the $3 \times 3$ force constant matrices $\Phi_{IJ}$ and atomic displacements $\Delta {\bf R}_I$. The Newton equations of motion therefore read

$$
\begin{align}
\left(
\begin{array}{c}
M_{1} \Delta \ddot{\mathbf{R}}_{1} \\
M_{2} \Delta \mathbf{R}_{2} \\
\vdots \\
M_{N} \Delta \ddot{\mathbf{R}}_{N}
\end{array}\right)=\left(\begin{array}{cccc}
\Phi_{11} & \Phi_{12} & \cdots & \Phi_{1 N} \\
\Phi_{21} & \Phi_{22} & \cdots & \Phi_{2 N} \\
\vdots & \vdots & \ddots & \vdots \\
\Phi_{N 1} & \Phi_{N 2} & \cdots & \Phi_{N N}
\end{array}\right)\left(\begin{array}{c}
\Delta \mathbf{R}_{1} \\
\Delta \mathbf{R}_{2} \\
\vdots \\
\Delta \mathbf{R}_{N}
\end{array}\right)~,
\end{align}
$$

and can be solved analytically, yielding

$$
\begin{align}
\Delta {\bf R}_I (t) 
	= \frac{1}{\sqrt{M_I}} \sum_s {\bf e}_{sI} A_s \sin (\omega_s t + \phi_s)~,
\label{eq:R2(t)}
\end{align}
$$

where ${\bf e}_{sI}$ denotes the eigenvector of the dynamical matrix $D_{IJ} = \Phi_{IJ} / \sqrt{M_I M_J}$ corresponding to atom $I$ and phonon mode $s$ and $\omega_s$ the respective frequency, cf. [the phonon tutorial](2_phonopy_intro.md#phonons-harmonic-vibrations-in-solids). The amplitudes $A_s$ and phases $\phi_s$ are fixed by the initial condition $\{{\bf R} (t_0), {\bf P} (t_0) \}$, where in thermal equilibrium [[Dove1993]](references.md#Dove1993)

$$
\begin{align}
\left\langle A_{s}\right\rangle
	=
	\sqrt{
	\frac{\hbar}{\omega_{s}}\left(n_{\mathrm{B}}(\omega_s, T)+\frac{1}{2}\right) 
	}
	~\stackrel{k_{\rm B} T \,\gg\, \hbar \omega_s}{\longrightarrow}~
	\frac{\sqrt{2 k_{\rm B} T}}{\omega_s}~.
	\label{eq:As}
\end{align}
$$

Noting that $\sin (\omega_s t + \phi_s)$ is essentially a random number in $[-1, 1]$, and mimicking thermal fluctuations, samples can be generated by

$$
\begin{align}
\Delta {\bf R}_{I} (n)
	=\frac{1}{\sqrt{M_{I}}} \sum_{s} \zeta_{s} (n) \left\langle A_{s}\right\rangle {\bf e}_{s I}~,
\label{eq:samples}
\end{align}
$$

where $\langle A_s \rangle$ is given by Eq. $\eqref{eq:As}$ and each $\zeta_s (n)$ is a normally distributed random number, where  $n$ labels the sample.

## Example: LDA-Si at 300K

### Obtain force constants

We will re-use the calculation from [the previous tutorial on phonon calculations](2_phonopy.md) in the conventional cell with 8 atoms. If you didn't run these calculations, you find the respective calculations in [our reference repository](https://gitlab.com/vibes-developers/vibes-tutorial-files/-/tree/master/2_phonopy/sc_8).

As start, [perform the postprocess](2_phonopy.md#basic-postprocessing) in the `phonopy` folder:

```
cd phonopy
vibes output phonopy -bs
```

Check the bandstructure in `output/bandstructure.pdf` for plausibility.

#### Remap the force constants to the supercell

The force constants are written to [`output/FORCE_CONSTANTS`](../Documentation/output_files.md#force_constants) which is a condensed representation in $(N_{\rm prim}, N, 3, 3)$ shape, where $N_{\rm prim}$ is the number of atoms in the primitive cell and $N$ is the number of atoms in the supercell. For creating samples, they need to be mapped to a full $3 N \times 3N$ shape. This can be done with the CLI tool `vibes utils fc remap`,

```
cd output
vibes utils fc remap
```

which will create [`FORCE_CONSTANTS_remapped`](../Documentation/output_files.md#force_constants_remapped) as a plain, `numpy` readable text file.

The current folder `…/phonopy/output` should now contain the following files:

```
>>> ls
bandstructure.pdf  band.yaml  FORCE_CONSTANTS  FORCE_CONSTANTS_remapped  geometry.in.primitive  geometry.in.supercell
```

We are now ready to create samples according to Eq. $\eqref{eq:samples}$.

### Create samples

Samples can be created with the CLI tool `utils create-samples`. We now create 10 samples according to Eq. $\eqref{eq:samples}$ at a temperature of $300\,{\rm K}$ by running

```
vibes utils create-samples geometry.in.supercell -fc FORCE_CONSTANTS_remapped -n 10 -T 300
```

The geometry files are written as `geometry.in.supercell.0300K.???` to the  current directory. Let's move them to a folder called `samples_20K`:

```
mkdir samples_300K
mv geometry.in.supercell.* samples_300K
```

Create a new working directory where you

- copy the  files `geometry.in.primitive` and `geometry.in.supercell`,
- copy your samples-folder `samples_300K`.

Move to that directory. _Remark: copying the `geometry.in.primitive` is optional. However, you should always know from where you started, it might save some hours of unnecesary work at some point :-)._

### Compute the samples

In order to compute pressure in all samples, prepare an `aims.in` like this:

```
[files]
geometries:                    samples_300K/geometry.in.*
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell

[calculator]
name:                          aims
socketio:                      true
workdir:                       aims_300K

[calculator.parameters]
xc:                            pw-lda
compute_analytical_stress:     true

[calculator.kpoints]
density:                       2

[calculator.basissets]
default:                       light
```

This will run a calculation for all structures found in `samples_300K/geometry.in.*` [as discussed earlier in the tutorial on singlepoint calculations](0_singlepoint.md), and write the results to a `trajectory.son` file in the directory `aims_300K`.

Run the calculation with

```
vibes run singlepoint aims.in | tee log.aims
```

This will take a few minutes depending on your machine.

### Compare sampling vs. MD

[Similar to MD](3_md_ab_initio.md#postprocess) you can create a trajectory dataset from `trajectory.son` by running `vibes output md`:

```
cd  aims_300K
vibes output md
```

Since this is not an MD, the time axis in the dataset will simply label the samples instead of corresponding to an actual simulation time.

The dataset in `trajectory.nc` can be inspected similar to [the MD case](3_md_postprocess.md). For example, run

```python
import xarray as xr
from ase.units import GPa


ds = xr.load_dataset("trajectory.nc")

p = ds.pressure_potential.to_series() / GPa

ax = p.plot(marker="x", lw=0)

p.expanding().mean().plot(ax=ax, color="k")

ax.set_xlabel('Sample number')
```

??? info "Pressure plot"
	![image](assets/mc_pressure.png)
	
The [mean pressure and standard error](3_md_postprocess.md#expectation-value-and-convergence-estimation)[^footnote1] are

$$
\begin{align*}
\langle p_{\rm Pot} (300\,{\rm K}) \rangle^{(2)} = (-0.11 \pm 0.05)\,{\rm GPa}~,
\end{align*}
$$

which is converged with a precision of $\sim 20\,\%$. **However**, the superscript $\langle \cdot \rangle^{(2)}$ reminds us that we used the harmonic approximation to create the samples. Indeed, the pressure found by harmonic sampling is about $40\,\%$ larger than the [MD reference computed in the previous tutorial](3_md_postprocess.md#expectation-value-and-convergence-estimation). Yet, the  two values coincide within their error margins. Further sampling on both sides would be needed to compute the actual difference between the two sampling techniques.

## Take Home Messages

The sampling trick can be used to rapidly estimate the expectation value of a static observable and typically converges 1-2 orders of magnitude faster than a MD simulation[^footnote2].

The sampling trick only works if the harmonic approximation is a good reference for the system of interest. This might or might not be the case and is difficult to tell _a priori_. We present a scheme to estimate the importance of anharmonic effects in the next chapter of the tutorial.



[^footnote1]: The samples are uncorrelated by construction, so $\tau = 1$ and $\tilde N_t = N_t$ in this case.

[^footnote2]: Remember that the correlation length was about 200 steps in the MD case while the samples created from the harmonic distribution are uncorrelated per construction.
# Singlepoint Calculations

!!! info
	We assume you are familiar with `FHI-aims` calculations and are here to learn how to perform them with `FHI-vibes`. The tutorial is however transferable to any other calculator supported by ASE.

In this tutorial, you will learn how to perform singlepoint calculations with `FHI-vibes`.

## Silicon in Equilibrium

As discussed in the introduction, this requires two files. The geometry of the system given by

??? info "geometry.in"

    ```
    lattice_vector 0.0000000000000000 2.72 2.72
    lattice_vector 2.72 0.0000000000000000 2.72
    lattice_vector 2.72 2.72 0.0000000000000000
    atom_frac 0. 0. 0. Si
    atom_frac .25 .25 .25 Si
    ```

and a file desribing the computational taks, which contains the
calculator settings discussed in the introduction. Additionally,
it contains a section `files` pointing to the geometry that shall
be calculated.

??? info "aims.in"

    ```
    [files]
    geometry:                      geometry.in
    
    [calculator]
    name:                          aims
    socketio:                      true
    
    [calculator.parameters]
    xc:                            pw-lda
    compute_forces:                true
    
    [calculator.kpoints]
    density:                       3
    
    [calculator.basissets]
    Si:                            light
    ```

You can run the calculation interactively via

```
vibes run singlepoint aims.in | tee log.aims
```

or by submitting it to a queue on a computing cluster.

??? info "Example `submit.sh` for `slurm` queue manager"
    ```
    #!/bin/bash -l

    #SBATCH -J md|vibes
    #SBATCH -o log/md.%j
    #SBATCH -e log/md.%j
    #SBATCH --mail-type=all
    #SBATCH --mail-user=your@mail.com
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=32
    #SBATCH --ntasks-per-core=1
    #SBATCH -t 24:0:00
    
    # Make sure that the correct Python environment is set up, e.g.
    module load miniconda/3/4.5.4
    source activate vibes
    
    vibes run singlepoint aims.in &> log.aims
    ```

The calculation should only take a few seconds and yield the following output:

??? info "`log.aims`"

    ```
    [vibes.run]    run singlepoint calculations with settings from aims.in
    
    [calculator]   Update aims k_grid with kpt density of 3 to [8, 8, 8]
    [calculator]   .. add `sc_accuracy_rho: 1e-06` to parameters (default)
    [calculator]   .. add `relativistic: atomic_zora scalar` to parameters (default)
    [calculator]   .. add `compensate_multipole_errors: False` to parameters (default)
    [calculator]   .. add `output_level: MD_light` to parameters (default)
    [calculator]   Add basisset `light` for atom `Si` to basissets folder.
    [calculator]   Calculator: aims
    [calculator]   settings:
    [calculator]     xc: pw-lda
    [calculator]     compute_forces: True
    [calculator]     k_grid: [8, 8, 8]
    [calculator]     sc_accuracy_rho: 1e-06
    [calculator]     relativistic: atomic_zora scalar
    [calculator]     compensate_multipole_errors: False
    [calculator]     output_level: MD_light
    [calculator]     use_pimd_wrapper: ('localhost', 10011)
    [calculator]     aims_command: /u/christia/Codes/vibes_v2/run_aims.sh
    [calculator]     species_dir: /draco/u/christia/Codes/vibes_v2/tutorials/SPC_2/aims/basissets
    [socketio]     Use SocketIO with host localhost and port 10011
    [backup]       /draco/u/christia/Codes/vibes_v2/tutorials/SPC_2/aims/calculations does not exists, nothing to back up.
    [vibes]        Compute structure 1 of 1: working
    [vibes]        Compute structure 1 of 1: finished.
    ```

The calculation will create a working directory called `aims`, the traditional input files for FHI-aims (control.in and geometry.in) as well as the output file
`aims.out` can be found in `aims/calculations/`. Additionally, vibes produces a trajectory file `aims/trajectory.son`, which contains all salient information for
postprocessing and is particularly useful for the adavanced tasks tackled in the next tutorials. The data from this file can be extracted and stored
in an `xarray.Dataset` in `trajectory.nc`, see [the documentation on output files](../Documentation/output_files.md). For this purpose, run

```
vibes output trajectory aims/trajectory.son
```

which yields the output
??? info "Output of `vibes output trajectory`"
    ```
    Extract Trajectory dataset from <Command trajectory>
    [trajectory]   Parse `aims/trajectory.son`
    [son] read file:  aims/trajectory.son
    [son] process:    |||||||||||||||||||||||||||||||||||||  2/2
    [trajectory]   .. time elapsed: 0.000s
    [trajectory]   .. create atoms
    [trajectory]   .. time elapsed: 0.256s
    [trajectory]   .. done in 0.257s
    [trajectory]   Get positions from trajectory
    * Message from file vibes/trajectory/trajectory.py, line 192, function times:
    --> time unit not found in trajectory metadata, use ase.units.fs

    ** Warning from file vibes/trajectory/trajectory.py, line 198, function times:
    --> no time steps found, return time as index
    
    /u/christia/.local/lib/python3.7/site-packages/numpy/core/fromnumeric.py:3373: RuntimeWarning: Mean of empty slice.
      out=out, **kwargs)
    /u/christia/.local/lib/python3.7/site-packages/numpy/core/_methods.py:170: RuntimeWarning: invalid value encountered in double_scalars
      ret = ret.dtype.type(ret / rcount)
    [trajectory]   .. time elapsed: 0.071s
    [trajectory]   Get velocities from trajectory
    [trajectory]   .. time elapsed: 0.001s
    ** Warning from file vibes/trajectory/trajectory.py, line 540, function set_displacements:
    --> SUPERCELL NOT SET, compute w.r.t to reference atoms
    
    [trajectory]   Compute displacements
    [trajectory]   .. time elapsed: 0.002s
    [trajectory]   Get pressure from trajectory
    [trajectory]   .. time elapsed: 0.001s
    Trajectory dataset written to trajectory.nc
    ```


## Multiple Singlepoint Calculations in One Run

`FHI-vibes` offers the possibility to run a set of related calculations in a single run, where "related calculations" means that the input geometries are allowed to differ in their positions and/or lattice.
The stoichometry, the number of atoms, as well as the computational settings must be the same for all systems. For instance, one can run singlepoint calculations for the following two geometries

??? info "geometry.in.000"

    ```
    lattice_vector 0.0000000000000000 2.72 2.72
    lattice_vector 2.72 0.0000000000000000 2.72
    lattice_vector 2.72 2.72 0.0000000000000000
    atom_frac 0.01 0. 0. Si
    atom_frac .25 .25 .25 Si
    ```

??? info "geometry.in.001"

    ```
    lattice_vector 0.0000000000000000 2.72 2.72
    lattice_vector 2.72 0.0000000000000000 2.72
    lattice_vector 2.72 2.72 0.0000000000000000
    atom_frac 0.02 0. 0. Si
    atom_frac .25 .25 .25 Si
    ```

which only differ by in the first fractional coordinate of the first atom, by using the
`geometries` tag in the file `aims.in`


??? info "aims.in"

    ```
    [files]
    geometries:                    geometry.in.???
    
    [calculator]
    name:                          aims
    socketio:                      true
    
    [calculator.parameters]
    xc:                            pw-lda
    compute_forces:                true
    
    [calculator.kpoints]
    density:                       3
    
    [calculator.basissets]
    Si:                            light
    ```

Note, that `FHI-vibes` supports wildcards to read input files and will sort the input files found by this wildcard alphabetically.


Again, you can run the calculation interactively via

```
vibes run singlepoint aims.in | tee log.aims
```

or by submitting to queue as above.

## Submit calculation on a cluster

To efficiently perform _ab initio_ calculations for systems larger than a few atoms, you will need a workstation or access to a computing cluster. To submit a `vibes` simulation to your cluster, follow these steps:

1. [Install `FHI-vibes` on your cluster](../../#installation),
2. set up a calculation as you have done earlier on your laptop,
3. submit the `vibes run` command to the queue.

??? info "Example `submit.sh` for `slurm` queue manager"

    ```
    #!/bin/bash -l
    
    #SBATCH -J md|vibes
    #SBATCH -o log/md.%j
    #SBATCH -e log/md.%j
    #SBATCH --mail-type=all
    #SBATCH --mail-user=your@mail.com
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=32
    #SBATCH --ntasks-per-core=1
    #SBATCH -t 24:0:00
    
    vibes run singlepoint aims.in
    ```

The log file will be written to a `log` folder.

For default tasks, `vibes` can write and submit the `submit.sh` file in a single command, see [documentation of slurm submission](../Documentation/input_files_slurm.md)

## Restart a calculation

If your calculation does not fit into a walltime or stops for another reason before the total number of simulation steps is reached, you can simply resubmit `vibes run`. It will restart the calculation from the last completed step. This also holds for `vibes run phonopy` and `vibes run md` as explained later.

### Automatic restarts

Optionally, `vibes` can restart the job by itself using a `[restart]` section in the input file, e.g., `aims.in`. To this end, add

```
[restart]
command = sbatch submit.sh
```

to your `aims.in`, where `sbatch submit.sh` is the command you use to submit the calculation to the queue. `vibes` will run this command shortly before the walltime is over to restart the job.# Acknowledgments

This tutorial is based on tutorials prepared for the annual "DFT and Beyond" Workshop series by Maja Lenz, Christian Carbogno, Martin Fuchs, Felix Hanke, Jörg Meyer, Karsten Rasim, Manuel Schöttler, Amrita Bhattacharya, Honghui Shang, and Johannes Hoja.
# Anharmonicity Quantification

!!! warning
	The tutorial assumes you are familiar with performing [phonon calculations](2_phonopy.md) and [molecular dynamics simulations](3_md_ab_initio.md).

## Background

As detailed [in our paper](https://arxiv.org/abs/2006.14672), we define the anharmonic contribution to the potential energy $\mathcal V ({\bf R})$ as

$$
\begin{align}
	\mathcal{V}^{\rm A}(\mathbf{R}) \equiv \mathcal{V}(\mathbf{R})-\mathcal{V}^{(2)}(\mathbf{R})~,
	\label{eq:VA}
\end{align}
$$

where $\mathcal{V}^{(2)}(\mathbf{R})$ at a given atomic configuration $\bf R$ is given by

$$
\begin{align}
	\mathcal{V}^{(2)}\left(\mathbf{R}=\mathbf{R}^{0}+\Delta \mathbf{R}\right)
	=\frac{1}{2} \sum_{I, J} \Phi_{\alpha \beta}^{I, J} \Delta R_{I}^{\alpha} \Delta R_{J}^{\beta}~,
\label{eq:V2}
\end{align}
$$

with the [harmonic force constants $\Phi^{IJ}$](2_phonopy_intro.md) obtained at the equilibrium configuration ${\bf R}^0$ as

$$
\begin{align}
	\Phi_{\alpha, \beta}^{I, J}
	=\left.\frac{\partial^{2} \mathcal{V}}{\partial R_{I}^{\alpha} \partial R_{J}^{\beta}}\right|_{\mathbf{R}^{0}}~.
	\label{eq:Phi}
\end{align}
$$

Likewise, we define the anharmonic contribution to the force components $F_{I, \alpha} ({\bf R})$ as

$$
\begin{align}
	F_{I, \alpha}^{\mathrm{A}}(\mathbf{R})
	&=
	F_{I, \alpha}(\mathbf{R})-F_{t, \alpha}^{(2)}(\mathbf{R})~,\text{ with} \label{eq:FA} \\
	F_{I, \alpha}^{(2)}
	&=
	-\sum_{J, \beta} \Phi_{\alpha, \beta}^{I, J} \Delta R_{J}^{\beta}
	\label{eq:F2}
\end{align}
$$

This is a depiction of Eq. $\eqref{eq:VA}$ and $\eqref{eq:FA}$ for a one-dimensional toy potential:

![image](assets/PES_sketch.png)

In order to estimate the strength of anharmonic effects in a material, we define the _anharmonicity measure_

$$
\begin{align}
\sigma^{\mathrm{A}}(T) \equiv \frac{\sigma\left[F^{\mathrm{A}}\right]_{T}}{\sigma[F]_{T}}=\sqrt{\frac{\sum_{I, \alpha}\left\langle\left(F_{I, \alpha}^{\mathrm{A}}\right)^{2}\right\rangle_{T}}{\sum_{I, \alpha}\left\langle\left(F_{I, \alpha}\right)^{2}\right\rangle_{T}}}~,
\label{eq:sigmaA}
\end{align}
$$

where $\langle \cdot \rangle_T$ denotes an [expectation value at a given temperature](3_md_postprocess.md#expectation-value-and-convergence-estimation),

$$
\begin{align}
	\left\langle O \right\rangle
	= \lim _{N_{\mathrm{t}} \rightarrow \infty}
	\frac{1}{N_{\mathrm{t}}} \sum_{n}^{N_{\mathrm{t}}} \left(t_{n}\right)~.
	\label{eq:meanO}
\end{align}
$$

$F_{I, \alpha} (t) \equiv F_{I, \alpha} [{\bf R} (t)]$ is the force component $\alpha$ on atom $I$ at time $t$, and $F^{\rm A}_{I, \alpha}$ is given by Eq. $\eqref{eq:FA}$. $\sigma^{\rm A} (T)$ therefore quantifies the _average strength of anharmonic force components $F_{I, \alpha}^{\rm A}$, normalized by the average strength of forces $F_{I, \alpha}$, observed at temperature $T$_.



## Evaluating anharmonicity with `FHI-vibes`

The necessary ingredient to evaluate Eq. $\eqref{eq:sigmaA}$ are:

- Atomic forces $F_{I, \alpha}$,
- harmonic force constants $\Phi^{IJ}$ to compute $F^{(2)}$ according to Eq. $\eqref{eq:F2}$ for evaluating $F^{\rm A}$ according to Eq. $\eqref{eq:FA}$ , and
- thermodynamic expectation values according to Eq. $\eqref{eq:meanO}$.

These ingredients can be obtained with `FHI-vibes` with the following workflow:

- Take the materials of interest and generate a reference structure, i.e., a primitive cell and a supercell.
- Obtain force constants for the supercell as introduced in the [phonons tutorial](2_phonopy.md).
- Run an MD simulation for the supercell as introduced in the [MD tutorial](3_md_ab_initio.md).

### Example: LDA-Silicon at room temperature

Assuming that  you performed the previous tutorials for LDA-Silicon in an 8-atoms supercell, we already have all the necessary ingredients available to evaluate $\sigma^{\rm A}$ for this system!

In a new working directory, copy your `trajectory.nc` dataset from the the [MD tutorial](3_md_ab_initio.md) and your force constants from the [phonopy tutorial](2_phonopy.md), i.e., the file `phonopy/output/FORCE_CONSTANTS`. You can attach the force constants to the trajectory dataset with the CLI tool `utils trajectory update`:

```
vibes utils trajectory update trajectory.nc -fc FORCE_CONSTANTS
```

This will attach read the force constants from `FORCE_CONSTANTS` and attach them to the trajectory dataset.

To evaluate Eq. $\eqref{eq:sigmaA}$, you can use the CLI tool `info anharmonicity`:

```
vibes info anharmonicity trajectory.nc
```

which will give you the total $\sigma^{\rm A}$ value (`sigma`), as well as an individual value for each atom species. The output should be

```
DataFrame:
       sigma  sigma [Si]  sigma_atom_mean  sigma_mode
Si  0.156109    0.156109         0.156109    0.156026
```

This tells you that the average magnitude of anharmonic contributions to the forces, $F^{\rm A}$, in LDA-Silicon at $300\,{\rm K}$ is about $16\,\%$.

## Mode resolved anharmonicity

To obtain a mode-resolved $\sigma^{\rm A}_s$ similar to the analysis of Fig. 8 in [our paper](https://arxiv.org/pdf/2006.14672.pdf), you can run

```
vibes info anharmonicity trajectory.nc --per_mode
```

which will produce a `.csv` file containing mode frequencies $\omega_s$ in THz and the respective mode-resolved anharmonicity $\sigma^{\rm A}_s$.

You can plot the file e.g. via

```python
import pandas as pd

s = pd.read_csv("sigmaA_mode_Si.csv", index_col=0)

ax = s.plot(marker=".", lw=0)

ax.set_xlim(0, 20)
ax.set_ylim(0, 0.5)

ax.set_xlabel("$\omega_s$ (THz)")
ax.set_ylabel(r"$\sigma^{\rm A}_s$")
```

??? info "Plot of $\sigma^{\rm A}_s$"
	![image](assets/sigma_mode_Si.png)

The plot won't look too impressive because we're using a small supercell and the anharmonicity in silicon is overall quite weak. But you should be good to go to investigate the anharmonicity of your material of choice by now -- Happy Computing 💪<a name="2_Phonopy"></a>

!!! info "Prerequisites"

	- For vibrational studies, it is crucial to use structures that are accurately  relaxed. Before starting with actual phonon calculations, make sure you are familiar with [geometry optimization](1_geometry_optimization.md).
	- Create a new working directory and copy over the `geometry.in.next_step` file you obtained from the previous geometry optimization as your new `geometry.in` file.
	- We assume you already have some working knowledge for working with [`phonopy`](https://phonopy.github.io/phonopy/index.html) and understand [the underlying method](https://phonopy.github.io/phonopy/formulation.html).



## Perform a phonon calculation

Setting up a `phonopy` calculation is similar to settings up a `relaxation` (or any other workflow supported by `FHI-vibes`). To ensure that our physical settings don't change, we will copy the `relaxation.in` obtained in the [previous part of the tutorial](1_geometry_optimization.md) to the new working directory and rename it to `phonopy.in`. Please delete the `relaxation` specific sections `[relaxation]` and `[relaxation.kwargs]` and add settings for a phonopy calculation by running

```
vibes template phonopy >> phonopy.in
```

??? info "`phonopy.in`"
	```
	[calculator]
    name:                          aims

    [calculator.parameters]
    xc:                            pw-lda

    [calculator.kpoints]
    density:                       2

    [calculator.basissets]
    default:                       light

    [calculator.socketio]
    port:                          12345

    [phonopy]
    supercell_matrix:              [1, 1, 1]
    displacement:                  0.01
    is_diagonal:                   False
    is_plusminus:                  auto
    symprec:                       1e-05
    q_mesh:                        [45, 45, 45]
    workdir:                       phonopy
    ```

Obviously the most important section in the `phonopy.in` input file is `[phonopy]` which containts information about how the supercells with displacements should be set up to compute the force constants from the [finite-differences method](0_intro.md#Phonons). An explanation for the full list of keywords is found in the [documentation](../Documentation/phonopy.md). The most important two are explaned in the following:

### Supercell Matrix (`supercell_matrix`)

The supercell matrix $M_{\rm S}$ given as `supercell_matrix` will be used to [generate the lattice of the supercell from the lattice of the primitive unitcell by matrix multiplication:](https://phonopy.github.io/phonopy/phonopy-module.html#supercell-matrix)

$$
\begin{align}
	\require{mediawiki-texvc}
	\def\t#1{\text{#1}}
	\begin{pmatrix}
		\mathbf a_\t{S}^\t{t} \\ \mathbf b_\t{S}^\t{t} \\ \mathbf c_\t{S}^\t{t}
	\end{pmatrix}
	=
	M_\t{S} \cdot
	\begin{pmatrix}
	\mathbf a_\t{u}^\t{t} \\ \mathbf b_\t{u}^\t{t} \\ \mathbf c_\t{u}^\t{t}
	\end{pmatrix}
	 ~.
	 \label{eq:smatrix}
\end{align}
$$

Here, $\mathbf a_\t{u}^\t{t}$, $\mathbf b_\t{u}^\t{t}$, $\mathbf c_\t{u}^\t{t}$ are the transposed lattice vectors (row-vectors) of the (primitive) unit cell and $\mathbf a_\t{S}$, $\mathbf b_\t{S}$, $\mathbf c_\t{S}$ label the lattice vectors of the supercell respectively. `supercell_matrix` can be given in any shape that lets itself transform trivially to a $3 \times 3$-matrix. For example, `[1, 1, 1]` gets transformed to the $3 \times 3$ unit matrix.

### Displacement (`displacement`)

The `displacement` tag will set the amplitude of the finite displacement in $\AA$. The same parameter is called [`DISPLACEMENT_DISTANCE` in `phonopy`](https://phonopy.github.io/phonopy/setting-tags.html#displacement-distance). In principle, this is a numerical parameter that needs to be optimized. A smaller `discplacement` results in a better approximation to the true second derivative of the potential. However, a too small displacement generates too small forces that can be severely affected by other sources of computational noise, e.g., finite grids etc. For production purposes, the default value of $d = 0.01\,\AA$ usually works quite well and with a properly set up force calculator, there is no need to increase the displacement further.

### Run the calculation

Let's stick to the default settings in `phonopy.in` for the moment and run the calculation with

```
vibes run phonopy | tee log.phonopy
```

The calculation should take only a few seconds (depending on you computer).

### Postprocessing

The `vibes run` command takes care that all _ab initio_ calculations are performed, but some additional,  postprocessing is needed to obtain the phonon-related quantities. The postprocessing  itself can be performed interactivley with

```
vibes output phonopy phonopy/trajectory.son --full
```

??? note "Terminal output"
	```
    [phonopy.postprocess] Start phonopy postprocess:
    [trajectory]   Parse `phonopy/trajectory.son`
    [son] read file:  phonopy/trajectory.son
    [son] process:    |||||||||||||||||||||||||||||||||||||  2/2
    [trajectory]   .. create atoms
    [progress]        |||||||||||||||||||||||||||||||||||||  1/1
    [trajectory]   .. done in 0.001s
    [phonopy.postprocess] .. done in 0.034s
    [phonopy.postprocess]
    Extract phonopy results:
    [phonopy.postprocess] .. q_mesh:   [45, 45, 45]
    [phonopy.postprocess] .. write force constants
    [phonopy.postprocess] Extract basic results:
    [phonopy.postprocess] .. write primitive cell
    [phonopy.postprocess] .. write supercell
    [phonopy.postprocess] .. write force constants to FORCE_CONSTANTS
    [phonopy.postprocess] Extract bandstructure
    [phonopy.postprocess] .. write yaml
    [phonopy.postprocess] .. plot
    [phonopy.postprocess] .. all files written to phonopy/output in 1.113s
    * Message from file vibes/phonopy/postprocess.py, line 123, function check_negative_frequencies:
        --> Negative frequencies found at G = [0 0 0]:

    # Mode   Frequency
        1 -1.62419e-07 THz
        2 -1.15231e-07 THz
    [phonopy.postprocess]
    Frequencies at Gamma point:
    q = [0. 0. 0.] (weight= 1)
    # Mode   Frequency
        1   -0.0000002 THz
        2   -0.0000001 THz
        3    0.0000002 THz
        4   15.7649077 THz
        5   15.7649078 THz
        6   15.7649079 THz
    ```
This will:

- Compute the phonon bandstructure along [high symmetry paths in the Brillouin zone](https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#high-symmetry-paths) and save it in `phonopy/output/bandstructure.pdf`.
- Compute the density of states using a $45 \times 45 \times 45$ $\bf q$ point grid and the Tetrahedron method.
  The density of states will be plotted alongside the bandstructure to a file `output/bandstructure_dos.pdf`, and written to a data file [`total_dos.dat`](https://phonopy.github.io/phonopy/output-files.html#total-dos-dat-and-projected-dos-dat).
  The q-grid can be adjusted by specifying it with an additional flag `--q_mesh`.
- Compute the harmonic free energy $F^{\rm ha}$ and the harmonic heat capacity at constant volume, $C_V$, i.e., the thermal properties accessible in the harmonic approximation using the DOS and it q-point settings.
  An overview plot is saved to `output/thermal_properties.pdf` and the detailed output is written to [`output/thermal_properties.yaml`](https://phonopy.github.io/phonopy/output-files.html#thermal-properties-yaml).
- Create animation files for visualization with [`v_sim`](http://www.mem-lab.fr/en/Pages/L_SIM/Softwares/V_Sim.aspx).
- Write a `phonopy.yaml` for [loading a `Phonopy` object directly within `python`](https://phonopy.github.io/phonopy/phonopy-module.html#shortcut-to-load-input-files-phonopy-load).

??? info "Bandstructure"
	![image](bandstructure.png)

**Congratulations!** You have just performed a full (but not yet converged!) _ab initio_ phonon bandstructure calculation.

Note that the CLI also allows to only run a subset of the postprocessing, e.g.,
```
vibes output phonopy phonopy/trajectory.son -v -bs
```
only outputs the bandstructure.

## Choosing a supercell size

!!! info
	The ideal supercell size and shape depends on your problem at hand and it is difficult to give definite advice. In practice, the supercell size needs to be converged until the target property of interest is not changing anymore.
        To facilitate this, there is a CLI tool that can help you creating supercells of different sizes.

There is a [CLI utility](../Documentation/cli.md#vibes-utils)  in`FHI-vibes` that can help you to find supercells of different sizes:

```
vibes utils make-supercell
```

For example

```
vibes utils make-supercell geometry.in -n 8
```

will find the conventional, cubic cell of silicon with 8 atoms:

```
...
Settings:
  Target number of atoms: 8

Supercell matrix:
 python:  [-1,  1,  1,  1, -1,  1,  1,  1, -1]
 cmdline: -1 1 1 1 -1 1 1 1 -1
 2d:
[[-1, 1, 1],
 [1, -1, 1],
 [1, 1, -1]]

Superlattice:
[[5.42906529 0.         0.        ]
 [0.         5.42906529 0.        ]
 [0.         0.         5.42906529]]

Number of atoms:  8
  Cubicness:         1.000 (1.000)
  Largest Cutoff:    2.715 AA
  Number of displacements: 1 (1)

Supercell written to geometry.in.supercell_8
```

It will tell you the supercell matrix that you can use in `phonopy.in` (`python:  [-1,  1,  1,  1, -1,  1,  1,  1, -1]`), the generated superlattice, a "cubicness" score based on the filling ratio of the largest sphere fitting into the cell, the largest cutoff in which any neighbor is not a periodic image of a closer neighbor to estimate boundary effects, and the number of supercells with displacements that  `phonopy` will create. It will also write the structure to `geometry.in.supercell_8` which you can inspect, e.g., with `jmol`.

To run a calculation for such a supercell, one just needs to replace the `supercell` keyword in the `[phonopy]` section is the required value:

```
...
[phonopy]
supercell_matrix: [-1,  1,  1,  1, -1,  1,  1,  1, -1]
...
```

Remember to use a new working directory to not mess up your previous results!

!!! info
	The force constants obtained for the `[-1,  1,  1,  1, -1,  1,  1,  1, -1]` supercell with 8 atoms are re-used in later tutorials, please don't delete them.

### Practical guideline

In practice, the convergence with supercell size needs always to be checked carefully, since it depends on the range of the interactions present in your system, e.g., long ranged unscreened van-der-Waals interactions require larger supercells than short-ranged covalent ones as here in Si. Along the same lines, the acceptable supercell size depends also on the properties you are interested in. Free energies and specific heats converge faster than individual frquenecies.  Using a cubic-as-possible supercell shape and playing around with `vibes utils make-supercell` and a little bit of experience will do the job. For the example at hand, it might for instance be instructive to check the convergence of different properties for larger supercell sizes. You find calculations for supercell sizes up 1728 atoms in our [tutorial repository](https://gitlab.com/vibes-developers/vibes-tutorial-files/-/tree/master/2_phonopy). When do you consider the calculations as converged?# Rapid Prototyping with Empirical Force-Fields

The aim of this tutorial is to learn how to use force-fields within FHI-vibes. For this purpose, we use [the Lennard-Jones Argon test case](0_intro.md#test-systems) at $20\,{\rm K}$
and perform the exact same calculation steps discussed for LDA-Silicon in the previous tutorial. Since force evaluations for such a toy system are order of magnitudes faster compared to _ab initio_
methods, this allows to quickly test and illustrate the influence of varios computational parameter on the MD. For instance, we suggest to test the workflow below for various supercell
sizes, temperatures, timesteps, etc.

!!! warning

	- This tutorial mimics the essential steps for performing MD simulations in bulk systems. How you obtain initial structures in your project is, of course, highly dependent on the systems you aim to study etc.
	- This tutorial needs ASE 3.20, which is not yet released. You can install the most recent ASE version with `pip install https://gitlab.com/ase/ase/-/archive/master/ase-master.tar.gz`.

##  Structure preparation

### Generate a structure

Copy the [the Argon structure](0_intro.md#test-systems) to a file called `geometry.in.primitive`. From this primitive cell, we can use the CLI to generate a supercell of about 100 atoms:

```
vibes utils make_supercell geometry.in.primitive -n 100
```

this will try to find a cubic-as-possible supercell with roughly 100 atoms and write it to an input file.

??? info "Output  of `vibes utils make_supercell geometry.in.primitive -n 100`"
    ```
    Find supercell for
    [vibes]        Geometry info
      input geometry:    Ar
      Symmetry prec.:    1e-05
      Number of atoms:   1
      Species:           Ar (1)
      Periodicity:       [ True  True  True]

      Spacegroup:          Fm-3m (225)
      Wyckoff positions:   1*a
      Equivalent atoms:    1*0
    
    Settings:
      Target number of atoms: 100
    
    Supercell matrix:
     python:  [-3,  3,  3,  3, -3,  3,  3,  3, -3]
     cmdline: -3 3 3 3 -3 3 3 3 -3
     2d:
    [[-3, 3, 3],
     [3, -3, 3],
     [3, 3, -3]]
    
    Superlattice:
    [[15.78  0.    0.  ]
     [ 0.   15.78  0.  ]
     [ 0.    0.   15.78]]
    
    Number of atoms:  108
      Cubicness:         1.000 (1.000)
      Largest Cutoff:    7.890 AA
    
    Supercell written to geometry.in.primitive.supercell_108
    ```

In this case it should find a perfectly cubic supercell with 108 atoms and write it to `geometry.in.primitive.supercell_108`. You can obtain detailed information about the structure by running

```
vibes info geometry geometry.in.primitive.supercell_108
```

??? info "Output of `vibes info geometry geometry.in.primitive.supercell_108`"
    ```
    [vibes]        Geometry info
      input geometry:    Ar
      Symmetry prec.:    1e-05
      Number of atoms:   108
      Species:           Ar (108)
      Periodicity:       [ True  True  True]
      Lattice:
        [15.78  0.    0.  ]
        [ 0.   15.78  0.  ]
        [ 0.    0.   15.78]
      Cubicness:         1.000 (1.000)
      Largest Cutoff:    7.890 AA

      Spacegroup:          Fm-3m (225)
      Wyckoff positions:   108*a
      Equivalent atoms:    108*0
    
    Cell lengths and angles [Å, °]:
      a, b, c:     15.7800     15.7800     15.7800
      α, β, γ:     90.0000     90.0000     90.0000
      Volume:             3929.3526 Å**3
      Volume per atom:      36.3829 Å**3
    ```

Additionally you can inspect the generated structure with the structure viewer of your choice, e.g., with [`jmol`](http://jmol.sourceforge.net/).

??? info "`jmol geometry.in.primitive.supercell_108`"
	![image](assets/geometry.in.supercell.png)


Assuming that we are happy with this structure, we save it as our supercell:

```
mv geometry.in.primitive.supercell_108 geometry.in.supercell
```

### Pre-thermalize the structure

To speed up the thermalization, we can pre-thermalize the system by giving momenta to the atoms according to a Maxwell-Boltzmann distribution at our target temperature of $20\,{\rm K}$. This can be done with the CLI utility `create-samples`:

```
vibes utils create-samples geometry.in.supercell -T 20
```

??? info "Output of `vibes utils create-samples geometry.in.supercell -T 20`"
    ```
    vibes CLI: create_samples
    [vibes]        Geometry info
      input geometry:    Ar
      Symmetry prec.:    1e-05
      Number of atoms:   108
      Species:           Ar (108)
      Periodicity:       [ True  True  True]

      Spacegroup:          Fm-3m (225)
      Wyckoff positions:   108*a
      Equivalent atoms:    108*0
    [vibes]        Geometry info
      input geometry:    Ar
      Symmetry prec.:    1e-05
      Number of atoms:   108
      Species:           Ar (108)
      Periodicity:       [ True  True  True]
    
      Spacegroup:          Fm-3m (225)
      Wyckoff positions:   108*a
      Equivalent atoms:    108*0
    [vibes]        Use Maxwell Boltzamnn to set up samples
    [vibes]        Sample   0:
    [vibes]        .. temperature before cleaning:    18.119K
    [vibes]        .. remove net momentum from sample and force temperature
    [vibes]        .. temperature in sample 0:        20.000K
    [vibes]        Sample   0:
    [vibes]        .. temperature in sample 0:        20.000K
    [vibes]        .. written to geometry.in.supercell.0020K
    ```
The geometry written to `geometry.in.supercell.0020K` will now include the appropriate velocities.

We will use this structure and the chosen velocities as the initial structure for the MD run. We suggest to rename this file to `geometry.in` accordingly.

## Run MD

### Prepare `md.in`

Before we can run the MD, we need to create an input file. To this end, copy the calculator section [for the Lennard-Jones calculator](0_intro.md#lj-argon) to a file called `md.in`. Next, we use the CLI command `template` to add settings for performing a NVT simulation:

```
vibes template md --nvt >> md.in
```

??? info "The generated `md.in`"

    ```
    [calculator]
    name:                          lj
    
    [calculator.parameters]
    # parameters for LJ Argon
    sigma:    3.405
    epsilon:  0.010325
    rc:       13.0


    [md]
    driver:                        Langevin
    timestep:                      1
    maxsteps:                      1000
    compute_stresses:              False
    workdir:                       md
    
    [md.kwargs]
    temperature:                   300
    friction:                      0.02
    logfile:                       md.log
    ```

We suggest to add and/or adjust the following  keywords:

```
[md]
timestep:                      4
maxsteps:                      7500

[md.kwargs]
temperature:                   20

[files]
geometry:                      geometry.in
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell
```

The `timestep` can be increased to $4\,{\rm fs}$ for Argon at $20\,{\rm K}$. With `maxsteps: 7500` we will run a total of 7500 MD steps, i.e., $30\,{\rm ps}$ simulation time. **The total simulation time depends on the system and the quantitiy of interest!**`temperature` should be set to $20\,{\rm K}$, our target temperature.
Adding `primitive: geometry.in.primitive` and `supercell: geometry.in.supercell` in the `[files]` section is not necessary to run the calculation. However, `vibes` will automatically attach this information to the trajectory so that it cannot get lost. This also makes life easer when post processing. For example, the displacements $\Delta {\bf R}_I(t)$ can only be properly calculated, when the reference supercell is known.

The final `md.in` should look like this:

```
[calculator]
name:                          lj

[calculator.parameters]
# parameters for LJ Argon
sigma:    3.405
epsilon:  0.010325
rc:       8.0

[md]
driver:                        Langevin
timestep:                      4
maxsteps:                      7500
compute_stresses:              False
workdir:                       md

[md.kwargs]
temperature:                   20
friction:                      0.02
logfile:                       md.log

[files]
geometry:                      geometry.in
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell
```

We are now ready to  run the simulation!

### Run the calculation

You can run this calculation with the CLI command `run`. We recommend to save it's output, e.g., with `tee`:

```
vibes run md | tee md.log
```

Depending on you computer, the calculation will take a few minutes.

## Create trajectory dataset and inspect the simulation

### Process the calculation

The data obtained at each time step will be written to the trajectory file `md/trajectory.son`. The CLI provides a tool to process the trajectory and create an `xarray.Dataset` from it. To this end, run

```
vibes output md md/trajectory.son
```

This command will create `trajectory.nc`, a dataset representation of the data contained in the MD trajectory saved as an `xarray.Dataset` to a NetCDF file. The included data can be viewed with

```
vibes info netcdf trajectory.nc
```

??? info "Output of `vibes info netcdf trajectory.nc`"

    ```
    <xarray.Dataset>
    Dimensions:              (I: 108, a: 3, b: 3, time: 2501)
    Coordinates:
      * time                 (time) float64 0.0 4.0 8.0 ... 9.996e+03 1e+04
    Dimensions without coordinates: I, a, b
    Data variables:
        positions            (time, I, a) float64 ...
        displacements        (time, I, a) float64 ...
        velocities           (time, I, a) float64 ...
        momenta              (time, I, a) float64 ...
        forces               (time, I, a) float64 ...
        energy_kinetic       (time) float64 ...
        energy_potential     (time) float64 ...
        stress               (time, a, b) float64 ...
        stress_kinetic       (time, a, b) float64 ...
        stress_potential     (time, a, b) float64 ...
        temperature          (time) float64 ...
        cell                 (time, a, b) float64 ...
        positions_reference  (I, a) float64 ...
        lattice_reference    (a, b) float64 ...
        pressure             (time) float64 ...
        pressure_kinetic     (time) float64 ...
        pressure_potential   (time) float64 ...
    Attributes:
        name:             trajectory
        system_name:      Ar
        natoms:           108
        time_unit:        fs
        timestep:         4.000000000000006
        nsteps:           2500
        symbols:          ['Ar', 'Ar', 'Ar', 'Ar', 'Ar', 'Ar', 'Ar', 'Ar', 'Ar', ...
        masses:           [39.948 39.948 39.948 39.948 39.948 39.948 39.948 39.94...
        atoms_reference:  {"pbc": [true, true, true],\n"cell": \n[[ 1.57800000000...
        volume:           3929.352552000002
        raw_metadata:     {"MD": {\n  "type": "molecular-dynamics",\n  "md-type":...
        hash:             0eff05aa63cd4019927c42af74bb0ff0a0e21009
    ```

### View simulation statistics

To get information about the simulation, you can use the CLI command `info md`, which summarizes the simulation and can produce an overview plot:

```
vibes info md trajectory.nc -p
```

This command should tell you, among other things, that the temperature is indeed thermalized to approximately $ 20\,{\rm K}$:

```
...
[info]         Summarize Temperature
Simulation time:            30.000 ps (7501 steps)
Temperature:                    19.826 +/-       1.7902 K
Temperature (1st 1/3):          19.402 +/-       2.2027 K
Temperature (2st 1/3):          20.426 +/-       1.4744 K
Temperature (3st 1/3):          19.651 +/-       1.4218 K
Temperature (last 1/2):         19.863 +/-       1.4301 K
...
```

The pdf file `md_summary.pdf` provides a visualization of the simulated properties for quick sanity checking that  the simulation went according to plan:

??? info "`md_summary.pdf`"
	![image](assets/md_summary.png)

### Visualize trajectory
The trajectory can be exported to an `xyz` file for visualizing the atomic motion, e.g., with [`VMD`](https://www.ks.uiuc.edu/Research/vmd/). To this end, run

```
vibes utils trajectory 2xyz trajectory.nc
vmd trajectory.xyz
```

??? info "`vmd trajectory.xyz`"
	![image](assets/LJ-Argon.gif)

# Phonons

!!! info
	We assume that you are familiar with the basics of phonon theory and you are here to learn how to perform them with `FHI-vibes`. We give a short recap on the background below.

## <a name="Phonons"></a> Phonons: Harmonic vibrations in solids

To determine the vibrations in a solid, we approximate the potential energy surface
for the nuclei by performing a Taylor expansion of the total energy $\mathcal V({\bf R})$ around the equilibrium positions:

$$
\require{cancel}
\begin{aligned}
\def\vec#1{{\mathbf{#1}}}
\def\t#1{\text{#1}}
\mathcal V \left(\{\vec{R}^0 + \Delta \vec{R}\}\right)
& \approx
\mathcal V\left(\{\vec{R}^0\}\right) \\
& + \cancel{ \sum\limits_{I} \left.\frac{\partial \mathcal V}{\partial \vec{R}_I}\right\vert_{\vec{R}^0} \Delta\vec{R}_{I} } \\
& + \frac{1}{2} \sum\limits_{I,J} \left.\frac{\partial^2 \mathcal V}{\partial \vec{R}_I\partial \vec{R}_J}\right\vert_{\vec{R}^0} \Delta\vec{R}_{I}\Delta\vec{R}_{J} \\
& + \mathcal{O}(\Delta\vec{R}^3)
\end{aligned}
$$

The linear term vanishes, since no forces $\vec{F} = - \nabla \mathcal V$ are acting on the system in equilibrium $\vec{R}^0$.
Assessing the Hessian $\Phi_{IJ} = \frac{\partial^2 \mathcal V}{\partial \vec{R}_I\partial \vec{R}_J}$ involves some additional
complications: In contrast to the forces $\vec{F}$, which only depend on the density, the Hessian $\Phi_{IJ}$ also depends
on its derivative with respect to the nuclear coordinates, i.e., on its _response_ to nuclear displacements. One can either
use _Density Functional Perturbation Theory (DFPT)_ [[Baroni2001](references.md#baroni2001)] to compute the response
or one can circumvent this problem by performing the second order derivative _numerically by finite differences_

$$
\begin{align}
\Phi_{IJ}
= \left.\frac{\partial^2 \mathcal V}{\partial \vec{R}_I\partial \vec{R}_J}\right\vert_{\vec{R}^0}
= - \left.\frac{\partial }{\partial \vec{R}_I} \vec{F}_J\right\vert_{\vec{R}^0}
= - \lim_{\epsilon \rightarrow 0} \frac{ \vec{F}_J(\vec{R}_I^0 + \epsilon \,\vec{d}_I)}{\epsilon}~,
\label{eq:FinDiff}
\end{align}
$$

as we will do in this tutorial.
The definition in Eq.$~\eqref{eq:FinDiff}$ is helpful to realize that the Hessian describes a coupling between different atoms, i.e., how the force acting on an atom $\vec{R}_J$ changes
if we displace atom $\vec{R}_I$, as you have already learned in tutorial 1. However, an additional complexity arises in the case of _periodic boundary conditions_,
since beside the atoms in the unit cell $\vec{R}_J$ we also need to account for the periodic images $\vec{R}_{J'}$. Accordingly, the Hessian is in principle a matrix
of infinite size. In non-ionic crystals, however, the interaction between two atoms$~I$ and $J$ quickly decays with their distance$~\vec{R}_{IJ}$, so that we can compute the Hessian from
finite supercells, the size convergence of which must be accurately inspected.

Once the real-space representation of the Hessian is computed, we can determine the _dynamical matrix_ by adding up the contributions from all periodic images$~J'$ in the mass-scaled Fourier transform of the Hessian:

$$
\begin{align}
D_{IJ}({\vec{q}}) = \sum\limits_{J'}
\frac{\t e^{\t{im} \left({\vec{q}}\cdot{\vec{R}_{JJ'}}\right)}}{\sqrt{M_I M_J}}
\;\Phi_{IJ'}
\quad .
\label{DynMat}
\end{align}
$$

In reciprocal space [[AshcroftMermin](references.md#AshcroftMermin)],
this _dynamical matrix_ determines the equation of motion for such a periodic array of harmonic
atoms for each reciprocal vector$~\vec{q}$:

$$
\begin{align}
D(\vec{q}) \, \vec e_s (\vec{q}) = \omega_s^2(\vec{q}) \, \vec e_s (\vec{q})
\; .
\label{eq:eigenproblem}
\end{align}
$$

The dynamical matrix has dimension $3N_\t{A} \times 3N_\t{A}$, where $N_\t{A}$
is the number of
atoms in the *primitive* unit cell. Equation$~\eqref{eq:eigenproblem}$ thus
constitutes
an eigenvalue problem with $3N_\t{A}$ solutions at each $\vec q$ point. The
solutions
are labelled by $s$ and are denoted as *phonon branches*. The lowest three
branches are commonly called the *acoustic branches*, whereas in solids
with more than one atom in the primitive unit cell, the remaining $(3N_\t{A} -
3)$ branches are denoted as *optical branches*.

The eigenvalues$~\omega_s^2(\vec{q})$ and eigenvectors$~\vec e_s(\vec{q})$ of the
dynamical matrix$~D(\vec{q})$
completely describe the dynamics of the system (in the harmonic approximation), which is nothing else than a superposition
of harmonic oscillators, one for each mode, i.e., for each
eigenvalue$~\omega_s (\vec{q})$.

From the set of eigenvalues $\{ \omega_s (\vec{q}) \}$, also denoted as
spectrum, the *density of states* (DOS) can be obtained by
calculating the number of states in an infinitesimal energy window
$[\omega, \omega + \t d \omega]$:

$$
\begin{align}
g(\omega) = \sum_s\int\frac{\t d \vec{q}}{(2\pi)^3}\delta(\omega -
\omega(\vec{q})) = \sum_s\int\limits_{\omega(\vec{q}) =
\omega}\frac{\t{d}S}{(2\pi)^3}\frac{1}{\vert\nabla\omega(\vec{q})\vert}~.
\label{DOS}
\end{align}
$$

The DOS is a very useful quantity, since it allows to determine any integrals
(the integrand of which only depends on $\omega$) by a
simple integration over a one-dimensional variable$~\omega$ rather than a three-dimensional variable$~\vec{q}$. This is much
easier to handle both in numerical and in analytical models. For instance, we can compute the associated thermodynamic
potential[^footnote1], i.e., the (harmonic) *Helmholtz free
energy*[^footnote2]

$$
\begin{equation}
F^{\mathrm{ha}}(T,V)  = \int \t d \omega\; g(\omega) \left (
\frac{\hbar\omega}{2} + k_\t{B}\, T~
\ln\left(1-\t{e}^{-\frac{\hbar\omega}{k_\t{B}\,T}}\right)
\right)\;.
\label{HFVT}
\end{equation}
$$

In turn, this allows to calculate the heat capacity at constant
volume

$$
\begin{equation}
C_V = - T \left(\frac{\partial^2 F^{\mathrm{ha}}(T,V)}{\partial T^2} \right)_V \;.
\label{CV}
\end{equation}
$$

For a comprehensive introduction to the field of lattice dynamics and its
foundations, we refer to [[BornHuang](references.md#BornHuang)].

To compute the quantities introduced above, we will use
`FHI-vibes`, which uses the package _phonopy_ [[Togo2015](references.md#Togo2015)] as a backend to compute vibrational properties via the finite-displacements method as outlined above. Please note that
_phonopy_ makes extensive use of symmetry
analysis [[Parlinski1997](references.md#Parlinski1997)], which allows to reduce numerical noise and to speed up the calculations considerably. Our system of choice will be fcc-diamond Silicon (you can run the tutorial as well with [LJ-Argon](0_intro.md#test-systems-for-the-tutorials)).

!!! warning
    In the following exercises, the computational settings, in particular the reciprocal space grid (tag `k_grid`), the basisset and supercell sizes, have been chosen to allow a rapid computation of the exercises. In a _real_ production calculation, the reciprocal space grid, the basis set, and the supercells would all have to be converged with much more care,  although the qualitative trends hold already with the present settings.

## Recap on solid state physics

### Brillouin Zone

As you know, the periodicity of the atomic positions
in a lattice is reflected by the fact that it is sufficient to look at $\bf q$
values within the first Brillouin zone of the reciprocal lattice: A wave vector
$\bf q$ that lies outside the first Brillouin zone corresponds to a wave whose
wavelength is _shorter_ than the distance between periodic images of the
same atom. It can thus be represented equally well by a wave with longer
wavelength, i.e. a smaller wave vector ${\bf q}'$ taken from within the
first [Brillouin zone](https://en.wikipedia.org/wiki/Phonon\#Crystal_momentum).

### High Symmetry Points
The Brillouin zone of our Silicon fcc diamond structure is displayed below

??? info "Plot of Brillouin zone if fcc lattice"
	![image](assets/BZ_fcc.png)

The labelled points correspond to $\bf q$
values of high symmetry. This means that there are symmetry operations in the
point group of the lattice that leave this point invariant (up to a reciprocal space vector)
[p. 218 in Dresselhaus].

You can list the high symmetry points of the lattice of your geometry with `vibes` by
running

```
vibes info geometry geometry.in -v
```

??? info "list  of high symmetry points"
    ```
    ...
    Special k points:
    G: [0. 0. 0.]
    K: [0.375 0.375 0.75 ]
    L: [0.5 0.5 0.5]
    U: [0.625 0.25  0.625]
    W: [0.5  0.25 0.75]
    X: [0.5 0.  0.5]
    ```

Please note that the list of points is given in
fractional coordinates as coefficients of the _reciprocal_ lattice. For
the meaning of the Symbols $\Gamma$ (G), $X$, etc., you can take a look at [the Wikipedia article](https://en.wikipedia.org/wiki/Brillouin_zone#Critical_points).

### Bandstructure
If we connect two or more $\bf q$ points from the Brillouin zone, solve the eigenvalue
problem for any $\bf q$ point in between, and plot the
obtained dispersions $\omega (q)$ versus $q$, we obtain the so-called phonon bandstructure. The bandstructure is typically computed for a path in the Brillouin
zone that connects several or all of the high symmetry points.



[^footnote1]: Given that the _Bose-Einstein distribution_ is used for
the derivation of the harmonic free energy in this case, we get the correct quantum-mechanical result including zero-point effects by this means.

[^footnote2]: A derivation of Eq.$\,\eqref{HFVT}$ will be presented in Tutorial 7 on
Molecular Dynamics Simulations.

# _ab initio_ Molecular Dynamics

!!! info
	We will now introduce _ab initio_ Molecular Dynamics simulations where the atomic forces come from a first-principles code. We will use `FHI-aims` as this calculator in the following. We assume:

	- You are familiar with running MD simulations in a different context and are here to learn how to run MD with `FHI-vibes`.
	- You are familiar with running `FHI-aims` calculations.
	- Optionally: You are familiar with running `FHI-aims` calculations on a workstation or computing cluster.

_Ab initio_ molecular dynamics simulations are MD simulations where the forces are computed from first principles, e.g., using density functional theory (DFT) with LDA or GGA xc-functional in the Born-Oppenheimer approximation. Thus,

$$
\begin{align}
\mathcal V ({\bf R}) = E_{\rm tot}^{\rm DFT} ({\bf R})~,
\label{eq:PES}
\end{align}
$$

where the potential energy $\mathcal V({\bf R})$ of a given atomic configuration $\bf R$ is given by the total energy  $E_{\rm tot}^{\rm DFT} ({\bf R})$ of the electronic and nuclear system computed for this structure[^footnote1].
Given that a single DFT calculations is required for each step of the trajectory, MD calculations are computationally more expensive than the calculations performed in the other tutorials.

## Setting up ab initio MD

We will use an 8 atoms supercell of silicon with LDA xc-functional as previously in the tutorial on geometry optimization and phonon calculations. You can re-use the structure from [the tutorial on geometry optimization](1_geometry_optimization.md), as well as the calculator setup.

```
cp ../path/to/relaxation/geometry.in.next_step geometry.in.primitive
```

Create a supercell with 8 atoms:

```
vibes utils make-supercell geometry.in.primitive -n 8
mv geometry.in.primitive.supercell_8 geometry.in.supercell
```

To speed up the thermalization, we pre-thermalize the supercell by giving the kinetic energy corresponding to a temperature of $300\,{\rm K}$ with the following command:
```
vibes utils create-samples geometry.in.supercell -T 300
mv geometry.in.supercell.0300K geometry.in
```

### Prepare `md.in`

Before we can run the MD, we need to create an input file. To this end, copy the calculator section [for LDA-Si](0_intro.md#lda-silicon) to a file called `md.in`. Next, we use the CLI command `template` to add settings for performing a NVT simulation:

```
vibes template md --nvt >> md.in
```

??? info "The generated `md.in`"

    ```
    [calculator]
    name:                          aims

    [calculator.parameters]
    xc:                            pw-lda

    [calculator.kpoints]
    density:                       2

    [calculator.basissets]
    default:                       light

    [calculator.socketio]
    port:                          12345

    [md]
    driver:                        Langevin
    timestep:                      1
    maxsteps:                      1000
    compute_stresses:              False
    workdir:                       md

    [md.kwargs]
    temperature:                   300
    friction:                      0.02
    logfile:                       md.log
    ```

We suggest to add and/or adjust the following  keywords:

```
[md]
timestep:                      4
maxsteps:                      2500
compute_stresses:              10

[files]
geometry:                      geometry.in
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell
```

The `timestep` can be increased to $4\,{\rm fs}$ for Silicon at $300\,{\rm K}$. With `maxsteps: 2500` we will run a total of 2500 MD steps, i.e., $10\,{\rm ps}$ simulation time. **The total simulation time depends on the system and the quantitiy of interest!**`temperature` should be set to $300\,{\rm K}$, our target temperature.
The flag `compute_stresses` in the section `[md]` will make FHI-aims compute the _ab initio_ stress every 10 steps during the MD simulation. This will provide access to pressure.[^footnote2]
by inspecting the Adding `primitive: geometry.in.primitive` and `supercell: geometry.in.supercell` in the `[files]` section is not necessary to run the calculation. However, `vibes` will automatically attach this information to the trajectory so that it cannot get lost. This also makes life easer when post processing. For example, the displacements $\Delta {\bf R}_I(t)$ can only be properly calculated, when the reference supercell is known.

The final `md.in` should look like this:

```
[calculator]
name:                          aims

[calculator.parameters]
xc:                            pw-lda

[calculator.kpoints]
density:                       2

[calculator.basissets]
default:                       light

[calculator.socketio]
port:                          12345

[md]
driver =           Langevin
timestep =         4
temperature =      300
friction =         0.02
maxsteps =         2500
compute_stresses = 10

[files]
geometry:                      geometry.in
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell
```

We are now ready to  run the simulation!

## Run a calculation

This step is similar to [before](2_phonopy.md#run-the-calculation), i.e., you run

```
vibes run md >> log.md &
```

??? info "`log.md"

    ```
    [vibes.run]    run MD workflow with settings from md.in

    [md]           driver: Langevin
    [md]           settings:
    [md]             type: molecular-dynamics
    [md]             md-type: Langevin
    [md]             timestep: 0.39290779153856253
    [md]             temperature: 0.02585199101165164
    [md]             friction: 0.02
    [md]             fix-cm: True
    [md]           ** /scratch/usr/becflokn/vibes/tutorial/3_md/ab_initio/si_8/md/trajectory.son does not exist, nothing to prepare
    [calculator]   Update aims k_grid with kpt density of 2 to [4, 4, 4]
    [calculator]   .. add `sc_accuracy_rho: 1e-06` to parameters (default)
    [calculator]   .. add `relativistic: atomic_zora scalar` to parameters (default)
    [calculator]   .. add `compensate_multipole_errors: False` to parameters (default)
    [calculator]   .. add `output_level: MD_light` to parameters (default)
    [calculator]   Add basisset `light` for atom `Si` to basissets folder.
    [calculator]   Calculator: aims
    [calculator]   settings:
    [calculator]     xc: pw-lda
    [calculator]     k_grid: [4, 4, 4]
    [calculator]     sc_accuracy_rho: 1e-06
    [calculator]     relativistic: atomic_zora scalar
    [calculator]     compensate_multipole_errors: False
    [calculator]     output_level: MD_light
    [calculator]     compute_forces: True
    [calculator]     compute_heat_flux: True
    [calculator]     use_pimd_wrapper: ('localhost', 12345)
    [calculator]     aims_command: run_aims
    [calculator]     species_dir: /scratch/usr/becflokn/vibes/tutorial/3_md/ab_initio/si_8/md/basissets
    [socketio]     Use SocketIO with host localhost and port 12345
    [backup]       /scratch/usr/becflokn/vibes/tutorial/3_md/ab_initio/si_8/md/calculations does not exists, nothing to back up.
    Module for Intel Parallel Studio XE Composer Edition (version 2019 Update 5) loaded.
    Module for Intel-MPI (version 2018.5) loaded.
    [md]           Step 1 finished, log.
    [md]           switch stresses computation off
    [md]           Step 2 finished, log.
    [md]           switch stresses computation off
    [md]           Step 3 finished, log.
    [md]           switch stresses computation off
    [md]           Step 4 finished, log.
    [md]           switch stresses computation off
    [md]           Step 5 finished, log.
    ...
    ```

For running on a cluster, see [additional remarks](0_singlepoint.md#submit-calculation-on-a-cluster).

## Postprocess

### Process the calculation

The data obtained at each time step will be written to the trajectory file `md/trajectory.son`. The CLI provides a tool to process the trajectory and create an `xarray.Dataset` from it. To this end, run

 ```
vibes output md md/trajectory.son
 ```

This command will create `trajectory.nc`, a dataset representation of the data contained in the MD trajectory saved as an `xarray.Dataset` to a NetCDF file. The included data can be viewed with

```
vibes info netcdf trajectory.nc
```

??? info "Output of `vibes info netcdf trajectory.nc`"

    ```
    <xarray.Dataset>
    Dimensions:              (I: 8, a: 3, b: 3, time: 2501)
    Coordinates:
      * time                 (time) float64 0.0 4.0 8.0 ... 9.996e+03 1e+04
    Dimensions without coordinates: I, a, b
    Data variables:
        positions            (time, I, a) float64 ...
        displacements        (time, I, a) float64 ...
        velocities           (time, I, a) float64 ...
        momenta              (time, I, a) float64 ...
        forces               (time, I, a) float64 ...
        energy_kinetic       (time) float64 ...
        energy_potential     (time) float64 ...
        stress               (time, a, b) float64 ...
        stress_kinetic       (time, a, b) float64 ...
        stress_potential     (time, a, b) float64 ...
        temperature          (time) float64 ...
        cell                 (time, a, b) float64 ...
        positions_reference  (I, a) float64 ...
        lattice_reference    (a, b) float64 ...
        pressure             (time) float64 ...
        pressure_kinetic     (time) float64 ...
        pressure_potential   (time) float64 ...
    Attributes:
        name:             trajectory
        system_name:      Si
        natoms:           8
        time_unit:        fs
        timestep:         4.000000000000006
        nsteps:           2500
        symbols:          ['Si', 'Si', 'Si', 'Si', 'Si', 'Si', 'Si', 'Si']
        masses:           [28.085 28.085 28.085 28.085 28.085 28.085 28.085 28.085]
        atoms_reference:  {"pbc": [true, true, true],\n"cell": \n[[ 5.42906529316...
        atoms_primitive:  {"pbc": [true, true, true],\n"cell": \n[[-0.00000000000...
        atoms_supercell:  {"pbc": [true, true, true],\n"cell": \n[[ 5.42906529316...
        volume:           160.02034201861315
        raw_metadata:     {"MD": {\n  "type": "molecular-dynamics",\n  "md-type":...
        hash:             ff1410ec05dc89c85cf148670ecb05947a0066c8
    ```

### Inspect results

You can perform postprocessing of the pressure by [inspecting the trajectory dataset in `trajectory.nc`](../Documentation/output_files.md#trajectorync). Be aware that the simulation time is shorter when discarding the thermalization period.

??? info "reference pressure"
	For 8 atoms LDA-Silicon, you should get a potential pressure of $-0.61 \pm 0.06 {}$

[A more detailed introduction to postprocessing including example scripts is given in the next chapter.](3_md_postprocess.md)

## References
Running the calculation will take some time depending on the computer your working with. You find references [in our reference repository](https://gitlab.com/vibes-developers/vibes-tutorial-files/-/tree/master/3_molecular_dynamics/ab_initio). There you also find reference calculations for 64 and 216 atoms.

[^footnote1]: Note that a corect _ab initio_ MD requires also to choose the different numerical settings in the DFT calculation with care. For instance, the inherent incompleteness of the SCF cycle in Kohn-Sham DFT schemes can introduce a systematic error that introduces energy drifts and other unphysical effects. Choosing the correct convergence settings is an important aspect of performing _ab initio_ MD simulations and is highly materials specific. Devising a strategy on how to choose these settings goes beyond the scope of this tutorial.

[^footnote2]: We compute the stress only every 10th step because computing the stress is numerically more expensive in _ab initio_ than computing the atomic forces only. Since we know that consecutive samples generated during MD are highly correlated, we don't loose valuable information by computing this quantity not in every single step.
# Getting started with the high-throughput version of FHI-vibes
FireWorks is a high-throughput framework that can greatly enhance the efficiency of your work, but does require initial setup to get working.
Because of this, it is important to consider if running calculations in high-throughput is necessary for your work, or if using the command line interface tools vibes provides is sufficient.
If you only need to do an in-depth study of a few materials then high-throughput aspects of FHI-vibes is probably not necessary, but if you need to do a systematic study over many materials this could be a useful tool for your research.

This setup guide will only focus on the basic aspects that you need do in order to run our FireWorks workflows, for more advanced usage please look at the [FireWorks online documentation](https://materialsproject.github.io/fireworks/).


## Setting up MongoDB
Before running or installing FireWorks you first need to have access to a MongoDB instance to act as the job management database (LaunchPad in FireWorks terminology).
To do this there are two options: set up your own database by installing MongoDB locally or on an accessible server, or by using a cloud provider.
If you want to install MongoDB yourself please follow [their instructions](https://www.mongodb.com/), for a cloud provider you can try [mlab](http://mlab.com/), [MongoDB Atlas](https://www.mongodb.com/cloud/atlas), or another option.

Once you have created a database you'll need to set up a user to access it.
Furthermore it is recommended that you also set up secure passwords for all your databases for data security reasons.
To set this up we suggest you follow this guide, but other methods are possible.
We have a separate admin/user structure in order to ensure that admin access to the database can't be obtained through a FireWorks configuration file.

### How to create a new MongoDB with authorization
*Note: This is for a locally managed database, if you are using a cloud service security procedures maybe different please consult their instructions*

First create a database with the correct base directory (\$BASE_DIR), port (\$PORT), and binding IP addresses (\$BIND_IP) with the following command
```
mongod --logpath $BASE_DIR/logs --dbpath $BASE_DIR/db --port $PORT --bind_ip $BIND_IP --fork
```
It is good practice to include 127.0.0.1 (localhost) as one of the binding IP's along with the one you will use to access the database from the outside.
The list of IP addresses can be stored in $BIND_IP with the following command
```
export BIND_IP=127.0.0.1,IP_1,IP_2,...
```
note it is a comma separated list.

From here you can access the database with
```
mongo --port $PORT --host $HOST_NAME
```
\$HOST_NAME should be one of the IP addresses in \$BIND_IP, if localhost is on the list and you are on that server, the host keyword is not necessary.
You should now be in a MongoDB terminal, from here type in the following commands
```
use admin
db.createUser(
  {
    user: "admin_user",
    pwd: "ADMIN_PASSWORD",
    roles: [ { role: "userAdminAnyDatabase", db: "admin" }, "readWriteAnyDatabase" ]
  }
)
exit
```
This will create an admin user that has global access to all databases on the server. Most importantly now you can restart your database requiring authorization to access it.
```
mongod --dbpath $BASE_DIR/db --shutdown
mongod --logpath $BASE_DIR/logs --dbpath $BASE_DIR/db --port $PORT --bind_ip $BIND_IP --fork --auth
```
Now that the database is running you can create additional users/databases inside the main one with the following commands.
To start launch the MongoDB terminal
```
mongo --port $PORT --host $HOST_NAME -u "admin_user" -p "ADMIN_PASSORD" --authenticationDatabase "admin"
```
then inside that terminal type the following commands
```
use FIREWORKS_DB_NAME
db.createUser(
  {
    user: "USER",
    pwd: "PASSWORD",
    roles: [ { role: "dbOwner", db: "FIREWORKS_DB_NAME" } ]
  }
)
exit
```
Please do not use the same username/password for the admin user and individual users.
Now that you have a database running you can store this information in a yaml file called my_launchpad.yaml with the following contents
```
host: HOST
port: PORT
name: FIREWORKS_DB_NAME
username: USER
password: PASSWORD
```
The host this time should be the IP address in $BIND_IP that you intend to use for external connections. If this is on the same machine you will be running the workflows you maybe able to use localhost, but check with your sysadmin on how the nodes are set up.
Now that you have a database set up now it is time to install FireWorks

## Installing FireWorks
The first step in using FireWorks is installing the python library.
By installing vibes with the FireWorks dependency with `pip install fhi-vibes[fireworks]` this is already included, but you may want to install your own version of FireWorks locally.
You can do this either by cloning the [git repository](https://github.com/materialsproject/fireworks) and using the python setup tools or via pip/conda.
In addition to FireWorks, if you want to use the remote clients/database you'll need to install paramiko and fabric or if you are using the NEWT queuing system you'll have to install requests.
To do all of this with pip simply type in
```
pip install FireWorks
pip install paramiko  # (only needed if using built-in remote file transfer!)
pip install fabric  # (only needed if using daemon mode of qlaunch!)
pip install requests  # (only needed if you want to use the NEWT queue adapter!)
```
The fireworks extension of vibes already includes paramiko and fabric, but requests would have to be installed via pip.
As a note to install fireworks using pip the Kerberos 5 develop package must be installed.
Please ensure this is installed for your system.

Once FireWorks is installed you can test your installation by attempting to connect to FireWorks read-only test database, by creating `my_launchpad_testing.yaml` with the following contents:
```
host: ds049170.mongolab.com
port: 49170
name: fireworks
username: test_user
password: testing123
```
Then you can run the following command

```
lpad -l my_launchpad_testing.yaml get_wflows
```
And you should get the following output
```
[
    {
        "name": "Tracker FW--1",
        "state": "READY",
        "states_list": "REA",
        "created_on": "2014-10-27T15:00:25.408000"
    },
    {
        "name": "Tracker FW--2",
        "state": "READY",
        "states_list": "REA",
        "created_on": "2014-10-27T15:00:25.775000"
    }
]
```

Now that FireWorks is installed properly it is time to set up your configuration.
Create a fireworks configuration directory (\$FW_CONFIG) to store all the configuration files.
We recommend you use `.fireworks/` for \$FW_CONFIG.
Move the my_launchpad.yaml file to \$FW_CONFIG.
If you are planning to use FireWorks with a queuing system also create a my_qadapter.yaml file in \$FW_CONFIG. To get a sample of what to do for your queue system go to the [FireWork's git hub](https://github.com/materialsproject/fireworks/tree/master/fw_tutorials/queue) page and download the correct one.
For example here is one for a SLURM System
```
_fw_name: CommonAdapter
_fw_q_type: SLURM
rocket_launch: vibes fireworks rlaunch singleshot
nodes: 1
ntasks_per_node: NUMBER OF CORES PER CPU
walltime: MAXIMUM WALL TIME
queue: BATCHING PARTITION
account: ACCOUNT TO TAKE CPU TIME FROM
job_name: DEFAULT JOB NAME
logdir: LOG_FILES
pre_rocket: COMMAND TO RUN BEFORE EACH JOB (if None use null)
post_rocket: COMMAND TO RUN AFTER EACH JOB (if None use null)
```
Finally you can also specify a FWorker by creating an my_fworker.yaml file in \$FW_CONFIG, for example:
```
name: my first fireworker
category: ''
query: '{}'
```
Once all files are in your \$FW_CONFIG folder edit PYTHON_SITE_PACKAGES_DIRECTORY/fireworks/fw_config.py to reflect where to find the correct FireWorks configuration files. The relevant portion changes to the file  should look like this with $FW_CONFIG replaced with the correct path
```
LAUNCHPAD_LOC = $FW_CONFIG/my_launchpad.yaml  # where to find the my_launchpad.yaml file
FWORKER_LOC = $FW_CONFIG/my_fworker.yaml  # where to find the my_fworker.yaml file
QUEUEADAPTER_LOC = $FW_CONFIG/my_qadapter.yaml  # where to find the my_qadapter.yaml file

CONFIG_FILE_DIR = $FW_CONFIG
```
To find where the site-packages file is located run `python -m site` and it should appear in the returned list.
```
[optimize_kgrid]
dfunc_min = 1e-3
```

## Sections

### `[optimize_kgrid]`

Section to optimize k-grid for all `FHI-aims` calculation

#### `dfunc_min`

`float`: Minimum change in total energy with changing k-point density```
[section.qadapter]
nodes = 1
ntasks_per_node = 20
walltime = 1-00:00:00
queue = partition_name
account = account_name
```

## Sections

### `[section.qadapter]`

A subsection of the high throughput workflows that specifies the queuing parameters for each job.

#### `nodes`

`int`: The number of nodes requested (Default defined in the `my_qadapter.yaml` file)

#### `ntasks_per_node`

`int`: The number of tasks per node (Default defined in the `my_qadapter.yaml` file)

#### `walltime`

`str`: The requested wall time in `DD-HH:MM:SS` format (Default defined in the `my_qadapter.yaml` file)

#### `queue`

`str`: The partition to submit the job to (Default defined in the `my_qadapter.yaml` file)

#### `account`

`str`: The account to charge for the calculation (Default defined in the `my_qadapter.yaml` file)
```
[relaxation]
use_aims_relax = True
method = trm
fmax = 1e-3
relax_unit_cell = full

[relaxation.1]
basis = light


[relaxation.2]
basis = intermediate

.
.
.

[relaxation.n]
basis = really_tight
```

## Sections

### `[relaxation]`

Sections to do relaxation of structures. This is a general definition for various step

#### `use_aims_relax`

`bool`: True if you want to use the aims relaxation instead of the one defined in the [`relaxation`](../../../Documentation/relaxation) portion of the documentation. For this documentation we will only show keywords if this is true. If this is false then use the keywords you'd normally use for a relaxation

#### `basis`

`str`: keyword for the basis set to use for the relaxation step

#### `method`

`str`: Relaxation method used for the calculation

#### `fmax`

`float`: Maximum residual force before ending the relaxation

#### `relax_unit_cell`

`str`: How to relax the unit cell within `FHI-aims` either `full`, `fixed_angles` or `none`

### `[relaxation.1]`

The first step of the relaxation. Only define parameters that are different from default parameters


### `[relaxation.2]`

The second step of the relaxation (parameters same as the first step)

### `[relaxation.n]`

The n<sup>th</sup> step of the relaxation

```
[phonopy]
supercell_matrix = [2, 2, 2]
walltime = 3500
serial = True
displacement = 0.01

[phonopy.convergence]
minimum_similarity_score = 0.80
sc_matrix_base = [1, 1, 1]

[gruneisen]
volume_factors = [0.99, 1.01]

```
## Sections

### `[phonopy]`

Parameters for phonopy calculations. Most keywords are the same. For full documentation of those see the `phonopy` section in the documentation

#### `serial`

`bool`: If True use serial calculations instead of parallel (calculate all supercells in one calculation v. separately) (Default is True)

#### `convergence`

`bool`: If True do phonon supercell convergence with the defaults defined in `phonopy.convergence` section (Default is False)

### `[phonopy.convergence]`

Section used to define the phonon convergence parameters. If both sets of defaults are desired then

#### `minimum_similarity_score`

`float`: Minimum Tanimoto similarity score to consider the phonon calculations converged with respect supercell size (Default is 0.80)

#### `sc_matrix_base`

`list(int)`: Base supercell matrix use to increase supercell size (Default is `phonopy.supercell_matrix`).

In the example above the next supercell matrix tested would be `[3, 3, 3]` without it the next would be `[4, 4, 4]`

### `[gruneisen]`

Set up to calculate the gruneisen parameters from finite difference using phonopy

#### `volume_factors`

`list(float)`: ratio of the equilibrium volume to calculate the phonons for.```
[md]
phonon_file = path/to/phonopy/trajectory.son
temperatures = [300, 600]
supercell_matrix = [1, 1, 1]
{The rest of the MD parameters described in [md]}
```

## Sections

### `[md]`

#### `phonon_file`

`str`: The trajectory file used for generating the thermally displaced structures

#### `supercell_matrix`:

`list(int)`: The supercell matrix for the calculation, if not given use the one from the phonopy calculation. If the supercell matrix is different from the one in the `phonon_file` the phonopy force constants will be remapped onto the new supercell.

#### `temperatures`:

`list(float)`: list of temperatures to calculate the anharmonicity at
```
[fireworks]
name = example_run
config_dir: "~/.fireworks"
tasks2queue = ["vibes.relaxation.bfgs.relax", "vibes.fireworks.tasks.calculate_wrapper.wrap_calc_socket", "vibes.k_grid.converge_kgrid.converge_kgrid", "vibes.fireworks.tasks.calculate_wrapper.wrap_calculate", "vibes.fireworks.tasks.md.run" ]

[fireworks.workdir]
remote = "test_run/"
local   = "test_run/"

[fireworks.remote]
host = ["remote.host.path"]
config_dir = ["/path/to/remote/home/.fireworks/"]
launch_dir = "."

[fireworks.remote.authorization]
user = remote_username
password = null

[fireworks.remote.launch]
njobs_queue = 0
njobs_block = 500
reserve = True
nlaunches = 0
sleep_time = 60
```

## Sections

### `[fireworks]`
General parameters for the FireWorks workflows

#### `name`
`str`: The name that will be perpended the workflow to better organize the LaunchPad

#### `config_dir`
`str`: Directory where FireWorks configuration file are located (Default set in `.fireworksrc` file)

#### `tasks2queue`
`list(str)`: List of functions to send to the queue (Default set in `.fireworksrc` file)

### `[fireworks.workdir]`

These are used to define the base working directory on remote and local machines

#### `local`

`str`: Base working directory on your local machine

#### `remote`

`str`: Base working directory on a remote directory (Default is `fireworks.workdir.local`)

### `[fireworks.remote]`

Parameters for remote FireWorks workers (Default defined in `.fireworksrc`)

#### `host`

`list(str)`: List of remote hosts to send jobs to (Default defined in `.fireworksrc`)

#### `config_dir`

`list(str)`: List of remote FireWorks configuration directories (Default defined in `.fireworksrc`)

#### `launch_dir`

`str`: Default launch directory on the remote host (Default defined in `.fireworksrc`)

### `[fireworks.remote.authorization]`

Parameters for authentication remote FireWorks workers (Default defined in `.fireworksrc`)

#### `user`

`str`: remote host username (Default defined in `.fireworksrc`)

#### `password`

`str`: remote host password (not recommended) (Default defined in `.fireworksrc`)

### `[fireworks.remote.launch]`

Parameters for launching jobs on FireWorks workers (Default defined in `.fireworksrc`)

#### `njobs_queue`

`int`: number of jobs to have on the queue at any given time (0 no limit) (Default defined in `.fireworksrc`)

#### `njobs_block`

`int`: number of launches to have in a single FireWorks block directory (Default defined in `.fireworksrc`)

#### `reserve`

`bool`:  If True run FireWorks in reservation mode (Default defined in `.fireworksrc`)

#### `nlaunches`

`int`: Maximum number of jobs to launch at any given (0 no limit) (Default defined in `.fireworksrc`)

#### `sleep_time`

`float`: Time to sleep in seconds between checking for jobs to run (Default defined in `.fireworksrc`)
```
[statistical_sampling]
phonon_file = path/to/phonopy/trajectory.son
supercell_matrix = [-1,1,1,1,-1,1,1,1,-1]
temperatures = [300, 600]
debye_temp_fact = [1.0]
serial = True
n_samples = 1
plus_minus = True
mc_rattle = False
quantum = True
deterministic = True
zacharias = True
gauge_eigenvectors = True
ignore_negative = False
failfast = True
random_seed = 13
propagate = False
```

## Sections

### `[statistical_sampling]`

Used for Monte Carlo sampling of a system for anharmonicity quantification

#### `phonon_file`

`str`: The trajectory file used for generating the thermally displaced structures

#### `supercell_matrix`:

`list(int)`: The supercell matrix for the calculation, if not given use the one from the phonopy calculation. If the supercell matrix is different from the one in the `phonon_file` the phonopy force constants will be remapped onto the new supercell.

#### `temperatures`:

`list(float)`: list of temperatures to calculate the anharmonicity at

#### `debye_temp_fact`:

`list(float)`: list of multipliers to add temperatures that are factors the materials Debye temperature

#### `serial`:

`bool`: If True then do this in serial (Default is True)

#### `n_samples`:

`int`: number of samples to calculate for each temperature (Default is 1)

#### `plus_minus`:

`bool`: Use the deterministic sampling regime from Zacharias, et al (Default is True)

#### `deterministic`:

`bool`: If True populate all phonon modes to k_B T energy (Default is True)

#### `gauge_eigenvectors`:

`bool`: If True use a plus minus gauge for the eigenmodes (Default is True)

#### `mc_rattle`:

`bool`: If True rattle the structures using a Monte Carlo rattling method (Default is False)

#### `quantum`:

`bool`: If True populate phonon modes according to a Bose-Einstein distribution (Default is False)

#### `ignore_negative`:

`bool`: If True ignore all imaginary modes (Default is False)

#### `failfast`:

`bool`: If True If True fail if any imaginary modes are present or acustic modes are not near zero at Gamma (Default is True)

#### `random_seed`:

`int`: The seed for random number generator (Default is Random number)

#### `propagate`:

`bool`: If True propagate the structure forward in time somewhat with ASE (Default is False)This section defines the high-throughput definitions for all the keywords for the following tasks:
- [Installing FireWorks Dependency](../../Installation/0_setup)
- [Setting up a high-throughput workflow](../1_general_high_throughput)
- [qadapters](../2_qadapter)
- [K-Point Density Convergence](../3_optimize_k_grid)
- [Relaxation](../4_relaxation)
- [Phonon Calculations](../5_phonons)
- [Monte Carlo Sampling](../6_statistical_sampling)
- [Molecular Dynamics](../7_md)

At the moment only FHI-aims is supported for high-throughput calculations. It is possible to extend it to other calculators, but certain transitions between steps assume the calculator is for FHI-aims and those would have to change.<a name="Running Multiple Phonopy Calculations"></a>

??? info "Prerequisite"
    Please complete the [high-throughput phonon tutorial](1_phonopy.md) before completing this one

## Summary
In this section we will learn how to set up and run multiple multi-step workflows using FHI-vibes and FireWorks.

We will start by learning how to converge the `k_grid` setting for FHI-aims and optimize the structure of Si before performing phonopy calculations.
We then go on to define how we can perform harmonic and molecular dynamics based sampling, and use the harmonic sampling quantify a material's anharmonicity.

## Optimizing the k-point Density and Geometry

Now that basic phonopy calculations have been explained, let's make the workflows more useful by adding two preprocessing tasks to the workflow before the phonopy calculations: k-grid convergence and geometry optimization.

### k-grid Optimization
```
[optimize_kgrid]
dfunc_min:                     1e-3

[optimize_kgrid.qadapter]
nodes:                         1
walltime:                      01:00:00
```
Because a single k-point density may not give converged yet efficient  results for every material in a test set, we provide the option to converge the k-point density with respect to the total energy for all materials.
This task will increase the k-point density defined in `calculator.kpoints.density` by 1.0 until successive increases of the density results in a change of energy less than `optimize_kgrid.dfunc_min` eV.
It will then use that k-point density for all future calculations in the workflow.

In this example the convergence criteria is set to 1 meV to decrease computational time, but this can be changed to any desired level of accuracy.
The default value is 0.001 meV.

### Relaxation

```
[relaxation.1]
basis:                         light
driver:                        BFGS
fmax:                          0.001
unit_cell:                     True
decimals:                      12
maxstep:                       0.2
fix_symmetry:                  True

[relaxation.1.qadapter]
nodes:                         1
walltime:                      00:15:00
```
Optimizing the geometries for high-throughput calculations are largely similar to the command line interface, but with numerically numbered sections.
Instead of using `[relaxation]` and `[relaxation.kwargs]` as done in the command line interface both of these sections are combined into `[relaxation.1]`, `[relaxation.2]`, ..., `[relaxation.n]` sections.
This is done to allow for a multi-stage relaxation, where the basis set or unit cell restrictions can be changed to improve the performance of the algorithms.
To accommodate the potential change in basis each relaxation step can have its own basis set that is different from the main calculation, set by the `relaxation.n.basis` keyword.
Each of these relaxation steps will also have their own `qadapter` as they all may have different costs associated with them.

### The Full Workflow
Below is the complete workflow file for the calculations, we removed the phonon convergence and only studying Si to save time.

??? info "`workflow.in`"
    ```
    [files]
    geometry:                      Si/geometry.in

    [fireworks]
    name:                          example_multistep_calculations

    [fireworks.workdir]
    local:                         analysis/
    remote:                        run/

    [calculator]
    name:                          aims
    socketio:                      True

    [calculator.parameters]
    xc:                            pw-lda

    [calculator.kpoints]
    density:                       1

    [calculator.basissets]
    default:                       light

    [phonopy]
    supercell_matrix:              [-2, 2, 2, 2, -2, 2, 2, 2, -2]
    serial:                        True

    [phonopy.qadapter]
    nodes:                         1
    walltime:                      00-01:00:00

    [optimize_kgrid]
    dfunc_min:                     1e-3

    [optimize_kgrid.qadapter]
    nodes:                         1
    walltime:                      01:00:00

    [relaxation.1]
    basis:                         light
    driver:                        BFGS
    fmax:                          0.001
    unit_cell:                     True
    decimals:                      12
    maxstep:                       0.2

    [relaxation.1.qadapter]
    nodes:                         1
    walltime:                      00:15:00
    ```

### Running the Workflow
We can now run the calculations as we did previously in the [phonopy tutorial](1_phonopy.md).
```
vibes fireworks add_wf
vibes fireworks rlaunch rapidfire
```
Once this is completed, we can then compare the results from these workflows and those from using the wrong structure.

### Analyzing the Results
Looking inside the `run/` directory there are the following directories:
```
run/Si/49833c381b84708fbcd174c47a777478ab5dec26:
1_relax  kgrid_opt  sc_natoms_64
```
`kgrid_opt` and `1_relax` store information about the k-point density and the first (and only) step of geometry optimization respectively.
By looking at the logs we see that the k-point density is already converged at a density of 1.0, and Si needs only slight optimizations to its starting geometries.

The `analysis/` directory only has information about the `phonopy` calculations because the optimizations do not have any post-processing steps that need to be stored locally.
Compare the DOS and bandstructure to those calculated in the [previous tutorial](1_phonopy.md) by running

```
vibes output phonopy -bs --dos
```
These results look similar to what was seen previously, with only slight differences in the bandstructures and density of states because there was only a slight change in the geometry upon relaxation; however, it is always important to relax all geometries before running phonopy to get the correct harmonic model.
Here is what the bandstructure_dos.pdf should look like:
![Silicon Bandstructure and DOS](images/Si_ms.png)


## Using the Harmonic Model to Quantify Anharmonicity

Once the a converged harmonic model is generated, the anharmonicity of the material can be quantified using the methods published [here](https://arxiv.org/abs/2006.14672) and previously described [here](../../Tutorial/5_anharmonicity_quantification.md).

### Harmonic Sampling

The `[statistical_sampling]` section calculates the anharmonicity of a material from a harmonic sampling of a its harmonic vibrational potential energy surface.

This task uses the harmonic model to generate a series of thermally displaced supercells as discussed [here](../../Tutorial/4_statistical_sampling.md).

Here is the section will be used to calculate $\sigma^\text{A}$ for Si at 300 and 600 K.
```
[statistical_sampling]
phonon_file:                   analysis/Si/0df71cea3a5446b7104554b9bada4da6eb4a802a/sc_natoms_64/phonopy_analysis/trajectory.son
serial:                        True
temperatures:                  [300, 600]
supercell_matrix:              [-1, 1, 1, 1, -1, 1, 1, 1, -1]
n_samples:                     1
plus_minus:                    True
```
In order to not rerun the phonopy calculation of the relaxed structure, the previously calculated `trajectory.son` file is passed via the `phonon_file` keyword.
If the workflow has a `[phonopy]` section in it, then it will automatically use the `trajectory.son` file from that calculation for all materials in the workflow; however, if `phonon_file` needs to be explicitly passed then the workflow can only be used for that material.
Because the `supercell_matrix` used in this calculation is not the same as the one used to calculate the harmonic model, the force constants will be remapped onto the new supercell.
If `supercell_matrix` is not provided then it will default to the one in `phonon_file`.

In this example, `statistical_sampling.plus_minus` is `True` so the scheme developed by [Zacharias and Giustino](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.075125) is used to generate the displaced supercells.
This scheme replaces the random scaling by deterministic one that approximates a thermal ensemble, where $\mathbf{d_I}$ is

\begin{equation}
    d_{I, \alpha} = \sum_{s=3}^{n_{w}} \frac{(-1)^s d_{s, I, \alpha}}{\sqrt{m_I}}.
\end{equation}

If this scheme is not used, then `statistical_sampling.n_samples` should be large enough to ensure convergence for all materials tested (normally on the order of 10-30 samples).
For a complete description of possible harmonic sampling keywords, have a look at the [documentation](../Documentation/6_statistical_sampling.md).

### Molecular Dynamics (MD)

```
[md]
phonon_file:                   run/Si/0df71cea3a5446b7104554b9bada4da6eb4a802a/sc_natoms_64/phonopy/trajectory.son
supercell_matrix:              [-1, 1, 1, 1, -1, 1, 1, 1, -1]
driver:                        Langevin
timestep:                      1
temperatures:                  [300, 600]
friction:                      0.02
maxsteps:                      50
logfile:                       md.log
```

The settings in the `[md]` section used for high-throughput calculations are almost identical to those introduced in [the molecular dynamics tutorial](../../Tutorial/3_md_ab_initio.md), but the initial structures are calculated automatically from phonopy calculation in the previous step or from the trajectory file defined in `md.phonon_file`.
If you are not chaining the `[phonopy]` and `[md]` steps together please ensure the `phonon_file` is on the resource that you intend to run the calculation on (most likely the remote host not your local machine).

As is the case with the `[statistical_sampling]` section, if the requested `supercell_matrix` is not the same as the one in `phonon_file` the force constants will be remapped onto the new supercell, and if no `supercell_matrix` is given it will default to the one in `phonon_file`.

If you want to run MD for multiple temperatures replace `md.temperature` with `md.temperatures` and give a list of temperatures and multiple Molecular Dynamics tasks will be generated, as is done above for the harmonic sampling example.

### The Full Workflow
Below is the complete workflow file for the calculations.
The molecular dynamics calculations are not run here, since they are too expensive
See [the tutorial](../../Tutorial/3_md_postprocess.md) for postprocessing options for MD.

??? info "`workflow_anharmonicty.in`"
    ```
    [files]
    geometry:                      Si/geometry.in

    [fireworks]
    name:                          example_anharmonicity_calculations

    [fireworks.workdir]
    local:                         analysis/
    remote:                        run/

    [calculator]
    name:                          aims
    socketio:                      True

    [calculator.parameters]
    xc:                            pw-lda

    [calculator.kpoints]
    density:                       1

    [calculator.basissets]
    default:                       light


    [statistical_sampling]
    phonon_file:                   analysis/Si/0df71cea3a5446b7104554b9bada4da6eb4a802a/sc_natoms_64/phonopy_analysis/trajectory.son
    serial:                        True
    temperatures:                  [300, 600]
    supercell_matrix:              [-1, 1, 1, 1, -1, 1, 1, 1, -1]
    n_samples:                     1
    plus_minus:                    True

    [statistical_sampling.qadapter]
    walltime:                      0:30:00
    nodes:                         1

    [md]
    phonon_file:                   run/Si/0df71cea3a5446b7104554b9bada4da6eb4a802a/sc_natoms_64/phonopy/trajectory.son
    supercell_matrix:              [-1, 1, 1, 1, -1, 1, 1, 1, -1]
    driver:                        Langevin
    timestep:                      1
    temperatures:                  [300, 600]
    friction:                      0.02
    maxsteps:                      5
    logfile:                       md.log
    compute_stresses:              5

    [md.qadapter]
    walltime:                      01:00:00
    nodes:                         1
    ```

For both `[statistial_sampling]` and `[md]` the structure that is used is the one in their respective `phonon_file`, not the one defined in `files.geometry`, but it is important to keep the reference file there so the correct atoms hashed is used.

### Running the Workflow
We can now run the calculations as we did previously, but with passing the workflow file name in explicitly.
```
vibes fireworks add_wf -w workflow_anharmonicty.in
vibes fireworks rlaunch rapidfire
```
Once this is completed, we can see what the anharmonicity of Silicon is.

## Analyzing the Results
Looking inside the `run/` directory we now see the following directories:
```
run/Si/49833c381b84708fbcd174c47a777478ab5dec26:
1_relax  kgrid_opt  sc_natoms_64 statistical_sampling
```
where we now see the statistical_sampling calculation directories.

Inside the `analysis/` directory there are now two folders:
```
analysis/Si/49833c381b84708fbcd174c47a777478ab5dec26:
sc_natoms_64 statistical_sampling_analysis
```
In `statistical_sampling_analysis` there are two files: `trajectory.son` and `sigma.dat`.
The `sigma.dat` file has $\sigma^\text{A}$ value at the requested 300 and 600 K.
Looking at the file with, e.g., `cat`, we see that for Si $\sigma^\text{A}$ is quite low, confirming that it is very harmonic
```
cat analysis/Si/0df71cea3a5446b7104554b9bada4da6eb4a802a/statistical_sampling_analysis/sigma.dat
300.0, 0.1312671429298908
600.0, 0.18100681527473686
```
# Configuring FHI-vibes for interaction with FireWorks
Once FireWorks is setup and running  a few additional steps are needed to configure FHI-vibes to work with FireWorks.
If you are having trouble installing FireWorks see [their documentation](https://materialsproject.github.io/fireworks/) or [our installation guide](../Installation/0_setup.md).

## Creating a `.fireworksrc` file
The first step in configuring FHI-vibes for use with FireWorks.
This file sets up default values for running FireWorks utilities that will be consistent throughout all calculation.
To make this file run:

`vibes template configuration fireworks > ~/.fireworksrc`

```
[fireworks]
config_dir: "~/.fireworks" # Directory containing the *yaml files for FireWorks
tasks2queue = ["vibes.relaxation.bfgs.relax", "vibes.fireworks.tasks.calculate_wrapper.wrap_calc_socket", "vibes.k_grid.converge_kgrid.converge_kgrid", "vibes.fireworks.tasks.calculate_wrapper.wrap_calculate", "vibes.fireworks.tasks.md.run" ] # DO NOT CHANGE

[fireworks.remote]
host = ["remote.host.path"] # List of remote host names
config_dir = ["/path/to/remote/home/.fireworks/"] # List of remote FireWorks configuration directories
launch_dir = "." # Default launch directory on the remote host

[fireworks.remote.authorization]
user = remote_username # remote host username
password = null # remote host password (not recommended try to use password-less login)

[fireworks.remote.launch]
njobs_queue = 0 # Number of jobs to have on the queue at any given time (0 no limit)
njobs_block = 500 # Number of launches to have in a single FireWorks block directory
reserve = True # If True run FireWorks in reservation mode
nlaunches = 0 # Maximum number of jobs to launch at any given (0 no limit)
sleep_time = 60 # Time to sleep in seconds between checking for jobs to run
```
For a complete description of each of these parameters see the [full documentation](../../Documentation/1_general_high_throughput).

## Testing if it works
Now that your FireWorks installation should be working properly test it with the vibes FireWorks test in `test/fireworks/test_fireworks.py`.
If the test runs successfully then when you run `lpad get_wflows` you should get the following output (created on should match today's date/time):
```
{
    "state": "COMPLETED",
    "name": "Ni_6d2a2be5a5c1c4549639b55c5403b438b3b0ccf7--1",
    "created_on": "2020-03-13T11:55:30.357000",
    "states_list": "C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C"
}

```
If you see that you have successfully set up the high-throughput portions of vibes. To use this on clusters you need to repeat the steps in Installing/Testing FireWorks on each machine you plan to use it on.
<a name="Running Multiple Phonopy Calculations"></a>

??? info "Prerequisite"
    FHI-vibes configured for use with FireWorks. See the [configuration guide](0_configuring_fw_for_vibes.md) for more information.


## Summary
In this section we will learn how to set up and run multiple phonopy calculations using FHI-vibes and FireWorks.
In the workflow we will also describe how we systematically determine if the supercell size is converged for a given material, so all materials will be calculated to the same level of precision.
As an example we will use Si and MgO for this tutorial, but it can be extended to any set of materials.

## Setup workflow.in file

Setting up a high-throughput workflow to perform multiple `phonopy` calculation is similar to setting up a single calculation, but with a few additional steps to ensure the calculations are done in a similar manner.
Because the high throughput workflows are designed to be flexible, there is no `vibes template` command to automatically generate them, but modifying the workflows you'll work with here and in [the multi-step tutorial](2_multistep.md) should be good guide on how to get started.

In this case only two materials (Si-diamond and MgO-rock salt) will be calculated, but the workflow can be used to generate a harmonic model for an arbitrary number of materials.
As before, we will not fully converge the results for these examples, so to allow for rapid execution and testing.

We start from already relaxed structures for silicon and magnesium oxide and store them in `Si/geometry.in` and `MgO/geometry.in`, respectively. *Note: Typically, one does not have relaxed structures already available. Information on how to perform a relaxation before running the phonon calculations (i.e. a multistep workflow) will be given in [the multistep tutorial](2_multistep.md).*

??? info "`Si/geometry.in`"
    ```
    lattice_vector      0.000      2.703      2.703
    lattice_vector      2.703      0.000      2.703
    lattice_vector      2.703      2.703      0.000
    atom_frac      0.00     0.00     0.00 Si
    atom_frac      0.25     0.25     0.25 Si

    ```

??? info "`MgO/geometry.in`"
    ```
    lattice_vector    0.000   2.104   2.104
    lattice_vector    2.104   0.000   2.104
    lattice_vector    2.104   2.104   0.000
    atom_frac    0.00    0.00    0.00  Mg
    atom_frac    0.50    0.50    0.50  O
    ```
Once the geometry files are added, create a phonopy workflow in `workflow.in` with the following contents.

??? info "`workflow.in`"
    ```
    [files]
    geometries:                    */geometry.in

    [fireworks]
    name:                          example_phonon_calculations
    
    [fireworks.workdir]
    local:                         analysis/
    remote:                        run/
    
    [calculator]
    name:                          aims
    
    [calculator.parameters]
    xc:                            pw-lda
    
    [calculator.kpoints]
    density:                       1
    
    [calculator.basissets]
    default:                       light
    
    [calculator.socketio]
    port:                          12345
    
    [phonopy]
    supercell_matrix:              [-2, 2, 2, 2, -2, 2, 2, 2, -2]
    displacement:                  0.01
    is_diagonal:                   False
    is_trigonal:                   False
    is_plusminus:                  auto
    symprec:                       1e-05
    q_mesh:                        [45, 45, 45]
    serial:                        True
    
    [phonopy.convergence]
    minimum_similarity_score:      0.05
    sc_matrix_base:                [-1, 1, 1, 1, -1, 1, 1, 1, -1]
    
    [phonopy.qadapter]
    nodes:                         1
    walltime:                      00-04:00:00
    ```

This `workflow.in` is very similar to the `phonopy.in` file from the [phonopy tutorial](../../../Tutorial/2_phonopy), but with a few extra sections and keywords.

### FireWorks sections

The largest difference between this workflow and the `phonopy.in` from the command line interface is the additional sections defining FireWorks specific parameters.
The purpose of these sections is to organize both the LaunchPad and file structure of the machines running the workflows.
`[fireworks.workdir]` specifies the base working directory for where the workflows will run (`fireworks.workdir.remote`) and where all the post-processing will happen (`fireworks.workdir.local`).
The naming convention was chosen to reflect that in most cases the electronic structure calculations will be done on a remote cluster, while the postprocessing will be done locally on a laptop or desktop computer.
If `fireworks.workdir.remote` is not set, both post-processing and running the jobs will be done in the base working directory defined in `fireworks.workdir.local`.
Furthermore to ensure the MgO and Si calculations do not overwrite each other these base working directories are appended with `{material_chemical_formula}/{atoms_hash}/` giving each material a unique working directory to store all of its data in.

### qadapter

Another new section in this workflow is `[phonopy.qadapter]`.
This section is used to define job/queue specific information (number of nodes and wallclock time limits in this case), for a more complete definition of what to include here see the [documentation](../../Documentation/2_qadapter).
Defaults for all of these parameters will be stored in the `my_qadapter.yaml` file.

### `[phonopy]`
In the `[phonopy]` section there are only two changes: no `workdir` keyword and the `serial` keyword.
For these workflows all working directories have been standardized, so each task does not require its own `workdir` keyword.
The `serial` keyword controls how FireWorks sets up the force evaluation calculations for phonopy.
If it is True then all calculations will be run serially on a single job, while if it is False then each force evaluation will be done as separate calculations and tasks.
We recommend this is always True to make use of ASE's [socketio calculators](https://wiki.fysik.dtu.dk/ase/ase/calculators/socketio/socketio.html), but it is not necessary to run the calculations.

### Phonopy Convergence

This section is used to converge the harmonic model with respect to the supercell size.
Because `phonopy` uses finite differences to calculate the harmonic approximation, if the supercell is too small an atomic displacement can interact with a periodic image of itself leading to errors.
The magnitude of these errors is material dependent, and is normally either ignored by using what is considered a sufficiently large supercell or checked via a manual inspection of the phonon density of states and band structures.
For high throughput applications manual inspection is impractical so we developed an automatic metric for determining supercell convergence using the [Tanimoto similarity score](https://en.wikipedia.org/wiki/Jaccard_index#Other_definitions_of_Tanimoto_distance) of two successively larger supercells to check convergence.
The score is calculated from the DOS on a 45x45x45 q-grid of the smaller, $\mathbf{d}_{\rm small}$, and larger, $\mathbf{d}_{\rm large}$, supercell

\begin{equation}
    score = \frac{\mathbf{d}_{\rm small} \cdot \mathbf{d}_{\rm large}}{\left|\mathbf{d}_{\rm small}\right|^2 + \left|\mathbf{d}_{\rm large}\right|^2 - \mathbf{d}_{\rm small} \cdot \mathbf{d}_{\rm large} }.
    \label{eq:tanimoto}
\end{equation}

If the similarity score is larger than `phonopy.convergence.minimum_similarity_score` then the harmonic model is considered converged and that part of the workflow ends.
If it is smaller than the threshold then the supercell is increased to
\begin{equation}
    M_\text{S, new} = \left(n\right) M_\text{S, base} + M_\text{0}
    \label{eq:update_phonon}
\end{equation}
where $M_\text{S, base}$ is defined in [`phonopy.convergence.sc_matrix_base`](../../Documentation/5_phonons/#sc_matrix_base), $n$ the current phonopy iteration, and $M_\text{0}$ is the original supercell defined by `phonopy.supercell_matrix`.
$M_\text{0}$ must be an integer scalar value of $M_\text{S, base}$, or the workflow will not be added to the `LaunchPad`.

A score of 0.80 is considered to be a good balance between getting fully converged results and not going to very large supercell sizes; however, you may want to increase it if very accurate results are needed or lower it if the unitcell of a material is already very large.
Here a significantly lower minimum of 0.05 is used so that none of the supercells get above 64 atoms and the workflows can be easily run on a laptop or desktop computer.
Because this supercell matrix used for this phonon model is `[-2, 2, 2, 2, -2, 2, 2, 2, -2]`, the forceconstants for $\mathbf{d}_{\rm small}$ do not have to be calculated explicitly, but are obtained from remapping the force constants of the larger supercell onto the smaller one.
In most cases if a 200 atom supercell is used in this scheme then a converged phonon model can be calculated from a single phonopy calculation.

## Running the Calculations

Unlike what is done with the command line interface running these calculations is done in two steps: first adding the workflow to the LaunchPad and then second running all jobs in the LaunchPad.
To add a calculation to the LaunchPad run:
```
vibes fireworks add_wf
```

and you should get the following output
??? info "Adding Workflow Output"
    ```
    * Message from file vibes/context.py, line 57, function workdir:
    --> workdir not set, return `workdir``

    * Message from file vibes/context.py, line 57, function workdir:
    --> workdir not set, return `workdir``
    
    [calculator]   Update aims k_grid with kpt density of 1 to [4, 4, 4]
    [calculator]   .. add `sc_accuracy_rho: 1e-06` to parameters (default)
    [calculator]   .. add `relativistic: atomic_zora scalar` to parameters (default)
    [calculator]   .. add `output_level: MD_light` to parameters (default)
    [calculator]   Calculator: aims
    [calculator]   settings:
    [calculator]     xc: pw-lda
    [calculator]     k_grid: [4, 4, 4]
    [calculator]     sc_accuracy_rho: 1e-06
    [calculator]     relativistic: atomic_zora scalar
    [calculator]     output_level: MD_light
    [calculator]     compute_forces: True
    [calculator]     use_pimd_wrapper: ('localhost', 12345)
    [calculator]     aims_command: mpiexec -n 1 /home/purcell/git/fhi_aims/bin/ipi.aims.190906.scalapack.mpi.x
    [calculator]     species_dir: /home/purcell/git/fhi_aims/species_defaults/light
    [fireworks]    Generating workflow for MgO
    * Message from file vibes/context.py, line 57, function workdir:
    --> workdir not set, return `workdir``
    
    * Message from file vibes/context.py, line 57, function workdir:
    --> workdir not set, return `workdir``
    
    [calculator]   Update aims k_grid with kpt density of 1 to [4, 4, 4]
    [calculator]   .. add `sc_accuracy_rho: 1e-06` to parameters (default)
    [calculator]   .. add `relativistic: atomic_zora scalar` to parameters (default)
    [calculator]   .. add `output_level: MD_light` to parameters (default)
    [calculator]   Calculator: aims
    [calculator]   settings:
    [calculator]     xc: pw-lda
    [calculator]     k_grid: [4, 4, 4]
    [calculator]     sc_accuracy_rho: 1e-06
    [calculator]     relativistic: atomic_zora scalar
    [calculator]     output_level: MD_light
    [calculator]     compute_forces: True
    [calculator]     use_pimd_wrapper: ('localhost', 12345)
    [calculator]     aims_command: mpiexec -n 1 /home/purcell/git/fhi_aims/bin/ipi.aims.190906.scalapack.mpi.x
    [calculator]     species_dir: /home/purcell/git/fhi_aims/species_defaults/light
    [fireworks]    Generating workflow for Si
    ```

If you are using a non-default launchpad file (normally `my_launchpad.yaml`) or a workflow not defined in `workflow.in` use the `-l` and `-w` flags respectively.

To confirm that both the workflows are in the LaunchPad run
```
lpad get_wflows
```
and you should see:
```
[
    {
        "state": "READY",
        "name": "example_phonon_calculations_Si_92b4ed53a4691c7621606216cab915fa8d5ac311--3",
        "created_on": "YYYY-MM-DDTHH:MM:SS",
        "states_list": "W-REA"
    },
    {
        "state": "READY",
        "name": "example_phonon_calculations_MgO_294b61fcf3f0b3b0c067796c380cade4e2153b30--1",
        "created_on": "YYYY-MM-DDTHH:MM:SS",
        "states_list": "W-REA"
    }
]
```
These workflows each have two jobs: one for setting up the calculation and one for performing the analysis.
The additional analysis step is done for creating a trajectory file if the calculations are running in parallel, and creating a trajectory locally if running with a combined launcher.

Once the workflows have been added to the LaunchPad you have multiple options to run it:

- `vibes fireworks claunch`: Run electronic structure calculations on clusters and everything else locally
- `vibes fireworks qlaunch`: Run all jobs on clusters using the queuing system
- `vibes fireworks rlaunch`: Run all jobs locally
- The FireWorks utilities:

A more detailed description for each running option can be found in their respective documentation.

For this example we'll use `rlaunch`. After running
```
vibes fireworks rlaunch rapidfire
```
the output from
```
lpad get_wflows
```
should now be
```
[
    {
        "state": "COMPLETED",
        "name": "example_phonon_calculations_MgO_294b61fcf3f0b3b0c067796c380cade4e2153b30--1",
        "created_on": "2020-04-10T11:08:43.213000",
        "states_list": "C-C-C-C-C-C"
    },
    {
        "state": "COMPLETED",
        "name": "example_phonon_calculations_Si_92b4ed53a4691c7621606216cab915fa8d5ac311--3",
        "created_on": "2020-04-10T11:08:43.292000",
        "states_list": "C-C-C-C-C-C"
    }

]
```
Four additional tasks have now been added to each workflow corresponding to the force evaluation for the original supercell and three more tasks for the larger supercell for convergence tests.
Now all the workflows have been completed, let's analyze the results.

## Analyzing the Results

The first step in analyzing the results is understanding the file structure of the directories.
First looking at the run directory there are the following folders:
```
run/MgO/294b61fcf3f0b3b0c067796c380cade4e2153b30:
sc_natoms_64

run/Si/92b4ed53a4691c7621606216cab915fa8d5ac311:
sc_natoms_64
```
The  `sc_natoms_64` folders each contain a `phonopy` directory where all the files generated from respectively running $n=2$ `phonopy` iterations.

Looking in the analysis folder there are the following folders
```
analysis/MgO/294b61fcf3f0b3b0c067796c380cade4e2153b30:
converged  sc_natoms_64

analysis/Si/92b4ed53a4691c7621606216cab915fa8d5ac311:
converged  sc_natoms_64
```
The `sc_natoms_64` folders each contain a `phonopy_analysis` directory that only contains the final phonopy `trajectory.son` file for those iterations.
Additionally the last `phonopy` iteration are stored in `converged` for easy access to the converged phonon calculations.

From here you can perform any analysis that is possible within `phonopy` on all the materials, and get results to a similar level of precision for all of them.
For example you can see the bandstructure and DOS of both materials by running
```
vibes output phonopy -bs --dos
```
in each of the `converged` folders.

For Si bandstructure and DOS should look like this
![Silicon Bandstructure and DOS](images/Si_phonopy.png)

And for MgO it should look like this
![Magnesium Oxide Bandstructure and DOS](images/MgO_phonopy.png)

Because of the large variety of possible analysis steps, there is no automated phonopy output scripts in the workflow, but the file structure can be used to easily make bash or python scripts to do all post-processing.

While we can now get converged phonon results, these workflows are incomplete because they require pre-relaxed structures to get physically relevant results.
It would be possible to separately relax all the structures and then use those in this workflow, but the easier solution would be to use [multi-step workflows](2_multistep.md) in the next tutorial.
