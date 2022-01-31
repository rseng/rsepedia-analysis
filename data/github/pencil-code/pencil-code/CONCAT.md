The Pencil Code
---------------

The Pencil Code moved to
[GitHub](https://github.com/pencil-code/pencil-code)
on 19 April 2015. It was previously hosted at
[Google Code](https://code.google.com/p/pencil-code/).

In order to checkout the code with
[Subversion](https://subversion.apache.org), use the command
```sh
svn checkout https://github.com/pencil-code/pencil-code/trunk pencil-code --username <github-username>
```
where `<github-username>` is your GitHub username.

To get started, run one of the samples:
```sh
unix>  cd pencil-code
unix>  source sourceme.csh  [or . sourceme.sh]
unix>  cd samples/conv-slab
unix>  mkdir data
```
To set up the symbolic links and compile the code:
```sh
unix>  pc_setupsrc
unix>  pc_build  [ -f /path/to/config/file.conf ]
```
To create the initial condition and run the code:
```sh
unix>  pc_start  [ -f /path/to/config/file.conf ]
unix>  pc_run    [ -f /path/to/config/file.conf ]
```

See `pencil-code/config/hosts/*/*.conf` for sample config files. For more
details, see the manual in the `doc/` directory (also available
[here](http://pencil-code.nordita.org/)).

-----------------------------------------------------------------------------

If you are using bash and you do not want to "source sourceme.sh" on each
session, you can insert the following into your .bashrc and/or .bash_profile:
```sh
export PENCIL_HOME=$HOME/pencil-code  [or wherever you have the code]
_sourceme_quiet=1; . $PENCIL_HOME/sourceme.sh; unset _sourceme_quiet
```
If you are using csh insert the following into your .cshrc:
```sh
setenv PENCIL_HOME $HOME/pencil-code  [or wherever you have the code]
source $PENCIL_HOME/sourceme.csh
```

## Documentation

* A new documentation webpage has been created with the goal of gather all the documentation together and create auto-documentation for the code:
[Pencil Code documentation](https://pencil-code.readthedocs.io/en/latest/index.html).
* The [manual][manual] is the main source of information.
* There is also a [quick start][quick_start] to help getting started.
* Information about [Python with the Pencil Code][PythonForPencil] and the
  [Python Coding Style][PythonCodingStyle] can be found on the [wiki][wiki].
* Updates to the community are provided through the [newsletter][newsletter].
* Talk to use during the Pencil Code Office Hours every second Friday
  of the month at 13:00 CET via zoom (641-599-5185).
  You may also want to inspect the activity of [Pencil Code Steering Committee][PCSC].
  Please contact them with your ideas.
* The [Pencil Code User Meeting][meetings] will be held every year.
  The 2021 meeting was organized from Switzerland by Jennifer Schober online: [Pencil Code User Meeting 2021](https://www.epfl.ch/labs/lastro/meetings/pcum2021/).
* See the [Scientific Usage of the Pencil Code][citations] for papers using or discussing the code.

## List of Contributors

* Around 100 people have contributed to various extent during the
  nearly 20 years of Pencil Code history.
* The current [list of contributors][contributors] shows the temporal
  check-in activity of the those who stayed connected with the code
  over the various host changes (Nordita 2001-2007, Google Code 2007-2015,
  and Github since 2015).
  Some additional contributors are also listed in the [manual][manual].

## How to contribute to the Pencil Code

* For all changes to the code, make sure the auto-test still runs

* If you have write access: check in your changes and make sure you can
  fix possible problems emerging on [travis-ci.com][travis] as well as the
  minutely, hourly, and daily [auto-tests][auto-tests].

* If you have only read access: fork this repository and use pull requests to contribute.

## Code of Conduct

* The Pencil Code community adheres to the [Contributor Covenant Code of Conduct][conduct].
  Please familiarize yourself with its details.

## License

* The Pencil Code is under the [GNU public license agreement][license].

[travis]: https://www.travis-ci.com/github/pencil-code/pencil-code
[auto-tests]: http://pencil-code.nordita.org/tests.php
[conduct]: https://github.com/pencil-code/pencil-code/blob/master/license/CODE_OF_CONDUCT.md
[manual]: https://github.com/pencil-code/website/raw/master/doc/manual.pdf
[quick_start]: https://github.com/pencil-code/website/raw/master/doc/quick_start.pdf
[license]: https://github.com/pencil-code/pencil-code/blob/master/license/GNU_public_license.txt
[contributors]: https://github.com/pencil-code/pencil-code/graphs/contributors
[wiki]: https://github.com/pencil-code/pencil-code/wiki
[PythonCodingStyle]: https://github.com/pencil-code/pencil-code/wiki/PythonCodingStyle
[PythonForPencil]: https://github.com/pencil-code/pencil-code/wiki/PythonForPencil
[newsletter]: http://www.nordita.org/~brandenb/pencil-code/newsletter
[citations]: https://github.com/pencil-code/website/raw/master/doc/citations.pdf
[PCSC]: https://www.nordita.org/~brandenb/pencil-code/PCSC/
[meetings]: http://pencil-code.nordita.org/meetings.php

# Julia module for the Pencil Code

**Author**:   Daniel Carrera (danielc@astro.lu.se)

**Date**:   Last modified on July 2014.

## Introduction

Julia is a high-level, high-performance dynamic programming language
for technical computing, with syntax that is familiar to users of other
technical computing environments. You can obtain Julia from [http://julialang.org],
or if you are using Ubuntu/Debian, you can install it with `apt-get install julia`.

This is the documentation for the Julia module for the Pencil Code. This module
contains convenience functions for post-processing data files from the Pencil
Code. To use this module, add this to your `~/juliarc.jl` so that Julia can find
the module.

    push!(LOAD_PATH,  ENV["PENCIL_HOME"] * "/julia")


At this point you can load the `Pencil` module:

	ubuntu$ julia        # Start the Julia interpreter.
	...
	julia> using Pencil  # Load the module.

**NOTE:** At the present time, you also need to add the `push!` line at the top of stand-alone programs.


## Plotting and graphics

Julia has several plotting packages. The one I like is
[PyPlot](https://github.com/stevengj/PyPlot.jl) which uses Python's
[Matplotlib](http://matplotlib.org/index.html) library. To install
`PyPlot` run the following:

```
	ubuntu$ sudo apt-get install python-matplotlib
	
	ubuntu$ julia
	...
	julia> Pkg.add("PyPlot")
```

Brief `PyPlot` tutorial:

````python
	julia> using PyPlot
	
	julia> x = linspace(0,2*pi,1000);
	julia> y = sin(3*x + 4*cos(2*x));
	
	julia> plot(x, y, color="red", linewidth=2.0, linestyle="--")

	julia> title("A sinusoidally modulated sinusoid")
	
	#
	# Save figure as a PNG:
	#
	julia> savefig("myfigure.png")
	
	#
	# LaTeX in labels and titles.
	#
	julia> title(L"Plot of $\Gamma_3(x)$")  # L for LaTeX.
	
	#
	# Colour mesh: plot a 2D grid.
	#
	julia> y = [1:128] * ones(128)';  # Col vector x Row vector
	
	julia> r2 = (y - 64).^2 + (y' - 64).^2;
	
	julia> pcolormesh(r2)
	
	julia> axis([0,128,0,128]) # [xmin, xmax, ymin, ymax]
	
	julia> savefig("colmesh.png")
	
	#
	# 3D plotting
	#
	julia> surf(r2)
	
	julia> mesh(r2)
````

The `PencilPlot` module provides a color map called "density" that may
be useful when plotting a density variable like `rhopmxz`. `PencilPlot`
loads `PyPlot`. Example usage:

```
	ubuntu$ cd /my/simulation/dir
	ubuntu$ julia
	...
	julia> using Pencil
	
	julia> using PencilPlot
	
	julia> rhopmxz = read_yaver(it=10)["rhopmxz"];
	
	julia> axis([0, 128, 0, 128])
	
	#
	# PencilPlot defines the "density" color map.
	#
	julia> pcolormesh(rhopmxz, cmap=ColorMap("density") )
	
	julia> savefig("rhopmxz.png")
```




$PENCIL_HOME/python/tests
=========================

Test Pencil Code Python modules, so we can feel better when refactoring the Python code.

```sh
$PENCIL_HOME/python/tests/test-python-modules.py
```

These tests are best run with the [_Proboscis_](https://pythonhosted.org/proboscis/) test runner[^1]:
```sh
pip3 install proboscis
```
but will fall back on a minimal mockup implementation of _Proboscis_ if necessary.


[^1]: The main reason for using _Probioscis_ is test discovery by decorator, as opposed to the standard practice of disovering by name.
python/stubs
============

Stub files (providing type hints) for some external modules.

Use with

```sh
  export MYPYPATH=${PENCIL_HOME}/python/stubs
```

See [the mypy docs](https://mypy.readthedocs.io/en/stable/stubs.html).


 Helical MHD turbulence
==========================

## Maintainer:

Axel Brandenburg <brandenb/nordita[dot]org>

## Added:

08-Jun-2002

## Status:

succeeds

## Recommended resolution:

32x32x32 for nu=eta=5e-3 is fine. For comparison with higher
resolution runs see Brandenburg (2001, ApJ 550, 824), except that
there a forcing wavenumber of kf=5 was used. The forcing function
works with preselected wavevectors that were computed with:
${PENCIL_HOME}/samples/helical-MHDturb/idl/generate_kvectors.pro

## Comments:

This is a helical MHD turbulence run. After about 600 time units
a large scale magnetic field develops from an initially random
magnetic field (if initaa=0 is set).

## Links:
* https://www.nordita.org/~brandenb/projects/LShelicityspec/
* http://pencil-code.nordita.org/samples/turbulence/helical-MHDturb32-4procs/

## References:

*  Brandenburg, A.: 2001 ``The inverse cascade and nonlinear alpha-effect in
   simulations of isotropic helical hydromagnetic turbulence,''
   *Astrophys. J.* **550**, 824-840 |
   [arXiv](http://arXiv.org/abs/astro-ph/0006186) |
   [ADS](http://esoads.eso.org/cgi-bin/nph-bib_query?bibcode=2001ApJ...550..824B)

*  Candelaresi, S., & Brandenburg, A.: 2013, ``Kinetic helicity needed to drive
   large-scale dynamos'' Phys. Rev. E 87, 043104 |
   [arXiv](https://arxiv.org/abs/1208.4529) |
   [ADS](http://adsabs.harvard.edu/abs/2013PhRvE..87d3104C)

 Gravitational Waves from switching on a Beltrami field
=======================================================

## Maintainer:

Axel Brandenburg <brandenb/nordita[dot]org>

## Added:

11-Jun-2018

## Status:

succeeds

## Recommended resolution:

16x16x16 for nu=eta=0 is fine for short times.
For longer runs, one can put nu=eta=1e-3.
To see several oscillations, one should for 400 steps (0.02*400=8).
5-7 us/step/pt on one processor nl6 (laptop), depending on output.

## Comments:

* The gravitational waves are positively circularly polarized.
  Therefore, hel_GWs=GWs and hel_GWh=GWh.
* For EEM=1/2, the amplitudes are EEGW=2pi/4=pi/2; see (28) of Ref.[1]
  and hrms=8pi/2=4pi; see (26) of Ref.[1].
  The spectra denote spec_GWs = S_hdot(k,t), spec_GWh = S_h(k,t),
  and both have amplitude (4pi)^2.
* The time step is fixed to be 0.02 for accuracy reasons (i.e., c*dt/dx=0.051).
  Based on stability alone, the timestep could be 10 times longer [1].
* Alternatively, can use SPECIAL=special/gravitational_waves_hTXk,
  which solves the GW equatin exactly.

## Links:
* https://www.nordita.org/~brandenb/projects/GW

## Reference:
[1] Roper Pol, A., Brandenburg, A., Kahniashvili, T., Kosowsky, A.,
    Mandal, S.: 2018, ``The timestep constraint in solving the
    gravitational wave equations sourced by hydromagnetic turbulence,''
    Geophys. Astrophys. Fluid Dyn., submitted (arXiv:1807.05479v1)

---
title: 'The Pencil Code, a modular MPI code for partial differential equations and particles: multipurpose and multiuser-maintained'

# The Pencil Code is used and developed by the 37 authors, who define the Pencil Code Collaboration.
# About half of the currently 560 papers that acknowledge the code are by others who picked up the code, which has been public since 2001.

tags:
 - Fortran90
 - fluid dynamics
 - magnetohydrodynamics
 - Python
 - IDL
 - astrophysics
 - radiation
 - inertial particles
 - combustion
authors:
 - name: The Pencil Code Collaboration
   affiliation: 1
 - name: Axel Brandenburg
   affiliation: "1, 2, 3"
   orcid: 0000-0002-7304-021X
 - name: Anders Johansen
   affiliation: 4
   orcid: 0000-0002-5893-6165
 - name: Philippe A. Bourdin
   affiliation: "5, 6"
   orcid: 0000-0002-6793-601X
 - name: Wolfgang Dobler
   affiliation: 7
 - name: Wladimir Lyra
   affiliation: 8
   orcid: 0000-0002-3768-7542
 - name: Matthias Rheinhardt
   affiliation: 9
 - name: Sven Bingert
   affiliation: 10
   orcid: 0000-0001-9547-1582
 - name: Nils Erland L. Haugen
   affiliation: "11, 12, 1"
   orcid: 0000-0002-9184-8722
 - name: Antony Mee
   affiliation: 13
 - name: Frederick Gent
   affiliation: "9, 14"
   orcid: 0000-0002-1331-2260
 - name: Natalia Babkovskaia
   affiliation: 15
 - name: Chao-Chin Yang
   affiliation: 16
   orcid: 0000-0003-2589-5034
 - name: Tobias Heinemann
   affiliation: 17
   orcid: 0000-0002-2476-9733
 - name: Boris Dintrans
   affiliation: 18
 - name: Dhrubaditya Mitra
   affiliation: 1
   orcid: 0000-0003-4861-8152
 - name: Simon Candelaresi
   affiliation: 19
   orcid: 0000-0002-7666-8504
 - name: Jörn Warnecke
   affiliation: 20
   orcid: 0000-0002-9292-4600
 - name: Petri J. Käpylä
   affiliation: 21
   orcid: 0000-0001-9619-0053
 - name: Andreas Schreiber
   affiliation: 15
 - name: Piyali Chatterjee
   affiliation: 22
   orcid: 0000-0002-0181-2495
 - name: Maarit J. Käpylä
   affiliation: "9, 20"
   orcid: 0000-0002-9614-2200
 - name: Xiang-Yu Li
   affiliation: 1
   orcid: 0000-0002-5722-0018 
 - name: Jonas Krüger
   affiliation: "11, 12"
   orcid: 0000-0001-8036-0695
 - name: Jørgen R. Aarnes
   affiliation: 12
   orcid: 0000-0002-5899-2597
 - name: Graeme R. Sarson
   affiliation: 14
   orcid: 0000-0001-6774-9372
 - name: Jeffrey S. Oishi
   affiliation: 23
   orcid: 0000-0001-8531-6570
 - name: Jennifer Schober
   affiliation: 24
   orcid: 0000-0001-7888-6671
 - name: Raphaël Plasson
   affiliation: 25
   orcid: 0000-0003-2319-1463
 - name: Christer Sandin
   affiliation: 1
   orcid: 0000-0002-6370-5505
 - name: Ewa Karchniwy
   affiliation: "12, 26"
   orcid: 0000-0001-6709-1160
 - name: Luiz Felippe S. Rodrigues
   affiliation: "14, 27"
   orcid: 0000-0002-3860-0525
 - name: Alexander Hubbard
   affiliation: 28
 - name: Gustavo Guerrero
   affiliation: 29
   orcid: 0000-0002-2671-8796
 - name: Andrew Snodin
   affiliation: 14
 - name: Illa R. Losada
   affiliation: 1
   orcid: 0000-0002-0416-7516
 - name: Johannes Pekkilä
   affiliation: 9
   orcid: 0000-0002-1974-7150
 - name: Chengeng Qian
   affiliation: 30
   orcid: 0000-0002-5560-5475

affiliations:
 - name: Nordita, KTH Royal Institute of Technology and Stockholm University, Sweden
   index: 1
 - name: Department of Astronomy, Stockholm University, Sweden
   index: 2
 - name: McWilliams Center for Cosmology & Department of Physics, Carnegie Mellon University, PA, USA
   index: 3
 - name: GLOBE Institute, University of Copenhagen, Denmark
   index: 4
 - name: Space Research Institute, Graz, Austria
   index: 5
 - name: Institute of Physics, University of Graz, Graz, Austria
   index: 6
 - name: Bruker, Potsdam, Germany
   index: 7
 - name: New Mexico State University, Department of Astronomy, Las Cruces, NM, USA
   index: 8
 - name: Astroinformatics, Department of Computer Science, Aalto University, Finland
   index: 9
 - name: Gesellschaft für wissenschaftliche Datenverarbeitung mbH Göttingen, Germany
   index: 10
 - name: SINTEF Energy Research, Trondheim, Norway
   index: 11
 - name: Norwegian University of Science and Technology, Norway
   index: 12
 - name: Bank of America Merrill Lynch, London, UK
   index: 13
 - name: School of Mathematics, Statistics and Physics, Newcastle University, UK
   index: 14
 - name: No current affiliation
   index: 15
 - name: University of Nevada, Las Vegas, USA
   index: 16
 - name: Niels Bohr International Academy, Denmark
   index: 17
 - name: CINES, Montpellier, France
   index: 18
 - name: School of Mathematics and Statistics, University of Glasgow, UK
   index: 19
 - name: Max Planck Institute for Solar System Research, Germany
   index: 20
 - name: Institute for Astrophysics, University of Göttinge, Germany
   index: 21
 - name: Indian Institute of Astrophysics, Bengaluru, India
   index: 22
 - name: Department of Physics & Astronomy, Bates College, ME, USA
   index: 23
 - name: Laboratoire d'Astrophysique, EPFL, Sauverny, Switzerland
   index: 24
 - name: Avignon Université, France
   index: 25
 - name: Institute of Thermal Technology, Silesian University of Technology, Poland
   index: 26
 - name: Radboud University, Netherlands
   index: 27
 - name: Department of Astrophysics, American Museum of Natural History, NY, USA
   index: 28
 - name: Physics Department, Universidade Federal de Minas Gerais, Belo Horizonte, Brazil
   index: 29
 - name: State Key Laboratory of Explosion Science and Technology, Beijing Institute of Technology, China
   index: 30

date: 17 September 2020
bibliography: paper.bib
---

# Summary

The Pencil Code is a highly modular physics-oriented simulation code
that can be adapted to a wide range of applications.
It is primarily designed to solve partial differential equations (PDEs)
of compressible hydrodynamics and has lots of add-ons ranging from
astrophysical magnetohydrodynamics (MHD) [@2010ascl.soft10060B] to
meteorological cloud microphysics [@2017JAMES.9.1116L] and engineering
applications in combustion [@2011JCoPh.230.1B].
Nevertheless, the framework is general and can also be applied to
situations not related to hydrodynamics or even PDEs, for example when
just the message passing interface or input/output strategies of the
code are to be used.
The code can also evolve Lagrangian (inertial and noninertial)
particles, their coagulation and condensation, as well as their
interaction with the fluid.
A related module has also been adapted to perform ray tracing
and to solve the eikonal equation.

The code is being used for Cartesian, cylindrical, and spherical geometries,
but further extensions are possible.
One can choose between different time stepping schemes and different
spatial derivative operators.
High-order first and second derivatives are used to deal with weakly 
compressible turbulent flows.
There are also different diffusion operators to allow for both direct numerical
simulations (DNS) and various types of large-eddy simulations (LES).

# High-level functionality

An idea about the range of available modules can be obtained by inspecting
the examples under pencil-code/samples/.
Those are low resolution versions related to applications published in the literature.
Some of the run directories of actual production runs are published through Zenodo.
Below a list of method papers that describe the various applications and tests:

* Coagulation and condensation in turbulence [@2008A&A.486.597J; @2017JAMES.9.1116L],
* Radiative transfer [@2006A&A.448.731H; @2014A&A.571A.68B; @2020GApFD.114.162B],
* Chiral magnetic effect in relativistic plasmas [@2018ApJ.858.124S],
* Primordial gravitational waves [@2020GApFD.114.130R],
* Modeling homochirality at the origin of life [@2004IJAsB.3.209B; @2019OLEB.49.49B],
* Modeling of patterned photochemical systems [@2012ChemEurJ],
* Gaseous combustion and detonation [@2011JCoPh.230.1B; @Zhang_etal_2020comb; @2017CNF.185a160],
* Burning particles, resolved or unresolved [@2020GApFD.114.58Q],
* Flows around immersed solid objects [@2019IJCFD.33.43A; @2020GApFD.114.35A; @2010JFM.661a239],
* Test-field method for turbulent MHD transport [@2010A&A.520A.28R; @2010PhST.142a4028B; @2018A&A.609A.51W],
* Mean-field MHD [@2013SoPh.287.293K; @2013A&A.556A.106J],
* Spherical shell dynamos and convection [@2009ApJ.697.923M; @2020GApFD.114.8K],
* Boris correction for coronal physics [@2020GApFD.114.213C],
* Thermal instability and mixing [@2012ApJ.758.48Y],
* Implicit solver for temperature [@2008A&A.484.29G],
* Dust-gas dynamics with mutual drag interaction [@2007ApJ.662.613Y; @2016ApJS.224.39Y],
* Boundary conditions for the solar atmosphere and HDF5 format [@2020GApFD.114.235B].

# Statement of need and purpose of software

The code is an easily adaptable tool for solving both standard
MHD equations as well as others, such as the test-field equations.
Significant amounts of runtime diagnostics 
as well as Python and IDL libraries for post-processing are available.

Among the currently 83 developers with check-in permission, there are
currently 18 owners who can give others check-in permission.
Of the developers, 35 have done more than 35 commits.
Users have access to the latest development version and can ask to
join the circle of developers by contacting one of the owners.

Every revision on GitHub is verified on 9 tests on travis-ci.com.
The current version is also automatically being tested on 59 hourly
tests and on 79 daily tests.
Continuous progress on the code is driven by the research of
individual developers.

Further developments and interactions between developers and users are
being promoted through annual user meetings since 2004 and a newsletters
since 2020.
Since 2016, a steering committee of five elected owners reviews the
progress and can take decisions of general concern to the Pencil Code
community.

# Ongoing research using the Pencil Code

Current research includes topics from stellar physics, interstellar and intercluster medium, the early universe,
as well as from meteorology and engineering:
small-scale dynamos and reconnection;
primordial magnetic fields and decaying turbulence;
gravitational waves from turbulent sources;
planet formation and inertial particles;
accretion discs and shear flows;
coronal heating and coronal mass ejections;
helical dynamos, helical turbulence, and catastrophic quenching;
helioseismology;
strongly stratified MHD turbulence and negative effective magnetic pressure instability;
convection in Cartesian domains;
global convection and dynamo simulations;
turbulent transport and test-field methods;
hydrodynamic and MHD instabilities and turbulence;
chiral MHD;
turbulent gaseous and solid combustion, particle clustering and deposition on solid walls,
front propagation, radiation & ionization.
As of July 2020, 564 papers have been published that acknowledge use of
the Pencil Code [@zenodo.3466444].

# Key references

The Pencil Code is unique in two ways:
the high level of flexibility and modularity, and the way it is organized
(open source, distributed ownership, openness of development version).

Other software addressing related needs include: 
Athena, CO5BOLD, ENZO, MuRAM, NIRVANA, Stagger, ZEUS, Snoopy, and several other LES codes.
There are also several other engineering DNS codes such as
Sandia-3-Dimensional (S3D), a high-order compressible code,
optimized for combustion, which is not open source, however.
Another example is SpECTRE, a task-based discontinuous Galerkin code for
relativistic astrophysics [@Kidder_2017].
In addition, there are frameworks like Dedalus or Cactus,
which allow one to program the equations in symbolic form.

Some recent research areas that made use of the Pencil Code, as
evidenced by the aforementioned document listing all those papers
[@zenodo.3466444], include:

* Flows around immersed solid objects [@2010JFM.661a239],
* Particle clustering in supersonic and subsonic turbulence [@2019MNRAS.483.5623M; @Karchniwy_etal_2019],
* Cloud microphysics [@2017JAMES.9.1116L],
* Planet and planetesimal formation [@2007Natur.448.1022J; @2007ApJ.670.805O; @2009A&A.497.869L],
* Global simulations of debris disks [@2013Natur.499.184L],
* Stratified shearing box simulations, also with dust [@2011ApJ.740.18O; @2018ApJ.861.47S; @2018ApJ.868.27Y],
* Supernova-driven turbulence [@2013MNRAS.432.1396G],
* Solar dynamo and sunspots [@2005ApJ.625.539B; @2007ApJ.669.1390H],
* Solar corona above active regions [@2011A&A.530A.112B; @2013A&A.555A.123B; @2016PhRvL.116j1101C],
* Fully convective star in a box [@2006ApJ.638.336D],
* Dynamo wave in spherical shell convection [@2012ApJ.755L.22K; @2014ApJ.796L.12W],
* Convection with Kramers opacity law [@2017ApJ.845.23K; @2019A&A.631.122K; @2020GApFD.114.8K],
* MHD turbulence and cascades [@2004PhRvE.70a6308H],
* Turbulent diffusivity quenching with test fields [@2008ApJ.676.740B; @2014ApJ.795.16K].

# Acknowledgements

We acknowledge contributions from all submitters and their supporting
funding agencies.
In particular, we mention the ERC Advanced Grant on Astrophysical Dynamos
(No 227952), the Swedish Research Council,
grants 2012-5797, 2013-03992, 2017-03865, and 2019-04234,
the National Science Foundation under the grant AAG-1615100,
the FRINATEK grant 231444 under the Research Council of Norway, SeRC,
the grant "Bottlenecks for particle growth in turbulent aerosols"
from the Knut and Alice Wallenberg Foundation, Dnr.\ KAW 2014.0048,
the ReSoLVE Centre of Excellence (grant number 307411),
the research project ‘Gaspro’, financed by the Research Council of
Norway (267916),
the European Research Council (ERC) under the European Union's
Horizon 2020 research and innovation programme (Project UniSDyn,
grant agreement n:o 818665), and the Deutsche Forschungsgemeinschaft
(DFG) Heisenberg programme grant KA 4825/2-1.

# References

Pencil Code Implementation Notes
--------------------------------

This directory is contains implementation notes for special
solvers, algorithms, boundary conditions, etc.
These notes tend to be more specialized than what is documented
in the manual.
As always, to have definitive information on what is really done in
the code, it is necessary to check on the actual source code and the
comments therein.

?to remove sample (hyper check Nils...)
##################################
The Pencil Code documentation
##################################

.. admonition:: Welcome!

   This is the new homepage of The Pencil Code documentation!
   
   Explore the page hierarchy below (or in the sidebar), and get started with
   :ref:`contributing your own documentation<Contributing to the documentation>`!

The Pencil Code is primarily designed to deal with weakly compressible turbulent flows,  which is why we use high-order first and second derivatives. To achieve good parallelization, we use explicit (as opposed to compact) finite differences. Typical scientific targets include driven MHD turbulence in a periodic box, convection in a slab with non-periodic upper and lower boundaries, a convective star embedded in a fully nonperiodic box, accretion disc turbulence in the shearing sheet approximation, self-gravity, non-local radiation transfer, dust particle evolution with feedback on the gas, etc. A range of artificial viscosity and diffusion schemes can be invoked to deal with supersonic flows. For direct simulations regular viscosity and diffusion is being used.

Please find `more details on our website <http://pencil-code.nordita.org/>`_.


.. toctree::
   :caption: Introduction
   :maxdepth: 2

   intro/getting_started
   intro/usingrst
   intro/links
   intro/discussion

.. toctree::
   :caption: User manuals
   :maxdepth: 2

   Quick Guide <manuals/quick-guide>
   

.. toctree::
   :caption: Tutorials
   :maxdepth: 2

   tutorials/pencil/tutpencil
   tutorials/python/tutpython
   tutorials/mathematica/tutmathematica


.. toctree::
   :caption: Code documentation
   :maxdepth: 2

   toc/modpython
   toc/modidl
   toc/modfortran

  
 


Discussion groups
======================

The best way to keep up-to-date with the code updates is to join any of this discussion/work groups:

* Pencil code commits: https://groups.google.com/u/1/g/pencil-code-commits
* Pencil code discuss: https://groups.google.com/u/1/g/pencil-code-discuss
* Python for pencil: https://groups.google.com/u/1/g/pencil-code-python
* Documentation for the Pencil Code: https://groups.google.com/g/pencil-code-doc


***************
Getting started
***************

Welcome to The Pencil Code documentation page! To help you get started, here are a few tips for using and contributing to this space:

Using the documentation
=======================

The purpose of this space is to bring together The Pencil Code documentation
in a unified format and organized in a logical, hierarchical structure.
The ultimate goal is to make the process of locating and reading documentation
as easy and enjoyable as possible.

Contributing to the documentation
=================================

In order to contribute to the documentation, you will need to clone the
GitHub repository containing the source of the documentation, edit the
necessary files, and then push your changes back to the repository. The
documentation will then be built and deployed automatically using  the ReadTheDocs platform.

The following instructions are for Ubuntu systems (64-bit, 16.04+):

Cloning the repository
----------------------


Go to a directory and type in a terminal:

.. code:: bash

   git clone git@github.com:pencil-code/pencil-code.git

If you have  a github username 'MY_GITHUB_NAME' and like to submit changes you can use: 

.. code:: bash

   git clone http://MY_GITHUB_NAME@github.com/pencil-code/pencil-code.git
   git config --global credential.helper 'cache --timeout=3600'
   git config --global branch.autosetuprebase always

After cloning the repository, you can access all the documentation files in the directory:

.. code:: bash

   cd pencil-code/doc/readthedocs

How to build locally (fast)
---------------------------

Follow these instructions to build the reStructuredText documentation (i.e., manuals), but not
the auto-generated code documentation. The build is very fast (few seconds).

1. Make sure you have ``sphinx`` installed and that ``sphinx-build`` is in your ``PATH``.

2. Make sure you have the following ``python3`` scripts installed (e.g., with ``pip3``)::

      sphinx-rtd-theme
      sphinxcontrib-images
      sphinx-fortran
      sphinx-git

3. Build:

   .. code:: bash

      # Fast: do not build auto-generated code documentation
      make fast

   The html files will be built into *_build/html*.

How to build locally (slow)
---------------------------

Follow these instructions to build the complete documentation, including
the auto-generated code documentation. The build is slow (several minutes).

.. warning:: 

   Sphinx imports the entire *Pencil* package in order to generate the documentation (right now, only the *Pencil Python* module. More to come!)

1. Make sure you have ``sphinx`` installed and that ``sphinx-build`` is in your ``PATH``.

2. Make sure you have the following ``python3`` scripts installed (e.g., with ``pip3``)::

      sphinx-rtd-theme
      sphinxcontrib-images
      sphinx-js

      astropy
      numpy
      scipy

3. Update your local copy of the repository (in order to have freshly autogenerated documentation) and build:

   .. code:: bash

      git pull --rebase

4. Build:

   .. code:: bash

      # Slow: build auto-generated code documentation
      make html

The html files will be built into *_build/html*.

Tips for the Python documentation (numpy style)
-----------------------------------------------

The *Pencil Python* documentation follows the numpy style docstring convention.

For a thorough example please see `the napoleon extension website <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_.



.. tip::

   To make sure sphinx will be successful in generating the documentation, go to
   the  python directory

   .. code:: bash

      # from the directory containing conf.py
      cd ../../python 
      python 
   
   and try to import
   the ``pencil`` package. If the import succeeds, it is likely
   that sphinx will also succeed.



Tips for the IDL documentation
------------------------------

Not yet available.


Tips for the Fortran documentation
----------------------------------

Not yet available.
The Pencil Code links
======================

These are the essential project links:

* Project homepage: http://pencil-code.nordita.org/
* Code repository: https://github.com/pencil-code/pencil-code
* Wiki: https://github.com/pencil-code/pencil-code/wiki

Using reStructuredText
======================

All documentation available on this page is written in the reStructuredText
(reST) markup language.

About reST
----------

reST is a simple markup language for plain text files. It is used to
semantically mark parts of a document (e.g., section headings, bold text, lists,
source code, etc.) for further processing and uniform rendering.

Sphinx is a documentation generator that, in our case, achieves two goals:

  - processes the reST documents and renders them to HTML and PDF;
  - autogenerates the code documentation for Python, IDL and Fortran projects
    (i.e., it extracts and formats lists of classes, functions, etc.,
    each with descriptions based on comments found in the source code).

While the `reST <https://docutils.sourceforge.io/rst.html>`_
(and `Sphinx <https://www.sphinx-doc.org/en/master/contents.html>`_)
official documentation pages are exhaustive, they are perhaps not recommended
for a beginner, as they necessarily contain a lot of information that is
not relevant for our documentation page.
We suggest starting with https://rest-sphinx-memo.readthedocs.io/en/latest/ReST.html,
which is a quick reference for reST and Sphinx that was specifically created
to cover a small subset of features that are likely to be used on a daily basis.

Converting existing documents
-----------------------------

The utility *pandoc* can be used to convert a variety of formats, including
Microsoft Word (*doc*, *docx*), Libre Office (*odt*), LaTeX,
HTML, XML, and if all else fails even PDF, to reST.

The syntax of the command is:

.. code:: bash

   pandoc input.doc -o output.rst

where *input.doc* is your input document, in any format other than *rst*.

Style guideline
---------------

Headings
~~~~~~~~

In reST, headings are marked by underlining them with the same character:

.. code:: rst

   This is a heading
   =================

In the The Pencil Code, the following markers should be used, in this order:

.. code:: rst

   Title of your page
   ==================

   Section
   -------

   Subsection
   ~~~~~~~~~~

   Sub-subsection
   ++++++++++++++

You should not use further levels of headings, as it would prevent optimal
rendering of the table of contents in the left-hand sidebar. You can structure
your document further by using the ``.. rubric::`` directive.

Do not use ``###`` and ``***``, as they are already used for higher-level headings
(e.g., on the main landing page).

Admonitions
~~~~~~~~~~~

The use of admonition directives can greatly enhance the user experience by
presenting tips, warnings, important notes, etc. in a way that stands out from
the rest of the document.

The following admonition directives are available for the Pencil Code: *attention*,
*caution*, *danger*, *error*, *hint*, *important*, *note*, *tip*, *todo*,
*warning*.

Any of the previous values can be used as follows:

.. code:: rst

   .. note::

      This is a note.

producing the following output:

.. note::

   This is a note.

Keep in mind that overuse of admonitions will detract from the
document flow too much, and consequently worsen the user experience.
**Use them sparingly.**


Images
~~~~~~


Three different directives allow for the addition images in the documentation.
Please, see `this guide <https://docutils.sourceforge.io/docs/ref/rst/directives.html#images>`_ 
for a full description.

#. The simplest one is the ``image`` directive:

   .. code:: rst

      .. image:: pics/myimage.png

   Accepted options for the directive are the  width and alternative text for screen readers:

   .. code:: rst

      .. image:: pics/myimage.png
         :width: 400
         :height: 100px
         :scale: 50 %
         :alt: alternate text
         :align: right
      

#. The ``figure`` directive supports all the options of the ``image`` directive and  allows for adding a caption to the figure:
   
   .. code:: rst

      .. figure:: pics/myimage.png
         :scale: 50 %
         :alt: Flow patterns in the Sun

         This is the caption of the figure (a simple paragraph).

         This is the legend of the figure, which can include a table:

         +-----------------------+-----------------------+
         | Symbol                | Meaning               |
         +=======================+=======================+
         | .. image:: arrow.png   | Magnetic field lines |
         +-----------------------+-----------------------+
         | .. image:: lines.png  | Velocity lines        |
         +-----------------------+-----------------------+
   
   There must be blank lines before the caption paragraph and before the legend. 
   To specify a legend without a caption, use an empty comment (“..”) in place of the caption.
 

#. The ``thumbnail`` directive allows you expand the image by clicking on it:

   .. code:: rst
   
      .. thumbnail:: pics/myimage.png
         :width: 500px


Videos
~~~~~~

You can add short movies to your documentation by using the ``.. video::``
directive. Any video that works inside an HTML5 *video* tag can be used (i.e.,
mp4, webm, ogg). Follow these steps to add your video:

- Add the ``.. video:: <video_url>`` directive in your rst file,
  where you want the video to be rendered.
- It is not necessary to specify any options (height, width, etc.), but if
  you want to have a look at the documentation of the extension:
  https://github.com/sphinx-contrib/video

This is the recommended way of adding videos, since they should not
be committed to the *ingdoc* git repository, but rather stored on a
separate server.

However, if you absolutely need to store the video with the documentation,
follow these steps instead:

- Copy the video file to the directory ``_static``. This is necessary at the
  moment, since we have not found a way (yet) for Sphinx to deploy the file
  otherwise.
- Add the ``.. video:: <relative_path_to_video>`` directive in your rst file,
  where you want the video to be rendered. The path is relative to your rst file,
  so it will probably look similar to ``../_static/video.mp4``.

{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

   .. rubric:: {{ _('Revision history') }}
   .. git_changelog::
      :fullname: {{ fullname }}

   {% if attributes %}
   .. rubric:: {{ _('List of attributes') }}
   .. autosummary::
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}

   {% if exceptions %}
   .. rubric:: {{ _('List of exceptions') }}
   .. autosummary::
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}

   {% if functions %}
   .. rubric:: {{ _('List of functions') }}
   .. autosummary::
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}

   {% if classes %}
   .. rubric:: {{ _('List of classes') }}
   .. autosummary::
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   {% for item in attributes %}
   .. autodata:: {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: {{ _('Exceptions') }}

   {% for item in exceptions %}
   .. autoexception:: {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions %}
   {% if functions %}
   .. rubric:: {{ _('Functions') }}

   {% for item in functions %}
   .. autofunction:: {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: {{ _('Classes') }}

   {% for item in classes %}
   .. autoclass:: {{ item }}
      :special-members: __init__
      :show-inheritance:
      :members:
      :undoc-members:
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% block modules %}
{% if modules %}
.. rubric:: Modules

.. autosummary::
   :toctree:
   :template: custom-module-template.rst
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}*****************
Quick start guide
*****************

Required software
=================

Linux
-----

A Fortran and a C compiler are needed to compile the code. Both
compilers should belong to the same distribution package and version
(e.g. GNU GCC or Intel).

MacOS X
-------

For Mac, you first need to install Xcode from the website
http://developer.apple.com/, where you have to register as a member.
Alternatively, an easy to install ``gfortran`` binary package can be
found at the website http://gcc.gnu.org/wiki/GFortranBinaries. Just
download the archive and use the installer contained therein. It
installs into ‘``/usr/local/gfortran``’ with a symbolic link in
‘``/usr/local/bin/gfortran``’. It might be necessary to add the
following line to the “``.cshrc``”-file in your ‘``/home``’ folder:

::

     setenv PATH /usr/local/bin:\$PATH

Download the Pencil Code
========================

The Pencil Code is an open source code written mainly in Fortran and
available under GPL. General information can be found at our official
homepage:

http://pencil-code.nordita.org/.

The latest version of the code can be downloaded with ``svn``. In the
directory where you want to put the code, type:

::

     svn checkout https://github.com/pencil-code/pencil-code/trunk/ pencil-code

Alternatively, you may also use ``git``:

::

     git clone https://github.com/pencil-code/pencil-code.git

More details on download options can be found here:
http://pencil-code.nordita.org/download.php

The downloaded ‘``pencil-code``’ directory contains several
sub-directories:

#. ‘``doc``’: you may build the latest manual as PDF by issuing the
   command ``make`` inside this directory

#. ‘``samples``’: contains many sample problems

#. ‘``config``’: has all the configuration files

#. ‘``src``’: the actual source code

#. ‘``bin``’ and ‘``lib``’: supplemental scripts

#. ‘``idl``’, ‘``python``’, ‘``julia``’, etc.: data processing for
   diverse languages

Configure the shell environment
===============================

You need to load some environment variables into your shell. Please
change to the freshly downloaded directory:

::

     cd pencil-code

Probably you use a ``sh``-compatible shell (like the Linux default shell
``bash``), there you just type:

::

     . sourceme.sh

(In a ``csh``-compatible shell, like ``tcsh``, use this alternative:
``source sourceme.csh`` )

Your first simulation run
=========================

Create a new run-directory
--------------------------

Now create a run-directory and clone the input and configuration files
from one of the samples that fits you best to get started quickly (here
from ‘``pencil-code/samples/1d-tests/jeans-x``’):

::

     mkdir -p /data/myuser/myrun/src
     cd /data/myuser/myrun
     cp $PENCIL_HOME/samples/1d-tests/jeans-x/*.in ./
     cp $PENCIL_HOME/samples/1d-tests/jeans-x/src/*.local src/

Your run should be put outside of your ‘``/home``’ directory, if you
expect to generate a lot of data and you have a tight storage quota in
your ‘``/home``’.

Linking to the sources
----------------------

One command sets up all needed symbolic links to the original Pencil
Code directory:

::

     pc_setupsrc

Makefile and parameters
-----------------------

Two basic configuration files define a simulation setup:
“``src/Makefile.local``” contains a list of modules that are being used,
and “``src/cparam.local``” defines the grid size and the number of
processors to be used. Take a quick look at these files...

Single-processor
~~~~~~~~~~~~~~~~

An example “``src/Makefile.local``” using the module for only one
processor would look like:

::

     MPICOMM=nompicomm

For most modules there is also a “``no``”-variant which switches that
functionality off.

In “``src/cparam.local``” the number of processors needs to be set to
``1`` accordingly:

::

     integer, parameter :: ncpus=1,nprocx=1,nprocy=1,nprocz=ncpus/(nprocx*nprocy)
     integer, parameter :: nxgrid=128,nygrid=1,nzgrid=128

Multi-processor
~~~~~~~~~~~~~~~

If you like to use MPI for multi-processor simulations, be sure that you
have a MPI library installed and change “``src/Makefile.local``” to use
MPI:

::

     MPICOMM=mpicomm

Change the ``ncpus`` setting in “``src/cparam.local``”. Think about how
you want to distribute the volume on the processors — usually, you
should have 128 grid points in the x-direction to take advantage of the
SIMD processor unit. For compilation, you have to use a configuration
file that includes the “``_MPI``” suffix, see below.

Compiling...
------------

In order to compile the code, you can use a pre-defined configuration
file corresponding to your compiler package. E.g. the default compilers
are ``gfortran`` together with ``gcc`` and the code is being built with
default options (not using MPI) by issuing the command:

::

     pc_build

Alternatively, for multi-processor runs (still using the default GNU-GCC
compilers):

::

     pc_build -f GNU-GCC_MPI

Using a different compiler (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you prefer to use a different compiler package (e.g. with MPI support
or using ``ifort``), you may try:

::

     pc_build -f Intel
     pc_build -f Intel_MPI
     pc_build -f Cray
     pc_build -f Cray_MPI

More pre-defined configurations are found in the directory
“``pencil-code/config/compilers/*.conf``”.

Changing compiler options (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Of course you can also create a configuration file in any subdirectory
of ‘``pencil-code/config/hosts/``’. By default, ``pc_build`` looks for a
config file that is based on your ``host-ID``, which you may see with
the command:

::

     pc_build -i

You may add your modified configuration with the filename
“``host-ID.conf``”, where you can change compiler options according to
the Pencil Code manual. A good host configuration example, that you may
clone and adapt according to your needs, is
“``pencil-code/config/hosts/IWF/host-andromeda-GNU_Linux-Linux.conf``”.

Running...
----------

The initial conditions are set in “``start.in``” and the parameters for
the main simulation run can be found in “``run.in``”. In “``print.in``”
you can choose which quantities are written to the file
“``data/time_series.dat``”.

Be sure you have created an empty ‘``data``’ directory.

::

     mkdir data

It is now time to run the code:

::

     pc_run

If everything worked well, your output should contain the line

::

     start.x has completed successfully

after initializing everything successfully. It would then start running,
printing in the console the quantities specified in “``print.in``”, for
instance,

::

   ---it--------t-------dt------rhom------urms------uxpt-----uypt-----uzpt-----
          0      0.00 4.9E-03 1.000E+00  1.414E+00 2.00E+00 0.00E+00 0.00E+00
         10      0.05 4.9E-03 1.000E+00  1.401E+00 1.98E+00 0.00E+00 0.00E+00
         20      0.10 4.9E-03 1.000E+00  1.361E+00 1.88E+00 0.00E+00 0.00E+00 
         .......

ending with

::

     Simulation finished after        xxxx  time-steps
     .....
     Wall clock time/timestep/meshpoint [microsec] = ...

An empty file called “``COMPLETED``” will appear in your run directory
once the run is finished.

If you work with one of the samples or an identical setup in a new
working directory, you can verify the correctness of the results by
checking against reference data, delivered with each sample:

::

     diff reference.out data/time_series.dat

Welcome to the world of Pencil Code!

Troubleshooting...
------------------

If compiling fails, please try the following — with or without the
optional ``_MPI`` for MPI runs:

::

     pc_build --cleanall
     pc_build -f GNU-GCC_MPI

If some step still fails, you may report to our mailing list:
http://pencil-code.nordita.org/contact.php. In your report, please state
the exact point in this quick start guide that fails for you (including
the full error message) — and be sure you precisely followed all
non-optional instructions from the beginning.

In addition to that, please report your operating system (if not
Linux-based) and the shell you use (if not ``bash``). Also please give
the full output of these commands:

::

     bash
     cd path/to/your/pencil-code/
     source sourceme.sh
     echo $PENCIL_HOME
     ls -la $PENCIL_HOME/bin
     cd samples/1d-tests/jeans-x/
     gcc --version
     gfortran --version
     pc_build --cleanall
     pc_build -d

If you plan to use MPI, please also provide the full output of:

::

     mpicc --version
     mpif90 --version
     mpiexec --version

Data post-processing
====================

IDL visualization (optional,)
-----------------------------------------

GUI-based visualization (recommended for quick inspection)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most simple approach to visualize a Cartesian grid setup is to run
the Pencil Code GUI and to select the files and physical quantities you
want to see:

::

   IDL> .r pc_gui

If you miss some physical quantities, you might want to extend the two
IDL routines ``pc_get_quantity`` and ``pc_check_quantities``. Anything
implemented there will be available in the GUI, too.

Command-line based processing of “big data”
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please check the documentation inside these files:

+------------------------------------------------+--------------------------------------------+
|``pencil-code/idl/read/pc_read_var_raw.pro``    | efficient reading of raw data              |
+------------------------------------------------+--------------------------------------------+
|``pencil-code/idl/read/pc_read_subvol_raw.pro`` | reading of sub-volumes                     |
+------------------------------------------------+--------------------------------------------+
|``pencil-code/idl/read/pc_read_slice_raw.pro``  | reading of any 2D slice from 3D snapshots  |
+------------------------------------------------+--------------------------------------------+
|``pencil-code/idl/pc_get_quantity.pro``         | compute physical quantities out of raw data|
+------------------------------------------------+--------------------------------------------+
|``pencil-code/idl/pc_check_quantities.pro``     | dependency checking of physical quantities |
+------------------------------------------------+--------------------------------------------+


in order to read data efficiently and compute quantities in physical
units.

Command-line based data analysis (may be inefficient)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several idl-procedures have been written (see in ‘``pencil-code/idl``’)
to facilitate inspecting the data that can be found in raw format in
‘``jeans-x/data``’. For example, let us inspect the time series data

.. code:: idl 

   IDL> pc_read_ts, obj=ts

The structure ``ts`` contains several variables that can be inspected by

.. code:: idl

   IDL> help, ts, /structure
   ** Structure <911fa8>, 4 tags, length=320, data length=320, refs=1:
      IT              LONG      Array[20]
      T               FLOAT     Array[20]
      UMAX            FLOAT     Array[20]
      RHOMAX          FLOAT     Array[20]

The diagnostic ``UMAX``, the maximal velocity, is available since it was
set in “``jeans-x/print.in``”. Please check the manual for more
information about the input files.

We plot now the evolution of ``UMAX`` after the initial perturbation
that is defined in “``start.in``”:

.. code:: idl 

   IDL> plot, ts.t, alog(ts.umax)

The complete state of the simulation is saved as snapshot files in
“``jeans-x/data/proc0/VAR*``” every ``dsnap`` time units, as defined in
“``jeans-x/run.in``”. These snapshots, for example “``VAR5``”, can be
loaded with:

.. code:: idl

   IDL> pc_read_var, obj=ff, varfile="VAR5", /trimall

Similarly ``tag_names`` will provide us with the available variables:

.. code:: idl

   IDL> print, tag_names(ff)
   T X Y Z DX DY DZ UU LNRHO POTSELF

The logarithm of the density can be inspected by using a GUI:

.. code:: idl

   IDL> cslice, ff.lnrho

Of course, for scripting one might use any quantity from the ``ff``
structure, like calculating the average density:

.. code:: idl

   IDL> print, mean(exp(ff.lnrho))

Python visualization (optional)
-------------------------------

Be advised that the Python support is still not complete or as
feature-rich as for IDL. Furthermore, we move to Python3 in 2020, and
not all the routines have been updated yet.

Python module requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example we use the modules: ``numpy`` and ``matplotlib``. A
complete list of required module is included in
“``pencil-code/python/pencil/README``”.

Using the ’pencil’ module
~~~~~~~~~~~~~~~~~~~~~~~~~

After sourcing the “``sourceme.sh``” script (see above), you should be
able to import the ``pencil`` module:

::

   import pencil as pc

Some useful functions:

===============================   ======
``pc.read.ts``                    read “``time_series.dat``” file. Parameters are added as members of the class
``pc.read.slices``                read 2D slice files and return two arrays: (nslices,vsize,hsize) and (time)
``pc.visu.animate_interactive``   assemble a 2D animation from a 3D array
===============================   ======


Some examples of postprocessing with Python can be found in the
:ref:` python documentation <modpython>` and in the :ref:` python tutorials <tutpython>`... tutpython:

***********************
Pencil Python Tutorials
***********************

Here you can find some tutorials on how to modify/contribute to the Python Code 
using the Coding style :ref:`pythonstyle` and how to use the code for post-processing :ref:`pythongeneral`.



.. _pythonstyle: 

Python Coding Style
===================

Good coding style greatly improves the readability of the code. Similar
to the guidelines for the Fortran routines, it is strongly recommended
to follow some basic style rules for the python routines. These are some
recommendations extracted from `PEP 008 <https://www.python.org/dev/peps/pep-0008/>`_ and 
`Google Python Style Guide
<https://google-styleguide.googlecode.com/svn/trunk/pyguide.html>`_.


General Style Guide for Python
------------------------------

Indentation and Spaces
~~~~~~~~~~~~~~~~~~~~~~

-  Use 4 spaces per indentation level.
-  Use hanging indent for function calls over multiple lines:

   .. code:: python


        # Aligned with opening delimiter.
        foo = long_function_name(var_one, var_two, var_three, var_four)


-  Surround top-level function and class definitions with two blank lines.

-  Wildcard imports ( from import \* ) should be avoided, as they make
   it unclear which names are present in the namespace, confusing both
   readers and many automated tools.
-  More than one space around an assignment (or other) operator to align
   it with another should be avoided. **No**:

   .. code:: python

      x             = 1
      y             = 2
      long_variable = 3

-  Always surround these binary operators with a single space on either
   side: assignment ( = ), augmented assignment ( += , -= etc.),
   comparisons ( == , < , > , != , <> , <= , >= , in , not in , is , is
   not ), Booleans ( and , or , not ).
-  If operators with different priorities are used, consider adding
   whitespace around the operators with the lowest priority(ies).
   
   **Yes**:

   .. code:: python

      i = i + 1
      submitted += 1
      x = x\*2 - 1

   **No**:

   .. code:: python

      
      i=i+1
      submitted +=1
      x = x * 2 - 1
      
-  Don’t use spaces around the = sign when used to indicate a keyword
   argument or a default parameter value. 
   
   **Yes**:

   .. code:: python

      def complex(real, imag=0.0):
            return magic(r=real, i=imag)
      

   **No**:

   .. code:: python

      def complex(real, imag = 0.0):
            return magic(r = real, i = imag)
     
Comments
~~~~~~~~

-  Comments should be complete sentences.
-  Block comments generally apply to some (or all) code that follows
   them, and are indented to the same level as that code. Each line of a
   block comment starts with a # and a single space (unless it is
   indented text inside the comment). Paragraphs inside a block comment
   are separated by a line containing a single # .

Docstrings
~~~~~~~~~~

Always use docstrings for classes and functions which can be accessed by
the user. 

We are now working with read the docs and sphinx to create automatic documentation for the code, hence we have updated the style guide for creating docstrings.

We are using Numpy docstring style, and require the following fields in the docstring:

- General description of the Class/function
- Signature: how the function can be called
- Parameters: list of parameters of the class/function
- Returns: type of variable the function returns
- Examples: at least one example of usage
- Notes (ptional): any further comments to the function


.. code:: python

   def complex(real=0.0, imag=0.0):
        """
        Form a complex number.

        Signature
        ---------
        complex(8,7)

        Parameters
        ----------
         *real*: float
             the real part (default 0.0)
         *imag*: float
             the imaginary part (default 0.0)

        Returns
        -------
        complex number with real and imaginary part

        Examples 
        --------
        Define two complex numbers:
        >>> a = complex(3,5)
        >>> b = complex(4,7)
        >>> print(a)
        (3+5j)
        >>> a + b
        (7+12j)
        """
  
Naming Convention
~~~~~~~~~~~~~~~~~

module_name, package_name, ClassName, method_name, ExceptionName,
function_name, GLOBAL_CONSTANT_NAME, global_var_name, instance_var_name,
function_parameter_name, local_var_name

Exceptions for >our< code: datadir, varfile, varfiles, …

pylint
~~~~~~

Run pylint over your code. pylint is a tool for finding bugs and style
problems in Python source code. It finds problems that are typically
caught by a compiler for less dynamic languages like C and C++.

Default Function Arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~

Do not use mutable objects as default values in the function or method
definition. 

**Yes**:

.. code:: python

   def foo(a, b=None):
           if b is None:
               b = []

**No**: 

.. code:: python

        def foo(a, b=[]):


Private Methods
~~~~~~~~~~~~~~~

Python does not know any private methods or class member. In order to
somewhat hide such methods use two underscores in the function
definition: ``def __magicAttributes(self, param):``.

Others
~~~~~~

-  Use ``''.startswith()`` and ``''.endswith()`` instead of string
   slicing to check for prefixes or suffixes. startswith() and
   endswith() are cleaner and less error prone. For example: **Yes**:
   ``if foo.startswith('bar'):`` **No**: ``if foo[:3] == 'bar':``
-  For sequences, (strings, lists, tuples), use the fact that empty
   sequences are false. 

   **Yes**:

   .. code:: python
     
      if not seq:
      if seq:
      

   **No**:

   .. code:: python
      
      if len(seq)
      if not len(seq)
      

-  Don’t compare boolean values to True or False using == . 
**Yes**: ``if greeting:`` **No**: ``if greeting == True:``
-  Check if a variable has a particular type by using ``isinstance``,
   e.g.: ``isinstance(my_variable, list)``.


Pencil Code Specific Style
--------------------------

Classes/Objects
~~~~~~~~~~~~~~~

Use classes as much as possible. When you write a function try to embed
it into a class as **init** function which should return the desired
result. This has the advantage of adding methods to the returned object
which can modify the data. Read-methods always give back objects
containing the whole information (container philosophy). Therefore we
use classes if possible.

Data Directory
~~~~~~~~~~~~~~

The default data directory is always ‘./data’ and not ‘data’.

File Headers
~~~~~~~~~~~~

Start each file with the file ID and  a short
description of the routines.
(The authors' list is no longer required since it can be easily accesed through git history.)

.. code:: python

   
   # varfile.py
   #
   # Read VAR files. Based on the read_var.pro IDL script.
   #
   # NB: the f array returned is C-ordered: f[nvar,nz,ny,nx]
   #     NOT Fortran as in Pencil (& IDL):  f[nx,ny,nz,nvar]
   
  

Import Libraries
~~~~~~~~~~~~~~~~

-  Import numpy as *np* instead of *N*.
-  Import pylab as *plt* instead of *P*.

If you need to access libraries in some routines in your module, import
them in the routine, rather than the head of the module. That way they
are not visible by the user.

**Yes**:

.. code:: python

        # my_module.py

   class MyClass(object):
       """
       Some documentation.
       """

       def __init__(self):
           import numpy as np

           self.pi = np.pi

**No**:

.. code:: python

        # my_module.py
        import numpy as np

        class MyClass(object):
        """
        Some documentation.
        """

        def __init__(self):
                self.pi = np.pi</pre>


Further Reading
---------------

`<https://www.python.org/dev/peps/pep-0008/#tabs-or-spaces>`_

`<https://google-styleguide.googlecode.com/svn/trunk/pyguide.html>`_



.. _pythongeneral: 

Pencil Code Commands in General
===============================

For a list of all Pencil Code commands start IPython and type ``pc. <TAB>`` (as with auto completion).
To access the help of any command just type the command followed by a '?' (no spaces), e.g.:

.. code:: 

        pc.math.dot?
        Type:       function
        String Form:<function dot at 0x7f9d96cb0cf8>
        File:       ~/pencil-code/python/pencil/math/vector_multiplication.py
        Definition: pc.math.dot(a, b)
        Docstring:
        take dot product of two pencil-code vectors a & b with shape

        a.shape = (3, mz, my, mx)
        
You can also use ``help(pc.math.dot)`` for a more complete documentation of the command.

There are various reading routines for the Pencil Code data. All of them return an object with the data. To store the data into a user defined variable type e.g.

.. code:: python

        ts = pc.read.ts()

Most commands take some arguments. For most of them there is a default value, e.g.

.. code:: python

        pc.read.ts(file_name='time_series.dat', datadir='data')

You can change the values by simply typing e.g.


.. code:: python

        pc.read.ts(datadir='other_run/data')


Reading and Plotting Time Series
================================

Reading the time series file is very easy. Simply type

.. code:: python

        ts = pc.read.ts()

and python stores the data in the variable ``ts``. 
The physical quantities are members of the object ``ts`` and can be accessed accordingly, e.g. ``ts.t, ts.emag``. 
To check which other variables are stored simply do the tab auto completion ``ts. <TAB>``.

 Plot the data with the matplotlib commands:

.. code:: python

        plt.plot(ts.t, ts.emag)


The standard plots are not perfect and need a little polishing. See further down about making pretty plots.
You can save the plot into a file using the GUI or with

.. code:: python

        plt.savefig('plot.eps')

Reading and Plotting VAR files and slice files
==============================================

Read var files:

.. code:: python

        var = pc.read.var()

Read slice files:

.. code:: python

        slices = pc.read.slices(field='bb1', extension='xy')

This returns an object ``slices`` with members ``t`` and ``xy``. 
The last contains the additional member ``xy``.


If you want to plot e.g. the x-component of the magnetic field at the central plane simply type:

.. code:: python
        
        plt.imshow(var.bb[0, 128, :, :].T, origin='lower', extent=[-4, 4, -4, 4], interpolation='nearest', cmap='hot')

For a complete list of arguments of ``plt.imshow`` refer to its documentation.

For a more interactive function plot use:

.. code:: python

        pc.visu.animate_interactive(slices.xy.bb, slices.t)

.. warning::

        arrays from the reading routines are ordered ``f[nvar, mz, my, mx]``, i.e. reversed to IDL. 
        This affects reading var files and slice files.

Create a custom VAR0 or var.dat
===============================

With the functionality of writing snapshots directly into ``VAR*`` or ``var.dat`` the user can now generate an initial condition directly from a numpy array or modify the last snapshot and continue running. The function to be used is in ``python/pencil/io/snapshot.py`` and is called ``write_snapshot``. Here we outline how to generate an initial condition. For modifying the ``var.dat`` only the last steps are necessary.

First we need an empty run. For this let us use ``samples/kin-dynamo``


.. code:: python

        cd pencil-code/samples/kin-dynamo
        pc_setupsrc

In principle we can use any initial condition, as we are going to over write it. But it is cleaner to use

.. code::

        INITIAL_CONDITION = noinitial_condition

in ``src/Makefile.local``. Compile and start:

.. code:: bash

        make
        pc_start

This generates a ``VAR0`` and ``var.dat`` in every proc directory.

Our snapshot writing routine needs to know the cpu structure. Furthermore, we need to know the indices of the primary variables. The first can be obtained from ``src/cparam.local``, while the latter can be read from the newly generated ``data/index.pro``. The numpy arrays that are written need to have the shape [nvar, nz, ny, nz] with the correct order of variables and no ghost zones. Optionally, the number of ghost zones, which is usually 3, can be specified.

Putting it all together our python routine would look something like this:

.. code:: python

        import numpy as np
        import pencil as pc

        # Read the data to obtain the shape of the arrays, rather than the actual data.
        var = pc.read.var(trimall=True)

        # Modify the data.
        var.aa += np.random.random(var.aa.shape)

        # Write the new VAR0 and var.dat files.
        pc.io.write_snapshot(var.aa, file_name='VAR0', nprocx=1, nprocy=1, nprocz=1)
        pc.io.write_snapshot(var.aa, file_name='var.dat', nprocx=1, nprocy=1, nprocz=1)


Examples
========

Standard plots with any plotting library are not the prettiest ones. The same is true for matplotlib. Here are a few pretty examples of plots where the default style is changed. You can add your commands into a script e.g. ``plot_results.py`` and execute it in IPython with ``execfile('plot_results.py')``.

Simple plot:

.. code:: python

        import pencil as pc
        import numpy as np
        import pylab as plt

        # Read the time_series.dat.
        ts = pc.read.ts()

        # Prepare the plot.
        # Set the size and margins.
        width = 8
        height = 6
        plt.rc("figure.subplot", left=0.2)
        plt.rc("figure.subplot", right=0.95)
        plt.rc("figure.subplot", bottom=0.15)
        plt.rc("figure.subplot", top=0.90)
        figure = plt.figure(figsize=(width, height))
        axes = plt.subplot(111)

        # Make the actual plot.
        plt.semilogy(ts.t, ts.brms/ts.brms[0], linestyle='-', linewidth=2, color='black', label=r'$\langle\bar{B}\rangle/\langle\bar{B}\rangle(0)$')
        plt.semilogy(ts.t, ts.jrms/ts.jrms[0], linestyle='--', linewidth=2, color='blue', label=r'$\langle\bar{J}\rangle/\langle\bar{J}\rangle(0)$')
        plt.semilogy(ts.t, ts.jmax/ts.jmax[0], linestyle=':', linewidth=2, color='red', label=r'$J_{\rm max}/J_{\rm max}(0)$')

        plt.xlabel(r'$t$', fontsize=25)
        plt.ylabel(r'$\langle\bar{B}\rangle, \langle\bar{J}\rangle, J_{\rm max}$', fontsize=25)
        plt.title('various quantities', fontsize=25, family='serif')

        # Prepare the legend.
        plt.legend(loc=1, shadow=False, fancybox=False, numpoints=1)
        leg = plt.gca().get_legend()
        # Change the font size of the legend.
        ltext = leg.get_texts() # all the text.Text instance in the legend
        for k in range(len(ltext)):
                legLine = ltext[k]
                legLine.set_fontsize(25)
        frame = leg.get_frame()
        frame.set_facecolor('1.0')
        leg.draw_frame(False)

        # Make plot pretty.
        plt.xticks(fontsize=20, family='serif')
        plt.yticks(fontsize=20, family='serif')
        axes.tick_params(axis='both', which='major', length=8)
        axes.tick_params(axis='both', which='minor', length=4)

        # Create an offset between the xylabels and the axes.
        for label in axes.xaxis.get_ticklabels():
                label.set_position((0, -0.03))
        for label in axes.yaxis.get_ticklabels():
                label.set_position((-0.03, 0))


Simple 2d plot:

.. code:: python

        import pencil as pc
        import numpy as np
        import pylab as plt

        # Read the slices.
        slices = pc.read.slices(field='bb1', extension='xy')

        # Read the grid size.
        grid = pc.read.grid()
        x0 = grid.x[3]
        x1 = grid.x[-4]
        y0 = grid.y[3]
        y1 = grid.y[-4]

        # Prepare the plot.
        # Set the size and margins.
        width = 8
        height = 6
        plt.rc("figure.subplot", left=0.15)
        plt.rc("figure.subplot", right=0.95)
        plt.rc("figure.subplot", bottom=0.15)
        plt.rc("figure.subplot", top=0.95)
        figure = plt.figure(figsize=(width, height))
        axes = plt.subplot(111)

        # Make the actual plot.
        plt.imshow(zip(*slices.xy.bb1[0, :, :]), origin='lower', interpolation='nearest', cmap='hot', extent=[x0, x1, y0, y1])
        plt.xlabel(r'$x$', fontsize=25)
        plt.ylabel(r'$y$', fontsize=25)

        # Set the colorbar.
        cb = plt.colorbar()
        cb.set_label(r'$B_{x}(x,y,z=0)$', fontsize=25)
        cbytick_obj = plt.getp(cb.ax.axes, 'yticklabels')
        plt.setp(cbytick_obj, fontsize=15, family='serif')

        # Make plot pretty.
        plt.xticks(fontsize=20, family='serif')
        plt.yticks(fontsize=20, family='serif')
        axes.tick_params(axis='both', which='major', length=8)
        axes.tick_params(axis='both', which='minor', length=4)

        # Create an offset between the xylabels and the axes.
        for label in axes.xaxis.get_ticklabels():
                label.set_position((0, -0.03))
        for label in axes.yaxis.get_ticklabels():
                label.set_position((-0.03, 0))


IDL to Python guide
===================

A large array of idl scripts have been developed over the years, and many of them served their purpose at the time, but there are many others
of general purpose. Below is a small selection of examples of idl call sequences along with their python counterparts.

Here are the links to a few potentially useful sites:

1. `IDL to Python bridge <https://www.l3harrisgeospatial.com/docs/IDLToPython.html>`_

2. `IDL commands in numerical Python <http://mathesaurus.sourceforge.net/idl-python-xref.pdf>`_

===============================   ======
IDL                               Python
===============================   ======
pc_read_var,obj=var,/trimall      var = pc.read.var(var_file = 'var.dat', trimall = True, sim = SIM)    
help,var                          help(var)       
pc_read_param,obj=param           pc.read.param()
===============================   ======
.. tutmathematica:

****************************
Pencil Mathematica Tutorials
****************************

Here you can find some tutorials on using Mathematica for post-processing.


Loading the package
===================

There are two ways of telling Mathematica where the package is:

1. Modifying ``init.m`` so that the path to the package is automatically added to the ``$Path`` variable in Mathematica.
First, type

.. code::

  FileNameJoin[{$UserBaseDirectory, "Kernel", "init.m"}]

in Mathematica to locate this ``init.m`` file.
Then, add line

.. code::

  AppendTo[$Path, "your/pencil/home/mathematica"]

in this file and save it. In general the path will be ``$PENCIL_HOME/mathematica/``, but of course you may put it somewhere else.
Mathematica will not search in subdirectories, so make sure the package in right in the folder.

After updating ``init.m``, restart the Mathematica kernel (``Evaluation`` -> ``Quit Kernel``).
To use the package, call ``Needs["pc`"]`` in a notebook or a script.

2. Alternatively, if you don't want to modify ``init.m``, you may also call

.. code::

 Needs["pc`","path/to/this/package"]

each time.

.. admonition:: Note:

        To run the package on subkernels you may need to do something like:

        .. code::

          LaunchKernels[];
          AppendTo[$Path, "your/pencil/home/mathematica"]//ParallelEvaluate;
          Needs["pc`"]
          ParallelNeeds["pc`"]

        Then you can do things like ``ParallelTable[readTS[...],...]``.
        Note that both ``Needs`` on the master kernel and ``ParallelNeeds`` on subkernels are needed.
        See also the discussions `here <https://mathematica.stackexchange.com/questions/11595/package-found-with-needs-but-not-with-parallelneeds>`_, and the 'Possible issues' section
        `here <https://reference.wolfram.com/language/ref/ParallelNeeds.html>`_.



Pencil Code Commands in General
===============================

For a list of all Pencil Code commands, load the package and type ``pcFunction[]``.
To access the help of any command just type '?' followed by the command, e.g. ``?readTS``.
You can also check the full definition of the command by typing '??' followed by the command.


Reading and Plotting Time Series
================================

To read the time series, type

.. code::

  data = readTS[sim,var1,var2,...]

where ``var1``, ``var2`` etc. are entries in the time series.
The return of the right side is a ``List`` object, with its elements corresponding to ``var1``, ``var2`` etc.
You can then access, for example, the time series of ``var2`` through ``data[[2]]`` (indexing in Mathematica starts from ``1``).

Alternatively, you may also put a ``List`` object on the left side so that ``var1``, ``var2`` etc. will be assigned to each of its elements.
For example,

.. code ::

  {t,urms} = readTS[sim,"t","urms"]

Make sure that ``Length`` of the left side is equal to the number of ``var``; otherwise Mathematica will complain.

To plot the data, you can say

.. code ::

  fig = ListPlot[Transpose[{t,urms}],Joined->True]

or, in a one-line command,

.. code ::

  (** same as ListPlot[Transpose[readTS[sim,"t","urms"]]] **)
  fig = readTS[sim,"t","urms"]//Transpose//ListPlot

A few options for some internal plotting functions have been reset by the package.
For details check ``??pcLabelStyle`` and ``??pcPlotStyle``.

To export the figure in ``.eps`` format,

.. code ::

  Export["directory/to/export/figure.eps",fig]


Reading VAR files
================================

VAR files can be read using

.. code ::

  data = readVARN[sim,iVAR]

Here ``iVAR`` is the index of the VAR file and starts from 0.

By default, ghost zones will not be trimmed.
You can do it using the option ``"ltrim"->True``, or ``"ltrim"->More``;
the latter will trim ``2*nghost`` cells on each boundary.

To compute "magic" variables, you can use

.. code ::

  data = readVARN[sim,iVAR,{"oo","bb","jj"}]

Here "oo", "bb", "jj" refer to vorticity, magnetic, and current fields, respectively.

The return of ``readVARN`` is an ``Association`` object (i.e., ``Head[data]=Association``).
You can obtain all of its keys by ``Keys[data]``. Here is an example of its return:

.. code ::

  data = readVARN[sim,iVAR,{"oo","bb","jj"}];
  Keys[data]
  (* {"t", "dx", "dy", "dz", "deltay", "lx", "ly", "lz", "x", "y", "z",
     "uu1", "uu2", "uu3", "lnrho", "ooo1", "ooo2", "ooo3", "bbb1", "bbb2",
     "bbb3", "jjj1", "jjj2", "jjj3"} *)

Magic variables are named using triple characters, to avoid shadowing the auxilliary ones
written by the code (which will be "oo1" etc.).

The ``x`` coordinates of the mesh points is then ``data["x"]``, which will have length
``(16+6)^3`` if the resolutoin is ``16^3`` and ``nghost=3``.
One can form a three-dimensional map of ``uu1`` using

.. code ::

  uu1 = Transpose[ data/@{"x","y","z","uu1"} ];
  (* {{x1,y1,z1,f1},{x2,y2,z2,f2},...} *)

Sometimes the following method is also useful:

.. code ::

  Clear[uu1]
  grid = Transpose[ data/@{"x","y","z"} ];
  uu1 = Association[ Thread[ grid->data["uu1"] ] ];

Then ``uu1`` becomes a "function" and its value at ``{x1,y1,z1}`` is simply ``uu1[{x1,y1,z1}]``.

Visualizing slices from VAR files
================================

A quick way to make a density plot from ``data`` is

.. code ::

  showSlice[data, "uu1", {"z", 8}]

Here ``{"z",8}`` instructs to plot the 8th slice in the ``z`` direction.

For vector fields one can also use

.. code ::

  showSliceVector[data, "uu", {"z", 8}]

Notice the second argument is just ``"uu"`` with no index.
The function then makes a density plot of the out-of-plane component of (here ``"uu3"``),
and a superposed vector plot of the in-plane components (here ``"uu1"`` and ``"uu2"``).

Reading video files
================================

To read video or slice files, one uses

.. code ::

  {slices,times,position}=readSlice[sim,"uu1","xy2"]

The returned ``slices`` variable is a ``List`` of all slices at different times, and can
be visualized by, say, ``DensityPlot[ slices[[1]] ]``.
``position`` tells you the spatial coordinate of the slices.

Here is an example to make a video:

.. code ::

  Clear[makeFrame]
  makeFrame[ slice_,time_ ] := DensityPlot[ slice, PlotLabel->"t="<>ToString@time]
  frames = MapThread[ makeFrame, {slices,times} ];
  (* to view the video in the notebook; can be slow if too many frames*)
  ListAnimate[ frame, AnimationRunning->False ]
  (* output to a movie file *)
  Export[ "your/output/directory/video.mov", frames, FrameRate->24 ]

One can also visualize variables in a 3D box.
For more information see the comments of ``makeBox`` and ``makeBoxes``.














.. tutpencil:

***********************
Pencil Code Tutorials
***********************

The Pencil Code is written in Fortran90, and is hosted at github 
`<https://github.com/pencil-code/pencil-code>`_.



Modifying the Pencil Code
=========================

.. note::

        Adapted from this `github wiki page <https://github.com/pencil-code/pencil-code/wiki/>`_

Commit rights for the Pencil Code are given out quite liberally, but
they come with responsibility: - only commit code that is meaningful and
necessary - avoid committing code that breaks any of the auto tests -
discuss major changes with other developers first (e.g. at the
pencil-code-discuss mailing list) - follow the existing :ref:`pencilstyleguide`.

When developing the Pencil Code, we often communicate via commit
messages. A typical example is:

.. code::

        Fixed a bug in hydro.do_something().
        @Paul, you wrote the original subroutine, can you check whether my
        fix is OK?

and ideally a few commits down the line, we would have

.. code::

        Improved Peter's fix of the hydro.do_something() subroutine.
        Thanks for finding and analyzing this.

For this mode of communication to work, we **must** be able to rely on
our co-developers to read the commit messages. If they do not, their
contributions have a strongly reduced value for the code, and they
probably should not have commit rights in the first place.

By far the easiest way of reading the commit messages is to subscribe to
``pencil-code-commits@googlegroups.com`` and read at least superficially
through the commit messages as they appear in the mailbox. Individual
developers may prefer other ways of keeping up to date, but you should
be aware that if you do not follow what is going on with the code, it is
likely that parts of the code that you contributed will eventually get
moved around, altered, or even removed.

If you want to become a member of any of the groups without being a
committer, and if your name is not already known to us (or by googling)
we would appreciate a brief email to us explaining your interest. This
is to prevent spammers entering our lists. The pencil-code-core list, on
the other hand, is reserved for project owners only.

.. _pencilstyleguide:

Coding Style guide
==================

.. note::

        Adapted from this `github wiki page <https://github.com/pencil-code/pencil-code/wiki/CodingStyle>`_

We describe here the Pencil Code coding style and best practice to write
code in form of a checklist that should applied to each of the given
items. Of course, no rules can be hammered in stone and always need some
reality check to ensure these rules help improving the code and not the
opposite.

Is a module…
------------

-  named well in the sense that its name describes its purpose?
-  abstract enough to form a separate module?
-  consistent within the existing modules scheme?
-  interface obviously revealing how the module should be used?
-  interface abstract enough so that the module can be used without
   thinking about how the services are implemented?


Is a subroutine…
----------------

-  name revealing exactly what this subroutine is doing?
-  implementing only one task that is well defined?
-  containing only code parts that would not better be put in a separate
   subroutine?
-  interface obviously revealing how the subroutine should be used?
-  interface abstract enough so that the subroutine can be used without
   thinking about how it is implemented in detail?


Is a data type…
---------------

-  named well so that its name describes its content type?
-  descriptive so that it helps to document its variable declaration?
-  simple so that it minimizes complexity?
-  that needs to be complex operated only through access subroutines
   (set_XY/get_XY)?

Is a variable…
--------------

-  really necessary to be variable (and not a constant)?
-  given all possible attributes (like “save”, “parameter”,
   “intent(in/out)”)?
-  not some kind of a “magic number” or “magic string” that should be
   converted to a named constant?
-  not redundant and distinct to all other variables that are otherwise
   available?
-  used only for the single purpose that it was intended for by its
   name?
-  additionally defined and used, if the clarity of the code is
   significantly improved?

Is a variable name…
-------------------

-  chosen well so that it describes its content (i.p. not its data
   type)?
-  readable and following the style
   “lower_case_words_connected_with_underscores”?
-  of a boolean indicating its flag-like behavior (currently “lflag” for
   a flag)?
-  of a loop counter more informative than just i, j, k, l, m, n? ([i,j]
   should be reserved only for general matrix and vector operations,
   while [l,m,n] are reserved for the m-n-loop, for handling the f-array
   or p-pencil indices).

Is a logical statement…
-----------------------

-  using only simple boolean expressions?
-  better stored it into an additional boolean variable or put into a
   boolean function, if it is to be reused?
-  not using double negations? (oops!)

Is an ``if``-``else``-construct…
--------------------------------

-  consisting of small blocks where both blocks are of similar size?
-  written so that the “normal” case appears first?
-  used to minimize complexity?

Is a loop…
----------

-  performing exactly one well-defined function, as a subroutine would?
-  implemented using the best matching type:
   ``do``/``while``/``repeat``?
-  nested in another loop only if necessary?

Is the code…
------------

-  representing its own logical structure?
-  nominal path or calling sequence clear and easy to follow?
-  organized so that related statements are grouped together?
-  free of relatively independent code blocks that could stay in
   subroutines?
-  hiding implementation details to the greatest extent?
-  written in respect to the problem solution and not in terms of
   programming-language requirements?
-  initializing all variables outside any conditional statements?
   (e.g. code inside ``if``-``elseif``-``elseif``-constructs might never
   be executed!)
-  compiling without any compiler warnings? (yes, they do have a serious
   background even if it is sometimes not obvious!)

Is a commit message…
--------------------

-  not “new settings” or “some corrections” or “minor changes” or
   similar?
-  telling which code block is affected?
-  containing all relevant major changes?
-  informative about the consequences of previous bugs that are now
   fixed?
-  giving hints what to search or where to start reading if anyone is
   interested in details?

FORTRAN formatting
------------------

-  use the Fortran95 language standard (F95)
-  use spaces for indentation, two spaces represent one indentation
   level
-  use spaces for formatting of tabular data or comments
-  no spaces at the end of a line
-  one empty line is enough to split code blocks
-  two empty lines can be used to split distinct parts of code
-  in-code comments should be indented together with the code
-  block-like comments (e.g. function headers) start at the beginning of
   a line
-  use spaces around operators, where applicable
-  line-breaks are required after 130 characters (F95)
-  line-breaks can be used to significantly improve readability of long
   code lines

Typical rule-breaker and its solution
-------------------------------------

-  ``goto`` => implement a loop or an ``if``-``else``-construct
-  ``entry`` => implement an interface or split into distinct
   subroutines
-  ``format`` => put the format string inside each ``write`` statement
-  hard-coded file units => use a named constant
-  hard-coded string length => use pre-defined global constants

Recommended further reading
---------------------------

-  Kernighan, Brian, and Plauger: “The Elements of Programming Style”,
   2nd ed., McGraw-Hill, New York, 1978
-  Kernighan, Brian, and Pike: “The Practice of Programming”, Addison
   Wesley, Reading (Massachusetts), 1999
-  McConnell: “Code Complete”, 2nd ed., Microsoft Press, Redmont
   (Washington), 2004
-  Hunt and Thomas: “The Pragmatic Programmer”, Addison Wesley, Reading
   (Massachusetts), 1999


Social Rules?
=============

.. note::

        Adapted from this `github wiki page <https://github.com/pencil-code/pencil-code/wiki/SocialRules>`_

Hi guys,

I discussed with a co-developer of a code that is developed pretty much
like Pencil, open-source, a team, version control, etc. Talking to him
about code development, I asked if there were cases of flame fights or
heated arguments in the code community. He mentioned a couple of cases,
and pointed me to **books** on open source development where such stuff
is discussed. Not surprisingly, it is quite a common occurrence.

Chapter 6 of the first link, from 102 on (“Difficult People”), is
particularly relevant.

`<http://producingoss.com/>`_ 

`<http://artofcommunityonline.org/>`_ 

Wlad.

Good electronic communication
-----------------------------

For an electronic discussion, there is no such thing as a meta-level of
information transfer. Therefore, every good electronic communicator just
stays with the facts. And *if* an interpretation needs to be done, one
chooses the interpretation that assumes *best* motives of your opponent.
Only then, one has a chance to understand the opponent right. And
without understanding an opponent *fully*, one has no right to answer.
(Philippe)


Code of Conduct
---------------

Although spaces may feel informal at times, we want to remind ourselves
that this is a professional space. As such, the Pencil Code community
adheres to a code of conduct adapted from the Contributor Covenant
(`<https://www.contributor-covenant.org/>`_) code of conduct. All
contributors will be required to confirm they have read our 
`code of conduct <https://github.com/pencil-code/pencil-code/blob/master/license/CODE_OF_CONDUCT.md>`_,
and are expected to adhere to it in all Pencil Code spaces and
associated interactions.

**************
Fortran module
**************

Testing...


Module ascalar
--------------

.. f:automodule:: ascalar


Module Geometrical Types
-------------------------

.. f:automodule:: geometrical_types*************************
IDL code documentation
*************************


Soon!.. _modpython:

********************************
Pencil python code documentation
********************************



:mod:`pencil`: Pencil package
-----------------------------

.. .. autosummary::
   :toctree: ../code/sourcePython
   :template: custom-module-template.rst
   :recursive:

   pencil
