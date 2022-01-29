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
