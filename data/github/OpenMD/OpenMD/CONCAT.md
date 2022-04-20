# What is OpenMD?

[![build](https://github.com/OpenMD/OpenMD/workflows/build/badge.svg)](https://github.com/OpenMD/OpenMD/actions?query=workflow%3Abuild)

OpenMD is an open source molecular dynamics engine which is capable of
efficiently simulating liquids, proteins, nanoparticles, interfaces,
and other complex systems using atom types with orientational degrees
of freedom (e.g. "sticky" atoms, point dipoles, and coarse-grained
assemblies). Proteins, zeolites, lipids, transition metals (bulk, flat
interfaces, and nanoparticles) have all been simulated using force
fields included with the code. OpenMD works on parallel computers
using the Message Passing Interface (MPI), and comes with a number of
analysis and utility programs that are easy to use and modify. An
OpenMD simulation is specified using a very simple meta-data language
that is easy to learn.

## Getting Started

Simulations are started in OpenMD using a single Molecular Dynamics (.omd)
file. These files must start with the `<OpenMD>` tag and must have two
sections:

  1) a C-based `<MetaData>` section, and

  2) a `<Snapshot>` block for initial coordinate and velocity information.

Detailed descriptions of the structures of these files are available
in the `doc` directory. Sample simulations are available in the
`samples` directory.

## Requirements

 1) A good C++11-compliant compiler. We've built and tested OpenMD on the
    following architecture & compiler combinations:

| Architecture                   |  CXX | Notes                                     |
|--------------------------------|:----:|-------------------------------------------|
| macOS 10.15 (intel)            | clang| (Xcode 12, Open MPI 3.1)                  |
| macOS 10.15 (intel)            | g++  | (Xcode 12, OpenMPI 3.1)                   |
| Linux (Ubuntu 20.04 - x86\_64) | clang| (clang version 7.0.0 Open MPI 3.1)        |
| Linux (Ubuntu 20.04 - x86\_64) | g++  | (GNU version 9.3.0, Open MPI 3.1)         |
| Linux (RHEL 7.6 - x86\_64)     | icpc | (Intel version 18.0.5, Open MPI 3.0.0)    |

  OpenMD uses features in the C++ standard library. Most (but
  not all) C++ compilers support these features.

 2) CMake, a cross-platform build system which is available at
    [cmake.org](http://www.cmake.org). Most Linux and some Unix distributions
    provide CMake as a standard package. If not, please download it,
    and make sure you get a recent version. Mac OS X users can either
    download the CMake installer or install it from the command line
    using macports.

 3) An implementation of MPI-2 is optional for the single processor
    version of OpenMD, but is required if you want OpenMD to run in
    parallel. We like OpenMPI. Other implementations of MPI-2 might
    work, but we haven't tried. You can get Open MPI here:
    [open-mpi.org](http://www.open-mpi.org/)

 4) Other optional (but quite useful) libraries that will unlock some
    features of OpenMD:

      + Open Babel:  [openbabel.org](http://openbabel.org)
      + Qhull:       [www.qhull.org](http://www.qhull.org)
      + FFTW:        [www.fftw.org](http://www.fftw.org)
      + Doxygen:     [www.doxygen.org](http://www.doxygen.org)

 5) Some of the utility scripts depend on Python, NumPy, SciPy, and Perl.  
    These interpreters are common installations on most flavors of Unix and
    Mac OS X.

## Instructions

 1) Get, build, and test the required pieces above.
 2) mkdir build
 3) cd build
 4) cmake ..
 5) make
 6) umask 0022; sudo make install

That's it.
# Water models in OpenMD

OpenMD supports a large number of water models, only some of which
have examples shown here.

### Single site models
+ SSD, SSD-E, SSD-RF
+ SSDQ, SSDQO

### Three site models
+ SPC, SPC-HW
+ SPC/E
+ TIP3P

### Four site models
+ TIP4P, TIP4P-Ew, TIP4P/2005, TIP4P/Ice

### Five site models
+ TIP5P, TIP5P-E

### Six site models
+ NE6



One of the powerful functionalities of OpenMD is the separation
between declaration of *components* and *stuntDoubles*. In every
(.omd) file, a component block is used to specifiy what types of atoms
or molecules the stuntDoubles correspond to.

```
	component{
		type = "SPCE";
		nMol = 256;
	}

```

```
    <StuntDoubles>
         0    pvqj         -13.882269          -6.595441         -10.074898  8.730000e-04  7.358000e-03  2.248000e-03  7.861640e-01 -3.843370e-01 -3.681500e-01  3.141600e-01 -1.691600e-02  1.558000e-03  1.901000e-02
         1    pvqj         -10.800233          -3.862809          -9.917968 -1.088300e-02  6.540000e-04  4.554000e-03  4.137720e-01 -7.080340e-01 -4.812670e-01 -3.096160e-01 -7.502000e-03  1.201400e-02  9.878000e-03
         2    pvqj         -13.359152         -13.264782          -5.753301  3.626000e-03  3.221000e-03  1.990000e-04  5.155220e-01 -1.470710e-01 -5.520980e-01  6.385880e-01 -3.590000e-03 -1.001200e-02 -2.305800e-02
 ```

Due to this separation, a simulation of water with one model can be
easily transformed into a simulation with a different water model
(regardless of the number of potential sites) by simply changing the
component block.

```
	component{
		type = "TIP4P";
		nMol = 256;
	}

```

Now, the same stuntDoubles will be treated as *TIP4P* water instead of
*SPC/E* water. Of course, due to the change in potential you may have
to re-equilibrate the system. However, no further changes need to be
made.

In addition to the ease of transferring between potentials, OpenMD
also has a water system builder which generates OpenMD (.omd) files
for you, **waterBoxer**.

## waterBoxer

**waterBoxer** is a perl script which generates FCC lattices of water
molecules at user-defined densities and system sizes. The user can
specifiy which water model to generate a system of, or as described
above, change the component block definition from the default once
generated.  **waterBoxer** prints a helpful discription of how its use
and functionalities when passed *-h*.


# References

| Water Model| Number of Sites | DOI  |
| ----------:|:---------------:|-----:|
| SSD        | 1 |[10.1016/S0009-2614(03)01044-3](https://doi.org/10.1016/S0009-2614(03)01044-3) |
| SSD/RF     | 1 |[10.1063/1.1697381](https://doi.org/10.1063/1.1697381)     |
| SSD/E      | 1 |[10.1063/1.1697381](https://doi.org/10.1063/1.1697381)     |
| SPC        | 3 |[10.1021/j100308a038](https://doi.org/10.1021/j100308a038) |
| SPC-HW     | 3 |[10.1063/1.1359183](https://doi.org/10.1063/1.1359183)     |
| SPC/E      | 3 |[10.1021/j100308a038](https://doi.org/10.1021/j100308a038) |
| TIP3P      | 3 |[10.1063/1.445869](https://doi.org/10.1063/1.445869)       |
| TIP4P      | 4 |[10.1063/1.445869](https://doi.org/10.1063/1.445869)       |
| TIP4P-Ew   | 4 |[10.1063/1.1683075](https://doi.org/10.1063/1.1683075)     |
| TIP4P/2005 | 4 |[10.1063/1.2121687](https://doi.org/10.1063/1.2121687)     |
| TIP4P/Ice  | 4 |[10.1063/1.1931662](https://doi.org/10.1063/1.1931662)     |
| TIP5P      | 5 |[10.1063/1.481505](https://doi.org/10.1063/1.481505)       |
| TIP5P-E    | 5 |[10.1063/1.1652434](https://doi.org/10.1063/1.1652434)     |
| NE6        | 6 |[10.1063/1.1562610](https://doi.org/10.1063/1.1562610)     |
# Structures for simulating Ice-Ih

This directory contains unit cell and larger structures for
simulations of ice crystals. Unless otherwise stated, these crystals
are oriented such that the basal face is normal to the z-axis.

## Proton-disordered ice crystals

The following proton-disordered ice-Ih crystals,

+ 3x3x2-C2.omd
+ 3x3x2-e.omd
+ 3x3x2.omd
+ 5x3x3.omd
+ 6x4x4.omd
+ 3x3x2-CH.omd
+ 3x3x2-h.omd
+ 4x3x2.omd
+ 6x3x3.omd
+ 9x5x1.omd

were taken from the supporting information of "Unit cells for
hexagonal ice" by J. A. Hayward and J. R. Reimers, *J. Chem. Phys.*
**106**, 1518 (1997).
DOI: [10.1063/1.473300](https://doi.org/10.1063/1.473300)

## Proton-ordered ice unit cells

The following structures are unit cells for proton-ordered ice-Ih
crystals.

+ HO-struct1.omd
+ HO-struct2.omd
+ HO-struct3.omd
+ HO-struct4.omd
+ HO-struct5.omd
+ HO-struct6.omd
+ HO-struct7.omd
+ HO-struct8.omd
+ HO-struct9.omd
+ HO-struct10.omd
+ HO-struct11.omd
+ HO-struct12.omd
+ HO-struct3.omd
+ HO-struct14.omd
+ HO-struct15.omd
+ HO-struct16.omd

These structures come from Table 1. in "Quantum-Chemical and
Force-Field Investigations of Ice Ih" by Thomas K. Hirsch and Lars
Ojamae, *J. Phys. Chem. B* **108**, 15856-15864 (2004).
DOI: [10.1021/jp048434u](https://doi.org/10.1021/jp048434u)

NOTE: HO-struct1.omd	is actually ice XI.

When replicated, HO-struct6.omd and HO-struct7.omd create proton
stripes on the basal surfaces.

## Creating large ice crystals from these structures 

In order to generate larger crystals from these structures, use
omd2omd with the -x -y -z flags to replicate these unit cells in the
x, y, and z dimensions.

```
omd2omd -i HO-struct1.omd -o bigCrystal.omd -x 5 -y 3 -z 5
```

Also, while you are unable to cleave the crystals with the current
OpenMD software, you are able to rotate these crystals exposing the
prismatic and secondary prismatic facets using the -p -q -r
functionality of omd2omd.

```
omd2omd -i bigCrystal.omd -o prismFace.omd -p 90 -q 90 -r 0
```

## Sample equilibration scheme

NOTE: These structures are ideal ice crystals, and should be *gently*
equilibrated with whichever water model you choose. Depending on the
model, these starting structures may be more or less favorable. A
sample equilibration scheme might be:

1. Short NPTxyz run at a low temperature, approximately 10 to 50 K,
   with resetTime set to a small time, approximately 10 to 50 fs.
2. Once the pressure tensor elements are nearly zero and the volume is
   oscillating around some average, use affineScale so scale the
   simulation cell to the average volume. Turn resetTime off.
3. Perform an NVT simulation with the targetTemperature set to your
   desired temperature.
4. When the temperature has reached the target, and the total energy
   is oscillating around some average, use thermalizer to scale the
   simulation energy to this average energy.
5. You can now perform NVE simulations.
# Sample zeolite structures

OpenMD supports the **CLAYFF** forcefield described in: "Molecular
Models of Hydroxide, Oxyhydroxide, and Clay Phases and the Development
of a General Force Field," Randall T. Cygan, Jian-Jie Liang, and
Andrey G. Kalinichev, *J. Phys. Chem. B* **108**, pp. 1255-1266 (2004)
DOI: [10.1021/jp0363287](https://doi.org/10.1021/jp0363287), and
provides an example of a solvated ZSM-5 zeolite structure.
# Sample metal surfaces

Contained here are examples of (111) cut surfaces of the coinage
metals, as well as a few of the catalytically active metals.

## slabBuilder

OpenMD also has a utility script which makes creation of these types
of systems trivial. **slabBuilder** is a python script which generates
*sc*, *bcc*, and *fcc* lattices and cleaves the crystals along a
desired *(hkl)* plane. The systems are then reoriented such that the
cleaved facet is presented to the z-dimension of the simulation
box. **slabBuilder** comes with a help message when passed *-h* or
*--help*.

# Air models in OpenMD

The files here can help set up some simple gas-phase simulations of
the common components of air. Three site rigid body models for
N<sub>2</sub>, O<sub>2</sub>, and CO<sub>2</sub> are based on the
TraPPE force field, while the SPC/E water model (also rigid) is
included for simulating various humidity levels. Lennard-Jones
parameters for various noble gases (Ar, He, Ne, Kr) are also included.

Note that if you want heat capacities at high temperatures, you would
need to include vibrational motion for these molecules as well (not
included in these parameters).

Air has a density at sea level and at 15C of roughly 1.225 kg / m<sup>3</sup>
(0.001225 g / cm<sup>3</sup>).  The components of dry air are

| Gas            |  Fraction by Volume |
|----------------|---------------------|
| N<sub>2</sub>  | 0.7809              |
| O<sub>2</sub>  | 0.2095              |
| Ar             | 0.00933             |
| CO<sub>2</sub> | 0.0003              |
| Ne             | 0.000018            |
| He             | 0.000005            |
| Kr             | 0.000001            |


### Rigid Linear Models
+ N2
+ O2
+ CO2

### Rigid three site models
+ SPC/E

### Building a box of air

The `air.inc` and `Air.frc` files contain parameters for most of the
simple components of air.  To build a small box, one might start with
the `mix.omd` which declares the three most prevalent components
(N<sub>2</sub>, O<sub>2</sub>, and Ar).

With all 3 of these files in this sample directory, you can make an air
mixture using randomBuilder:
```
randomBuilder mix.omd -o test.omd --density=0.001225 --nx=5 --ny=5 --nz=5 --molFraction=0.78084 --molFraction=0.20946
```
This creates a `test.omd` structure containing 500 molecules in
roughly the correct proportions (and density).

To warm this mixture up to 15C (and assign initial velocities):
```
thermalizer -t 288 -i test.omd -o warm.omd
```
The simulation is relatively short:
```
openmd warm.omd
```
Then `Dump2XYZ` can output the base atom types and map back to the simulation box:
```
Dump2XYZ -i warm.dump -b -m
```

# References

| Molecular Model| Number of Sites | DOI  |
| ----------:|:---------------:|-----:|
| N<sub>2</sub>  | 3 |[10.1002/aic.690470719](https://doi.org/10.1002/aic.690470719) |
| O<sub>2</sub>  | 3 |[10.1007/s00214-005-0073-1](https://doi.org/10.1007/s00214-005-0073-1) |
| CO<sub>2</sub> | 3 |[10.1002/aic.690470719](https://doi.org/10.1002/aic.690470719)     |
| SPC/E | 3 |[10.1021/j100308a038](https://doi.org/10.1021/j100308a038) |
| CH<sub>4</sub> | 1 |[10.1021/jp972543+](https://doi.org/10.1021/jp972543+) |
| He, Ne, Ar, Kr | 1 | Maitland, G.C., Rigby, M., Smith, E.B., and Wakeham, W.A. (1981) *Intermolecular forces: their origin and determination*. Clarendon Press, Oxford. |

# Sample (.omd) files using Reverse Non-Equilibrium Molecular Dynamics (RNEMD)

## Momentum transport in bulk fluids

The file **shearWater.omd** is a box of 1500 SPC/E water molecules
which has the momentum flux functionality of RNEMD turned on. Notice
in the RNEMD block of the (.omd) file that

```
	fluxType = "Px";
	momentumFlux = 6.0e-7;
```

With these parameters, the applied flux will be an x-axis momentum
flux transferred across the z-dimension of the box. Using this
functionality, one is able to measure the *shear viscosity* of the
liquid at the simulated temperature by relating the imposed momentum
flux to the system's gradient response of the velocity.

## Thermal transport in bulk materials

It is also possible to measure the *thermal conductivity* of a
material using the RNEMD functionality in OpenMD. As an example,
**graphene.omd** is an (.omd) file where two sheets of graphene have a
thermal flux applied accross the long axis of the sheet. In the (.omd)
file the fluxType has been set to a kinetic energy flux, and also that
the kineticFlux is defined.

```
	fluxType = "KE";
	kineticFlux = -6.55e-11;
```
	
The system responds to the thermal flux by developing a temperature
gradient across the z-axis of the system. The thermal conductivity can
be computed by relating the resulting thermal gradient to the imposed
kinetic energy flux.

## Thermal transport across an interface

While computing the thermal conductivity of a bulk material is
certainly of interest, one may also want to investigate interfacial
thermal conductance across an interface. The RNEMD functionality of
OpenMD easily allows for this, and the RNEMD definition block in the
(.omd) file only needs a minor tweaks.  Before our selection, was of
only one component, now we need to allow for more than one atom or
molecule type to be selected.

**gold_water_interface.omd** is an example system where we can compute
the thermal conductivity across an interface, here, a gold / water
interface. Notice in the RNEMD declaration block in the (.omd) file
that the objectSelection is now,

```
	objectSelection = "select SPCE_RB_0 or Au";
```

It is important to make sure your simulation cell is constructed
properly for these kinds of simulations. Since the two RNEMD exchange
regions are defined along the z-dimension at the middle of the
simulation cell and at the far edges (wrapping about the periodic
box), the gold and water need to be properly distributed throughout
the box or else your computation will not give you what you want.


## Simultaneous shearing and thermal transport of bulk materials

The file **2744_shear.omd** is a box of 2744 Argon atoms which has a
simultaneous momentum and kinetic energy flux through the box. Notice
in the RNEMD block of the (.omd) file that

```
	fluxType = "KE+Pvector";
```

and both *kineticFlux* and *momentumFluxVector* are defined.

```
	kineticFlux = -5.0e-6;
	momentumFluxVector = (-2e-7, 0, 0);
```

Application of simultaneous momentum and kinetic energy fluxes result
in both a velocity and thermal gradient response of the system,
allowing for measurement of the shear viscosity of the fluid at a
large number of temperature domains with one simulation.


## Thermal transport in non-periodic systems

OpenMD can perform non-periodic simulations using the Langevin Hull
along with the RNEMD functionality, allowing for computation of
thermal conductance across solvated nanoparticle interfaces. OpenMD
hosts a large number of builder utility scripts which aide in the
construction of these nanoparticles (nanospheres, icosohedra,
cuboctahedra), which can be found in *samples/builders/*.

**NP20_hex_KEflux.omd** is an example of a Gold nanosphere solvated in
  hexane, with a kinetic energy flux that moves thermal energy from
  the solvent into the particle. The RNEMD declaration block in the
  (.omd) file is only slightly different than above. The syntax is
  described in detail in the OpenMD manual.

```
	useRNEMD = "true";
	objectSelection = "select Au or Hexane";
	sphereAradius = 10;
	sphereBradius = 41;
	method = "VSS";
	fluxType = "KE";
	kineticFlux = 1E-5;
	exchangeTime = 10;
	outputBins = 60;

```
The only notable change to the RNEMD declaration block is the addition
of *sphereAradius* and *sphereBradius*, which define the two exchange
regions for the RNEMD moves.
# Metal Oxides in OpenMD
This Sample Library contains Metal Oxide compounds utilizing metals defined in the DR-EAM forcefield and oxygen.  This library contains a variety of compounds, both real and computed, that simulate multiple oxidation states for many of the metals.

All .omd samples contained were obtained from CIF files found on [http://www.crystallography.net/cod/](http://www.crystallography.net/cod/).
# Madelung Energy Sample Calculation

The sample in this directory provides a way of checking the value of
the Madelung Energy for a perfect crystal of NaCl.  The relevant
quantities are:
```
M (Madelung constant) = 1.74756
a (lattice constant)  = 5.65 Angstroms
q^2 / (4 pi e0 a)     = 58.77233   kcal / mol
M q^2 / (4 pi e0 a)   = 102.708173 kcal / mol
```

The file NaCl.omd contains 8000 ions, so the total electrostatic energy
of the perfect crystal in this file should be:
```
V_electrostatic = -821665.38  kcal / mol
```

Using different electrostatic calculation methods, we can get quite
close to this value.

For example, with :
```
cutoffMethod = "shifted_force";
electrostaticScreeningMethod = "damped";
cutoffRadius = 28;
dampingAlpha = 0.14159292;
```

The resultant electrostatic potential is:
```
V_electrostatic = -821667.68 kcal / mol
```

To obtain values for the electrostatic potential in OpenMD, we add the
`ELECTROSTATIC_POTENTIAL` keyword to the end of the statFileFormat:
```
statFileFormat = "TIME|TOTAL_ENERGY|POTENTIAL_ENERGY|KINETIC_ENERGY|TEMPERATURE|PRESSURE|VOLUME|CONSERVED_QUANTITY|ELECTROSTATIC_POTENTIAL";
```
# Sample configurations of liquid alkanes

OpenMD supports the Dipolar Unified-atom Force Field (DUFF) with many
of the parameters coming from the
[TraPPE-UA](http://chem-siepmann.oit.umn.edu/siepmann/trappe/index.html)
forcefield from Siepmann's group.
This sample illustrates a simulation of united-atom (UA) propylene
monomers confined between and all-atom (AA) representation of graphene
sheets.  The force field is based partially on OPLS-AA, and partially
on TraPPE-UA.

Propylene parameters from: C.D. Wick, M.G. Martin, and J.I. Siepmann,
"Transferable potentials for phase equilibria. 4. United-atom
description of linear and branched alkenes and of alkylbenzenes,"
*J. Phys. Chem. B*, **104**, pp. 8008-8016 (2000).
DOI: [10.1021/jp001044x](https://doi.org/10.1021/jp001044x)

Note that the molecule definition in `graphene.inc` includes bonds that span the
box boundaries.  If you want to extend to larger boxes, start with the molecule
definition in `graphene.raw.inc` and then add bonds that span the boundaries of
the new box.

Files included:

1. `graphene.frc`: Force field file
2. `graphene.inc`: Molecule definition include file for a graphene
   sheet (includes cross-box bonding, so assumes a particular box
   geometry)
3. `graphene.raw.inc`: Molecule definition without the additional
   cross-box bonds. Contains unterminated sp<sup>2</sup> carbon atoms.
4. `graphene.omd`: initial OpenMD file for starting a simulation
5. `propylene.omd`: Molecule definition include file for the propylene monomers
6. `propylene.xyz`: A skeletal propylene xyz file for use with Packmol
7. `system.pack`: A [packmol](http://www.ime.unicamp.br/~martinez/packmol/)
   input script that places propylene molecules inside bounds of the
   graphene sheets.

To create the initial configurations, we typically run:

`packmol < system.pack`

which will create a `system.xyz` file containing all of the propylene moleucles.
To get this into a reasonable OpenMD starting configuration, we would run:

`atom2omd -ixyz system.xyz`

This creates `system.omd`, which has periodic box guessed from the bounding box
of the solvent molecules, and this box is generally not the same size as the
box containing the graphene sheets.  To successfully combine the monomers with
the graphene sheets, the `system.omd` file must be edited to modify the Hmat
line to read:

~~~~
        Hmat: {{ 22.23, 0, 0 }, { 0, 42.78, 0 }, { 0, 0, 50 }}
~~~~

Then to combine the propylene monomers with the graphene sheets, eliminating
molecules that overlap, we would run:

`omd-solvator -u graphene.omd -v system.omd -r 3.5 -o combined.omd -n 360 -p 3`

Following this, we typically edit the `combined.omd` file to include the
molecule definitions, set the force field, and set various simulation
parameters:

~~~~
#include "graphene.inc"
#include "propylene.omd"

component{
  type = graphene;
  nMol = 2;
}
component{
  type = propylene;
nMol = 23;
}

forceField = "graphene";
ensemble = NVT;
cutoffMethod = "shifted_force";
electrostaticScreeningMethod = "damped";
cutoffRadius = 9;
dampingAlpha = 0.2;
targetTemp = 300;
tauThermostat = 1000;
dt = 1.0;
runTime = 1e3;
sampleTime = 100;
statusTime = 10;
~~~~
# Utility scripts to aide in construction of OpenMD (.omd) files

Much of the magic outlined below is described in greater detail in the
file **runMe.in** in this directory.

Here we outline a variety of types of systems able to be constructed
using **simpleBuilder**.

+ fcc lattices of a given unit cell
+ nanospheres
+ spherically-capped nanorods
+ pentagonal nanorods
+ icosohedra
+ cuboctahedra
+ truncated cube particles

Note: to use simpleBuilder, you must have the *<MetaData>* section of
an OpenMD (.omd) file as input. Examples of such *<MetaData>* sections
are **one_component.omd**, **three_component.omd**, **gold.omd**, and
**bimetallic.omd**.

## Examples of builder commands shown in **runMe.in** are:

1. Builds an FCC lattice from the <MetaData> block in one_component.omd
2. Builds an FCC lattice from the <MetaData> block in three_component.omd
3. Builds a spherical nanoparticle (FCC) from the <MetaData> block in gold.omd
4. Builds a random alloy spherical nanoparticle (FCC) from the <MetaData>
5. Builds a Au(core)-Ag(shell) spherical nanoparticle (FCC) from the
   <MetaData> block in bimetallic.omd
6. Reverses example 5 by building a Ag(core)-Au(shell) spherical
   nanoparticle. Uses the same <MetaData> block from bimetallic.omd
7. Builds a Au(core)-Ag(shell) spherical nanoparticle (FCC) from the
   <MetaData> block in bimetallic.omd
8. Builds a random alloy spherical nanoparticle with 30% vacancies
   using the <MetaData> block in bimetallic.omd
9. Builds a spherically-capped nanorod (FCC) from the <MetaData> block in gold.omd
10. Builds a pentagonal nanorod from the <MetaData> block in gold.omd
11. Builds a Mackay icosahedral nanoparticle from the <MetaData> block in gold.omd
12. Builds a regular decahedral nanoparticle from the <MetaData> block in gold.omd
13. Builds a ino-decahedral nanorod from the <MetaData> block in gold.omd
14. Builds a cuboctahedral particle from the <MetaData> block in gold.omd
15. Builds a truncated cube particle from the <MetaData> block in gold.omd
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Custom issue template
about: Describe this issue template's purpose here.
title: ''
labels: ''
assignees: ''

---
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:

1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**

- OS: [e.g. iOS]
- Browser [e.g. chrome, safari]
- Version [e.g. 22]

**Smartphone (please complete the following information):**

- Device: [e.g. iPhone6]
- OS: [e.g. iOS8.1]
- Browser [e.g. stock browser, safari]
- Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
