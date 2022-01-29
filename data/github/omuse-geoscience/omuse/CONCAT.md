# OMUSE #

OMUSE stands for Oceanographic MUltipurpose Software Environment. It is a 
package to conduct numerical experiments in oceanography and other Earth 
sciences.

### What is this repository for? ###

This repository contains the source tree for OMUSE. 

### How do I get set up? ###

Easiest way is to use a pip developer install:

- setup a python environment, e.g. using virtualenv, and activate it.
- optionally (see instructions [below](https://github.com/omuse-geoscience/omuse/blob/master/README.md#amuse-developer-install)) first do a develop install of [AMUSE](http://www.amusecode.org), 
- clone the OMUSE repository: `git clone https://github.com/omuse-geoscience/omuse`,
- go into the source directory `cd omuse` and set the environment variable `DOWNLOAD_CODES`, e.g. `export DOWNLOAD_CODES=latest`.
- now, do `pip install -e .` from the root of the package
- type `python setup.py build_codes --inplace` to build the codes. 
- the file `build.log` will report any errors in the build process.

The community codes of OMUSE can be build manually by going into each directory:

 + src/omuse/community/adcirc
 + src/omuse/community/swan
 + etc

and typing: first `make download` (for some) and then `make`

OMUSE has been tested on OSX and linux machines, with ifort and gfortran 
compilers, on desktop machines and on the Carthesius supercomputer.

In addition to the AMUSE dependencies, OMUSE needs/ can use the following 
packages:

 + matplotlib basemap
 + netCDF and netCDF for fortran and the python bindings
 + GRIB_API

#### AMUSE developer install ####

OMUSE depends on the "amuse-framework" package, and the instructions above will install this automatically from pypi. 
If you want to also have a developer install (from the repository source) for AMUSE you should take care that you install the amuse-framework package:

- clone the AMUSE [repository](https://github.com/amusecode/amuse): `git clone https://github.com/amusecode/amuse`
- go into the amuse-framework package directory: `cd amuse/packages/amuse-framework`
- do the developer install from here: `pip install -e .` 

### Documentation ###

Documentation can be found [here](https://omuse.readthedocs.io). In addition the base  [AMUSE documentation](https://amuse.readthedocs.io) can be consulted.

### Reporting issues ###

Issues can be reported at the OMUSE issue tracker; for framework issues, 
report them at the AMUSE [repository](https://github.com/amusecode/amuse).

### Contribution guidelines ###

Contributions are welcome. Note that most framework development happens at 
the AMUSE [repository](https://github.com/amusecode/amuse) A primer for 
writing code interfaces and other documentation can be found on the amuse 
website (www.amusecode.org).
### DALES python wrapper test suite

* To run the test suite type `nosetests -v` within this directory (using your python environment where OMUSE is installed)
* To test the dales low-level fortran interface to OMUSE, run the `run_dales` binary for an appropriate test case (see `dales-repo/cases`)## examples ##

* bubble-notebook.ipynb: Jupyter notebook containing a DALES convective bubble experiment

The Dales code directory contains a number of examples for its use (src/omuse/community/dales/example).

More examples can be found at: https://github.com/omuse-geoscience/omuse-examples
This package installs the i-emic community code for OMUSE.
This package installs the ERA5 interface for OMUSE.

In order to use the ERA5 interface, UID and API key from the CDS portal is 
needed (see the documentation for [CDSAPI](https://pypi.org/project/cdsapi/) 
for more details).
This package installs the SWAN community code for OMUSE.

SWAN is a third-generation wave model, developed at Delft University of Technology, that computes random, short-crested wind-generated waves in coastal regions and inland waters.

The latest version of SWAN can be found at [here](http://swanmodel.sourceforge.net)
This package installs the base OMUSE framework.
This package installs the qgmodel community code for OMUSE.
This package installs all available OMUSE packages.
