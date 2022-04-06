[![Homepage](https://img.shields.io/badge/Home-plumed.org-green.svg)](http://www.plumed.org)
[![Homepage](https://img.shields.io/badge/Google_group-plumed--users-green.svg)](http://groups.google.com/forum/#!forum/plumed-users)
[![codecov](https://codecov.io/gh/plumed/plumed2/branch/master/graph/badge.svg)](https://codecov.io/gh/plumed/plumed2)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/plumed/plumed2.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/plumed/plumed2/context:python)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/plumed/plumed2.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/plumed/plumed2/context:cpp)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](http://www.gnu.org/licenses/lgpl-3.0)
[![Github Releases](https://img.shields.io/github/release/plumed/plumed2.svg)](https://github.com/plumed/plumed2/releases)
[![MacPorts package](https://repology.org/badge/version-for-repo/macports/plumed.svg)](https://repology.org/project/plumed/versions)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/plumed/badges/version.svg)](https://anaconda.org/conda-forge/plumed)
[![AUR package](https://repology.org/badge/version-for-repo/aur/plumed.svg)](https://repology.org/project/plumed/versions)
[![DPorts package](https://repology.org/badge/version-for-repo/dports/plumed.svg)](https://repology.org/project/plumed/versions)
[![FreeBSD port](https://repology.org/badge/version-for-repo/freebsd/plumed.svg)](https://repology.org/project/plumed/versions)
[![Spack package](https://repology.org/badge/version-for-repo/spack/plumed.svg)](https://repology.org/project/plumed/versions)
[![Twitter Follow](https://img.shields.io/twitter/follow/plumed_org.svg?style=social&label=Follow)](https://twitter.com/plumed_org)

Branches and releases
---------------------

Several branches and tags are stored on the git repository.

Branches named `v2.X` correspond to release branches.

Master branch may contain non tested features and is not expected to be used by non-developers.
It typically contains features that will be available on the next release.

Tags named `v2.XbY` correspond to beta releases, use it with care.
Tags named `v2.X.Y` correspond to official releases, use the latest available.

In addition, the repository contains a number of other branches related to specific features.
Please contact the developers that are committing on those branches before basing your work
there, since they might contain temporary work and might be rebased later.
For instance, branch `testdoc` is setup so as to push a test copy of the manual
and is often force pushed.

To report problems found on beta or official releases, use the normal
[plumed-users@googlegroups.com](mailto:plumed-users@googlegroups.com)
mailing list. Please state exactly which version you are using.
To report problems found on `master` branch, use the
[plumed2-git@googlegroups.com](plumed2-git@googlegroups.com) mailing list.
This is also the correct place for discussions about new features etc.
When reporting please provide the git hash (you can obtain it with `git rev-parse HEAD`).

Status
------

Below you find the status on [GitHub Actions](https://github.com/plumed/plumed2/actions) for the release branches.

| Branch   |      Status   | First stable release (year) | Still supported |
|:--------:|:-------------:|:--------:|:------:|
| master   | [![CI](https://github.com/plumed/plumed2/workflows/CI/badge.svg?branch=master)](https://github.com/plumed/plumed2/actions) | 2022 (expected) | / |
| v2.8     | [![CI](https://github.com/plumed/plumed2/workflows/CI/badge.svg?branch=v2.8)](https://github.com/plumed/plumed2/actions)   | 2021 | yes |
| v2.7     | [![CI](https://github.com/plumed/plumed2/workflows/CI/badge.svg?branch=v2.7)](https://github.com/plumed/plumed2/actions)   | 2020 | yes |
| v2.6     | [![CI](https://github.com/plumed/plumed2/workflows/CI/badge.svg?branch=v2.6)](https://github.com/plumed/plumed2/actions)   | 2019 | yes |
| v2.5     | [![CI](https://github.com/plumed/plumed2/workflows/CI/badge.svg?branch=v2.5)](https://github.com/plumed/plumed2/actions)   | 2018 | no |
| v2.4     | [![CI](https://github.com/plumed/plumed2/workflows/CI/badge.svg?branch=v2.4)](https://github.com/plumed/plumed2/actions)   | 2017 | no |
| v2.3     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.3)](https://travis-ci.org/plumed/plumed2)   | 2016 | no |
| v2.2     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.2)](https://travis-ci.org/plumed/plumed2)   | 2015 | no |
| v2.1     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.1)](https://travis-ci.org/plumed/plumed2)   | 2014 | no |
| v2.0     | Not available | 2013 | no |

Content
-------

Here's a description of the content of each file and directory in the root PLUMED directory.

    CHANGES          : change log
    COPYING.LESSER   : license
    Makefile         : makefile
    Makefile.conf.in : template configuration makefile
    PEOPLE           : list of authors
    README.md        : this file
    VERSION          : version file
    astyle           : a local version of astyle, used to format code
    configure        : configuration script
    configure.ac     : configuration script (autoconf)
    developer-doc    : developer documentation
    docker           : directory where Docker is generated
    macports         : directory where Portfiles are generated
    patches          : patch scripts
    python           : python stuff
    regtest          : regression tests, including reference results
    release.sh       : developer utility to publish releases
    scripts          : shell tools
    sourceme.sh.in   : template configuration script
    src              : source code
    test             : examples
    user-doc         : user documentation
    vim              : directory where vim syntax is generated

Required software
-----------------

Required software:

* GNU make.
* C/c++ compiler (c++11 support is required as of version 2.4).
* A modern version of the `patch` command line tool.
* Support for POSIX library `dirent.h`.

Suggested software (libraries are checked by `./configure` and enabled if available):

* MPI library to run parallel simulations. It should be the same library used by your MD code.
* Optimized blas and lapack libraries. They are automatically replaced by an internal version if not available.
* [VMD molfile plugins](http://www.ks.uiuc.edu/Research/vmd/plugins) to read arbitrary file formats. They are automatically replaced by an internal version supporting a few formats if not available.
* [Zlib library](http://zlib.net/) to use compressed data files.
* [Doxygen](http:://www.doxygen.org) to build user manual. Doxygen might need the following packages:
  * Latex to build the pdf user manual.
  * [Graphviz](http://www.graphviz.org) to show class hierarchy in
    developer manual.

Quick compilation instructions
------------------------------

Extensive installation instructions are in the [user documentation](http://www.plumed.org/doc).
Quick instructions:

    ./configure --prefix=$HOME/opt
    make
    make doc # optional
    make test # optional

User documentation can be found at `user-doc/html/index.html`.
Developer documentation can be found at `developer-doc/html/index.html`.
[Pre-compiled documentation](http://www.plumed.org/doc) is available online, so this is only required
if you are working with a modified version of the code!

In order to run PLUMED without installing it you should type `source sourceme.sh`. However,
we recommend installing PLUMED. 
To install it in `$HOME/opt` (directory should be set during `./configure`):

    umask 022
    make install
    
Now you will be able to run plumed using e.g.

    plumed help

If you compiled your own documentation, paths to the installed documentation can be found with command `plumed info --user-doc`.

A sample modulefile with environment variable will be placed in
`$HOME/opt/lib/plumed/src/lib/modulefile`. This can be useful if you want to
install multiple PLUMED versions side by side and select them with env modules.

\defgroup internal-lapack Internal LAPACK

Internal implementation of LAPACK, imported from GROMACS.

The module in src/lapack contains an internal implementation
of LAPACK routines which is automatically imported from GROMACS
using the src/lapack/import.sh script. This set of routines
is compiled when __PLUMED_HAS_EXTERNAL_BLAS is not defined
and allow PLUMED to be used when installed LAPACK libraries
are not available. Notice that the import script
creates a lapack.h file with function declarations which
are used also when installed lapack are employed. This is
done because there are lapack installation written in FORTRAN
that do not provide header files.

Since files are automatically generated, do not edit them directly.
In case you need PLUMED specific modifications
please do it by modifying the import script.

Within the PLUMED doxygen (this page) the available
macros are listed but not documented. Have a look
at the corresponding documentation at http://www.netlib.org/lapack

ANN (Artificial Neural Network) function for plumed
====================

This is plumed ANN function (annfunc) module.  It implements `ANN` class, which is a subclass of `Function` class.  `ANN` class takes multi-dimensional arrays as inputs for a fully-connected feedforward neural network with specified neural network weights and generates corresponding outputs.  The `ANN` outputs can be used as collective variables, inputs for other collective variables, or inputs for data analysis tools.  

## Installation

Enable compilation by adding the `--enable-modules=annfunc` to the configure command.

## Usage

It is used in a similar way to [other plumed functions](https://www.plumed.org/doc-v2.5/user-doc/html/_function.html).  To define an `ANN` function object, we need to define following keywords:

- `ARG` (string array): input variable names for the fully-connected feedforward neural network

- `NUM_LAYERS` (int): number of layers for the neural network

- `NUM_NODES` (int array): number of nodes in all layers of the neural network

- `ACTIVATIONS` (string array): types of activation functions of layers, currently we have implemented "Linear", "Tanh", "Circular" layers, it should be straightforward to add other types as well

- `WEIGHTS` (numbered keyword, double array): this is a numbered keyword, `WEIGHTS0` represents flattened weight array connecting layer 0 and layer 1, `WEIGHTS1` represents flattened weight array connecting layer 1 and layer 2, ...  An example is given in the next section.

- `BIASES` (numbered keyword, double array): this is a numbered keyword, BIASES0 represents bias array for layer 1, BIASES1 represents bias array for layer 2, ...

Assuming we have an `ANN` function object named `ann`, we use `ann.node-0, ann.node-1, ...` to access component 0, 1, ... of its outputs (used as collective variables, inputs for other collective variables, or data analysis tools).

## Examples

Assume we have an ANN with numbers of nodes being [2, 3, 1], and weights connecting layer 0 and 1 are

```
[[1,2],
[3,4],
[5,6]]
```

weights connecting layer 1 and 2 are

```
[[7,8,9]]
```

Bias for layer 1 and 2 are

```
[10, 11, 12]
```

and 

```
[13]
```

respectively.

All activation functions are `Tanh`.

Then if input variables are `l_0_out_0, l_0_out_1`, the corresponding `ANN` function object can be defined using following plumed script: 

```
ann: ANN ARG=l_0_out_0,l_0_out_1 NUM_LAYERS=3 NUM_NODES=2,3,1 ACTIVATIONS=Tanh,Tanh  WEIGHTS0=1,2,3,4,5,6 WEIGHTS1=7,8,9  BIASES0=10,11,12 BIASES1=13
```

This plumed script can be generated with function `Plumed_helper.get_ANN_expression()` in [this](https://github.com/weiHelloWorld/plumed_helper/blob/master/plumed_helper.py) repository.  Following is the Python code using this function to generate the script above:

```Python
from plumed_helper import Plumed_helper
ANN_weights = [np.array([1,2,3,4,5,6]), np.array([7,8,9])]
ANN_bias = [np.array([10, 11, 12]), np.array([13])]
Plumed_helper.get_ANN_expression('ANN', node_num=[2, 3, 1], 
                                 ANN_weights=ANN_weights, ANN_bias=ANN_bias,
                                 activation_list=['Tanh', 'Tanh'])
```

## Authors

Wei Chen (UIUC, weichen9@illinois.edu) and Andrew Ferguson (University of Chicago, andrewferguson@uchicago.edu)

## Copyright

See ./COPYRIGHT
Experiment Directed Simulation (EDS)
====================================


Install
------------------------------------
Enable compilation by adding the `--enable-modules=+eds`
to the configure command.


Documentation
------------------------------------
See the generated documentation for information on
using


Authors
------------------------------------
Glen Hocky (University of Chicago) <hockyg@uchicago.edu>
Andrew White (University of Rochester) <andrew.white@rochester.edu>


Copyright
------------------------------------
See ./COPYRIGHT
Permutation Invariat Vector (PIV)
====================================


Install
------------------------------------
Enable compilation by adding the `--enable-modules=+piv`
to the configure command.


Documentation
------------------------------------
See the generated documentation for information on
using


Authors
------------------------------------
Silvio Pipolo (University of Lille) <silvio.pipolo@univ-lille.fr>
Fabio Pietrucci (Sorbonne University, Paris) <fabio.pietrucci@impmc.upmc.fr>


Copyright
------------------------------------
See ./COPYRIGHT
Welcome to the plumed2 wiki!

This fork of PLUMED has eABF/DRR implementation.

**Requirements**

Boost::serialization and C++11 compiler

**Compiling instruction:**

After clone this repository, please cd to the plumed2 directory and run:

1. `autoconf`
2. `./configure --enable-boost_serialization --enable-modules=drr`
3. Modify your Makefile.conf and add `-lboost_serialization` in `DYNAMIC_LIBS=`
4. `make` and `sudo make install`

**Usage**

Run `make doc` and the usage is in `user-doc/html/_d_r_r.html`

**Authors**

Chen Haochuan (All files in drr module except colvar_UIestimator.h)

Fu Haohao (colvar_UIestimator.h)

**COPYRIGHT**

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


\defgroup internal-blas Internal BLAS

Internal implementation of BLAS, imported from GROMACS.

The module in src/blas contains an internal implementation
of BLAS routines which is automatically imported from GROMACS
using the src/blas/import.sh script. This set of routines
is compiled when __PLUMED_HAS_EXTERNAL_BLAS is not defined
and allow PLUMED to be used when installed BLAS libraries
are not available. Notice that the import script
creates a blas.h file with function declarations which
are used also when installed blas are employed. This is
done because there are blas installation written in FORTRAN
that do not provide header files.

Since files are automatically generated, do not edit them directly.
In case you need PLUMED specific modifications
please do it by modifying the import script.

Within the PLUMED doxygen (this page) the available
macros are listed but not documented. Have a look
at the corresponding documentation at http://www.netlib.org/blas

maze
================================================================================
This version of PLUMED2 has the maze module implemented.

Install
--------------------------------------------------------------------------------
Enable the compilation of maze by adding the `--enable-modules=maze` to the 
configure command.

Page
--------------------------------------------------------------------------------
See [this link](http://maze-code.github.io) for further information.

Documentation
--------------------------------------------------------------------------------
Run `make doc`; the documentation should be in `user-doc/html/_m_a_z_e.html`.

Author
--------------------------------------------------------------------------------
Jakub Rydzewski (Nicolaus Copernicus University) <jr@fizyka.umk.pl>

Copyright
--------------------------------------------------------------------------------
This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU Lesser General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any 
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE.  

See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along 
with this program.  If not, see <http://www.gnu.org/licenses/>.
Instructions for using Artistic Style are included in the *doc* directory.

The file **install.html** contains instructions for compiling and
installing Artistic Style.

The file **astyle.html**' contains information on using Artistic Style.

The files **news.html** and **notes.html** contain information on changes
made to the various releases.
@page CHANGES-2-3 Version 2.3

## Version 2.3 (Dec 12, 2016)

Version 2.3 contains several improvements with respect to 2.2. Users currently working with 2.2
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files.

Below you find a list of all the changes with respect to version 2.2.
Notice that version 2.3 includes already all the fixes in branch 2.2 up to 2.2.3 indicated in \ref CHANGES-2-2 .

Changes from version 2.2 which are relevant for users:
- Changes leading to incompatible behavior:
  - \ref COMMITTOR can now be used to define multiple basins, but the syntax has been changed
  - Syntax for \ref SPRINT and \ref DFSCLUSTERING has changed.
    We have separated the Actions that calculate the contact matrix from these actions.  These actions thus now take a contact
    matrix as input.  This means that we these actions can be used with contact matrices that measures whether or not a pair of atoms
    are hydrogen bonded.  For more details on this see \ref contactmatrix.  For clustering the output can now be passed to the actions
    \ref CLUSTER_PROPERTIES, \ref CLUSTER_DIAMETER, \ref CLUSTER_NATOMS, \ref OUTPUT_CLUSTER and \ref CLUSTER_DISTRIBUTION.  These
    provide various different kinds of information about the connected components found by clustering 
  - In \ref driver masses and charges are set by default to NaN.
    This makes it less likely to do mistakes trying to compute centers of mass or electrostatic-dependent variables
    when masses or charges were not set. To compute these variables from the driver you are now forced to use
    `--pdb` or `--mc`.
  - In rational switching functions, by default MM is twice NN. This is valid both in \ref switchingfunction with expanded
    syntax and when specifying MM on e.g. \ref COORDINATION
  - Patch script `plumed patch` now patches by default with `--shared`. This should make the procedure more robust (see \issue{186}).
  - Faster \ref GYRATION but new default behavior is not mass weighted
  - When using \ref HISTOGRAM you now output the accumulated grid using \ref DUMPGRID or \ref DUMPCUBE to get the free energy you use
    the method \ref CONVERT_TO_FES.  These changes allow one to use grids calculated within PLUMED in a work flow of tasks similarly to 
    the way that you can currently use Values.
  - The way that reweighting is performed is now different.  There are three separate actions \ref REWEIGHT_BIAS, \ref REWEIGHT_TEMP and
    \ref REWEIGHT_METAD.  These actions calculate the quantities that were calculated using the keywords REWEIGHT_BIAS and REWEIGHT_TEMP that
    used to appear in the old HISTOGRAM method.  Now those these methods can be used in any methods that calculate ensemble averages for
    example \ref HISTOGRAM and \ref AVERAGE
  - Manual is now build with locally compiled plumed
  - Removed CH3SHIFT
  - \ref CS2BACKBONE is now native in PLUMED removing the need to link ALMOST, small syntax differences
  - \ref CS2BACKBONE, \ref NOE, \ref RDC, removed the keyword ENSEMBLE: now ensemble averages can only be calculated using \ref ENSEMBLE
  - \ref RDC, syntax changes
  - It is not possible anymore to select modules using `modulename.on` and `modulename.off` files. Use `./configure --enable-modules` instead.
  - Removed IMD modules. In case someone is interested in restoring it, please contact the PLUMED developers.
- New actions:
  - \ref FIXEDATOM
  - \ref HBOND_MATRIX
  - \ref CLUSTER_PROPERTIES
  - \ref CLUSTER_DIAMETER
  - \ref CLUSTER_NATOMS
  - \ref OUTPUT_CLUSTER
  - \ref CLUSTER_DISTRIBUTION
  - \ref ROWSUMS
  - \ref COLUMNSUMS
  - \ref UPDATE_IF
  - \ref DUMPGRID
  - \ref DUMPCUBE
  - \ref CONVERT_TO_FES
  - \ref INTERPOLATE_GRID
  - \ref FIND_CONTOUR
  - \ref FIND_SPHERICAL_CONTOUR
  - \ref FIND_CONTOUR_SURFACE
  - \ref AVERAGE
  - \ref REWEIGHT_BIAS
  - \ref REWEIGHT_TEMP
  - \ref REWEIGHT_METAD
  - \ref PCA
  - \ref PRE
  - \ref STATS
  - \ref METAINFERENCE
  - \ref LOCALENSEMBLE
  - \ref FRET
  - \ref RESET_CELL
  - \ref JCOUPLING
  - \ref ERMSD
- New features in MD patches (require re-patch):
  - Patch for amber 14 now passes charges with appropriate units (fixes \issue{165}). Notice that
    the patch is still backward compatible with older PLUMED version, but the charges will only be passed
    when using PLUMED 2.3 or later.
  - Patch for GROMACS 5.1 incorporates Hamiltonian replica exchange, see \ref hrex
  - Gromacs 2016, 5.1.x, 5.0.x, flush the plumed output files upon checkpointing
  - Added patch for Gromacs 2016.1
  - gromacs 5.1.x patch updated to 5.1.4
  - Removed the patch for Gromacs 4.6.x 
  - LAMMPS patch updated to support multiple walkers and report plumed bias to LAMMPS (thanks to Pablo Piaggi).
- New features for existing actions:
  - The SPECIES and SPECIESA keyword in MultiColvars can now take a multicolvar as input.  This allows one
    to calculate quantities such as the Q4 parameters for those atoms that have a coordination number greater
    than x.
  - Added MATHEVAL type in \ref switchingfunction
  - Added Q type native contacts in \ref switchingfunction (thanks to Jan Domanski).
  - \ref COMMITTOR can now be used to define multiple basins
  - The number of atoms admitted in \ref BRIDGE has been significantly increased, see \issue{185}.
  - \ref driver now allows --trajectory-stride to be set to zero when reading with --ixtc/--itrr. In this case, step number is read from the trajectory file.
  - \ref METAD and \ref PBMETAD can now be restarted from a GRID 
  - Added keywords TARGET and DAMPFACTOR in \ref METAD
  - When using \ref METAD with file-based multiple walkers and parallel jobs (i.e. mpirun) extra suffix is not added (thanks to Marco De La Pierre).
  - \ref ENSEMBLE added keywords for weighted averages, and calculation of higher momenta
  - \ref MOLINFO now allows single atoms to be picked by name.
  - \ref FIT_TO_TEMPLATE now supports optimal alignment.
  - \ref CONSTANT added the possibility of storing more values as components with or without derivatives
  - \ref PUCKERING now supports 6 membered rings.
  - Extended checkpoint infrastructure, now \ref METAD and \ref PBMETAD will write GRIDS also on checkpoint step (only the GROMACS patch
    is currently using the checkpointing interface)
- Other features:
  - Added a plumed-config command line tool. Can be used to inspect configuration also when cross compiling.
  - Added a `--mpi` option to `plumed`, symmetric to `--no-mpi`. Currently, it has no effect (MPI is initialized by default when available).
  - PLUMED now generate a VIM syntax file, see \ref VimSyntax
  - The backward cycle is now parallelized in MPI/OpenMP in case many collective variables are used.
  - GSL library is now searched by default during `./configure`.
  - Tutorials have been (partially) updated to reflect some of the changes in the syntax
  - Parser now reports errors when passing numbers that cannot be parsed instead of silently replacing their default value. See \issue{104}.
  - More and more documentation
- Bug fixes:
- Fixed a bug in \ref PBMETAD that was preventing the writing of GRIDS if a hill was not added in that same step 

For developers:
- IMPORTANT: BIAS can now be BIASED as well, this changes can lead to some incompatibility: now the "bias" component is always defined automatically
  by the constructor of Bias as a componentWithDerivatives, derivatives are automatically obtained by forces. The main change is that you don't have to define
  the bias component anymore in your constructor and that you can use setBias(value) to set the value of the bias component in calculate. 
- Added new strings for plumed cmd: setMDMassUnits, setMDChargeUnits, readInputLine, performCalcNoUpdate, update and doCheckPoint.
- Easier to add actions with multiple arguments
- New functions to access local quantities in domain decomposition
- Active modules to enable regtests are chosen using `plumed config`.
- A script is available to check if source code complies plumed standard. Notice that this script is run together with cppcheck on travis-ci.
- Cppcheck on travis-ci has been updated to 1.75. Several small issues triggering errors on 1.75 were fixed (e.g. structures passed by value
   are now passed by const ref) and false positives marked as such.
- Added coverage scan.

## Version 2.3.1 (Mar 31, 2017)

- Fix to FIT_TO_TEMPLATE as in 2.2.5. Notice that in 2.3.0 also the case with TYPE=OPTIMAL was affected. This is fixed now.
- small change in \ref CS2BACKBONE to symmetrize the ring current contribution with respect to ring rotations (also faster)
- fixed `plumed-config` that was not working.
- log file points to the `config.txt` files to allow users to check which features were available in that compiled version.
- `make clean` in root dir now also cleans `vim` sub-directory.
- Updated gromacs patch to version 2016.3 

For developers:
- Cppcheck on travis-ci has been updated to 1.77.
- Doxygen on travis-ci has been updated to 1.8.13

## Version 2.3.2 (Jun 12, 2017)

See branch \branch{v2.3} on git repository.

- Resolved problem with nan in \ref SMAC with SPECIESA and SPECIESB involving molecules that are the same
- PDB reader is now able to read files with dos newlines (see \issue{223}).
- Fixed bug in \ref CS2BACKBONE (v2.3.1) related to ring currents of HIS and TRP
- Fixed bug in if condition in \ref PCAVARS so that you can run with only one eigenvector defined in input 
- Fixed bug with timers in \ref sum_hills \issue{194}.
- Fixed bug when using \ref MOVINGRESTRAINT with periodic variables such as \ref TORSION \issue{225}.
- Fixed bug in \ref HBOND_MATRIX that used to appear when you used DONORS and ACCEPTORS with same numbers of atoms 
- Fixed bug in \ref DISTANCES that appears when using BETWEEN and link cells.
- Prevented users from causing segfaults by storing derivatives without LOWMEM flag.  In these cases PLUMED crashes with meaningful errors.
- Fixed bug in \ref HISTOGRAM that causes NaNs when using KERNEL=DISCRETE option
- Fixed a bug in the parser related to braces, see \issue{229}
- Fixed a bug that appeared when using \ref Q3, \ref Q4 and \ref Q6 with LOWEST or HIGHEST flag
- Fixed a bug that appears when you use \ref MFILTER_LESS as input to \ref COORDINATIONNUMBER with SPECIESA and SPECIESB flags
- Fixed a bug that was making flushing when gromacs checkpoints not functional (thanks to Summer Snow).
- Fixed a bug affecting \ref EXTENDED_LAGRANGIAN and \ref METAD with ADAPT=DIFF when using an argument
  with periodicity (min,max) such that min is different from -max.
  This does not affect normal \ref TORSION, but would affect \ref PUCKERING component phi
  with 6-membered rings. In addition, it would affect any variable that is created by the user with a periodicity
  domain not symmetric around zero. See \issue{235} (thanks to Summer Snow for reporting this bug).
- Fixed numerical issue leading to simulations stuck (LatticeReduction problem) with intel compiler and
  large simulation cells.
- Fixed a bug affecting \ref LOCAL_AVERAGE and outputting all multicolvars calculated by \ref Q6 with \ref DUMPMULTICOLVAR
- `plumed info --user-doc` and `plumed info --developer-doc` now fall back to online manual when local doc is not installed,
  see \issue{240}.

For developers:
- IMPORTANT: we started to enforce code formatting using astyle. Check the developer documentation to learn how to
  take care of not-yet-formatted branches.
- plumedcheck validation has been made stricter. All the checks are now described in the developer manual.
- New flag `--disable-libsearch` for `configure`, allowing an easier control of linked libraries when installing PLUMED
  with a package manager such as MacPorts.
- Added `--disable-static-patch` to `./configure` to disable tests related to static patching. It can be used
  when static patching is not needed to make sure a wrong c++ library is not linked by mistake.
- Using `install_name_tool` to fix the name of the installed library on OSX. Allows linking the PLUMED
  shared library without explicitly setting `DYLD_LIBRARY_PATH`.
- Added environment variable `PLUMED_ASYNC_SHARE` to enforce synchronous/asynchronous atom sharing (mostly for debug purpose).
- On travis-ci, using ccache to speedup builds.
- On travis-ci, added a regtest using Docker with gcc6 and MPI.
- On travis-ci, docs for unofficial or unsupported branches are set not to be indexed by search engines (see \issue{239})
- Cppcheck on travis-ci has been updated to 1.79.

## Version 2.3.3 (Oct 3, 2017)

For users:
- Fixed a bug in \ref switchingfunction MATHEVAL, leading to inconsistent results when using OpenMP with multiple threads (see \issue{249}).
- \ref FIT_TO_TEMPLATE now reports when it is used with a reference file with zero weights.
- Fixed logging of \ref UNITS (thanks to Omar Valsson).
- Fixed a possible bug with \ref EFFECTIVE_ENERGY_DRIFT and domain decomposition with a domain containing zero atoms.


For developers:
- Fixed a bug in `./configure --disable-libsearch` when searching for molfile plugins.
- Cppcheck on travis-ci has been updated to 1.80.
- Configure script now has a list of better alternatives to find a working `ld -r -o` tool to merge object files.
  This solves linking issues on some peculiar systems (see \issue{291}, thanks to Massimiliano Culpo). 
- Using `install_name_tool` also on non-installed libraries. This makes it possible to link them and later
  find them without explicitly setting `DYLD_LIBRARY_PATH`. This should also make the `DYLD_LIBRARY_PATH` irrelevant.
  Notice that `DYLD_LIBRARY_PATH` is not well behaved in OSX El Capitan.

## Version 2.3.4 (Dec 15, 2017)

For users:
- GROMACS patch updated to gromacs-2016.4. This patch was also fixed in order to properly work with \ref ENERGY (see \issue{316})
  and to implement `-hrex` option (see \issue{197}).
- Patch for GROMACS 5.1.4 updated to fix an error with \ref ENERGY (see \issue{316}).
- Solved a bug in \ref ERMSD leading to incorrect results when using non-default length units (e.g. with `UNITS LENGTH=A`).

For developers:
- Regtest script also reports when exitcode different from zero is returned.
- Patch script reports errors returning a nonzero exit code.
- cppcheck update to 1.81
- Solved small bug in stored PLUMED_ROOT directory as obtained from statically patched MD codes.
  Namely, the compilation directory was stored rather than the installation one.

## Version 2.3.5 (Mar 2, 2018)

For users:
- Fixed `plumed partial_tempering` to agree with GROMACS conventions for the choice of dihedral angles (see \issue{337}).
  Should be irrelevant for the vast majority of cases.
- Fixed small bug in regexp parser - the part outside the parentheses was just ignored.

For developers:
- Doxygen on travis-ci has been updated to 1.8.14.
- Embedded astyle updated to 3.1.
- `make clean` now correctly removes the `src/lib/plumed` executable.

## Version 2.3.6 (Jul 2, 2018)

For users:
- Fixed a problem leading to NaN derivatives of \ref switchingfunction `Q` when distance between two atoms is large.
- GROMACS patch updated to gromacs-2016.5.
- `./configure` crashes if prefix is set to present working directory (notice that this choice was already leading to issues).
- \ref DUMPATOMS reports an error when trying to write xtc/xdr files without the xdrfile library installed.
- Fixed a bug appearing when using \ref PATH or \ref GPROPERTYMAP with virtual atoms without simultaneously using the same
  atoms in a different action.
- Fixed incorrect format of the pdb file written by \ref PCA (see \issue{363}).
- Fixed behavior of natural units. When an MD code asks for natural units, it is not necessary to also set units within PLUMED using \ref UNITS (see \issue{364}).

For developers:
- Fixed small issue in debug options of \ref driver (see \issue{245}).
- `plumed patch -e` now accepts a name closely matching the patch name (e.g. `plumed patch -e gromacs2016.5` will try to patch
  even if the stored patch is for `gromacs-2016.4`). This simplifies managing Portfiles. Nothing changes when picking the patch
  from the interactive menu.
- Install newer ccache on travis-ci, build faster.
- Small fix in provided env modules (`PLUMED_VIMPATH` is set also when shared libraries are disabled).

## Version 2.3.7 (Oct 5, 2018)

For users:
- Fixed flag DETAILED_TIMERS in \ref DEBUG (flag was ignored and detailed timers always written).
- Small fix in \ref DUMPMASSCHARGE (atoms are now correctly requested only at first step).

## Version 2.3.8 (Dec 19, 2018)

\plumednotmaintained

For users:
- Fixed some openMP regression (some related to the whole codes and some specifics for Coordination and Multicolvar), this were compiler dependent so not all users may have experienced them
- Fixed an issue with \ref CS2BACKBONE when more than 2 chains were used
- Fixed memory leak in \ref RDC.
- Fixed segmentation fault with more than two CVs in reweighting \ref METAD (see \issue{399}, thanks to Fiskissimo).

For developers:
- Small fix in LDFLAGS when enabling coverage.
- Fixed order of flags in tests for static linking done by configure (see \issue{407}).
- Fixed the way paths are hard-coded so as to facilitate conda packaging (see \issue{416}).


*/
@page CHANGES-2-7 Version 2.7
  
## Version 2.7.0 (Dec 23, 2020)

Changes from version 2.6 which are relevant for users:

- Changes leading to differences with previous versions
  - The definition of the omega angle has been modified to adhere to the IUPAC standard (i.e. with the previous amino acid)

- New contributed modules:
  - A new Funnel module by Stefano Raniolo and Vittorio Limongelli 
     - \ref FUNNEL_PS 
     - \ref FUNNEL 
  - A new Infinite Switch Simulated Tempering in Force module by Glen Hocky
     - \ref FISST
  - A new OPES module by Michele Invernizzi
     - \ref OPES_METAD

- New actions:
  - \ref ENVIRONMENTSIMILARITY from Pablo Piaggi
  - \ref PROJECTION_ON_AXIS from
  - \ref FUNCPATHGENERAL from

- Other improvements:
  - \ref MOLINFO action can now be used multiple times. Every action doing a search will use the latest
    appearance. See \issue{134}.
  - Neighbor lists are now OpenMP and MPI parallel so improving the scalability of all actions employing them
  - It is now possible to pass pdb files with all weights set to zero. Instead of reporting an error,
    PLUMED will now assume they are all equal to 1/n, where n is the number of atoms (see \issue{608}).
  - All the examples in the manual are now displayed with contextual help and regularly tested for correctness.
  - A tool to build PLUMED input directly within a python script has been added (see \issue{611}
    and documentation for class `plumed.InputBuilder()`).
  - Python function `plumed.read_as_pandas()` now also accepts an argument `index_col`.
  - Lepton arithmetics can be used also when reading integers (e.g., `METAD PACE=2*5`, see \issue{614}).

- GROMACS:
  - When using `-hrex` flag, the neighbor lists are update automatically at every exchange step.
    This relaxes the requirements on the choice of `-replex` stride (see \issue{579}, thanks to Chang Junhan).

- Changes in the DRR module
  - Support multi-time stepping. Now the STRIDE keyword should work with DRR correctly.
  - Support reflecting boundary conditions, which should be a better solution to the boundary effect of eABF in non-periodic cases. You can use REFLECTINGWALL to enable it.
  - Stop the simulation when the temperature is not passed from the MD engine to PLUMED. In this case, users should set the temperature by the TEMP keyword.

- Changes in the ISDB module
  - There is a new option for OPTSIGMAMEAN, SEM_MAX that allows to automatically determine an optimal value for SIGMA_MAX

- Changes in the VES module
  - Small changes to TD_MULTICANONICAL and TD_MULTITHERMAL_MULTIBARIC. Bug fix concerning the calculation of the logarithm of the target distribution. Added the keyword EPSILON to avoid dealing with regions of zero target probability.

For developers:
- small fix in `Plumed.h` too avoid unique global symbols (see \issue{549})
- Added `cmd("readInputLines")` to allow reading input from a buffer with comments and continuation lines (see \issue{571}).
- fixed error when the install prefix contained unicode characters

## Version 2.7.1 (Apr 16, 2021)

- Includes all fixes up to 2.6.3
- In python interface, fixed usage of python arrays to allow compatibility with PyPy.
- New/updated patches:
  - updated patch for gromacs-2020.5
  - new patch for gromacs-2021
    - this should work with multiple-time stepping (plumed forces are integrated with the smallest time step, plumed can internally implement a multiple-time step if needed). 
    - Modular simulator is still not supported
    - hrex, lambda cv and replica-exchange are not yet tested
 
## Version 2.7.2 (Jul 27, 2021)

- Includes all fixes up to 2.6.4
- Fixed a bug in the `-hrex` implementation for GROMACS 2020 and 2021 (see #691, thanks to Chang Junhan).
- Changes in the OPES module
  - the CALC_WORK option now outputs the accumulated work, as in METAD, instead of the work done in the last bias update

## Version 2.7.3 (Dec 1, 2021)

- Includes all fixes up to 2.6.5
- GROMACS patches now take a note of the used PLUMED version in the GROMACS log (see \issue{737})
- GROMACS 2021 patch renamed to 2021.4 for consistency.

## Version 2.7.4 (Feb 22, 2022)

- Includes all fixes up to 2.6.6

## Version 2.7.5

- Minor fixes in error reporting.
- Fix in building python package with MacPorts and MacOS 11.

@page CHANGES-2-2 Version 2.2

Version 2.2 (Oct 13, 2015)
----------------------------

Version 2.2 contains several improvements with respect to 2.1. Users currently working with 2.1
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files. In 2.2 we restored more features of 1.3
that were missing in 2.1, so users still working with 1.3 could opt for an upgrade.
A tutorial explaining how to move from PLUMED 1 to PLUMED 2 is available (see \ref moving).

Below you find a list of all the changes with respect to version 2.1.
Notice that version 2.2 includes already all the fixes in branch 2.1 up to 2.1.4 indicated in \ref CHANGES-2-1 .

Changes from version 2.1 which are relevant for users:
- Changes leading to incompatible behavior:
  - Labels of quantities calculates by \ref SPRINT have changed from <i>label</i>.coord_<i>num</i> to <i>label</i>.coord-<i>num</i>
  - \ref METAD with WALKERS_MPI now writes a single hills file, without suffixes
  - removed the ./configure.sh script of v2.0.x, now plumed can only be configured using autotools (./configure)
  - \ref COM, \ref CENTER, and \ref GYRATION now automatically make molecules whole. In case you do not want them to do it, use NOPBC flag,
    which recovers plumed 2.1 behavior
  - Some MD code could now automatically trigger restart (e.g. gromacs when starting from cpt files). This can be overwritten using
    \ref RESTART NO.
  - Replica suffixes are now added by PLUMED *before* extension (e.g. use plumed.0.dat instead of plumed.dat.0)
  - When using \ref switchingfunction the STRETCH keyword is now implicit. NOSTRETCH is available to enforce the old behavior.
- Module activation can now be controlled during configure with `--enable-modules` option.
- Almost complete refactoring of installation procedure. Now
  DESTDIR and other standard autoconf directories (e.g. bindir) are completely supported.
  Additionally, everything should work properly also when directory names include spaces (\issue{157}).
  Finally, compiler is not invoked on install unless path are explicitly changed (\issue{107}).
- Related to installation refactoring, upon install a previously installed PLUMED is not removed.
  This is to avoid data loss if prefix variable is not properly set
- Several changes have been made in the Makefile.conf that makes it not compatible with those
  packaged with plumed 2.0/2.1. Please use ./configure to generate a new configuration file.
- Added partial OpenMP parallelization, see \ref Openmp
- Added multiple time step integration for bias potentials, see \ref MTS
- Link cells are now used in all multicolvars that involve \ref switchingfunction.  The link cell cutoff is
  set equal to 2.*\f$d_{\textrm{max}}\f$.  Where \f$d_{\textrm{max}}\f$ is the (user-specified) point at which
  the switching function goes to zero. Users should always set this parameter when using a switching function
  in order to achieve optimal performance.
- DHENERGY option is no longer possible within \ref DISTANCES.  You can still calculate the DHENERGY colvar by using \ref DHENERGY
- Reweighting in the manner described in \cite Tiwary_jp504920s is now possible using a combination of the \ref METAD and \ref HISTOGRAM actions.  The relevant keywords in \ref METAD are REWEIGHTING_NGRID and REWEIGHTING_NHILLS.  The \f$c(t)\f$ and the appropriate weight to apply to the configurations are given by the values labeled rct and rbias. 
- News in configure and install:
  - ./configure now allows external BLAS to be used with internal LAPACK. This is done automatically if only BLAS are available,
    and can be enforced with --disable-external-lapack.
  - ./configure supports --program-prefix, --program-suffix, and --program-transform-name.
  - make install supports DESTDIR and prefix.
  - Environment variables PLUMED_LIBSUFFIX and PLUMED_PREFIX are deprecated and will be removed in a later version.
- New actions
  - \ref DUMPMASSCHARGE to dump a file with mass and charges during MD.
  - \ref EFFECTIVE_ENERGY_DRIFT to check that plumed forces are not screwing the MD integrator.
  - \ref EXTENDED_LAGRANGIAN : in combination with  \ref METAD it implements metadynamics with Extended Lagrangian; standalone it implements TAMD/dAFED.
  - \ref DFSCLUSTERING calculate the size of clusters 
  - \ref DUMPMULTICOLVAR print out a multicolvar
  - \ref MFILTER_LESS filter multicolvar by the value of the colvar
  - \ref MFILTER_MORE 
  - \ref MFILTER_BETWEEN
  - \ref PCARMSD PCA collective variables using OPTIMAL rmsd measure
  - \ref PCAVARS PCA collective variables using any one of the measures in reference
  - \ref GRADIENT can be used to calculate the gradient of a quantity.  Used to drive nucleation
  - \ref CAVITY
  - \ref PUCKERING implemented for 5-membered rings (thanks to Alejandro Gil-Ley).
  - \ref WRAPAROUND to fix periodic boundary conditions.
- New features for existing actions:
  - Keywords UPDATE_FROM and UPDATE_UNTIL to limit update step in a defined time window, available only for actions where it would be useful.
  - Keyword UNNORMALIZED for \ref HISTOGRAM.
  - Possibility to use Tiwary-Parrinello reweighting for \ref METAD
  - Keywords for \ref GROUP (REMOVE, SORT, UNIQUE) to allow more flexible editing of groups.
  - \ref DUMPATOMS now supports dumping xtc and trr files (requires xdrfile library).
  - \ref driver can now read xtc and trr files also with xdrfile library.
  - \ref driver accepts a --mc flag to read charges and masses from a file produced during
    molecular dynamics with \ref DUMPMASSCHARGE
  - Possibility to enable or disable \ref RESTART on a per action basis, available only for actions where it would be useful.
  - \ref MOLINFO now supports many more special names for rna and dna (thanks to Alejandro Gil-Ley).
  - VMEAN and VSUM allow one to calculate the sum of a set of vectors calculated by VectorMultiColvar.  Note these
  can also be used in tandem with \ref AROUND or \ref MFILTER_MORE to calculate the average vector within a particular
  part of the cell or the average vector among those that have a magnitude greater than some tolerance
  - New way of calculating the minimum value in multicolvars (ALT_MIN). This is less susceptible to overflow for certain 
    values of \f$\beta\f$.  
  - New keywords for calculating the LOWEST and HIGHEST colvar calculated by a multicolvar
  - Added components to \ref DIPOLE (\issue{160}).
- Other changes:
  - File reader now supports dos newlines as well as files with no endline at the end.

For developers:

- In order to be able to use openMP parallelism within multicolvar, secondarystructure, manyrestraints and crystallisation
we had to make some substantial changes to the code that underlies these routines that is contained within vesselbase. In 
particular we needed to get rid of the derivatives and buffer private variables in the class ActionWithVessel.  As a consequence
the derivatives calculated in the various performTask methods are stored in an object of type MultiValue.  Within multicolvar
this is contained within an object of type AtomValuePack, which stores information on the atom indices.  If you have implemented
a new multicolvar it should be relatively straightforward to translate them so they can exploit this new version of the code.  Look 
at what has been done to the other multicolvars in there for guidance.  Sorry for any inconvenience caused.
- Changed the logic of several PLUMED ifdef macros so as to make them consistent.
  Now every feature based on external libraries is identified by a __PLUMED_HAS_* macro.

Version 2.2.1 (Jan 18, 2016)
---------------------------------------------

For users:
- \ref PBMETAD implement the new Parallel Bias Metadynamics flavor of the Metadynamics sampling method.
- PLUMED now reports an error when using \ref HISTOGRAM with UNNORMALIZED without USE_ALL_DATA. See \issue{175}
- Fixed a bug in configure together with --enable-almost. The check for lbz2 library was not working properly.
- Fixed a bug in install procedure that was introducing an error in linking with CP2K.
- Fixed a bug that sometimes was preventing the printing of a useful error message.

For developers:
- Vector and Tensor now support direct output with `<<`.
- Added some missing matmul operation Vector and Tensor.
- ./configure is automatically relaunched when changing ./configure or Makefile.conf. This makes it more robust
  to switch between branches.

Version 2.2.2 (Apr 13, 2016)
----------------------------------------------

For users:
- \ref MOLINFO for RNA accepts more residue names, see \issue{180}.
- added two mpi barries (one was missing in PBMetaD for multiple walkers) to help synchronized initialisation
- Fixed a bug in internal stopwatches that was making \ref DEBUG logRequestedAtoms not working
- Some multicolvars (including \ref BRIDGE, \ref ANGLES, and \ref INPLANEDISTANCES) now crashes if one
  asks for too many atoms, see \issue{185}.
- Optimisations (activation of the dependencies, secondary structures, DRMSD)
- Fixed a performance regression with RMSD=OPTIMAL-FAST
- Fixed a bug in the normalization of kernel functions (relevant for \ref HISTOGRAM).
- Fixed a regression introduced in v2.2 that was making \ref METAD with non-MPI multiple walkers crash
  if reading frequently. See \issue{190}
- Updated patch for gromacs 5.x. Patches for gromacs 5.0 and 5.1 have been fixed so as to allow
  patching in runtime mode.
- Possibility to control manual generation (including pdf) from ./configure. Pdf manual is now off
  by default. Notice that on travis CI it is still generated.

For developers:
- Fixed a bug in the interpretation of cmd strings. Namely, an erroneous string was not triggering an error.
  This is harmless for MD codes properly patched, but could have introduced problems in MD codes with typoes
  in cmd strings.
- ./configure is not automatically relaunched anymore when doing `make clean`.

Version 2.2.3 (Jun 30, 2016)
----------------------------------------------

For users:
- Updated patches for gromacs 5.1.x and 5.0.x to fix a problem when plumed was trying to write to an already
  closed gromacs log file.
- When looking for a value outside the GRID now the error include the name of the responsible 
  collective variable
- Numerical check in LatticeReduction made less picky. This should solve some of the internal errors reported
  by `LatticeReduction.cpp` when using aggressive compilers.
- Files are now flushed at the correct step. Before this fix, they were flushed at the step before the requested one
  (e.g. with \ref FLUSH STRIDE=100 at step 99, 199, etc).
- In \ref METAD, INTERVAL with periodic variables now report an error.
- \ref LOAD now works also when plumed is installed with a suffix.
- Added `--md-root` option to `plumed patch` which allows it to be run from a directory different from the one
  where the md code is located.
- Wham script in \ref munster tutorial now writes weights in scientific notation.

For developers:
- `./configure` checks if dependencies can be generated. If not, they are disabled.
- Added --disable-dependency-tracking to ./configure
- Added a make target `all_plus_doc` that builds both code and docs.
- Added possibility to set a default location for plumed library in runtime binding.
  If the plumed wrapped is compiled with `-D__PLUMED_DEFAULT_KERNEL=/path/libplumedKernel.so`,
  then if the env var PLUMED_KERNEL is undefined or empty PLUMED will look in the path at compile time.
- Tentative port files are now available at [this link](http://github.com/plumed/ports). 
  They can be used to install PLUMED using MacPorts.

Version 2.2.4 (Dec 12, 2016)
-------------

For users:
- Fix a bug in \ref PBMETAD when biasing periodic and not periodic collective variables at the same time 
- GSL library is now treated by `./configure` in the same way as other libraries, that is `-lgsl -lgslcblas` are only
  added if necessary.
- Fix a bug in \ref METAD when using INTERVAL and ADAPTIVE gaussians at the same time
- Updated gromacs patch for 5.1.x to 5.1.4
- Fix a performance regression in the calculate loop where derivatives and forces were set to zero even if an action
  was not active, this is relevant for postprocessing and for the on-the-fly analysis
- Torsion calculation has been made slightly faster and improved so as to provide correct
  derivatives even for special angles (e.g. +pi/2 and -pi/2).

For developers:
- Macports portile is now tested on travis at every plumed push.

Version 2.2.5 (Mar 31, 2017)
-------------

\plumednotmaintained

For users:
- Fixed a problem with large step numbers in driver (see \issue{209}).
- Fixed a problem leading to crashes when using switching functions without cutoff with some compiler (see \issue{210}).
- Fixed a bug when using \ref FIT_TO_TEMPLATE and domain decomposition (see \issue{214}).
- Added an automatic flush of HILLS files when using \ref METAD with file-based multiple walkers.
- Root dir is logged to allow easier debugging of problems.

@page CHANGES-2-4 Version 2.4

## Version 2.4 (Dec 15, 2017)

Version 2.4 contains several improvements with respect to 2.3. Users currently working with 2.3
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files.
Notice that version 2.4 includes already all the fixes in branch 2.3 up to 2.3.3 indicated in \ref CHANGES-2-3 .

Changes from version 2.3 which are relevant for users:
- Changes leading to incompatible behavior:
  - A c++11 compliant compiler is required (see \issue{212}). This should mean:
    - gcc 4.8
    - clang 3.3
    - intel 15
    Since the number of c++11 features that we use is limited, older compilers might work as well.
  - The meaning of `BIASFACTOR=1` in \ref METAD has been modified and can now be used to indicate unbiased
    simulations. Non-well-tempered metadynamics is BIASFACTOR=-1, which is the new default value.
    Notice that this has an implication on the bias factor written in the HILLS file when doing
    non-well-tempered metadynamics.
  - Due to a change in \ref COMMITTOR, the format of its output file has been slightly changed.
  - \ref HISTOGRAM : When using weights default is now to output histogram divided by number of frames from which data was taken.  In addition the 
    UNORMALIZED flag has been replaced with the keyword `NORMALIZATION`, which can be set equal to true, false or ndata.
  - All switching functions are now stretched by default, also when using the "simple syntax" (e.g. `COORDINATION NN=6`).
    Switching functions were already stretched by default when using the advanced syntax (e.g. `COORDINATION SWITCH={}`)
    since version 2.2.  Notice that this will introduce small numerical differences in the computed switching functions.
- New modules:
  - A new PLUMED-ISDB module have been included, this module includes a number of CVs to calculate experimental data with the internal ability
    to also calculate a \ref METAINFERENCE score.
    - New actions include:
      - \ref EMMI
      - \ref SAXS
      - \ref RESCALE, \ref SELECT, \ref SELECTOR
    - Updated actions include:
      - \ref CS2BACKBONE
      - \ref FRET
      - \ref JCOUPLING
      - \ref METAINFERENCE
      - \ref NOE
      - \ref PRE
      - \ref RDC, \ref PCS
      - \ref PBMETAD
  - A new EDS module have been included, contributed by Glen Hocky and Andrew White.
    This module implements the following methods:
    - \ref EDS
  - A new DRR module have been included, contributed by Haochuan Chen and Haohao Fu.
    This module implements the following methods:
    - \ref DRR
    - \ref drr_tool
  - A new VES module have been included, contributed by Omar Valsson.
    This module implements the following methods:
    - \ref BF_CHEBYSHEV
    - \ref BF_COMBINED
    - \ref BF_COSINE
    - \ref BF_CUSTOM
    - \ref BF_FOURIER
    - \ref BF_LEGENDRE
    - \ref BF_POWERS
    - \ref BF_SINE
    - \ref OPT_AVERAGED_SGD
    - \ref OPT_DUMMY
    - \ref TD_CHI
    - \ref TD_CHISQUARED
    - \ref TD_CUSTOM
    - \ref TD_EXPONENTIAL
    - \ref TD_EXPONENTIALLY_MODIFIED_GAUSSIAN
    - \ref TD_GAUSSIAN
    - \ref TD_GENERALIZED_EXTREME_VALUE
    - \ref TD_GENERALIZED_NORMAL
    - \ref TD_GRID
    - \ref TD_LINEAR_COMBINATION
    - \ref TD_PRODUCT_COMBINATION
    - \ref TD_PRODUCT_DISTRIBUTION
    - \ref TD_UNIFORM
    - \ref TD_VONMISES
    - \ref TD_WELLTEMPERED
    - \ref VES_LINEAR_EXPANSION
    - \ref VES_OUTPUT_BASISFUNCTIONS
    - \ref VES_OUTPUT_FES
    - \ref VES_OUTPUT_TARGET_DISTRIBUTION
    - \ref ves_md_linearexpansion
- New collective variables:
  - \ref DIMER (thanks to Marco Nava).
  - \ref EEFSOLV : EEF1 implicit solvent solvation energy
  - \ref ADAPTIVE_PATH : Adaptive path variables using the method from \cite BerndAdaptivePath
- New actions:
  - \ref INENVELOPE
  - \ref TOPOLOGY_MATRIX
  - \ref BOND_DIRECTIONS
  - \ref DUMPGRAPH
  - \ref GRID_TO_XYZ
  - \ref INTEGRATE_GRID
  - \ref LWALLS
  - \ref MAXENT
  - \ref MCOLV_COMBINE
  - \ref MCOLV_PRODUCT
  - \ref POLYMER_ANGLES
  - \ref XANGLES , \ref YANGLES , \ref ZANGLES
  - \ref XYTORSIONS , \ref XZTORSIONS , \ref YXTORSIONS , \ref YZTORSIONS , \ref ZXTORSIONS , and \ref ZYTORSIONS
- New command line tools:
  - \ref pesmd : Tool for performing Langevin dynamics on an energy landscape that is specified using a PLUMED input file
  - \ref pathtools 
- Other changes:
  - Sharing coordinates and applying force is now faster (in some cases these can result in much better scaling of the performances in parallel).
  - \ref COMMITTOR : new flag to use committor to keep track of the visited basins without stopping the simulation
  - \ref PBMETAD : multiple walkers using files (thanks to Marco De La Pierre).
  - \ref PBMETAD : adaptive Gaussian kernels
  - \ref PBMETAD : default names for `GRID` and `FILE` (useful with many collective variables) 
  - \ref METAD : BIASFACTOR=1 is allowed and performs unbiased sampling. HILLS file can be used
    to recover free energy also in this case.
  - \ref METAD : a RECT option is available that allows setting an array of bias factors, one for each replica.
  - \ref METAD : added options to perform Transition Tempered Metadynamics (thanks to James Dama)
  - \ref PATHMSD and \ref PROPERTYMAP now support alignment to a close structure (thanks to Jana Pazurikova)
  - PDB files with more than 100k atoms can now be read using [hybrid 36](http://cci.lbl.gov/hybrid_36/) format,
    see \issue{226}.
  - Added lepton support. Set env var `export PLUMED_USE_LEPTON=yes` to activate lepton as a matheval replacement
    in \ref MATHEVAL, \ref CUSTOM, and \ref switchingfunction "MATHEVAL switching function".
    Notice that in v2.5 matheval support will be dropped and all these keywords will use lepton.
    See \issue{244}.
  - When parsing constants, PLUMED uses lepton library. This allows to pass
    arguments such as `HEIGHT=exp(0.5)` (see \ref parsing-constants).
  - \ref CUSTOM function has been added as an alias to \ref MATHEVAL .
  - Trajectories read in \ref driver also support the usual replica convention, that is if
    trajectory with replica suffix is not found the driver will look for a trajectory without the replica suffix.
  - A new syntax (`@replicas:`) can be used to specify different arguments for different replicas (see \ref special-replica-syntax).
  - Internal molfile implementation has been updated to VMD 1.9.3.
  - Examples in the documentation now have syntax highlighting and links to the documentation of used actions.
  - \ref COORDINATIONNUMBER : Added option to have pairwise distance moments of coordination number in the multicolvar module
  - GROMACS patch updated to gromacs-2016.4
  - Implemented HREX for gromacs-2016.4.
  - Added patch for Quantum ESPRESSO 6.2 (thanks to Ralf Meyer).
  - Fixed a bug in \ref LOCAL_AVERAGE which appears when you use `SPECIESA` and `SPECIESB` keywords instead of just `SPECIES`
  - Added possibility to pass `--kt` from \ref driver.

Changes from version 2.3 which are relevant for developers:
  - A few fixes has been made to improve exception safety. Although we still cannot declare
    PLUMED totally exception safe (there are still many non-safe pointers around),
    this made it possible to add a regtest that actually tests erroneous cmd strings
    and erroneous inputs.
  - Due to the required c++11 support, travis-ci test on Ubuntu Precise has been removed.
  - `gettimeofdate` and `gettime` have been replaced with portable `chrono` classes introduced in c++11.
  - C++ exceptions are enabled by default.
  - A large number of loops have been changed to use the `auto` keyword in order to improve code readability.
  - Stack trace is not written upon error anymore, unless environment variable `PLUMED_STACK_TRACE` is set at runtime.
  - Fixed a potential bug using single precision system BLAS on a mac (notice that currently plumed only uses
    double precision, so it is harmless).
  - Added `--enable-rpath` option for autoconf (off by default).
  - Files related to changelog are now stored as `.md` files. This makes
    it possible to navigate them from github.
  - `configure.ac` has been simplified and improved in order to more easily probe C++ libraries.
  - added `plumed_custom_skip` function to regtests in order to skip specific tests based on specific conditions (e.g. OS).
  - environment variable `LDSO` has been renamed to `LDSHARED`, which is standard in the python community.
  - a `libplumedWrapper.a` library is installed as well, that is used in `--runtime` patching.
  - pkgconfig files are installed.
  - `plumed config makefile_conf` can be used to retrieve `Makefile.conf` file a posteriori.
  - Store `MPIEXEC` variable at configure time and use it later for running regtests. Notice that in case
    `MPIEXEC` is not specified regtests will be run using the command stored in env var `PLUMED_MPIRUN` or, if this is
    also not defined, using `mpirun`.
  - Added canonical Makefile targets `check` and `installcheck`. Notice that `check` runs checks with
    non-installed plumed whereas `installcheck` uses the installed one, including its correct program name if it
    was personalized (e.g. with suffixes). Notice that this modifies the previously available `check` target.
    
    
## Version 2.4.1 (Mar 2, 2018)

For users:
  - Fixed an important bug affecting RMSD calculations with compilers supporting OpenMP 4 (e.g.: intel compiler). Notice that this bug might potentially affect not only
    \ref RMSD variable, but also \ref PATHMSD variables using RMSD, \ref FIT_TO_TEMPLATE, \ref PCAVARS, and possibly other variables based on RMSD calculations and optimal alignments
    (see \issue{343}). Results might depend on the exact architecture and on how aggressive is the compiler. The bug is a consequence of some erroneous SIMD directives introduced in 2.4.0, so it does not affect PLUMED 2.3.x. 
  - Resolved a problem with \ref CS2BACKBONE and glycine atom names.
  - Module VES: Fixed a bug with basis functions that have a constant function different from 1 (e.g. scaled version of the Legendre basis functions, \ref BF_LEGENDRE) that was causing a time-dependent shift in the bias potential.
  - Module VES: In optimizers (\ref OPT_AVERAGED_SGD and \ref OPT_DUMMY) the output of quantities related to the instantaneous gradients are now off by default as these quantities are generally not useful for normal users, their output can instead by re-enabled by using the `MONITOR_INSTANTANEOUS_GRADIENT` keyword. Also added an keyword `MONITOR_AVERAGE_GRADIENT` that allows to monitor the averaged gradient and output quantities related to it. 
  - \ref RMSD variable and other collective variables using reference PDB files now crash when zero weights are passed (see \issue{247}).
  - Using \ref COM with \ref driver without passing masses now triggers an error instead of reporting NaNs (see \issue{251}).

For developers:
  - `plumed patch -p` command can be used twice without triggering an error. This will allow e.g. building again
    on MacPorts in cases where the build was interrupted. Notice that this only works for patches without special
    after/before patch/revert functions.

## Version 2.4.2 (Jul 2, 2018)

For users:
  - All fixes done in version 2.3.6. Notice that \issue{363} in version 2.4 also applies to \ref pathtools.
  - Additional residue names (without the prefix `D`) are now supported by \ref MOLINFO for DNA. See \issue{367}.
  - Solved an important bug appearing in NAMD interface. Notice that the bug was a regression introduced in 2.4.0. As consequence, versions <= 2.3 and versions >=2.4.2
    are expected to work correctly. See \issue{254}.
  - GROMACS patch for gromacs-2018.1.
  - \ref VimSyntax now highlights `__FILL__` strings.
  - \ref METAD and \ref PBMETAD give a warning when one restarts a simulation and the old hills file is not found. See \issue{366}.

For developers:
  - `LDSHARED` is now correctly taken into account when launching `./configure`.
  - Fixed installation with `--disable-shared`.
  - Cppcheck upgraded to 1.84.

## Version 2.4.3 (Oct 5, 2018)

For users:
  - All fixes done in version 2.3.7.
  - Module VES: Fixed a bug in `TD_GRID` for 2D grids where the grid spacing is not the same for both dimensions.
  - GROMACS patch for gromacs-2018.3.

## Version 2.4.4 (Dec 19, 2018)

For users:
  - Fixed some performances regression issue with OpenMP
  - Updated NAMD patches to version 2.12 and 2.13. Old patches have been removed.
  - GROMACS patch for gromacs-2018.4.
  - Fixed a thread safety issue using forces on \ref HISTOGRAM 
  - Fixed error message suggesting wrong actions (see \issue{421}).

For developers:
  - All fixed done in version 2.3.8
  - Cppcheck updated to 1.85

## Version 2.4.5 (Apr 1, 2019)

For users:
  - Fixed an inconsistency in parsing of braces.
    It is now possible to pass individual options
    including spaces (e.g. with `FILE={/path with space/file}`). Notice 
    that this invalidates syntax such as `ATOMS={1}{2}{3}{4}`. See more
    at \issue{434}.
  - Fixed \ref simplemd so as to call "runFinalJobs" at the end of the simulation.
  - GROMACS patch for gromacs-2016.6.
  - GROMACS patch for gromacs-2018.6.
  - Added aliases for some actions/options containing dashes (`-`) in their name. This will improve
    backward compatibility when these actions/options will be removed (see \issue{449}).

## Version 2.4.6 (Jul 19, 2019)

For users:
  - Fixed a bug in \ref COORDINATIONNUMBER where derivatives were wrong when using `R_POWER` > 2, thanks to `@MoleOrbitalHybridAnalyst` for spotting and fixing
  - Fixed a bug in library search, possibly affecting linked blas/lapack on OSX (see \issue{476}).
  - Fixed a bug in \ref METAD with `TARGET` and `GRID_SPARSE` (see \issue{467}).

## Version 2.4.7 (Jan 27, 2020)

For users:
  - Fixed a bug with \ref CONVERT_TO_FES and periodic variables, see \issue{441} (backported from v2.5.3).
  - More robust backup for output files when running over multiple processes
  - Fixed a regression in the performances of `GEOMETRY` based flexible hills in \ref METAD and \ref PBMETAD
  - Fixed \issue{538}.
  - Fixed potential issue with VMD plugins from 1.9.4 (\issue{545}, thanks to Lixin Sun).
  - Module VES: Fixed an off-by-one bug in the output of target distribution averages. The bug only affects the output and does not affect results. The bug also affected the output of coefficients when using a bias cutoff. 
  - Module VES: Made sure that all relevant output files are written out at the final step when shutting down the simulation. This solves issues reported by @PabloPiaggi with restarting when there is a mismatch been the output of files and the number of MD steps. 

## Version 2.4.8 (Jul 8, 2020)

\plumednotmaintained

For users:
   - Take into account `UNITS` when using MD codes that feeds one line at a time to PLUMED (e.g., OpenMM). See \issue{582}.
   - Fix PDB parser for non justified atom numbers. See \issue{592}.
   - Fix in \ref INPLANEDISTANCES. See \issue{595}.

For developers:
   - Tests and doc building moved from Travis-CI to [GitHub Actions](https://github.com/plumed/plumed2/actions) (see \issue{634}).

@page CHANGES-2-5 Version 2.5

## Version 2.5 (Dec 19, 2018)

This page contains changes that will end up in 2.5

Changes from version 2.4 which are relevant for users:
- Changes leading to incompatible behavior:
  - \ref RMSD, \ref MULTI-RMSD, \ref PATHMSD, \ref PROPERTYMAP, \ref PCAVARS, \ref PCARMSD, \ref FIT_TO_TEMPLATE,
    \ref DIPOLE, \ref ALPHARMSD, \ref ANTIBETARMSD, and \ref PARABETARMSD now automatically make molecules whole.
    In case you do not want them to do it, use NOPBC flag,
  - There is some subtle change in the installation layout (see below). There should be no visible effect, however it is now compulsory
    to set correctly the `LD_LIBRARY_PATH` variable for the linux executable to work correctly. The procedure has been tested well on OSX and Linux,
    but could give problems on other platform. Please report possible problems on the mailing list.
  - \ref driver now stops correctly when using \ref COMMITTOR. If you want to continue the analysis, use the `NOSTOP` flag in \ref COMMITTOR.
  - \ref METAD the calculation of the reweighting factor is now activated by CALC_RCT instead of REWEIGHTING_NGRID and REWEIGHTING_NHILLS, the frequency of update can be set 
    by RCT_USTRIDE, the default value is 1 and should be OK for most of the cases
  - Fixed sign in Cartesian components of \ref PUCKERING with 6 membered rings (thanks to Carol Simoes and Javi Iglesias).

- New actions:
  - \ref COLLECT_FRAMES
  - \ref EUCLIDEAN_DISSIMILARITIES
  - \ref HBPAMM_MATRIX
  - \ref HBPAMM_SH
  - \ref LANDMARK_SELECT_FPS
  - \ref LANDMARK_SELECT_RANDOM
  - \ref LANDMARK_SELECT_STAGED
  - \ref LANDMARK_SELECT_STRIDE
  - \ref OUTPUT_ANALYSIS_DATA_TO_COLVAR
  - \ref OUTPUT_ANALYSIS_DATA_TO_PDB
  - \ref OUTPUT_PCA_PROJECTION
  - \ref PAMM
  - \ref PLUMED
  - \ref PRINT_DISSIMILARITY_MATRIX
  - \ref PROJECT_ALL_ANALYSIS_DATA
  - \ref READ_DISSIMILARITY_MATRIX
  - \ref RESELECT_LANDMARKS
  - \ref REWEIGHT_WHAM
  - \ref SKETCHMAP_CONJGRAD
  - \ref SKETCHMAP_POINTWISE
  - \ref SKETCHMAP_READ
  - \ref SKETCHMAP_SMACOF
  - \ref SKETCH_MAP
  - \ref SMACOF_MDS
  - \ref WHAM_HISTOGRAM
  - \ref WHAM_WEIGHTS

- New command line tools:
  - \ref completion (used to generate command line completion scripts).
  - \ref pdbrenumber (see \issue{371}).

- New modules:
  - A new PIV module has been included, contributed by Silvio Pipolo and Fabio Pietrucci.
    This module implements the following collective variable:
    - \ref PIV
  - A new LOGMFD module has been included, contributed by Tetsuya Morishita.
    This module implements the following bias:
    - \ref LOGMFD

- Changes in the ISDB module
  - \ref CS2BACKBONE is now mpi parallelized in particular with DOSCORE and CAMSHIFT
  - \ref SAXS has an additional implementation based on Bessel functions that can be faster for large systems (new keyword BESSEL)
  - \ref SAXS keyword SCEXP has been renamed into SCALEINT
  - \ref SAXS includes the MARTINI bead structure factors for Proteins and Nucleic Acids
  - \ref SAXS includes a GPU implementation based on ArrayFire (need to be linked at compile time) that can be activated with GPU
  - \ref METAINFERENCE and all related methods has a new keyword REGRES_ZERO to scale data using a linear scale fit
  - \ref CALIBER new bias to perform Maximum Caliber replica-averaged restrained simulations 

- Changes in the eABF/DRR module (contributed by Haochuan Chen and Haohao Fu):
  - \ref DRR now supports the extended generalized ABF(egABF) method.
  - \ref DRR accepts different GRID options for CVs and extended variables.
  - The MAXFACTOR option is added in \ref DRR to control the factor of biasing force.
  - \ref drr_tool can calculate the divergence of gradients now. (Maybe useful for future pABF)
  - Fixed conflicts of output files in multiple replicas.

- Changes in the EDS module:
  - \ref EDS implements Levenberg-Marquardt optimization in addition to previous gradient descent. 
  - \ref EDS no longer automatically increases prefactor for bias parameter updates. This results in more stable optimization for the cases tested.
  - \ref EDS now has a larger default RANGE parameter to go with these other changes.

- Other changes:
  - \ref METAD there is a new FLYING_GAUSSIAN keyword to activate the flying gaussian methods by Spiwok (contributed by Spiwok and Hozzova)
  - \ref EXTERNAL can now SCALE the input grid. This allows for more flexibility without modifying the grid file.
  - \ref ALPHABETA can now combine dihedral angles with different coefficients
  - \ref INCLUDE can now be used also before setup actions.
  - \ref CENTER can now be computed using trigonometric functions (PHASES) to simplify its calculation with periodic boundary conditions.
  - Libmatheval is not used anymore. \ref MATHEVAL (and \ref CUSTOM) are still available
    but employ an internal implementation of the lepton library.
    Functions available in libmatheval and absent in the original lepton library have been added so as to have backward compatibility.
    `atan2(y,x)` function has also been added.
    Notice that MATHEVAL (and CUSTOM) \ref switchingfunction "switching functions"
    using the lepton library have been further optimized with respect to PLUMED 2.4.
    Finally, notice that it is possible to use asmjit to optimize performance (see \ref Lepton).
  - Implemented bash autocompletion, see \ref BashAutocompletion.
  - \ref MOLINFO now allows selecting atoms from chains with a numeric ID (see \issue{320}).
  - Removed the patch for GMX 5.1.4
  - LAMMPS patch has been finally removed. Notice that LAMMPS has native support for PLUMED now.
  - AMBER patch has been finally removed. Notice that AMBER (sander module) has native support for PLUMED starting from version 15.
  - \ref RMSD calculation has been optimized. This should positively affect the performances of CVs where
     many RMSD values are computed on small groups of atoms, such as secondary structure variables.
  - In \ref METAD, when using a bias factor equal to one (no bias) the `rct` component is set to zero rather than to one.
  - New shortcuts are available for selecting atoms: `@allatoms` and `@mdatoms` (see \ref atomSpecs).
  - When using \ref MOLINFO, also the following shortcuts are available for selecting atoms: `@nucleic`, `@protein`, `@water`, `@ions`, `@hydrogens`, `@nonhydrogens`.
  - When using \ref MOLINFO, individual atoms can be chosen also from water molecules (e.g. `@OW-100`).
  - Additional switching function COSINUS contributed by Michael King
  - added API to set the number of used openMP threads from the linked code, updated gromacs 2018.3 patch to use it

Changes from version 2.4 which are relevant for developers:
- Code has been cleanup up replacing a number of pointers with `std::unique_ptr`. All `delete` statements
  in the core parts of the code have been eliminated.
- Exceptions cannot be disabled (`--disable-cxx-exceptions` option has been removed from `./configure`).
- Every exception thrown in PLUMED now also writes its message on PLUMED log.
- Runtime loader in `Plumed.c` now works also when linked without `-rdynamic` (that is, 
  its names are not exported). Notice that all the combinations are expected to
  work, that is: `Plumed.c` from <=2.4 or >=2.5 combined with libplumedKernel
  from <=2.4 or >=2.5. In order to achieve this the following changes are implemented:
  - libplumedKernel does not depend anymore on `Plumed.c`. This allows loading it even
    in cases where names in the loader are not visible. The relevant function needed
    to be compatible with `Plumed.c` <=2.4 are found using `dlsym`.
  - `Plumed.c` does not need anymore libplumedKernel to register itself, but rather
    searches the relevant functions using `dlsym`. In addition, if it is not able to
    load `libplumedKernel` since the latter is <=2.4 and needs `Plumed.c` to be visible,
    it just uses as a fallback `libplumed`, which should load properly.
- In addition to the capability mentioned above, the MD-code interface has been significantly
  improved and allows for:
  - Translation of exception (allowing to mix PLUMED and an MD-code linked against a different C++ library).
  - Possibility to choose the path to the PLUMED kernel while instantiating a Plumed object.
  See the developer documentation for more information.
- The installation layout of shared libraries has been modified. In particular,
  both `libplumed.so` and `plumed` links to `libplumedKernel.so`.
  This reduces considerably the size of the installed package. In addition, it allows
  using two-level namespace on OSX. Notice that this implies that on Linux one should
  always set the `LD_LIBRARY_PATH` flag to have a working executable.
- A smaller number of header files is installed. In particular, all the files that were historically generated in subdirectories
  (such as `plumed/core/tools/Vector.h', just including `plumed/tools/Vector.h`) are not installed and the related include
  statements are fixed. This makes the installed package smaller.
- List of preferred compilers (used when `CXX` or `CC` are not set) has been changed. On OSX, `./configure` will try `clang++/clang` as first choices.
- Added `--enable-static-archive` to `./configure` to build a `libplumed.a` static library (yes by default).
- Stop setting `DYLD_LIBRARY_PATH` in `sourceme.sh` and in modulefile. Notice that as of PLUMED v2.3.3
  it should not be needed.
- Coverage scan is not anymore contained in developer manual. It can be found in a separate repository
  `github.com/coverage-branch` (see \issue{348}). In addition, coverage for third-party libraries included in PLUMED
  is reported as well.
- It is not possible anymore to use `make install prefix=/path`. Prefix can only be changed during `./configure` (see \issue{332}).
- Exception class has been rewritten to allow more extensive messages. Now also function name is shown.
- On linux, library is linked with `-Bsymbolic`.
- When launching `plumed`, flags `--no-mpi` and `--mpi` can appear multiple times. The last appearance is the effective one.
- Internal BLAS and LAPACK libraries updated to gromacs 2018.
- Choosing `./configure --prefix=$PWD` does not lead anymore to deletion of all header files.
- A copy of `plumed-runtime` is installed in `prefix/lib/plumed` and can be used for testing.
- Absolute/relative soname/install_name can be configured on linux/OSX. This feature is only
  for testing, the default choice is the typical one used on the respective operating system.
- On OSX, `plumed` and `libplumed.dylib` will find `libplumedKernel.dylib` using `@loader_path`.
- Using CXX compiler to link the main program.
- plumed can be compiled with ArrayFire to enable for gpu code. \ref SAXS collective variable is available as part of the isdb module to provide an example of a gpu implementation for a CV


## Version 2.5.1 (Apr 1, 2019)

For users:
- in \ref SAXS the keyword ADDEXP is removed. Furthemore, SAXS intensities are automatically normalised for I(0)=1, in case experimental data are provided, the intensity is rescaled with the intensity of the lowest q provided. As a consequence SCALEINT is only needed for additional adjustments.
- gromacs patch updated to gromacs 2018.5
- Fixed a bug in gromacs patch that was resulting in incorrect number of threads (0) set when not explicitly using `-ntomp` on the 
  command line or setting `OMP_NUM_THREADS` (see \issue{446}). To apply this fix you need to re-patch gromacs.
  Notice that setting the number of threads to zero might lead to inconsistent results when using secondary structure variables
  or other multicolvars.
- Fixed PLUMED so that when zero threads are selected from gromacs (see previous fix) the number of used threads is set to 1.
  This fix allows to use a GROMACS executable patched with PLUMED 2.5.0 and linked at runtime with PLUMED 2.5.1 without introducing
  errors. However, re-patching is preferred since it selectes the correct number of threads.
- Python wrappers:
  - Fixed building of python interface on MacOS Mojave (see \issue{445}, thanks to Omar Valsson).
  - Numpy is not required anymore at build time (though it is required at runtime for our tests).
  - Raw python arrays can be passed as an alternative to Numpy ndarrays.

## Version 2.5.2 (Jul 19, 2019)

For users:
- New shortcuts are available for selecting protein atoms: `@chi2-#`, `@chi3-#`,`@chi4-#` and `@chi5-#`
- Fixed performance of \ref CUSTOM when having zero derivatives with respect to some arguments.
- New --parse-only option in \ref driver to check the validity of a plumed input file
- New patch for GROMACS 2019.2
- Module VES: Fixed performance of \ref BF_CUSTOM for basis functions with linear terms (e.g. having zero derivatives). 
- Python wrappers:
  - Python module is now always named `plumed` irrespectively of program prefix and suffix. Notice 
    that python module is installed inside the `lib/program_name` directory and thus it is not necessary to
    use `program_name` in order to install multiple modules side by side.
  - Python module can be compiled without compiling PLUMED first.
  - `Plumed` object can be explicitly finalized using `finalize()`. Can be used to make sure all files are closed,
    but it is not necessary if the `Plumed` object gets correctly collected by Python.
  - `Plumed` object can be used in context managers (e.g. `with plumed.Plumed() as p:`).
- Precompiled binaries are available on Anaconda cloud on the [conda-forge channel](https://anaconda.org/conda-forge/plumed).

## Version 2.5.3 (Oct 11, 2019)

For users:
- Fixed a bug with \ref CONVERT_TO_FES and periodic variables, see \issue{441}
- Fixed a bug with \ref FOURIER_TRANSFORM 
- Updated patch for GROMACS 2019.4
- Updated patch for GROMACS 2018.8
- Python module:
  - Fixed building with clang-8.
  - Set `language_level` for cython to the actually used language level.
  - Force using cython when compiling from source. Still using the pre-generated cpp file
    when installing from PyPI, to avoid cython dependency.
  - Using python 2 to create the cpp file uploaded on PyPI (this will change to python 3 in 2.6, see \issue{502}).
- Module VES: Fixed a bug in updating of bias potential in \ref VES_LINEAR_EXPANSION that is present for certain integrators that call the calculation of the bias multiple times (see [here](https://groups.google.com/d/msg/plumed-users/kPZu_tNZtgk/LrkS0EqrCQAJ)) and replica exchange.

## Version 2.5.4 (Jan 27, 2020)

For users:
- Includes all fixes up to 2.4.7

## Version 2.5.5 (Jul 8, 2020)

For users:
- Includes all fixes up to 2.4.8

For developers:
- Small fix to avoid unique global symbols (see \issue{549})

## Version 2.5.6 (Oct 26, 2020)

For users:
- Report an error when using all weights set to zero in reference PDB files. Same as \issue{247} fixed in 2.4, but now the check is done in more cases.
- Fixed overflow in reweighting factor when using non well-tempered metadynamics (thanks to Michele Invernizzi, see \issue{599}).
- Fixed \ref READ with EVERY when reading a trajectory file (see \issue{619}).

For developers:
- Fixed a warning in `wrapper/Plumed.h` appearing with recent clang versions.

## Version 2.5.7 (Apr 16, 2021)

\plumednotmaintained

For users:
- Fixed handling of periodic variables in \ref ves_md_linearexpansion (see \issue{649}).
- Small fix that might affect performance (backport of a fix needed for master branch, see \issue{680}).

@page CHANGES-2-9 Version 2.9
  
## Version 2.9 (under development)

This page contains changes that will end up in 2.9

@page CHANGES-2-1 Version 2.1

Version 2.1.0 (September 15, 2014)
----------------------------

Version 2.1 contains several improvements with respect to 2.0. Users currently working with 2.0 
should have a look at the section "Changes leading to incompatible behavior" below and
might need tiny adjustments in their input files. In 2.1 we restored more features of 1.3
that were missing in 2.0, so users still working with 1.3 could opt for an upgrade.
A tutorial explaining how to move from PLUMED 1 to PLUMED 2 is available (see \ref moving).

Below you find a list of all the changes with respect to version 2.0.
Notice that version 2.1 includes already all the fixes in branch 2.0 up to 2.0.4.

Changes from version 2.0 which are relevant for users:
- Changes leading to incompatible behavior:
  - \ref COORDINATION now skips pairs of one atom with itself.
  - Labels of quantities calculated by \ref BIASVALUE have changed from <i>label</i>.bias.<i>argname</i> to <i>label</i>.<i>argname</i>_bias, which is more consistent with steered MD
  - Labels of quantities calculated by \ref ABMD have change from <i>label</i>.min_<i>argname</i> to <i>label</i>.<i>argname</i>_min, which is more consistent with steered MD
  - Labels of quantities calculated by \ref PIECEWISE have change from <i>label</i>.<i>argnumber</i> to <i>label</i>.<i>argname</i>_pfunc, which is more consistent with steered MD
  - For multicolvars components calculated with LESS_THAN and MORE_THAN keywords are now labelled lessthan and morethan. This change is necessary as the underscore
character now has a special usage in component names.
  - In \ref CONTACTMAP components are now labelled <i>label</i>.contact-\f$n\f$.
  - The command SPHERE has been replaced by \ref UWALLS.
- New configuration system based on autoconf (use ./configure from root directory).
  Optional packages are detected at compile time and correctly
  enabled or disabled. An internal version of LAPACK and BLAS will be used
  if these libraries are not installed.
- New actions:
  - \ref SPRINT topological collective variables.
  - CH3SHIFTS collective variable.
  - \ref POSITION collective variable.
  - \ref FIT_TO_TEMPLATE.
  - \ref COMMITTOR analysis.
  - \ref LOCAL_AVERAGE.
  - \ref NLINKS.
  - \ref DIHCOR.
  - \ref NOE.
  - \ref RDC. 
  - \ref CLASSICAL_MDS.
  - \ref XDISTANCES.
  - \ref YDISTANCES.
  - \ref ZDISTANCES.
  - \ref DUMPMULTICOLVAR.
  - Crystallization module, including \ref Q3, \ref LOCAL_Q3, \ref Q4, \ref Q6, \ref LOCAL_Q4, \ref LOCAL_Q6, \ref MOLECULES, \ref SIMPLECUBIC, \ref TETRAHEDRAL and \ref FCCUBIC.
  - \ref ENSEMBLE to perform Replica-Averaging on any collective variable.
- New features for existing actions:
  - \ref METAD : WALKERS_MPI flag (multiple walkers in a mpi-based multi-replica framework),
    ACCELERATION flag (calculate on the fly the Metadynamics acceleration factor),
    TAU option (alternative way to set Gaussian height in well-tempered metadynamics),
    GRID_SPACING (alternative to GRID_BIN to set grid spacing).
    Notice that now one can also omit GRID_BIN and GRID_SPACING when using
    fixed size Gaussian, and the grid spacing will be automatically set.
  - \ref DISTANCE : added SCALED_COMPONENTS
  - \ref COORDINATION : if a single group is provided, it avoids permuted atom indexes and runs
    at twice the speed.
  - \ref DUMPATOMS : PRECISION option to set number of digits in output file.
  - \ref GROUP : NDX_FILE and NDX_GROUP options to import atom lists from ndx (gromacs) files.
  - In many multicolvars, MIN and MAX options can be used.
  - \ref HISTOGRAM : GRID_SPACING (alternative to GRID_BIN to set grid spacing),
    FREE-ENERGY flags in addition to standard probability density,
    additional option for KERNEL=DISCRETE to accumulate standard histograms. 
  - \ref sum_hills : added options --spacing (alternative to --bin to set grid spacing)
    and --setmintozero to translate the minimum of the output files to zero.
  - \ref CONTACTMAP : parallelized and added weights.
- New features in MD patches (require re-patch):
  - New patch for Gromacs 5.0
  - Gromacs 4.6.X patch updated to 4.6.7
  - Gromacs 4.6.7 supports \ref COMMITTOR analysis; can be now be used to perform energy minimization;
     now passes temperature to PLUMED (this allows temperature to be omitted in some actions,
     namely \ref METAD and analysis actions).
  .
  Notice that if you use runtime binding it is not compulsory to re-patch,
  and that all combinations should work correctly
  (new/old PLUMED with re-patched/non-re-patched MD code).
- Other new features:
  - \ref driver can now read trajectories in many formats using VMD molfile plugin
    (requires VMD plugins to be compiled and installed). In case VMD plugins are not installed,
    the configuration system falls back to an internal version which implements a minimal list
    of plugins (gromacs and dcd) (kindly provided by T. Giorgino).
  - \ref switchingfunction : added STRETCH flag.
  - Negative strides in atom ranges (e.g. ATOMS=10-1:-3 is expanded to ATOMS=10,7,4,1).
  - \ref COORDINATION and \ref DHENERGY with NLIST now work correctly in replica exchange simulations.
  - Multicolvars with neighbor lists now work correctly in replica exchange simulations.
  - Improved multicolvar neighbor lists.
- Optimization:
  - Root-mean-square deviations with align weights different from displace weights
    are now considerably faster. This will affect \ref RMSD calculations plus
    other variables based on RMSD.
  - \ref WHOLEMOLECULES is slightly faster.
  - \ref COORDINATION is slightly faster when NN and MM are even and D_0=0.
  - Atom scattering with domain decomposition is slightly faster.
  - Link cells are now exploited in some multicolvars.
  - Derivatives are not calculated unless they are specifically required, because for instance you are adding
    a bias.
- Documentation:
  - All tutorial material from the recent plumed meeting in Belfast is now in the manual
  - Improvements to documentation, including lists of quantities that are output by each action that can be referenced 
  - Manual has been re-organized following suggestions received at the plumed meeting.
  - An experimental PDF version of the manual is now provided (a link can be found in the documentation homepage).

Changes from version 2.0 which are relevant for developers:
- Added regtests for plumed as a library (e.g. basic/rt-make-0). plumed command has an additional
  flag (--is-installed) to probe if running from a compilation directory or from a fully installed copy
  (this is needed for regtests to work properly).
- Improved class Communicator. Many operations can now be done directly on Vectors, Tensors, std::vector and PLMD::Matrix.
- Modified class RMSD.
- Patches for GPL codes (Quantum Espresso and Gromacs) now also include
  original code so as to simplify their modification.
- Fixed dependencies among actions such that it is now possible (and reliable)
  to use MPI calls inside Action::prepare()
- colvar/CoordinationBase.cpp has been changed to make it faster. If you devised a class which inherits from here,
  consider that CoordinationBase::pairing now needs _squared_ distance instead of distance
- It is possible to run "make install" from sub-directories (e.g. from src/colvar)
- There is a small script which disables/enables all optional modules (make mod-light/mod-heavy/mod-reset)
- Added "-q" option to plumed patch
- You can now create new metrics to measure distances from a reference configurations. If you do so such
  metrics can then be used in paths straightforwardly
- You can now use multicolvars in tandem with manyrestraints in order to add a large numbers of restraints.
- Can now do multicolvar like things in which each colvar is a vector rather than a scalar.
- Updated script that generated header files so that they properly show years. Notice that the script
  should new be run from within a git repository

This list is likely incomplete, if you are developing in PLUMED you are encouraged to follow changes on github.

Version 2.1.1 (December 15, 2014)
----------------------------------------------

This release includes all the fixes available in branch 2.0 until 2.0.5.

For users:
- New patch for AMBER 14 (sander module only). This patch should be compatible
  with any PLUMED 2 version (including 2.0). It includes most PLUMED features
  with the notable exception of multi-replica framework.
- Changed definition in arbitrary phase of eigenvectors. This will change the result of some
  analysis method where the phase does matter (e.g. \ref CLASSICAL_MDS) and make
  some regression test better reproducible.
- Fixed a portability issue in BG/P where gettimeofday is not implemented.
  Notice that this fix implies that one should execute again ./configure to have
  plumed timing working correctly.
- CS2Backbone: fixed a bug that resulted in only a fraction of the chemical shifts being printed with WRITE_CS and 
  parallel simulations (requires to get the last almost updated from SVN)
- NOE: fixed a bug in the replica-averaging
- Fixed a linking issue with ALMOST, where bz2 was always used to link ALMOST to PLUMED even if it is not compulsory 
  to build ALMOST.
- Fixed a wrong include in the GMX5 patch.
- \ref FUNCPATHMSD can now be used together with \ref CONTACTMAP to define pathways in contact-map space
- Configuration is more verbose, a warning is given if a default option cannot be enabled and an error is given if 
  an option explicitly enabled cannot be enabled.
- Compilation is less verbose (use "make VERBOSE=1" to have old behavior)
- Small fixes in documentation.

For developers:
- Tests are now performed at every single push on travis-ci.org
- Manual is built and pushed to the online server from travis-ci.org (see developer doc)
- Fixes in developer doc.

Version 2.1.2 (Mar 16, 2015)
----------------------------------------------

For users:
- Added two new short tutorials to the manual ( \ref cambridge and \ref munster ).
- Fixed a severe bug on \ref DRMSD - cutoff values were ignored by PLUMED.
  Notice that this bug was introduced in 2.1.0, so that it should not affect the 2.0.x series.
- Fixed a bug affecting LAMMPS patch used with a single processor. Notice that
  the fix is inside PLUMED, thus it does not necessarily requires re-patching.
- Sander patch now works with multiple replica (no replica exchange yet). It also contains
  some fix from J. Swails.
- GMX5 patch was not working for bias-exchange like cases
- Patching system now checks for the availability of shared/static/runtime version of plumed before
  patching
- Configure now check better if compiler flag are accepted by the compiler. This makes
  configure on bluegene more robust.
- Sourceme.sh now sets proper library path in linux also.


Version 2.1.3 (June 30, 2015)
----------------------------------------------

For users:
- Fixed bug in \ref ENSEMBLE derivatives when more than 1 argument was provided
- Fixed bug in \ref GHOST : virial is now computed correctly.
- Fixed a serious bug in virial communicated from plumed to gromacs, for both gromacs versions 4.6 and 5.0.
  See \issue{132}.
  This fix requires gromacs to be re-patched and could be very important if you run biased simulations in the NPT ensemble.
- Fixed a bug in the virial computed with \ref FIT_TO_TEMPLATE when the reference pdb had center non located at the origin.
- Fixed a bug in the the forces computed with \ref FIT_TO_TEMPLATE when used in combination with \ref COM, \ref CENTER, or \ref GHOST
- Fixed a bug that could lead plumed to be stuck with domain decomposition in some extreme case (one domain with all atoms, other domains empty).
- Fixed a bug when \ref COMBINE or \ref MATHEVAL are used with PERIODIC keyword. Now when PERIODIC keyword is used the result
  of the calculation is brought within the periodicity domain. See \issue{139}.
- Fixed a bug related to \ref RANDOM_EXCHANGES followed by \ref INCLUDE
- Fixed bug in derivatives of histogram bead with triangular kernels
- Updated gromacs patch 4.5.5 to 4.5.7
- Updated internal molfile plugins to VMD 1.9.2.
- Included crd and crdbox formats to internal molfile.
- Added --natoms to \ref driver . This is required to read coordinate
  files with VMD plugins when number of atoms is not present (e.g. amber
  crd files)
- Added the checks in the driver to detect cases where molinfo does not provide box information
  (e.g. pdb).
- Added support for readdir_r when available, which makes opening files thread safe.
- CFLAGS now include -fPIC by default
- Added a warning when using \ref METAD without grids with a large number of hills.
- Fixes in user documentation.

For developers:
- Allow external VMD plugins to be detected with --has-external-molfile. This
  is required to enable some regtest with amber files.
- Added --dump-full-virial to \ref driver
- Allow definition of variables where some of the components have derivatives and some haven't (\issue{131}).
- Improved travis tests with more debug options.
- Improved some regtest to check out-of-diagonal virial components
- Improved make cppcheck options.
- Fixes in developer documentation.

Version 2.1.4 (Oct 13, 2015)
-----------------------------

For users:
- Fixed NAMD patch. Masses and charges were not passed correctly, thus resulting in wrong
  \ref COM or \ref CENTER with MASS.
  This fix required re-patching NAMD.
  Notice that this bug was present also in v2.0 but in a different form.
  More information here (\issue{162}), including a workaround that allows masses to be fixed
  without re-patching.
- When installing with PLUMED_LIBSUFFIX an underscore is used as separator instead of a dash.
  E.g. `make install PLUMED_LIBSUFFIX=2.1` will result in an executable named `plumed_v2.1`.
  This fix a potential problem (see \ref Installation).
- Fixed erroneously reported message about MPI at the end of ./configure.
- Changed warning message about undocumented components.
- PLUMED now says in the log file if it was compiled from a dirty git repository.
- Fixed a problem leading to rare random crashes when using \ref METAD with WALKERS_MPI and multiple
  processors per replica.
- Small change in numerical accuracy of lattice reduction. Should be more
  robust when running with highly optimizing compilers.
- Fixed a bug in normalization of kernel functions.  This affects \ref HISTOGRAM
  If these actions were used with previous versions of the code care should be taken when analyzing the 
  results.
- Fixed a bug in derivatives of kernel functions with non-diagonal covariance matrices. This affects the 
  derivatives output by \ref sum_hills

Version 2.1.5 (Jan 18, 2016)
---------------------------------------------

\plumednotmaintained

For users:
- PLUMED now reports an error when using \ref HISTOGRAM with FREE-ENERGY without USE_ALL_DATA. See \issue{175}
- Fixed a bug in configure together with --enable-almost. The check for lbz2 library was not working properly.

@page CHANGES-2-0 Version 2.0

Version 2.0.0 (September 27, 2013)
----------------------------

Version 2.0 is a complete rewrite, so there is no way to write a complete set of difference
with respect to plumed 1.3. Here is a possibly incomplete summary of the difference:
- The input is simpler, more flexible, and more error proof.
  Many checks are now performed and in this way common errors are avoided. 
- The units are now the same for all MD codes.
  If you want to use a different unit than the default you set it in the input file. 
- The analysis tools are now much more flexible.
  As an example of this it is now possible to write different collective variables with different frequencies.
- Many complex collective variables are considerably faster than they were in plumed1.
  In particular, all variables based on RMSD distances. 
- Centers of mass can be used as if they were atoms.
  Hence, unlike plumed 1.3, you can use center of mass positions in ALL collective variables.
- The virial contribution is now computed and passed to the MD code.
  Plumed can thus now be used to perform biased NPT simulations.
- Variables can be dumped on different files, and are
  computed only when this is necessary.
- PLUMED is now compiled as a separate library. This simplifies the patching
  procedure, but might require some extra work to configure PLUMED properly.
  Since PLUMED can be loaded as a shared library, it is possible to setup
  everything such that PLUMED and MD codes can be updated independently from each
  other.

In addition, it is now much easier to contribute new functionality to the code because: 
- There is a much simpler interface between plumed and the base MD codes.
  This makes it much easier to add plumed to a new MD code. Hopefully, in the future,
  interfaces with MD codes will be maintained by the developers of the MD codes
  independently from PLUMED developers. This will allow more MD codes
  to be compatible with PLUMED.
- There is C++ object oriented programming and full compatibility with the C++ standard library 
- A modular structure.
- New collective variables and methods can be released independently.
- There is an extensive developer documentation.
- User documentation is provided together inside the implementation files.

Caveats:
- PLUMED 2 input file (plumed.dat) has a syntax which is not
  compatible with PLUMED 1.
  Transition should be easy, but cannot
  be done just using the new version with the old input file.
- PLUMED 2 is written in C++, thus requires a C++ compiler
- PLUMED 2 may not include all the features that were available
  in PLUMED 1.

A tutorial explaining how to move from PLUMED 1 to PLUMED 2 is available (see \ref moving).

Version 2.0.1 (Nov 14, 2013)
----------------------------

For users:
- Fixed a bug in \ref HISTOGRAM with REWEIGHT_BIAS. Reweighting was only done when also temperature-reweighting was enabled.
- Fixed a bug that was sometime crashing code with domain decomposition and
  non-dense simulation boxes (e.g. implicit solvent).
- Performance improvements for \ref GYRATION.
- Flush all files every 10000 steps by default, without need to use \ref FLUSH
- Errors when writing input for \ref switchingfunction are now properly
  recognized.
- Added message when \ref simplemd is used on a non-existing file.
- Fixed `plumed mklib` such that it deletes the target shared library in case
  of compilation error.
- Several small fixes in documentation and log file.

For developers:
- Added possibility to setup replica exchange from MD codes in Fortran (commands "GREX setMPIFIntercomm" and "GREX setMPIFIntracomm").
- cmd("setStopFlag") should now be called after PLUMED initialization.
- Several small fixes in documentation.

Version 2.0.2 (Feb 11, 2014)
----------------------------

For users:
- Fixed bug with \ref METAD with INTERVAL and replica exchange, including bias exchange.
  Now the bias is correctly computed outside the boundaries. Notice that this is different
  from what was done in PLUMED 1.3. Also notice that INTERVAL now works
  correctly with grids and splines.
- Fixed bug with \ref READ and periodic variables.
- Fixed bug with \ref HISTOGRAM (option USE_ALL_DATA was not working properly).
- Gromacs patch updated to 4.6.5.
- Gromacs patch for 4.6 has been modified to allow for better load balancing when
  using GPUs.
- Added option 'plumed info --long-version' and 'plumed info --git-version'.
- Added full reference (page/number) to published paper in doc and log.
- Fixed a bug in file backups (only affecting Windows version - thanks to T. Giorgino).
- Added possibility to search in the documentation.
- Several small fixes in documentation and log file.

For developers:
- Fixed Makefile dependencies in some auxiliary files in src/lib (*cmake and *inc).
- Changed way modules are linked in src/.
  E.g. src/colvar/tools/ is not anymore a symlink to src/colvar but a real directory.
  (Notice that this introduces a regression: when using plumed as an external library
  some include files could not work - this only applies when plumed is installed;
  also notice that this is fixed in 2.0.3)
- Patch for gromacs 4.6 now also include original code so as to simplify its modification.
- Added option 'plumed patch --save-originals'.
- Fixed regtest regtest/secondarystructure/rt32 to avoid problems with NUMERICAL_DERIVATIVES.
- Removed include graphs in the documentation (too large).
- Several small fixes in documentation.

Version 2.0.3 (June 30, 2014)
----------------------------

For users:
- Now compiles on Blue Gene Q with IBM compilers.
- Fixed bug in \ref CENTER where default WEIGHTS were missing. 
- Fixed broken \ref CONTACTMAP with SUM
- Fixed \ref DUMPATOMS with gro file and more than 100k atoms.
- Added CMDIST in \ref CONTACTMAP to emulate plumed1 CMAP.
- Several small fixes in documentation and log file.

For developers:
- Fixed cmd("getBias") to retrieve bias. It was not working with
  single precision codes and it was not converting units properly.
- Fixed a regression in 2.0.2 concerning include files from installed plumed
  (see commit 562d5ea9dfc3).
- Small fix in tools/Random.cpp that allows Random objects to be
  declared as static.
- Small fix in user-doc compilation, so that if plumed is not found
  the sourceme.sh file is sourced
- Fixed non-ANSI syntax in a few points and a non-important memory leakage.
- Split cltools/Driver.cpp to make parallel compilation faster.

Version 2.0.4 (September 15, 2014)
----------------------------------------------

For users:
- Fixed a bug in \ref BIASVALUE that could produce wrong acceptance with replica exchange simulations.
- Fixed a few innocuous memory leaks.
- Fixed reader for xyz files, that now correctly detects missing columns. Also a related regtest has
  been changed.
- Several small fixes in documentation and log file.

For developers:
- Renamed Value.cpp to BiasValue.cpp

Version 2.0.5 (December 15, 2014)
----------------------------------------------

\plumednotmaintained

For users:
- Fixed a bug in replica exchange with different Hamiltonians (either lambda-dynamics
  or plumed XX-hrex branch) possibly occurring when using charge or mass dependent
  variables.
- Fixed a bug in analysis (e.g. \ref HISTOGRAM) leading to wrong accumulation
  of statistics when running a replica exchange simulation.
- Fixed a bug in the calculation of derivatives in histograms. This should
  be harmless since people usually only consider the value in histograms
  and not the derivatives.
- Fixed an issue in Makefile that could results in problems when
  patching an MD code with --shared option (pointed out by Abhi Acharya).
  This fixes a regression introduced in 2.0.2.
- Small fixes in documentation.

For developers:
- Added warning when performing regtests using an instance of plumed from
  a different directory

@page CHANGES-2-8 Version 2.8
  
## Version 2.8 (Feb 22, 2022)

This page contains changes that will end up in 2.8

- Changes leading to differences with previous versions
  - in \ref METAD and \ref PBMETAD, Gaussians are now stretched rather than truncated, making the energy a continuous function
    of the collective variable. See \issue{420}.
  - \ref sum_hills is now aware of stretched Gaussians. This change also fixes a minor bug in the set of grid points
    where Gaussians were different from zero that is still present up to version 2.7.
  - it is possible to restart from a HILLS file produced with PLUMED < 2.8, but Gaussians will be reinterpreted as stretched
    and a warning will be written in the log file. This might lead to small numerical changes in bias potentials.
  - in \ref METAD if possible the root walker in WALKERS_MPI will set the folder from which reading the GRID/HILLS file upon restart
  - in \ref METAD work is not calculated by default anymore, if needed it can be obtained using CALC_WORK
  - in \ref METAD an error will be thrown if, when restarting from FILE, the file is not found
  - the parser is more strict. Specifically, the explicitly crashes when a string cannot be parsed correctly.
    This was true only in a limited number of cases until v2.7 and might lead to errors when reading incorrectly
    formatted files. See \issue{717}.

- New actions:
  - \ref GHBFIX to compute generalized hydrogen-bond fixes

- New contributed module:
  - A new SASA module by Andrea Arsiccio
     - \ref SASA_HASEL
     - \ref SASA_LCPO
  - A new S2 contact model module by Omar Valsson 
     - \ref S2CM

- Fixed patches:
  - A bug in using GROMACS with expanded ensemble in combination with PLUMED has been fixed (version 2020.6 and 2021.4, see \issue{793}).
    Notice that this fix requires PLUMED 2.8, so it won't be backward compatible.

- Other improvements
  - in \ref METAD a new keyword NLIST has been added to use a neighbor list for bias evaluation, this should be faster than grids with many CVs
  - in \ref METAD there are more checks that a restart of WALKERS_MPI is working consistently among walkers
  - in \ref driver there is a flag `--restart` that can be used to enforce restart (similar to using \ref RESTART in the PLUMED input file).
  - Added configure option `--enable-cxx`. Can be used to select C++14 with `--enable-cxx=14`. Required to compile against libraries
    whose header files need C++14.

- Changes in the OPES module
  - new action \ref OPES_EXPANDED
  - various new actions of type \ref EXPANSION_CV to be used with \ref OPES_EXPANDED
  - new action \ref OPES_METAD_EXPLORE
  - new option EXTRA_BIAS in \ref OPES_METAD, to sample custom target distributions
  - new option EXCLUDED_REGION in \ref OPES_METAD, to define a region where no kernels are deposited

- Changes in the VES module
  - New localized basis functions: Wavelets (\ref BF_WAVELETS), Gaussians (\ref BF_GAUSSIANS), and cubic splines (\ref BF_CUBIC_B_SPLINES). In particular, symmetric wavelets (symlets) have shown the best performance and are recommended of the localized basis functions. Furthermore, symlets have been shown to perform better than delocalized Chebyshev and Legendre polynomials.  
  - New optimizer based on Adam (\ref OPT_ADAM). Still experimental, and restarting with it does not work yet. 
  - New optimizer based on classical Robbins Monro stochastic gradient descent (\ref OPT_ROBBINS_MONRO_SGD). Only included for reference and not recommended for usage in simulations. 
  - Fixed a bug in \ref VES_LINEAR_EXPANSION for multidimensional bias potential if one (or more) of the CVs is outside the range of the bias potential. Previously, there was a force acting on the CVs if this happened. Now, there is no biasing force acting on the CVs if one (or more) of the CVs is outside the bias potential range. 

- Changes in the DRR module
  - Added a new option MERGEHISTORYFILES to output a single history file instead of many .drrstate files.

- For developers:
  - The C++ interface now performs type checking (see https://github.com/plumed/plumed2/pull/653).
    This should require no change for MD codes that were calling PLUMED with correct arguments.
    Checks could be disabled at runtime with `export PLUMED_TYPESAFE_IGNORE=yes`.
  - Two new Fortran modules have been added. One of them provides explicit interfaces for the already available wrappers.
    With no change in calling code, just by including this module, one could perform runtime type/shape checking.
    In addition, a novel object oriented Fortran interface has been designed which allow to better manipulate
    PLUMED instances from Fortran.
    Both interfaces were written with a significant help from Balint Aradi.
  - The C interface (`plumed_cmd`) also performs type checking and allows overload-like syntax to pass
    additional size and shape information. This is obtained redefining `plumed_cmd` to a macro that calls the C++ interface,
    when using a C++ compiler, or using C11 `_Generic`, if the C compiler supports it.
    This feature is not supported if used a pre-C11 C compiler (pre-C++11 C++ compilers are ok instead).
  - `xxd` replaced by a awk script. This removed the build dependence on vim.
  - Lepton has been updated with OpenMM 7.6.0
  - Asmjit is now enabled by default on supported architectures.
  - Xdrfile library is now embedded and always available.
  - `--enable-rpath` now also includes the path where `libplumedKernel.so` is installed (see \issue{767}).
@page CHANGES-2-6 Version 2.6
  
## Version 2.6 (Jan 27, 2020)

Changes from version 2.5 which are relevant for users:
- Changes leading to incompatible behavior:
  - PLUMED input file parsing is now case insensitive that is that all directives can be written using uppercase characters (compatible with former versions) as well as lowercase characters (not compatible) internally PLUMED still uses uppercase definitions
  - `plumed partial_tempering` now uses `gawk` instead of `awk`. You might need to install `gawk` for it to work correctly.

- Other changes:
  - Asmjit is now embedded into PLUMED. In order to enable it, it is sufficient to configure with `--enable-asmjit`. See \ref Lepton "this page".
  - Fixed grids so as to decrease memory footprint of derivatives (see \issue{465}).
  - Added option `--idlp4` to \ref driver to read DLPOLY4 HISTORY files (see \issue{478}, thanks to Alin Marin Elena).
  - Added atom selectors using mdtraj/MDAnalysis/VMD syntax, see \ref MOLINFO and \issue{448}.
  - \ref EEFSOLV is now faster in scalar and also mpi/openmp parallel
  - New shortcuts are available for selecting protein atoms: `@sidechain-#`, `@back-#`
  - VIM syntax highlight is now case insensitive. Notice that autocompletion still only works with upper case commands.

- New contributed modules:
  - A new Maze module by Jakub Rydzewski
     - \ref MAZE_LOSS
     - \ref MAZE_MEMETIC_SAMPLING
     - \ref MAZE_RANDOM_ACCELERATION_MD
     - \ref MAZE_RANDOM_WALK
     - \ref MAZE_SIMULATED_ANNEALING
     - \ref MAZE_STEERED_MD
     - \ref MAZE_OPTIMIZER_BIAS
  - A new ANN module by Wei Chen and Andrew Ferguson
     - \ref ANN

- New patches:
  - added support for AMBER PMEMD 18 (contributed by Viktor Drobot, see \issue{486}).

- Changes in the VES module
  - new \ref VES_DELTA_F bias.
  - ves_md_linearexpansion now outputs one-dimensional free energy projections of the potential energy landscape. 

- Changes in the DRR module
  - The MAXFACTOR option now is tunable for each CV in multidimensional cases.
  - Output .zcount file (the same as .czar.count) for compatibility with newer abf_integrate.
  - The citation of DRR module has been updated.

- Changes in the ISDB module
  - in \ref METAINFERENCE we removed the MC_STRIDE keyword
  - in \ref METAINFERENCE the bias value (metainference score) now includes the Jeffrey's prior (values are different, but forces are equal)
  - components were previously named using _ but now they abide to the standard is -
  - removed ADDEXP keywords for \ref JCOUPLING \ref NOE \ref PRE \ref RDC
  - \ref METAINFERENCE performs more check on the input and restart files to ensure a consistent setup
  - \ref SAXS is slightly faster and scales better, removed BESSEL options

- Python module:
  - Removed compatibility with Python 2.
  - Added capability to read and write pandas dataset from PLUMED files (see \issue{496}).

Changes from version 2.5 which are relevant for developers:
  - Components documentation is now enforced
  - `readdir_r` is deprecated and is thus not used by default (can be enabled with `./configure --enable-readdir-r`).

## Version 2.6.1 (Jul 8, 2020)

For users:
- Includes all fixes up to 2.5.5
- New patches:
  - added gromacs 2019.6 
  - added gromacs 2020.2 (experimental) 
- Fixed handling of truncated octahedron box in Amber (see \issue{584}).
  Notice that the fix is for the PMEMD patch to be used with Amber 18.
  Amber 20 has been fixed upstream, both in PMEMD and Sander code.

For developers:
- Small fix to avoid unique global symbols (see \issue{549})

## Version 2.6.2 (Oct 26, 2020)

For users:
- Includes all fixes up to 2.5.6
- Updated patches:
  - added gromacs 2020.4 (experimental: it does not yet support modular simulator) 

## Version 2.6.3 (Apr 16, 2021)

For users:
- Includes all fixes up to 2.5.7

## Version 2.6.4 (Jul 27, 2021)

For users:
- Fixed `plumed partial_tempering` so as to correctly process `[ pairs ]` sections.
  The incorrect script was leading to unscaled 14 interactions with Glycam force field.
  (reported by Isabell Grothaus).

For developers:
- Added integer macros `PLUMED_VERSION_MAJOR` `PLUMED_VERSION_MINOR` and `PLUMED_VERSION_PATCH` to `config/version.h`.
  Can be used to write \ref LOAD -able source code portable across multiple versions.
- Fix for compilation with GCC 11 (reported by Axel Kohlmeyer, see \issue{693}).

## Version 2.6.5 (Dec 1, 2021)

For users:
- Fixed configure problem on XL compiler (see \issue{731}).
- Fixed a bug in \ref METAINFERENCE where the score was not properly updated upon multiple MC moves in the same MD step

For developers:
- Fixed several regtests decreasing their numeric precision.

## Version 2.6.6 (Feb 22, 2022)

For users:
- Fixed some incorrectly formatted output

For developers:
- Several fixes to improve portability on Debian and FreeBSD
<!--
  Feel free to delete not relevant sections below
-->

##### Description

<!-- describe here your pull request -->

##### Target release

<!-- please tell us where you would like your code to appear (e.g. v2.4): -->
I would like my code to appear in release __XXXXX__

##### Type of contribution

<!--
  Please select the type of your contribution among these:
  (Change [ ] to [X] to tick an option)
-->
- [ ] changes to code or doc authored by PLUMED developers, or additions of code in the core or within the default modules
- [ ] changes to a module not authored by you
- [ ] new module contribution or edit of a module authored by you

##### Copyright

<!--
  In case you picked one of the first two choices
  MAKE SURE TO TICK ALSO THE FOLLOWING BOX
-->

- [ ] I agree to transfer the copyright of the code I have written to the PLUMED developers or to the author of the code I am modifying.

<!--
  In case you picked the third choice (new module authored by you)
  MAKE SURE TO TICK ALSO THE FOLLOWING BOX
-->

- [ ] the module I added or modified contains a `COPYRIGHT` file with the correct license information. Code should be released under an open source license. I also used the command `cd src && ./header.sh mymodulename` in order to make sure the headers of the module are correct. 

##### Tests

<!--
  Make sure these boxes are checked. For Travis-CI tests, you can wait for them
  to be completed monitoring this page after your pull request has been submitted:
  http://travis-ci.org/plumed/plumed2/pull_requests
-->

- [ ] I added a new regtest or modified an existing regtest to validate my changes.
- [ ] I verified that all regtests are passed successfully on [GitHub Actions](https://github.com/plumed/plumed2/actions).

<!--
  After your branch has been merged to the desired branch and then to plumed2/master, and after the
  plumed official manual has been updated, please check out the coverage scan at
  http://www.plumed.org/coverage-master
  In case your new features are not well covered, please try to add more complete regtests.
-->

You are welcome to contribute to making PLUMED a better tool.

If you want to submit a small change in the code, such as fixing a bug
or adding a small feature, the best way to do it is to open a pull request.
If it is a bug fix, please open a merge request on the oldest maitained branch
where the bug appears. If you are submitting a new feature, please open a merge request on the
latest developer branch (master). Please ensure that your change does not break any of the
existing regression tests on Travis. In addition, add regression tests for any new features that you have added.  If you do not do 
this it is very likely that someone else will break the code for your new feature it in the future.
We (the developers) need to be certain that we will be able to maintain
your code in the future.  Consequently, please  be ready to answer to specific questions about your changes.
Finally, we are very happy to accept contributions to the documentation.

Notice that when you open a pull request you
*implictly agree to transfer the copyright of your code to the PLUMED developers* (or to the authors
of the code that you are modifying).
We understand that you might think this unfair.  However, we want to be 100% sure that in the
future we can make drastic changes to the code, including changes to the  license and that we will not have to 
contact all the developers that contributed a small number of lines when doing so.

If you want to contribute some large change, please consider adding a new module.
Documentation about adding new modules is still limited, but you can get inspiration
from the existing ones. This procedure will allow you to keep the ownership on your code.
On the other hand, we expect that you will maintain it in the future.
In order to incorporate a new module into the main repository, we ask contributors to declare that
the module is available with an open source license.

Finally, notice that you can always share modified versions of PLUMED with your changes.
We are happy if you want to host on github a fork of PLUMED with additional features.
\page tutorials Tutorials 

The following pages describe how to perform a variety of tasks using PLUMED

@TUTORIALS@

In addition, the following websites contain resources that might be helpful  

@WEBSITES@

Some older tutorials are also available here (some of them cover topics not covered by the recent tutorial but syntax may be outdated): 

- \subpage oldtutorials

\page oldtutorials Old Tutorials

@OLDTUTORIALS@


\page Miscellaneous Miscellaneous

- \subpage comments
- \subpage ContinuationLines
- \subpage VimSyntax
- \subpage BashAutocompletion
- \subpage includes
- \subpage load
- \subpage embed
- \subpage degub
- \subpage exchange-patterns
- \subpage mymodules
- \subpage special-replica-syntax
- \subpage parsing-constants
- \subpage misc

\page comments Comments

If you are an organized sort of person who likes to remember what the hell you were trying to do when you ran a 
particular simulation you might find it useful to put comments in your input file.  In PLUMED you can do this as 
comments can be added using a # sign.  On any given line everything after the # sign is ignored so 
erm... yes add lines of comments or trailing comments to your hearts content as shown below (using Shakespeare is optional):

\plumedfile
# This is the distance between two atoms:
d1: DISTANCE ATOMS=1,2 
UPPER_WALLS ARG=d1 AT=3.0 KAPPA=3.0 LABEL=Snout # In this same interlude it doth befall.
# That I, one Snout by name, present a wall.
\endplumedfile
(see \ref DISTANCE and \ref UPPER_WALLS)

An alternative to including comments in this way is to use the command \subpage ENDPLUMED.  Everything in the PLUMED input after this
keyword will be ignored.

\page ContinuationLines Continuation lines

If your input lines get very long then editing them using vi and other such text editors becomes a massive pain in the arse.  
We at PLUMED are aware of this fact and thus have provided a way of doing line continuations so as to make your life that much 
easier - aren't we kind?  Well no not really, we have to use this code too.  Anyway, you can do continuations by using the "..." syntax
as this makes this: 

\plumedfile
DISTANCES ATOMS1=1,300 ATOMS2=1,400 ATOMS3=1,500 LABEL=dist
\endplumedfile
(see \ref DISTANCES)

equivalent to this:

\plumedfile
DISTANCES ...
  LABEL=dist
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed

... DISTANCES
\endplumedfile

Notice that the closing `...` is followed by the word `DISTANCES`. This is optional, but might be
useful to find more easily which is the matching start of the statement. The following is equally correct
\plumedfile
dist: DISTANCES ...
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed

...
\endplumedfile

Notice that PLUMED makes a check that the word following the closing `...` is actually identical to
the first word in the line with the first `...`. If not, it will throw an error.
Also notice that you might put more than one word in the first line. E.g.
\plumedfile
DISTANCES LABEL=dist ...
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed
...
\endplumedfile
or, equivalently,
\plumedfile
dist: DISTANCES ...
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed
...  
\endplumedfile

\page BashAutocompletion Using bash autocompletion

When possible, PLUMED tries to install bash autocompletion so that
you do not have to do anything. Just use the `<TAB>` key to complete
plumed commands (e.g. `plumed dr<TAB>`) or even options (e.g. `plumed driver --i<TAB>`).
In case this does not work, you might have to add the following lines to your .bashrc file:
\verbatim
_plumed() { eval "$(plumed --no-mpi completion 2>/dev/null)";}
complete -F _plumed -o default plumed
\endverbatim

\par Effect

When typing on the the shell you should observe the following behavior.
\verbatim
> plumed <TAB>
\endverbatim
will autocomplete with the names of the available PLUMED commands (e.g. `driver`, `help`, etc).
\verbatim
> plumed -<TAB>
\endverbatim
will autocomplete with the available PLUMED options (e.g. `--no-mpi`, etc).

PLUMED also knows which are the options available for each command
(e.g. `plumed driver --natoms`). So, the following
\verbatim
> plumed driver -<TAB>
\endverbatim
(notice the `-`) will autocomplete to the options of `plumed driver`. On the contrary
\verbatim
> plumed driver --ixtc <TAB>
\endverbatim
(notice the there is no `-` before `<TAB>`) will autocomplete to the files in the current directory.

Also notice that every time you use the `<TAB>` key to autocomplete the command `plumed` will be invoked.
This should allow the correct commands and options to be reported depending on the exact `plumed` command
in the current execution path. For instance, if you have multiple PLUMED versions installed with
env modules, you should be able to see the commands available in the currently loaded version.
Clearly, this feature will only be available if `plumed` can run on this machine (that is: will not work
if you are cross compiling). This is not a problem since you are not expecting to run the `plumed` command
in this specific case.

\par Technicalities

At configure time if the variable `BASH_COMPLETION_DIR` is defined it will be used to decide where PLUMED
autocompletion should be installed. Otherwise, configure will look for the presence of the `bash-completion` package
and, in case it is installed on the same prefix as PLUMED, also PLUMED autocompletion will be installed.
Finally, if none of these two conditions are satisfied, autocompletion will not be enabled. You will
have to change your bashrc file once adding the following lines:
\verbatim
_plumed() { eval "$(plumed --no-mpi completion 2>/dev/null)";}
complete -F _plumed -o default plumed
\endverbatim
The command `plumed completion` just writes on its standard output the body of a bash function that is
then used by bash to construct the autocompletion.
The `--no-mpi` flag makes it more likely that the command can be executed correctly e.g. when you are on the login node of a cluster and
PLUMED was compiled with MPI but the login node does not support MPI. In other cases, it is harmless.
The `-o default` options will make sure that if `plumed --no-mpi completion` returns an error the default bash completion
will be used. This is what will happen if you load an older PLUMED version for which the `completion` command is not available yet.
In future PLUMED versions the `plumed completion` command might return more sophisticated functions. You should
be able to benefit of these features without ever changing your bash configuration file again.

\par Multiple versions and suffixes

In case you have multiple versions of PLUMED installed in separate env modules there is nothing more to do.
However, if you have have multiple versions of PLUMED installed with different suffixes you should
consistently add more lines to your profile file. For instance, if you installed two executables named
`plumed` and `plumed_mpi` your configuration file should look like:
\verbatim
_plumed() { eval "$(plumed --no-mpi completion 2>/dev/null)";}
complete -F _plumed -o default plumed
_plumed_mpi() { eval "$(plumed_mpi --no-mpi completion 2>/dev/null)";}
complete -F _plumed_mpi -o default plumed_mpi
\endverbatim

\page VimSyntax Using VIM syntax file

For the impatient use:
- Add the following to your .vimrc file:
\verbatim
" Enable syntax
:syntax on
" This allows including the proper PLUMED syntax file:
:let &runtimepath.=','.$PLUMED_VIMPATH
" The former command requires PLUMED_VIMPATH to be set. Alternatively, use this:
" let &runtimepath.=',/usr/local/lib/plumed/vim'
" properly adjusted to the path where PLUMED is installed.
" This makes autocompletion work in the expected way:
:set completeopt=longest,menuone
" This enables bindings of F2/F3/F4 to plumed specific commands:
:let plumed_shortcuts=1
\endverbatim
- When you open a PLUMED input file, you can enable syntax highlighting with:
\verbatim
:set ft=plumed
\endverbatim
This will also enable autocompletion. Use `<CTRL-X><CTRL-O>` to autocomplete a word.
- If you want to fold multi-line statements, type
\verbatim
:setlocal foldmethod=syntax
\endverbatim
- While editing a plumed input file, you can use command `:PHelp` (or shortcut `<F2>`)
  to show in a split window a short help about the action defined in the line where the cursor is.
  Typing `:PHelp` again (or pushing `<F2>`) you will
  close that window. With `<CTRL-W><CTRL-W>` you go back and forth between the two windows.
- When you open a file starting with `#! FIELDS`, VIM will automatically understand it
  is a PLUMED output file (VIM filetype = plumedf) and will color fields and data columns with
  alternating colors. Typing `:PPlus` and `:PMinus` (or pushing `<F3>` and `<F4>`)
  you can move a highlighted column.

See below for more detailed instructions.

\par Configuration

When PLUMED is compiled, directories `help` and `syntax` will appear in `builddir/vim`.
They contain a VIM plugin that can be used to highlight proper PLUMED instructions in a PLUMED
input file and to quickly retrieve help.
There is also a file `builddir/vim/scripts.vim` that helps VIM in recognizing PLUMED output files.

\warning
Notice that these file do not appear if you are cross compiling.
In this case, you must copy the plugin files from another machine.

To make VIM aware of these files, you should copy them to your `$HOME/.vim` directory.
Later you can
enable plumed syntax with the command
\verbatim
:set ft=plumed
\endverbatim

If you work in an environment where several PLUMED versions are installed (e.g. using env modules),
we recommend the following procedure:
- Install PLUMED
- Add to your `.vimrc` file the following line:
\verbatim
:let &runtimepath.=','.$PLUMED_VIMPATH
\endverbatim

The modulefile provided with PLUMED should set the PLUMED_VIMPATH environment variable
to the proper path.
Thus, when working with a given PLUMED module loaded, you should be able to
enable to proper syntax by just typing
\verbatim
:set ft=plumed
\endverbatim
in VIM.
Notice that the variable `PLUMED_VIMPATH` is also set in the `sourceme.sh` script in the build directory.
Thus, if you modify your `.vimrc` file as suggested, you will be able to use the correct syntax both
when using an installed PLUMED and when running from a just compiled copy.
Finally, in case you have both a pre-installed PLUMED **and** you have your development version
the following command would give you the optimal flexibility:
\verbatim
:let &runtimepath.=','.$PLUMED_VIMPATH.',/opt/local/lib/plumed/vim/'
\endverbatim
The environment variable `PLUMED_VIMPATH`, if set, will take the precedence.
Otherwise, vim will resort to the hard coded path.
In this case we assumed that there is a PLUMED installed in `/opt/local/` (e.g. using MacPorts),
but you can override it sourcing a `sourceme.sh` file in the compilation directory
or loading a PLUMED module with `module load plumed`.

If you are tired of typing `:set ft=plumed`, you can use a modeline.
Add to your `.vimrc` file the following commands
\verbatim
:set modeline
:set modelines=5
\endverbatim
Then, at the beginning of your PLUMED input file, put the following comment:
\plumedfile
# vim:ft=plumed
d: DISTANCE ATOMS=1,2
RESTRAINT ARG=d AT=0.0 KAPPA=1.0
\endplumedfile
Now, every time you open this file, you will see it highlighted.

\par Syntax highlighting

The syntax file contains a definition of all possible PLUMED actions and keywords.
It is designed to allow for a quick validation of the PLUMED input file before running it.
As such, all the meaningful words in the input should be highlighted:
- Valid action names (such as `METAD`) and labels (such as `m:` or `LABEL=m`) will be
  highlighted in the brightest way (`Type` in VIM). Those are the most important words.
- Keyword and flag names (such as `ATOMS=` or `COMPONENTS` when part of the action \ref DISTANCE) will be highlighted with a different color
  (`Statement` in VIM).
- Values provided by users (such as the number of the atoms following `ATOMS=`) will be highlighted with a different color
  (`String` in VIM).
- Comments (see \ref comments) will be highlighted as comments (`Comment` in VIM).
- String `__FILL__` (extensively used in tutorials to indicate parts to be completed) is highlighted (`Todo` in VIM).

If you see something that is not highlighted and appears in black, this is likely going to result in an error at runtime.
Think of this as a sort of preliminary spell-check.
For this checks to be effective, we recommend to use a syntax file generated with
exactly the same version of PLUMED that you are using.
In case you find that parts of an input file that is valid are not highlighted, then please
report it as a bug.
On the contrary, you cannot expect the VIM syntax file to recognize all possible errors
in a PLUMED input. Thus, a file for  which the highlighting looks correct might still contain errors.

\par Multi-line folding

Notice that syntax highlighting also allow VIM to properly fold multi-line actions.
Try to do the following:
- Open a PLUMED input file
- Enable PLUMED syntax
\verbatim
:set ft=plumed
\endverbatim
- Enable syntax-based folding
\verbatim
:setlocal foldmethod=syntax
\endverbatim

Now look at what happened to all the multi-line statements in PLUMED (i.e. those using
\ref ContinuationLines).  As you can see, they will be folded into single lines.
Folded lines can be expanded with `zo` and folded with `zc`. Look at VIM documentation to
learn more.
In case you want to use this feature, we suggest you to put both label
and action type on the first line of multi-line statements. E.g.
\plumedfile
d: DISTANCE ATOMS=1,2

m: METAD ...
  ARG=d
  HEIGHT=1.0
  SIGMA=0.5
  PACE=100
...
\endplumedfile
will be folded to
\verbatim
d: DISTANCE ATOMS=1,2

+--  6 lines: m: METAD ...------------------------------------------------------
\endverbatim
and
\plumedfile
d: DISTANCE ATOMS=1,2

METAD LABEL=m ...
  ARG=d
  HEIGHT=1.0
  SIGMA=0.5
  PACE=100
...
\endplumedfile
will be folded to
\verbatim
d: DISTANCE ATOMS=1,2

+--  6 lines: METAD LABEL=m ...-------------------------------------------------
\endverbatim
This will allow you to easily identify the folded lines by seeing the most important information,
that is the action type (`METAD`) and its label (`m`). This feature is convenient if
you want to browse files that contain a lot of actions defined on multiple lines.

\par Autocompletion

Another VIM feature that comes when you load PLUMED syntax is autocompletion of PLUMED
actions and keywords. Open your favorite PLUMED input file and set it to PLUMED syntax highlighting with
\verbatim
:set ft=plumed
\endverbatim
Now go into insert mode pressing `i` and type `DU` followed by `<CTRL+X><CTRL+O>`.
Here `<CTRL+X>` stands for autocompletion and `<CTRL+O>` for omnifunc autocompletion. You will see a short
menu listing the following actions
\verbatim
DUMPATOMS       
DUMPDERIVATIVES 
DUMPFORCES      
DUMPMASSCHARGE  
DUMPMULTICOLVAR 
DUMPPROJECTIONS 
\endverbatim
That is, all the actions starting with `DU`.
You can navigate it with up and down arrows so as to choose the
best match.

Notice that the default behavior of VIM is to use the first match by default.
In the first example (`DU<CTRL+X><CTRL+O`), it would be `DUMPATOMS`.
The following settings make it work as most of the people expect:
\verbatim
:set completeopt=longest,menuone
\endverbatim
With these settings, in the first example (`DU<CTRL+X><CTRL+O`) VIM will only complete up to the longest common part (`DUMP`).

As you can imagine,
if you use autocompletion after you have typed the word `DISTANCE` followed by a space you will see
a menu listing `LABEL=`, `COMPONENTS`, etc. Basically, all the keywords that are possibly used within a `DISTANCE` line
will be shown. This is very useful if you do not remember the exact name of the keywords associated with
a given action.

\par Quick help

You can also retrieve quick explanation of the input options for a specific action.
Try to do the following. Enable plumed syntax:
\verbatim
:set ft=plumed
\endverbatim
Then add the following line
\verbatim
DISTANCE
\endverbatim
Now, in normal mode, go with the cursor on the `DISTANCE` line and type
\verbatim
:PHelp
\endverbatim
A new split window should appear containing some documentation about the \ref DISTANCE collective variable.
You can go back and forth between the two windows with `<CTRL+W><CTRL+W>`, as usually in vim.
Notice that if you are in the help window and type `:PHelp` this window will be closed.

To make the navigation easier, you can add a shortcut in your .vimrc file. For example, adding:
\verbatim
: nmap <F2> : PHelp<CR>
\endverbatim
you should be able to open and close the manual hitting the F2 key.
This is done automatically in the PLUMED syntax file if you add `let plumed_shortcuts=1` to your
vimrc file.

\par Displaying output files

Most of the PLUMED output files look like this
\verbatim
#! FIELDS A B C
1 2 3
\endverbatim
This is useful since in the header you can see the name of the quantities that are printed in the
data lines. However, when you have an output file with many columns it might be a bit error prone to count them.
To simplify this, when PLUMED syntax for VIM is configured properly VIM should be able to:
- Detect that this file is a PLUMED output file with fields, automatically setting its type to `plumedf`. If not, just type
  `:set ft=plumedf`.
- Show this file with syntax highlighting to increase its readability.

Notice that the syntax file for the output files (`plumedf.vim`) is not the same one that is used for the PLUMED
input file (`plumed.vim`).

To make output files more readable, vim will show `FIELDS` and `SET` words in a different color,
and data columns with alternating colors (e.g. dark/light/dark/light).
The colors in the columns are consistent with those shown in the FIELD line.
In the example above, 1, 2, and 3 will be of the same color as A, B, and C respectively.
This should make it much easier to find which columns correspond to a given quantity.

It is also possible to highlight a specific field of the file. Typing
\verbatim
:5PCol
\endverbatim
you will highlight the fifth field. Notice that in the `FIELDS` line (the first line of the file)
the seventh word of the line will be highlighted, which is the one containing the name of the field.
This allows for easy matching of values shown
in the file and tags provided in the `FIELDS` line.
The highlighted column can be moved back and forth using `:PPlus` and `:PMinus`.
Adding a count to the command will move the highlighted column more. E.g. `:2PPlus` will move
the column to the right twice.

If you have a long output file, it might be convenient to split it with
`:split` so that one of the two windows will only show the header. The other window
can be used to navigate the file.


To make the navigation easier, you can add a shortcut in your .vimrc file. For example, adding:
\verbatim
: map <F3> :PMinus<CR>
: map <F4> :PPlus<CR>
\endverbatim
you should be able to move the highlight column using F3 and F4 buttons.
This is done automatically in the PLUMED syntax file if you add `let plumed_shortcuts=1` to your
vimrc file.

\page includes Including other files 

If, for some reason, you want to spread your PLUMED input over a number of files you can use \subpage INCLUDE as shown below:

\plumedfile
INCLUDE FILE=filename
\endplumedfile

So, for example, a single "plumed.dat" file:

\plumedfile
DISTANCE ATOMS=1,2 LABEL=dist
RESTRAINT ARG=dist AT=2.0 KAPPA=1.0
\endplumedfile

could be split up into two files as shown below:
 
\plumedfile
DISTANCE ATOMS=1,2 LABEL=dist
INCLUDE FILE=toBeIncluded.inc
\endplumedfile
plus a "toBeIncluded.inc" file
\plumedfile
#SETTINGS FILENAME=toBeIncluded.inc
# this is toBeIncluded.inc
RESTRAINT ARG=dist AT=2.0 KAPPA=1.0
\endplumedfile

However, when you do this it is important to recognize that \ref INCLUDE is a real directive that is only resolved
after all the \ref comments have been stripped and the \ref ContinuationLines have been unrolled.  This means it
is not possible to do things like:

\plumedfile
# this is wrong:
DISTANCE INCLUDE FILE=options.dat
RESTRAINT ARG=dist AT=2.0 KAPPA=1.0
\endplumedfile

\page load Loading shared libraries

You can introduce new functionality into PLUMED by placing it directly into the src directory and recompiling the 
PLUMED libraries.  Alternatively, if you want to keep your code independent from the rest of PLUMED (perhaps
so you can release it independently - we won't be offended), then you can create your own dynamic library.  To use this 
in conjunction with PLUMED you can then load it at runtime by using the \subpage LOAD keyword as shown below:

\plumedfile
LOAD FILE=library.so
\endplumedfile
 
N.B.  If your system uses a different suffix for dynamic libraries (e.g. macs use .dylib) then PLUMED will try to 
automatically adjust the suffix accordingly.

\page embed Embed a separate PLUMED instance

\subpage PLUMED

\page degub Debugging the code

The \subpage DEBUG action provides some functionality for debugging the code that may be useful if you are doing 
very intensive development of the code of if you are running on a computer with a strange architecture.

\page exchange-patterns Changing exchange patterns in replica exchange

Using the \subpage RANDOM_EXCHANGES keyword it is possible to make exchanges between randomly
chosen replicas. This is useful e.g. for bias exchange metadynamics \cite piana.

\page special-replica-syntax Special replica syntax

(this part of the manual is based on \ref trieste-5-replica-special-syntax).

In many cases, we need to run multiple replicas with almost identical PLUMED files.
These files might be prepared with cut-and-paste, which is very error prone,
or could be set up with some smart bash or python script. Additionally,
one can take advantage of the \ref INCLUDE keyword so as to have a shared input
file with common definitions and specific input files with replica-dependent keywords.
However, as of PLUMED 2.4, we introduced a simpler manner to manipulate multiple replica
inputs with tiny differences. Look at the following example:

\plumedfile 
#SETTINGS NREPLICAS=3
# Compute a distance
d: DISTANCE ATOMS=1,2

# Apply a restraint.
RESTRAINT ARG=d AT=@replicas:1.0,1.1,1.2 KAPPA=1.0
# On replica 0, this means:
#   RESTRAINT ARG=d AT=1.0 KAPPA=1.0
# On replica 1, this means:
#   RESTRAINT ARG=d AT=1.1 KAPPA=1.0
# On replica 2, this means:
#   RESTRAINT ARG=d AT=1.2 KAPPA=1.0
\endplumedfile

If you prepare a single `plumed.dat` file like this one and feeds it to PLUMED while using 3 replicas,
the 3 replicas will see the very same input except for the `AT` keyword, that sets the position of the restraint.
Replica 0 will see a restraint centered at 1.0, replica 1 centered at 1.1, and replica 2 centered at 1.2.

The `@replicas:` keyword is not special for \ref RESTRAINT or for the `AT` keyword. Any keyword in PLUMED can accept that syntax.
For instance, the following single input file can be used to setup a bias exchange metadynamics \cite piana simulations:
\plumedfile
#SETTINGS NREPLICAS=2
# Compute distance between atoms 1 and 2
d: DISTANCE ATOMS=1,2

# Compute a torsional angle
t: TORSION ATOMS=30,31,32,33

# Metadynamics.
METAD ...
  ARG=@replicas:d,t
  HEIGHT=1.0
  PACE=100
  SIGMA=@replicas:0.1,0.3
  GRID_MIN=@replicas:0.0,-pi
  GRID_MAX=@replicas:2.0,pi
...
# On replica 0, this means:
#  METAD ARG=d HEIGHT=1.0 PACE=100 SIGMA=0.1 GRID_MIN=0.0 GRID_MAX=2.0
# On replica 1, this means:
#  METAD ARG=t HEIGHT=1.0 PACE=100 SIGMA=0.3 GRID_MIN=-pi GRID_MAX=pi
\endplumedfile

This would be a typical setup for a bias exchange simulation.
Notice that even though variables `d` and `t` are both read in both replicas,
`d` is only computed on replica 0 (and `t` is only computed on replica 1).
This is because variables that are defined but not used are never actually calculated by PLUMED.

If the value that should be provided for each replica is a vector, you should use curly braces as delimiters.
For instance, if the restraint acts on two variables, you can use the following input:

\plumedfile
#SETTINGS NREPLICAS=3
# Compute distance between atoms 1 and 2
d: DISTANCE ATOMS=10,20

# Compute a torsional angle
t: TORSION ATOMS=30,31,32,33

# Apply a restraint:
RESTRAINT ...
  ARG=d,t
  AT=@replicas:{{1.0,2.0} {3.0,4.0} {5.0,6.0}}
  KAPPA=1.0,3.0
...
# On replica 0 this means:
#  RESTRAINT ARG=d AT=1.0,2.0 KAPPA=1.0,3.0
# On replica 1 this means:
#  RESTRAINT ARG=d AT=3.0,4.0 KAPPA=1.0,3.0
# On replica 2 this means:
#  RESTRAINT ARG=d AT=5.0,6.0 KAPPA=1.0,3.0
\endplumedfile

Notice the double curly braces. The outer ones are used by PLUMED to know there the argument of the `AT` keyword ends,
whereas the inner ones are used to group the values corresponding to each replica.
Also notice that the last example can be split in multiple lines exploiting the fact that
within multi-line statements (enclosed by pairs of `...`) newlines are replaced with simple spaces:
\plumedfile
#SETTINGS NREPLICAS=3
d: DISTANCE ATOMS=10,20
t: TORSION ATOMS=30,31,32,33
RESTRAINT ...
  ARG=d,t
# indentation is not required (this is not python!)
# but makes the input easier to read
  AT=@replicas:{
    {1.0,2.0}
    {3.0,4.0}
    {5.0,6.0}
  }
  KAPPA=1.0,3.0
...
\endplumedfile

In short, whenever there are keywords that should vary across replicas, you should set them using the `@replicas:` keyword.
As mentioned above, you can always use the old syntax with separate input file, and this is recommended when the
number of keywords that are different is large.

\page parsing-constants Parsing constants

You might have noticed that from time to time constants are specified using strings rather than numbers.
An example is the following

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt65/AA.pdb
MOLINFO STRUCTURE=AA.pdb  MOLTYPE=rna
e1: TORSION ATOMS=@epsilon-1
t: METAD ARG=e1 SIGMA=0.15 PACE=10 HEIGHT=2 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=200
\endplumedfile

Notice that the boundaries for `GRID_MIN` and `GRID_MAX` are `-pi` and `pi`. Until PLUMED 2.3,
we used a very dummy parses that could recognize only `pi` as a special string, plus strings such
as `0.5pi` and `-pi`. However, as of version 2.4, we use the Lepton library in order to parse every constant
that we read. This means that you can also employ more complicated expressions such as `1+2` or `exp(10)`:

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt65/AA.pdb
MOLINFO STRUCTURE=AA.pdb  MOLTYPE=rna
e1: TORSION ATOMS=@epsilon-1
RESTRAINT ARG=e1 AT=1+0.5
\endplumedfile

Notice that this applies to any quantity read by plumed as a real number, but does not apply
yet to integer numbers (e.g.: the PACE argument of \ref METAD).

\page misc Frequently used tools

@DICTIONARY@
<TABLE ALIGN="center" FRAME="void" WIDTH="95%%" CELLPADDING="5%%">
<TR>
<TD WIDTH="5%"> 
\subpage Regex </TD><TD> </TD><TD> POSIX regular expressions can be used to select multiple actions when using ARG (i.e. \ref PRINT).
</TD>
</TR>
<TR>
<TD WIDTH="5%"> 
\subpage Files </TD><TD> </TD><TD> Dealing with Input/Output
</TD>
</TR>
</TABLE>

\page Regex Regular Expressions

When you use need to pass many arguments to a PLUMED action, being them
components of a few collective variables or also multiple collective variables,
you might find it convenient to use [regular expressions](https://en.wikipedia.org/wiki/Regular_expression).

Since version 2.1, plumed takes advantage of a configuration scripts that
detects libraries installed on your system. If regex library is found,
then you will be able to use regular expressions to refer to collective variables
or function names.

Regular expressions are enclosed in round braces and must not contain spaces (the components 
names have no spaces indeed, so why use them?).

As an example the command:
\plumedfile
d1: DISTANCE ATOMS=1,2 COMPONENTS
PRINT ARG=(d1\.[xy])   STRIDE=100 FILE=colvar FMT=%8.4f
\endplumedfile
will cause both the d1.x and d1.y components of the DISTANCE action to be printed.

Notice that selection does not happen in alphabetic order, nor in the order in which `[xy]` are listed, but rather in the order in which
the two variables have been created by PLUMED.
Also notice that the
`.` character must be escaped as `\.` in order to interpret it as a literal `.`. An un-escaped dot is a wildcard which is matched by any character,
So as an example
\plumedfile
d1: DISTANCE ATOMS=1,2 COMPONENTS
dxy: DISTANCE ATOMS=1,3

# this will match d1.x,d1.y,dxy
PRINT ARG=(d1.[xy])   STRIDE=100 FILE=colvar FMT=%8.4f

# while this will match d1.x,d1.y only
PRINT ARG=(d1\.[xy])   STRIDE=100 FILE=colvar FMT=%8.4f
\endplumedfile

You can concatenate more than one regular expression by using comma separated regular expressions.
The resulting matches will be concatenated:
\plumedfile
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
d1: DISTANCE ATOMS=7,17 COMPONENTS

# The first expression matches d1.x and d1.y
# The second expression matches t1 and t2
PRINT ARG=(d1\.[xy]),(t[0-9]) STRIDE=100 FILE=colvar FMT=%8.4f
# Thus this is the same as ARG=d1.x,d1.y,t1,t2
\endplumedfile

Be aware that if you have overlapping selections they will be duplicated.
As an alternative you could use the "or" operator `|`:
\plumedfile
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
d1: DISTANCE ATOMS=7,17 COMPONENTS

# Here is a single regular expression
PRINT ARG=(d1\.[xy]|t[0-9]) STRIDE=100 FILE=colvar FMT=%8.4f
# Thus this is the same as ARG=t1,t2,d1.x,d1.y
\endplumedfile

this selects the same set of arguments as the previous example.

\note
Be careful you do not confuse regular expressions, which are triggered by the parenthesis `()` and only available when
PLUMED has been compiled with the regex library, with the capability of PLUMED to use `*` as a wildcard in arguments:
\plumedfile
d1: DISTANCE ATOMS=1,2 COMPONENTS
# this is a regular expression that selects all components of d1
# i.e. d1.x d1.y and d1.z
PRINT ARG=(d1\..*)   STRIDE=100 FILE=colvar_reg FMT=%8.4f

# this is a wildcard that selects all the components of d1 as well
PRINT ARG=d1.*   STRIDE=100 FILE=colvar_wild FMT=%8.4f
\endplumedfile
Regular expressions are way more flexible than wildcards!

You can check the log to see whether or not your regular expression is picking the set of components you desire.

For more information on regular expressions visit http://www.regular-expressions.info/reference.html.

\page EDSMOD Experiment Directed Simulation

<!-- 
description: Methods for incorporating additional information about CVs into MD simulations by adaptively determined linear bias parameters
authors: Glen Hocky, Andrew White
reference: \cite white2014efficient \cite hocky2017cgds \cite Amirkulova2019Recent
-->

## Overview

This Experiment Directed Simulation module contains methods for adaptively determining linear bias parameters such that each biased CV samples a new target mean value. This module implements the stochastic gradient descent algorithm in the original EDS paper \cite white2014efficient as well as additional minimization algorithms for Coarse-Grained Directed Simulation \cite hocky2017cgds.
There is a recent review on the method and its applications here: \cite Amirkulova2019Recent.

Notice that a similar method is available as \ref MAXENT, although with different features and using a different optimization algorithm.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=eds' to your './configure' command when building PLUMED to enable these features.

## Usage
Currently, all features of the EDS module are included in a single EDS bias function: \ref EDS

A tutorial using EDS specifically for biasing coordination number can be found on <a href="http://thewhitelab.org/Blog/tutorial/2017/05/10/lammps-coordination-number-tutorial/">Andrew White's webpage</a>.

## Module Contents
- \subpage EDSMODBias

\page EDSMODBias Biases Documentation

The following list contains descriptions of biases developed for the PLUMED-EDS module. They can be used in combination with other biases outside of the EDS module.

@EDSMOD_BIAS@
\page ISDB PLUMED-ISDB

<!-- 
description: Integrative Structural and Dynamical Biology with PLUMED
authors: Max Bonomi and Carlo Camilloni
reference: \cite Bonomi:2017cc 
-->

Here are listed the collective variables, functions and biases originally developed for the Integrative Structural and Dynamical Biology module of PLUMED. They are related but not limited to the interpretation and modelling of experimental data in molecular modelling.

- \subpage ISDBColvar
- \subpage ISDBFunction
- \subpage ISDBGeneric
- \subpage ISDBBias

Additional tutorials focused on the ISDB module are included in the following and are meant as advanced tutorials.

- \subpage ISDBTutorial

\page ISDBColvar CVs Documentation

The following list contains descriptions of a number of the colvars that are currently implemented in the PLUMED-ISDB module.
These collective variables are related to the definitions of models to interpret experimental observables. They can be used in combination with any other collective variable, function or bias also outside the ISDB module.

@ISDB_COLVAR@

\page ISDBFunction Functions Documentation

The following list contains descriptions of functions originally developed for the PLUMED-ISDB module. They can be used in combination with any other collective variable, function or bias also outside the ISDB module.

@ISDB_FUNCTION@

\page ISDBGeneric General Actions Documentation

The following list contains descriptions of actions originally developed for the PLUMED-ISDB module. They can be used in combination with any other collective variable, function or bias also outside the ISDB module. 

Using \ref SELECTOR it is possible to define a variable inside the PLUMED code that can be used and modified by other actions. For example, a \ref SELECTOR can be used in combination with \ref RESCALE to activate a simulated-tempering like approach.

@ISDB_GENERIC@

\page ISDBBias Biases Documentation

The following list contains descriptions of biases originally developed for the PLUMED-ISDB module. They can be used in combination with any other collective variable, function or bias also outside the ISDB module.

@ISDB_BIAS@

\page ISDBTutorial Tutorials

The following are tutorials meant to learn how to use the different methods implemented in the ISDB module.

@ISDB_TUTORIALS@


\mainpage Introduction

PLUMED is a plugin that works with a large number of molecular dynamics codes (\ref codes ). 
It can be used to analyze features of the dynamics on-the-fly or to perform a wide variety of free energy methods.
PLUMED can also work as a \ref tools to perform analysis on trajectories saved in most of the
existing formats. If PLUMED is useful for your work please read and cite \cite plumed2, if you are interested in 
the PLUMED 1 original publication please read and cite \cite plumed1 .

To follow the development of PLUMED 2, you can look at the detailed \ref ChangeLog .

To install PLUMED, see this page: \ref Installation , while in \ref Syntax you can find a brief introduction on how to write your first PLUMED input file.

\ref tutorials are available to introduce basic as well as more advanced features of PLUMED.
 
\section AboutManual About this manual

@VERSION@

This is the user manual -  if you want to modify PLUMED or to understand how it works internally, have a look at the 
<a href="../../developer-doc/html/index.html"> developer manual </a>.

@PDFMANUAL@

\section codes Codes interfaced with PLUMED 

PLUMED can be incorporated into an MD code and used to analyze or bias a molecular dynamics run on the fly.
Some MD code could already include calls to the PLUMED library
and be PLUMED-ready in its original distribution.
As far as we know, the following MD codes can be used with PLUMED out of the box:
- [Amber](http://ambermd.org/), pmemd module, since version 20.
- [AmberTools](http://ambermd.org/), sander module, since version 15.
- [CP2K](http://www.cp2k.org), since Feb 2015.
- [ESPResSo](http://espressomd.org), in a version that has been patched with PLUMED can be found
  [here](http://davidebr.github.io/espresso/).
- [PINY-MD](http://github.com/TuckermanGroup/PINY), in its plumed branch.
- [IPHIGENIE](http://sourceforge.net/projects/iphigenie/).
- [AceMD](http://www.multiscalelab.org/acemd/), see [this link](https://github.com/tonigi/ACEMD-PLUMED).
- [OpenMM](http://openmm.org), using the [openmm-plumed plugin](http://github.com/peastman/openmm-plumed).
- [DL_POLY4](https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx).
- [VNL-ATK](https://quantumwise.com), see [this link](https://docs.quantumwise.com/tutorials/metadynamics_with_plumed/metadynamics_with_plumed.html).
- [ABIN](https://github.com/PHOTOX/ABIN).
- [i-pi](https://github.com/i-pi/i-pi).
- [LAMMPS](https://lammps.sandia.gov/) since Nov 2018.
- [Yaff](https://github.com/molmod/yaff), since Jul 2019.
- [DFTB+](https://www.dftbplus.org/), since release 20.1.
- [Metalwalls](https://gitlab.com/ampere2/metalwalls)

Please refer to the documentation of the MD code to know how to use it with the latest PLUMED release.
If you maintain another MD code that is PLUMED-ready let us know and we will add it to this list.

Additionally, we provide patching procedures for the following codes:

@CODESL@

Alternatively, one
can use PLUMED as a \ref tools for post processing the results from molecular dynamics 
or enhanced sampling calculations.  Notice that PLUMED can be used as an analysis tool
also from the following packages:
- [PLUMED-GUI](http://github.com/tonigi/vmd_plumed) is a [VMD](http://www.ks.uiuc.edu/Research/vmd/) plugin that computes PLUMED collective variables.
- [HTMD](http://www.htmd.org/) can use PLUMED collective variables for analysis.
- [OpenPathSampling](http://openpathsampling.org/), using the [PLUMED Wrapper for OpenPathSampling](https://e-cam.readthedocs.io/en/latest/Classical-MD-Modules/modules/OpenPathSampling/ops_plumed_wrapper/readme.html).

\page tools Command Line Tools

PLUMED contains a number of simple command line tools.  To use one of these tools 
you issue a command something like:

\verbatim
plumed <toolname> <list of input flags for that tool>
\endverbatim

The following is a list of the various standalone tools that PLUMED contains.

@TOOLS@

For all these tools and to use PLUMED as a plugin in an MD calculation you will need an input file.

\page ChangeLog Change Log

Here you can find a history of changes across different PLUMED versions.
The future releases are expected to follow more or less the pace
of the old release. This means:
- Approximately once per year, after summer, a new release (2.X). These releases
  typically group together all the features that were contributed during the
  year.
- Approximately every three month, we announce a patch (e.g. 2.2.X).
  This typically contains bug fixes, and could occasionally contain a new feature.

A few months before each new release we provide a beta release.
We typically maintain release branches until the fifth patch release (2.X.5),
which should come out approximately 15 month after the original release (2.X).
After that, branches are not supported anymore.

Notice that occasionally we publish patches on the mailing list.
These patches are always included in the following release, but we encourage
users that want to be up to date to follow the mailing list.

Below you can find change logs for all the published releases.
We mostly add new features without breaking existing ones.
However, some of the changes lead to incompatible behavior.
In the Change Log we try to give as much visibility as possible to these changes
to avoid surprises.

We also log changes that are relevant if you are developing the code. These
change lists are however not complete, and if you want to put your hands in the code
and maintain your own collective variables we suggest you to follow the development
on github.

@CHANGES@

\page LOGMFDMOD Logarithmic Mean Force Dynamics

<!-- 
description: Method for enhanced sampling and for free energy calculations along collective variables
authors: Tetsuya Morishita, Naoki Watanabe
reference: \cite MorishitaLogMFD \cite MorishitaLogPD \cite MorishitaVsLogMFD
-->

## Overview

The LOGMFD module contains the LogMFD/LogPD method for enhanced sampling in a CV space and for on-the-fly free energy reconstruction along the CVs. This module implements the multiple-replica algorithm (LogPD \cite MorishitaLogPD) as well as the single-replica algorithm (LogMFD \cite MorishitaLogMFD), the former invoking the Crooks-Jarzynski non-equilibrium work relation. In addition, TAMD/d-AFED \cite AbramsJ2008 can also be implemented by this module.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=logmfd' to your './configure' command when building PLUMED to enable these features.

## Usage
Currently, all features of the LOGMFD module are included in a single LOGMFD bias function: \ref LOGMFD

## Module Contents
- \subpage LOGMFDMODBias

\page LOGMFDMODBias Biases Documentation

The following list contains descriptions of biases developed for the PLUMED-LOGMFD module. They can be used in combination with other biases outside of the LOGMFD module.

@LOGMFDMOD_BIAS@
\page S2CMMOD S2 contact model collective variable

<!-- 
description: S2 contact model collective variable (S2CM) 
authors: Omar Valsson
reference: \cite Palazzesi_s2_2017  
-->

## Overview

S2 contact model CV used in \cite Palazzesi_s2_2017, based on NH order parameter from \cite Zhang_s2_2002 and methyl order parameter from \cite Ming_s2_2004.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=s2cm' to your './configure' command when building PLUMED to enable these features.

## Usage
Currently, all features of the S2 contact model module are included in a single S2 contact model collective variable: \ref S2CM

## Module Contents
- \subpage S2CMMODColvar 

\page S2CMMODColvar CVs Documentation

@S2CMMOD_COLVAR@
\page Performances Performances 

In this page we collect hints on how to use the features available in PLUMED to speed
up your calculations. Please note that PLUMED performs many different tasks, it
can calculate a number of different collective variables, functions of collective 
variables, bias, on-the-fly analysis, etc in a way that is compatible with a number of
different molecular dynamics codes. This means that there cannot be a single 
strategy to speed up all the possible calculations. 

\ref performance-optimization "Here" 
you can find a step-by-step tutorial on optimizing PLUMED performances, discussing some of the topics
below in more detail and using practical examples.

PLUMED makes use of MPI and OpenMP to parallelize some of its functions, try to always
compile it with these features enabled. Furthermore, newer compilers with proper optimization 
flags can provide a dramatic boost to performances.

PLUMED collects atoms from an external code and sends back forces, so it is key to minimize
the effect of PLUMED on highly parallel calculations to keep to the minimum the number of atoms 
used by PLUMED at every calculation step. The less is the number of atoms you need to send 
to PLUMED the less will be the overhead in the communication between PLUMED and the code.

In the following you can find specific strategies for specific calculations, these could
help in taking the most by using PLUMED for your simulations.

- \subpage GMXGPU 
- \subpage Metadyn
- \subpage MTS
- \subpage Multicolvar 
- \subpage Neighbour 
- \subpage Openmp
- \subpage Secondary
- \subpage Time
- \subpage Lepton

\page GMXGPU GROMACS and PLUMED with GPU

Since version 4.6.x GROMACS can run in an hybrid mode making use of both
your CPU and your GPU (either using CUDA or OpenCL for newer versions of
GROMACS). The calculation of the short-range non-bonded interactions is 
performed on the GPU while long-range and bonded interactions are at the
same time calculated on the CPU. By varying the cut-off for short-range
interactions GROMACS can optimize the balance between GPU/CPU loading 
and obtain amazing performances.

GROMACS patched with PLUMED takes into account PLUMED in its load-balancing, 
adding the PLUMED timings to the one resulting from bonded interactions and long-
range interactions. This means that the CPU/GPU balance will be optimized 
automatically to take into account PLUMED!  

It is important to notice that the optimal setup to use GROMACS alone
on the GPU or GROMACS + PLUMED can be different, try to change the number
of MPI/OpenMP processes (\ref Openmp) used by GROMACS and PLUMED to find
optimal performances. Remember that in GROMACS multiple MPI threads
can use the same GPU:

i.e. if you have 4 cores and 2 GPU you can:

- use 2 MPI/2GPU/2OPENMP:

\verbatim
export PLUMED_NUM_THREADS=2
mpiexec -np 2 gmx_mpi mdrun -nb gpu -ntomp 2 -pin on -gpu_id 01
\endverbatim

- use 4 MPI/2GPU:

\verbatim
export PLUMED_NUM_THREADS=1
mpiexec -np 4 gmx_mpi mdrun -nb gpu -ntomp 1 -pin on -gpu_id 0011
\endverbatim

Of notice that since plumed 2.5 and gromacs 2018.3 the number of openMP threads can automatically set by gromacs (so PLUMED_NUM_THREADS is not needed, and the number of OpenMP threads used by plumed is set by -ntomp)

\verbatim
mpiexec -np 2 gmx_mpi mdrun -nb gpu -ntomp 2 -pin on -gpu_id 01
\endverbatim


\page Metadyn Metadynamics

Metadynamics can be sped up significantly using grids,
which are activated setting the GRID_MIN and GRID_MAX keywords of \ref METAD.
This makes addition of a hill to the list a bit slower (since
the Gaussian has to be evaluated for many grid points)
but the evaluation of the potential very fast. Since
the former is usually done every few hundred steps, whereas the latter 
typically at every step, using grids will make the simulation
 faster in particular for long runs.

Notice that when restarting a simulation the history is read  by default
from a file and hills are added again to the grid.
This allows one to change the grid boundaries upon restart. However,
the first step after restart is usually very slow.
Since PLUMED 2.3 you can also store the grid on a file
and read it upon restart. This can be particularly
useful if you perform many restarts and if your hills are large.

For the precise syntax, see \ref METAD

\page MTS Multiple time stepping

By setting a STRIDE different from 1, you change how frequently
an action is calculated. In the case of actions such as \ref PRINT, this just
means how frequently you dump some quantity on the disk.
Notice that variables are only computed when necessary. Thus,
if a variable is only appearing as the argument of a \ref PRINT statement with
STRIDE=10, it will be computed every 10 steps.

In a similar fashion, the STRIDE keyword can be used in a bias potential
so as to apply the bias potential every few steps.
In this case, forces from this bias potential are scaled up by
a factor equal to STRIDE.

This technique can allow your simulation to run faster if you need
the apply a bias potential on some very expensive collective variable.
Consider the following input:
\plumedfile
c1: COM ATOMS=1-1000
c2: COM ATOMS=1001-2000
d:  DISTANCE ATOMS=c1,c2
METAD ARG=d HEIGHT=1 SIGMA=0.1 BIASFACTOR=5 PACE=500
\endplumedfile
This performs a \ref METAD simulation biasing the distance between two
centers of mass. Since computing these centers requires a lot of atoms
to be imported from the MD engine, it could slow down significantly the
simulation. Notice that whereas the bias is changed every PACE=500 steps,
it is applied every STRIDE step, where STRIDE=1 by default.
The following input could lead to a significantly faster simulation at the price
of a negligible systematic error
\plumedfile
c1: COM ATOMS=1-1000
c2: COM ATOMS=1001-2000
d:  DISTANCE ATOMS=c1,c2
METAD ARG=d HEIGHT=1 SIGMA=0.1 BIASFACTOR=5 PACE=500 STRIDE=2
\endplumedfile
Similarly, the STRIDE keyword can be used with other biases (e.g. \ref RESTRAINT).

The technique is discussed in details here \cite Ferrarotti2015.
See also \subpage EFFECTIVE_ENERGY_DRIFT.

\page Multicolvar Multicolvar

Whenever you have a multicolvar action such as:

\plumedfile
COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1. D_MAX=3.0} MORE_THAN={RATIONAL R_0=6.0 NN=6 MM=12 D_0=0}
\endplumedfile

You will get a colossal speedup by specifying the D_MAX keyword in all switching functions that act on distances.
D_MAX tells PLUMED that the switching function is strictly zero if the distance is greater than this value.  As a result
PLUMED knows that it does not need to calculate these zero terms in what are essentially sums with a very large number of terms.
In fact when D_MAX is set PLUMED uses linked lists when calculating these coordination numbers, which is what 
gives you such a dramatic increase in performance.

\page Neighbour Neighbor Lists

Collective variables that can be speed up making us of neighbor lists:
- \ref COORDINATION
- \ref DHENERGY
- \ref PATHMSD

By tuning the cut-off for the neighbor list and the frequency for the recalculation of the list it is
possible to balance between accuracy and performances.

Notice that for \ref COORDINATION and \ref DHENERGY using a neighbor list could imply that a smaller
number of atoms are requested to the host MD engine. This is typically true when considering
\ref COORDINATION of a small number of atoms (e.g. a ligand) again many atoms (e.g. water).
When the neighbor list is used, only the water atoms close to the ligand will be requested at each step.

\warning
Notice that the calculation of the neighbor list is not not parallelized for \ref COORDINATION and \ref DHENERGY.
As a consequence, if you run
with many processors and/or OpenMP threads, the neighbor list might even make the calculation slower.


\page Openmp OpenMP

PLUMED is partly parallelized using OpenMP.
This should be enabled by default if your compiler supports it,
and can be disabled with `--disable-openmp`..
At runtime, you should set the environment variable
PLUMED_NUM_THREADS to the number of threads you wish to use with PLUMED.
The number of OpenMP threads can be set either by the MD code, if implemented in the patch, or generally by setting PLUMED_NUM_THREADS.
If they are not set openmp will be disabled at runtime. 

E.g., to run with gromacs you can do:
\verbatim
export PLUMED_NUM_THREADS=8
mdrun -plumed
\endverbatim

or as well

\verbatim
mdrun -plumed -ntomp 8
\endverbatim

In the first case the number of OpenMP threads used by plumed is 8 while the one used by gromacs can be 1 or something else, this is usually sub optimal.
In the second case GROMACS and plumed will use the same number of OpenMP threads.

Notice that:
- This option is likely to improve the performance, but could also slow down
  the code in some case.
- Results could be slightly different because of numerical round off and
  different order in summations. This should be harmless.
- The optimum number of threads is not necessary "all of them", nor should be
  equal to the number of threads used to parallelize MD.
- Only a few CVs are parallelized with openMP (currently, \ref COORDINATION and
  \ref DHENERGY).
- You might want to tune also the environmental variable PLUMED_CACHELINE_SIZE,
  by default 512, to set the size of cache lines on your machine. This is used
  by PLUMED to decrease the number of threads to be used in each loop so as to
  avoid clashes in memory access. This variable is expected to affect
  performance only, not results.


\page Secondary Secondary Structure

Secondary Structure collective variables (\ref ALPHARMSD, \ref PARABETARMSD and \ref ANTIBETARMSD)
can be particularly demanding if you want to calculate them for all the residues of a protein. 
This is particularly true for the calculation of beta structures.

The FIRST thing to speed up \ref PARABETARMSD and \ref ANTIBETARMSD is to use the keyword
STRANDS_CUTOFF (i.e. STRANDS_CUTOFF=1), in this way only a subset of possible fragments, the one
less than 1. nm apart, are used in the calculation.

The metric used to calculate the distance from ideal secondary structure elements can also influence 
the performances, try to use TYPE=OPTIMAL or TYPE=OPTIMAL-FAST instead of TYPE=DRMSD.

At last, try to reduce the number of residues in the calculation.

\page Lepton Making lepton library faster

In case you are using a lot of \ref CUSTOM functions or \ref switchingfunction "switching functions",
notice that these commands depend on the lepton library that is included in PLUMED.
This library replaces libmatheval since PLUMED 2.5, and by itself it is significantly faster than libmatheval.
However, you can make it even faster using a [just-in-time compiler](https://github.com/asmjit/asmjit.git).
As of PLUMED 2.6, the correct version of ASMJIT is embedded in PLUMED.
As of PLUMED 2.8, ASMJIT is enabled by default on supported architectures (X86/X64).
You can disable it at runtime setting the environment variable
`PLUMED_USE_ASMJIT`:
\verbatim
export PLUMED_USE_ASMJIT=no
\endverbatim

In some case using a custom expression is almost as fast as using a hard-coded
function. For instance, with an input that contained the following lines:
\plumedfile
c: COORDINATION GROUPA=1-108 GROUPB=1-108 R_0=1
d_fast: COORDINATION GROUPA=1-108 GROUPB=1-108 SWITCH={CUSTOM FUNC=1/(1+x2^3) R_0=1}
\endplumedfile
I (GB) obtained the following timings (on a Macbook laptop):
\verbatim
...
PLUMED: 4A  1 c                                          108     0.126592     0.001172     0.000701     0.002532
PLUMED: 4A  2 d_fast                                      108     0.135210     0.001252     0.000755     0.002623
...
\endverbatim

Notice the usage of `x2` as a variable for the switching function (see \ref switchingfunction), which
avoids an unnecessary square root calculation (this is done automatically by the hard-coded switching functions
when you use only even powers). The asmjit calculation (`d_fast`) takes less than 10% more than the hard-coded
one (`c`).

\page Time Time your Input

Once you have prepared your plumed input file you can run a test simulation, or use driver, 
to see which collective variable, function, bias or analysis is consuming more time and can 
thus be the target for a different definition (use less atoms, change relevant parameters,
or just use something else)

To have an accurate timing of your input you can use the \ref DEBUG DETAILED_TIMERS.
  
\page FUNNELMOD Funnel-Metadynamics (FM)

<!-- 
description: a collective variable and a bias action necessary to perform Funnel-Metadynamics on Molecular Dynamics simulations
authors: Stefano Raniolo, Vittorio Limongelli
reference: \cite limongelli2013funnel \cite raniolo2020ligand
-->

## Overview
FM is a combination of Metadynamics bias potential \cite metad with a funnel-shape restraint potential applied to the target structure of a binding interaction. 
The latter is composed of a cone restraint, which covers the ligand binding site, and a cylindrical one that heads towards the solvent \cite limongelli2013funnel. 
When inside the funnel volume, the ligand does not feel any restraint potential, proceeding as regular Metadynamics.
Upon reaching the boundaries of the funnel, a repulsive bias is applied forcing the ligand to remain in the allowed funnel space. 
The result is an acceleration in the sampling of the binding/unbinding process, leading to a swift convergence of the calculation and a well-defined binding free-energy surface.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=funnel' to your './configure' command when building PLUMED to enable these features.

## Usage
This module is a direct evolution of the original FM \cite limongelli2013funnel since it incorporates an alignment function that removes the necessity to block the target macromolecule in the simulation box.

The user can follow a comprehensive protocol \cite raniolo2020ligand, which will help in all stages of the simulation, including pre- and post-processing.
An example of input file can be found on <a href="https://www.plumed-nest.org/eggs/19/039/">FUNNEL-NEST's webpage</a>

## Module Contents

The funnel module is composed of a collective variable that calculates the position of a ligand with respect to a line and a potential that creates a funnel-shape restraint centered on the line (\ref FUNNEL_PS and \ref FUNNEL, respectively).

- \subpage funnel_cv
- \subpage funnel_bias

\page funnel_cv CV documentation

The following list contains descriptions of the collective variables developed for the PLUMED-FUNNEL module. They should be used in combination with the funnel-shaped restraint potential and Metadynamics to enable Funnel-Metadynamics.

@FUNNELMOD_COLVAR@

\page funnel_bias Bias Documentation

The following list contains descriptions of biases developed for the PLUMED-FUNNEL module. They should be used in combination with the collective variable to calculate the position relative to the funnel-shape restraint potential and Metadynamics to enable Funnel-Metadynamics.

@FUNNELMOD_BIAS@
\page ANNMOD ANN (Artificial Neural Network) function

<!-- 
description: ANN (Artificial Neural Network) function
authors: Wei Chen and Andrew Ferguson 
reference:  
-->

## Overview 

This is plumed ANN function (annfunc) module.  It implements `ANN` class, which is a subclass of `Function` class.  `ANN` class takes multi-dimensional arrays as inputs for a fully-connected feedforward neural network with specified neural network weights and generates corresponding outputs.  The `ANN` outputs can be used as collective variables, inputs for other collective variables, or inputs for data analysis tools.  

## Installation

This module is not installed by default. Add '\-\-enable-modules=annfunc' to your './configure' command when building PLUMED to enable these features.

## Usage

Currently, all features of the ANNfunc module are included in a single ANNfunc collective variable: \ref ANN

## Module Contents
- \subpage ANNMODFunction

\page ANNMODFunction Functions Documentation

The following list contains descriptions of functions developed for the PLUMED-ANNfunc module. They can be used in combination with other actions outside of the ANNfunc module.

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage ANN </td> <td>Calculates the ANN-function.</td> </tr>
</table>
\page SASAMOD SASA collective variable

<!-- 
description: Solvent Accessible Surface Area collective variable (SASA)
authors: Andrea Arsiccio
reference: \cite Hasel1988 \cite Weiser1999 \cite Arsiccio-SASA-2021 
-->

## Overview

This SASA module contains methods for the calculation of the solvent accessible surface area (SASA) of proteins using either the fast algorithm by Hasel et al. \cite Hasel1988 or the LCPO algorithm \cite Weiser1999. This module can be used to include the SASA as a collective variable in metadynamics simulations, and also for implicit solvent simulations as described in \cite Arsiccio-SASA-2021.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=sasa' to your './configure' command when building PLUMED to enable these features.

## Usage
Currently, all features of the SASA module are included in two SASA functions: \ref SASA_HASEL \ref SASA_LCPO

## Module Contents
- \subpage SASAMODColvar

\page SASAMODColvar CVs Documentation

The following list contains descriptions of biases developed for the PLUMED-SASA module. They can be used in combination with other biases outside of the SASA module.

@SASAMOD_COLVAR@
\page Syntax Getting Started 

To run PLUMED you need to provide one input file.  In this file you specify what it
is that PLUMED should do during the course of the run.  Typically this will involve calculating 
one or more collective variables, perhaps calculating a function of these CVs
 and then doing some analysis of values of your collective variables/functions or running
some free energy method. A very brief introduction to the syntax used in the PLUMED input file
is provided in this <a href="http://www.youtube.com/watch?v=PxJP16qNCYs"> 10-minute video </a>.

Within this input file every line is an instruction for PLUMED to perform some particular action.  This could be
 the calculation of a colvar, an occasional analysis of the trajectory or a biasing of the dynamics.  The first
word in these lines specify what particular action is to be performed.  This is then followed by a number of keywords
which provide PLUMED with more details as to how the action is to be performed.  These keywords are either single words
(in which they tell PLUMED to do the calculation in a particular way - for example NOPBC tells PLUMED to not use the periodic
boundary conditions when calculating a particular colvar) or they can be words followed by an equals sign and a comma separated 
list _with no spaces_ of numbers or characters (so for example ATOMS=1,2,3,4 tells PLUMED to use atom numbers 1,2,3 and 4 in 
the calculation of a particular colvar).
The reason why spaces are not admitted is that PLUMED should be able to understand when the list of atoms
ended and a new keyword should be expected. 
Space separated lists can be used instead of comma separated list if the entire list
is enclosed in curly braces (e.g. ATOMS={1 2 3 4}).  Please note that you can split commands over multiple lines by using
\ref ContinuationLines. 

The most important of these keywords is the label keyword as it is only by using these labels that we can pass data 
from one action to another.  As an example if you do:

\plumedfile
DISTANCE ATOMS=1,2
\endplumedfile

Then PLUMED will do nothing other than read in your input file.  In contrast if you do:

\plumedfile
DISTANCE ATOMS=1,2 LABEL=d1
PRINT ARG=d1 FILE=colvar STRIDE=10
\endplumedfile

then PLUMED will print out the value of the distance between atoms 1 and 2 every 10 steps to the file colvar as you have told
PLUMED to take the value calculated by the action d1 and to print it. You can use any character string to label your actions
as long as it does not begin with the symbol \@.  Strings beginning with \@ are used by within PLUMED to reference special, 
code-generated groups of atoms and to give labels to any Actions for which the user does not provide a label in the input. 

Notice that if a word followed by a column is added at the beginning of the line (e.g. pippo:), PLUMED automatically
removes it and adds an equivalent label (LABEL=pippo).
Thus, a completely equivalent result can be obtained with the following shortcut:
\plumedfile
d1: DISTANCE ATOMS=1,2
PRINT ARG=d1 FILE=colvar STRIDE=10
\endplumedfile

Also notice that all the actions can be labeled, and that many actions besides normal collective variables can define
one or more value, which can be then referred using the corresponding label.

Actions can be referred also with POSIX regular expressions (see \ref Regex) if regex library is available on your system
and detected at configure time.
You can also add \ref comments to the input or set up your input over multiple files and then create a composite input by
\ref includes.

More information on the input syntax as well as details on the the various trajectory
analysis tools that come with PLUMED are given in: 

- \ref colvarintro tells you about the ways that you can calculate functions of the positions of the atoms.
- \ref Analysis tells you about the various forms of analysis you can run on trajectories using PLUMED.
- \ref Bias tells you about the methods that you can use to bias molecular dynamics simulations with PLUMED.

\section units Plumed units
By default the PLUMED inputs and outputs quantities in the following units:

- Energy - kJ/mol
- Length - nanometers
- Time - picoseconds

Unlike PLUMED 1 the units used are independent of the MD engine you are using.  If you want to change these units you can do this using the 
\subpage UNITS keyword. 

\page colvarintro Collective Variables

Chemical systems contain an enormous number atoms, which, in most cases makes it simply impossible for
us to understand anything by monitoring the atom positions directly.  Consequently,
we introduce Collective variables (CVs) that describe the chemical processes we are
interested in and monitor these simpler quantities instead.  These CVs are used in many of the methods
implemented in PLUMED - there values can be monitored using \ref PRINT, \ref Function of them can be calculated
or they can be analyzed or biased using the \ref Analysis and \ref Bias "Biasing" methods implemented in PLUMED.
Before doing any of these things however we first have to tell PLUMED how to calculate them.

The simplest collective variables that are implemented in PLUMED take in a
set of atomic positions and output one or multiple scalar CV values.  Information on these variables is given on the page entitled 
\ref Colvar while information as to how sets of atoms can be selected
can be found in the pages on \ref Group.  Please be aware that PLUMED contains implementations of many other collective variables 
but that the input for these variables may be less transparent when it is first encountered.
In particular, the page on \ref dists describes the various ways that you can calculate the distance from a particular reference
configuration.  So you will find instructions on how to calculate the RMSD distance from the folded state of a protein here.
Meanwhile, the page on \ref Function describes the various functions of collective variables that can be used in the
code.  This is a very powerful feature of PLUMED as you can use the \ref Function commands to calculate any function or 
combination of the simple collective variables listed on the page \ref Colvar.  Lastly the page on \ref mcolv describes MultiColvars.  
MultiColvars allow you to use many different colvars and allow us to
implement all these collective variables without a large amount of code.  For some things (e.g.
\ref DISTANCES GROUPA=1 GROUPB=2-100 LESS_THAN={RATIONAL R_0=3}) there are more computationally efficient options available in plumed
(e.g. \ref COORDINATION).  However, MultiColvars are worth investigating as they provide a flexible syntax for many quite-complex CVs.

- \subpage Group
- \subpage Colvar
- \subpage dists
- \subpage Function
- \subpage mcolv
- \subpage contactmatrix

\page Colvar CV Documentation

The following list contains descriptions of a number of the colvars that are currently implemented in PLUMED.

@COLVAR@

\page dists Distances from reference configurations

One colvar that has been shown to be very successful in studying protein folding is the distance between the instantaneous configuration
and a reference configuration - often the structure of the folded state.  When the free energy of a protein is shown as a function
of this collective variable there is a minima for low values of the CV, which is due to the folded state of the protein.  There is 
then a second minima at higher values of the CV, which is the minima corresponding to the unfolded state.

A slight problem with this sort of collective variable is that there are many different ways of calculating the distance from a 
particular reference structure.  The simplest - adding together the distances by which each of the atoms has been translated in
going from the reference configuration to the instantaneous configuration - is not particularly sensible.  A distance calculated
in this way does not neglect translation of the center of mass of the molecule and rotation of the frame of reference.  A common practice
is thus to remove these components by calculating the \ref RMSD distance between the reference and instantaneous configurations.
This is not the only way to calculate the distance, however.  One could also calculate the total amount by which a large number 
of collective variables change in moving from the reference to the instantaneous configurations.  One could even combine RMSD distances
with the amount the collective variables change.  A full list of the ways distances can be measured in PLUMED is given below:

@DCOLVAR@

These options for calculating distances are re-used in a number of places in the code.  For instance they are used in some of the 
analysis algorithms that are implemented in PLUMED and in \ref PATH collective variables. 
Notice that most of these actions read the reference configuration from a PDB file. Be sure
you understand how to format properly a PDB file to use used in PLUMED (see \ref pdbreader).

\page mcolv MultiColvar 

Oftentimes, when you do not need one of the collective variables described elsewhere in the manual, what you want instead is a 
function of a distribution of collective variables of a particular type.  In other words, you would like to calculate a
function something like this:
\f[
s = \sum_i g[f(\{X\}_i)]
\f]
In this expression \f$g\f$ is a function that takes in one argument and \f$f\f$ is a function that takes a set of atomic positions
as argument. The symbol \f$\{X\}_i\f$ is used to indicate the fact that the function \f$f\f$ is evaluated for a number of different
sets of atoms.  If you would just like to output the values of all the various \f$f\f$ functions you should use the command \ref DUMPMULTICOLVAR

This functionality is useful if you need to calculate a minimum distance or the number of coordination numbers greater than a 3.0.  
To avoid duplicating the code to calculate an angle or distance many times and to make it easier to implement very complex collective 
variables PLUMED provides these sort of collective variables using so-called MultiColvars.  MultiColvars are named in this way because a single
PLUMED action can be used to calculate a number of different collective variables.  For instance the \ref DISTANCES
action can be used to calculate the minimum distance, the number of distances less than a certain value, the number of
distances within a certain range... A more detailed introduction to multicolvars is provided in this 
<a href="http://www.youtube.com/watch?v=iDvZmbWE5ps">10-minute video</a>. Descriptions of the various multicolvars
that are implemented in PLUMED 2 are given below: 

@MCOLVAR@  

To instruct PLUMED to calculate a multicolvar you give an instruction that looks something like this:

\verbatim
NAME <atoms involved> <parameters> <what am I calculating> TOL=0.001 LABEL=label
\endverbatim

Oftentimes the simplest way to specify the atoms involved is to use multiple instances of the ATOMS keyword 
i.e. ATOMS1, ATOMS2, ATOMS3,...  Separate instances of the quantity specified by NAME are then calculated for 
each of the sets of atoms.  For example if the command issued contains the following:

\plumedfile
DISTANCES ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
\endplumedfile

The distances between atoms 1 and 2, atoms 3 and 4, and atoms 5 and 6 are calculated. Obviously, generating 
this sort of input is rather tedious so short cuts are also available many of the collective variables. 
These are described on the manual pages for the actions.
 
After specifying the atoms involved you sometimes need to specify some parameters that required in the 
calculation.  For instance, for \ref COORDINATIONNUMBER - the number of atoms in the first coordination
sphere of each of the atoms in the system - you need to specify the parameters for a \ref switchingfunction
that will tell us whether or not an atom is in the first coordination sphere.  Details as to how to do this
are provided on the manual pages.  

One of the most important keywords for multicolvars is the TOL keyword.  This specifies that terms in sums 
that contribute less than a certain value can be ignored.  In addition, it is assumed that the derivative
with respect to these terms are essentially zero.  By increasing the TOL parameter you can increase the speed 
of the calculation.  Be aware, however, that this increase in speed is only possible because you are lowering the 
accuracy with which you are computing the quantity of interest.

Once you have specified the base quantities that are to be calculated from the atoms involved and any parameters
you need to specify what function of these base quantities is to be calculated.  For most multicolvars you can calculate 
the minimum, the number less than a target value, the number within a certain range, the number more than a target
value and the average value directly.  

\section multicolvarfunction MultiColvar functions

It is possible to use multicolvars to calculate complicated collective variables by exploiting the fact that the output
from one multicolvar can be used as input to a second multicolvar.  One simple way of exploiting this functionality is to
filter the atoms based on the value they have for a symmetry function.  For example you might want to consider only those
atoms that with a \ref COORDINATIONNUMBER higher that a certain threshold when calculating some particularly expensive symmetry
function such as \ref Q6.  The following methods can thus all be used to filter the values of multicolvars in this way:

@MFILTERS@

An alternative way of filtering atoms is to consider only those atoms in a particular part of the simulation box.  This can be
done by exploiting the following methods

@VOLUMES@

The idea with these methods is that function of the form:
\f[
s = \sum_i w(\{X\}_i) g[f(\{X\}_i)]
\f]
can be evaluated where once again \f$g\f$ is a function with one argument and \f$g\f$ is a function of a set of atomic positions.  
The difference from the more general function described earlier is that we now have a weight \f$w\f$ which is again a function of the
atomic positions.  This weight varies between zero and one and it is this weight that is calculated in the list of filtering methods
and volume methods described in the lists above.  

In addition to these volume and filtering methods it is also possible to calculate the local average of a quantities in the manner 
described in \cite dellago-q6 using the \ref LOCAL_AVERAGE method.  Furthermore, in many cases \ref Q6, \ref MOLECULES and 
\ref PLANES the symmetry function being evaluated is a vector.  You can thus construct a variety of novel collective variables by
taking dot products of vectors on adjacent atoms as described below: 

@MCOLVARF@ 

The final set of functions that you can apply on multicolvars are functions that transform all the colvars calculated using a 
multicolvar using a function.  This can be useful if you are calculating some complicated derived quantity of some simpler 
quantity.  It is also useful if you are calculating a Willard Chandler surface or a histogram.  The actions that you can use to 
perform these transforms are:

@MTRANSFORMS@

\section multicolvarbias MultiColvar bias

There may be occasions when you want add restraints on many collective variables. For instance if you are studying a cluster
you might want to add a wall on the distances between each of the atoms and the center of mass of the cluster in order to
prevent the cluster subliming.  Alternatively, you may wish to insist that a particular set of atoms in your system all have a 
coordination number greater than 2.  You can add these sorts of restraints by employing the following biases, which all act 
on the set of collective variable values calculated by a multicolvar.  So for example the following set of commands:

\plumedfile
COM ATOMS=1-20 LABEL=c1
DISTANCES GROUPA=c1 GROUPB=1-20 LABEL=d1
UWALLS DATA=d1 AT=2.5 KAPPA=0.2 LABEL=sr
\endplumedfile

creates the aforementioned set of restraints on the distances between the 20 atoms in a cluster and the center of mass of the cluster.

The list of biases of this type are as follows:

@MCOLVARB@

Notice that (in theory) you could also use this functionality to add additional terms to your force field or to implement your 
force field.

\section usingbase Extracting all the base quantities

There may be occasions where you want to get information on all the individual colvar values that you have calculated.
For example you might want to output the values of all the coordination numbers calculated by a \ref COORDINATIONNUMBER 
action.  You can thus use the following command to extract this sort of information, \ref DUMPMULTICOLVAR.

\page contactmatrix Exploiting contact matrices

A contact matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether or not the \f$i\f$th
and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  There are various ways of defining
whether a pair of atoms/molecules are adjacent or not.  For example we can say two atoms are adjacent if the distance between
them is less than some cutoff.  Alternatively, if we have a have a pair of molecules, we might state they are adjacent if their
centers of mass are within a certain cutoff and if the two molecules have the same orientation.  Two electronegative atoms
might be said to be adjacent if there is a hydrogen bond between them.  For these reasons then PLUMED contains all of the 
following methods for calculating an adjacency matrix 

@MATRIX@

Once you have calculated an adjacency matrix you can then perform any one of the following operations on this object in order
to reduce it to a scalar number or a set of connected components.

@MATRIXF@

If the function you have chosen reduces your contact matrix to a set of connected components you then need a method to convert 
these connected components into a scalar number or to output this information to a file.  The various things that you can do
with a set of connected components are listed below:

@CONCOMP@

\page EABFMOD Extended-System Adaptive Biasing Force 

<!--
description: Methods for performing eABF or DRR method to calculate PMF along CVs
authors: Haochuan Chen, Haohao Fu
reference: \cite Chen2018 \cite Lelievre2007 \cite Lesage2016 \cite Fu2016
-->

## Overview

This module contains the eABF/DRR method to do free energy calculation or enhance sampling along CVs.

## Installation

This module is not installed by default and depends on the boost serialization module. Please make sure the boost serialization library is compiled and installed in your system before trying to compile this module. Add '\-\-enable-modules=drr \-\-enable-boost_serialization' to your './configure' command when building PLUMED to enable these features.

## Usage

Please read \ref drr_tool and \ref DRR for more information.

## Module Contents

- \subpage EABFMODBias
- \subpage EABFMODCLTools

\page EABFMODBias Biases Documentation

The following list contains descriptions of biases developed for the eABF module. They can be used in combination with other biases outside of the eABF module.

@EABFMOD_BIAS@

\page EABFMODCLTools Command Line Tools

The following list contains the command line tools available in the eABF module.

@EABFMOD_TOOLS@
\page mymodules List of modules

The functionality in PLUMED 2 is divided into a small number of modules.  Some
users may only wish to use a subset of the functionality available within the 
code while others may wish to use some of the more complicated features that are available.
For this reason the plumed source code is divided into modules, which users can
activate or deactivate to their hearts content.  

You can activate a module at configure time using the keyword `--enable-modules`.
For example:
\verbatim
./configure --enable-modules=modulename
\endverbatim
will enable module called modulename. A module that is on by default can be disabled
using the following syntax
\verbatim
./configure --enable-modules=-modulename
\endverbatim
To enable or disable multiple modules one should provide them as a : separated
list. Notice that `+modulename` and `modulename` both activate the module, whereas
`-modulename` deactivates it. E.g.
\verbatim
./configure --enable-modules=+crystallization:-colvar
\endverbatim
will disable the colvar module and enable the crystallization module.
Also notice that `:` can be omitted when using `+` or `-`. Thus, the same can be obtained
with
\verbatim
./configure --enable-modules=+crystallization-colvar
\endverbatim

If you repeat the `--enable-modules` keyword only the last instance will be used. Thus
`./configure --enable-modules=crystallization --enable-modules=-colvar` will _not_ do what you expect!

There are also some shortcuts available:
- `./configure --enable-modules=all` to enable all optional modules. This includes the maximum number of features in PLUMED,
including modules that might not be properly functional.
- `./configure --enable-modules=none` or `./configure --disable-modules` to disable all optional modules. This produces a minimal
PLUMED which can be used as a library but which has no command line tools and no collective variables or biasing methods.
- `./configure --enable-modules=reset` or `./configure --enable-modules` to enable the default modules.

The two kinds of syntax can be combined and, for example, `./configure --enable-modules=none:colvar` will result
in a PLUMED with all the modules disabled with the exception of the colvar module.

Some modules are active by default in the version of PLUMED 2 that you download from 
the website while others are inactive.  The following lists all of the modules that
are available in plumed and tells you whether or not they are active by default.

@MODULES@

Until PLUMED 2.2, it was also possible to switch on or off modules by adding files
in the plumed2/src directory. Since PLUMED 2.3 this is discouraged, since any choice made
in this manner will be overwritten next time `./configure` is used.

\page VES Variationally Enhanced Sampling (VES code) 

<!-- 
description: Module that implements enhanced sampling methods based on Variationally Enhanced Sampling
authors: Omar Valsson
reference: \cite Valsson-PRL-2014
-->

The VES code is a module for PLUMED that implements enhanced sampling methods
based on _Variationally Enhanced Sampling_ (VES) \cite Valsson-PRL-2014.
The VES code is developed by [Omar Valsson](http://www.valsson.info), 
see the [homepage of the VES code](http://www.ves-code.org) for further information.

The VES code is an optional module that needs to be enabled when configuring the
compilation of PLUMED by using the '\-\-enable-modules=ves' 
(or '\-\-enable-modules=all') flag when running the 'configure' script. 

In the \ref ves_tutorials "tutorials" you can learn how to use the methods 
implemented in the VES code.

The various components of the VES code module are listed and described in the following sections 

- \subpage ves_biases
- \subpage ves_basisf
- \subpage ves_targetdist
- \subpage ves_optimizer
- \subpage ves_utils
- \subpage ves_cltools
- \subpage ves_tutorials


\page ves_biases Biases

The following list contains the biases available in the VES code.

@VES_BIAS@

\page ves_basisf Basis functions

The following list contains the one-dimensional basis functions available in the VES code.

@VES_BASISF@


\page ves_targetdist Target Distributions

The following list contains the target distributions available in the VES code.

@VES_TARGETDIST@


\page ves_optimizer Optimizers

The following list contains the optimizers available in the VES code.

@VES_OPTIMIZER@


\page ves_utils Utilities

The following list contains various utilities available in the VES code. 

@VES_UTILS@



\page ves_cltools Command Line Tools

The following list contains the command line tools available in the VES code.

@VES_TOOLS@



\page ves_tutorials Tutorials

The following tutorials are available for the VES code. 

\subpage ves_tutorial_lugano_2017

@VES_TUTORIALS@




\page ves_tutorial_lugano_2017 MARVEL-VES School February 2017

\image html ves-lugano2017-logo.png  width=800px

Tutorials from the [MARVEL School on Variationally Enhanced Sampling]
(https://sites.google.com/site/vesschool2017/home) that was held in
Lugano, February 14-17, 2017.

\par Suggested readings

Metadynamics:

[Enhancing Important Fluctuations: Rare Events and Metadynamics from a Conceptual Viewpoint](https://doi.org/10.1146/annurev-physchem-040215-112229), Annual Reviews in Physical Chemistry 2016



Variationally Enhanced Sampling:

[Variational Approach to Enhanced Sampling and Free Energy Calculations](https://doi.org/10.1103/PhysRevLett.113.090601), Physical Review Letters 2014

[Variationally Optimized Free-Energy Flooding for Rate Calculation](https://doi.org/10.1103/PhysRevLett.115.070601), Physical Review Letters 2015



\par Tuesday February 14

\ref marvel-1 "Tutorial 1": Introduction to PLUMED and analyzing molecular simulations

\par Wednesday February 15

\ref ves-lugano2017-metad "Tutorial 2": Biasing with metadynamics

\ref ves-lugano2017-ves1 "Tutorial 3": Biasing with variationally enhanced sampling

\par Thursday February 16

\ref ves-lugano2017-ves2 "Tutorial 4": Further on variationally enhanced sampling

Tutorial 5: Advanced collective variables
- \ref marvel-2 "Path CVs"
- \ref belfast-10 "Multicolvar"
- \ref belfast-3 "Dimensionality reduction"

\par Friday February 17

\ref ves-lugano2017-kinetics "Tutorial 6": Obtaining kinetics from molecular simulations
\page Analysis Analysis

PLUMED can be used to analyze trajectories either on the fly during an MD run or via
post processing a trajectory using \ref driver.  A molecular dynamics trajectory is in essence an ordered 
set of configurations of atoms.  Trajectory analysis algorithms are methods that allow us to extract meaningful 
information from this extremely high-dimensionality information.  In extracting this information much of the 
information in the trajectory will be discarded and assumed to be irrelevant to the problem at hand.  For example, 
when we calculate a histogram from a trajectory we throw away all information on the order the frames were visited during the
trajectory.  We instead opt to display a time average that shows the parts of configuration space that were  
visited most frequently.  There are many situations in which this is a reasonable thing to do as we know that
time averages are equivalent to ensemble averages in the long timescale limit and that these average probabilities
of being in different parts of configuration space, \f$P(s)\f$, are thus related to the underlying free
energy, \f$F(s)\f$, via:
\f[
F(s) = - k_B T \ln P(s)
\f]
About the simplest form of analysis 
that PLUMED can perform involves printing information to a file.  PLUMED can output
various different kinds of information to files as described below:

@PRINTANALYSIS@  

The \ref UPDATE_IF action allows you to do more complex things using the above print
commands. As detailed in the documentation for \ref UPDATE_IF when you put any of the above 
actions within an UPDATE_IF block then data will only be output to the file if colvars
are within particular ranges.  In other words, the above printing commands, in tandem 
with \ref UPDATE_IF, allow you to identify the frames in your trajectory that satisfy
some particular criteria and output information on those frames only.

Another useful command is the \ref COMMITTOR command. 
As detailed in the documentation for \ref COMMITTOR this command tells PLUMED (and the underlying 
MD code) to stop the calculation one some criteria is satisfied, alternatively one can use it to keep
track of the number of times a criteria is satisfied.

A number of more complicated forms of analysis can be performed that take a number of frames from 
the trajectory as input.  In all these commands the STRIDE keyword is used to tell PLUMED how 
frequently to collect data from the trajectory.  In all these methods the output from the analysis
is a form of ensemble average.  If you are running with a bias it is thus likely that you may want 
to reweight the trajectory frames in order to remove the effect the bias has on the static behavior
of the system.  The following methods can thus be used to calculate weights for the various trajectory
frames so that the final ensemble average is an average for the canonical ensemble at the appropriate 
temperature.

\section analysisbias Reweighting and Averaging

@REWEIGHTING@

You can then calculate ensemble averages using the following actions.

@GRIDCALC@

For many of the above commands data is accumulated on the grids.  These grids can be further 
analyzed using one of the actions detailed below at some time.  

@GRIDANALYSIS@

As an example the following set of commands instructs PLUMED to calculate the distance between 
atoms 1 and 2 for every fifth frame in the trajectory and to accumulate a histogram from this data
which will be output every 100 steps (i.e. when 20 distances have been added to the histogram).

\plumedfile
x: DISTANCE ATOMS=1,2
h: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1 STRIDE=5
DUMPGRID GRID=h FILE=histo STRIDE=100 
\endplumedfile

It is important to note when using commands such as the above the first frame in the trajectory is assumed 
to be the initial configuration that was input to the MD code. It is thus ignored.  Furthermore, if you are 
running with driver and you would like to analyze the whole trajectory (without specifying its length) 
and then print the result you simply call \ref DUMPGRID (or any of the commands above) without a STRIDE 
keyword as shown in the example below. 

\plumedfile
x: DISTANCE ATOMS=1,2
h: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1 STRIDE=5
DUMPGRID GRID=h FILE=histo 
\endplumedfile

Please note that even with this calculation the first frame in the trajectory is ignored when computing the 
histogram.

Notice that all the commands for calculating smooth functions described above calculate some sort of 
average.  There are two ways that you may wish to average the data in your trajectory:

- You might want to calculate block averages in which the first \f$N\f$N frames in your trajectory are
averaged separately to the second block of \f$N\f$ frames.  If this is the case you should use the 
keyword CLEAR in the input to the action that calculates the smooth function.  This keyword is used to 
specify how frequently you are clearing the stored data.

- You might want to calculate an accumulate an average over the whole trajectory and output the average
accumulated at step \f$N\f$, step \f$2N\f$...  This is what PLUMED does by default so you do not need to 
use CLEAR in this case.

\section diag Diagnostic tools

PLUMED has a number of diagnostic tools that can be used to check that new Actions are working correctly: 

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \ref DUMPFORCES </td> <td>Dump the force acting on one of a values in a file.  </td> </tr>
<tr> <td width=5%> \ref DUMPDERIVATIVES </td> <td>Dump the derivatives with respect to the input parameters for one or more objects (generally CVs, functions or biases).</td> </tr>
<tr> <td width=5%> \ref DUMPMASSCHARGE </td> <td>Dump masses and charges on a selected file.</td> </tr>
<tr> <td width=5%> \ref DUMPPROJECTIONS </td> <td>Dump the derivatives with respect to the input parameters for one or more objects (generally CVs, functions or biases).</td> </tr>
</table>

These commands allow you to test that derivatives and forces are calculated correctly
within colvars and functions.  One place where this is very useful is when you are testing whether or
not you have implemented the derivatives of a new collective variables correctly.  So for example if
we wanted to do such a test on the distance CV we would employ an input file something like this:

\plumedfile
d1: DISTANCE ATOMS=1,2
d1n: DISTANCE ATOMS=1,2 NUMERICAL_DERIVATIVES
DUMPDERIVATIVES ARG=d1,d1n FILE=derivatives
\endplumedfile

The first of these two distance commands calculates the analytical derivatives of the distance
while the second calculates these derivatives numerically.  Obviously, if your CV is implemented
correctly these two sets of quantities should be nearly identical.

\section storing Storing data for analysis

All the analysis methods described in previous sections accumulate averages or output diagnostic information on the fly.
That is to say these methods calculate something given the instantaneous positions of the atoms or the instantaneous 
values of a set of collective variables.  Many methods (e.g. dimensionality reduction and clustering) will not work like 
this, however, as information from multiple trajectory frames is required at the point when the analysis is performed.  In other
words the output from these types of analysis cannot be accumulated one frame at time.  When using these methods you must therefore
store trajectory frames for later analysis.  You can do this storing by using the following action:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage COLLECT_FRAMES </td> <td> Collect and store trajectory frames for later analysis with one of the methods detailed below. </td> </tr>
</table>  

\section dissimilaritym Calculating dissimilarity matrices

One of the simplest things that we can do if we have stored a set of trajectory frames using \ref COLLECT_FRAMES is we can calculate the dissimilarity between 
every pair of frames we have stored.  When using the \ref dimred "dimensionality reduction" algorithms described in 
the sections that follow the first step is to calculate this matrix.  Consequently, within PLUMED the following 
command will collect the trajectory data as your simulation progressed and calculate the dissimilarities: 

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage EUCLIDEAN_DISSIMILARITIES </td> <td> Calculate the matrix of dissimilarities between a trajectory of atomic configurations. </td> </tr>
</table>

By exploiting the functionality described in \ref dists you can calculate these dissimilarities in
a wide variety of different ways (e.g. you can use \ref RMSD, or you can use a collection of collective variable
values see \ref TARGET).  If you wish to view this dissimilarity information you can print these quantities 
to a file using:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage PRINT_DISSIMILARITY_MATRIX </td> <td> Print the matrix of dissimilarities between a trajectory of atomic configurations. </td> </tr>
</table>

In addition, if PLUMED does not calculate the dissimilarities you need you can read this information from an 
external file

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage READ_DISSIMILARITY_MATRIX </td> <td> Read a matrix of dissimilarities between a trajectory of atomic configurations from a file. </td> </tr>
</table>
 
N.B. You can only use the two commands above when you are doing post-processing.  

\section landmarks Landmark Selection

Many of the techniques described in the following sections are very computationally expensive to run on large trajectories.
A common strategy is thus to use a landmark selection algorithm to pick a particularly-representative subset of trajectory
frames and to only apply the expensive analysis algorithm on these configurations.  The various landmark selection algorithms
that are available in PLUMED are as follows

@LANDMARKS@

In general most of these landmark selection algorithms must be used in tandem with a \ref dissimilaritym "dissimilarity matrix" object as as follows:

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
data: COLLECT_FRAMES ARG=d1,d2,d3 STRIDE=1
ss1: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=data 
ll2: LANDMARK_SELECT_FPS USE_OUTPUT_DATA_FROM=ss1 NLANDMARKS=300
OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=ll2 FILE=mylandmarks
\endplumedfile

When landmark selection is performed in this way a weight is ascribed to each of the landmark configurations.  This weight is
calculated by summing the weights of all the trajectory frames in each of the landmarks Voronoi polyhedron 
(https://en.wikipedia.org/wiki/Voronoi_diagram).  The weight of each trajectory frame is one unless you are reweighting using the
formula described in the \ref analysisbias to counteract the fact of a simulation bias or an elevated temperature.  If you are reweighting
using these formula the weight of each of the points is equal to the exponential term in the numerator of these expressions.

\section dimred Dimensionality Reduction

Many dimensionality reduction algorithms work in a manner similar to the way we use when we make maps. You start with distances 
between London, Belfast, Paris and Dublin and then you try to arrange points on a piece of paper so that the (suitably transformed) 
distances between the points in your map representing each of those cities are related to the true distances between the cities.  
Stating this more mathematically MDS endeavors to find an <a href="http://en.wikipedia.org/wiki/Isometry">isometry</a> 
between points distributed in a high-dimensional space and a set of points distributed in a low-dimensional plane.  
In other words, if we have \f$M\f$ \f$D\f$-dimensional points, \f$\mathbf{X}\f$, 
and we can calculate dissimilarities between pairs them, \f$D_{ij}\f$, we can, with an MDS calculation, try to create \f$M\f$ projections, 
\f$\mathbf{x}\f$, of the high dimensionality points in a \f$d\f$-dimensional linear space by trying to arrange the projections so that the 
Euclidean distances between pairs of them, \f$d_{ij}\f$, resemble the dissimilarities between the high dimensional points.  In short we minimize:

\f[
\chi^2 = \sum_{i \ne j} w_i w_j \left( F(D_{ij}) - f(d_{ij}) \right)^2
\f]

where \f$F(D_{ij})\f$ is some transformation of the distance between point \f$X^{i}\f$ and point \f$X^{j}\f$ and \f$f(d_{ij})\f$ is some transformation
of the distance between the projection of \f$X^{i}\f$, \f$x^i\f$, and the projection of \f$X^{j}\f$, \f$x^j\f$.  \f$w_i\f$ and \f$w_j\f$ are the weights
of configurations \f$X^i\f$ and \f$^j\f$ respectively.  These weights are calculated using the reweighting and Voronoi polyhedron approaches described in
previous sections.  A tutorial on dimensionality reduction and how it can be used to analyze simulations can be found in the tutorial \ref lugano-5 and in 
the following <a href="https://www.youtube.com/watch?v=ofC2qz0_9_A&feature=youtu.be" > short video.</a>

Within PLUMED running an input to run a dimensionality reduction algorithm can be as simple as:

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
data: COLLECT_FRAMES STRIDE=1 ARG=d1,d2,d3
ss1: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=data 
mds: CLASSICAL_MDS USE_OUTPUT_DATA_FROM=ss1 NLOW_DIM=2
\endplumedfile

Where we have to use the \ref EUCLIDEAN_DISSIMILARITIES action here in order to calculate the matrix of dissimilarities between trajectory frames.
We can even throw some landmark selection into this procedure and perform

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
data: COLLECT_FRAMES STRIDE=1 ARG=d1,d2,d3
matrix: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=data
ll2: LANDMARK_SELECT_FPS USE_OUTPUT_DATA_FROM=matrix NLANDMARKS=300
mds: CLASSICAL_MDS USE_OUTPUT_DATA_FROM=ll2 NLOW_DIM=2
osample: PROJECT_ALL_ANALYSIS_DATA USE_OUTPUT_DATA_FROM=matrix PROJECTION=mds
\endplumedfile

Notice here that the final command allows us to calculate the projections of all the non-landmark points that were collected by the action with
label matrix.

Dimensionality can be more complicated, however, because the stress function that calculates \f$\chi^2\f$ has to optimized rather carefully using
a number of different algorithms.  The various algorithms that can be used to optimize this function are described below

@DIMRED@

\section output Outputting the results from analysis algorithms

The following methods are available for printing the result output by the various analysis algorithms:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage OUTPUT_ANALYSIS_DATA_TO_COLVAR </td> <td> Output the results from an analysis using the PLUMED colvar file format. </td> </tr>
<tr> <td width=5%> \subpage OUTPUT_ANALYSIS_DATA_TO_PDB </td> <td> Output the results from an analysis using the PDB file format.</td> </tr>
</table>

Using the above commands to output the data from any form of analysis is important as <b> the STRIDE with which you output the data to a COLVAR or PDB file
controls how frequently the analysis is performed on the collected data </b>.  If you specified no stride on the output lines then PLUMED assumes
you want to perform analysis on the entire trajectory.

If you use the above commands to output data from one of the \ref landmarks algorithms then only the second will give you information on the 
atomic positions in your landmark configurations and their associated weights.  The first of these commands will give the values of the colvars
in the landmark configurations only.  If you use the above commands to output data from one of the \ref dimred algorithms then 
\ref OUTPUT_ANALYSIS_DATA_TO_COLVAR will give you an output file that contains the projection for each of your input points.  \ref OUTPUT_ANALYSIS_DATA_TO_PDB
will give you a PDB that contains the position of the input point, the projections and the weight of the configuration.

A nice feature of plumed is that when you use \ref landmarks algorithms or \ref dimred algorithms the output information is just a vector of 
variables.  As such you can use \ref HISTOGRAM to construct a histogram of the information generated by these algorithms.

\page Function Functions

When performing biased dynamics or analyzing a trajectory you may wish to analyze/bias the value of
some function of a set of collective variables rather than the values of the collective variables
directly.  You can do this with PLUMED by using any one of the following list of functions.

Notice that in many functions you should explicitly say to PLUMED whether the result
is a periodic variable or not using the keyword `PERIODIC`.
This is crucial to allow a variable to be properly based.
To know if a function is periodic
of not you should answer to the following question:

- Can my function change with a discontinuity when I move my atoms in a continuous manner?

In case the answer is no, than you should use `PERIODIC=NO`. In case the answer is yes, then you should
consider the following question:

- Are the values of the function at the discontinuity always the same or do they change?

In case the answer is that they are the same, you should use `PERIODIC=A,B` where `A`
is the smallest value and `B` is the largest value. In case the answer is that the
values at the discontinuity are not always the same, then you cannot construct a variable that
can be biased with PLUMED. Consider the following examples:

\plumedfile
t: TORSION ATOMS=1,2,3,4
# When atoms are moved, t could jump suddenly from -pi to +pi

c: MATHEVAL ARG=t FUNC=x*x*x PERIODIC=-31.0062766802998,31.0062766802998
# When atoms are moved, c could jump suddenly from -pi**3 to +pi**3

# equivalently, we could have used:
# c: COMBINE ARG=t POWERS=3 PERIODIC=-31.0062766802998,31.0062766802998

# compute x/y/z components of the distance between atoms 1 and 10
d: DISTANCE ATOMS=1,10 COMPONENTS

# make a new variable equal to d.z but with the correct periodicity
dz: COMBINE ARG=d.z PERIODIC=-10,10
# here we assumed the system is in a orthorhombic box with z side = 20
\endplumedfile

@FUNCTION@

\page Installation Installation

In this page you can learn how to
\link ConfiguringPlumed configure\endlink,
\link CompilingPlumed compile\endlink,
and \link InstallingPlumed install\endlink
PLUMED.
For those of you who are impatient, the following might do the job:
\verbatim
> ./configure --prefix=/usr/local
> make -j 4
> make doc # this is optional and requires proper doxygen version installed
> make install
\endverbatim
Notice that `make install` is not strictly necessary  as plumed can be used from the compilation
directory. This is very useful so as to quickly test the implementation of new features.
However, we strongly recommend to perform a full install.

Once the above is completed the `plumed` executable should be in your execution path
and you will be able to use PLUMED to analyze
existing trajectories or play with the Lennard-Jones code that is included.
However, because PLUMED is mostly used to bias on the fly simulations
performed with serious molecular dynamics packages, 
you can find instructions about how to
\link Patching patch \endlink
your favorite MD code so that it can be combined with PLUMED below.
Again, if you are impatient, something like this will do the job:
\verbatim
> cd /md/root/dir
> plumed patch -p
\endverbatim
Then compile your MD code.
For some MD codes these instructions are insufficient.
It is thus recommended that you read the instructions
at the end of this page.
Notice that MD codes could in principle be "PLUMED ready"
in their official distribution. If your favorite MD code is available "PLUMED ready" 
you will have to compile PLUMED first, then (optionally) install it, then check the MD codes' manual to
discover how to link it.

\section SupportedCompilers Supported compilers

As of PLUMED 2.4, we require a compiler that supports C++11. The following compilers
(or later versions) should be sufficient:

- gcc 4.8.1
- clang 3.3
- intel 15

Notice that the `./configure` script verifies that your compiler supports C++11.
Some compilers do not declare full support, but implement anyway a number of C++11 features
sufficient to compile PLUMED (this is the case for instance of intel 15 compiler).
In case you see a warning about C++11 support during `./configure`
please make sure that PLUMED compiles correctly and, if possible, execute the regtests
(using `make regtest`). Notice that we regularly test a number of compilers on travis-ci,
and at least those compilers are guaranteed to be able to compile PLUMED correctly.

\section ConfiguringPlumed Configuring PLUMED

The `./configure` command 
just generates a Makefile.conf file and a sourceme.sh file.
In PLUMED 2.0 these files were prepared and stored in the 
directory configurations/. The new ones generated by ./configure
are similar to the old ones but are not completely compatible.
In particular, some of the -D options have been changed in version 2.2,
and several new variables so as  to specify the installation directories have been added. For this reason,
you now should run `./configure` again.
Anyway, it should be easy to enforce a similar setup with autoconf by passing
the proper arguments on the command line.
If you have problems on your architecture, please
report them to the mailing list.

Useful command line options for ./configure can be found by typing
\verbatim
> ./configure --help
\endverbatim
PLUMED is made up of modules. Some of them are on by default, some others aren't.
Since version 2.3, the activation of modules should be made during configuration using the `--enable-modules`
option (see \ref mymodules).

Notice that some of the methods within PLUMED depend on external
libraries which are looked for by configure. You can typically
avoid looking for a library using the "disable" syntax, e.g.
\verbatim
> ./configure --disable-mpi --disable-gsl
\endverbatim

Notice that when MPI search is enabled (by default) compilers
such as "mpic++" and "mpicxx" are searched for first. On the other hand,
if MPI search is disabled ("./configure --disable-mpi") non-mpi
compilers are searched for. Notice that only a few of the
possible compiler name are searched. Thus, compilers such as
"g++-mp-4.8" should be explicitly requested with the CXX option.

You can better control which compiler is used by setting the
variables CXX and CC. E.g., to use Intel compilers use the following command:
\verbatim
> ./configure CXX=icpc CC=icc
\endverbatim
Notice that we are using icpc in this example, which is not an MPI compiler as a 
result MPI will not be enabled. Also consider that this is different with respect
to what some other configure script does in that variables such as MPICXX are
completely ignored here. In case you work on a machine where CXX is
set to a serial compiler and MPICXX to a MPI compiler, to
compile with MPI you should use
\verbatim
> ./configure CXX="$MPICXX"
\endverbatim

\warning
This procedure could be somehow confusing since many other programs behave in a different way.
The flag `--enable-mpi` is perfectly valid but is not needed here.
Autoconf will check if a code containing MPI calls can be compiled,
and if so it will enable it. `--disable-mpi` could be used if you are
using a compiler that supports MPI but you don't want PLUMED to  be compiled
with MPI support. Thus the correct way to enable MPI is to pass
to ./configure  the name of a C++ compiler that implements MPI using the CXX option.
In this way, MPI library is treated similarly to all the other libraries
that PLUMED tries to link by default.

To tune the compilation options you can use the CXXFLAGS variable:
\verbatim
> ./configure CXXFLAGS=-O3
\endverbatim

If you are implementing new functionality and want to build with debug flags 
in place so as to do some checking you can use
\verbatim
> ./configure --enable-debug
\endverbatim
This will perform some extra check during execution (possibly slowing down PLUMED)
and write full symbol tables in the executable (making the final executable  much larger).

The main goal of the automatic configure is to find the libraries.
When they are stored in unconventional places it is thus sensible to tell autoconf where 
to look! To do this there are some environment variable that can be used to instruct the linker
which directories it should search for libraries inside. These variables are compiler dependent,
but could have been set by the system administrator so that libraries are found
without any extra flag. Our suggested procedure is to 
first try to configure without any additional flags and to then check the log so as to see whether
or not the libraries were properly detected.

If a library is not found during configuration, you can try to use options to modify the
search path.
For example if your gsl libraries is in /opt/local (this is where MacPorts put it)
and configure is not able to find it you can try
\verbatim
> ./configure LDFLAGS=-L/opt/local/lib CPPFLAGS=-I/opt/local/include
\endverbatim
Notice that PLUMED will first try to link a routine from say gsl
without any additional flag, and then in case of failure will retry adding
"-lgsl" to the LIBS options.
If also this does not work, the gsl library will be
disabled and some features will not be available.
This procedure allows you to use libraries
with custom names. So, if
your gsl library is called /opt/local/lib/libmygsl.so you can 
link it with
\verbatim
> ./configure LDFLAGS=-L/opt/local/lib CPPFLAGS=-I/opt/local/include LIBS=-lmygsl
\endverbatim
In this example, the linker will directly try to link `/opt/local/lib/libmygsl.so`.
This rule is true for all the libraries, so that you will always be able to link
a specific version of a library by specifying it using the LIBS variable.

Since version 2.3.2, the search for the library functions passing to the linker a flag with the standard library name (in the gsl example,
it would be `-lgsl`) can be skipped by using the option `--disable-libsearch`.
Notice that in this manner only libraries that are explicitly passed using the `LIBS` option will be linked. For instance
\verbatim
> ./configure --disable-libsearch LIBS=-lgsl
\endverbatim
will make sure that only gsl is linked and, for instance, BLAS and LAPACK libraries are not.
This might be useful when installing PLUMED within package managers such as MacPorts to
make sure that only desired libraries are linked and thus to avoid to introduce spurious
dependencies. The only exception to this rule is `-ldl`, which is anyway a system library on Linux.

\warning On OSX it is common practice to hard code the full path
to libraries in the libraries themselves. This means that, after having linked
a shared library, that specific shared library will be searched in the same
place (we do the same for the `libplumed.dylib` library, which has an install name hard coded).
On the other hand, on Linux it is common practice not to hard code the full path.
This means that if you use the `LDFLAGS` option to specify the path
to the libraries you want to link to PLUMED (e.g. `./configure LDFLAGS="-L/path"`)
these libraries might not be found later.
The visible symptom is that `src/lib/plumed-shared` will not be linked correctly.
Although the file 'src/lib/plumed-shared' is not necessary, being
able to produce it means that it will be possible to link PLUMED dynamically
with MD codes later.
The easiest solution is to hard code the library search path in this way:
\verbatim
> ./configure LDFLAGS="-L/path -Wl,-rpath,/path"
\endverbatim
Notice that as of PLUMED v2.4 it is possible to use the configure option `--enable-rpath`
to automatically hard code the path defined in `LIBRARY_PATH`:
\verbatim
> ./configure LIBRARY_PATH=/path --enable-rpath
\endverbatim
In this way, the search path used at link time (`LIBRARY_PATH`) and the one saved in the `libplumed.so`
library will be consistent by construction.
In a typical environment configured using module framework (http://modules.sourceforge.net),
`LIBRARY_PATH` will be a variable containing the path to all the modules loaded at compilation
time.

PLUMED needs BLAS and LAPACK. These are treated slightly different from
other libraries. The search is done in the usual way (i.e., first
look for them without any link flag, then add "-lblas" and "-llapack", respectively).
As such  if you want to use a specific version of BLAS or LAPACK
you can make them available to configure by using
\verbatim
> ./configure LDFLAGS=-L/path/to/blas/lib LIBS=-lnameoflib
\endverbatim
If the functions of these libraries are not found, the compiler
looks for a version with a final underscore added.
Finally, since BLAS and LAPACK are compulsory in PLUMED,
you can use a internal version of these libraries that comes as part of PLUMED.
If all else fails the internal version of BLAS and LAPACK are the ones that will be
used by PLUMED.
If you wish to disable any search for external libraries
(e.g. because the system libraries have problems) this can be done with
\verbatim
> ./configure --disable-external-blas
\endverbatim
Notice that you can also disable external LAPACK only, that is use internal LAPACK with external BLAS
using
\verbatim
> ./configure --disable-external-lapack
\endverbatim
Since typically it is the BLAS library that can be heavily optimized, this configuration
should not provide significant slowing down and could be used on systems where
native LAPACK libraries have problems.


As a final resort, you can also edit the resulting Makefile.conf file.
Notable variables in this file include:
- DYNAMIC_LIB : these are the libraries needed to compile the PLUMED
library (e.g. -L/path/to/gsl -lgsl etc). Notice that for the
PLUMED shared library to be compiled properly these should be dynamic
libraries. Also notice that PLUMED preferentially requires BLAS and LAPACK library;
see \ref BlasAndLapack for further info. Notice that the variables 
that you supply with `configure LIBS=something` will end up in this
variable. This is a bit misleading but is required to keep the configuration
files compatible with PLUMED 2.0.
- LIBS : these are the libraries needed when patching an MD code; typically only "-ldl" (needed to have functions for dynamic loading).
- CPPFLAGS : add here definition needed to enable specific optional functions;
e.g. use -D__PLUMED_HAS_GSL to enable the gsl library
- SOEXT : this gives the extension for shared libraries in your system, typically
"so" on UNIX, "dylib" on mac; If your system does not support dynamic libraries or, for some other reason, you would like a static executable you can
just set this variable to a blank ("SOEXT=").

\subsection BlasAndLapack BLAS and LAPACK

We tried to keep PLUMED as independent as possible from external libraries and as such those features
that require external libraries are optional. However, to have a properly working version
of plumed PLUMED you need BLAS and LAPACK libraries.  We would strongly recommend you download these libraries and 
install them separately so as to have the most efficient possible implementations of the functions contained within 
them.  However, if you cannot install BLAS and LAPACK, you can use the internal ones.
Since version 2.1, PLUMED uses a configure script to detect libraries. In case system LAPACK or BLAS
are not found on your system, PLUMED will use the internal replacement.

We have had a number of emails (and have struggled ourselves) with ensuring that PLUMED 
can link BLAS and LAPACK.  The following describes some of the pitfalls that you can fall
into and a set of sensible steps by which you can check whether or not you have set up the configuration
correctly.

Notice first of all that the DYNAMIC_LIB variable in the Makefile.conf
should contain the flag necessary to load the BLAS and LAPACK libraries.  Typically this will be
-llapack -lblas, in some case followed by -lgfortran.  Full path specification with -L may be necessary 
and on some machines the BLAS and LAPACK libraries may not be called -llapack and -lblas.
Everything will depend on your system configuration.

Some simple to fix further problems include:
- If the linker complains and suggests recompiling LAPACK with -fPIC, it means that you have static LAPACK libraries. Either install dynamic LAPACK libraries
or switch to static compilation of PLUMED by stopping to set the SOEXT variable
in the configuration file.
- If the linker complains about other missing functions (typically starting with
  "for_" prefix) then you should also link some Fortran libraries. PLUMED
  is written in C++ and often C++ linkers do not include Fortran libraries by default.
  These libraries are required for LAPACK and BLAS to work. Please check the documentation of your compiler.
- If the linker complains that dsyevr_ cannot be found, try adding
  -DF77_NO_UNDERSCORE to CPPFLAGS
  Notice that "./configure" should automatically try this solution.

\subsection installation-vmdplugins VMD trajectory plugins

PLUMED source code already includes a few selected VMD molfile plugins so as to read a small number
of additional trajectory formats (e.g., dcd, gromacs files, pdb, and amber files).
If you configure PLUMED with the full set of VMD plugins you will be able to read
many more trajectory formats, basically all of those supported by VMD.
To this aim,
you need to download the SOURCE of VMD, which contains
a plugins directory. Adapt build.sh and compile it. At
the end, you should get the molfile plugins compiled as a static
library `libmolfile_plugin.a`. Locate said file and `libmolfile_plugin.h`,
they should be in a directory called `/pathtovmdplugins/ARCH/molfile`
(e.g. `/pathtovmdplugins/MACOSXX86_64/molfile`). Also locate file `molfile_plugin.h`,
which should be in `/pathtovmdplugins/include`.
Then customize the configure command with something along the lines of:

\verbatim
> ./configure LDFLAGS="-L/pathtovmdplugins/ARCH/molfile" CPPFLAGS="-I/pathtovmdplugins/include -I/pathtovmdplugins/ARCH/molfile"
\endverbatim

Notice that it might be necessary to add to `LDFLAGS` the path to your TCL interpreter, e.g.

\verbatim
> ./configure LDFLAGS="-ltcl8.5 -L/mypathtotcl -L/pathtovmdplugins/ARCH/molfile" \
            CPPFLAGS="-I/pathtovmdplugins/include -I/pathtovmdplugins/ARCH/molfile"
\endverbatim

Then, rebuild plumed.

\subsection additional-modules Additional Modules

PLUMED includes some additional modules that by default are not compiled, but can be enabled during configuration.
You can use the option `--enable-modules` to activate some of them, e.g.

\verbatim
> ./configure --enable-modules=module1name+module2name
\endverbatim

For more information on modules see \ref mymodules.

\section CompilingPlumed Compiling PLUMED

Once configured, PLUMED can be compiled using the following command:
\verbatim
> make -j 4
\endverbatim
This will compile the entire code and produce a number of files
in the 'src/lib' directory, including the executable 
'src/lib/plumed'. When shared libraries are enabled,
a shared libraries called 'src/lib/libKernel.so' should also be present.
Notice that the extension could be '.dylib' on a Mac.

In case you want to run PLUMED *without installing it* (i.e. from the compilation
directory), you can use the file 'sourceme.sh' that has been created by
the configure script in the main PLUMED directory.
This file can be "sourced" (presently only working for bash shell)
if you want to use PLUMED *before installing it* (i.e. from the compilation
directory). It is a good idea to source it now, so that you can play with the just compiled
PLUMED:
\verbatim
> source sourceme.sh
\endverbatim

Now
a "plumed" executable should be in your path. Try to type
\verbatim
> plumed -h
\endverbatim

\warning If you are cross compiling, the plumed executable
will not work. As a consequence, you won't be able to run regtests
or compile the manual. This is not a problem.

You can also check if PLUMED is correctly compiled by performing our regression tests.
Be warned that some of them fail because of the different numerical accuracy on different machines.
As of version 2.4, in order to test the `plumed` executable that you just compiled
(prior to installing it) you can use the following command
\verbatim
> make check
\endverbatim
On the other hand, in order to test the `plumed` executable that you just installed (see \ref InstallingPlumed)
you should type
\verbatim
> make installcheck
\endverbatim
In addition, similarly to previous versions of PLUMED, you can test the `plumed` executable
that is in your current path with
\verbatim
> cd regtest
> make
\endverbatim 
You can check the exact version they will use by using the command
\verbatim
> which plumed
\endverbatim
Thus, you can easily run the test suite using a different version of PLUMED
(maybe an earlier version that you already installed), just making sure that it can be 
found in the path. Clearly, if you test a given
version of PLUMED with a test suite from a different version you can expect two
possible kinds of innocuous errors:
- If `plumed` executable is older than the test suite, the tests might fail since they rely on
  some feature introduced in PLUMED in a newer version.
- If `plumed` executable is newer than the test suite, the tests might fail since some
  non-backward compatible change was made in PLUMED. We try to keep the number
  of non-backward compatible changes small, but as you can see in the \ref ChangeLog there
  are typically a few of them at every new major release.

\attention
Even though we regularly perform tests on [Travis-CI](http://travis-ci.org/plumed/plumed2),
it is possible that aggressive optimization or even architecture dependent features
trigger bugs that did not show up on travis. So please always perform the regtests when you install
PLUMED.

Notice that the compiled executable, which now sits in 'src/lib/plumed', relies
on other resource files present in the compilation directory.
This directory should thus stay in the correct place. One should thus not
rename or delete it. In fact the path to the PLUMED root directory is 
hard coded in the plumed executable as can be verified using 
\verbatim
> plumed info --root
\endverbatim
In case you try to use the plumed executable without the compilation
directory in place (e.g. you move away the src/lib/plumed static executable
and delete or rename the compilation directory) PLUMED will 
not work correctly and will give you an
error message
\verbatim
> plumed help
ERROR: I cannot find /xxx/yyy/patches directory
\endverbatim
You can force plumed to run anyway by using the option --standalone-executable:
\verbatim
> plumed --standalone-executable help
\endverbatim
Many features will not be available if you run in this way. However, 
this is currently the only way to use the PLUMED static executable on Windows.

\section InstallingPlumed Installing PLUMED

It is strongly suggested to
install PLUMED in a predefined location.
This is done using
\verbatim
> make install
\endverbatim
This will allow you to remove the original compilation directory,
or to recompile a different PLUMED version in the same place.

To install PLUMED one should first decide the location:
\verbatim
> ./configure --prefix=$HOME/opt
> make
> make install
\endverbatim
As of PLUMED 2.5 you cannot anymore change the location during install.
If you didn't specify the `--prefix` option during configure PLUMED will be installed in /usr/local.
The install command should be executed with root permissions (e.g. "sudo make install")
if you want to install PLUMED on a system directory.

Notice that upon installation PLUMED might need to relink a library.
This was always true until version 2.1, but in version 2.2 libraries should
only be relinked if one changes the install prefix during when typing `make install`.
If root user does not have access to compilers, "sudo -E make install" might solve
the issue.

Upon install, the executable is copied to $prefix/bin, libraries to $prefix/lib,
include files to $prefix/include, and
documentation to $prefix/shared/doc/plumed. Additionally, a directory
$prefix/lib/plumed is created containing several other files, including
patch files, object files (for static patches), etc.
Notice also that these path can be further customized using standard autoconf
directories (e.g. `./configure --bindir=/usr/bin64`).

One should then set the environment properly. We suggest to do it using
the module framework (http://modules.sourceforge.net). An ad hoc generated
module file for PLUMED can be found in $prefix/lib/plumed/src/lib/modulefile
Just edit it as you wish and put it in your modulefile directory.
This will also allow you to install multiple PLUMED versions on your machine and to
switch among them. If you do not want to use modules, you can 
still have a look at the modulefile we did so as to know which
environment variables should be set for PLUMED to work correctly.

If the environment is properly configured one should be able to do
the following things:
- use the "plumed" executable from the command line. This is also possible before installing.
- link against the PLUMED library using the "-lplumed" flag for the linker. This allows
  one to use PLUMED library in general purpose programs
- use PLUMED internal functionality (C++ classes) including
  header files such as "#include <plumed/tools/Vector.h>". This is useful as it may be expedient to
  exploit the PLUMED library in general purpose programs

As a final note, if you want to install several PLUMED versions without using modules then you
should provide a different suffix and/or prefix at configure time:
\verbatim
> ./configure prefix=$HOME/opt --program-suffix=_2.2 --program-prefix=mpi-
> make install
\endverbatim
This will install a plumed executable named "mpi-plumed_2.2". All the other files will be renamed similarly,
e.g. the PLUMED library will be loaded with "-lmpi-plumed_2.2" and the PLUMED header files
will be included with "#include <mpi-plumed_2.2/tools/Vector.h>".
Notice that you can also use arbitrary scripts to edit the name of the executable
with the option --program-transform-name=PROGRAM
(see <a href="http://www.gnu.org/software/autoconf/manual/autoconf-2.69/html_node/Transformation-Examples.html#Transformation-Examples"> autoconf documentation </a> for more info).
These options are useful if you
do not want to set up modules, but we believe that using modules as described above is more flexible.

\section Patching Patching your MD code

A growing number of MD codes can use PLUMED without any modification.
If you are using one of these codes, refer to its manual to know how to activate PLUMED.
In case your MD code is not supporting PLUMED already, you should modify it.
We provide scripts to adjust some of the most popular MD codes
so as to provide PLUMED support.
At the present times we support patching the following list of codes:

@CODESL@

In the section \subpage CodeSpecificNotes you can find information specific for each MD code.

To patch your MD code, you should have already installed PLUMED properly.
This is necessary as you need to have the command "plumed" in your execution
path.  As described above this executable will be in your paths if plumed was 
installed or if you have run sourceme.sh

Once you have a compiled and working version of plumed, follow these steps to add it to
an MD code
- Configure and compile your MD engine (look for the instructions in its documentation).
- Test if the MD code is working properly.
- Go to the root directory for the source code of the MD engine.
- Patch with PLUMED using:
\verbatim
> plumed patch -p
\endverbatim
The script will interactively ask which MD engine you are patching.
- Once you have patched recompile the MD code (if dependencies are set up properly in the MD engine,
  only modified files will be recompiled)

There are different options available when patching. You can check all of them using 
\verbatim
> plumed patch --help
\endverbatim
Particularly interesting options include:
- --static just link PLUMED as a collection of object files. This is only suggested if for external reasons you
  absolutely need a static executable. Notice that with this setting it is often more complicated to configure
  properly the MD code, since all the libraries that PLUMED depends on should be properly specified. The `./configure` script
  does its best in this sense, but sometime it cannot solve the problem. Additionally, this patching mode has been reported
  not to work properly on OSX.
- --shared (default) allows you to link PLUMED as a shared library. As a result when PLUMED is updated, there will be no need to recompile the MD code.
  This is way better than --static since the libraries that PLUMED depends on should be automatically linked.
  Notice that if you later remove the directory where PLUMED is installed also the MD code will not run anymore.
- --runtime allows you to choose the location of the PLUMED library at runtime by setting the variable PLUMED_KERNEL.
  This is probably the most flexible option, and we encourage system administrators to use this option when installing
  PLUMED on shared facilities. Indeed, using this setting it will be possible to update separately the PLUMED library
  and the MD code, leaving to the user the possibility to combine different versions at will. 
  We also recommend to use the provided modulefile (see above) to properly set the runtime environment.

Notice that with PLUMED version <2.5 there was no possibility to link PLUMED as a static library (something like 'libplumed.a').
However, starting with PLUMED 2.5, the `./configure` script will try to set up the system so that a `libplumed.a` file is produced.
Patching an MD code with `--static` with try to link against this static library.
Creation of the `libplumed.a` library can be avoided with `./configure --disable-static-archive`.

If your MD code is not supported, you may want to implement an interface for
it. Refer to the <a href="../../developer-doc/html/index.html"> developer
manual </a>.

\section CrossCompiling Cross compiling

If you are compiling an executable from a different machine, then
`plumed` executable will not be available in the compilation environment.
This means that you won't be able to perform regtests on the machine
nor to compile the manual.
You can try to run the regtests on the computing nodes, but this might require some tweak
since often machines where people do cross compiling have architectures with limited capabilities
on the compute nodes. Also notice that many of the `plumed` options (e.g. patch) are implemented
as shell scripts launched from within the `plumed` executable. If the compute nodes have some limitation
(e.g. they do not allow to fork new processes) these options will not work. Anyway, the PLUMED library
in combination with an MD software should work if both PLUMED and the MD software have been properly compiled.

Also notice that it will not be possible
to use the command `plumed patch` on the machine where you are compiling.
You should thus use `plumed-patch` instead of `plumed patch` (notice that it should be written as a single word).

Try e.g.:
\verbatim
> plumed-patch --help
\endverbatim
This script provides a "shell only" implementation of `plumed patch` that will skip the launch of the `plumed` executable.

Notice that other command line tools will be available in the directory `prefix/lib/progname/`. If configuring with
default values this would be `/usr/local/lib/plumed/plumed-*`. These files are not included in the execution path (prefix/bin)
to avoid clashes, but can be executed also when plumed is cross compiled and the main plumed executable cannot be
launched.

\section Installation-macports Installing PLUMED with MacPorts

If you are using a Mac, notice that you can take advantage of a MacPorts package.
Installing a working plumed should be as easy as:
- Install [MacPorts](https://www.macports.org/)
- Type `sudo port install plumed`

Notice that plumed comes with many variants that can be inspected with the command

    > sudo port info plumed

Plumed uses variants to support different compilers.
For instance, you can install plumed with mpich using

    > sudo port install plumed +mpich

Using more recent clang instead of native compilers is recommended so as to
take advantage of openMP

    > sudo port install plumed +mpich +clang50

Notice that support for c++11 with gcc compilers is someway problematic within MacPorts
due to impossibility to use the system c++ library. For this reason, only clang compilers are supported
(see also [this discussion](https://github.com/macports/macports-ports/pull/1252)).

Variants can be also used to compile with debug flags (`+debug`), to pick a linear algebra library
(e.g. `+openblas`) and to enable all optional modules (`+allmodules`).
Notice that the default variant installed with `sudo port install plumed` is shipped as a compiled
binary, which is significantly faster to install.

In addition, we provide a developer version (typically: a later version not yet considered as stable)
under the subport `plumed-devel` that can be installed with

    > sudo port install plumed-devel

`plumed-devel` also supports the same variants as `plumed` in order to customize the compilation.
`plumed-devel` and `plumed` cannot be installed at the same time.

It is also possible to install a plumed-patched version of gromacs.
For instance, you can use the following command to install
gromacs patched with plumed with clang-5.0 compiler and mpich:

    > sudo port install plumed +mpich +clang50
    > sudo port install gromacs-plumed +mpich +clang50

In case you want to combine gromacs with the unstable version of plumed, use this instead:

    > sudo port install plumed-devel +mpich +clang50
    > sudo port install gromacs-plumed +mpich +clang50

Notice that gromacs should be compiled using the same compiler
variant as plumed (in this example `+mpich +clang50`). In case this is not
true, compilation will fail.

Also notice that gromacs is patched with plumed in runtime mode
but that the path of libplumedKernel.dylib in the MacPorts tree
is hard coded. As a consequence:

- If gromacs is run with `PLUMED_KERNEL` environment variable unset (or set to empty),
  then the MacPorts plumed is used.

- If gromacs is run with `PLUMED_KERNEL` environment variable pointing to another instance
  of the plumed library, the other instance is used.

This is especially useful if you are developing PLUMED since you will be able to install
gromacs once for all and combine it with your working version of PLUMED.

\section Installation-conda Installing PLUMED with conda

If you use the conda package manager you can install a pre-compiled PLUMED binary using the following command:
\verbatim
> conda install -c conda-forge plumed
\endverbatim
Similarly, the python wrappers can be installed with
\verbatim
> conda install -c conda-forge py-plumed
\endverbatim

These packages are part of [conda-forge](https://anaconda.org/conda-forge) and as such should be binary compatible
with other codes from the same distribution. Notice that it should also be possible to combine the installed
plumed kernel with an MD code compiled outside of conda (or within a different conda environment)
if plumed is linked in runtime mode.
The only variable that you need to set in order to access to the installed plumed kernel is
`PLUMED_KERNEL` (e.g., `export PLUMED_KERNEL=/conda/prefix/lib/libplumedKernel.so`).

Notice that binaries are only available for Linux and MacOS and that they have a limited number of features.
In particular, they do not support MPI and do not include optional modules.
However, they can be used to quickly install a working PLUMED version without the need to have a compiler.

Notice that there are additional conda packages on the [plumed](https://anaconda.org/plumed/plumed) channel.
Those packages are for testing only.

\section installingonacluster Installing PLUMED on a cluster

If you are installing PLUMED on a cluster and you want several users to take advantage of it
consider the following suggestions.

First of all, we highly recommend using the module file that PLUMED provides to set up the environment.
Just edit it as necessary to make it suitable for your environment.

Notice that PLUMED can take advantage of many additional features if specific libraries are available upon
compiling it.

Try to patch all MD codes with the `--runtime` option. This will allow independent update of PLUMED and MD codes.
  Users will be able to combine any of the installed gromacs/amber/etc versions with any of the installed PLUMED versions.
Notice that it is sometime claimed that statically linked codes are faster. In our experience, this is not true.
In case you absolutely need a static executable, be ready to face non trivial linking issues. PLUMED is written in C++,
thus required the appropriate C++ library to be linked, and might require additional libraries (e.g. libgsl).

Sometime we make small fixes on the patches. For this reason, keep track of which version of PLUMED you used
to patch each of the MD code. Perhaps you can call the MD code modules with names such as `gromacs/4.6.7p1`,
`gromacs/4.6.7p2` and write somewhere in the module file which version of PLUMED you used. Alternatively, call them
something like `gromacs/4.6.7p2.2.0`. In this way, when we report a bug on the mailing list, users will know if the version
they are using is affected by it.

Usually it is not necessary to install both a MPI and a non-MPI PLUMED version. PLUMED library only calls MPI functions
when the MD code is compiled with MPI. PLUMED executable calls MPI functions only when it is invoked without `--no-mpi`.
In many machines it is thus sufficient to run the plumed executable on the login node as
\verbatim
> plumed --no-mpi
\endverbatim
even though PLUMED was compiled with MPI and the login node does not support MPI.
The only case where you might need two different PLUMED installation for compute
and login node is when you are cross compiling.

PLUMED needs to be well optimized to run efficiently.
If you need a single PLUMED binary to run efficiency on machines with different levels of hardware (e.g.: some
of your workstations support AVX and some do not), with intel compiler you can use something like
\verbatim
> ./configure CXX=mpicxx CXXFLAGS="-O3 -axSSE2,AVX"
\endverbatim
It will take more time to compile but it will allow you to use a single module. Otherwise, you should install two
PLUMED version with different optimization levels.

Using modules, it is not necessary to make the PLUMED module explicitly dependent on the used library. Imagine a
scenario where you first installed a module `libgsl`, then load it while you compile PLUMED. If you
provide the following option to configure `--enable-rpath`, the PLUMED executable and
library will remember where libgsl is, without the need to load libgsl module at runtime.
Notice that this trick often does not work for fundamental libraries such as C++ and MPI library. As a consequence,
usually the PLUMED module should load the compiler and MPI modules.

\attention
In case you found out how to compile PLUMED on some fancy architecture please share your tricks! You can
either post it in your blog, send it to the mailing list, or ask as to update this paragraph in the manual, we will
be happy to do so.

\section installingpython Installing Python wrappers

As of PLUMED 2.5 it is possible to use the PLUMED library through Python wrappers. Notice that this is not something for end users but rather for developers. The interface is very similar to the one used in MD codes linked with PLUMED.

There are two ways to install Python wrappers.

\subsection installingpython-inside Installing Python wrappers within PLUMED

If `./configure` finds a `python` executable that also has the `cython` module available, Python wrappers will be installed within `/prefix/lib/plumed/python`. In order to access them, you should add this directory to the environment variable `PYTHONPATH`. Notice that if your python interpreter has a different name you might have to pass it to `./configure` with `PYTHON_BIN=python3.6`. The whole thing would then be:

````
./configure PYTHON_BIN=python3.6 --prefix=$HOME/opt
make && make install
export PYTHONPATH="$HOME/opt/lib/plumed/python:$PYTHONPATH"
python3.6
>>> import plumed
````

Notice that in this manner you will have to commit to a specific python version **before** installing PLUMED.

\subsection installingpython-outside Installing Python wrappers outside PLUMED

If you use multiple python versions, you might find it easier to install the Python wrappers separately from PLUMED. The simplest way is to do it with `pip`:

````
pip3.6 install --user plumed
````

Here the `--user` flag allows you to install the packages on your home. Notice that you don't even need to download PLUMED in order to install the wrappers, but you will need PLUMED in order to use them. You can tell the wrappers where PLUMED is by setting the `PLUMED_KERNEL` environment variable:

````
export PLUMED_KERNEL=$HOME/opt/lib/libplumedKernel.so
python3.6
>>> import plumed
````

Notice that by installing the wrappers in this manner you will download those that are packaged on [Pypi](https://pypi.org/project/plumed/).
If you want to install using pip the development version of the wrappers you should download the PLUMED repository and use
the following commands:

````
pip3.6 install --user cython # cython is required in this case
cd plumed2/python
make pip
pip3.6 install --user .
````

If you want to install the development version it is recommended to use a virtualenv so that it will not interfere with the released packages.

\section installinghints Other hints

We here collect a list of suggestions that might be useful on particular
machines.

- On Blue Gene Q (likely on AIX) the prelinking made with `ld -r` is not
  working properly. There is no easy way to detect this at configure time.
  If during `make` you receive an error in the form
\verbatim
ld: TOC section size exceeds 64k
\endverbatim
  please configure plumed again with the following flag
\verbatim
> ./configure --disable-ld-r
\endverbatim
- On Cray machines, you might have to set the following environment variable
  before configuring and building both PLUMED and the MD code that you want
  to patch with PLUMED (kindly reported by Marco De La Pierre):
\verbatim
> export CRAYPE_LINK_TYPE=dynamic
\endverbatim
- Intel MPI seems to require the flags `-lmpi_mt -mt_mpi` for compiling and linking and the flag `-DMPICH_IGNORE_CXX_SEEK` for compiling
  (kindly reported by Abhishek Acharya).
  You might want to try to configure using
\verbatim
> ./configure LDFLAGS=-lmpi_mt CXXFLAGS="-DMPICH_IGNORE_CXX_SEEK -mt_mpi" STATIC_LIBS=-mt_mpi
\endverbatim
  Adding libraries to `STATIC_LIBS` uses them for all the linking steps, whereas those in `LIBS` are only used when linking the PLUMED kernel library.
  See more at [this thread](https://groups.google.com/d/msgid/plumed-users/CAB1aw3y0m%3D5qwzsZY4ZB-aBevsL5iuS%3DmQuSWK_cw527zCMqzg%40mail.gmail.com?utm_medium=email&utm_source=footer).

\page CodeSpecificNotes Code specific notes

Here you can find instructions that are specific for patching each of the supported MD codes.
Notice that MD codes with native PLUMED support are not listed here.

@CODES@

\page Files Files

We tried to design PLUMED in such a manner that input/output is done consistently
irrespective of the file type. Most of the files written or read by PLUMED thus follow
the very same conventions discussed below. 

\section Restart

Whenever the \ref RESTART option is used, all the files written by PLUMED are appended.
This makes it easy to analyze results of simulations performed as a chain of several sub-runs.
Notice that most of the PLUMED textual files have a header. The header is repeated at every
restart. Additionally, several files have time in the first column. PLUMED just takes the value
of the physical time from the MD engine. As such, you could have that time starts again from zero
upon restart or not.

An exception from this behavior is given by files which are not growing as the simulation proceeds.
For example, grids written with \ref METAD with GRID_WFILE are overwritten by default during the simulation.
As such, when restarting, there is no point in appending the file. Internally, PLUMED opens the file in append
mode but then rewinds it every time a new grid is dumped.

\section Backup

Whenever the \ref RESTART option is not used, PLUMED tries to write new files. If an old file
is found in the way, PLUMED takes a backup named "bck.X.filename" where X is a progressive number.
Notice that by default PLUMED only allows a maximum of 100 backup copies for a file.
This behavior can be changed by setting the environment variable PLUMED_MAXBACKUP to the desired number
of copies. E.g. export PLUMED_MAXBACKUP=10 will fail after 10 copies. PLUMED_MAXBACKUP=-1 will never fail - be careful
since your disk might fill up quickly with this setting.

\section Replica-Suffix Replica suffix

When running with multiple replicas (e.g., with GROMACS, -multi option) PLUMED adds the replica index as a suffix to
all the files. The following command
will thus print files named COLVAR.0, COLVAR.1, etc for the different replicas.
\plumedfile
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=COLVAR
\endplumedfile

When reading a file, PLUMED will try to add the suffix. If the file is not found, it will fall back to
the name without suffix. The most important case is the reading of the plumed input file.
If you provide a file for each replica (e.g. plumed.0.dat, plumed.1.dat, etc) you will be able to
setup plumed differently on each replica. 
On the other hand, using a single plumed.dat will make all the replicas read the same file.

\warning This rule is true for almost all the files read by PLUMED. As of
  PLUMED version 2.4, the only exception is PDB files, where the replica suffix is not added.

Notice that when PLUMED adds the replica suffix, it recognizes the file extension and add the suffix _before_ the
extension. Before PLUMED 2.2, the only recognized suffix was ".gz". Since 2.2, any suffix with length
less or equal to five letters is recognized.

This means that using in a multi-replica context an input such as
\plumedfile
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=COLVAR.gz
METAD ARG=d FILE=test.HILLS SIGMA=0.1 HEIGHT=0.1 PACE=100
\endplumedfile
PLUMED will write files named COLVAR.0.gz, COLVAR.1.gz, test.0.HILLS, test.1.HILLS, etc
etc. This is useful since the preserved extension makes it easy
to process the files later.

\page Bias Bias

PLUMED allows you to run a number of enhanced sampling algorithms.
The list of enhanced sampling algorithms contained in PLUMED is as follows: 

@BIAS@

Methods, such as \ref METAD or \ref PBMETAD, that work by introducing a history dependent bias can be restarted 
using the \subpage RESTART keyword

\page Group Groups and Virtual Atoms 

\section atomSpecs Specifying Atoms

The vast majority of the CVs implemented in PLUMED are calculated from a list of atom positions.  Within PLUMED
atoms are specified using their numerical indices in the molecular dynamics input file. 

In PLUMED lists of atoms can be either provided directly inside the definition of each collective variable, or
predefined as a \subpage GROUP that can be reused multiple times. Lists of atoms can be written as:

- comma separated lists of numbers (`GROUP ATOMS=10,11,15,20 LABEL=g1`)
- numerical ranges.  So `GROUP ATOMS=10-20 LABEL=g2` is equivalent to `GROUP ATOMS=10,11,12,13,14,15,16,17,18,19,20 LABEL=g2`
- numerical ranges with a stride. So `GROUP ATOMS=10-100:10 LABEL=g3 is equivalent to `GROUP ATOMS=10,20,30,40,50,60,70,80,90,100 LABEL=g3`
- atoms ranges with a negative stride. So `GROUP ATOMS=100-10:-10 LABEL=g4 is equivalent to `GROUP ATOMS=100,90,80,70,60,50,40,30,20,10 LABEL=g4`

In addition, there are a few shortcuts that can be used:

- `@mdatoms` indicate all the physical atoms present in the MD engine (e.g. `DUMPATOMS ATOMS=@mdatoms`).
- `@allatoms` indicates all atoms, including \ref vatoms "those defined only in PLUMED" (e.g. `DUMPATOMS ATOMS=@allatoms`).

The list of the virtual atoms defined in PLUMED can be obtained by using the command `GROUP ATOMS=@allatoms REMOVE=@mdatoms`.

Other shortcuts are available if you loaded the structure of the molecule using the \ref MOLINFO command.

All the above methods can be combined just putting one name after the other separated by a comma:
\plumedfile
DUMPATOMS ATOMS=1,2,10-20,40-60:5,100-70:-2 LABEL=g5 FILE=test.xyz
\endplumedfile

Some collective variable must accept a fixed number of atoms, for example a \ref DISTANCE is calculated
using two atoms only, an \ref ANGLE is calculated using either 3 or 4 atoms and \ref TORSION is calculated using 4 atoms.

Additional material and examples can be also found in the tutorial \ref belfast-1. 

\subsection mols Molecules

In addition, for certain colvars, pdb files can be read in using the following keywords and used to select ATOMS:

@TOPOLOGY@

\subsection pbc Broken Molecules and PBC 

PLUMED is designed so that for the majority of the CVs implemented the periodic boundary conditions are treated 
in the same manner as they would be treated in the host code.  In some codes this can be problematic when the colvars
you are using involve some property of a molecule.  These codes allow the atoms in the molecules to become separated by 
periodic boundaries, a fact which PLUMED could only deal with were the topology passed from the MD code to PLUMED.  Making this
work would involve a lot laborious programming and goes against our original aim of having a general patch that can be implemented 
in a wide variety of MD codes.  Consequentially, we have implemented a more pragmatic solution to this problem - the user specifies
in input any molecules (or parts of molecules) that must be kept in tact throughout the simulation run.  In PLUMED 1 this was done
using the ALIGN_ATOMS keyword.  In PLUMED 2 the same effect can be achieved using the \subpage WHOLEMOLECULES command.

The following input computes the end-to-end distance for a polymer of 100 atoms and keeps it at a value around 5.

\plumedfile
WHOLEMOLECULES ENTITY0=1-100
e2e: DISTANCE ATOMS=1,100 NOPBC
RESTRAINT ARG=e2e KAPPA=1 AT=5
\endplumedfile

Notice that NOPBC is used to be sure in \ref DISTANCE that if the end-to-end distance is larger than half the simulation box the distance 
is compute properly. Also notice that, since many MD codes break molecules across cell boundary, it might be necessary to use the 
\ref WHOLEMOLECULES keyword (also notice that it should be before distance).

Notice that most expressions are invariant with respect to a change in the order of the atoms,
but some of them depend on that order. E.g., with \ref WHOLEMOLECULES it could be useful to
specify atom lists in a reversed order.

\plumedfile
# to see the effect, one could dump the atoms as they were before molecule reconstruction:
# DUMPATOMS FILE=dump-broken.xyz ATOMS=1-20
WHOLEMOLECULES STRIDE=1 ENTITY0=1-20
DUMPATOMS FILE=dump.xyz ATOMS=1-20
\endplumedfile

Notice that there are other ways to manipulate the coordinates stored within PLUMED:
- Using the \subpage FIT_TO_TEMPLATE they can be aligned to a template structure.
- Using \subpage WRAPAROUND you can bring a set of atom as close as possible to another set of
  atoms.
- Using \subpage RESET_CELL you can rotate the periodic cell.

\section vatoms Virtual Atoms

Sometimes, when calculating a colvar, you may not want to use the positions of a number of atoms directly.  Instead
 you may wish to use the position of a virtual atom whose position is generated based on the positions of a collection 
of other atoms.  For example you might want to use the center of mass of a group of atoms.  Plumed has a number of routines
for calculating the positions of these virtual atoms from lists of atoms:

@VATOM@

To specify to a colvar that you want to use the position of a virtual atom to calculate a colvar rather than one of the atoms
in your system you simply use the label for your virtual atom in place of the usual numerical index.  Virtual
atoms and normal atoms can be mixed together in the input to colvars as shown below:

\plumedfile
COM ATOMS=1,10 LABEL=com1
DISTANCE ATOMS=11,com1
\endplumedfile

If you don't want to calculate CVs from the virtual atom.  That is to say you just want to monitor the position of a virtual atom 
(or any set of atoms) over the course of your trajectory you can do this using \ref DUMPATOMS.

\page MAZE MAZE

<!-- 
description: Module that implements enhanced sampling methods for ligand unbinding from protein tunnels
authors: Jakub Rydzewski
reference: \cite RydzewskiMaze 
-->

maze is a module for PLUMED2, which implements enhanced sampling methods for 
ligand unbinding from protein tunnels. The maze module is developed and 
maintained by [Jakub Rydzewski](http://www.fizyka.umk.pl/~jr) at the Institute 
of Physics, Nicolaus Copernicus University, Torun, Poland. See this 
[link](https://www.fizyka.umk.pl/~jr/maze.html) for additional information.

The maze module is an optional module for PLUMED2 that needs to be enabled when 
configuring the compilation of PLUMED2. You can either pass a flag
'\-\-enable-modules=maze' or a '\-\-enable-modules=all' when running the 
configure script. 

See the following sections for further information:

- \subpage maze_loss
- \subpage maze_optimizer
- \subpage maze_bias

\page maze_loss Loss

The following list contains the loss functions available in the maze module.

@MAZE_LOSS@

\page maze_optimizer Optimizers

The following list contains the optimizers available in the maze module.

@MAZE_OPTIMIZER@

\page maze_bias Biases

The following list contains the biases available in the maze module.

@MAZE_BIAS@
\page AddMod Additional Modules

Here is collected the documentation for the additional modules contributed to PLUMED.
For information on how to enable modules see \ref mymodules.

<table align=center frame=void cellpadding=5%%>
<tr><td><b>Module</b></td><td><b>Description</b></td><td><b>Authors</b></td><td><b>References</b></td></tr>

@ADDITIONALMODULES@

</table>
\page FISSTMOD FISST (Infinite Switch Simulated Tempering in Force)

<!-- 
description: Infinite Switch Simulated Tempering in Force (FISST)
authors: Glen Hocky
reference: \cite Hartmann-FISST-2019
-->

## Overview

This FISST module contains methods for adaptively determining weight parameters to construct a bias function that represents the Infinite Switch limit of Simulated Tempering for a linear bias coefficient of a CV, as described in \cite Hartmann-FISST-2019.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=fisst' to your './configure' command when building PLUMED to enable these features.

## Usage
Currently, all features of the FISST module are included in a single FISST bias function: \ref FISST

## Module Contents
- \subpage FISSTMODBias

\page FISSTMODBias Biases Documentation

The following list contains descriptions of biases developed for the PLUMED-FISST module. They can be used in combination with other biases outside of the FISST module.

@FISSTMOD_BIAS@
\page PYTORCH PYTORCH (Machine Learning Collective Variables)

<!-- 
description: Machine Learning Collective Variables with PyTorch (pytorch)
authors: Luigi Bonati
reference: \cite bonati2020data
-->

## Overview 

The PYTORCH module is an interface between PyTorch machine learning library and PLUMED. It implements the \ref PYTORCH_MODEL class, which is a subclass of `Function` class. `PYTORCH_MODEL` provide the ability to load models defined in Pytorch and compiled with <a href="https://pytorch.org/docs/stable/jit.html#"> TorchScript</a>. 

For instance, this allows one to use the outputs of a neural network as collective variables, as done in \cite bonati2020data and in \cite bonati2021deep. Furthermore, the \ref PYTORCH_MODEL outputs can also be used as inputs for other collective variables and for data analysis tools. 

## Installation

This module is not installed by default. It requires the PyTorch C++ APIs (LibTorch) to be linked against PLUMED. Since these APIs are still in beta phase regarding stability, it is strongly suggested to use the 1.8.* LTS version of both PyTorch and LibTorch. Please note that if you want to link a different version it might be necessary to specify the required libraries within LIBS in configure. 
Furthermore, note that if the versions of PyTorch and LibTorch do not match it might not be possible to correctly load the model. 

### Download LibTorch C++ API library

You can download the pre-built LibTorch library from their <a href="https://pytorch.org/get-started/locally/"> website</a>. The following script downloads the recommended 1.8.2 LTS version (CPU, with C++11 ABI compatibility).

\verbatim
> wget https://download.pytorch.org/libtorch/lts/1.8/cpu/libtorch-cxx11-abi-shared-with-deps-1.8.2%2Bcpu.zip 
> unzip libtorch-cxx11-abi-shared-with-deps-1.8.2+cpu.zip 
> rm libtorch-cxx11-abi-shared-with-deps-1.8.2+cpu.zip
> LIBTORCH=${PWD}/libtorch
\endverbatim

The location of the include and library files need to be exported in the environment. For convenience, we can save them in a file `sourceme.sh` inside the libtorch folder:

\verbatim
> echo "export CPATH=${LIBTORCH}/include/torch/csrc/api/include/:${LIBTORCH}/include/:${LIBTORCH}/include/torch:$CPATH" >> ${LIBTORCH}/sourceme.sh
> echo "export INCLUDE=${LIBTORCH}/include/torch/csrc/api/include/:${LIBTORCH}/include/:${LIBTORCH}/include/torch:$INCLUDE" >> ${LIBTORCH}/sourceme.sh
> echo "export LIBRARY_PATH=${LIBTORCH}/lib:$LIBRARY_PATH" >> ${LIBTORCH}/sourceme.sh
> echo "export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH" >> ${LIBTORCH}/sourceme.sh
> . ${LIBTORCH}/sourceme.sh
\endverbatim

Remember to source the `sourceme.sh` file in your `~/.bashrc` or  `~/.bash_profile` file. 

### Configure PLUMED

In order to install the `PYTORCH` module when compiling PLUMED we need to (1) specify to look for libtorch (`--enable-libtorch`) and (2) enable the related module (`--enable-modules=pytorch` or also `--enable-modules=all`):

\verbatim
> ./configure --enable-libtorch --enable-modules=pytorch  
\endverbatim

### Notes about the linking of LibTorch

- A compiler with C++14 support is required. 
- If you want to use the pre-cxx11 ABI LibTorch binaries the following flag should be added to the configure: `CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"`.


## Usage

Currently, all features of the PYTORCH module are included in a single function: \ref PYTORCH_MODEL

## Training CVs with PyTorch: the mlcvs package

<a href="https://mlcvs.readthedocs.io/"> `mlcvs` </a> is a Python package for the design of different kinds of neural-networks based CVs, e.g. that discriminate between states \ref \cite bonati2020data or that approximate the slow dynamical modes of the system \cite bonati2021deep. The CVs are optimized in Python and the resulting model is compiled with TorchScript, in order to allowed the models to be employed without Python dependencies.

## Module Contents
- \subpage PYTORCHFunction

\page PYTORCHFunction Functions Documentation

The following list contains descriptions of functions developed for the PYTORCH module. They can be used in combination with other actions outside of the PYTORCH module.

@PYTORCH_FUNCTION@\page PIVMOD PIV collective variable

<!-- 
description: Permutation invariant collective variable (PIV) 
authors: S. Pipolo, F. Pietrucci 
reference: \cite gallet2013structural \cite pipolo2017navigating 
-->

## Overview

To be completed

## Installation 
This module is not installed by default. Add '\-\-enable-modules=piv' to your './configure' command when building PLUMED to enable these features.

## Usage
Currently, all features of the PIV module are included in a single PIV collective variable: \ref PIV

## Module Contents
- \subpage PIVMODColvar

\page PIVMODColvar CVs Documentation

The following list contains descriptions of biases developed for the PLUMED-PIV module. They can be used in combination with other actions outside of the PIV module.

@PIVMOD_COLVAR@
\page glossary Index of Actions 

The following page contains an alphabetically ordered list of all the Actions and command line tools 
that are available in PLUMED 2. For lists of Actions classified in accordance with the 
particular tasks that are being performed see:

- \ref colvarintro tells you about the ways that you can calculate functions of the positions of the atoms.
- \ref Analysis tells you about the various forms of analysis you can run on trajectories using PLUMED.
- \ref Bias tells you about the methods that you can use to bias molecular dynamics simulations with PLUMED.

\section ActionList Full list of actions

@GLOSSARY@

\page OPES OPES (On-the-fly Probability Enhanced Sampling)

<!-- 
description: On-the-fly Probability Enhanced Sampling (OPES)
authors: Michele Invernizzi
reference: \cite Invernizzi2020rethinking \cite Invernizzi2020unified
-->

## Overview

The OPES module contains the implementation of the on-the-fly probability enhanced sampling mehtod (OPES) \cite Invernizzi2020rethinking \cite Invernizzi2022explore \cite Invernizzi2020unified.

The OPES method aims at sampling a given target distribution over the configuration space, \f$p^{\text{tg}}(\mathbf{x})\f$,
different from the equilibrium Boltzmann distribution, \f$P(\mathbf{x})\propto e^{-\beta U(\mathbf{x})}\f$.
To do so, it incrementally builds a bias potential \f$V(\mathbf{x})\f$, by estimating on-the-fly the needed probability distributions:
\f[
V(\mathbf{x}) = -\frac{1}{\beta}\log\frac{p^{\text{tg}}(\mathbf{x})}{P(\mathbf{x})}\, .
\f]
The bias quickly becomes quasi-static and the desired properties, such as the free energy, can be calculated with a simple reweighting \ref REWEIGHT_BIAS.

Depending on the kind of target distribution one wishes to sample, different \ref OPES_BIAS "OPES biases" can be used.

## Installation 
This module is not installed by default. Add '\-\-enable-modules=opes' to your './configure' command when building PLUMED to enable these features. See also \ref mymodules.

## Usage
The OPES module contains three bias actions, \ref OPES_METAD and \ref OPES_METAD_EXPLORE that sample metadynamics-like target distributions (e.g. the well-tempered one), and \ref OPES_EXPANDED that samples expanded ensembles target distributions (replica-exchange-like).
It also contains various expansion collective variables (ECVs) to define such expanded targets.

## Module Contents
- \subpage OPES_bias
- \subpage expansion_CV
- \subpage OPES_tutorial

\page OPES_bias Biases

The following list contains the biases available in the OPES module.

@OPES_BIAS@

\page expansion_CV Expansion Collective Variables

Expansion collective variables (ECVs) are used to define the expanded ensemble to be sampled by \ref OPES_EXPANDED.
See Ref.\cite Invernizzi2020unified for details.
ECVs take as argument some underlying colvar and have as output components the same colvars.

The following list contains the expansion CVs available in the OPES module.

@OPES_EXPANSION_CV@

\page OPES_tutorial Tutorials

The following list contains the tutorials available in the OPES module.

@OPES_TUTORIALS@

\page AddingAMetric Implementing methods for calculating the distances between pairs of configurations 

To implement a new method for calculating the distance between a pair of trajectory frames you will need to work with the
PLMD::reference module.  This module is used in many parts of PLUMED (e.g. path collective variables, field-cvs and analysis methods).
Consequently, if you implement your distance measure using the functionality in this class you will be able to use it in a 
wide variety of difffernt contexts.  As always if you would like us to incorporate your measure in the release version of 
PLUMED you will need to write at least one regression test for it.  Details on how to write regression tests are provided 
here: \ref regtests  

The classes in PLMD::reference allow one to calculate the distance between pairs of trajectory frames.  The reason that this is 
a module rather than a single method is that there are a wide variety of different ways of calculating the distance between
two frames.  One could, for example, calculate the RMSD distance between the two frames.  Alternatively, one can calculate 
a large number of collective variables and compare the values these variables have in the two configurations.  Lastly, one
could somehow combine some element of RMSD calculation with the calculation of some set of collective variables.  As with
so many things the way in which one wishes to calculate the distance between two configurations will depend on the problem
one is endeavoring to solve.  The aim of the PLMD::reference module is thus to provide the user unlimited flexibility in the
way that the matrix of distances can be calculated in analysis methods such as sketch-map or in biasing methods such as path
cvs.  I say unlimited because, although the particular distance measure the user needs many not currently be implemented in 
PLUMED, he/she always has the option to implement this new distance metric in the reference module and then access it in the
full range of analysis and biasing tools that are already available in the code.  The following provides instructions for 
implementing a new way of calculating the distance between a pair of trajectory frames in the code using the PLMD::reference 
module.  As always once your new method is implemented you can access it in all the places where pairwise distances are employed
because of the way that PLUMED exploits inheritance and polymorphism.   

\section adding Creating your measure

To create a new way of measuring the distances between pairs of atoms you must write a new class in the reference directory.
An example declaration for such a class is given below:

\verbatim
class OptimalRMSD : public RMSDBase {
private:
  bool fast;
  RMSD myrmsd;
public:
  OptimalRMSD(const ReferenceConfigurationOptions& ro);
  void read( const PDB& );
  double calc( const std::vector<Vector>& pos, const bool& squared );
};
\endverbatim

In this case we are inheriting from PLMD::RMSDBase but the particular class you will want to inherit from will depend on the
manner in which the distance is calculated.  To then ensure that your new measure can be used throughout the code you need
to include the MetricRegister.h file and the following line:

\verbatim
PLUMED_REGISTER_METRIC(OptimalRMSD,"OPTIMAL")
\endverbatim

Once again dynamic polymorphism is exploited again here.  With this example the command above ensures that PLMD::OptimalRMSD objects
are used whenever the command METRIC=OPTIMAL is found in the PLUMED input.

Your new class must contain a constructor, a method to read in the configuration from a pdb file (read) and a method to calculate
the distance between the reference configuration and the input position (calc).  Please be aware that the number of arguments
to this method will change depending on which base class you inherit from in creating your new measure object.

The inheritance structure of these routines is rather complicated looking.  In essence, however, the base class underlying all
these classes in PLMD::ReferenceConfiguration.  There are then two classes PLMD::ReferenceArguments and PLMD::ReferenceAtoms that
can be multiply inherited in any derived classes.  PLMD::ReferenceArguments provides tools for dealing with distance measures that
involve colvars. PLMD::ReferenceAtoms provides tools for dealing with distance measures that involve the positions of atoms.
Base classes such as PLMD::SingleDomainRMSD and PLMD::RMSDBase are there essentially so that we can use particular sets of metrics
in secondary structure variables and the RMSD class respectively.  Essentially within these two objects the RMSD is calculated by
calling the calc member of the abstract base class PLMD::SingleDomainRMSD and PLMD::RMSDBase.  This allows one to use the
appropriate set of measures within these particular functions.

\section args Dealing with colvars

There are essentially three ways in which you might wish to calculate the distance between two sets of reference colvar values:

- You will calculate the euclidean distance using pythagoras theorem
- You calculate the normalised euclidean distance in which pythagoras theorem is again used but each pair of components is given a separate weight
- You calculate the Mahalonobis distance in which the distance is calculated as \f$x^T M x\f$ where \f$x\f$ is the displacement vector and \f$M\f$ is a matrix.

These three methods are all implemented within PLUMED and are in the classes PLMD::EuclideanDistance, PLMD::NormalizedEuclideanDistance and PLMD::MahalanobisDistance
respectively.  If you look in these three classes you will note that there is a very small amount of code in each of them.  Essentially the calculation of the
distance and the reading in of the PDB file are looked after by the methods PLMD::ReferenceArguments::calculateArgumentDistance and
PLMD::ReferenceArguments::readArgumentsFromPDB respectively.  To reuse these functionalities you need to add the command:

- hasmetric=true in the constructor if you want to use the Mahalonobis distance in your new metric
- hasweights=true in the constructor if you want to use the normalised euclidean distance in your new metric

If you want to use the euclidean distance you can use PLMD::ReferenceArguments::calculateArgumentDistance and PLMD::ReferenceArguments::readArgumentsFromPDB
without any further instructions.  Notice these methods will still work if you use a combination of atom positions and colvars in your measure.  They will
give the part of the measure due to the arguments - you will be required to do some further calculation to get the bit involving the atoms.

\section atoms Dealing with atoms

All the distance measures that we work with involve the positions of the atoms in the two configurations.  If we work with colvars we
just do the calculation of the distance from the positions of the atoms in an indirect way - we calculate some intermediary quantities
from the positions of the atoms and then calculate the set of difference between these intermediary quantities.  There are cases,
however, such as when we are working with RMSD distances where it is useful to calculate the distance directly from the set of atomic
positions.  That is to say there are cases where these intermediary quantities are not useful.  If you are implementing such a measure
you will need to write a class that inherits from PLMD::ReferenceAtoms either directly or indirectly.

Within the read method you can read the atoms in the PDB file by using the method PLMD::ReferenceAtoms::readAtomsFromPDB.  Within calc you can then
access the positions of these read in atoms by using the method PLMD::ReferenceAtoms::getReferencePositions() or by using PLMD::ReferenceAtoms::getReferencePosition.
It is useful to think carefully about how you can reuse other parts of the code when implementing new reference methods.  As an example notice
how PLMD::OptimalRMSD and PLMD::SimpleRMSD make extensive use of the PLMD::RMSD class.  Similarly notice how the PLMD::OptimalRMSD,
PLMD::SimpleRMSD and PLMD::DRMSD classes are reused in PLMD::MultiDomainRMSD.

\section AddingAMeasureDocs Adding documentation for your measure

To test whether you have implemented you new measure correctly you should implement a small Action that calculates the
distance between the instantaneous positions of the atoms and some read in reference configuration.  The methods
PLMD::colvar::RMSD, PLMD::colvar::DRMSD, PLMD::colvar::MultiRMSD and PLMD::function::Target perform these functions for
the set of measures that are currently implemented within the reference module.  You will need to do something similar
in your test action.

You will notice that the documentation in these files starts with the line:

\verbatim
//+PLUMEDOC DCOLVAR TARGET
\endverbatim

The DCOLVAR tag is important here as including this tag ensures that the documentation for these objects appears in the
appropriate part of the manual.  This page of the manual is linked to from all of the pages that exploit the reference
module's functionality to provide multiple methods for calculating the distance between two trajectory frames.  Any measure
that you implement should thus have one of these wrapper actions associated with it and the documentation for the measure
should be included in the wrapper code's source code file.

\page regtests Adding regressions tests

When you write new functionality for PLUMED it is important that you document these features AND that you provide suitable
regression tests so as to ensures that developers working on the code in the future do not break any existing features.
The regression tests that are currently within PLUMED are all in the regtests directory separated.  The tests are contained 
in a number of separate directories that (roughly) correspond to the modules of PLUMED.  These tests are run whenever the 
user executes the

\verbatim
make tests
\endverbatim

command.  Furthermore, they are also executed every time one of the core developers pushes commits to the main plumed2 fork
on github.  In fact these tests are executed across multiple architectures whenever a push is performed.  To do this we use 
a web-app called travis-ci.  The current status of the code on travis-ci can be viewed here: https://travis-ci.org/plumed/plumed2
<b> We take code testing various seriously and will not consider including contributions from developers unless appropriate 
regression tests are included. </b> 

\section addingregtests Creating a regression test

The first step in creating a regression test is to create a directory for it somewhere in the regtest directory of PLUMED.  If you
are implementing your own module for PLUMED then it is best to keep your regression tests in a subdirectory named after your module.
If you are not doing this then you can create a directory in one of the already available module directories.  Your regression test
module <b> must </b> be named:

\verbatim
rt-something
\endverbatim

where the something should be replaced with a suitable description of the feature that is being tested within the module.  The directory 
name must begin with rt-*, however, as this is how the various Makefiles than run the regression tests identify the directories that 
contain tests that must be run.

\subsection regconfig Creating the config file and the Makefile

Once you have created the directory for your regression test you will need to create two files: a Makefile and a config file.  These files
can be copied from one of the regtest directories that are already there.  The Makefile does not need to be changed as it essentially 
just includes a script from elsewhere in PLUMED.  You may have to edit the config file, however, so in order to understand what the various 
lines in this file do here is an example config file: 

\verbatim
mpiprocs=2
type=driver
plumed_modules=adjmat
arg="--plumed plumed.dat --trajectory-stride 50 --timestep 0.005 --ixyz diala_traj_nm.xyz --dump-forces forces --dump-forces-fmt=%10.6f"
extra_files="../../trajectories/diala_traj_nm.xyz  ../trajectories/path_msd/all.pdb"
\endverbatim

- The first line in this file - the one containing the string mpiprocs=2 - tells PLUMED that when this regtest should be run using MPI on two nodes.
If mpiprocs is not specified then PLUMED will, by default, run the test one node.

- The second line tells PLUMED that the test is to be performed using the command line tool <a href="../../user-doc/html/driver.html">driver</a>.  
There are four options that you can use here:
type=driver tells PLUMED to run the test by reading in a trajectory and analysing it using <a href="../../user-doc/html/driver.html">driver</a>, 
type=simplemd tells PLUMED that the test will involve
running some Lennard Jones MD using <a href="../../user-doc/html/simplemd.html">simplemd</a>, 
type=sum_hills tells PLUMED that the test will involve summing gaussians using the 
<a href="../../user-doc/html/sum_hills.html">sum_hills</a> untility
 and type=make tells PLUMED that a main function (called main.cpp) is included in the regtest directory and that the test should compile this function, link in
the PLUMED library and run the resulting executible. The vast majority of features are best tested using type=driver here although the type=make function
can be useful for testing tools that are used in many places in the code (see the code in regtest/basic/rt-make-4 for an example of how use type=make to test
PLUMED's internal matrix class).  If you need to test PLUMED in a way that is not included in these various types please contact the PLUMED developers before
making any changes.

- The third line - the one containing the string plumed_modules=adjmat - tells PLUMED that this regtest needs the optional module adjmat to be installed in 
order for it to be run.  Before trying to run this test PLUMED will thus check whether or not the module is installed.  If it is not installed then it will 
abandon the attempt without running the test.  <b> You must incorporate a line like this if you are testing some feature that is included in an optional module </b>

- The fourth line specifies the command line arguments that need to be passed to the PLUMED tool.  In the above example the final command the test will 
run is thus as follows:

\verbatim
plumed driver --plumed plumed.dat --trajectory-stride 50 --timestep 0.005 --ixyz diala_traj_nm.xyz --dump-forces forces --dump-forces-fmt=%10.6f
\endverbatim    

- The final line tells PLUMED that a few files need to be copied from elsewhere in the PLUMED directory hierarchy in order to run the test.  This copying of files
is useful as it ensures that our repository does not grow large because the same trajectory is duplicated in many regtest directory.  However, if all the files you 
need to run the test  are contained in the test directory this line is not necessary.  Having said that, however, it may be useful to reuse the trajectories contained 
in the directory regtest/trajectories in creating your own regression tests.  

\subsection reffiles Creating the reference files

The next step in creating your regression test will be to copy suitable trajectories to run the calculation into the test directory and to write a plumed input file.
Obviously, we cannot tell you what to do here (although there are some suggestions in the next section) as the tests you will need to write will depend on the feature
you are incorporating.  <b>In general though it should be possible it should be possible to run the test in a few seconds and the output files should be small.</b>
Once you have created this input you can run the test by simply executing 

\verbatim
make
\endverbatim

inside the directory you have created for the regression test.  When this process completes a new directory called tmp is created, which contains the output from your
calculation.  Regression tests work by comparing this output to a reference output and you should thus create this reference now.  This is simple to do you simply copy 
the output file named <i>name</i> from the tmp directory to the regtest directory and rename it <i>name</i>.reference.  As an example suppose we were testing the following
PLUMED input by analysing a trajectory using driver:

\verbatim
d1: DISTANCE ATOMS=1,2
PRINT ARG=d1 FMT=%8.4f FILE=colvar
\endverbatim 

This input creates an output file called colvar so in future we can check the code is running in the same way by comparing with a reference version of the colvar file.  To
create this reference file we execute the following command:

\verbatim
cp tmp/colvar colvar.reference
\endverbatim

Immediately after our first successful run of the test.

\subsection regtrepo Adding the regtest to PLUMED's git repository

Suppose you use the instructions outlined in the previous section to create a regression test in a directory called rt-mynewtest.  You can add this test to the git 
repository by running the commands:

\verbatim
git add rt-mynewtest
git commit
\endverbatim

in the directory containing rt-mynewtest. By adding the test in this way you will ensure that your feature is not broken in future versions of the code.

\section reg-suggestions-1 Some suggestions on how to test new collective variables

If the new feature you are implementing is a new collective variable or a new function of a collective variable and you inherited from either 
\ref PLMD::Colvar, \ref PLMD::function::Function or \ref PLMD::multicolvar::MultiColvarBase then there is a reasonably well established way of testing
that your implementing is correct.  In essence you have to ensure that the numerical derivatives and the analytical derivatives of the variable 
are equal.  This sort of test is easy to setup with PLUMED.  To explain how to setup the test let's suppose for a moment that you want to test 
\ref DISTANCE.  You would do this using an input file something like the one below:

\verbatim
d1: DISTANCE ATOMS=1,2
d1n: DISTANCE ATOMS=1,2 NUMERICAL_DERIVATIVES
PRINT ARG=d1 FILE=colvar FMT=%8.4f 
DUMPDERIVATIVES ARG=d1,d1n FILE=deriv FMT=%8.4f
\endverbatim 

The input above outputs the analytical and numerical derivatives of the distance between atom 1 and 2 to a file called deriv.  The analytical 
derivatives will be in the third column of this file while the numerical derivatives are in the fourth.  If you had run a regtest with this 
input you can thus ensure that the two sets of derivatives were equal using the command:

\verbatim
awk '{print $3-$4}' tmp/deriv
\endverbatim

Lets suppose the derivatives are equal. You might thus conclude the creation of your regtest by copying the deriv and colvar file to files called 
deriv.reference and colvar.reference and checking in the regtest directory.  We would urge you not to do the following, however, as we have found 
that numerical derivatives calculated on different architectures vary, which makes it difficult to determine if the regtests are failing for real 
or if it is just small numerical errors in the calculated numerical derivatives.  On top of this calculating numerical derivatives is expensive.
We would thus ask that once you have ensured that the numerical and analytical derivatives match that you <b> check in deriv files that contain the 
analytical derivatives only. </b>  Also notice the use of the FMT keywords in the above input.  This is important - there will be small numerical 
differences when the code is run by different architectures.  These small differences are not a problem, however, if the format to use 
for real numbers in the output is stated explicitly in the input using FMT.

Notice that this method of testing will also work if you are testing a \ref PLMD::Function.  As an example we can test \ref COMBINE
as follows:

\verbatim
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
c1: COMBINE ARG=d1,d2 PERIODIC=NO 
c1n: COMBINE ARG=d1,d2 PERIODIC=NO NUMERICAL_DERIVATIVES
PRINT ARG=c1 FILE=colvar FMT=%8.4f
DUMPDERIVATIVES ARG=c1,c1n FILE=deriv FMT=%8.4f 
\endverbatim

Now, however, the numerical and analytical derivatives in deriv are the derivatives of the value c1 with respect to the distances d1 and d2.  In
other words the file deriv does not contain derivatives with respect to atomic positions.

\section reg-suggestions-2 More complicated tests

If your feature is more complicated than those covered in the previous section - if for example you are using multiple lines in the PLUMED input
to calculate a final CV or if you have had to write a new apply method in your method - then you will need to do more to test your implementation.
It is probably still a good idea to test analytical derivatives against numerical derivatives in these complicated cases, however.  You can do this
by adding a bias to the system and by exploiting the --debug-forces feature that is contained in <a href="../../user-doc/html/driver.html">driver</a>.  
To understand how this works suppose you have a PLUMED input that looks like this:

\verbatim
d1: DISTANCE ATOMS=1,2
RESTRAINT ARG=d1 AT=1 KAPPA=10
\endverbatim

This input applies a harmonic restraint on the value of the distance between atoms 1 and 2.  If you use the above as input to a driver command such as the one
below:

\verbatim
plumed driver --ixyz trajectory.xyz --debug-forces forces.num
\endverbatim

Then a file called forces.num will be output that contains three columns. The first column is a numerical indexr. The second contains the forces on the atoms
in your system calculated using the analytical equations for the forces that are contained within PLUMED.  The final column also contains the forces on the atoms
in your system but this time these forces are calculated numerically using finite differences.  Calculating the numerical derivatives using --debug-forces is 
very expensive so you definitely should not include regtests that exploit this option in PLUMED.   Instead you can just output the analytical forces using a command
such as:

\verbatim
plumed driver --ixyz trajectory.xyz --dump-forces forces --dump-forces-fmt %8.4f 
\endverbatim

Notice, however, that with both these commands, even thought there are only forces on the 6 coordinates that are used to specify the positions of 1 and 2 
and the cell vectors, the forces with respect to all the atomic positions are output.  In other words, the output files produced usig --debug-forces and 
-dump-foces will contain a lot of zeros. 

\page HowToContributeToPlumed How to contribute new functionality to PLUMED

We welcome researchers to contribute new functionality to the PLUMED code.  In fact, we would 
argue that new biasing methods, cvs and analysis tools should be implemented in common software 
libaries and that the community should adopt shared standards when it comes to the format of input files
and output files in order to make collaboration and communication more straightforward.  There are 
some (not particularly onerous) caveats, however. We would thus ask any person who is 
considering contributing some functionality to PLUMED to read the following page carefully before 
commencing. 

\section maintanence Our phillosophy on sharing code

Writing programs is not so difficult.  Writing software that other people can use and that can
be maintaned is a lot of work, however.   In fact much of the work that we (the core developers)
do to maintain PLUMED does not involve coding.  It involves ensuring that the code has a sufficiently
large test suite, ensuring that old features are not broken when new methods are added to the code,
updating the manual, answering users questions on the mail list and organising PLUMED
meetings and tutorials.  None of the core developers are employed to do these things full-time and as 
such we would like to minimise the amount of time we spend on this sort of maintence.  As such it is
very important to us that we <b>ensure that we understand who wrote every file within PLUMED</b>.  There are
two reasons why we believe this is important:

- It is important that credit is given where credit is due.  In other words, we want users to know who
wrote the functionalities they are using and not to simply assume that we (the core developers) wrote 
everything in PLUMED.  
- The author of any feature is responsible for ensuring their feature is maintained in perpetuity.  
That is to say we (the core developers) will not take responsibility for maintaining features that we
did not ourselves write.  Furthermore, this maintenance involves answering users questions on the mail
lists as well as the occasional fix to the relevant cpp files.

In addition, to these two requirements we do not want the core-developer team to expand as we believe that
if it does expand our meetings will become unweildy.  

We decided on these priorities early during the 
development of PLUMED 2 as we learnt a lot of hard lessons about software management based on our experience 
with PLUMED 1.  We have thus always thought of PLUMED 2 as consisting a small set of <b>core</b> functionalities
together with various extensions.  This led us to write the code so that developers <b>can implement new 
functionalities without changing any other files within PLUMED</b>. 

We believe this distinction between core code and extensions divides the community of PLUMED developers into two
distinct groups.   There is a small group of core-developers who work to ensure that the core of the code is maintained
and a second larger group of contributors who work on new features.  We should be clear, however, that we do not
make this distinction because we believe that what we (the core developers) do is more valuable than the work 
of contributors.  The core developers are only called the core developers because they maintain the core parts of the code
that are used in every PLUMED calculation.  This maintence work is not even software development per say as
<b>the core of the code has not changed much since the publication of the original paper.</b>  In fact we consider almost every
addition we have made since publication and even some of the functionalities that were described in the original papers to 
be extensions upon the original core code.

With all the above in mind we ask that developers who wish to contribute features to PLUMED put any new cpp files in a 
separate module (see what follows for an explanation of how to create a module).  Authors of these modules can use 
different copyright information at the top of their source code files to make it clear that they (and not the core 
developers of PLUMED) wrote these functions.  Furthermore, information about the modules that have been contributed to 
PLUMED will be put here.  On this page information on the authors of the module, the relevant papers and a graphical
image to illustrate the purpose of the module will be provided for the PLUMED users to see.  

Ultimately we would like PLUMED to be a community code that serves the requirements of all the users and developers in this field.
We feel that the best model for achieving this is to have a code that is composed of a number of semi-autonomous modules with their 
own individual identities that are developed in separate research groups.  These modules should (as much as possible) make use of 
a common input/output syntax and a common interface with the various large MD codes.  Furthermore, it should be possible to use 
functionalities from different modules concurrently.  Dividing the code into a core and extensions is what will allow us to achieve 
this aim.       

\section cmodule Creating a module 

In the following sections a set of step-by-step instructions for creating a new module for incorporating some new functionality 
within PLUMED is provided.

\subsection forking Fork the PLUMED repository

Git is amazing!  If you learn to use it properly it will make your life much easier.  It is perhaps hard to get started but it is
well worth the effort and there are some excellent tutorials online e.g. https://try.github.io/levels/1/challenges/1.  The first 
things you should do when you start working on developing PLUMED is to learn a bit about git, create an account on github 
https://github.com and create your own fork of plumed on github (see \ref https://help.github.com/articles/fork-a-repo/).  
By forking plumed on github you are creating your own independent repository on a remote server.  You can change this repository
without breaking the main plumed2 repository as it is your own personal repository.  Furthermore, because it is on a remote server
it is easy for you transfer your code between your laptop and desktop computers.  

The main advantages of having your own fork are that you can merge changes to the main
plumed2 repository into your own repository so that you have all the new features added by other developers using the following
instructions (https://help.github.com/articles/syncing-a-fork/).  Furthermore, as we shall see in later section working in a fork makes it 
straightforward to merge your changes (once they are ready) into the main plumed2 repository and the release version of the code. 

\subsection cdirectory Create a directory for the module source code

Once you have your own fork of PLUMED you can begin to add your new features.  Ideally when you do so you should not need to modify
any of the cpp files that are already part of PLUMED.  In other words, all your features should be implemented in new cpp files and new
header files.  It is a good idea to bundle all these files together into a single directory as has been done for other PLUMED modules 
such as crystallisation, adjmat and metainference.  Within your cpp files it is also a good idea to put all your new code into its own
separate namespace within the main PLMD namespace.  By doing so you prevent naming conflicts with other developers.  

Notice that conflicts may still happen if you pick a name for your collective variable that collides with something
else existing in the code. In this respect, before merging your contribution, the core developers may ask you to change
the name of some of the keywords that you added.

The first step in writing your new feature will thus be to create a sub-directory within src in which to hold your new feature.  Before you 
write any c++ code you will need to create two files within this directory. The first of these files will be called module.type and will 
just contain the following text:

\verbatim
default-off
\endverbatim

By setting this file up this way you ensure that your module is not compiled unless the user specifically asks for it to be compiled at
configure time.  Explanations on how to configure PLUMED to include your module will appear in the documentation for the module automatically.
Furthermore, the procedure for compiling with your code enabled is relatively straightforward.

The second file you will need to create is the Makefile.  An example module Makefile is shown below:

\verbatim
USE=core tools vesselbase multicolvar
# generic makefile
include ../maketools/make.module
\endverbatim

You should only need to modify the first line of this file - the line starting with USE=.  This line is used to tell PLUMED at compile time
what other modules are required in order for this module to function.  This module relies on functionality that is contained in the core, tools, 
vesselbase and multicolvar modules and so these modules are all required during compilation of this particular module.  For your new module you 
will most likely always need to use core and tools.  The other modules that are required will depend on what you are implementing.

The last thing you will need to do before you start programming is that you will need to modify the .gitignore file in the src directory in order 
to stop git from ignoring your new module directory.   If you look in the .gitignore file you will see that it reads something like this:

\verbatim
/*

# Only track modules that are part of the plumed 
!/Makefile
!/README
!/adjmat
!/analysis
\endverbatim

If your new module is called newmodule then you need to add a !/newmodule in the .gitignore file in order to prevent git from ignoring the directory.

\subsection cwrite Writing your source code

Obviously, the source code you will write in your modules directory will depend on the particular feature you are implementing
and there is thus little generic advice that we can give at this stage.  We would ask that you observe a few rules, however.  
In particular:

- Please ensure that you fully document all the new features that you add to the code and that you include examples in your documentation.
There is information on how to write documentation for PLUMED on this page: \ref usingDoxygen.  Note that whenever you make a commit that 
changes the documentation for the code you need to add the string [makedoc] somewhere in the commit message to update the website.  
- We ask that you maintain the portability of plumed by only using the STL library and lapack in modifications.
If you need to use any less standard library (e.g. Boost, Sockets) please ensure that your functionality is not 
installed during a default compilation.  Instead add new options to the configure script in order to make it search for 
libraries so as to make compilation straightforward for users who are using/not using your new feature.  There is information 
on how to add complilation options on this page: \ref UsingExternalLibs.  N.B. you should only need to modify the configure
scripts if you are using modules.  Flags for activating your module at configure time will be generated automatically.
- We ask you not to include C++ features that are too new and not supported by the vast majority of compilers.
Until PLUMED v2.3, we were not using any C++11 feature. Since PLUMED v2.4, a C++11-compliant compiler is explicitly
requested and so you can use C++11 features.

\subsection ctesting Writing regression tests 

It is really important for us (the core developers) to be able to tell if a changes to the code affect the way PLUMED behaves and the answers that 
it produces.  It is thus really important that you write suitable regression tests for new features.  In fact you should be writing regtests as
you code - every time you implement a new feature write a regression tests immediately after you are done.  Testing is an important part of 
developing code and you shouldn't leave it to the end as an afterthought - doing so is horrendously bad software development practise.  

There are instructions as to how to create regression tests on this 
page: \ref regtests.  If you put all your new rt* directories in a single subdirectory in the regtest directory and if you give the directory you
create the same name as your module then that makes our life easier.  Just remember, if you do it this way, that you need to copy the Makefile 
from one of the other module folders to your new directory and that you need to modify the .gitignore file in the regtest folder to make it not
ignore the files in your new folder.  If you look in the .gitignore file you will see it reads soemthing like this:

\verbatim
/*

# Only track modules that are part of the plumed 
!/Makefile
!/README
!/adjmat
!/analysis
!/bias
\endverbatim 

If your module is called newmodule then you need to add a !/newmodule in the .gitignore file in order to prevent git from ignoring the directory.

\subsection pull-request Do a git pull request

Once you have finished coding and once you have written your regression tests submit a git pull request.  There is an explanation of how 
to do a pull request here: https://help.github.com/articles/using-pull-requests/  By doing the pull request you will make the core developers aware
that you want to add a new feature to the release version of PLUMED.  We will be able to see all the changes that you would like to make 
and any new files that you would like to add.  We may come back with a few questions or suggestions at this stage but we will eventually merge your
code into the main PLUMED repository.  Once we have merged your code it is very important that you <b> do not delete your fork of the plumed repository</b>.
If you need to make changes to your source code files you will first need to change them in your own forked repository.  Once you have made these changes
you will then need to <b> submit another pull request</b> in order to change the release version of the code.  Remember also that you are responsible for
maintaining your source code in perpetuity and that we may therefore contact you and ask you to fix something in response to some change we have made.  
In case you will not have time to do it, we will be forced to delete your module from the main PLUMED repository. People
will still be able to download it from your fork, but obviosuly your module won't benefit anymore from the further
enhancement we will add to the main PLUMED repository.
As always this fix should be done on your fork of the repository first and then merged into the release version of PLUMED using a pull request. 

\page ABriefIntroduction A brief introduction to the plumed core

Plumed 2, unlike its predecessor plumed 1, which was written in plain C, is written in C++.  
C++, unlike C, fortran and many of the other languages that are commonly used in the 
atomistic simulation community, is an object oriented programming language.  As such the way things
are done in plumed may feel unfamiliar to developers in the scientific community who,
if our experience is anything to go by, are more used to programming in non-object
oriented languages.  For this reason we have tried in what follows to explain how
we have used the features of object oriented programming in plumed 2.  We
hope that this guide is helpful and appologize in advance to any developers
who feel patronized.     

\section intro Object oriented programming

The main objective in object oriented programming is to write code that is more 
resilient to bugs.  There are two ways that object oriented programing allows
us to acchieve these aims:

- In object oriented programs one generally needs fewer lines of code
- Object oriented programming allows us to use the compiler to do many more checks of the code 

To be clear though object oriented programming does not allow us to do things that 
would be impossible with other programming languages.  All programs perform some set of 
mathematical operations and any programming language can be used to implement these mathematical 
operations. The only advantage of C++ is that the advanced, object-oriented features of the 
language make implementing things more straighforward. 

As you are no doubt aware, in C one can create structures to order the variables in your code.  
A naive way of thinking about the objects, or more correctly classes, that one uses in C++ is
that these are structures that also contain functions.  This is useful for making neat header files
as the parameters are kept near the functions.  However, at this level of thinking the C++ way 
of doing things:

\verbatim
class fclass {
bool param0;
int param1,param2;
double param3;
double f( double x );
}; 
\endverbatim

is not much better than the C way of doing things:

\verbatim
struct fparams {
bool param0;
int param1,param2;
double param3;
}; 
double f( double x, struct fparams myparams ); 
\endverbatim

<b>
Nevertheless for reasons that will hopefully become clear as you read this document every bias, 
colvar and function that is implemented in plumed 2 is inside its own separate class.
</b>

\section publicprivate Public and private members

The variables in a C struct can be accessed anywhere in the code - any function in the code
can copy the information from a structure or change the values of the variables in the structure.  This
was a particularly fun feature of plumed 1.0, every function in the old code could change any of the variables
in the code!  Obviously this causes problems as new functions can accidentally change the 
values of some variable in a widely used structure that should never have been changed.  As one can imagine this can cause 
drastic problems.  To prevent errors like this C++ provides a set of functionalities to 
allow one to specify what any given function can do to the members of a class.  This is not possible in 
C and it was the ability to use this functionality to create flexible, easily-extendable code that motivated
our choice of C++ for plumed 2.  The example class below shows how this is done in practise:

\verbatim
class myclass{
private:
  bool b;  //--This can only be accessed when you are in one of the methods in the class
public:
  int i;  //--This can be acessed by anywhere in the code
};

\section constructors Constructors

As someone who learnt to program in fortran it was this aspect of C++, more than any other, that confused me 
the most. In actually though it is rather simple and I don't really know why I was confused. In essence every 
class must contain some method for creating it. This class should set the initial values of all the variables
inside the class. Obviously the functions (or more correctly methods) that the class
contains cannot be used until an instance of the class has been created using the constructor.

An example of how all this works in practise is given below:

\verbatim
class myclass{
private:
  bool b;
  int i;
  double d;
public:
  myclass( bool bb, int ii, double dd ) : b(bb), i(ii), d(dd) {}
  static double g(int j);
  double f( double x );
};

// Here there are currently no instances of myclass and so I cannot run f
// I can run g however as it is a static member of the class - I run it using
double d=myclass::g(j);

// Now I create an instance of the class using the constructor
myclass thisIsTheInstance(false, 3, 6.5);
// so that I can run the method called f 
double s=thisIsTheInstance.f(4.0); 
\endverbatim         

<b>
In plumed 2 all the lines in the input file are read in inside the constructors.  This ensures
that the parameters inside any given method are set correctly from the outset.
</b> 

\section operators Operators

Addition, subtraction, multiplication, division and so on are all functions (they are obviously
not variables).  We usually don't think of them as functions however because we use 
these operations all the time.  C++ recognizes that the short cuts of +, -, *, / and so on
 are very useful.  It thus allows one to define operators in our new classes that explain
to the compiler what a given symbol means for a given class.  Among other things we can define:

- How to perform the logical operators !=, ==, etc
- How to perform arithmatic for a class: +, -, /, *, +=, -= etc
- What brackets mean i.e. the meanings of (), [] 

We do not use this extensively in plumed 2 but it does occasionally appear.

\section inclusion Including the functionality of one class in a new class 1: Inclusion

There are various ways that one can include the functionality of one class inside a second class.
By far the simplest is to create an instance of class 1 inside class 2 as shown below:

\verbatim
class class1 {
private:
  double d1,d2,d3;
public:
  class1();
  void f(double x);
};

class class2 {
private:
  class1 myclass;
public:
  class2();
  // The methods of class 2 here
};
\endverbatim

This is really simple one includes a class in this way in exactly the same way that one includes a 
double, int or whatever variable.

<b>
This kind of inclusion is used extensively in plumed 1.0 and there are a huge number of classes that you
can re-use in this way to create new colvars, functions or biases.  For a full list of the classes that are 
available see \subpage TOOLBOX.
</b>

\section inheritance Including the functionality of one class in a second class 2: Inheritance

There is an alternate way of reusing the functionality from one class in a second class that is available
in C++.  This method is called inheritance and it offers some advantages over simply including class A
inside class B as described above.  To create a class called B that inherits from A one writes:

\verbatim
class B : public A {
// Contents of class B
};  
\endverbatim

One advantage of this method over inclusion is that I can use protected members to more closely control
what members of A class B is allowed to access.  Hence, rather than simply having private and public members
I now have:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td> <b> public </b> </td> <td> These members can be accessed by anyone </td>
</tr> <tr> 
<td> <b> protected </b> </td> <td> These members can only be accessed by the methods of class A and class B </td>
</tr> <tr>
<td> <b> private </b> </td> <td> These members can only by accessed by the methods of class A (not by class B) </td>
</tr>
</table>

In addition, I can use inheritance to treat pointers to objects of class B as if they were pointers to objects of 
class A.  In other words, if I create an object of type B I can convert it to an object of type A using dynamic_cast
as shown below:

\verbatim
B* mynewB=new B();   // This is a special way of calling the constructor so you get a pointer
A* BpretendingToBeA=dynamic_cast<A*>(mynewB); 
\endverbatim

All the colvars and free energy methods of plumed use inheritance. In fact all these methods are
built on a single base class called PLMD::Action. This class contains all the functionality for
reading stuff from input, the stuff for controlling the dependencies Actions and a set of controls
that decide which actions are performed when. All the functionality for the different methods is
then built on this root. As you can see (PLMD::Action) the inheritance tree for the code is
quite complicated.  However, in practise if you are implementing a CV, function, bias or virtual atom
the correct start point is with one of the classes listed on this page \ref INHERIT all of which contain
detailed descriptions of how to use them. 

\section minheritance Including the functionality of one class in a second class 3: Multiple inheritance

Immediately above the PLMD::Action root of the inheritance tree in plumed there is a very complicated looking
layer in the inheritance structure of the code.  This layer looks ugly because in this layer
we are using multiple inheritance - the classes in the layer above inherit from multiple classes simultaneously.
This way of incorporating functionality from classes is unique to C++ and brings with it a special set of
 difficulties in programming.  Its great advantage is though that one can 
create classes that incorporate bring a set of particular attributes.  This will perhaps be most clear
if you look at what each of the classes in the multiple inheritance layer is doing (see \ref MULTIINHERIT) and see how these
functionalities are used in Colvars, Functions and Biases.  Please be aware that, unless you are doing something really wacky, 
you should be able to implement whatever you need to implement without writing classes that take advantage of multiple inheritance. Furthermore,
you should not need to touch the classes in this region of the code.  The information here is there for
completeness only.  If you feel you really must change something in this part of the code please contact the
developers before doing anything.     

\section static-poly Static Polymorphism

Polymorhpism is a way of using the same code to do many different things.  As an example consider a 
Matrix.  The elements of a Matrix can be ints, doubles, floats or even some fancy new class
but we would still want the operator (i,j) to return the element in row i and column j.  That is to
say the operator (const int i, const int j) of a matrix is independent of what is actually inside the
matrix.  Using C++ we can use so called template classes to implement thse kinds of things and can then 
re-use them to do an enormous variety of different operations.  To see how this works in practise take
a look at PLMD::Matrix, which is a working version of our Matrix example.  Be aware that all the 
routines in a template class must be inside the header file.  To use a template within the code
you declare it as follows:

\verbatim
Matrix<double> mat;   // This is a matrix of doubles
\endverbatim

<b>
The most common way we use this kind of functionality in plumed 2 is when we take advantage of the
features that are available in the C++ standard library.  For more details on the standard library
visit http://www.cplusplus.com/reference/
</b>  

\section dynamic-poly Dynamic Polymorhpism

When you run a calculation with plumed the code calculates a number of CVs.  The bias and the forces
due to the bias are then calculated and in the final step these forces are propegated back onto the 
atoms using the chain rule.  For example PLMD::colvar::Distance
contains the function that calculates a distance between atoms, while PLMD::bias::MetaD contains
the function for doing metadynamics.  What may thus seem remarkable to the programmer unfamiliar with
C++ is that the class that calls the functions that calculate the CVs, biases and so on only uses 
PLMD::Action.  To make that clear it looks like the code can calculate the distances between
atoms without ever calling any of the routines from PLMD::colvar::Distance!      

We can program in this way because we take advantage of dynamic polymorhpism.  If you look at the
documenation for PLMD::Action you will see that the method PLMD::Action::calculate is declare inside
PLMD::Action as:

\verbatim
virtual void calculate()=0;
\endverbatim

This kind of declaration promises two things to a class:

- That the class will only ever be used in derived classes.  No PLMD::Action class is ever 
constructed in the code. The functionality in PLMD::Action is only ever used in the derived classes 
that inherit PLMD::Action. Classes like PLMD::Action are called abstract base classes.
- That in one of the classes that inherits from PLMD::Action a method called calculate will be defined.

The great advantage of declaring calculate() inside PLMD::Action in this way is that the calculate 
routine that we declare in the derived class is a member of PLMD::Action.  We thus can thus write
a class for doing all the business of plumed in the manner described previously.

\section ForwardDeclaration Forward declaration

One problem of including classes inside other classes in C++ is that this enforces one
to include one .h file into another one, thus leading to a large set of objects
needing to be recompiled just because a single .h file was touched. In some cases
this is not avoidable, e.g. when classes inherits from other classes. However,
when only a pointer (or a reference) to another class is used, it might be better
to just use a forward declaration as in this example:
\verbatim
/////////////////////////////////////////////
// This is file A.h
namespace PLMD{

class A{
  int pippo;
};

}

/////////////////////////////////////////////
// This is file B-bad.h
// it has to include A.h
#include "A.h"
namespace PLMD{

class B{
public:
// notice that here we only use a reference to class A
  int do_something(A&a);
};

}

/////////////////////////////////////////////
// This is file B-good.h
namespace PLMD{

// this command just instructs the compiler that A is a class:
class A;
// no inclusion of A.h is required!

class B{
public:
// notice that here we only use a reference to class A
  int do_something(A&a);
};

}

\endverbatim

This trick however does not work is a class is including an instance of another class.
E.g., if B _contains_ an instance of A one should know exactly the A declaration
to build a B object.
In this case, a similar effect can be obtained at the price of adding some
more lines of code in constructor and destructor of B as in the following example

\verbatim

/////////////////////////////////////////////
// This is file B-bad.h
// it has to include A.h
#include "A.h"
namespace PLMD{

class B{
  A content;
};

}

/////////////////////////////////////////////
// This is file B-good.h
namespace PLMD{

// this command just instructs the compiler that A is a class:
class A;
// no inclusion of A.h is required!

class B{
// As "content" is a reference, it can be used exactly as if it was a normal object
// However, it is represents by the compiler as a pointer.
  A& content;
public:
  B();
  ~B();
};

}

/////////////////////////////////////////////
// Using B-good.h enforces to add something in B-good.cpp

#include "A.h"
#include "B-good.h"

using namespace PLMD;

B::B():
// now "content" needs to be explicitly allocated ...
  content(*new A){
}

B::~B(){
// ... and deallocated
  delete &content;
}

\endverbatim

Notice that this trick cannot be always be applied, e.g., if the constructor to be A needs parameter,
or if object "content" is to be accessed by inline methods of B for efficiency. Another example where
this does not work is when inline methods are used because of template expressions.

This trick is extensively used in plumed so as to avoid too many indirect dependencies
among .h files.

\section cxx11features C++11 Features

Since PLUMED 2.4 we systematically use C++11 features. Some of the most important ones are discussed here.

\subsection cxx11features-auto Using auto

Auto allows to implicitly declare the type of a variable.
This is particularly handy when using iterators from STL containers.
For instance, you can replace the following:
\verbatim
   map<string,vector<AtomNumber> >::const_iterator m=atoms.groups.find(strings[i]);
\endverbatim
with:
\verbatim
   const auto m=atoms.groups.find(strings[i]);
\endverbatim
Notice that the syntax is significantly simpler, especially if you do not
remember which exact type of map is the variable `groups`.

Iterators are often used in loops. Thus, you can now replace
\verbatim
  for(std::map<AtomNumber,Tensor>::const_iterator p=gradients.begin();p!=gradients.end();++p){
    a+=(*p).second;
  }
\endverbatim
with
\verbatim
  for(auto p=gradients.begin();p!=gradients.end();++p){
    a+=(*p).second;
  }
\endverbatim
However, in cases where you do not need to explicitly use `p` as an iterator, you might find
even more convenient to use range-based loops:
\verbatim
  for(const auto & p : gradients){
    a+=p.second;
  }
\endverbatim
Notice that now `p` is a constant reference, so it is not anymore necessary to use the `*` operator.

\subsection cxx11features-smart-pointers Using smart pointers

There are many resources on the web about this topic. Have a look at <a href="https://mbevin.wordpress.com/2012/11/18/smart-pointers/"> this link </a> for a
concise introduction.

Smart pointers can be most of the time used in place of regular pointers so as to better manage memory.
In PLUMED you can find many times sections of code such as
\verbatim
  object* obj;
  if(...) obj=new type1;
  else obj=new type2;

  ...

  obj->method();

  ..

  delete obj;
\endverbatim
Here we use a pointer to allow dynamic polymorphism.

In this case, the object pointed by `obj` is not transferred anywhere else.
In other words, you can think that the `obj` pointer owns the object itself.
You can replace it with a `std::unique_ptr` as follows:
\verbatim
  std::unique_ptr<object> obj;
  if(...) obj.reset(new type1);
  else obj.reset(new type2);

  ...
  
  obj->method();

  ..
\endverbatim

Notice that instead of assigning it with `=` you should assign it with `reset()`. This is because
the `std::unique_ptr` cannot be copied and so does not understand the assignment operator.
More importantly, notice that the delete command has disappeared. Indeed, when `obj` goes
out of scope, the pointee is automatically deleted.

You can also use vectors of pointers. Consider the following example
\verbatim
  std::vector<object*> objs;
  for(unsigned i=0;i<10;i++) objs.push_back(new object);

  ...

  for(unsigned i=0;i<10;i++) delete objs[i];
\endverbatim
This can be replaced with
\verbatim
  std::vector<std::unique_ptr<object>> objs;
  for(unsigned i=0;i<10;i++) objs.emplace_back(new object);

  ...
\endverbatim

Notice that instead of using `push_back()` we used `emplace_back()`. The reason is that the latter move
the pointer instead of copying it. More importantly, notice that the delete command has disappeared.
When the vector is cleared, also the contained objects are deleted.

Notice that `emplace_back` needs to be fed with a so-called rvalue. In case you created the
unique_ptr in advance, you should insert it with the following syntax
\verbatim
  std::vector<std::unique_ptr<object>> objs;
  std::unique_ptr<object> to_insert;
  to_insert.reset(new object);
  objs.emplace_back(std::move(to_insert));
\endverbatim

\subsection cxx11features-forward Forward declarations using C++11

Notice that also forward declarations discussed above are a bit simpler
to implement using C++11 syntax. This can be done using a std::unique_ptr
or, even better, using the class ForwardDecl, which is a small utility class that
only implement two methods:
- A constructor, which takes an arbitrary number of parameters and use them
  to construct an internally stored `std::unique_ptr`.
- A `operator *`, which returns a pointer to the object.

An example usage is below:


\verbatim
/////////////////////////////////////////////
// This is file B-good.h
#include "tools/ForwardDecl.h"

namespace PLMD{

// this command just instructs the compiler that A is a class:
class A;
// no inclusion of A.h is required!

class B{
  ForwardDecl<A> content_fwd;
  A& content1=*content1_fwd;
  ForwardDecl<A> content_fwd;
  A& content2=*content2_fwd;
public:
  B();
  ~B();
};

}

/////////////////////////////////////////////
// Using B-good.h enforces to add something in B-good.cpp

#include "A.h"
#include "B-good.h"

using namespace PLMD;

B::B():
// constructors that need no argument can be omitted
  content2_fwd(argument)
{
}

B::~B(){
// empty destructor
}

\endverbatim

Notice that it is necessary to add a destructor, even though it is empty.
The reason is that if the compiler tries to construct an inline destructor for this class
it will not be able to create it (the class is not completely defined in `B.h`.
However, the advantage is that objects are deallocated in the correct order as if they were
normal members of class B, that is the inverse of the initialization order.

\section conc Conclusion

The above is meant to give you some feel as to how plumed works.  If there is stuff you do not
understand it is not necessarily that important.  The great advantage of the code as it is currently
written is that you can implement new methods without touching the core of the code.  So in conclusion
give it a go and have fun!

\section Notes Notes

More information about C++
http://www.parashift.com/c++-faq-lite/

\page UsingExternalLibs Using external libraries in PLUMED

When implementing new CVs and methods you may at time find it useful to make use functionality that has
been implemented in some external library (eg <a href="http://www.boost.org"> boost </a>).  This is strongly 
encouraged and good practise - you introduce fewer bugs if you write less code.  Having said that we would 
prefer it if any functionality that is reliant on external libraries is <b> not </b> enabled by default.  The 
reason for this is that we would like to ensure that the code remains easy to compile.  We would rather not 
have users struggling to resolve lots of dependencies on external libraries just so that they can compile 
PLUMED to run metadynamics with an distance and a torsion.

The first step in ensuring that the code will compile even if your fancy library is not available is to thus 
put all the calls to the functions in these libraries inside an ifdef block.  So for example here is a block
of code that uses the boost graph library:

\verbatim
#ifdef __PLUMED_HAS_BOOST_GRAPH
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#endif
\endverbatim

Here the `#ifdef __PLUMED_HAS_BOOST_GRAPH` ensures that this part of the code is only compiled when the 
`-D__PLUMED_HAS_BOOST_GRAPH` flag is set at compile time.  

This example is a bit of a special case as the particular PLMD::Action I took this example from both when the 
library is available and when it is not.  Obviously, if your PLMD::Action does not work without the library
then you should place everything in the file inside the file (the definition of the class, the definitions of 
the methods and the PLUMED_REGISTER_ACTION command) an ifdef block like the one shown above. 

Once you have written the code and surrounded all the calls to the library by an ifdef block like the one above
you then need to adjust the configure scripts so that the user is provided with the option to compile the code
with your new functionality and the link to the external library available.  To do this you need to edit the
file configure.ac in the plumed2 directory.  The first edit you need to make to this file should read as follows:

\verbatim
PLUMED_CONFIG_ENABLE([library_name],[search for library_name],[no])
\endverbatim

Add this in the part of the file where there are similar names.  This command allows users to enable the package
using `--enable-library_name` when they run the configure command.  Obviously, library_name here should be replaced
with the name of the particular library you are using. A shell variable named `library_name` will be set
to `true` if the library has been requested. The last argument says that, by default, the library is not requested.
If you replace it with a `[yes]`, then you will be able to use `--disable--library_name` to disable the library.

The second edit you need to make to the configure.ac file is as follows:

\verbatim
if test $library_name == true ; then
  PLUMED_CHECK_PACKAGE([header_file.hpp],[function],[__PLUMED_HAS_LIBRARY_NAME],[library_name])
fi
\endverbatim 

This command checks if the library is available on the machine that PLUMED is being compiled on.  If it is the 
Makefiles are generated with the appropriate precompiler directives and lists of libraries to link.  If it is not
available a warning is issued to the user that tells him/her that PLUMED will not be compiled with the requested
features.  The PLUMED_CHECK_PACKAGE command that we use to do this check here takes in four arguments.  The first
is the name of one of the header files from the library that you are using.  You specify the location of this header
file in the same way as you specified its location within the code.  So for example if you had #include <boost/graph/adjacency_list.hpp> 
in your cpp code you would replace header_file.hpp in the above with boost/graph/adjacency_list.hpp.  The next 
argument to PLUMED_CHECK_PACKAGE is one of the functions that is within the library that you are trying to link.  If in doubt
you can use the exit function here, which should work even if you are using template functions.  The third argument is the precompiler
directive that appears around your library calling code and that we discussed earlier.  The last argument meanwhile is the name of the library - 
the part that appears after the -l.  In the example above the code would thus try and link -llibrary_name. 

Notice that in the for C++ libraries sometimes one has to check something more complicated than the presence of a single function.
In this case you can use a more advanced version of the command which allows you to write a short test. Look at this example:
\verbatim
if test $boost_serialization == true ; then
  PLUMED_CHECK_CXX_PACKAGE([boost serialization],[
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
int main() {
    std::ofstream ofs("filename");
    boost::archive::text_oarchive oa(ofs);
    return 0;
}
  ], [__PLUMED_HAS_BOOST_SERIALIZATION],[boost_serialization boost_serialization-mt])
fi
\endverbatim

The first argument here (`[boost serialization]`) is just used in the configure log. The second argument is a
complete C++ program that should compile and link. The last arguments are similar to those of the `PLUMED_CHECK_PACKAGE` macro.

Notice that for both macros (`PLUMED_CHECK_PACKAGE` and `PLUMED_CHECK_CXX_PACKAGE`) you can omit the final library, which
means that the function should be found without adding any extra `-l` option. In addition, you can put multiple libraries.
In this case, autoconf will scan them to find the appropriate one.

Once you have made these edits issue the command:

\verbatim
autoconf
\endverbatim

and the necessary changes will be made to the configure script.  <b>You should never edit the configure script directly</b>

\page AddingAnAnalysis Implementing analysis methods in PLUMED

Information on implementing methods that perform analysis on stored trajectory information eg dimensionality reduction

Implementing methods for analysing trajectory data is more complex than implementing collective variables.
Consequently it is difficult to write a step by step guide like those we have written on implementing \ref AddingAColvar "colvars" or
\ref AddingAFunction "functions".  Hence, this document tries to explain the things that we have considered and the way these 
have been incorporated into the PLMD::analysis::AnalysisBase and PLMD::analysis::AnalysisWithDataCollection abstract base classes.  
Hopefully this will provide some insight into our rationale in writing these parts of the code and will help you to understand how 
any new analysis method can be implemented in the PLUMED code in a way that exploits those features that are already there.

\section overview An overview of analysis methods in PLUMED

There are two distinct ways in which one may wish to perform some form of analysis on a molecular dynamics trajectory.  In the first method some quantity
is calculated for each of the atoms in the trajectory in each of the frames and this is then averaged over the whole trajectory.  Velocity 
autocorrelation functions or mean squared displacements are examples of forms of analysis of this type.   The methods implemented in PLMD::analysis are 
not of this type.  These methods are designed to collect set of snapshots of the trajectory and to perform some analysis of these snapshots.
These trajectory snapshots might be the values of a particular set of collective variables for each of the frames, they might be the 
instantaneous positions of the atoms in each frame or they might be some combination of the above.  The assumption then when running one of these
analysis methods is that a representation (or snapshot) will be collected intermittently from the trajectory and then once a sufficiently large
 collection of these snapshots are collected they will be analysed. 

\section antraj Analysis on the fly

It is important to remember that PLUMED is primarily a code for performing biased molecular dynamics.  The code's original purpose was
to be a plugin that could be easily added to a number of different molecular dynamics engines that allowed additional forces to be 
incorporated when integrating the equations of motion.  The command line tools that allows one to analyse trajectories during post
processing were added later as an afterthought.  This consideration is particularly important when considering analysis algorithms 
in this code because most analysis codes that are used in the community read in the trajectory and to do all their calculations 
during post processing.  The analysis functions that have been implemented in PLUMED can be used to post-process trajectories - you
simply make use of the command line tool driver - however, they can also be used to analyse trajectories on the fly.  We believe this 
is useful for a number of reasons:

- Computers are becoming more powerful so it is possible to run simulations for much longer.  At the same time, however, hard drives
space is at a premium and it is thus difficult to store these large trajectories for post-processing.  If trajectories can be analysed
on the fly this presents less of a problem.
- A number of free energy methods (e.g. Gaussian mixture umbrella sampling, adaptive umbrella sampling and reconnaissance metadynamics)
work by performing a sophistacated analysis of the trajectory and then using the result from this analysis to design a bias for the 
dynamics. 
- Analysis methods implemented in PLUMED can take advantage of the many different collective variables that we have implemented in 
this code and are thus extremely flexible implementations of these techniques.
- Analysis methods implemented in PLUMED make use of the PLUMED input syntax, which hopefully allows users already familiar with 
PLUMED to get to grips with using these tools more rapidly. 

\section genphil General Phillosopy

The aim in the PLMD::analysis and PLMD::dimred modules is to write code that is very flexible and that allows the user a great deal of 
flexibility in the input.  For example we would like to be able to write the code so that a user can collect data from the trajectory and 
then at analysis time they can:

- Select a subset of landmark points from the stored data
- Generate projections of these landmark points using sketch-map
- Project the remaining non-landmark points using the sketch-map projection generated and construct a histogram as a function of the sketch-map coordinates.

Furthermore, we would like to be able to do all the above with a minimum of syntax for the user and with a minimum amount of code to maintain.
This is why the analysis class is structured in the way it is.  The general idea is that one PLMD::analysis::AnalysisWithDataCollection object 
collects the trajectory data as the simulation progresses (for the example above this would be an object of type PLMD::analysis::EuclideanDissimilarityMatrix).
Then when it is time to analyse a chain of analysis objects are run on the data collected by the PLMD::analysis::AnalysisWithDataCollection.
There are thus two types of analysis actions:

- Those that inherit from PLMD::analysis::AnalysisWithDataCollection - these can collect data from a trajectory
- Those that inherit from PLMD::analysis::AnalysisBase - these cannot collect data from a trajectory.  They get their input from another PLMD::analysis::AnalysisBase Action

If you look at the code most of the routines in the PLMD::analysis and PLMD::dimred modules do not inherit from PLMD::analysis::AnalysisWithDataCollection.
In fact the only ones where this is an option that users really see in the manual are PLMD::analysis::Histogram and PLMD::analysis::EuclideanDissimilarityMatrix
(many of the other analysis actions that inherit from PLMD::analysis::AnalysisWithDataCollection only do so because this allows us to write straightforward 
regression tests for this part of the code).  The remaining analysis actions inherit from PLMD::analysis::AnalysisBase because in general they require some 
dissimilarity information in order to function and because this information is calculated within PLMD::analysis::EuclideanDissimilarityMatrix.  
Consequently, if you are writing a new analysis class it should probably not inherit from PLMD::analysis::AnalysisWithDataCollection.

\section storing Storing trajectory data for analysis

As discussed in the previous section the methods in PLMD::analysis all work by storing a snapshot of the trajectory every \f$N\f$ steps
and by then running an analysis method every \f$M\f$ steps, where \f$M\f$ is some substantial multiple of \f$N\f$.  This intermittent
storing of trajectory frames and occasional analysis of the trajectory is all looked after within the PLMD::analysis::AnalysisWithDataCollection
abstract base class.  Users can then set of a chain of analysis actions on the stored data by using actions that inherit from
PLMD::analysis::AnalysisBase.   Any class inheriting from PLMD::analysis::AnalysisBase must have a method within it named performAnalysis(),
which will actually perform the analysis on the stored trajectory frames.  When implementing a new analysis method the majority of your
development time will probably be spent implementing some part of this performAnalysis method.

\section reweight Reweighting trajectories

As discussed in previous sections PLUMED is primarily a code for doing biased molecular dynamics simulations.  This bias is used to force rare events
to occur in the short time scales that are accessible within a molecular dynamics simulation.  Oftentimes when analysing a trajectory we would like
to remove the contribution of the bias and to reweight the simulation so as to get the underlying free energy surface.  When performing any analysis of
the trajectory one may similarly wish to remove the effect of the bias and to consider what weight each sampled point would have had if it had been sampled
in accordance with the underlying canonical probability distribution function.  This process of reweighting points - of ascribing a weight to each snapshot
that discounts the effect of the simulation bias - is again looked after within PLMD::analysis::AnalysisWithDataCollection.  If you wish to take these
weights into account in your analysis method you should use the method PLMD::analysis::AnalysisBase::getWeight to access them.  Obviously, if you have
no simulation bias on the system then each point will have a weight of one and this will be the weight returned by PLMD::analysis::AnalysisBase::getWeight.

\section output Outputting data files

The fact that PLMD::analysis::AnalysisWithDataCollection can be used to run trajectory analysis in either post processing or on the fly during a trajectory means
that this class must also look after a number of things.  For example one might wish to perform multiple analyses of the trajectory 
during a simulation.  Obviously, you do not want to overwrite the output file from your first analysis when you perform the second 
analysis of the trajectory.  In addition, you do not want to overwrite files from earlier runs if you choose to rerun your analysis 
in a directory where you had already run an earlier calculation.  For these reasons whenever you wish to read in the name of an output file 
you should use the following code to make sure that any old files are backed up on restart:

\verbatim
if( !getRestart() ){ OFile ofile; ofile.link(*this); ofile.setBackupString("analysis"); ofile.backupAllFiles(fname); }
\endverbatim 

where fname is the name of your output file. On top of this when you open an output file in your analysis method you should use the following 
set of commands:

\verbatim
OFile gfile; gfile.link(*this);
gfile.setBackupString("analysis");
gfile.open( ofilename.c_str() ); 
\endverbatim

The second line ensures that files are named analysis.0.ofilename, analysis.1.ofilename and so on. Having said that it is probably best 
to not write routines to output data in analysis classes and to simply ensure that you can pass any output data from your method to the 
PLMD::analysis::OutputColvarFile and PLMD::analysis::OutputPDBFile methods.  If you have done everything properly these classes should be
able to interact with your new class through methods that are in PLMD::analysis::AnalysisBase.

\section metrics Calculating dissimilarity matrices

One of the most basic ways in which you can analyse a collection of trajectory snapshots is to calculate the matrix of dissimilarities between each of the pairs
of trajectories frames.  In plumed this is done in by the class PLMD::analysis::EuclideanDissimilarityMatrix.  Notice as well that this class makes full use
of the PLMD::reference module that is discussed in \ref AddingAMetric and so this class alone allows you to calculate the dissimilairity between any pair of 
trajectory frames in a wide variety of different ways.  In addition, you can use PLMD::analysis::ReadDissimilarityMatrix to get the dissimilarities from an 
input file.  There should thus be little need to implement new objects for calculating dissimilarity matrices -- if you really feel you need something other than
what is already there or that is not implementable by \ref AddingAMetric then you are doing a calculation that is very unusual.

\section landmarks Landmark selection algorithms

Many analyses methods scale poorly with the number of trajectory frames that you wish to analyse.  This happens in part because you need to compute the matrix
of pairwise disimiarities (a square matrix in which the number of columns is equal to the number of trajectory frames) but also because you then have to 
do some algebra involving this matrix.  To alleviate these problems a common strategy is to perform the analysis on a set of so-called landmark frames and to 
then project the non-landmark snapshots from the trajectory using some out-of-sample extension of your analysis algorithm.  Classes that inherit from
PLMD::analysis::LandmarkSelectionBase are implementations of the various landmark selection algorithms that are commonly employed.  If your favourite landmark
selection algorithm is not there you may choose to implement a new landmark selection algorithm by writing a new PLMD::Action that inherits from this class.
Within this new class you need to only define a registerKeywords function, a constructor and a method that actually selects the landmarks that will be a function
that must be called selectLandmarks.  Once you are satisfied that you want frame \f$k\f$ in your landmark set then you can select this using 
PLMD::analysis::LandmarkSelectionBase::selectFrame.  Please note that you can get the number of landmarks to select by calling:

\verbatim
getNumberOfDataPoints()
\endverbatim

If you would like to get the total number of frames from which you can get your subset of landmarks you should call:

\verbatim
mydata->getNumberOfDataPoints()
\endverbatim

Lastly, if you are using a method, which like PLMD::analysis::FarthestPointSampling uses the distances between your input points in some way, you should 
add something akin to:

\verbatim
if( !dissimilaritiesWereSet() ) error("dissimilarities have not been calcualted in input action");
\endverbatim

in your constructor as this will ensure that users are told what they are doing wrong if they do not specify how to calculate distances between points in the
input.  To then get the dissimilirity between input point \f$i\f$ and input point \f$j\f$ use:

\verbatim
mydata->getDissimilarity( landmarks[i], k );
\endverbatim

Calling PLMD::analysis::AnalysisBase::getDissimilarity will give you the distance between a pair of landmarks, which is not what you need.

\section dimred Dimensionality reduction

The aim when writing any large code such as PLUMED is to minise the number of lines of code as fewer lines of code means fewer bugs on average.
Hence, as explained in other sections of this developer manual, all the object oriented programming, inheritance and polymorphism.  Given this
consider how we would go about implementing a library of dimensionality reduction algorithms.  In LLE, ISOMAP, sketch-map or MDS the aim is to
generate a low-dimensional projection of some set of high-dimensional data points.  For all these methods we can use the same code to to store 
the high and low dimensional points and to output this data to a file.  In fact the only things that differ in these various different methods are
the ways in which the dissimilarities between the high-dimensional points are calculated and the manner in which the low-dimensional projections
are generated.  We have already discussed how PLUMED calculates matrices of dissimilarities between points using PLMD::analysis::EuclideanDissimilarityMatrix
and how one can think about introduce new methods of calculating dissimilarities.  Priting to files meanwhile can be looked after by PLMD::analysis::OutputColvarFile 
and PLMD::analysis::OutputPDBFile.  Furthermore, if dimensionality reduction classes are written properly it is even possible to pass the projections
generated by them to PLMD::analysis::Histogram and to output histograms and free energy surfaces as a function of the low-dimensional coordinates.
As such in PLMD::dimred::DimensionalityReductionBase and its daughter classes we are thus solely concerned with calculating projections of data points.  
The dissimilarities between the input high dimensional frames will always be calculated by some method akin to PLMD::analysis::EuclideanDissimilarityMatrix.  
PLMD::dimred::DimensionalityReductionBase inherits from PLMD::analysis::AnalysisBase and expects another PLMD::analysis::AnalysisBase object as input.  
Furthermore, with the exception of PLMD::dimred::ClassicalMultiDimensionalScaling the input PLMD::AnalysisBase must be an 
PLMD::dimred::DimensionalityReductionBase as initial guesses must be suggested when using an iterative optimization algorithm such as PLMD::dimred::SMACOF.

Much of the work of dimensionality reduction done in the base PLMD::dimred::DimensionalityReductionBase class.  When implementing any new
dimensionality reduction algorithm your focus will be on writing a PLMD::dimred::calculateProjections routines.  This function will take as
input the matrix of pairwise dissimilarities between points (be they landmark points or otherwise) and is expected to return a matrix 
containing the projections of the high-dimensional data points.  Depending on the stress function you minimise to find projections you 
may also have to implement PLMD::dimred::DimensionalityReductionBase::setTargetDistance and PLMD::dimred::DimensionalityReductionBase::calculateStress
functions.  This is necessary with PLMD::dimred::SketchMapBase for example because sketch-map uses transformed versions of the dissimilarities
and distances in its stress function.  These two routines are used within PLMD::dimred::ProjectNonLandmarkPoints which will product optimal projections of 
points that were not selected as landmarks.   

\section cluster Clustering trajectory data

There are currently no clustering methods implemented in the PLMD::analysis module.  This section is thus here to explain how I (Gareth Tribello) imagined one
might go about implementing these methods.  Again there are many commonalities between methods such as kmeans, Gaussian mixture models and so on, which should 
be thought about when constructing an abstract base class.  Furthermore, this abstract base class should (like PLMD::analysis::DimensionalityReductionBase) 
inherit from PLMD::analysis::AnalysisBase and be implemented in a way that allows one to exploit the use of landmark selection algorithms that inherit from 
PLMD::analysis::LandmarkSelectionBase and the measure actions such as PLMD::analysis::EuclideanDissimilarityMatrix.  It should also be able to work with the 
weights you get by reweiting the trajectories with the bias and so on.  I don't think that you should have this inheriting from 
PLMD::analysis::AnalysisWithDataCollection as I believe you want the clustering to work with projections of data generated by dimensionality reduction algorithms.
Obviously, if you are thinking of adding methods to cluster trajectory frames within PLUMED please feel free to get in touch with me (gareth.tribello\@gmail.com).  
I will be more than happy to discuss these ideas with you.

\page parsing Parsing functionality

By now you are probably familiar with the way that plumed2 input looks:

\verbatim
DISTANCE ATOMS=0,300 LABEL=dist NOPBC
RESTRAINT ARG=dist KAPPA=1.0 AT=1.0
\endverbatim

Within the code this information will be read by either parseVector or parseFlag.  parseFlag is called using:

\verbatim
parseFlag("NOPBC",nopbc)
\endverbatim

This will then read the list of action objects you passed to the constructor and look for the keyword NOPBC.  If the keyword is found then it is deleted from the list of action objects, while the boolian nopbc is returned as true otherwise the boolian is returned as false.  parseVector is called using:

\verbatim
std::vector<double> vec;
parseVector("KAPPA",vec);
\endverbatim

This routine will then read the list of action objects you passed to the constructor and look for the keyword KAPPA.  This keyword will be followed by an equals sign and a list of comma separated numbers.  These numbers are read into the vector vec and passed back to the main code.  (N.B.  The size of the vector is worked out automitically by parseVector from the input.  In addition the vector can be a vector of int or a vector of real.)  Much like parseFlag, parseVector will delete the keyword and the list from the list of actionObjects once it has completed.  When you have finished reading all your arguments you should call checkRead() - this routine checks that everything in the list of ActionOptions taken from the input has been read in. 

Please note when you are implementing functionality to read the plumed input that you never need to implement anything to read ARGS and LABEL as these keywords are read elsewhere in the code. 

\page usingDoxygen Creating plumed documentation

To create the plumed manual you should go to the <b> user-doc </b> directory and type <b> make </b>. 
This command works because user documentation for all the PLMD::Action is inside the source code.  If
you look at the documentation page for any of the actions that are implemented in plumed you will
see that it is composed of three pars:

- A short introduction which describes what the method does.
- A description of the various keywords for the calculation.
- An example/some examples of how the PLMD::Action can be used.

Furthermore, you will also have noticed that if you make an error in the input for any PLMD::Action the 
descriptions of all the keywords from the manual appears in the log file.  This is possible because manual
pages for PLMD::Action are inside the code and because the manual is created using the 
following packages:

- Doxygen:   http://www.doxygen.org
- Graphviz:  http://www.graphviz.org/ 

In addition a special class, PLMD::Keywords, is used to store the descriptions of the syntax for any given
action so that this data can be produced in the output when the user makes a mistake in input.  In the following
a step-by-step explanaition as to how to use the documentation prodcuing functionality of plumed is provided.  When
you have finished creating your new documentation they should be incorporated automatically in the manual the next 
time you make the manual.  Furthermore, a number of checks of your manual are performed when you make the manual
so any errors should be straightforward to find. 

The plumed manual does not only contain descriptions of the various actions that have been implemented and their 
keywords. There are also pages that place these actions in context and which attempt to explain the interplay 
between colvars, functions, biases, analysis methods and so on.  More importantly, the plumed manual contains
instructive examples that explain to users how to use particular features of the code.  We would encourage all
users of the code to submit these sort of things to the plumed repository particularly if they have found it difficult to 
start using a particular method.  Instructions as to how to go about writing a material that will appear in the How-tos 
section of the manual can be found \ref tutorials here

\section registerkeys Registering Keywords

When you implement any new PLMD::Action in plumed you must first create the documentation for the keywords that you
will use to read the input for your new method.  In fact you cannot read in undocumented keywords using plumed.  The
documentation for keywords is created in a static method of the action called registerKeywords.  This method
should be declared in the definition of the class as follows:

\verbatim
static void registerKeywords(Keywords& keys);
\endverbatim

The static attribute allows one to use the registerKeywords routine when no instances of the class have been created.
This is essential as keywordRegistration in plumed is done before the list of PLMD::Action is created.  This means
that when the keywords are created plumed has no understanding of the hierarchy of inherited Actions.  Hence, before
adding your own Keywords you must ensure that the keywords for the class from which your new class inherits have been
added.  In pracise this is done by calling PLMD::Colvar::registerKeywords, PLMD::function::Function::registerKeywords or 
PLMD::Bias::registerKeywords for a new PLMD::Colvar, PLMD::function::Function or PLMD::bias::Bias respectively.  To be clear these
functions will ensure that generic keywords such as LABEL, NUMERICAL_DERIVATIVES or ARG are registered and explained in the
manual. If your method requires the derivatives of some value and you have no way of implementing the analytical derivatives
you should also call PLMD::ActionWithValue::noAnalyticalDerivatives.  This routine will ensure that plumed's numerical
derivatives routines are used to calculation your derivatives automatically and will ensure that a message is put in the plumed
manual so that other users are aware that numerical derivatives are being used. 

Once you have called the reigsterKeywords routine for the PLMD::Action above yours in the hierarchy you can begin to add
the keywords for your new method.  These keywords will have one of 5 attributes:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td width=5%%> <b> compulsory </b> </td> <td> These are the quantities that must be defined in order to perform your action </td> 
</tr> <tr>
<td> <b> optional </b> </td> <td> If there is some alternate way of performing your calculation that requires numerical input you should declare your keyword as optional </td>
</tr> <tr>
<td> <b> flag </b> </td> <td> This is used to declare keywords such as NOPBC that tell plumed to turn on/off some feature of the calculation </td>
</tr> <tr>
<td> <b> atoms </b> </td> <td> If you are reading a list of atoms after the keyword then you should use this keyword. You can easily specify if there are multiple ways of defining the atoms involved in the action in the manual.  To register the keywords for this first method of specifying the atoms using atoms-1. Then register the keywords for the second way of specifying the atoms using atoms-2 and so on.  A manual that states that these keywords can be used in an either or fashion, much that for <a href="../../user-doc/html/_t_o_r_s_i_o_n.html"> TORSION </a>, will then be generated.  </td>
</tr> <tr>
<td> <b> numbered </b> </td> <td> If you need to read in a list of similar keywords such as keyword0, keyword1, keyword2... then you must use this option.  These keywords will be assumed to be optional.  However, you can set them to be atoms or whatever by using reset_style(keyword,newstyle).  </td>
</table>

All keywords (other than flags) are added using the add method of PLMD::Keywords.  This command has the following syntax:

\verbatim
keys.add( attribute, keyword, explanation );
\endverbatim

where <i> attribute </i> is one of the options from the above table, <i> keyword </i> is the word that appears on the input line and <i> explanation </i> is an explantion
of what the keyword does.  If your keyword is compulsory it can also be added using:

\verbatim
keys.add( attribute, keyword, default, explanation );
\endverbatim

where <i> default </i> is a string containing the default value to use for the quantity.

Flags are added using the add flag method, this has syntax:

\verbatim
keys.addFlag( keyword, default, explantion );   
\endverbatim

where default is a bool that tells plumed if by default this option is/is not in use.  

\section reading Reading the input keywords

Keywords are read in using either PLMD::Action::parse, PLMD::Action::parseVector, PLMD::Action::parseNumberedVector or PLMD::Action::parseFlag.  
These routines will use the information provided during keyword registration to check the sanity of any input.  For instance if you declare a 
compulsory keyword and do not specify a default value then the code will automatically complain if the particular keyword is missing from input.  
In addition, if the vector you pass to PLMD::Action::parseVector and PLMD::Action::parseNumberedVector has a
size greater than 0 plumed will assume that the input should contain a vector of this size and will complain if it finds a different sized vector.

\section components Registering components

In plumed 2.1 we will also begin registering all the components that are calculated by your action.
In plumed 2.2 registering components will become compulsory and your features will not work if this is not done.
This registering of components means that in the registerKeywords method for commands such as:

\verbatim
d1: DISTANCE COMPONENTS
\endverbatim

which calculates quantities that can be referenced in the input as d1.x, d1.y and d1.z, you will have to provide
documentation for the manual that describes what information is stored in the x, y and z components. As an example
for the distances components this documentation takes the following form:

\verbatim
keys.addOutputComponent("x","COMPONENTS","the x-component of the vector connecting the two atoms");
keys.addOutputComponent("y","COMPONENTS","the y-component of the vector connecting the two atoms");
keys.addOutputComponent("z","COMPONENTS","the z-component of the vector connecting the two atoms");
\endverbatim

As you can see this new feature works in a similar manner to the feature by which the keywords are registered.  Both
these features serve to keep the manual up to date with only a relatively small amount of effort on the part of the
developer.

All components registered by plumed should be added using a variant of the addOutputComponent method of PLMD::Keywords.
This command has the following syntax:

\verbatim
keys.addOutputComponent( name, keyword, explanation )
\endverbatim

where <i> name </i> is the name the component will be given i.e. the quantity will be referencable in input as <i> label.name </i>.
<i> keyword </i> is the name of any keyword that must be added to the input for the action in order to calculate this quantity
and <i> explanation </i> is an explanation of the quantity that will appear in the manual.

If your Action always generates a particular set of components then the form of this command changes slightly. That is to say
if it makes no sense to use reference the isolated label for your command then this command should read:

\verbatim
componentsAreNotOptional(keys);
keys.addOutputComponent( name, "default", explanation )
\endverbatim

So for example in RESTRAINT, which always generates a component called bias and a component called force2, the components are registered
using the following code:

\verbatim
componentsAreNotOptional(keys);
keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
\endverbatim

Lastly, if you have some method that takes in arguments and gives back one component per argument then these components should be
labelled <i>argument-name</i>_<i>description</i>.  The MOVINGRESTRAINT command below gives an example of how this is done in practise.

\verbatim
DISTANCE ATOMS=1,2 LABEL=d1
DISTANCE ATOMS=3,4 LABEL=d2
MOVINGRESTRAINT ARG=d1,d2 AT0=2,2 AT1=6,6 STEP0=0 STEP1=100 KAPPA=1
\endverbatim

This command has components called d1_steer, d2_steer, d1_cntr and d2_cntr.  These components describe the work done in moving the 
system along the d1 and d2 axis and the instantaneous positions of the harmonic potential on the d1 and d2 axis respectively.

The various components in MOVINGRESTRAINT are registered using:

\verbatim
componentsAreNotOptional(keys);
keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
keys.addOutputComponent("_cntr","default","one or multiple instances of this quantity will be refereceable elsewhere in the input file. "
                                          "these quantities will named with  the arguments of the bias followed by "
                                          "the character string _cntr. These quantities give the instantaneous position "
                                          "of the center of the harmonic potential.");
keys.addOutputComponent("_work","default","one or multiple instances of this quantity will be refereceable elsewhere in the input file. "
                                          "These quantities will named with the arguments of the bias followed by "
                                          "the character string _work. These quantities tell the user how much work has "
                                          "been done by the potential in dragging the system along the various colvar axis.");
\endverbatim

Appropriately labelled components can be created using the following command:

\verbatim
comp=getPntrToArgument(i)->getName()+"_cntr";
addComponent(comp); componentIsNotPeriodic(comp);
\endverbatim

They can then be set by using something like:

\verbatim
getPntrToComponent(getPntrToArgument(i)->getName()+"_work")->set(val[i]);
\endverbatim

\section reserved Reserved Keywords

To maintain some consistency for end users of the code certain keywords (e.g. ARG, PERIODIC) are reserved.
The reserved keywords for PLMD::Colvar, PLMD::function::Function and PLMD::bias::Bias are explained inside the documentation for
these actions.  To use one of the registered keywords you shold insert the following command into the registerKeywords
method of your new function.

\verbatim
keys.use( keyword );
\endverbatim

where <i> keyword </i> is a string that tells the method which reserved keyword you wish to use.  To be clear
when you use a reserved keyword all the parsing and registering for it is looked after automatically. In addition,
if new components are generated in the output of the action when the keyword is present the registering of those
components is also looked after elsewhere.  So, for example, if you do something like:

\verbatim
keys.use("MIN")
\endverbatim

In a new PLMD::multicolvar::MultiColvar then there is no need to add the command:

\verbatim
keys.addOutputComponent("min","MIN","the minimum value. This is calculated using the formula described in the description of the "
                                    "keyword so as to make it continuous.");
\endverbatim 

As this documentation will already have been created for you elsewhere in the code.

\section errors Generating errors

You may need to check for other mistakes in input.  When you find these mistakes you can report them to users using PLMD::Action::error.  This routine will 
output a description of how the input to your Action should appear.  It takes as input a string that describes the particular nature of the error that the user has made.

\section manual Creating the rest of the manual

The remainder of the manual - the detailed description of your action and some examples of how the PLMD::Action can be used - is created 
from the comments at the top of the cpp file that contains the various subroutines that your PLMD::Action performs.  This is converted
to manual using Doxygen as this allows one to incorporate equations, bibliographic information and bits and pieces of HTML.  At the
start of the block of manual information the following lines should appear:

\verbatim
//+PLUMEDOC TYPE ACTIONNAME 
/*
\endverbatim

ACTIONAME is the first word that appears on the input line - i.e. it is the command that a user would use in input in order to make
use of your particular PLMD::Action.  TYPE, meanwhile, tells Doxygen where in the manual the Docuementation should be placed.  TYPE
should be one of the following:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> 
<td width=5%%> <b> COLVAR </b> </td> <td> This is used if your PLMD::Action is the calculation of a CV </td>
</tr> <tr>
<td width=5%%> <b> DCOLVAR </b> </td> <td> This is used if your PLMD::Action is a colvar that measures the distance from a reference frame in some metric </td>
</tr> <tr>
<td width=5%%> <b> MCOLVAR </b> </td> <td> This is used if your PLMD::Action calculates some Function of a distribution of CVs. (Action inherits from PLMD::multicolvar::MultiColvar)  </td>
</tr> <tr>
<td width=5%%> <b> MCOLVARF </b> </td> <td> This is used if your PLMD::Action calculates a paricularly complicated function of a distribution of CVs. (Action inherits from PLMD::multicolvar::MultiColvarFunction) </td>
</tr> <tr>
<td width=5%%> <b> FUNCTION </b> </td> <td> This is used if your PLMD::Action calculates some Function of a set of CVs </td>
</tr> <tr>
<td width=5%%> <b> VATOM </b> </td> <td> This is used if your PLMD::Action calculates the position of a new atom e.g. for COM </td>
</tr> <tr>
<td width=5%%> <b> ANALYSIS </b> </td> <td> This is used if your PLMD::Action does some analysis of the trajectory </td>
</tr> <tr>
<td width=5%%> <b> BIAS </b> </td> <td> This is used if your PLMD::Action is a bias that adds supplemental forces to the potential in order to enhance sampling </td> 
</tr> <tr>
<td width=5%%> <b> TOPOLOGY </b> </td> <td> PLMD::setup::MolInfo has documentation of this type.  This command is used to provide information about the chemistry of the system under study.  For MolInfo this is what constitute the backbone atoms of the protin, what the residues are etc. </td>
</tr> <tr>
<td width=5%%> <b> GENERIC </b> </td> <td> This should be used if you want to specify manually where in the manual your documentation should appear.  If you feel this really is the correct way to incorporate your new feature please contact the core developers so as to discuss it. </td>
</tr> 
</table>

Immediately after the start comment symbol you should place a single line that describes in a sentence or two what it is your PLMD::Action does.  This information will
appear beside the link to your more detailed manual page in the general pages of the user manual.  The code will use everything up to the first blank
line in input to create this brief description.  You can then write a longer description of your PLMD::Action to appear at the start of its
particular page in the manual.  As described below this description can incorporate equations and bibliographic information.

\subsection Equations

You can add formulae in latex using:

\verbatim
This is an inline equation \f$s=y+x\f$ but this is an equation:

\f[
r = \sqrt{ \mathbf{s}^T \mathbf{C}^{-1} \mathbf{s} }
\f]

And this is an equation array:

\f{eqnarray*}{
 f &=& \frac{1}{2} \\
 g &=& \frac{2}{3}
\f}
\endverbatim

In the manual this will be translated into:

This is an inline equation \f$s=y+x\f$ but this is an equation:
 
\f[
r = \sqrt{ \mathbf{s}^T \mathbf{C}^{-1} \mathbf{s} }
\f]

And this is an equation array:

\f{eqnarray*}{
 f &=& \frac{1}{2} \\
 g &=& \frac{2}{3}
\f} 

\subsection Lists

You can create lists of data using: 

\verbatim
- First item in list
- Second item in list
\endverbatim

which becomes:

- First item in list
- Second item in list

\subsection Formatting 

You can create a new section in your documentation using:

\verbatim
\section manual Creating the rest of the manual
\endverbatim

In fact I used this very command earlier in writing this page.  I can therefore reference it here (\ref manual) by using:

\verbatim
\ref manual 
\endverbatim

You can also reference external webpages by typing web addresses directly in the documentation.

\subsection Citations

You can create citations using:

\verbatim
\cite bibtex-tag
\endverbatim

This command uses an interface between Doxygen and bibtex to create bibliographic data.  Inside
the user-doc directory you will find a bibtex file called bibliography.bib that contains all
the references that are included in the user documentation for plumed.  To add your reference
you should add bibliographic data for the article you want to cite in this file.

\section Examples

Manual entries for actions and tutorials <b>must</b> contain some examples.  The most basic way to include these is as follows:

\verbatim
\par Example

The following input tells plumed to print the distance between atoms 3 and 5,
the distance between atoms 2 and 4 and the x component of the distance between atoms 2 and 4.
\plumedfile
DISTANCE ATOMS=3,5             LABEL=d1
DISTANCE ATOMS=2,4 COMPONENTS  LABEL=d2
PRINT ARG=d1,d2,d2.x
\ endplumedfile /*** But with no space between the \ and the endplumedfile
\endverbatim 

In the manual this will be converted to:

\par Example

The following input tells plumed to print the distance between atoms 3 and 5,
the distance between atoms 2 and 4 and the x component of the distance between atoms 2 and 4.
<pre class="fragment">
<a href="../../user-doc/html/_d_i_s_t_a_n_c_e.html" style="color:green">DISTANCE</a> ATOMS=3,5             LABEL=d1
<a href="../../user-doc/html/_d_i_s_t_a_n_c_e.html" style="color:green">DISTANCE</a> ATOMS=2,4 COMPONENTS  LABEL=d2
<a href="../../user-doc/html/_p_r_i_n_t.html" style="color:green">PRINT</a> ARG=d1,d2,d2.x
</pre>

Please be aware of the blank line between after the title of the paragraph.  If this line is not present your manual will look ugly.  
Also be aware that your Examples section <b> must </b> be called Examples and not Example because of a perculiarity in the 
script that generates the manual.

By including the example input in a plumedfile environment you ensure two things:

- That the action names are converted to links to the relevant pages in the manual when the manual is constructed.
- That the code to construct the user manual will test to see if your example input can be parsed by PLUMED whenever the user manual is built.

To achieve the second of these objectives with the input shown above it is sufficient to include the example input in a plumedfile environment.
As detailed in the following sections, however, there are some cases where things are a little more complicated.

\subsection multirepeg Including example inputs for multiple replica simulations

If you have an input for a simulation that is to be run with three replicas such as the one below:

<pre class="fragment">
<span style="color:blue"># Compute a distance</span>
d: <a href="../../user-doc/html/_d_i_s_t_a_n_c_e.html" style="color:green">DISTANCE</a> ATOMS=1,2
<span style="color:blue"># Apply a restraint.</span>
<a href="../../user-doc/html/_r_e_s_t_r_a_i_n_t.html" style="color:green">RESTRAINT</a> ARG=d AT=@replicas:1.0,1.1,1.2 KAPPA=1.0
</pre>

Then you must specify that the input is to be run on three replicas in the first (SETTINGS) line of the input file as shown below: 

\verbatim
\plumedfile{3}
#SETTINGS NREPLICAS=3
# Compute a distance
d: DISTANCE ATOMS=1,2
# Apply a restraint.
RESTRAINT ARG=d AT=@replicas:1.0,1.1,1.2 KAPPA=1.0
\ endplumedfile /*** But with no space between the \ and the endplumedfile
\endverbatim

Notice that there should not be a space between the hash sign at the start of this line and word settings. 

\subsection auxfileeg Including example inputs that require an auxiliary file

Suppose that you have an input such as the one below:

<pre class="fragment">
<a href="../../user-doc/html/_r_m_s_d.html" style="color:green">RMSD</a> REFERENCE=file.pdb TYPE=OPTIMAL
</pre>

As RMSD has been used here you are also required to provide an input file which in this case would be called file.pdb.  You can include 
this input in an auxfile environment as shown below:

\verbatim
\auxfile{file.pdb}
ATOM      1  CL  ALA     1      -3.171   0.295   2.045  1.00  1.00
ATOM      5  CLP ALA     1      -1.819  -0.143   1.679  1.00  1.00
ATOM      6  OL  ALA     1      -1.177  -0.889   2.401  1.00  1.00
ATOM      7  NL  ALA     1      -1.313   0.341   0.529  1.00  1.00
ATOM      8  HL  ALA     1      -1.845   0.961  -0.011  1.00  1.00
END
\ endauxfile /*** But with no space between the \ and the endauxfile
\endverbatim

Obviously, the file.pdb inside the curly braces in the top line here indicates that the auxiliary file to be constructed from this data should be named 
file.pdb.  Files input in this way can be given any name but:

- If two auxfiles are used on the same page they must be given different names (if they are on different pages it does not matter)
- auxfiles should not be named *.dat as the script that builds the user manual assumes that all *.dat files are plumed input files. 

\subsection incfileeg Using INCLUDE in your example input files

Suppose that you have split your input by using an INCLUDE file as shown below:

<pre class="fragment">
<a href="../../user-doc/html/_d_i_s_t_a_n_c_e.html" style="color:green">DISTANCE</a> ATOMS=1,2 LABEL=dist
<a href="../../user-doc/html/_i_n_c_l_u_d_e.html" style="color:green">INCLUDE</a> FILE=toBeIncluded.inc
</pre>

<pre class="fragment">
<span style="color:blue"># this is toBeIncluded.inc</span>
<a href="../../user-doc/html/_r_e_s_t_r_a_i_n_t.html" style="color:green">RESTRAINT</a> ARG=dist AT=2.0 KAPPA=1.0
</pre>

To include an input like this in the manul you would write the following:

\verbatim
\plumedfile
DISTANCE ATOMS=1,2 LABEL=dist
INCLUDE FILE=toBeIncluded.inc
\ endplumedfile    /*** But with no space between the \ and the endplumedfile

\plumedfile
#SETTINGS FILENAME=toBeIncluded.inc  
RESTRAINT ARG=dist AT=2.0 KAPPA=1.0
\ endplumedfile   /*** But with no space between the \ and the endplumedincludefile
\endverbatim

By including the FILENAME attribute on the SETTINGS line you can set the name of the plumed input file that is generated when the input is tested.
Also notice that if, as in the example above, the included file is not (by itself) a valid plumed input it CANNOT be called *.dat as the script that 
checks the input will complain.  

\subsection molfileeg Using MOLFILE in your example input files

If you use have used a \ref MOLINFO command in the example input that you specified as has been done here:

<pre class="fragment">
<a href="./_m_o_l_i_n_f_o.html" style="color:green">MOLINFO</a> STRUCTURE=helix.pdb
<a href="./_w_h_o_l_e_m_o_l_e_c_u_l_e_s.html" style="color:green">WHOLEMOLECULES</a> ENTITY0=1-100
alpha: <a href="./_a_l_p_h_a_r_m_s_d.html" style="color:green">ALPHARMSD</a> RESIDUES=all TYPE=OPTIMAL R_0=0.1
</pre> 

Then you must provide information on the location from whence PLUMED can the reference input so that the example checking script can copy the input
for the MOLINFO.   The above input would thus be included in the manual as shown below:

\verbatim
\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=helix.pdb
WHOLEMOLECULES ENTITY0=1-100
alpha: ALPHARMSD RESIDUES=all TYPE=OPTIMAL R_0=0.1
\ endplumedfile    /*** But with no space between the \ and the endplumedfile
\endverbatim

\subsection otherfiles Other actions requiring external files/folder

Other actions in plumed may require reading input files, examples include reading gromacs .ndx files in \ref GROUP, reading chemical shifts in \ref CS2BACKBONE, etc.
To make these example work correctly in the manual you can use the keywords AUXFILE and AUXFOLDER as in the following:

\verbatim
\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt77/peptide.pdb
MOLINFO MOLTYPE=protein STRUCTURE=peptide.pdb
WHOLEMOLECULES ENTITY0=1-111

# This allows us to select only non-hydrogen atoms
#SETTINGS AUXFILE=regtest/basic/rt77/index.ndx
protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H

# We extend the cutoff by 0.1 nm and update the neighbor list every 40 steps
solv: EEFSOLV ATOMS=protein-h

# Here we actually add our calculated energy back to the potential
bias: BIASVALUE ARG=solv

PRINT ARG=solv FILE=SOLV

#SETTINGS AUXFOLDER=regtest/isdb/rt-cs2backbone/data NREPLICAS=2
cs: CS2BACKBONE ATOMS=1-174 DATADIR=data/
encs: ENSEMBLE ARG=(cs\.hn-.*),(cs\.nh-.*)
stcs: STATS ARG=encs.* SQDEVSUM PARARG=(cs\.exphn-.*),(cs\.expnh-.*)
RESTRAINT ARG=stcs.sqdevsum AT=0 KAPPA=0 SLOPE=24

PRINT ARG=(cs\.hn-.*),(cs\.nh-.*) FILE=RESTRAINT STRIDE=100
\ endplumedfile    /*** But with no space between the \ and the endplumedfile
\endverbatim

\section tutorials Writing how-to instructions

On every page of the plumed user manaul there are three tabs: Main-page, Glossary and How-tos.  Here we are going to describe how to
go about writing something that will appear on the How-tos page.  On that page you will find that there are two kinds of resources.  
The first set of resources are a set of links to shortist descriptions of how to particular things with plumed.  The second set of 
resources are links to a set of external websites that we think might be of interest to users of plumed.

\subsection websites Adding a link to your website 

If you have a website that you think it would be useful for us to link from the plumed How-tos page please send us a file that contains
the following information in the following format:

\verbatim
link: http://en.wikipedia.org/wiki/Metadynamics

description: A wikipedia article on metadynamics
\endverbatim 

In case it isn't abundantly clear the line that starts "link:" contains the hyperlink to your page, while the second line (the one starting "description:")
contains the description of the page that will appear beside the link on the how-tos page.  If this file is placed in the user-doc/tutorials directory of the
plumed source and if it is given a name that ends in .site then the link will appear on the how-tos page.

\subsection tute Writing a how-to

Lets say you now want to write a set of how-to instructions.  You will first need to create an additional file in the user-doc/tutorials directory of the plumed
source.  Much like the rest of plumed manual this file is written in Doxygen so you should read the instructions about how to go about using \ref Equations, \ref Lists, 
\ref Formatting, \ref Citations and \ref Examples.  You will also need to ensure that the text that you want to appear on the your page is contained between
the special doxygen comment symbols as shown below:

\verbatim
This will not appear on the final how-to page.

/**
This will appear on the final how-to page.
*/

This will also not appear on the final how-to page.
\endverbatim 

To ensure that a link is created on the main How-tos page you need to include the following instructions after the closing doxygen comment as shown below:

\verbatim
/**
\page filename My explanation of something in plumed

Text of how-to page

*/
link: @subpage filename

description: The description of what is described on your page
\endverbatim

For this particular example the above should be contained in a file inside the user-doc/tutorials directory called filename.txt.  If the user comes to the 
How-tos page they will see a link with the text "My explanation of something in plumed," which they can click on to get to your page.  Beside that link the 
description "The description of what is described on your page" will appear.

One final thing, if your How-to is a tutorial - that is to say if you have a set of exercises for users to work through - then it may be useful to 
provide them with some example files that they can download.  This is possible you add the files to our repository and you can put a link to download them
somewhere on your page.
To keep things manageable there should only be one (or a few) tar-ball file per tutorial.
As such if you need to provide users with multiple files please put them in a single directory. The script that 
build the documentation will make one or more tar-balls from there.
please put them in a zipped tar ball.
The simplest way is thus to add a directory
in the user-doc/tutorials directory and ideally given a name so that it can be identified
with your corresponding *.txt file.
For example if the name of the tutorial is `mytute` you should make a directory named
`user-doc/tutorials/mytute`. 
You can then download it by including the following into your *.txt file.

\verbatim
/**
\page mytute A tutorial that helps users learn to do something with plumed

Please download the following <a href="tutorial-resources/mytute.tar.gz" download="mytute.tar.gz"> file </a>.

*/

link: @subpage mytute

description: A tutorial for users to work through

additional-files: mytute
\endverbatim 

In this case the tar ball you add is called mytute.tar.gz.  The user can download this file by clicking on the word file.

\section updating-web-manuals Updating web manuals

Precompiled versions of PLUMED manuals can be found on github at an address such as http://www.plumed.org/doc-v2.1/user-doc/html/index.html
(replace v2.1 with the actual version number). These manuals take advantage of a nice github feature: any branch named gh-pages
is shown as a webpage. In this example, the repository is located at http://github.com/plumed/doc-v2.1 .
Before version 2.1.1 it was necessary to upload the precompiled manual by hand. Since version 2.1.1, this is done
from Travis CI automatically whenever a commit containing in its log the tag [makedoc] is pushed into the plumed2 github repository.
Since version 2.3.3, manual is always updated, and tag [makedoc] is ignored.
Notice that Travis CI will try to push the manual on a repository named http://github.com/plumed/doc-NAMEOFTHEBRANCH , so that 
this should work for all the future release branches as long as an appropriate repository is created on the github.com/plumed
organization.
We could even easily create repositories to host the documentation of temporary branches.
Also notice that these repositories (plumed/doc-xxx) need to give write access to a dummy github account (PlumedBot). A token
for that user enabling html access is stored in the environment variable GIT_TOKEN which is saved (not visible) on travis-ci.org.
In this way, any commit made in the plumed repository by one of the developers will have access to the variable and will trigger
manual build and push. Conversely, pull requests by external users should not be able to
access the token and won't update manual changes.
Starting with PLUMED 2.3 the web manual also contains a coverage scan that shows which parts of the code are actually
covered by the regtests. As of PLUMED 2.5, the coverage scan ends up on a separate repository named 
http://github.com/plumed/coverage-NAMEOFTHEBRANCH.

Notice that to solve [this issue](https://github.com/plumed/plumed2/issues/239) as of PLUMED 2.3.2 the script that
pushes the documentation to travis-ci adds special information to remove from search engine results pages from
unofficial or unsupported branch (see .ci/push script).

Bottom line: manual will always be updated after a commit that can pass the tests.
Twenty minutes or so after your push the manual should be up to date, remember to double check on the web
and to revert the commit if there are errors!

It is possible to generate PLUMED manuals for your own personal forks 
using a similar procedure as described above. 
For this to work you need to enable Travis CI for your forked repository 
and define appropriately the environment variables on Travis CI. 
The github account used to automatically push the generated manuals 
should be defined using the `GIT_BOT` variable, 
preferably this should be a dummy account. A github token
enabling html access for that account should be defined using the `GIT_TOKEN` variable. 
Furthermore, you need to define an email address associated to the account using the `GIT_BOT_EMAIL` variable. 
It is better to make all these environment variable hidden such that they are 
not shown in the public logs on travis-ci.org. 
To generate a manual for a specific branch you need to create a repository 
`USERNAME/doc-NAMEOFTHEBRANCH` and give write access to the account given in 
`GIT_BOT`. The generated manuals will be accessible on 
https://USERNAME.github.io/doc-NAMEOFTHEBRANCH. Note that manuals generated in 
this way will always be labeled as unofficial and not shown in search engine results.
Starting with PLUMED 2.5, if you want to show the results of the coverage scan you should
similarly create arepository `USERNAME/coverage-NAMEOFTHEBRANCH`.


\page HowToPlumedYourMD How to add plumed to an MD code

\brief Learn how to use plumed in a not yet supported MD code

Plumed ships with scripts that can be used to add it to many of the standard MD packages.  Obviously though, if no patch is provided for the MD code you use then you will have to write one yourself.  Plumed has been designed so that it can be added to an MD code
either statically (i.e. as a collection of objects) or as a dynamic library.
For technical reasons it
is NOT possible to pack all the plumed objects as a static library and link the library
(if you really want to know it: when static librearies are linked the linker
typically discards all the objects which are not used; plumed is using an
automatic registration process which is not compatible with this choice).
Furthermore, unlike the previous version of plumed, the plumed source code is not compiled at the same time as your MD code is compiled.  Instead plumed now has its own makefile and is compiled separately to the MD engines that rely on it.  This makes plumed patches considerably simpler as now they only do two things:

- Modify the makefile so that the plumed is linked to the MD code
- Modify the source code to add all the required calls to plumed

\section makefile Modifying your makefile

Once the plumed source code has been compiled a file src/Plumed.inc is
generated. This file should be included in your MD code's makefile as it
informs the code where to find all the plumed source. There are three possible
ways to link PLUMED: static (linking the .o files directly), shared (linking
the libplumed.so library) or runtime (which links only a wrapper to plumed,
whereas the location of the remaining part on plumed, i.e. the libplumedKernel.so,
can be specified at runtime). The relevant variables are
- \$(PLUMED_STATIC_LOAD) the options to be added to the linker for static
  linking, i.e. all the plumed objects with their full path,
   plus special linking options, plus al the libraries used by plumed (e.g.
  matheval)
- \$(PLUMED_SHARED_LOAD) the options to be added to the linker for shared
  linking, i.e. full plumed shared library (libplumed.so) with its full path
  plus special linking options.
- \$(PLUMED_RUNTIME_LOAD) the options to be added to the linker for runtime
  linking, i.e. the wrapper object with its full path plus special linking options.

The libplumedKernel.so basically contains the same objects as the
libplumed.so except for the wrapper.

The simplest approach is to take advantage of the three variants
src/Plumed.inc.static , src/Plumed.inc.shared and src/Plumed.inc.runtime .
Including one of those, it is sufficient to add to the linker line the macro
\$(PLUMED_LOAD).

The best way to patch your MD file is:
- go into the root directory of the MD code and type
\verbatim
> plumed patch --new name-of-the-code
\endverbatim
  this will create a new (empty) patch file for your MD code
- patch the code with the empty patch
\verbatim
> plumed patch --patch --engine name-of-the-code
\endverbatim
  this will produce a symbolic link Plumed.inc in the root directory of the MD
  code
- make a backup copy of the Makefile of your MD code, naming it
  Makefile.preplumed
\verbatim 
> cp src/Makefile src/Makefile.preplumed
\endverbatim
- edit the Makefile including the Plumed.inc file
  (e.g. "include ../Plumed.inc") and adding "\$(PLUMED_LOAD)" to the linking
command
- save your modification into the plumed patch
\verbatim
> plumed patch --save
\endverbatim

Now you can patch/unpatch your MD code, automatically linking PLUMED, as if it
was an officially released interface (even though at this stage you do not
have any functionality of plumed yet in your code - see below how to add it).

Note that after you have this done, PLUMED creates a file in
$PLUMED_ROOT/patches/ that is generally named 
name-of-the-code.diff ( if your code is named name-of-the-code in the "save"
phase). You can also manually create, in the same directory,  a file  
called name-of-the-code.config so that you can use to perform five actions (totally optional).
The first one can be used to do some heuristic check to verify that the patch
is for the correct MD engine and that the MD source tree is in the correct state
(this typically means that it hass been already configured).
Then there are two functions that can
do some modifications to some
files that might vary considerable (as the Makefiles produced by the config
procedure that change with the machine) and do some actions after the
patching procedure is performed (a sanity check or a makedependence script to
be rerun).
The last two ones should revert the effect of the previous ones in such
a way that it is possible to smoothly revert the patch. Try to keep
them consistent so as the revert feature works correctly.
This file is in bash and typically looks like this 
\verbatim
function plumed_preliminary_test(){
# Just check of the code is the right one: no error means it is ok 
  grep -q  name-of-the-code Somefile 
}

function plumed_before_patch(){
# typically do something on the Makefile
}

function plumed_after_patch(){
# typically call some makedepends if necessary
}

function plumed_after_revert(){
# revert exacty the effect of plumed_before_patch
# typically undoes what was done on the Makefile
# additionally, one could add a makedepends
}

function plumed_before_revert(){
# revert exactly the effect of plumed_after_patch
# however, a makedepends script should always go in plumed_after_revert
}

\endverbatim

\attention Be careful with these scripts. You have access to a lot of flexibility
(adding files, modifying existing ones arbitrarily etc) but you should always 
take care of the revertibility with "patch -r".




\section language Language dependence

The next step is to add proper calls to plumed into the MD source code.

The interface between plumed and the MD codes is written is plain C so it is compatible with codes that are written in plain C, C++ and fortran.
To use this interface one should first create a plumed object (plumed_create()), then sent to it
messages (using plumed_cmd()) and deallocate it (plumed_finalize()) at the end
of the simulation. All of these routines have C++ and FORTRAN equivalents.
Notice that in C the plumed object is represented using a struct (plumed), in
C++ it is represented using a class (PLMD::Plumed) and in FORTRAN using a
string of 32 characters (CHARACTER(LEN=32)).
Also notice that whenever passing strings from FORTRAN you should add to them
"char(0)" as
a string terminator to be sure that the string is properly interpreted
(e.g. "pass-this-string"//char(0)).

As people might get confused by
the idea of "plumed being an object", we also added the possibility of using
a "generic global plumed instance", which can be accessed
using routines with a _g_ in the name (plumed_g_create(),
plumed_g_cmd() and plumed_g_finalize()).

This interface is very general, will not change in future releases, and is
fully contained in Plumed.h and Plumed.c files. For a reference to it, see
\ref ReferencePlumedH


\section mdImplementation A practical example

We now describe how to pass data to plumed in from an MD code that is written in C/C++.
First save all the files that you intend to modify into a .preplumed file (for
example if you want to modify myfile.c save it first in its original version
into myfile.c.preplumed): this will tell the patching procedure that this  is
a file that will enter the set of the files to be modified. 

In C or C++ files containing calls to plumed you will have to include the Plumed.h file.
In addition, you will also have to define a plumed obect that is visible in all these routines
and you will probably like to define some sort of plumedswitch,
which can be read in from the input to your code and used to tell the code whether or not this is to be a run with plumed.
Finally, you might like to include something in input so that you specify the name of the plumed input file.
How these things are best done will depend on your code and so we leave it to your discretion.  

Plumed must perform three actions inside your MD code:

- It must be initialized before the main MD loop starts so that the plumed input files are read in.
- It must be called every MD step so that the forces from the bias can be computed.
- It must be finalized when the simulation is completed

The only routine which can be used to send message to plumed is plumed_cmd
(or, equivalently, Plumed::cmd in C++ and plumed_f_cmd in FORTRAN).

Notice that as we PLUMED evolves new commands could be available.
Thus, if you want your interface to be compatible with multiple PLUMED versions, you should first
check that the PLUMED version is new enough to accept a given command.
To know which is the
currently linked PLUMED version see \ref apiversion .

The various calls that can be used during initialization are as follows:

\verbatim
plumed plumedmain; plumedmain=plumed_create();                 // Create the plumed object

// Calls to pass data to plumed
plumed_cmd(plumedmain,"setRealPrecision",&real_precision);     // Pass a pointer to an integer containing the size of a real number (4 or 8)
plumed_cmd(plumedmain,"setMDEnergyUnits",&energyUnits);        // Pass a pointer to the conversion factor between the energy unit used in your code and kJ mol-1
plumed_cmd(plumedmain,"setMDLengthUnits",&lengthUnits);        // Pass a pointer to the conversion factor between the length unit used in your code and nm 
plumed_cmd(plumedmain,"setMDTimeUnits",&timeUnits);            // Pass a pointer to the conversion factor between the time unit used in your code and ps

// This is valid only if API VERSION > 3
plumed_cmd(plumedmain,"setMDChargeUnits",&chargeUnits);        // Pass a pointer to the conversion factor between the charge unit used in your code and e

// This is valid only if API VERSION > 3
plumed_cmd(plumedmain,"setMDMassUnits",&massUnits);            // Pass a pointer to the conversion factor between the mass unit used in your code and amu

plumed_cmd(plumedmain,"setPlumedDat",&plumedInput);            // Pass the name of the plumed input file from the md code to plumed
plumed_cmd(plumedmain,"setMPIComm",&MPI_COMM_WORLD);           // Pass a pointer to the MPI communicator to plumed
// notice that from fortran the command "setMPIFComm" should be used instead
plumed_cmd(plumedmain,"setNatoms",&natoms);                    // Pass a pointer to the number of atoms in the system to plumed
plumed_cmd(plumedmain,"setMDEngine","gromacs");                // Pass the name of your md engine to plumed (now it is just a label) 
plumed_cmd(plumedmain,"setLog",fplog);                         // Pass the file on which to write out the plumed log (if the file is already open)
plumed_cmd(plumedmain,"setLogFile",fplog);		       // Pass the file  on which to write out the plumed log (to be created)
plumed_cmd(plumedmain,"setTimestep",&delta_t);                 // Pass a pointer to the molecular dynamics timestep to plumed

// This is valid only if API VERSION > 1
plumed_cmd(plumedmain,"setKbT",&kbT);                          // Pointer to a real containing the value of kbT

// This is valid only if API VERSION > 2
plumed_cmd(plumedmain,"setRestart",&res);                      // Pointer to an integer saying if we are restarting (zero means no, one means yes)

// This is valid only if API VERSION > 3
plumed_cmd(plumedmain,"readInputLine","d: DISTANCE ATOMS=1,2");// Read a single input line directly from a string

// This is valid only if API VERSION > 7
plumed_cmd(plumedmain,"readInputLines","d: DISTANCE ATOMS=1,2\n"
                                       "PRINT ARG=d");         // Read a multiple lines directly from a string. Allows comments and continuation lines.

// Calls to do the actual initialization (all the above commands must appear before this call)
plumed_cmd(plumedmain,"init",NULL);                            // Do all the initialization of plumed
plumed_cmd(plumedmain,"read",read);                            // Read the plumed input.  N.B. This is called during init and so this call is only required in special cases. 
\endverbatim

Please note that if your code is in FORTRAN you should append a "char(0)"
token to every string. Also, remember that FORTRAN is by default passing
arguments by reference, so that the "&" symbols which are required in C are
not necessary in FORTRAN.


The various calls that can be used pass data and calculate the forces due to the bias are as follows:

\verbatim
// Calls to pass data to plumed
plumed_cmd(plumedmain,"setStep",&step);                      // Pass a pointer to the current timestep to plumed
/ *** The way that you pass positions will depend on how they are stored in your code.  If the x, y and z position are all stored in a single array you may use:
plumed_cmd(plumedmain,"setPositions",&pos[0][0]);            // Pass a pointer to the first element in the atomic positions array to plumed  
                                                             // assuming they
                                                             // are stored in
                                                             // a
                                                             // x1,y1,z1,x2,y2,z2 ...
                                                             // kind of ordering
/ *** Othersize if you pass the three separate vectors of x, y and z positions using:
plumed_cmd(plumedmain,"setPositionX",&x[0]);                 // Pass a pointer to the first element in the array of x component of the atomic positions to plumed
plumed_cmd(plumedmain,"setPositionY",&y[0]);                 // Pass a pointer to the first element in the array of y component of the atomic positions to plumed
plumed_cmd(plumedmain,"setPositionZ",&z[0]);                 // Pass a pointer to the first element in the array of z component of the atomic positions to plumed
plumed_cmd(plumedmain,"setMasses",&mass[0]);                 // Pass a pointer to the first element in the masses array to plumed
plumed_cmd(plumedmain,"setCharges",&charge[0]);              // Pass a pointer to the first element in the charges array to plumed
plumed_cmd(plumedmain,"setBox",&box[0][0]);                  // Pass a pointer to the first element in the box share array to plumed
plumed_cmd(plumedmain,"setEnergy",&poteng);                  // Pass a pointer to the current value of the potential energy to plumed?
/ *** The way that you pass forces will depend on how they are stored in your code.  If the x, y and z force are all stored in a single array you may use:
plumed_cmd(plumedmain,"setForces",&f[0][0]);                 // Pass a pointer to the first element in the foces array to plumed
/ *** Othersize if you pass the three separate vectors of x, y and z forces using:
plumed_cmd(plumedmain,"setForcesX",&fx[0]);                  // Pass a pointer to the first element in the array of the x components of the atomic forces to plumed
plumed_cmd(plumedmain,"setForcesY",&fy[0]);                  // Pass a pointer to the first element in the array of the y components of the atomic forces to plumed
plumed_cmd(plumedmain,"setForcesZ",&fz[0]);                  // Pass a pointer to the first element in the array of the z components of the atomic forces to plumed
plumed_cmd(plumedmain,"setVirial",&force_vir[0][0]);         // Pass a pointer to the first element in the virial array to plumed

// Calls to do actual calculations
plumed_cmd(plumedmain,"calc",NULL);                          // Calculate and apply forces from the biases defined in the plumed input

// One can break up the "calc" command in two parts:
plumed_cmd(plumedmain,"prepareCalc",NULL);                   // Prepare to do a calculation by requesting all the atomic positions from the MD code
plumed_cmd(plumedmain,"performCalc",NULL);                   // Use the atomic positions collected during prepareCalc phase to calculate colvars and biases.
// The "performCalc" step can be further split into:
// This is valid only if API VERSION > 3
plumed_cmd(plumedmain,"performCalcNoUpdate",NULL);           // Same as "performCalc", skipping the update phase. Could be called multiple time per step
// This is valid only if API VERSION > 3
plumed_cmd(plumedmain,"update",NULL);                        // Only performs the update phase. Should be called once per step

// After the first part it will be possible to ask PLUMED e.g. if the energy is required with
plumed_cmd(plumedmain,"isEnergyNeeded,&flag);                // assuming flag is an int, that will be set to 0 if energy is not needed and 1 if it is needed

// The "prepareCalc" can be further split into:
plumed_cmd(plumedmain,"prepareDependencies",NULL);           // Work out what we are calculating during this MD step (this is the first step of prepareCalc)
plumed_cmd(plumedmain,"shareData",NULL);                     // Request all the atomic positions from the MD code (this is the second step of prepareCalc)
// This will allow to overlap sharing of atoms across multiple processors (sent in shareData) and calculation

// Some extra calls that might come in handy
plumed_cmd(plumedmain,"createFullList",&n);                  // Create a list containing of all the atoms plumed is using to do calculations (return the number of atoms in n)
plumed_cmd(plumedmain,"getFullList",&list);                  // Return a list (in list) containing all the indices plumed is using to do calculations. list should be `const int*`
plumed_cmd(plumedmain,"clearFullList",NULL);                 // Clear the list of all the atoms that plumed is using to do calculations
plumed_cmd(plumedmain,"clear",clear);                        // Clear and delete all the pointers inside plumed.
\endverbatim

The plumed calls for the finalization tasks is as follows:

\verbatim
plumed_finalize(plumedmain);          // Call the plumed destructor
\endverbatim

\section mpicodes Dealing with parallelism

Plumed has functionality to deal with parallel MD codes.  The particular form of the functionality used to do this (and the frequency with which you will have to call these routines) will depend on whether your code is parallelism using a domain decomposition or particle decomposition strategy.  The calls for required for using this functionality are as follows:

\verbatim
plumed_cmd(plumedmain,"setAtomsNlocal",&nlocal);            // Pass a pointer to the number of atoms on this node
plumed_cmd(plumedmain,"setAtomsGatindex",gatindex);         // Pass an array (from a code in c) containing the indices of all the atoms on this node (used for domain decomposition)
plumed_cmd(plumedmain,"setAtomsFGatindex",gatindex);        // Pass an array (from a code in fortran) containing the indices of all the atoms on this node (used for domain decomposition)
plumed_cmd(plumedmain,"setAtomsContiguous",&start);         // Number the atoms on this node from start to start+nlocal   (used for particle decomposition)
\endverbatim

\section apiversion Inquiring for the plumed version

New functionalities might be added in the future to plumed. The description of
all the possible "commands" sent to plumed is called its API (application
programming interface). If you want to know which API is supported by plumed
you can use this call:
\verbatim
plumed_cmd(plumedmain,"getApiVersion",&api);                 // Pass the api version that plumed is using
\endverbatim 
With the current version, this will set the api variable (an integer) to 2. As
we add new features, this number will be increased.

\section Saving the diffs

This is similar to plumed 1. All the files that you want to modify should be
first copied to .preplumed files. Then use "plumed patch -s" to save the diff.

\section debugging-patch Debugging your PLUMED+MD code

If you want to be sure that your MD code is going to work properly with PLUMED we strongly suggest to go through the checks below.

\subsection debugging-patch-positions Check positions

The first thing to do is to checkout if positions are properly passed to PLUMED.
Let's assume your system has 100 atoms.
Just run PLUMED with an input such as the following one
\verbatim
DUMPATOMS ATOMS=1-100 FILE=test.xyz
\endverbatim
File test.xyz should look like:
\verbatim
100
 1000.000000 1000.000000 1000.000000
X 0.2 0.9 1.2
X 0.7 2.9 3.6
....... other 98 lines ......
100
 1000.000000 1000.000000 1000.000000
X 0.3 0.8 1.1
X 0.8 2.8 3.5
....... and so on ......
\endverbatim
That is, foreach snapsnot: first line, number of atoms; second line, box; then for each atom name/x/y/z coordinates.
Notice that coordinates are expected to be in in nanometers (PLUMED units).
For orthorhombic boxes the box line will contain three numbers, namely box size in x, y, and z directions.
For non orthorombic boxes, the box line will contain nine numbers, namely (ax,ay,az,bx,by,bz,cx,cy,cz),
where ax is the x component of the first lattice vector.

The produced test.xyz should match the trajectory file produced by your MD code.
Check all the atoms (not just the first ones).
Also, check that results are consistent also when you run in parallel (with particle or domain decomposition).

\subsection debugging-timestep Check timestep

PLUMED relies on the timestep from the MD code for several calculations. To be sure it has been passed properly,
run PLUMED with
\verbatim
d: DISTANCE ATOMS=1,10
PRINT ARG=d FILE=colvar
\endverbatim
Check the first column of `colvar` file (time). It should increase by one timestep at every line. Notice
that time here should appear in picoseconds (PLUMED units).

\subsection debugging-pass-energy Check energy

(This only applies if your MD codes actually passes energy to PLUMED)

Just try the following
\verbatim
e: ENERGY
PRINT ARG=e FILE=colvar
\endverbatim
Check the second column of `colvar` file (energy). It should be equivalent to the energy computed
in the MD code. Notice that PLUMED will write energy in kj/mol.

\subsection debugging-mass-charges Check masses and charges

The best way to debug the masses and charges is to print their values from PLUMED and
compare it with the MD code. This will be possible with PLUMED 2.2.

Meanwhile, you can try to compute center of masses or dipole moments and compare with what you expect.

\subsection debugging-patch-forces Check passed forces

Then you should double check whether PLUMED is able to return proper forces to the MD engine.
A first check could be the following. Run a short MD with the following PLUMED input
\verbatim
d: DISTANCE ATOMS=1,10
PRINT ARG=d FILE=colvar
\endverbatim
The run again with
\verbatim
d: DISTANCE ATOMS=1,10
RESTRAINT ARG=d AT=0.0 SLOPE=10
PRINT ARG=d FILE=colvar-biased
\endverbatim
Now plot both files with
\verbatim
gnuplot> plot "colvar" u 1:2 , "colvar-bias" u 1:2
\endverbatim

The two lines should start at the same value, but the second one should be systematically below the first one.

Notice that this test is not quantitative. Since an error in the force units would just give a qualitative
change, it is better to do a more rigorous check.

The best way to do a quantitative check is to alternatively add a restraint with the MD code and
with PLUMED and checkout that the obtained results are equivalent.

\todo better explain here

If your code passes these tests, you can likely start to do biased MD simulations in the NVT ensemble.
If you need NPT ensemble or if you want to bias the total energy you should continue with further tests.

\subsection debugging-patch-virial Check virial contribution

Most of the people use plumed to bias a small number of coordinates. This makes the contribution of plumed forces to the total virial
very small and makes it particularly difficult to find errors in this part. In case you combine plumed with a new MD code we highly
suggest to use the following procedure to debug the virial.

First run a short simulation (1000 steps should be enough) at constant pressure with pressure=1bar and the following plumed input:
\verbatim
v: VOLUME
PRINT ARG=v FILE=volume
\endverbatim

Then run another short simulation starting from identical conditions (i.e. same positions, velocities, random seed, etc) at constant pressure
with pressure=1001bar and the following plumed input:
\verbatim
v: VOLUME 
# slope should be just 10 times the Avogadro constant:
RESTRAINT AT=0.0 ARG=v SLOPE=-60.2214129
PRINT ARG=v FILE=volume2
\endverbatim
In this way, the negative pressure provided by plumed should exactly compensate for the extra 1000bar set in the barostat.
Thus, the two files `volume` and `volume2` should be virtually identical. Notice that small differences
due to different precision in the storage of avogadro number and other issues related to difficulties in exactly reproducing
a simulation could make the two trajectory slightly different.

If you pass this test, you can safely run biased MD in the NPT ensemble. Otherwise, there could be some
issue in the way the virial is passed to the MD code.

\subsection debugging-patch-energy Check forces on energy

(This only applies if your MD codes actually passes energy to PLUMED)

First run a short simulation (1000 steps should be enough) at constant temperature
\f$T=300\f$K and a given integration timestep \f$\Delta t\f$.
and the following PLUMED input:
\verbatim
e: ENERGY
PRINT ARG=e FILE=energy1
\endverbatim

Then run another short simulation starting from identical conditions (i.e. same positions, velocities, random seed, etc) at constant temperature
\f$T'=\alpha^2 T\f$K where  \f$\alpha=1.1\f$ (that is \f$T=363\f$) and a shorter timestep \f$\Delta t'=\frac{\Delta t}{\alpha}\f$.
Notice that all absolute times provided in the MD input (e.g. relaxation time for the thermostat) should be consistently divided by \f$\alpha\f$.
Use the following PLUMED input:
\verbatim
e: ENERGY
# slope is such that 
PRINT ARG=e FILE=energy2
# slope should be (alpha-1)=0.1
RESTRAINT AT=0.0 ARG=e SLOPE=0.1
\endverbatim

The two files `energy1` and `energy2` should be virtually identical.

In case you were able to have the virial properly working (see previous section), then you can try the same with a constant temperarure-constant pressure
simulation. In this case, please also monitor the volume of the resulting trajectory.
\page CodeFormatting How to format code properly

Since version 2.3.2, we format code using <a href="http://astyle.sourceforge.net/"> astyle </a>.
As a convention, we use `astyle` version 3.00, with the options that are
reported in the file `.astyle.options` located in the root directory of PLUMED.
You might want to automatize the application of `astyle` using those options.
As of now, you can use the following command from root PLUMED directory:
\verbatim
> make astyle
> git commit
\endverbatim
Notice that this command will both apply `astyle` to all the relevant files as well as
add them to the next git commit. After having inspected the changes, you can commit them.
Also notice that running this command from PLUMED root directory will also
compile the `astyle` version that is distributed with PLUMED. We decided to distribute
`astyle` within the PLUMED repository to make sure that everyone is using exactly the same version.
In addition, you can run `make astyle` directly within a module directory
so as to only reformat that specific module.

Additional care must be used while merging branches. In this case, you should
make sure that both branches are formatted with `astyle` before merging them.
The procedure discussed below should be made once for each not-yet-formatted branch that you are maintaining.

Let's say that you are working on a branch `feature` that is not yet formatted, and you
would like to merge changes from branch `v2.3` that is already formatted.
You should use the following list of commands
\verbatim
# Bring in all the changes in v2.3 up to astyle formatting (excluded):
> git merge astyle-v2.3~1
# Notice that this will include the astyle scripts, but not
# the big commit formatting the whole code.

# Mark the following commit (astyle-v2.3) as merged, even though
# it is completely ignored:
> git merge -s ours astyle-v2.3
# This is necessary since this commit is too big to be really merged.

# Indeed, instead of merging it, we apply directly astyle to the
# current branch:
> make astyle

# Now the two branches are both formatted and can be merged one into the other.
# Merge all the newer commits on v2.3
> git merge v2.3
\endverbatim

Notice that here `astyle-v2.3` is a tag that refers to the commit where we introduced
formatting in v2.3 branch. 
In a similar way you can bring in any changes from another branch, just replace
`v2.3` with the proper name (e.g. `master`).
After merging, your branch will be also formatted correctly.
Notice that you cannot work it in the opposite direction (that is, unformat an already
formatted branch). Finally, consider that rebasing might be more complicated.
When all branches will be formatted this will not be an issue anymore.

\page intro-git A brief introduction to git

Clone the github repository:
\verbatim
> git clone git@github.com:plumed/plumed2.git
> cd plumed2
\endverbatim

Stay up to date:
\verbatim
> git pull
\endverbatim

Make a small fix (working locally):
\verbatim
> git add FILENAME
> git commit FILENAME
\endverbatim

Share it (needs internet connection):
\verbatim
> git pull # always check if you are up-to-date
> git push
\endverbatim

Look at what's happening:
\verbatim
> git log
\endverbatim
or better
\verbatim
> gitk --all
\endverbatim

Start working on a new feature, opening a new branch:
\verbatim
> git checkout -b new-feature
\endverbatim

List the available branches:
\verbatim
> git branch
\endverbatim

Check the available branches including the ones on the origin
\verbatim
> git branch -a
\endverbatim

Switch among branches
\verbatim
> git checkout master
> git checkout new-feature
\endverbatim

Do a commit on your new-feature branch
\verbatim
> git checkout new-feature
> # now edit NEWFILE
> git add NEWFILE
> git commit NEWFILE
\endverbatim

Merge recent work from the master branch, doing a "rebase"
\verbatim
> git checkout master
> git pull # to stay up-to-date with remote work
> git checkout new-feature
> git rebase master
\endverbatim
Notice that rebasing is only recommended if your new-feature branch
has not been shared yet with other people.

Collect the changes to the log and get rid of branches that have 
been deleted on the remote repo:
\verbatim
> git fectch --all --prune
\endverbatim

After several commits, your new feature is ready for merge.
\verbatim
# checkout the branch you want to merge your work on (e.g. master)
> git checkout master
> git pull # to stay up-to-date with remote work
> git checkout new-feature
# You can squeeze your commits with:
> git rebase -i master
# then interactively picking your commits (follow onscreen instructions)
# (Without -i, all the commits are retained)
# Merge your work into master branch
> git checkout master
> git merge new-feature
# Remove the branch
> git branch -d new-feature
# In case you want to remove a remote branch you should use:emove a remote branch
# > git push origin :crap-branch
# Analyze the results
> gitk --all
# If everything seems right, push your master branch to the central repository
> git push
\endverbatim

Checkout a remote branch and switch to it on your local repo 
\verbatim
> git checkout -b experimental origin/experimental
\endverbatim

Compare list of commits on two different branches (here branches are named master and experimental)
\verbatim
git log master..experimental
\endverbatim

All these things can be done with a GUI:
\verbatim
> git gui
\endverbatim

\page AddingAColvar How to add a new collective variable

To implement a CV you need to create a single cpp file called <i>ColvarName</i>.cpp in the directory src/colvar. 
In addition, if you would like us to incorporate your CV in the release version of PLUMED you will need to write at least
one regression test for your CV.  Details on how to write regression tests are provided here: \ref regtests

If you use the following template for your new <i>ColvarName</i>.cpp file then the manual and the calls to the CV will 
be looked after automatically.

\verbatim
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR NAME
/*
\endverbatim

At this point you provide the description of your CV that will appear in the manual along with an description of the input file syntax and an example.  Merging new features of the code into the plumed main branch without proper documentation is punishable by death!  Some instructions as to how to format this information is provided here: \ref usingDoxygen

\verbatim
*/
//+ENDPLUMEDOC

/**** We begin by declaring a class for your colvar.  This class inherits everything from the Colvar class.
      This ensures it has a label, a place to store its value, places to the store the values of the derivatives
      and that it can access the various atoms it will employ.

class ColvarNAME : public Colvar {
\endverbatim

Insert declarations for your colvar's parameters here using plumed's \ref parsing.

\verbatim
public:
 /---- This routine is used to create the descriptions of all the keywords used by your CV
 static void registerKeywords( Keywords& keys ); 
 /---- This is the constructor for your colvar.  It is this routine that will do all the reading.
       Hence it takes as input a line from the input file.
  ColvarNAME(const ActionOptions&);
 /---- This is the routine that will be used to calculate the value of the colvar, whenever its calculation is required.
       This routine and the constructor above must be present - if either of them are not the code will not compile.
  virtual void calculate();
};

 /------ The following command inserts your new colvar into plumed by inserting calls to your new
        routines into the parts of plumed where they are required.  This macro takes two arguments:
        The first is the name of your ColvarClass and the second is the keyword for your CV
        (the first word in the input line for your CV).
PLUMED_REGISTER_ACTION(ColvarNAME,"KEYWORD")

 /----- The following routine creates the documentation for the keyowrds used by your CV
void ColvarName::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
\endverbatim

In here you should add all your descriptions of the keywords used by your colvar as well as descriptions of any components
that you can use this colvar to calculate. Descriptions as to how to do this can be found here: \ref usingDoxygen

\verbatim
}

 /---- We now write the actual readin (constructor) and calculations routines for the colvar

ColvarName::ColvarName(const ActionOptions&ao):
 /------ This line sets up various things in the plumed core which colvars rely on.
PLUMED_COLVAR_INIT(ao)
{
 vector<int> atoms;  /----- You almost always have atoms -----/
\endverbatim

Insert code here to read the arguments of the CV here using plumed's parsing functionality.  N.B.  The label is read in already elsewhere.

\verbatim
  checkRead();     /--- This command checks that everything on the input line has been read properly

/--- The following two lines inform the plumed core that we require space to store the value
     of the CV and that the CV will act on a particular list of atoms.
  addValueWithDerivatives("");
  requestAtoms(atoms);

/ --- For a number of the free energy methods in plumed it is necessary to calculate the
      distance between two points in CV space.  Obviously, for periodic CVs one must take
      periodicities into account when calculating distances and use the minimum image
      convention in distance calculations.  Hence, we set the periodicity of the cv using
      the following two lines.
   getValue("")->setPeridodicity(true);  // Set this true if the CV is periodic otherwise set if false.
   getValue("")->setDomain(min,max);     // This routine is only required if the function is periodic.  It sets the minimum and maximum values of the colvar.
}

void ColvarName::calculate(){
/--- These are the things you must calculate for any cv ---/
  double cv_val;              /--- The value of the cv ----/
  Tensor boxDerivatives;      /--- The derivative of the cv with respect to the box vectors ----/
  vector<double> derivatives; /--- The derivative of the cv with respect to the atom positions ---/
\endverbatim

Insert the code to calculate your cv, its derivatives and its contribution to the virial here. Please use, where possible, the library of tools described in \ref TOOLBOX.

\verbatim
/---- Having calculated the cv, its derivative and the contribution to the virial you now
      transfer this information to the plumed core using the following three commands. 
  for(int i=0;i<derivatives.size();i++){ setAtomsDerivatives(i,derivatives[i]); }
  setBoxDerivatives(boxDerivatives);
  setValue(cv_val);
}
\endverbatim

\section multicvs Mult-component CVs

To avoid code duplication, and in some cases computational expense, plumed has functionality so that a single line in input can calculate be used to calculate multiple components for a CV.  For example, PATH computes the distance along the path,\f$s\f$, and the distance from the path, \f$z\f$.  Alternatively, a distance can give one the \f$x\f$, \f$y\f$ and \f$z\f$ components of the vector connecting the two atoms.  You can make use of this functionality in your own CVs as follows:

- In the constructor we create an additional value for the CV by adding the call PLMD::addValueWithDerivative("new") as well as PLMD::addValueWithDerivatives(””).  In addtion set any periodicity for our component using getValue("new")->setPeridicity() and getValue("new")->setDomain(min,max).  If this CV is called plum in our input file we can now use both plum and plum.new in any of the functions/methods in plumed.
- Obviously in calculate we now must provide functionality to calculate the values, boxDerivatives and the atom derivatives for both our original plum and its component plum.new. Furthermore, all of this data must be transferred to the plumed core.  This is done by the following code:

Here we transfer the value, box derivatives and atomic derivatives for plum.
\verbatim
for(int i=0;i<derivatives.size();i++){ setAtomsDerivatives(i,derivatives[i]); }
setBoxDerivatives(boxDerivatives);
setValue(cv_val);
\endverbatim
Here we transfer the value, box derivatives and atomic derivatives for plum.new.
\verbatim
Value* nvalue=getValue("new");
for(int i=0;i<nderivatives.size();i++){ setAtomsDerivatives(nvalue i,nderivatives[i]); }
setBoxDerivatives(nvalue,nboxDerivatives);
setValue(nvalue,ncv_val);
\endverbatim

Please only use this functionality for CVs that are VERY similar.
\page AddingACLTool How to add a new command-line tool

To implement a command line tool you need to create a single cpp file call CLToolNAME.cpp. You can, in a command line
tool, use functionality from plumed to perform simple post-processing tasks.  For example, sum_hills uses the 
functionality inside from inside the biasing PLMD::Action, PLMD::BiasMetaD to calculate free energy surfaces.
Regardless, of what you are endeavoring to do your CLToolNAME.cpp file should be formatted in accordance with
the following template:

\verbatim
#include "CLTool.h"
#include "CLToolRegister.h"
#include "PlumedConfig.h"
#include "ActionRegister.h"

using namespace std;

namespace PLMD {
/**
//+PLUMEDOC TOOLS name
\endverbatim

Insert the documentation for your new tool here

\verbatim
\par Examples
\endverbatim

Insert some examples of how to use your tool here 

\verbatim
*/
//+ENDPLUMEDOC

/******* This is how you should declare the class for your command line tool.  main() does
         all the analsysis you require. The constructor and the registerKeywords routine 
         only do anything if you are using one of the template forms of input described below.
         However, reigsterKeywords() must always be declared.

class CLToolNAME:
public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  CLToolNAME(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc);
  string description()const{
    return "a description of the particular kind of task CLToolNAME performs";
  }
};

PLUMED_REGISTER_CLTOOL(CLToolNAME,"name")
\endverbatim

Insert the code for main, registerKeywords and the constructor here.

\verbatim
}  /--- don't forget to close the namespace PLMD { at the end of the file
\endverbatim

\section input Input

The input stream is passed to main so you can create tools that have an input in any form.
However, there are certain standardized forms of input that can be used for command line
tools. If your tool takes its input in one of these forms we strongly encourage you to
use the code that is already present.  If you do so a great deal of manual generation will 
be looked after automatically.

There are two command line input types that are implemented in the base class.  The first
is for tools such as driver where the input is specified using a series of command line flags
e.g.

\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --dump-forces
\endverbatim

The other are tools like simplmd that take an input file that contains one directive per line
and a corresponding value for that directive.  For both these forms of input it is possible to
read in everything that is required to run the calculation prior to the actual running of the 
calculation.  In other words these the user is not prompted to provide input data once the main 
calculation has started running (N.B. you can do tools with input of this sort though as the input stream
is passed to main).  

If you wish to use driver-like or simple-md like input then you have to specify this in the constructor.
For driver-like input you would write the following for the constructor:

\verbatim
CLToolNAME::CLToolNAME(const CLToolOptions& co ):
CLTool(co)
{
 inputdata=commandline;
}
\endverbatim

For simplemd-like input you write the following for the constructor:

\verbatim
CLToolNAME( const CLToolOptions& co ) :
CLTool(co)
{
  inputdata=ifile;
}
\endverbatim

If you are not using one of these input forms then you don't need to write a constructor although
you may choose to for reasons not connected to the input of data.

Once you have created the constructor the actual readin and manual creation is done in much the same
manner as it is done for Actions in the main plumed code (\ref usingDoxygen). You write a
registerKeywords routine as follows:

\verbatim
void CLToolNAME::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys );
}  
\endverbatim

Inside this routine you add descriptions of all your various command-line flags (driver-like) or
input directives (simple-md-like) and these descriptions are used to create the manual. The code 
will automatically check if a particular flag is present and read any input directive connected 
with a particular flag (i.e. the data after the space). Within main you can recover the read in 
data using CLTool::parse and CLTool:parseFlag.   

\section getplumed Re-using plumed

To re-use the functionality that is present in plumed you use the same tools that are used to 
patch the various MD codes (\ref HowToPlumedYourMD).  Alternatively, if you want to create
an instance of a particular Action you can do so by issuing the following commands:

\verbatim
PlumedMain* plumed=new PlumedMain(); std::vector<std::string> words;
Action* action=actionRegister().create(ActionOptions(plumed,words));
delete plumed; delete action;
\endverbatim

Please be aware that words should contain everything that would be required in an input
line to make a valid instance of the Action you require.
\page AddingAMultiColvar How to add a new MultiColvar

As you are no doubt aware within plumed 2 you can calculate multiple 
instances of a collective coorinate from a single line in the input file.
One can then calculate functions such as the minimum, number less than,...
from the resulting distribution of collective variables.  To create these kinds
of collective variables we use the functionality implemented in PLMD::multicolvar::MultiColvar.
In fact in writing a single PLMD::multicolvar::MultiColvar you actually write many CVs at once
as the minimum of a distribution, the number less than and so on come with no 
additional effort on your part.

To better understand how to go about writing a new MultiColvar examine the interior
of one of the existing MultiColvars, e.g. PLMD::multicolvar::Distances or PLMD::multicolvar::CoordinationNumbers.
In fact a good way to start is to copy one of these files as certain features (e.g. the fact that you
have to include MultiColvar.h and ActionRegister.h and that all your code should be inside the namespace
PLMD) never change. 

The essential idea with multicolvar is that you can reuse the same functionality within many collective variables.
In particular, you can exploit the same parallelisation strategy in a wide variety of different collective variables.
Similarly one can use the same implementation of link cells in a wide variety of different cvs.  Lastly, one can
create functions of multicolvars and thus use the calculated cvs in a wide variety of different contexts.  If you are 
reading this you are probably aware of this and are (hopefully) impressed by the flexibility of this objects. 

The reason for this flexibility is that many CVs can be calculated using a variant on the following formula:

\f[
s = \sum_i g[f(\{X\}_i)]
\f]

where \f$g\f$ and \f$f\f$ are functions and the sum runs over a sets of atomic positions \f$\{X\}_i\f$.
In the above formula the sum and a set of likely \f$g\f$ functions are implemented within the base classes
base classes of multicolvar PLMD::multicolvar::MultiColvarBase, PLMD::multicolvar::MultiColvar, PLMD::multicolvar::MultiColvarFunction and 
PLMD::multicolvar::BridgedMultiColvarFunction.  When you implement a new multicolvar what you thus need to do is 
write code to calculate the function \f$f\f$ from a set of atomic positions.  One way to understand how to do this would be to look
through these classes and attempt to understand the code contained within them.  The following is an alternative.  It
provides a rather prescriptive description as to how to implement a new multicolvar.  I hope this helps but I am keenly
aware that any description is bound to be incomplete so if you have questions please do email the user list to ask.

\section AddingAMultColvarDocs Creating documentation  

The first thing you will have to change is the documentation.  As discussed on the \ref usingDoxygen page
of this manual the documentation is created using Doxygen.  You are implementing a cv so your PLMEDOC line
should read:

\verbatim
//+PLUMEDOC MCOLVAR MYCVKEYWORD 
\endverbatim 

Your documentation should contain a description of what your CV calculates and some examples.  You do not
need to write a description of the input syntax as that will be generated automatically.  For more information
on how to write documentation go to \ref usingDoxygen.

\section class Creating the class

The first step in writing the executable code for your CV is to create a class that calculates your CV.  The
declaration for your class should appear after the documentation and will look something like:

\verbatim
class MyNewMultiColvar : public MultiColvar {
private:
   // Declare all the variables you need here  
public:
  static void registerKeywords( Keywords& keys );
  MyNewMultiColvar(const ActionOptions&);
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms );
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(MyNewMultiColvar,"MYCVKEYWORD")
\endverbatim

This new class (MyNewMultiColvar) inherits from MultiColvar and so contains much of the functionality we 
require already.  Furthermore, by calling PLUMED_REGISTER_ACTION we have ensured that whenever the keyword
MYCVKEYWORD is found in the input file an object of type MyNewMultiColvar will be generated and hence that your 
new CV will be calculated wherever it is required.

The four functions that are defined in the above class (registerKeywords, the constructor, isPeriodic and compute) are mandatory.  
Without these functions the code will not compile and run.  Writing your new CV is simply a matter of writing these
four subroutines.

\section register RegisterKeywords

RegisterKeywords is the routine that is used by plumed to create the remainder of the documentation.  As much of this
documentation is created inside the MultiColvar class itself the first line of your new registerKeywords routine must
read:

\verbatim
MultiColvar::registerKeywords( keys )
\endverbatim

as this creates all the documentation for the features that are part of all PLMD::MultiColvar objects.  To see how to create 
documentation that describes the keywords that are specific to your particular CV please read \ref usingDoxygen.

\par Reading the atoms

Creating the lists of atoms involved in each of the colvars in a PLMD::MultiColvar is quite involved as you have to generate all the sets of 
atoms for which you want to calculate your \f$f\f$ function.  What is more this is
an area where we feel it is important to maintain some consistency in the input. Many of the multicolvars that are currently
implemented in the code thus use one or multiple of the following keywords: 

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td width=5%> ATOMS </td> <td> The atoms keyword specifies that one collective coordinate is to be calculated for each set of atoms specified.  
                               Hence, for MultiColvarDistance the command DISTANCES ATOMS1=1,2 ATOMS2=2,3 ATOMS3=3,4 specifies 
                               that three distances should be calculated. </td>
</tr> <tr>
<td width=5%> GROUP GROUPA GROUPB </td> <td> The GROUP keyword is used for quantities such as distances and angles.  A single GROUP specifies that
                                             a CV should be calculated for each distinct set of atoms that can be made from the group.  Hence,
                                             MUTIDISTANCE GROUP=1,2,3 specifies that three distances should be calculated (1,2), (1,3) and (2,3).
                                             If there is a GROUPA and GROUPB one CV is calculated for each set of atoms that includes at least one
                                             member from each group.  Thus MULTIDISTANCE GROUPA=1 GROUPB=2,3 calculates two distance (1,2) and (1,3). </td>
</tr> <tr>
<td width=5%> SPECIES SPECIESA SPECIESB </td> <td> The SPECIES keywords is used for quantities like coordination numbers.  The way this works is
                                                   best explained using an example.  Imagine a user working on NaCl wishes to calculate the average
                                                   coordination number of the sodium ions in his/her system.  To do this he/she uses COORDINATIONNUMBER SPECIES=1-100
                                                   which tells plumed that atoms 1-100 are the sodium ions.  To calculate the average coordination number
                                                   plumed calculates 100 coordination numbers (one for each sodium) and averages.  Obviously, each of these coordination
                                                   numbers involves the full set of 100 sodium atoms.  By contrast if the user wanted to calculate the coordination number
                                                   of Na with Cl he/she would do COORDINATIONNUMBER SPECIESA=1-100 SPECIESB=101-200, where obviously 101-200 are the 
                                                   chlorine ions.  Each of these 100 heteronuclear coordination numbers involves the full set of atoms speciefied using
                                                   SPECIESB and one of the ions speciefied by SPECIESA. </td>
</tr>
</table>

If you wish to use the ATOMS or SPECIES options above you can do so by adding:

\verbatim
keys.use("ATOMS");
\endverbatim

or 

\verbatim
keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB")
\endverbatim

The documentation for these keywords will then be generated automatically.  If you wish to use the 
something like the GROUP keyword above then you need to register GROUP.  Here is how this is done within 
PLMD::multicolvar::Distances

\verbatim
keys.add("atoms-1","GROUP","Calculate the distance between each distinct pair of atoms in the group");
keys.add("atoms-2","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
                            "the atoms in GROUPB. This must be used in conjunction with GROUPB.");
keys.add("atoms-2","GROUPB","Calculate the distances between all the atoms in GROUPA and all the atoms "
                            "in GROUPB. This must be used in conjunction with GROUPA.");
\endverbatim

To be clear, if these keywords are registered the actual readin of the atoms is looked after by the method
PLMD::multicolvar::MultiColvar::readAtoms.  However, there are a number of things I would like to point out 
about the registering syntax above.  Firstly, writing atoms-1 for GROUP and atoms-2 for GROUPA and GROUPB 
ensures nice formatting in the manual.  In particular it is made clear to the user that these are two 
distinct options for specificying the atoms in the input.   Notice also that we use similar commands to 
register keywords in PLMD::multicolvar::Angles:

\verbatim
keys.add("atoms-1","GROUP","Calculate angles for each distinct set of three atoms in the group");
keys.add("atoms-2","GROUPA","A group of central atoms about which angles should be calculated");
keys.add("atoms-2","GROUPB","When used in conjunction with GROUPA this keyword instructs plumed "
                            "to calculate all distinct angles involving one atom from GROUPA "
                            "and two atoms from GROUPB. The atom from GROUPA is the central atom.");
keys.add("atoms-3","GROUPC","This must be used in conjunction with GROUPA and GROUPB.  All angles "
                            "involving one atom from GROUPA, one atom from GROUPB and one atom from "
                            "GROUPC are calculated. The GROUPA atoms are assumed to be the central "
                            "atoms");
\endverbatim

Again though the actual reading of the GROUP, GROUPA, GROUPB and GROUPC keywords are looked after within
PLMD::multicolvar::MultiColvar::readAtoms

In some cases you may wish to use keywords other than ATOMS, GROUP, SPECIES etc.  In these cases the following
sets of routines may be useful:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td width=10%> PLMD::multicolvar::MultiColvar::readAtomsLikeKeyword </td> <td> This is used to read in the ATOMS keyword.  By using this function you can change the ATOMS keyword to whatever you wish. </td>
</tr> <tr>
<td width=10%> PLMD::multicolvar::MultiColvar::readTwoGroups </td> <td> This is used to read in the GROUPA and GROUPB keywords.  When you use this funciton you can change GROUPA and GROUPB to whatever you desire </td>
</tr> <tr>
<td width=10%> PLMD::multicolvar::MultiColvar::readThreeGroups </td> <td> This is used to read in the GROUPA, GROUPB and GROUPC keywords. 
                                                                          When you use this function you can change GROUPA, GROUPB and GROUPC to whatever you desire.
                                                                          For an example see PLMD::multicolvar::Bridge </td>
</tr> <tr>
<td width=10%> PLMD::multicolvar::MultiColvar::readSpeciesKeyword </td> <td> This is used to read in SPECIESA and SPECIESB keywords.  Once again you can use this funciton and use different keywords. </td>
</tr>
</table>

\section label4 The constructor

Here is the constructor for the COORDINATIONNUMBERS multicolvar:

\verbatim
CoordinationNumbers::CoordinationNumbers(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0){
     switchingFunction.set(sw,errors);
     if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else {
     double r_0=-1.0, d_0; int nn, mm;
     parse("NN",nn); parse("MM",mm);
     parse("R_0",r_0); parse("D_0",d_0);
     if( r_0<0.0 ) error("you must set a value for R_0");
     switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  coordination of central atom and those within %s\n",( switchingFunction.description() ).c_str() );
  // Set the link cell cutoff
  setLinkCellCutoff( switchingFunction.get_dmax() );
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();

  // Read in the atoms
  int natoms=2; readAtoms( natoms );
  // And setup the ActionWithVessel
  checkRead();
}
\endverbatim

This code does three things:

- It reads in all the keywords that are specific to your new CV.  In this case it is a switching function that is read by parse("SWITCH",sw).
- It reads in the atoms involved in the CV.  In this case this is all looked after by the base class method PLMD::multicolvar::MultiColvar::readAtoms
- It sets the cutoff for link cells by exploiting the method PLMD::multicolvar::MultiColvarBase::setLinkCellCutoff.  This is really important - your code will run much faster if you exploit link cells. 
- It checks that readin was successful.  This is the final call to checkRead

Within the constructor of a multicolvar there must be a call to PLMD::multicolvar::MultiColvar::readAtoms even if you choose to read in the atoms using
one of the methods described in the table above.  

As you can see from the above example the majority of the constructor will be concerned with reading keywords from the input directive. More information
and guidance on PLUMED's parsing functionality can be found in \ref parsing.

\section per Implementing isPeriodic

Some \f$f\f$ functions (e.g. TORSIONS) have a co-domain that is periodic the majority, however, do not.  If the \f$f\f$ function of your 
new multicolvar falls into this category your isPeriodic function will be as follows:

\verbatim
bool isPeriodic(){ return false; }
\endverbatim

If your \f$f\f$ function does have a periodic co-domain then you will need to have methods like the following:

\verbatim
bool isPeriodic(){ return true; }
void retrieveDomain( std::string& min, std::string& max ){ min="-pi"; max="pi"; }
\endverbatim

As is hopefully obvious the retrieveDomain should return the minimum and maximum values (as strings) that your function, \f$f\f$, can take. 
The above example is taken from the code in PLMD::multicolvar::Torsions. 

\section label5 Compute

Compute is the routine that actually calculates the function \f$f\f$ and this function's derivatives from the set of atomic positions \f$\{X\}\f$.
This function should have as its return value the value of \f$f\f$.  The atomic positions can be extracted by using the PLMD::multicolvar::AtomValuePack::getPosition method 
that belongs to the PLMD::multicolvar::AtomValuePack class that is passed to the compute method.  One can then calculate the distance vector connecting two atoms by using
PLMD::multicolvar::MultiColvarBase::getSeparation.  Once you have calculated the derivatives of your \f$f\f$ function you can store these within the 
PLMD::multicolvar::AtomValuePack by using the methods PLMD::multicolvar::AtomValuePack::addAtomsDerivatives and PLMD::multicolvar::AtomValuePack::addBoxDerivatives.
 
\section mregtests Writing regression tests

Once you have written all the above function you will need to write regression tests for your new feature if you would like
us to incorporate it in the release version of PLUMED.  Details on how to write regression tests are provided here: \ref regtests  

\page InstallationLayout Installation Layout

I write here some notes related to how plumed is installed.

As a first notice, plumed package is mostly designed to make these tools available:
- a `plumed` executable, that can be used to launch command line tools.
- a `libplumed.so` library, that can be linked to an externally supplied MD code.

These are the main entry points to plumed, but they require several other resources to be
properly located so as to work. Moreover, plumed is designed to be usable both
when "just compiled" (so as to allow for fast development) and when properly installed.
This results in the non trivial problem of knowing where the required resources are located.

Additionally, we provide shell-only alternatives to command line tools. E.g.,
to use `plumed patch` in a cross compiled environment, one can call `plumed-patch` which
does the same but bypass the `plumed` executable and directly executes as a bash script.

As a result, plumed routines and scripts can be entered in three ways:
- calling `plumed` from the command line.
- calling `plumed-patch` from the command line.
- entering the shared library from another code (typically an MD code).

This is achieved in the following way:
- `plumed-*` scripts contains hardcoded environment variables pointing at the correct paths when they are launched.
- `plumed` executable and `libplumed.so` have access to methods (in the config namespace)
  to access to the paths and locate resources.

\warning
Since paths are hardcoded, plumed executable and library are de facto not relocatable.
It is however possible to use them after relocation provided that some environment
variables are set as discussed below.

As an example, the `PLUMED_ROOT` variable is defined to tell to the plumed scripts where to find
most of the plumed-related files. Similarly, from C++ you can use \ref config::getPlumedRoot() to retrieve
the same path.

When a plumed command line tool implemented as script is invoked by the plumed executable,
thus transferring the control from C++ to
an external script, the environment should be consistently set. This is done in method \ref config::getEnvCommand()
which builds a string in the form `env PLUMED_ROOT=/path env PLUMED_INCLUDEDIR=/path ` etc.
In this ways, the scripts are run in an environment with the correct settings.

The following paths need to be set for plumed to work properly. Here they are listed
together with the name of the corresponding environment variables (that can be used in plumed
scripts) and the method that can be used to retrieve them from C++ code.
- Root of plumed: `$PLUMED_ROOT` \ref config::getPlumedRoot()
- Path to include files: `$PLUMED_INCLUDEDIR` \ref config::getPlumedIncludedir()
- Path to html files: `$PLUMED_HTMLDIR` \ref config::getPlumedHtmldir()
- Name of plumed program: `$PLUMED_PROGRAM_NAME` \ref config::getPlumedProgramName()

When using plumed from its build directory (without installing it) these paths will be set to the
value reported below:
- `PLUMED_ROOT=/build/directory`
- `PLUMED_INCLUDEDIR=$PLUMED_ROOT/src/include` (this works thanks to a symlink of `/build/directory/src` to `/build/directory/src/include/plumed`)
- `PLUMED_HTMLDIR=$PLUMED_ROOT`
- `PLUMED_PROGRAM_NAME=plumed`

These paths are hardcoded in `plumed` executable and in `plumed-*` scripts when they are compiled.
Notice that it is possible to set the `PLUMED_ROOT` variable before calling plumed overriding the hard code values.
E.g., you can compile plumed in directory `/build/directory1`, move it to `/build/directory2`, and launch it
with `PLUMED_ROOT=/build/directory2 /build/directory2/src/lib/plumed`. Notice however that although plumed will find all the
required resources in this way, it might not be possible to perform some task such as patching MD code. Also notice that
since the structure of the build directory is fixed the `PLUMED_ROOT` variable is sufficient to reconstruct the other paths.

When using plumed after it has been installed, these paths will be set to the value reported below:
- `PLUMED_ROOT=/usr/local/lib/plumed`
- `PLUMED_INCLUDEDIR=/usr/local/include`
- `PLUMED_HTMLDIR=/usr/local/share/doc/plumed`
- `PLUMED_PROGRAM_NAME=plumed`

These paths are hardcoded in `plumed` executable and in `plumed-*` scripts when they are installed.
Notice that these value can be customized at configure step using standard arguments to ./configure.
When using an installed copy of plumed one can override the hard code values by setting the variables
`PLUMED_ROOT`, `PLUMED_INCLUDEDIR` ,`PLUMED_HTMLDIR`, and `PLUMED_PROGRAM_NAME` before launching plumed.

Notice that to enforce a consistent behavior of scripts and plumed executable the same logic needed to
be implemented twice. One implementation is found in the `src/config/Config.inc.in` file, another implementation
is prependend to the installed scripts by the `src/lib/Makefile`.

Also consider that environment is inherited by subprocesses. That means that if you want to
launch another plumed version from a plumed script (crazy idea, perhaps nobody will ever do it)
you should unexport the relevant environment variables so that the second plumed executable
will find its paths correctly.

\section InstallationLayout-files Installed files

I here describe what's the content of the most important files installed by plumed.

`/usr/local/bin/plumed`: this is a static executable that can be used to launch plumed.
It is typically used to launch a command line tool (e.g. `plumed sum_hills`).
Notice that some command line tools are actually implemented as bash scripts (e.g. `plumed patch`).
Those scripts are located in `$PLUMED_ROOT/scripts/` with an extra `.sh`  suffix. E.g.
the `plumed patch` command will set properly the environment then call `$PLUMED_ROOT/scripts/patch.sh`.

`/usr/local/lib/libplumed.so`: this is a library containing all the plumed routines.
Notice that `/usr/local/bin/plumed` described above is equal to the combination of
`/usr/local/lib/libplumed.so` with a single object file compiled from `buildroot/src/main/main.cpp`.

`/usr/local/lib/libplumedKernel.so`: this is a library containing almost all the plumed routines,
with the exception of those called from MD engines.
Notice that `/usr/local/lib/libplumed.so` described above is equal to the combination of
`/usr/local/lib/libplumedKernel.so` with a single object file compiled from `buildroot/src/wrapper/PlumedStatic.cpp`

`/usr/local/lib/libplumedWrapper.a`: this is a static library containing exclusively the
object file compiled from `buildroot/src/wrapper/Plumed.cpp`

To summarize:
- `bin/plumed` = `buildroot/src/main/main.cpp` + `lib/libplumed.so`
- `lib/libplumed.so` = `buildroot/src/wrapper/PlumedStatic.cpp` + `lib/libplumedKernel.so`
- `lib/libplumedWrapper.a` = `buildroot/src/wrapper/Plumed.cpp`

The logic of this subdivision is that it is possible to either link the MD code to `/usr/local/lib/libplumed.so`
or to link it to a single object file (the one compiled from `buildroot/src/wrapper/Plumed.c` or the installed `libplumedWrapper.a`)
so as to allow linking at run time an a posteriori chosen plumed library. This is the trick behind the `--runtime` patching procedure.

Notice that the only differences between `buildroot/src/wrapper/PlumedStatic.cpp` and `buildroot/src/wrapper/Plumed.c` are that
runtime binding is disabled for the first one and that the second one is compiled as plain C.
This makes it less likely to do mistakes when linking lib/libplumed.so (by unintentionally using a different version
of plumed), and makes C++ library unnecessary if an external code is only interesting in linking the PLUMED
wrappers in `buildroot/src/wrapper/Plumed.c` or in `libplumedWrapper.a`.

We can then dissect more the content of `/usr/local/lib/libplumedKernel.so`.
This library puts together a large list of object files. The same object files will be located after install
in `/usr/local/lib/plumed/obj/k*.o`. I use a wildcard here because these might be many files (named `k0.o`, 'k1.o', etc) or
a single `kernel.o` file (when `ld -r -o` can be used to merge them together). The reason why we 
store object files in the installed directory is that this is the most portable way to link statically C++
objects to another executable. Indeed, merging them in a single .a file (such as libplumed.a) 
would require this library to be linked with special flags so as to allow dropping all the static constructors.
Whereas the special flags could be found by autoconf, it seems simpler to directly link `/usr/local/lib/plumed/obj/k*.o`.


Also notice that this library changes slighlty in the installed version (`/usr/local/lib/libplumedKernel.so`)
and in the pre-install version (`buildroot/src/lib/libplumedKernel.so`). Indeed, whereas the former
include the object file from `buildroot/src/config/ConfigInstall.cpp` the latter includes the object file from
`buildroot/src/config/Config.cpp`. This object file is the one containing the hardcoded paths discussed above,
and thus should include different strings in the installed and pre-install versions.

\note
New in PLUMED v2.5, the `./configure` script will check if it is possible to build a `/usr/local/lib/libplumed.a` library.
This library contains basically `buildroot/src/wrapper/PlumedStatic.cpp` and the single object obtained
merging all the objects in the kernel. When this library is linked, if at least one of the functions in the wrappers
is called (e.g. `plumed_cmd`) then all the objects are pulled in. In principle, this should solve the problem
with C++ static constructors. This feature can be disabled with `--disable-static-archive`.

\section InstallationLayout-installation Installation procedure

When `make` is invoked, several things are performed. First, all the source files are compiled.
The `plumed` executable and the library files are put in `buildroot/src/lib`.
Then, the "to be installed" versions
of the executable and library files are produced and located in `buildroot/src/lib/install`. These are different from
those located in `buildroot/src/lib` in that they include the `buildroot/src/config/ConfigInstall.o` object so as
to hardcode the proper paths.

When `make install` is invoked, the makefile checks if the objects in `buildroot/src/lib/install` should be updated
and, if necessary, recompiles them. If not, it just copies all the material in place. Notice that all the resulting
files are real files (no symlinks). This is a novelty with respect to PLUMED 2.1 and allows for a proper implementation
of the DESTDIR feature required by unix distributions.

Using the standard behavior explained in the autoconf documentation, it is possible to change the paths
for plumed install either during configure (with `--prefix`) or by setting `prefix` during `make install`.


\mainpage Introduction

This is the developer manual. Please first have a look at the <a href="../../user-doc/html/index.html"> user manual </a>.

Plumed 2 is written in C++ and uses many of the advanced, object-oriented features of this language.  This structure makes the implementation of collective coordinates and free energy methods straightforward.  In fact, it should be possible to implement methods and collective coordinates (CV) by creating a single file and without touching any other part of the code. Futhermore, to implement new methodology does not require one to be some sort of C++ wizzard. Rather, the code has been specifically redisigned to make the implementation of new CVs and new free energy methods straightforward so as to encourage people to implement whatever new functionality they require.  This document serves then to provide an introduction as to how to go about implementing new functionality in plumed. A good starting point is \ref INHERIT as this page contains links to parts of the manual where you can find information on how to go about implementing CV, functions and biases. Another useful page is the \subpage TOOLBOX page, which contains information on the many reusable objects that have been implemented in plumed.  

If you want to understand a little more about the code and the way that we use the various features of C++ before you start then we describe this briefly here:

\ref ABriefIntroduction 

And finally, for the developers of MD codes, we provide information as to how to incorperate plumed into your codes here:

\ref HowToPlumedYourMD

If you would like to contribute new functionalities to PLUMED please read the following guidance:

\ref HowToContributeToPlumed

We ask that contributors endeavor to maintain the portability of plumed by, as much as possible, by only using the STL library and lapack in modifications.  
If you need to use any less standard library (e.g. Boost, Sockets) please ensure that your functionality is not installed during a default compilation.  
However, do feel free to provide alternative compilation options that incorperate your functionality.

Information about C++
http://www.parashift.com/c++-faq-lite/

\par Code Coverage

This manual might  also contain a detailed analysis of which parts of the PLUMED code have been tested when compiling it.
In case so, you will find it at <a href="../coverage/index.html"> this link </a>.
If this manual was compiled on Travis-CI, notice that as of PLUMED 2.5 the coverage scan ends up on a separate 
repository named `github.com/plumed/coverage-branchname`.

\defgroup TOOLBOX Tool Box
@{
Classes providing basic tools in plumed.

Classes of this group are designed to be reusable and to incorporate all sorts of functionality in plumed.
We try to keep their documentation as complete and clear as possible so as to increase the
chance that they will be reused.

If you implement a new class that you think might be useful to others please add it to the 
list by including the following inside the header file.
\verbatim
\ingroup TOOLBOX
\endverbatim
@}

\defgroup MULTIINHERIT Classes for multiple inheritance
@{
Classes for multiple inheritance.

Each of these classes implements some special feature which can
be then used to compose complex Actions.
All of them are "public virtual" derivatives of PLMD::Action,
so that it is possible to build ad Action which is based on multiple
classes from this group. This is the only place in the Action hierarchy
where multiple inheritance should be used.

Multiple inheritance allows for immediate combination of these features,
but add some C++ subtleties. If you do not fully understand them don't worry
and directly inherits from classes in the \ref INHERIT group.

To add a class to this group, just put a
\verbatim
\ingroup MULTIINHERIT 
\endverbatim
statement somewhere inside the header.
@}

\defgroup INHERIT Base classes for CVs, functions, biases, etc.
@{
Classes which can be used to create CVs, functions, biases and so on.

The typical way to add a new feature to plumed is to create a new class
which just inherits from one of the classes of this group. For example,
a new collective variable can be created by inheriting from PLMD::Colvar.
Most of the \ref INPUTDIRECTIVES are based on classes from this group.

To add a class to this group, just put a
\verbatim
\ingroup INHERIT
\endverbatim
statement somewhere inside the header.
@}

\defgroup INPUTDIRECTIVES Classes providing input directives
@{
Classes which implement directive that we used in the plumed input file.

Most of these classes are only used to provide a new feature which will be
available from the plumed input file. As such, they are typically not reused in
other places of the code. For this reason, almost all of them are directly
provided into an implementation file (.cpp), and have no accociated header file (.h).
A notable exceptions is PLMD::SetupMolInfo, which needs to be accessed directly
from other classes.




Each of these classes provides one directive for the plumed input file.
This list is built automatically based on the PLUMED_REGISTER_ACTION macro.
@}


\page AddingAFunction How to add a new function

Many collective variables are a function of a number of some set of simpler collective variables. 
These sorts of collective variables should be implemented should be implemented in plumed as functions so as not to duplicate code.

Much like CVs you can implement a function by creating a single cpp file called FunctionNAME.cpp.  
If you would like us to incorporate your CV in the release version of PLUMED you will need to write at least
one regression test for your function.  Details on how to write regression tests are provided here: \ref regtests

If you use the following template for your new FunctionNAME.cpp file then the manual and the calls to the CV will be looked after automatically.

\verbatim
#include "ActionRegister.h"
#include "Function.h"

#include <cmath>
#include <cassert>
using namespace std;
namespace PLMD{

//+PLUMEDOC FUNCTION COMBINE
/*
\endverbatim

At this point you provide the description of your function that will appear in the manual along with a description of the input file syntax and an example.  Merging new features of the code into the plumed main branch without proper documentation is punishable by death!  Some instructions as to how to format this information is provided here: \ref usingDoxygen

\verbatim
*/
//+ENDPLUMEDOC

/**** We begin by declaring a class for your function.  This class inherits everything from the function class.
      This ensures it has a label, a place to store its value, places to the store the values of the derivatives
      and that it knows which of the colvar objects its value depends on.
class FunctionNAME :
  public Function
{
\endverbatim

Insert declarations for your function's parameters here.

\verbatim
  /---- This routine is used to create the descriptions of all the keywords used by your new function 
  static void registerKeywords( Keywords& keys );
  /---- This is the constructor for your function.  It is this routine that will do all the reading.
       Hence it takes as input a line from the input file.
  FunctionNAME(const ActionOptions&);
  /---- This is the routine that will be used to calculate the value of the function, whenever its calculation is required.
        This routine and the constructor above must be present - if either of them are not the code will not compile.
  void calculate();
};

  /------ The following command inserts your new function into plumed by inserting calls to your new
          routines into the parts of plumed where they are required.  This macro takes two arguments:
          The first is the name of your FunctionClass and the second is the keyword for your function
          (the first word in the input line for your function).
PLUMED_REGISTER_ACTION(FunctionNAME,"COMBINE")

/----- The following routine creates the documentation for the keyowrds used by your CV
void FunctionNAME::registerKeywords( Keywords& keys ){
  Function::registerKeywords(keys);
\endverbatim

In here you should add all your descriptions of the keywords used by your colvar. Descriptions as to how to
do this can be found here: \ref usingDoxygen

\verbatim
}

FunctionNAME::FunctionNAME(const ActionOptions&ao):
/--- These two lines set up various things in the plumed core which functions rely on.
Action(ao),
Function(ao)
{
\endverbatim

Insert code here to read the arguments of the function here using plumed's \ref parsing.  N.B. The label and arguments (i.e. the cvs on which the function depends are read in already elsewhere.

\verbatim

/---- For a number of the free energy methods in plumed it is necessary to calculate the
      distance between two points in CV space.  Obviously, for periodic CVs one must take
      periodicities into account when calculating distances and use the minimum image
      convention in distance calculations.  Functions too are used as cvs in these methods
      and thus it is necessary to provide periodicities for these objects too.  In theory it
      should be possible to determine the periodicity of a function from the periodicity of the
      underlying CVs.  However, in practise this is very difficult to do.  We therefore recommend
      that you include the following few lines of code so that the periodicity of functions can
      be specified by the user in input.
  vector<string> period;
  double min(0),max(0);
  parseVector("PERIODIC",period);
  if(period.size()==0){
  }else if(period.size()==1 && period[0]=="NO"){
    getValue("")->setPeriodicity(false);
  } else if(period.size()==2 && Tools::convert(period[0],min) && Tools::convert(period[1],max)){
    getValue("")->setPeriodicity(true);
    getValue("")->setDomain(min,max);
  }
  checkRead();    /--- This command checks that everything on the input line has been read properly

  /--- The following line informs the plumed core that we require space to store the
       value of the function and the derivatives. 
  addValueWithDerivatives("");
}

\verbatim
void FunctionCombine::calculate(){
/--- These are the things you must calculate for any function ---/
  double cv_val;              /--- The value of the function ----/
  vector<double> derivatives; /--- The derivative of the function with respect to the cvs ---/
\endverbatim

Insert code here to calculate your function and its derivatives with repsect to the underlying cvs here. Please use, where possible, the library of tools described in \ref TOOLBOX.

\verbatim
  /---- Having calculated the function and its derivatives you now transfer this information
        to the plumed core using the following two commands. 
  for(int i=0;i<derivatives.size();i++){ setAtomsDerivatives(i,derivatives[i]); }
  setValue(cv_val);
}

}
\endverbatim

\section multicvsf Multi-component functions

To avoid code duplication, and in some cases computational expense, plumed has functionality so that a single line in input can calculate be used to calculate multiple components for a function.  You can make use of this functionality in your own CVs as follows:

- In the constructor we create an additional value for the function by adding the call PLMD::addValueWithDerivative("new") as well as PLMD::addValueWithDerivatives(””).  Please note you must provide keywords in input to provide periodicity information for both your function and all of its components.  The periodicities of the components should then be set using getValue("new")->setPeridicity() and getValue("new")->setDomain(min,max). If we call this function flum in our input file we can now use both flum and flum.new in any of the functions/methods in plumed.
- Obviously in calculate we now must provide functionality to calculate the values and derivatives for both flum and its component flum.new. Furthermore, all of this data must be transferred to the plumed core.  This is done by the following code:

Here we transfer the value and derivatives for flum.
\verbatim
for(int i=0;i<derivatives.size();i++){ setDerivatives(i,derivatives[i]); }
setValue(cv_val);
\endverbatim
Here we transfer the value and derivatives for plum.new.
\verbatim
Value* nvalue=getValue("new");
for(int i=0;i<nderivatives.size();i++){ setDerivatives(nvalue i,nderivatives[i]); }
setValue(nvalue,ncv_val);
\endverbatim

Please only use this functionality for functions that are VERY similar.

Python wrappers for plumed
==========================

Install using the following command::

     python -m pip install plumed

WARNING: You will need to also build and install the plumed library (see http://www.plumed.org) and make sure the file `libplumedKernel.so` (or `libplumedKernel.dylib`) is available on your system.

You should then make sure the library is found setting the environment variable `PLUMED_KERNEL`::

     export PLUMED_KERNEL=/path/to/libplumedKernel.so
     python
     >>> import plumed
     >>> p=plumed.Plumed()

If you manage multiple plumed versions on your system using tcl environment modules, this should be taken care automatically
by the plumed module.

Alternatively, a pure python solution is::

    >>> import plumed
    >>> os.environ["PLUMED_KERNEL"]="/path/to/libplumedKernel.so"
    >>> p=plumed.Plumed()

Finally, notice that you can set the path to the plumed library directly when declaring a Plumed object::

    >>> import plumed
    >>> p=plumed.Plumed(kernel="/path/to/libplumedKernel.so")

This will allow you to mix different plumed versions in the same python script.

CHANGES: See the PLUMED documentation.
