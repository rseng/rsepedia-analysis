# DEVSIM Installation

## Introduction

Please visit the DEVSIM homepage at [https://devsim.org](https://devsim.org) for the latest information concerning the availability of the software.

Supported Platforms:

* Centos, Red Hat 7
* macOS
* Microsoft Windows (64 bit)

## Installing the Binary Package

### Prerequisites

#### Windows

This runtime may be necessary to run the software.

*  Microsoft Visual C++ Redistributable for Visual Studio 2019
   - [https://www.visualstudio.com/downloads](https://www.visualstudio.com/downloads)
   - [https://aka.ms/vs/16/release/VC_redist.x64.exe](https://aka.ms/vs/17/release/vc_redist.x64.exe)


### Install Anaconda Python

*Please note that most Python 3 environments should work.  The Anaconda distributions are recommended since they provide ready access to the Intel Math Kernel Library (MKL).  As of DEVSIM version 2.1.0, you can configure the software to load alternative libraries implementing the math library functions.*

Install [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) environment.  The following install scripts are available in the ``devsim/scripts`` directory to install the latest version of Miniconda into ``${HOME}/anaconda``. 

* ``scripts/install_miniconda_macos.sh``
* ``scripts/install_miniconda_linux.sh``

The Python environment can be loaded in a terminal using:

```
source ${HOME}/anaconda/bin/activate
```

From a system command prompt, with an active Python environment, install the MKL.

    conda install mkl numpy

### Install DEVSIM

Binary releases are available from [Zenodo](https://doi.org/10.5281/zenodo.1186952), or from [GitHub](https://github.com/devsim/devsim/releases).  Uncompress the file to the desired directory.  The go into the the directory and run the Python installer.

    cd devsim
    python install.py

This script will install DEVSIM into your environment.  If you should need to remove this environment, you can use ``pip``.  Alternatively, you can set the PYTHONPATH environment variable to the ``devsim/lib`` directory.

You should then be able to run an example in the ``testing`` directory.

    cd testing
    python cap2.py

### Troubleshooting

If you have any issues with this process.  Please contact us through the DEVSIM forum at [https://forum.devsim.org](https://forum.devsim.org).

At this time, the macOS version only targets Intel-based systems.  There has been a report that users may have to grant permission for code from "unknown" developers.

## Building from Source

***For the builds in this section, an Anaconda (or Miniconda) installation is required.  Please install using the instructions above.  Ensure that you are running the scripts below from a shell with an activated Python environment, having the ``conda`` command in your path.  Please note that the build scripts will create a conda environment named ``devsim_python3_build``.***

#### Building Centos 7 Linux

***The scripts will install many system packages.  It is highly recommended to run these scripts in a docker container.***

The Centos 7 version is used for binary Linux releases, as it is forward compatible with most Linux distributions and used the Intel MKL.  It will install all required dependencies, including Anaconda Python.  It is recommended to run this in a docker image.  ***There is no need to run ``install_miniconda_linux.sh`` if it is already installed***.

```
git clone https://github.com/devsim/devsim.git
cd devsim
git submodule init
git submodule update
bash scripts/build_centos_7.sh devsim_linux_version
```

where ``version`` is replaced with the desired version.  The compressed file is then ``devsim/dist/devsim_linux_version.tgz``

#### Building Other Linux Systems

***Please note that these scripts may attempt to install packages using the system installer on your system.  Please evaluate the scripts carefully to ensure you agree with these changes to your system.  It is highly recommended to run these scripts in a docker container.  Anaconda Python is not required for these builds.***

Additional build scripts are available for ``Ubuntu`` and ``Fedora`` in the ``scripts`` directory.  These version do not use the Intel MKL and replace the use of the Intel MKL Pardiso Solver with SuperLU.

### Building macOS

#### Release Version

***The builds scripts will make changes to an existing Homebrew installation.***

The macOS release version is built using a Travis build server (see ``.travis.yml``), but it can be built locally with the following prerequisites:

* Xcode 10.1
* Homebrew [https://brew.sh/](https://brew.sh/)
* GCC Compiler (Brew package gcc@9)

The rest of the components, including Anaconda Python and various build components are downloaded and installed in the scripts below.  ***There is no need to run ``install_miniconda_macos.sh`` if Anaconda Python is already installed.***.

```
git clone https://github.com/devsim/devsim.git
cd devsim
git submodule init
git submodule update
source scripts/install_miniconda_macos.sh
bash scripts/build_macos.sh gcc version
```

where ``version`` is replaced with the desired version.  The compressed file is then ``devsim/dist/devsim_macos_version.tgz``
Note that the resulting package file does not have any dependencies on Homebrew on the target system.

#### Using macOS system compiler

Please note that you will need to install the Xcode developer tools, including the C++ compiler for your system.  This version does not have 128-bit floating point precision support and is not provided as part of the regular release.  By default it uses SuperLU 4.3 as the direct solver, but the build scripts can be modified to support Intel MKL Pardiso.  The rest of the components, including Anaconda Python and various build components are downloaded and installed in the scripts below.  ***There is no need to run ``install_miniconda_macos.sh`` if Anaconda Python is already installed.***

```
git clone https://github.com/devsim/devsim.git
cd devsim
git submodule init
git submodule update
source scripts/install_miniconda_macos.sh
bash scripts/build_macos.sh gcc version
```

### Building on Windows

***Please note that the Windows build scripts will make changes to an existing MSYS 2 installation.  It will also create new Anaconda environments.***

#### MSYS version

Please see ``appveyor.yml`` in the source directory for more information about the automated build process.

To build manually, the prerequisites are to have installed:

* Anaconda Python (or Miniconda)
* MSYS 2 [https://www.msys2.org/](https://www.msys2.org/)

First download the source and initialize the submodule from [https://github.com/devsim/devsim.git](https://github.com/devsim/devsim.git).

From a Anaconda environment command prompt:

    scripts\build_appveyor.bat msys conda version

where ``version`` is replaced with the desired version.  The compressed file is then ``devsim/dist/devsim_msys_version.zip``.


#### Visual Studio 2019 version

Extended 128-bit floating point precision is not available in this version.  Please build the MSYS version if this is required.  It requires installing:

* Anaconda Python (or Miniconda)
* MSYS 2 [https://www.msys2.org/](https://www.msys2.org/)
* Microsoft Visual Studio 2019

as well as the prerequisites for the MSYS version above.

From a Anaconda environment command prompt:

    scripts\build_appveyor.bat x64 conda version

where ``version`` is replaced with the desired version.  The compressed file is then ``devsim/dist/devsim_win64_version.zip``



### Citing The Simulator

This software has been issued a DOI which is citable in publications.  Please see the instructions at [Zenodo](https://doi.org/10.5281/zenodo.1186952) for citing this version, or all versions of this simulator.

### Citing The Documentation

The documentation has been issued a DOI which is citable in publications.  Please see the instructions at [Zenodo](https://doi.org/10.5281/zenodo.4583208) for citing this version, or all versions of this simulator.

### Other References

Please visit this projects [Official Website](https://devsim.org/introduction.html#documentation) for additional references concerning this work.



# Contributing

Thank you for your interest in this project.

TCAD simulation is a broad topic utilizing skills in:

* Programming
* Semiconductor Modeling
* Numerics
* Simulation Methods
* Meshing
* Extraction
* Visualization
* Documentation
* Testing

If you would like to participate in this project, please join our email forum:
https://forum.devsim.org
and let us know how you would like to help.

Please see the [Testing](README.md#Testing) and [Related Projects](README.md#Related-Projects) section in [README.md](README.md) for additional places where help is needed.
[![Build Status](https://travis-ci.org/devsim/devsim.svg?branch=master)](https://travis-ci.org/devsim/devsim)
[![Build status](https://ci.appveyor.com/api/projects/status/github/devsim/devsim?branch=master&svg=true)](https://ci.appveyor.com/project/devsim/devsim)

# DEVSIM

## Introduction:
**DEVSIM** - TCAD Device Simulator

**DEVSIM** is a tool for TCAD Device Simulation, using finite volume methods.  The source code is provided by Devsim LLC.

## Website:

The official website is here:

[https://devsim.org](https://devsim.org)

## Installation:

Binary releases are available at [Zenodo](https://doi.org/10.5281/zenodo.1186952).  Please see [INSTALL.md](INSTALL.md) for installation instructions from binary or source.

## Citing This Work:

Please see [CITATION.md](CITATION.md).  *Please do not cite this Github Repository as it will be moving in the future.*

## License:

DEVSIM is licensed under the [Apache 2.0 License](https://www.apache.org/licenses/LICENSE-2.0.html).  Example scripts are also provided under this license.  Other files are subject to the license terms of their copyright owners.
Please see [LICENSE](LICENSE) and [NOTICE](NOTICE) for license terms.

## Documentation:

The PDF documentation is located in ``doc/devsim.pdf``.  An online version of the documentation is available at [https://devsim.net](https://devsim.net).  A list of documentation resources is available [online](https://devsim.org/introduction.html#documentation).

The repository for the documentation is at [https://github.com/devsim/devsim_documentation](https://github.com/devsim/devsim_documentation).

## Supported Platforms:
* macOS 10.13 (High Sierra)
* Microsoft Windows (64 bit)
* Red Hat 7 (Centos Compatible)


## Software Features:
* Python scripting 
* DC, small-signal AC, impedance field method, transient
* User specified partial differential equations (PDE).
* 1D, 2D, and 3D simulation
* 1D, 2D mesher
* Import 3D meshes.
* 2D cylindrical coordinate simulation
* ASCII file format with PDE embedded.

## Support:

For support and general discussion, please join our forum:
[https://forum.devsim.org](https://forum.devsim.org)

## Contributing:

Please see:
[Contribution guidelines for this project](CONTRIBUTING.md)

## Testing

If you wish to see simulation results, or run the regression tests yourself.  They are available from:

| Name | Description |
| --- | --- |
| [devsim_tests_msys](https://github.com/devsim/devsim_tests_msys) | Windows MSYS Build |
| [devsim_tests_win64](https://github.com/devsim/devsim_tests_win64) | Windows Visual Studio 2019 Build |
| [devsim_tests_macos_gcc_x86_64](https://github.com/devsim/devsim_tests_macos_gcc_x86_64) | macOS Build |
| [devsim_tests_linux_x86_64](https://github.com/devsim/devsim_tests_linux_x86_64) | Linux |

*The results are platform dependent due to differences in the compiler, operating system, and math libraries used on each platform.*

## Related Projects

### Used directly by the simulator
| Name | Description |
| --- | --- |
| [symdiff](https://github.com/devsim/symdiff) | Symbolic differentiation engine for the simulator |
| [devsim_documentation](https://github.com/devsim/devsim_documentation) | Documentation for the simulator |

### Extended examples
| Name | Description |
| --- | --- |
| [devsim_bjt_example](https://github.com/devsim/devsim_bjt_example) | Bipolar Junction Transistor example |
| [devsim_density_gradient](https://github.com/devsim/devsim_density_gradient) | Quantum Corrections to Drift Diffusion simulation |

# CHANGES

## Introduction

Please see the release notes in doc/devsim.pdf or at https://devsim.net for more detailed information about changes.

## Version 2.1.0

### Explicit math library loading

#### Introduction

Since the Intel Math Kernel Library started versioning the names of their dynamic link libraries, it has been difficult to maintain a proper Anaconda Python environment when the version has been updated.  With this release, it is possible to use any recent version of the Intel MKL.  In addition, the user is able to load alternative BLAS/LAPACK math libraries.

#### Intel MKL

From DEVSIM Version 2.1.0 onward, a specific version is not required when loading the Intel MKL.  If the Intel MKL is not found, the import of the ``devsim`` module will fail, and an error message will be printed.  This method is the default, and should work when using an Anaconda Python environment with the ``mkl`` package installed.

When using a different Python distribution, or having an installation in a different place, it is possible to specify the location by modifying the ``LD_LIBRARY_PATH`` environment variable on Linux, or using ``DYLD_LIBRARY_PATH`` on macOS.  The explicit path may be set to the MKL math libraries may be set using the method in the next section.

#### Loading other math libraries

It is possible to load alternative implementations of the BLAS/LAPACK used by the software.  The ``DEVSIM_MATH_LIBS`` environment variable may be used to set a ``:`` separated list of libraries.  These names may be based on relative or absolute paths.  The program will load the libraries in order, and stop when all of the necessary math symbols are supplied.  If symbols for the Intel MKL are detected, then the Pardiso direct solver will be enabled.

#### New CMAKE build option
For those building the software, the ``EXPLICIT_MATH_LOAD`` CMAKE option has been added to control the new explicit math loading feature.  An important benefit of this option is that it is possible to build a release version of the software, even if the Intel MKL has not been installed on the build computer.

### Direct solver selection

The direct solver may be selected by using the ``direct_solver`` parameter.
```
devsim.set_parameter(name='direct_solver', value='mkl_pardiso')
```

The following options are available:
- ``mkl_pardiso`` Intel MKL Pardiso
- ``superlu`` SuperLU 4.3

The default is ``mkl_pardiso`` when the Intel MKL is loaded.  Otherwise, the default will switch to ``superlu``.

### Kahan summation in extended precision mode

The ``kahan3`` and ``kahan4`` functions are now using the Kahan summation algorithm for extended precision model evaluation.  Previously, this algorithm was replaced with 128-bit floating point addition and subtraction in releases that support extended precision mode.  With this change, better than 128-bit floating precision is available when extended precision is enabled.
```
devsim.set_parameter(name = "extended_model", value=True)
```

The ``testing/kahan_float128.py`` test has been added.

### Visual Studio 2022

The Microsoft Windows``win64`` release version is now built using the Visual Studio 2022 compiler.  For users needing extended precision on the Windows platform, the ``msys`` build is recommended.

## Version 2.0.1

### Update documentation files

The following files were updated in the text documentation distributed with the software.
- ``CONTRIBUTING.md``
- ``INSTALL.md``
- ``README.md``

This was done to create a version to coincide with this paper in the Journal of Open Source Software.

Sanchez, J. E., (2022). DEVSIM: A TCAD Semiconductor Device Simulator. Journal of Open Source Software, 7(70), 3898, [https://doi.org/10.21105/joss.03898](https://doi.org/10.21105/joss.03898)

### Update MKL Version

The release version of this software is build against version 2 of the Intel MKL, which corresponds to the latest version of Anaconda Python.  If you have issues running ``DEVSIM`` with this new version, please contact us at [https://forum.devsim.org](https://forum.devsim.org) for assistance.

### Update SYMDIFF reference

The submodule reference to ``SYMDIFF`` was updated due some changes to its build files.  Some ``SYMDIFF`` build scripts in the ``DEVSIM`` repository were also updated.

## Version 2.0.0

### Versioned MKL DLL in release build

The Intel Math Kernel Library now uses versioned library names.  Binary releases are now updated against the latest versioned dll names from MKL available in the Anaconda Python distribution.

### Fixed issue in ramp function.

The ``rampbias`` function in the ``devsim.python_packages.ramp`` module has been fixed to properly reduce the bias when there is a convergence failure.

### Transient Simulation

Fixed bug with ``transient_tr`` (trapezoidal) time integration method where the wrong sign was used to integrate previous time steps.

Fixed bug in the charge error calculation, which calculates the simulation result with that a forward difference projection.

Added ``testing/transient_rc.py`` test which compares simulation with analytic result for RC circuit.

Added ``set_initial_condition`` command to provide initial transient conditions based on current solution.

### Create interface from node pairs

Added ``create_interface_from_nodes`` to make it possible to add interface from non-coincident pairs of nodes.

### Solver

#### Convergence Tests

The ``maximum_error`` and ``maximum_divergence`` options where added to the ``solve`` command.  If the absolute error of any iteration goes above ``maximum_error``, the simulation stops with a convergence failure.  The ``maximum_divergence`` is the maximum number of iterations that the simulator error may increase before stopping.

#### Verbosity

During the ``solve``, circuit node and circuit solution information is no longer printed to the screen for the default verbosity level.  In addition, the number of equations per device and region is no longer displayed at the start of the first iteration.

#### SuperLU

The code now supports newer versions of ``SuperLU``.  The release version is still using SuperLU 4.3 for the iterative solution method, and the Intel MKL Pardiso for the direct solve method.

#### Simulation Matrix

The ``get_matrix_and_rhs`` command was not properly accepting the ``format`` parameter, and was always returning the same type.

### Build Scripts

The build scripts have been updated on all platforms to be less dependent on specific Python 3 versions.

An updated fedora build script has been added.  It uses the system installed ``SuperLU`` as the direct solver.

### Documentation Files

Some out of date files (e.g. RELEASE, INSTALL, . . .) have been removed.  The [README.md](README.md) has been updated and the [INSTALL.md](INSTALL.md) have been updated.

### Command Options

The ``variable_name`` option is no longer recognized for the ``devsim.contact_equation`` and ``devsim.interface_equation`` as it was not being used.

## Version 1.6.0

### Array Type Input and Output

In most circumstances, the software now returns numerical data using the Python ``array`` class.  This is more efficient than using standard lists, as it encapsulates a contiguous block of memory.  More information about this class can be found at [https://docs.python.org/3/library/array.html](https://docs.python.org/3/library/array.html).  The representation can be easily converted to lists and ``numpy`` arrays for efficient manipulation.

When accepting user input involving lists of homogenous data, such as ``set_node_values`` the user may enter data using either a list, string of bytes, or the ``array`` class.  It may also be used to input ``numpy`` arrays or any other class with a ``tobytes`` method.

### Get Matrix and RHS for External Use

The ``get_matrix_and_rhs`` command has been added to assemble the static and dynamic matrices, as well as their right hand sides, based on the current state of the device being simulated.  The ``format`` option is used to specify the sparse matrix format, which may be either in the compressed column or compressed row formats, ``csc`` or ``csr``.

### Maximum Divergence Count

If the Newton iteration errors keep increasing for 20 iterations in a row, then the simulator stops.  This limit was previously 5.

### Mesh Visualization Element Orientation

Elements written to the ``tecplot`` format in 2d and 3d have node orderings compatible with the element connectivity in visualization formats.  Specifying the ``reorder=True`` option in ``get_element_node_list`` will result in node ordering compatible with meshing and visualization software.

## Version 1.5.1

### Math Functions

The following inverse functions and their derivatives are now available in the model interpreter.
- ``erf_inv`` Inverse Error Function
- ``erfc_inv`` Inverse Complimentary Error Function
- ``derf_invdx`` Derivative of Inverse Error Function
- ``derfc_invdx`` Derivative of Complimentary Inverse Error Function

The Gauss-Fermi Integral, using Paasch's equations are now implemented.
- ``gfi`` Gauss-Fermi Integral
- ``dgfidx`` Derivative of Gauss-Fermi Integral
- ``igfi`` Inverse Gauss-Fermi Integral
- ``digfidx`` Derivative of Inverse Gauss-Fermi Integral

Each of these functions take two arguments, ``zeta`` and ``s``.  The derivatives with respect to the first argument are provided.  Please see ``testing/GaussFermi.py`` for an example.

In extended precision mode, the following functions are now evaluated with full extended precision.
- ``Fermi``
- ``dFermidx``
- ``InvFermi``
- ``dInvFermidx``

The following double precision tests:

- ``testing/Fermi1.py`` Fermi Integral Test
- ``testing/GaussFermi.py`` Gauss Fermi Integral Test

Have extended precision variants:

- ``testing/Fermi1_float128.py``
- ``testing/GaussFermi_float128.py``

### Installation Script

A new installation script is in the base directory of the package.
It provides instructions of completing the installation to the ``python`` environment without having to set the ``PYTHONPATH`` environment variable.
It notifies the user of missing components to finish the installation within an ``Anaconda`` or ``Miniconda`` environment.


To use the script, use the following command inside of the ``devsim`` directory.

```
    python install.py
```

The install script will write a file named ``lib/setup.py``, which can be used to complete the installation using ``pip``.  The script provides instructions for the installation and deinstallation of ``devsim``.

```
    INFO: Writing setup.py
    INFO:
    INFO: Please type the following command to install devsim:
    INFO: pip install -e lib
    INFO:
    INFO: To remove the file, type:
    INFO: pip uninstall devsim
```

## Version 1.5.0

The ``custom_equation`` command has been modified to require a third return value.  This boolean value denotes whether the matrix entries should be row permutated or not.  For the bulk equations this value should be ``True``.  For interface and contact boundary conditions, this value should be ``False``.

It is now possible to replace an existing ``custom_equation``.

The file ``examples/diode/diode_1d_custom.py`` demonstrates custom matrix assembly and can be directly compared to ``examples/diode/diode_1d.py``.

The ``EdgeNodeVolume`` model is now available for the volume contained by an edge.

The ``contact_equation`` command now accepts 3 additional arguments.

- ``edge_volume_model``
- ``volume_node0_model``
- ``volume_node1_model``

These options provide the ability to do volume integration on contact nodes.

The ``equation`` command has replaced the ``volume_model`` option with:

- ``volume_node0_model``
- ``volume_node1_model``

so that nodal quantities can be more localized.

More details are in the manual.

## Version 1.4.14
### Platforms

Windows 32 bit is no longer supported.  Binary releases of the ``Visual Studio 2019`` ``MSYS2/Mingw-w64`` 64-bit builds are still available online.

On Linux, the releases are now on Centos 7, as Centos 6 has reached its end of life on November 30, 2020.

### C++ Standard

The C++ standard has been raised to C++17.

## Version 1.4.13

The node indexes with the maximum error for each equation will be printed when ``debug_level`` is ``verbose``.

```
devsim.set_parameter(name="debug_level", value="verbose")
```

These are printed as ``RelErrorNode`` and ``AbsErrorNode``:

```
    Region: "gate"	RelError: 5.21531e-14	AbsError: 4.91520e+04
      Equation: "ElectronContinuityEquation"	RelError: 4.91520e-16	AbsError: 4.91520e+04
	RelErrorNode: 129	AbsErrorNode: 129
```

This information is also returned when using the ``info=True`` option on the ``solve`` command for each equation on each region of a device.

If the ``info`` flag is set to ``True`` on the ``solve`` command, the iteration information will be returned, and an exception for convergence will no longer be thrown.  It is the responsibility of the caller to test the result of the ``solve`` command to see if the simulation converged.  Other types of exceptions, such as floating point errors, will still result in a Python exception that needs to be caught.


## Version 1.4.12

Element assembly for calculation of current and charges from the device into the circuit equation are fixed.  These tests are added:

- ``testing/cap_2d_edge.py``
- ``testing/cap_2d_element.py``
- ``testing/cap_3d_edge.py``
- ``testing/cap_3d_element.py``

The ``edge`` variant is using standard edge based assembly, and the ``element`` variant is using element-based assembly.

## Version 1.4.11

The ``element_pair_from_edge_model`` is available to calculate element edge components averaged onto each node of the element edge.  This makes it possible to create an edge weighting scheme different from those used in ``element_from_edge_model``.

Fixed issue where command option names where not always shown in the documentation.

The platform specific notes now clarify that any version of Python 3 (3.6 or higher) is supported.

- ``linux.txt``
- ``windows.txt``
- ``macos.txt``

## Version 1.4.10

Fixed crash when evaluating element edge model in 3D.

Fixed potential error using ``delete_node_model`` and similar deletion commands.

## Version 1.4.9

Support for loading mesh files containing element edge data.

## Version 1.4.8

In transient mode, the convergence test was flawed so that the ``charge_error`` was the only convergence check required for convergence.  The software now ensures all convergence criteria are met.

## Version 1.4.7

### Models

In the simple physics models, the sign for time-derivative terms was wrong for the electron and hole continuity equations.  This affects small-signal and noise simulations.  The example at ``examples/diode/ssac_diode.py`` was updated to reflect the change.

### Platforms

Fix build script issue for macOS on Travis CI, updated the compiler to ``g++-9``.

Update Centos 6 build from ``devtoolset-6`` to ``devtoolset-8``.

## Version 1.4.6

### Version Information

Parameter ``info`` can be queried for getting version information.  The file ``testing/info.py`` contains an example.

```
  python info.py
  {'copyright': 'Copyright © 2009-2020 DEVSIM LLC', 'direct_solver': 'mkl_pardiso', 'extended_precision': True, 'license': 'Apache License, Version 2.0', 'version': '1.4.6', 'website': 'https://devsim.org'}
```

### Extended Precision

The example ``examples/diode/gmsh_diode3d_float128.py`` provides an example where extended precision is enabled.

### Python Formatting

The Python scripts in the ``examples`` and ``testing`` directories have been reformatted to be more consistent with language standards.

### Platforms

Microsoft Windows 10 is supported and is now compiled using Microsoft Visual Studio 2019.

Microsoft Windows 7 is no longer supported, as Microsoft has dropped support as of January 14, 2020.

### External Meshing

Support for reading meshes from Genius Device Simulator has been completely removed from DEVSIM.


## Version 1.4.5

* Platform Support:
  * An MSYS2/Mingw-w64 build is available for 64-bit Windows.  This build, labeled ``devsim_msys_v1.4.5``, enables the use of the 128-bit floating point precision already available on the macOS and Linux platforms.

## Version 1.4.4

* Bug Fixes:
  * Intermittent crash on Windows 10 at the end of the program
* CHANGES.md containing version changes in markdown format.
* Internal changes:
  * Regression system script refactored to Python.
  * Refactor threading code using C++11 functions
  * Refactor timing functions for verbose mode using C++11 functions.
  * Refactor FPE detection code to C++11 standard.

---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Before opening an issue, please create a discussion concerning your request at** https://forum.devsim.org

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**This template is for bug reporting. For general questions about the simulator, please create a new discussion at** https://forum.devsim.org

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
