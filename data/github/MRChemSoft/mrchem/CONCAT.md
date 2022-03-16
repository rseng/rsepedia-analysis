## Versioning scheme

This project adheres to a [semantic versioning scheme](https://semver.org/),
with `MAJOR`.`MINOR`.`PATCH`-`LABEL`

- `MAJOR`: Introduces backward incompatible changes to the API
- `MINOR`: Introduces new features without breaking backward compatibility
- `PATCH`: Only simple bugfixes
- `LABEL`: Marks a pre-released version or release candidate

We will try to keep only a *single* line of development on the `master`
branch, and from this branch off to tag new releases. The type of release
(`MAJOR` or `MINOR`) is determined by the nature of the introduced changes.
`MAJOR` changes should not be introduced lightly, and could be held back at pull
request-level in anticipation of prior `MINOR` release(s). Once a new `MAJOR` is
released, the development on the old `MAJOR` is terminated and no `MINOR` pull
requests will be accepted on top of old `MAJOR`s, only simple `PATCH`es.

## Branching model

In the following, X, Y and Z should *always* be substituted with appropriate
numeric values, they should never appear as 'X', 'Y' or 'Z' in any version or
branch name.

The upstream MRChemSoft/mrchem repository should contain only two kinds of
branches: the `master` branch which represents the main development line, as
well as a separate `release/X.Y` branch for each `MINOR` release. All `PATCH`
releases are applied linearly on the `MINOR` release branches.

### The `master` branch:

<img src="doc/gfx/git-master.png" alt="drawing" width="600"/>

- This branch should *not* carry any release tags
- New features should *always* be directed to this branch
- Bugfixes may be directed to this branch
- Bugfixes may be *cherry-picked* from `release` branches, but *never* merged
- The VERSION file should point to the next *expected* release,
  and *always* carry the pre-release label `-alpha`
- When a new `release/X.Y` branch is created, the VERSION file on `master`
  is bumped to the next expected `MAJOR`/`MINOR`, i.e. `X.(Y+1).0-alpha` or
  `(X+1).0.0-alpha`

### The `release/X.Y` branches:

<img src="doc/gfx/git-release.png" alt="drawing" width="530"/>

- This branch should carry *all* release tags associated with the `X.Y` `MINOR`
  release, including `PATCH`es and release candidates `-alpha1`, `-alpha2`, etc
- The VERSION file should point to the *latest* release tag on this branch
- The creation of this branch marks a feature freeze for version `X.Y`
- New features should *never* be directed to this branch
- Bugfixes and release preparations may be directed to this branch
- Bugfixes may be *cherry-picked* from `master`, but *never* merged

## Contributing

All code changes are incorporated through the `fork -> pull request (PR) ->
code review` work flow, with the following principles:

- New features should *always* branch off `master`
- Bugfixes may branch off either `master` or `release/X.Y`
- PRs should *always* be directed back at its original base branch, `master` or
  `release/X.Y`
- Bugfixes should be small and specific (cherry-pickable)
- Cherry-picks between `master` and `release` branches will be handled by code
  administrators
- All version tagging and changes to the VERSION file will be handled by code
  administrators, and should *not* be part of any PR

### Contributing a new feature:

- Always start from latest `master`
- Branch off to a local feature branch
- Implement new feature
- Regularly incorporate the latest changes from `master`, by merging or
  (preferably) rebasing
- File PR from the local feature branch back to `master`

### Contributing a bugfix:

- Start from latest appropriate branch, `master` or `release/X.Y`
- Branch off to a local bugfix branch
- Implement bugfix
- Regularly incorporate the latest changes from original branch by *rebasing*
- File PR from the local bugfix branch back to original branch
- Evaluate whether the bugfix should be cherry-picked to other branches
  and communicate it to the administrators


## Automatic formatting

You can install Git hooks to keep in check formatting and licensing headers:

```
cd .git/hooks
ln -s ../../.githooks/* .
```

![MRChem logo](https://github.com/MRChemSoft/mrchem/raw/master/doc/gfx/logo_full.png)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3606658.svg)](https://doi.org/10.5281/zenodo.3606658)
[![License](https://img.shields.io/badge/license-%20LGPLv3-blue.svg)](../master/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/mrchem/badge/?version=latest)](http://mrchem.readthedocs.io/en/latest/?badge=latest)
![Build and test MRChem](https://github.com/MRChemSoft/mrchem/workflows/Build%20and%20test%20MRChem/badge.svg)
[![CircleCI](https://circleci.com/gh/MRChemSoft/mrchem/tree/master.svg?style=svg)](https://circleci.com/gh/MRChemSoft/mrchem/tree/master)
[![codecov](https://codecov.io/gh/MRChemSoft/mrchem/branch/master/graph/badge.svg)](https://codecov.io/gh/MRChemSoft/mrchem)

MRChem is a numerical real-space code for molecular electronic structure
calculations within the self-consistent field (SCF) approximations of quantum
chemistry (Hartree-Fock and Density Functional Theory).

The code is being developed at the Hylleraas Centre for Quantum Molecular
Sciences at UiT - The Arctic University of Norway.

### User support: [mrchem.slack.com](https://join.slack.com/t/mrchem/shared_invite/enQtNTI3MjMzNjM0NTk0LWNkODZjNTMwYmM4NmRmODExMjQzMDc3NThlMzNmNmIyNWQwM2YwOGY0OWY4NmNmNzE4ZmM2NzgxYzUzNDg3NDM)
### Documentation: [mrchem.readthedocs.io](http://mrchem.readthedocs.io)


## Installation

For optimal performance it is recommended to build from source, as the packaged
builds are quite generic without architecture specific optimizations.


### From source

To build MRChem from source with MPI+OpenMP parallelization:

    $ git clone https://github.com/MRChemSoft/mrchem.git
    $ cd mrchem
    $ ./setup --prefix=<install-dir> --omp --mpi --cxx=<mpi-compiler> <build-dir>
    $ cd <build-dir>
    $ make
    $ make test
    $ make install

All dependencies will be fetched at configure time, if not already available.
For more information on different kinds of builds, see
[installation instructions](http://mrchem.readthedocs.io/en/latest/installation.html).


### Using Conda

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mrchem/badges/version.svg)](https://anaconda.org/conda-forge/mrchem)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mrchem/badges/latest_release_date.svg)](https://anaconda.org/conda-forge/mrchem)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mrchem/badges/downloads.svg)](https://anaconda.org/conda-forge/mrchem)

To install MRChem in a Conda environment `myenv`:

    $ conda create -n myenv
    $ conda activate myenv
    $ conda install -c conda-forge mrchem               # latest version (OpenMP only)
    $ conda install -c conda-forge mrchem=1.0.0         # tagged version (OpenMP only)
    $ conda install -c conda-forge mrchem=*=*openmpi*   # latest version (MPI+OpenMP)
    $ conda install -c conda-forge mrchem=*=*mpich*     # latest version (MPI+OpenMP)

To list all available versions

    $ conda search -c conda-forge mrchem


### Using Spack

To install MRChem in a Spack environment `myenv`:

    $ spack env create myenv
    $ spack env activate myenv
    $ spack install mrchem                              # latest version (MPI+OpenMP)
    $ spack install mrchem @1.0.0                       # tagged version (MPI+OpenMP)
    $ spack install mrchem -mpi                         # latest version (OpenMP only)

For information on available Spack builds:

    $ spack info mrchem


### Using EasyBuild

To install MRChem in an EasyBuild/Lmod environment (only MPI+OpenMP version
available):

    $ eb MRChem-<version>-<toolchain> --fetch
    $ eb MRChem-<version>-<toolchain> --robot
    $ module load MRChem/<version>-<toolchain>

See
[EasyBuild](https://github.com/easybuilders/easybuild-easyconfigs/tree/develop/easybuild/easyconfigs/m/MRChem)
for available `<versions>` and `<toolchains>`.


### Using Singularity

Singularity recipe files are provided under `recipes/` for building local container images using
the current state of the source. Requires Singularity >= v3.2 as well as `sudo` rights on the
machine you are building on:

    $ sudo singularity build <image_name>.sif recipes/Singularity.<variant>

Recipes are provided for a pure OpenMP build (`recipes/Singularity.nompi`) and one MPI+OpenMP version,
using `OpenMPI-4.0` (`recipes/Singularity.openmpi4.0`).

Official MRChem images can also be downloaded from the GitHub Container Registry.

Latest `master` version (here OpenMP variant):

    $ singularity pull oras://ghcr.io/MRChemSoft/mrchem/mrchem-nompi:latest

Tagged version (here MRChem-v1.0.2 using OpenMPI-v4.0):

    $ singularity pull oras://ghcr.io/MRChemSoft/mrchem/mrchem-openmpi4.0:v1.0.2

Note that the MPI image requires that a compatible MPI library is installed and
available on the host. For information on how to launch the container:

    $ singularity run-help mrchem-mpi.sif

# Change log

## Version 1.0.0 2020-10-28

### Added

- Installation instructions to README

### Changed

- Updated MRCPP to v1.3.6

### Fixed

- Faulty MPI_INTEGER type

## Version 1.0.0-alpha3 2020-10-10

### Added

- Dagger terms in exchange hessian
- Screening of exchange contributions

### Changed

- Updated MRCPP to v1.3.5
- Updated XCFun to v2.1.0
- Updated parselglossy to v0.7.0
- Removed runtest dependency
- Removed parselglossy as runtime dependency
- Encapsulated MRChem and MRCPP parallelizations

### Fixed

- Installation path for SAD files
- Misc CMake fixes to enable packaging

## Version 1.0.0-alpha2 2020-06-29

### Added

- New JSON output

### Changed

- Improved mrchem launcher script
- Improved documentation
- Input keywords for v1 defined and fixed
- Updated MRCPP to v1.2.0
- Updated XCFun to v2.0.1
- Updated nlohmann/json to v3.6.1

## Version 1.0.0-alpha1 2020-05-05

### Features

- Hartree-Fock
- Kohn-Sham DFT (LDA/GGA/hybrid)
- Restricted (closed-shell) and unrestricted
- External electric field
- Ground state energy
- Dipole moment
- Quadrupole moment
- Polarizability (static/dynamic)
- Magnetizability
- NMR shielding
- Density/orbital plots

## Version 0.2.2 2019-11-20

### Fixed

- Updated MRCPP to v1.0.2
- OpenMPI error with complex data types

## Version 0.2.1 2019-08-12

### Fixed

- Updated MRCPP to v1.0.1
- Eigen installation

# How to update the input parser

Run:

``` bash
$ cd python
$ parselglossy generate --template template.yml --docfile user_ref.rst --doc-header="User input reference" --target="mrchem/input_parser"
```
Remember to also update the documentation:

``` bash
cp python/mrchem/input_parser/docs/user_ref.rst doc/users/user_ref.rst
```
This file was automatically generated by parselglossy on 2022-03-09
Editing is *STRONGLY DISCOURAGED*
The CMake infrastructure for this project is generated using [Autocmake]
by Radovan Bast, Roberto Di Remigio, Jonas Juselius and contributors.
The `update.py` Python script and the contents of the directories `autocmake` and `downloaded` are licensed
under the terms of the [BSD-3-Clause license], unless otherwise stated.

[Autocmake]: http://autocmake.org
[BSD-3-Clause license]: https://tldrlegal.com/license/bsd-3-clause-license-(revised)### Generate new recipes using HPC Container Maker (HPCCM)

    $ hpccm --recipe recipe_<variant>.py --format singularity --singularity-version=3.2 > Singularity.<variant>

### Build Singularity image locally (must be done from project root directory)

    $ sudo singularity build <image_name>.sif recipes/Singularity.<variant>

### Pull Singularity image from GitHub Container Registry

Latest `master` version (here OpenMP variant):

    $ singularity pull oras://ghcr.io/MRChemSoft/mrchem/mrchem-nompi:latest

Tagged version (here MRChem-v1.0.2 using OpenMPI-v4.0):

    $ singularity pull oras://ghcr.io/MRChemSoft/mrchem/mrchem-openmpi4.0:v1.0.2

### Run Singularity container (non MPI)

    $ singularity exec <image-name>.sif mrchem molecule

### Run Singularity container (MPI)

    $ singularity exec <image-name>.sif mrchem -D molecule
    $ mpirun singularity exec <image-name>.sif mrchem.x molecule.json
.. raw:: html

    <style> .red {color:#aa0060; font-weight:bold; font-size:18px} </style>

.. role:: red

.. This documentation was autogenerated using parselglossy. Editing by hand is not recommended.

====================
User input reference
====================

- Keywords without a default value are **required**.
- Default values are either explicit or computed from the value of other keywords in the input.
- Sections where all keywords have a default value can be omitted.
- Predicates, if present, are the functions run to validate user input.

:red:`Keywords`
 :world_prec: Overall relative precision in the calculation. 

  **Type** ``float``

  **Predicates**
    - ``1.0e-10 < value < 1.0``

 :world_size: Total size of computational domain given as 2**(``world_size``). Always cubic and symmetric around the origin. Negative value means it will be computed from the molecular geometry. 

  **Type** ``int``

  **Default** ``-1``

  **Predicates**
    - ``value <= 10``

 :world_unit: Length unit for *all* coordinates given in user input. Everything will be converted to atomic units (bohr) before the main executable is launched, so the JSON input is *always* given in bohrs. 

  **Type** ``str``

  **Default** ``bohr``

  **Predicates**
    - ``value.lower() in ["bohr", "angstrom"]``

 :world_origin: Global gauge origin of the calculation. 

  **Type** ``List[float]``

  **Default** ``[0.0, 0.0, 0.0]``

  **Predicates**
    - ``len(value) == 3``

:red:`Sections`
 :Precisions: Define specific precision parameters. 

  :red:`Keywords`
   :exchange_prec: Precision parameter used in construction of Exchange operators. Negative value means it will follow the dynamic precision in SCF. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :helmholtz_prec: Precision parameter used in construction of Helmholtz operators. Negative value means it will follow the dynamic precision in SCF. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :poisson_prec: Precision parameter used in construction of Poisson operators. 
  
    **Type** ``float``
  
    **Default** ``user['world_prec']``
  
    **Predicates**
      - ``1.0e-10 < value < 1.0``
  
   :nuclear_prec: Precision parameter used in smoothing and projection of nuclear potential. 
  
    **Type** ``float``
  
    **Default** ``user['world_prec']``
  
    **Predicates**
      - ``1.0e-10 < value < 1.0``
  
 :Printer: Define variables for printed output. 

  :red:`Keywords`
   :print_level: Level of detail in the written output. Level 0 for production calculations, negative level for complete silence. 
  
    **Type** ``int``
  
    **Default** ``0``
  
   :print_mpi: Write separate output from each MPI to file called ``<file_name>-<mpi-rank>.out``. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :print_prec: Number of digits in property output (energies will get twice this number of digits). 
  
    **Type** ``int``
  
    **Default** ``6``
  
    **Predicates**
      - ``0 < value < 10``
  
   :print_width: Line width of printed output (in number of characters). 
  
    **Type** ``int``
  
    **Default** ``75``
  
    **Predicates**
      - ``50 < value < 100``
  
 :Plotter: Give details regarding the density and orbital plots. Three types of plots are available, line, surface and cube, and the plotting ranges are defined by three vectors (A, B and C) and an origin (O): ``line``: plots on line spanned by A, starting from O. ``surf``: plots on surface spanned by A and B, starting from O. ``cube``: plots on volume spanned by A, B and C, starting from O. 

  :red:`Keywords`
   :path: File path to plot directory. 
  
    **Type** ``str``
  
    **Default** ``plots``
  
    **Predicates**
      - ``value[-1] != '/'``
  
   :type: Type of plot: line (1D), surface (2D) or cube (3D). 
  
    **Type** ``str``
  
    **Default** ``cube``
  
    **Predicates**
      - ``value.lower() in ['line', 'surf', 'cube']``
  
   :points: Number of points in each direction on the cube grid. 
  
    **Type** ``List[int]``
  
    **Default** ``[20, 20, 20]``
  
    **Predicates**
      - ``all(p > 0 for p in value)``
      - ``not (user['Plotter']['type'] == 'line' and len(value) < 1)``
      - ``not (user['Plotter']['type'] == 'surf' and len(value) < 2)``
      - ``not (user['Plotter']['type'] == 'cube' and len(value) < 3)``
  
   :O: Origin of plotting ranges. 
  
    **Type** ``List[float]``
  
    **Default** ``[0.0, 0.0, 0.0]``
  
    **Predicates**
      - ``len(value) == 3``
  
   :A: First boundary vector for plot. 
  
    **Type** ``List[float]``
  
    **Default** ``[1.0, 0.0, 0.0]``
  
    **Predicates**
      - ``len(value) == 3``
  
   :B: Second boundary vector for plot. 
  
    **Type** ``List[float]``
  
    **Default** ``[0.0, 1.0, 0.0]``
  
    **Predicates**
      - ``len(value) == 3``
  
   :C: Third boundary vector for plot. 
  
    **Type** ``List[float]``
  
    **Default** ``[0.0, 0.0, 1.0]``
  
    **Predicates**
      - ``len(value) == 3``
  
 :MPI: Define MPI related parameters. 

  :red:`Keywords`
   :numerically_exact: This will use MPI algorithms that guarantees that the output is invariant wrt the number of MPI processes. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :shared_memory_size: Size (MB) of the MPI shared memory blocks of each shared function. 
  
    **Type** ``int``
  
    **Default** ``10000``
  
   :share_nuclear_potential: This will use MPI shared memory for the nuclear potential. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :share_coulomb_potential: This will use MPI shared memory for the Coulomb potential. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :share_xc_potential: This will use MPI shared memory for the exchange-correlation potential. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :bank_size: Number of MPI processes exclusively dedicated to manage orbital bank. 
  
    **Type** ``int``
  
    **Default** ``-1``
  
 :Basis: Define polynomial basis. 

  :red:`Keywords`
   :order: Polynomial order of multiwavelet basis. Negative value means it will be set automatically based on the world precision. 
  
    **Type** ``int``
  
    **Default** ``-1``
  
   :type: Polynomial type of multiwavelet basis. 
  
    **Type** ``str``
  
    **Default** ``interpolating``
  
    **Predicates**
      - ``value.lower() in ['interpolating', 'legendre']``
  
 :Derivatives: Define various derivative operators used in the code. 

  :red:`Keywords`
   :kinetic: Derivative used in kinetic operator. 
  
    **Type** ``str``
  
    **Default** ``abgv_55``
  
   :h_b_dip: Derivative used in magnetic dipole operator. 
  
    **Type** ``str``
  
    **Default** ``abgv_00``
  
   :h_m_pso: Derivative used in paramagnetic spin-orbit operator. 
  
    **Type** ``str``
  
    **Default** ``abgv_00``
  
 :Molecule: Define molecule. 

  :red:`Keywords`
   :charge: Total charge of molecule. 
  
    **Type** ``int``
  
    **Default** ``0``
  
   :multiplicity: Spin multiplicity of molecule. 
  
    **Type** ``int``
  
    **Default** ``1``
  
   :translate: Translate coordinates such that center of mass coincides with the global gauge origin. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :coords: Coordinates in xyz format. Atoms can be given either using atom symbol or atom number 
  
    **Type** ``str``
  
 :WaveFunction: Define the wavefunction method. 

  :red:`Keywords`
   :method: Wavefunction method. See predicates for valid methods. ``hf``, ``hartreefock`` and ``hartree-fock`` all mean the same thing, while ``lda`` is an alias for ``svwn5``. You can set a non-standard DFT functional (e.g. varying the amount of exact exchange) by choosing ``dft`` and specifing the functional(s) in the ``DFT`` section below. 
  
    **Type** ``str``
  
    **Predicates**
      - ``value.lower() in ['core', 'hartree', 'hf', 'hartreefock', 'hartree-fock', 'dft', 'lda', 'svwn3', 'svwn5', 'pbe', 'pbe0', 'bpw91', 'bp86', 'b3p86', 'b3p86-g', 'blyp', 'b3lyp', 'b3lyp-g', 'olyp', 'kt1', 'kt2', 'kt3']``
  
   :restricted: Use spin restricted wavefunction. 
  
    **Type** ``bool``
  
    **Default** ``True``
  
 :DFT: Define the exchange-correlation functional in case of DFT. 

  :red:`Keywords`
   :density_cutoff: Hard cutoff for passing density values to XCFun. 
  
    **Type** ``float``
  
    **Default** ``0.0``
  
   :functionals: List of density functionals with numerical coefficient. E.g. for PBE0 ``EXX 0.25``, ``PBEX 0.75``, ``PBEC 1.0``, see XCFun documentation <https://xcfun.readthedocs.io/>_. 
  
    **Type** ``str``
  
    **Default** `` ``
  
   :spin: Use spin separated density functionals. 
  
    **Type** ``bool``
  
    **Default** ``not(user['WaveFunction']['restricted'])``
  
 :Properties: Provide a list of properties to compute (total SCF energy and orbital energies are always computed). 

  :red:`Keywords`
   :dipole_moment: Compute dipole moment. 
  
    **Type** ``bool``
  
    **Default** ``True``
  
   :quadrupole_moment: Compute quadrupole moment. Note: Gauge origin dependent, should be used with ``translate = true`` in Molecule. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :polarizability: Compute polarizability tensor. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :magnetizability: Compute magnetizability tensor. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :nmr_shielding: Compute NMR shielding tensor. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :geometric_derivative: Compute geometric derivative. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :plot_density: Plot converged electron density. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :plot_orbitals: Plot converged molecular orbitals from list of indices, negative index plots all orbitals. 
  
    **Type** ``List[int]``
  
    **Default** ``[]``
  
 :ExternalFields: Define external electromagnetic fields. 

  :red:`Keywords`
   :electric_field: Strength of external electric field. 
  
    **Type** ``List[float]``
  
    **Default** ``[]``
  
    **Predicates**
      - ``len(value) == 0 or len(value) == 3``
  
 :Polarizability: Give details regarding the polarizability calculation. 

  :red:`Keywords`
   :frequency: List of external field frequencies. 
  
    **Type** ``List[float]``
  
    **Default** ``[0.0]``
  
 :NMRShielding: Give details regarding the NMR shileding calculation. 

  :red:`Keywords`
   :nuclear_specific: Use nuclear specific perturbation operator (h_m_pso). 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :nucleus_k: List of nuclei to compute. Negative value computes all nuclei. 
  
    **Type** ``List[int]``
  
    **Default** ``[-1]``
  
 :Files: Defines file paths used for program input/output. Note: all paths must be given in quotes if they contain slashes "path/to/file". 

  :red:`Keywords`
   :guess_basis: File name for GTO basis set, used with ``gto`` guess. 
  
    **Type** ``str``
  
    **Default** ``initial_guess/mrchem.bas``
  
   :guess_gto_p: File name for paired orbitals, used with ``gto`` guess. 
  
    **Type** ``str``
  
    **Default** ``initial_guess/mrchem.mop``
  
   :guess_gto_a: File name for alpha orbitals, used with ``gto`` guess. 
  
    **Type** ``str``
  
    **Default** ``initial_guess/mrchem.moa``
  
   :guess_gto_b: File name for beta orbitals, used with ``gto`` guess. 
  
    **Type** ``str``
  
    **Default** ``initial_guess/mrchem.mob``
  
   :guess_phi_p: File name for paired orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/phi_p_scf_idx_<0...Np>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_p``
  
   :guess_phi_a: File name for alpha orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/phi_a_scf_idx_<0...Na>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_a``
  
   :guess_phi_b: File name for beta orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/phi_b_scf_idx_<0...Nb>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_b``
  
   :guess_x_p: File name for paired response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/x_p_rsp_idx_<0...Np>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/X_p``
  
   :guess_x_a: File name for alpha response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/x_a_rsp_idx_<0...Na>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/X_a``
  
   :guess_x_b: File name for beta response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/x_b_rsp_idx_<0...Nb>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/X_b``
  
   :guess_y_p: File name for paired response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/y_p_rsp_idx_<0...Np>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/Y_p``
  
   :guess_y_a: File name for alpha response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/y_a_rsp_idx_<0...Na>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/Y_a``
  
   :guess_y_b: File name for beta response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/y_b_rsp_idx_<0...Nb>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/Y_b``
  
   :guess_cube_p: File name for paired orbitals, used with ``cube`` guess. Expected path is ``<path_orbitals>/phi_p_scf_idx_<0...Np>_<re/im>.cube 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_p``
  
   :guess_cube_a: File name for alpha orbitals, used with ``cube`` guess. Expected path is ``<path_orbitals>/phi_a>_scf_idx_<0...Na>_<re/im>.cube 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_a``
  
   :guess_cube_b: File name for beta orbitals, used with ``cube`` guess. Expected path is ``<path_orbitals>/phi_b_scf_idx_<0...Nb>_<re/im>.cube 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_b``
  
   :cube_vectors: Directory where cube vectors are stored for mrchem calculation. 
  
    **Type** ``str``
  
    **Default** ``cube_vectors/``
  
 :SCF: Includes parameters related to the ground state SCF orbital optimization. 

  :red:`Keywords`
   :run: Run SCF solver. Otherwise properties are computed on the initial orbitals. 
  
    **Type** ``bool``
  
    **Default** ``True``
  
   :max_iter: Maximum number of SCF iterations. 
  
    **Type** ``int``
  
    **Default** ``100``
  
   :kain: Length of KAIN iterative history. 
  
    **Type** ``int``
  
    **Default** ``5``
  
   :rotation: Number of iterations between each diagonalization/localization. 
  
    **Type** ``int``
  
    **Default** ``0``
  
   :localize: Use canonical or localized orbitals. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :energy_thrs: Convergence threshold for SCF energy. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :guess_prec: Precision parameter used in construction of initial guess. 
  
    **Type** ``float``
  
    **Default** ``0.001``
  
    **Predicates**
      - ``1.0e-10 < value < 1.0``
  
   :guess_screen: Screening parameter used in GTO evaluations, in number of standard deviations. Every coordinate beyond N StdDev from the Gaussian center is evaluated to zero. Note that too aggressive screening is counter productive, because it leads to a sharp cutoff in the resulting function which requires higher grid refinement. Negative value means no screening. 
  
    **Type** ``float``
  
    **Default** ``12.0``
  
   :start_prec: Incremental precision in SCF iterations, initial value. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :final_prec: Incremental precision in SCF iterations, final value. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :guess_type: Type of initial guess for ground state orbitals. ``chk`` restarts a previous calculation which was dumped using the ``write_checkpoint`` keyword. This will load MRA and electron spin configuration directly from the checkpoint files, which are thus required to be identical in the two calculations. ``mw`` will start from final orbitals in a previous calculation written using the ``write_orbitals`` keyword. The orbitals will be re-projected into the new computational setup, which means that the electron spin configuration and MRA can be different in the two calculations. ``gto`` reads precomputed GTO orbitals (requires extra non-standard input files for basis set and MO coefficients). ``core`` and ``sad`` will diagonalize the Fock matrix in the given AO basis (SZ, DZ, TZ or QZ) using a Core or Superposition of Atomic Densities Hamiltonian, respectively. 
  
    **Type** ``str``
  
    **Default** ``sad_dz``
  
    **Predicates**
      - ``value.lower() in ['mw', 'chk', 'gto', 'core_sz', 'core_dz', 'core_tz', 'core_qz', 'sad_sz', 'sad_dz', 'sad_tz', 'sad_qz', 'sad_gto', 'cube']``
  
   :write_checkpoint: Write orbitals to disk in each iteration, file name ``<path_checkpoint>/phi_scf_idx_<0..N>``. Can be used as ``chk`` initial guess in subsequent calculations. Note: must be given in quotes if there are slashes in the path "path/to/checkpoint". 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :path_checkpoint: Path to checkpoint files during SCF, used with ``write_checkpoint`` and ``chk`` guess. 
  
    **Type** ``str``
  
    **Default** ``checkpoint``
  
    **Predicates**
      - ``value[-1] != '/'``
  
   :write_orbitals: Write final orbitals to disk, file name ``<path_orbitals>/phi_<p/a/b>_scf_idx_<0..Np/Na/Nb>``. Can be used as ``mw`` initial guess in subsequent calculations. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :path_orbitals: Path to where converged orbitals will be written in connection with the ``write_orbitals`` keyword. Note: must be given in quotes if there are slashes in the path "path/to/orbitals". 
  
    **Type** ``str``
  
    **Default** ``orbitals``
  
    **Predicates**
      - ``value[-1] != '/'``
  
   :orbital_thrs: Convergence threshold for orbital residuals. 
  
    **Type** ``float``
  
    **Default** ``10 * user['world_prec']``
  
 :Response: Includes parameters related to the response SCF optimization. 

  :red:`Keywords`
   :run: In which Cartesian directions to run response solver. 
  
    **Type** ``List[bool]``
  
    **Default** ``[True, True, True]``
  
   :max_iter: Maximum number of response iterations. 
  
    **Type** ``int``
  
    **Default** ``100``
  
   :kain: Length of KAIN iterative history. 
  
    **Type** ``int``
  
    **Default** ``5``
  
   :property_thrs: Convergence threshold for symmetric property. Symmetric meaning the property computed from the same operator as the response purturbation, e.g. for external magnetic field the symmetric property corresponds to the magnetizability (NMR shielding in non-symmetric, since one of the operators is external magnetic field, while the other is nuclear magnetic moment). 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :start_prec: Incremental precision in SCF iterations, initial value. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :final_prec: Incremental precision in SCF iterations, final value. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :guess_prec: Precision parameter used in construction of initial guess. 
  
    **Type** ``float``
  
    **Default** ``0.001``
  
    **Predicates**
      - ``1.0e-10 < value < 1.0``
  
   :guess_type: Type of initial guess for response. ``none`` will start from a zero guess for the response functions. ``chk`` restarts a previous calculation which was dumped using the ``write_checkpoint`` keyword. ``mw`` will start from final orbitals in a previous calculation written using the ``write_orbitals`` keyword. The orbitals will be re-projected into the new computational setup. 
  
    **Type** ``str``
  
    **Default** ``none``
  
    **Predicates**
      - ``value.lower() in ['none', 'chk', 'mw']``
  
   :write_checkpoint: Write perturbed orbitals to disk in each iteration, file name ``<path_checkpoint>/<X/Y>_rsp_<direction>_idx_<0..N>``. Can be used as ``chk`` initial guess in subsequent calculations. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :path_checkpoint: Path to checkpoint files during SCF, used with ``write_checkpoint`` and ``chk`` guess. 
  
    **Type** ``str``
  
    **Default** ``checkpoint``
  
    **Predicates**
      - ``value[-1] != '/'``
  
   :write_orbitals: Write final perturbed orbitals to disk, file name ``<path_orbitals>/<X/Y>_<p/a/b>_rsp_<direction>_idx_<0..Np/Na/Nb>``. Can be used as ``mw`` initial guess in subsequent calculations. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :path_orbitals: Path to where converged orbitals will be written in connection with the ``write_orbitals`` keyword. 
  
    **Type** ``str``
  
    **Default** ``orbitals``
  
    **Predicates**
      - ``value[-1] != '/'``
  
   :orbital_thrs: Convergence threshold for orbital residuals. 
  
    **Type** ``float``
  
    **Default** ``10 * user['world_prec']``
  
   :localize: Use canonical or localized unperturbed orbitals. 
  
    **Type** ``bool``
  
    **Default** ``user['SCF']['localize']``
  
 :Environment: Includes parameters related to the computation of the reaction field energy of a system in an environment. 

  :red:`Keywords`
   :max_iter: Max number of iterations allowed in the nested procedure. 
  
    **Type** ``int``
  
    **Default** ``100``
  
   :run_environment: Perform the reaction field calculation of the reaction potential of the interaction between environment and molecule.  
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :algorithm: What algorithm to use for the reaction field ``scrf`` runs a nested algorithm where the generalized Poisson equation is solved iterativelly until self consistency wrt. the convergence threshold. 
  
    **Type** ``str``
  
    **Default** ``scrf``
  
    **Predicates**
      - ``value.lower() in ['scrf']``
  
   :convergence_criterion: Adjust the convergence threshold for the nested procedure. ``dynamic`` Uses the absolute value of the latest orbital update as convergence threshold. When the orbitals are close to convergence (``mo_residual < world_prec*10``) the convergence threshold will be equal to ``world_prec``. ``static`` uses ``world_prec`` as convergence threshold. 
  
    **Type** ``str``
  
    **Default** ``dynamic``
  
    **Predicates**
      - ``value.lower() in ['dynamic', 'static']``
  
   :extrapolate_Vr: Extrapolate on the reaction potential if true, or on the surface charge distribution in the convergence acceleration. 
  
    **Type** ``bool``
  
    **Default** ``True``
  
   :density_type: What part of the total molecular charge density to use in the algorithm. ``total`` uses the total charge density. ``nuclear`` uses only the nuclear part of the total charge density. ``electronic`` uses only the electronic part of the total charge density. 
  
    **Type** ``str``
  
    **Default** ``total``
  
    **Predicates**
      - ``value.lower() in ['total', 'nuclear', 'electronic']``
  
   :kain: Number of previous reaction field iterates kept for convergence acceleration during the nested precedure. 
  
    **Type** ``int``
  
    **Default** ``user['SCF']['kain']``
  
  :red:`Sections`
   :Cavity: Define the interlocking spheres cavity. 
  
      :red:`Keywords`
       :spheres: Coordinates and radii  of the spheres written as $spheres x_0    y_0    z_0    R_0 ... x_N    y_N    z_N    R_N $end The units used are the same specified with the `world_unit` keyword. 
      
        **Type** ``str``
      
        **Default** ````
      
       :cavity_width: Width of cavity boundary 
      
        **Type** ``float``
      
        **Default** ``0.2``
      
   :Permittivity: Parameters for the permittivity function. 
  
      :red:`Keywords`
       :epsilon_in: Permittivity inside the cavity. 1.0 is the permittivity of free space, anything other than this is undefined behaviour. 
      
        **Type** ``float``
      
        **Default** ``1.0``
      
       :epsilon_out: Permittivity outside the cavity. This is characteristic of the solvent used. 
      
        **Type** ``float``
      
        **Default** ``2.0``
      
       :formulation: Formulation of the Permittivity function. Currently only the exponential is used. 
      
        **Type** ``str``
      
        **Default** ``exponential``
      
        **Predicates**
          - ``value.lower() in ['exponential']``
      ============
Installation
============

-------------------
Build prerequisites
-------------------

- Python-3.6 (or later)
- CMake-3.12 (or later)
- GNU-5.4 or Intel-17 (or later) compilers (C++14 standard)

.. hint::
    We have collected the recommended modules for the different Norwegian HPC
    systems under ``tools/<machine>.env``. These files can be sourced in order
    to get a working environment on the respective machines, and may also serve
    as a guide for other HPC systems.


C++ dependencies
----------------

The MRChem program depends on the following C++ libraries:

- Input handling: `nlohmann/json-3.6  <https://github.com/nlohmann/json>`_
- Multiwavelets: `MRCPP-1.3  <https://github.com/MRChemSoft/mrcpp>`_
- Linear algebra: `Eigen-3.3  <https://gitlab.com/libeigen/eigen>`_
- DFT functionals: `XCFun-2.0  <https://github.com/dftlibs/xcfun>`_

All these dependencies will be downloaded automatically at configure time by
CMake, but can also be linked manually by setting the variables::

    MRCPP_DIR=<path_to_mrcpp>/share/cmake/MRCPP
    XCFun_DIR=<path_to_xcfun>/share/cmake/XCFun
    Eigen3_DIR=<path_to_eigen3>/share/eigen3/cmake
    nlohmann_json_DIR=<path_to_nlohmann_json>


Python dependencies
-------------------

**Users** only need a Python3 interpreter, which is used for configuration
(``setup`` script) as well as launching the program (``mrchem`` script).

**Developers** will need some extra Python packages to update the input
parser and build the documentation locally with Sphinx.

We **strongly** suggest not to install these Python dependencies globally, but
rather to use a local virtual environment. We provide a ``Pipfile`` for
specifying the Python dependencies.
We recommend using `Pipenv <https://pipenv.readthedocs.io/en/latest/>`_, since
it manages virtual environment and package installation seamlessly.
After installing it with your package manager, run::

    $ pipenv install --dev

to create a virtual environment with all developer packages installed.

The environment can be activated with::

    $ pipenv shell

Alternatively, any Python command can be run within the virtual environment by
doing::

    $ pipenv run python -c "print('Hello, world')"


-------------------------------
Obtaining and building the code
-------------------------------

The latest development version of MRChem can be found on the ``master``
branch on GitHub::

    $ git clone https://github.com/MRChemSoft/mrchem.git

The released versions can be found from Git tags ``vX.Y.Z`` under the
``release/X.Y`` branches in the same repository, or a zip file can be
downloaded from `Zenodo <https://doi.org/10.5281/zenodo.3606658>`_.

By default, all dependencies will be **fetched** at configure time if they are
not already available.


Configure
---------

The ``setup`` script will create a directory called ``<build-dir>`` and run
CMake. There are several options available for the setup, the most
important being:

``--cxx=<CXX>``
  C++ compiler [default: g++]
``--omp``
  Enable OpenMP parallelization [default: False]
``--mpi``
  Enable MPI parallelization [default: False]
``--type=<TYPE>``
  Set the CMake build type (debug, release, relwithdebinfo, minsizerel) [default: release]
``--prefix=<PATH>``
  Set the install path for make install [default: '/usr/local']
``--cmake-options=<STRING>``
  Define options to CMake [default: '']
``-h --help``
  List all options

The code can be built with four levels of parallelization:

 - no parallelization
 - only shared memory (OpenMP)
 - only distributed memory (MPI)
 - hybrid OpenMP + MPI

.. note::
    In practice we recommend the **shared memory version** for running on your
    personal laptop/workstation, and the **hybrid version** for running on a
    HPC cluster. The serial and pure MPI versions are only useful for debugging.

The default build is *without* parallelization and using GNU compilers::

    $ ./setup --prefix=<install-dir> <build-dir>

To use Intel compilers you need to specify the ``--cxx`` option::

    $ ./setup --prefix=<install-dir> --cxx=icpc <build-dir>

To build the code with shared memory (OpenMP) parallelization,
add the ``--omp`` option::

    $ ./setup --prefix=<install-dir> --omp <build-dir>

To build the code with distributed memory (MPI) parallelization, add the
``--mpi`` option *and* change to the respective MPI compilers (``--cxx=mpicxx``
for GNU and ``--cxx=mpiicpc`` for Intel)::

    $ ./setup --prefix=<install-dir> --omp --mpi --cxx=mpicxx <build-dir>

When dependencies are fetched at configuration time, they will be downloaded
into ``<build-dir>/_deps``. For the example of MRCPP, sources are saved into
the folders ``<build-dir>/_deps/mrcpp_sources-src`` and built into
``<build-dir>/_deps/mrcpp_sources-build``.

.. note::
    If you compile the MRCPP library manually as a separate project, the level
    of parallelization **must be the same** for MRCPP and MRChem. Similar
    options apply for the MRCPP setup, see
    `mrcpp.readthedocs.io <https://mrcpp.readthedocs.io/en/latest/>`_.


Build
-----

If the CMake configuration is successful, the code is compiled with::

    $ cd <build-dir>
    $ make


Test
----

A test suite is provided to make sure that everything compiled properly.
To run a collection of small unit tests::

    $ cd <build-dir>
    $ ctest -L unit

To run a couple of more involved integration tests::

    $ cd <build-dir>
    $ ctest -L integration


Install
-------

After the build has been verified with the test suite, it can be installed with
the following command::

    $ cd <build-dir>
    $ make install

This will install *two* executables under the ``<install-path>``::

    <install-path>/bin/mrchem       # Python input parser and launcher
    <install-path>/bin/mrchem.x     # MRChem executable

Please refer to the :ref:`User's Manual` for instructions for how to run the program.

.. hint::
    We have collected scripts for configure and build of the hybrid OpenMP + MPI
    version on the different Norwegian HPC systems under ``tools/<machine>.sh``.
    These scripts will build the current version under ``build-${version}``,
    run the unit tests and install under ``install-${version}``, e.g. to build
    version v1.0.0 on Fram::

        $ cd mrchem
        $ git checkout v1.0.0
        $ tools/fram.sh

    The configure step requires internet access, so the scripts must be run on
    the login nodes, and it will run on a single core, so it might take some
    minutes to complete. The scripts will *not* install the :ref:`Python
    dependencies`, so this must be done manually in order to run the code.

.. MRChem documentation master file, created by
   sphinx-quickstart on Tue Jan 26 15:03:29 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==================================
Welcome to MRChem's documentation!
==================================

MRChem is a numerical real-space code for molecular
electronic structure calculations within the self-consistent field (SCF)
approximations of quantum chemistry (Hartree-Fock and Density Functional
Theory). The code is divided in two main parts: the `MultiResolution
Computation Program Package <https://mrcpp.readthedocs.io/en/latest>`_ (MRCPP),
which is a general purpose numerical
mathematics library based on multiresolution analysis and the multiwavelet
basis which provide low-scaling algorithms as well as rigorous error control
in numerical computations, and the MultiResolution Chemistry (MRChem) program
that uses the functionalities of MRCPP for computational chemistry applications.

The code is being developed at the `Hylleraas Centre for Quantum Molecular
Sciences <http://www.ctcc.no/>`_ at
`UiT - The Arctic University of Norway <http://en.uit.no>`_.

--------------------------------------------------------------------------------

The code is under active development, and the latest stable releases as well as
development versions can be found on `GitHub <https://github.com/MRChemSoft/mrchem>`_.

Features in MRChem-1.0.0:
-------------------------

* Wave functions:
    + Kohn-Sham DFT
        - Spin-polarized
        - Spin-unpolarized
        - LDA, GGA and hybrid functionals
    + Hartree-Fock
        - Restricted closed-shell
        - Unrestricted
    + Explicit external fields
        - Electric field
* Properties:
    + Ground state energy
    + Dipole moment
    + Quadrupole moment
    + Polarizability
    + Magnetizability
    + NMR shielding constant
    + Geometric derivative
* Parallel implementation:
    + Shared memory (OpenMP): ~20 cores
    + Distributed memory (MPI): ~100 procs
    + Hybrid scheme (MPI + OpenMP): ~1000 cores
* Current size limitations:
    + ~200 orbitals on ~50 medium-memory (64GB) compute nodes
    + ~100 orbitals on a single high-memory (1TB) compute node

Upcoming features:
------------------

* Wave functions:
    + Meta-GGAs
    + ZORA Hamiltonian
    + Solvent effects
    + Periodic Boundary Conditions
    + External magnetic field
* Properties:
    + Optical rotation
    + Spin-spin coupling constant
    + Hyperfine coupling constant
    + Magnetically induced currents
    + Hyperpolarizability
    + Geometry optimization
* Performance:
    + Improved parallel scalability
    + Improved exact exchange performance
    + More efficient memory distribution


.. toctree::
   :maxdepth: 2

   installation
   users/manual
   programmers/manual
Documentation
=============

This documentation is generated using `Sphinx <http://sphinx-doc.org/>`_ and
`Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_ The two softwares are
bridged by means of the `Breathe extension <https://breathe.readthedocs.org/>`_
The online version of this documentation is built and served by `Read The Docs
<https://readthedocs.org/>`_.  The webpage http://mrchem.readthedocs.io/ is
updated on each push to the public GitHub repository.


How and what to document
------------------------

Doxygen enables documenting the code in the source code files thus removing a
"barrier" for developers.  To avoid that the code degenerates into a Big Ball
of Mud, it is mandatory to document directly within the source code classes and
functions.  To document general programming principles, design choices,
maintenance etc. you can create a .rst file in the ``doc`` directory. Remember
to refer the new file inside the ``index.rst`` file (it won't be parsed
otherwise).  Sphing uses `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ and `Markdown
<https://daringfireball.net/projects/markdown/>`_. Support for Markdown is not
as extensive as for reStructuredText, see `these comments
<https://blog.readthedocs.com/adding-markdown-support/>`_ Follow the guidelines
in :cite:`Wilson2014` regarding what to document.

Documeting methods in derived classes
-------------------------------------

Virtual methods should only be documented in the base classes.
This avoids unnecessary verbosity and conforms to the principle: "Document
_what_, not _how_" :cite:`Wilson2014`
If you feel the _how_ needs to be explicitly documented, add some notes in the
appropriate ``.rst`` file.

How does this work?
-------------------

To have an offline version of the documentation just issue ``make doc`` in the
build directory.  The HTML will be stored in ``doc/html``. Open the
``doc/html/index.html`` file with your browser to see and browse the
documentation.

.. warning::

   Building the documentation requires Doxygen, Sphinx, Perl and the Python
   modules PyYAML, Breathe and Matplotlib.
CMake usage
===========

This is a brief guide to our CMake infrastructure which is managed
*via* `Autocmake <http://autocmake.readthedocs.org/en/latest/>`_

.. warning::

   The minimum required CMake version is 2.8.8

Adding new source subdirectories and/or files
---------------------------------------------

Developers **HAVE TO** manually list the sources in a given subdirectory
of the main source directory ``src/``.

New subdirectory
................

First of all, you will have to let CMake know that a new source-containing
subdirectory has been added to the source tree. Due to the hierarchical
approach CMake is based upon, you will need to modify the ``CMakeLists.txt`` in
the ``src`` directory and create a new one in your new subdirectory.  For the
first step:

   1. if your new subdirectory contains header files, add a line like
   the following to the ``CMakeLists.txt`` file contained in the ``src`` directory:

   .. code-block:: cmake

      ${CMAKE_CURRENT_LIST_DIR}/subdir_name

   to the command setting the list of directories containing headers.  This
   sets up the list of directories where CMake will look for headers with
   definitions of classes and functions. If your directory contains Fortran
   code you can skip this step;

   2. add a line like the following to the ``CMakeLists.txt`` file contained in the
   ``src`` directory:

   .. code-block:: cmake

      add_subdirectory(subdir_name)

   This will tell CMake to go look inside ``subdir_name`` for a ``CMakeLists.txt``
   containing more sets of instructions.  It is preferable to add these new
   lines in **alphabetic order**

Inside your new subdirectory you will need to add a ``CMakeLists.txt`` file containing
the set of instructions to build your cutting edge code. This is the second step.
Run the ``make_cmake_files.py`` Python script in the ``src/`` directory:

.. code-block:: bash

   python make_cmake_files.py --libname=cavity --lang=CXX

to generate a template ``CMakeLists.txt.try`` file:

.. code-block:: cmake

   # List of headers
   list(APPEND headers_list Cavity.hpp ICavity.hpp Element.hpp GePolCavity.hpp RegisterCavityToFactory.hpp RestartCavity.hpp)

   # List of sources
   list(APPEND sources_list ICavity.cpp Element.cpp GePolCavity.cpp RestartCavity.cpp)

   add_library(cavity OBJECT ${sources_list} ${headers_list})
   set_target_properties(cavity PROPERTIES POSITION_INDEPENDENT_CODE 1 )
   set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
   # Sets install directory for all the headers in the list
   foreach(_header ${headers_list})
      install(FILES ${_header} DESTINATION include/cavity)
   endforeach()

The template might need additional editing.
Each source subdirectory is the lowest possible in the CMake
hierarchy and it contains set of instructions for:

#. exporting a list of header files (.h or .hpp) to the upper level in the
   hierarchy, possibly excluding some of them
#. define install targets for the files in this subdirectory.

All the source files are compiled into the unique static library ``libpcm.a`` and unique
dynamic library ``libpcm.so``.
This library is the one the host QM program need to link.

Searching for libraries
.......................

In general, the use of the `find_package <http://www.cmake.org/cmake/help/v3.0/command/find_package.html>`_
macro is to be preferred, as it is standardized and ensured to work on any
platform.  Use of ``find_package`` requires that the package/library you want to
use has already a module inside the CMake distribution.  If that's not the
case, you should *never* use the following construct for third-party libraries:

.. code-block:: cmake

   target_link_libraries(myexe -lsomesystemlib)

If the library does not exist, the end result is a cryptic linker error. See
also `Jussi Pakkanen's blog <http://voices.canonical.com/jussi.pakkanen/2013/03/26/a-list-of-common-cmake-antipatterns/>`_
You will first need to find the library, using the macro
`find_library <http://www.cmake.org/cmake/help/v3.0/command/find_library.html>`_,
and then use the ``target_link_libraries`` command.
Maintenance
===========

Description and how-to for maintenance operations.
Some of the maintenance scripts have been moved to the `pcmsolvermeta
repository <https://gitlab.com/PCMSolver/pcmsolvermeta>`_

Branching Model and Release Process
-----------------------------------

.. warning::
   **Incomplete or outdated information!**

Releases in a ``X.Y.Z`` series are annotated tags on the corresponding branch.

Pull Request Requirements
-------------------------

The project is integrated with `Danger.Systems <http://danger.systems/ruby/>`_.
On each PR, one CI job will run the integration and a `bot <https://github.com/minazobot>`_ will
report which requirements are **not met** in your PR.
These reports can be _warnings_ and _errors_. You will discuss and solve both
of them with the reviewers.
The automatic rules are laid out in the ``Dangerfile`` and are used to enforce an
adequate level of testing, documentation and code quality.

Danger.Systems Warnings
=======================

- PRs classed as Work in Progress.
- Codebase was modified, but no tests were added.
- Nontrivial changes to the codebase, but no documentation added.
- Codebase was modified, but ``CHANGELOG.md`` was not updated.
- Source files were added or removed, but ``.gitattributes`` was not updated.

Danger.Systems Errors
=====================

- Commit message linting, based on some of `these recommendations <https://chris.beams.io/posts/git-commit/>`_:
  - Commit subject is more than one word.
  - Commit subject is no longer than 50 characters.
  - Commit subject and body are separated by an empty line.

- Clean commit history, without merge commits.

- Code style for ``.hpp``, ``.cpp``, ``.h`` files follows the conventions in
  ``.clang-format``.

Bump Version
------------

Version numbering follows the guidelines of `semantic versioning <http://semver.org/>`_
To update, change the relevant field in the ``README.md`` file.

Changelog
---------

We follow the guidelines of `Keep a CHANGELOG <http://keepachangelog.com/>`_
On all **but** the release branches, there is an ``Unreleased`` section
under which new additions should be listed.
To simplify perusal of the ``CHANGELOG.md``, use the following subsections:

1. ``Added`` for new features.
2. ``Changed`` for changes in existing functionality.
3. ``Deprecated`` for once-stable features removed in upcoming releases.
4. ``Removed`` for deprecated features removed in this release.
5. ``Fixed`` for any bug fixes.
6. ``Security`` to invite users to upgrade in case of vulnerabilities.

Updating Eigen Distribution
---------------------------

The C++ linear algebra library Eigen comes bundled with the module. To update
the distributed version one has to:

1. download the desired version of the library to a scratch location. Eigen's
   website is: http://eigen.tuxfamily.org/
2. unpack the downloaded archive;
3. go into the newly created directory and create a build directory;
4. go into the newly created build directory and type the following (remember
   to substitute @PROJECT_SOURCE_DIR@ with the actual path)

   .. code-block:: bash

    cmake .. -DCMAKE_INSTALL_PREFIX=@PROJECT_SOURCE_DIR@/external/eigen3

Remember to commit and push your modifications.

Git Pre-Commit Hooks
--------------------

[Git pre-commit hooks](https://git-scm.com/book/gr/v2/Customizing-Git-Git-Hooks) are used to
keep track of code style and license header in source files.
Code style is checked using ``clang-format``.

.. warning::
   **You need to install ``clang-format`` (v3.9 recommended) to run the code style validation hook!**

License headers are checked using the ``license_maintainer.py`` script and the
header templates for the different languages used in this project.
The Python script checks the ``.gitattributes`` file to determine which license
headers need to be maintained and in which files:

.. code-block:: bash

   src/pedra/pedra_dlapack.F90 !licensefile
   src/solver/*.hpp licensefile=.githooks/LICENSE-C++

The first line specifies that the file in ``src/pedra/pedra_dlapack.F90`` should
not be touched, while the second line states that all ``.hpp`` files in ``src/solver``
should get an header from the template in ``.githooks/LICENSE-C++``
Location of files in ``.gitattributes`` are always specified with respect
to the project root directory.

The hooks are located in the ``.githooks`` subdirectory and **have to be installed by hand**
whenever you clone the repository anew:

.. code-block:: bash

   cd .git/hooks
   cp --symbolic-link ../../.githooks/* .

Installed hooks will **always** be executed. Use ``git commit --no-verify`` to
bypass explicitly the hooks.
Coding standards
================

General Object-Oriented design principles you should try to follow:
  1. Identify the aspects of your application that vary and separate them from what stays the same;
  2. Program to an interface, not an implementation;
  3. Favor composition over inheritance;
  4. Strive for loosely coupled designs between objects that interact;
  5. Classes should be open for extension, but closed for modification;
  6. Depend upon abstractions. Do not depend upon concrete classes;
  7. Principle of Least Knowledge. Talk only to your immediate friends;

:cite:`Sutter2004,Cline1998,CppFAQs`

Including header files
----------------------

Do not include header files unnecessarily. Unnecessary include
directives and/or forward declarations introduce nasty
interdependencies among different parts of the code.  This reflects
mainly in longer compilation times, but also in uglier looking code
(see also the discussion in :cite:`Sutter1999`).

Follow these guidelines to decide whether to include or forward declare:
  1. class A makes no reference to class B. Neither include nor forward declare B;
  2. class A refers to class B as a friend. Neither include nor forward declare B;
  3. class A contains a pointer/reference to a class B object. Forward declare B;
  4. class A contains functions with a class B object (value/pointer/reference) as parameter/return value. Forward declare B;
  5. class A is derived from class B. include B;
  6. class A contains a class B object. include B.

.. code-block:: cpp

    #pragma once

    //==============================
    // Forward declared dependencies
    class Foo;
    class Bar;

    //==============================
    // Included dependencies
    #include <vector>
    #include "Parent.hpp"

    //==============================
    // The actual class
    class MyClass : public Parent // Parent object, so #include "Parent.h"
    {
      public:
        std::vector<int> avector; // vector object, so #include <vector>
        Foo * foo;                // Foo pointer, so forward declare
        void Func(Bar & bar);     // Bar reference as parameter, so forward declare

        friend class MyFriend;    // friend declaration is not a dependency
                                  //    don't do anything about MyFriend
    };


Proper overloading of `operator<<`
----------------------------------

Suppose we have an inheritance hierarchy made of an abstract base class, Base, and
two derived classes, Derived1 and Derived2.
In the Base class header file we will define a pure virtual private function printObject
and provide a public friend overload of operator<<:

.. code-block:: cpp

    #include <iosfwd>

    class Base
    {
      public:
        // All your other very fancy public members
        friend std::ostream & operator<<(std::ostream & os, Base & base)
        {
                return base.printObject(os);
        }
      protected:
        // All your other very fancy protected members
      private:
        // All your other very fancy private members
        virtual std::ostream & printObject(std::ostream & os) = 0;
    }

The printObject method can also be made (impure) virtual, it really depends on your class hierarchy.
Derived1 and Derived2 header files will provide a public friend overload of operator<< (friendliness
isn't inherited, transitive or reciprocal) and an override for the printObject method:

.. code-block:: cpp

    #include <iosfwd>

    #include "Base.hpp"

    class Derived1 : public Base
    {
      public:
        // All your other very fancy public members
        friend std::ostream & operator<<(std::ostream & os, Derived1 & derived)
        {
          return derived.printObject(os);
        }
      protected:
        // All your other very fancy protected members
      private:
        // All your other very fancy private members
        virtual std::ostream & printObject(std::ostream & os);
    }

    class Derived2 : public Base
    {
      public:
        // All your other very fancy public members
        friend std::ostream & operator<<(std::ostream & os, Derived2 & derived)
        {
          return derived.printObject(os);
        }
      protected:
        // All your other very fancy protected members
      private:
        // All your other very fancy private members
        virtual std::ostream & printObject(std::ostream & os);
    }

Code formatting
---------------

We conform to the so-called Linux (aka kernel) formatting style for C/C++ code
(see http://en.wikipedia.org/wiki/Indent_style#Kernel_style) with minimal
modifications.
Using `clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_ is the
preferred method to get the source code in the right format.
Formatting style is defined in the ``.clang-format`` file, kept at the root of the project.

.. note::
   We recommend using at least v3.9 of the program, which is the version used to
   generate the ``.clang-format`` file defining all formatting settings.

``clang-format`` can be `integrated with both
Emacs and Vim. <https://clang.llvm.org/docs/ClangFormat.html#vim-integration>`_
It is also possible to install the Git pre-commit hooks to perform the necessary code style
checks prior to committing changes:

.. code-block:: bash

   cd .git/hooks
   cp --symbolic-link ../../.githooks/* .
Testing
-------

We perform unit testing of our code. The unit testing framework used is
`Catch <https://github.com/philsquared/Catch>`_ The framework provides quite an
extensive set of macros to test various data types, it also provides facilities
for easily setting up test fixtures.  Usage is extremely simple and the
`documentation <https://github.com/philsquared/Catch/blob/master/docs/Readme.md>`_
is very well written.  For a quick primer on how to use Catch refer to:
https://github.com/philsquared/Catch/blob/master/docs/tutorial.md
The basic idea of unit testing is to test each building block of the code
separataly. In our case, the term "building block" is used to mean a class.

To add new tests for your class you have to:

#. create a new subdirectory inside tests/ and add a line like the following
   to the CMakeLists.txt

   .. code-block:: cmake

      add_subdirectory(new_subdir)

#. create a CMakeLists.txt inside your new subdirectory.
   This CMakeLists.txt adds the source for a given unit test to the global ``UnitTestsSources``
   property and notifies CTest that a test with given name is part of the test suite.
   The generation of the CMakeLists.txt can be managed by ``make_cmake_files.py`` Python script.
   This will take care of also setting up CTest labels. This helps in further grouping
   the tests for our convenience.
   Catch uses tags to index tests and tags are surrounded by square brackets. The Python script
   inspects the sources and extracts labels from Catch tags.
   The ``add_Catch_test`` CMake macro takes care of the rest.

   We require that each source file containing tests follows the naming convention
   new_subdir_testname and that testname gives some clue to what is being tested.
   Depending on the execution of tests in a different subdirectory is bad practice.
   A possible workaround is to add some kind of input file and create a text fixture
   that sets up the test environment. Have a look in the ``tests/input`` directory
   for an example
===================
Programmer's Manual
===================


Classes and functions reference
-------------------------------

.. toctree::

   code_reference/chemistry
   code_reference/environment
   code_reference/initial_guess
   code_reference/properties
   code_reference/qmfunctions
   code_reference/qmoperators
   code_reference/scf_solver
   
:orphan:

=============================
MRChem Programmers' Manual
=============================


.. toctree::

   general-structure
   coding-standards
   documentation
   cmake-usage
   maintenance
   testing
General Structure
=================

External libraries:

+ parts of the C++ `Boost <http://www.boost.org/>`_ libraries are used to provide
  various functionality, like timing and metaprogramming.
  The source for the 1.54.0 release is shipped with the
  module's source code. Some of the libraries used
  need to be compiled. Boost is released under the terms
  of the `Boost Software License, v1.0 <http://opensource.org/licenses/BSL-1.0>`_ (see also
  http://www.boost.org/users/license.html) We encourage the use of
  Boost whenever some functionality has already been coded within those
  libraries. However, consider **carefully** the introduction of functionality
  depending on compiler Boost libraries.
+ the `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ template
  library for linear algebra.  Almost every operation involving matrices and
  vectors is performed through Eigen.  Eigen provides convenient type
  definitions for vectors and matrices (of arbitrary dimensions) and the
  corresponding operations. Have a look
  `here <http://eigen.tuxfamily.org/dox/group__QuickRefPage.html>`_ for a quick
  reference guide to the API and
  at the `getting started guide <http://eigen.tuxfamily.org/dox/GettingStarted.html>`_ to get started.
  Eigen is distributed under the terms of the `Mozilla Public License, v2.0
  <http://opensource.org/licenses/MPL-2.0>`_
+ the `Getkw library <https://github.com/juselius/libgetkw>`_ by Jonas Juselius is
  used to manage input.  It is distributed under the terms of the `GNU General
  Public License, v2.0 <http://opensource.org/licenses/GPL-2.0>`_
  under the terms of the `MIT License <(http://opensource.org/licenses/MIT>`_.
+ the `XCFun library <https://xcfun.address.here/missing>`_ by Ulf
  Ekström under the terms of the `GNU General
  Public License, v2.0 <http://opensource.org/licenses/GPL-2.0>`_

 Documentation build:
  
+ the build of this documentation has been copied with minor
  adaptations from the `PCMSolver API
  <https://pcmsolver.link.here>`_. In particular the Sphinx configuration
  script (conf.py) and the Doxygen configuration file (Doxygen.in) and
  the corresponding cmake structure (FindSphinx, find_python_module).

  
Environment
===========

Classes for the solvent environment overlay

Cavity
------------

.. doxygenclass:: mrchem::Cavity
   :project: MRChem
   :members:
   :protected-members:
   :private-members:

Permittivity
------------

.. doxygenclass:: mrchem::Permittivity
   :project: MRChem
   :members:  
   :protected-members:
   :private-members: 

SCRF
------------

.. doxygenclass:: mrchem::SCRF
   :project: MRChem
   :members:  
   :protected-members:
   :private-members: 
Quantum Mechanical Functions
============================

Classes to handle quantum mechanical functions such as electronic
density, molecular orbitals.

Properties
==========

Classes for the calculation of molecular properties

Initial Guess
=============

Classes providing the initial guess of the orbitals

SCF Solver
==========

Classes for the resolution of the SCF equations of HF and DFT

QMOperators
===========

The classes that implement quantum mechanical operators

QMPotential
-----------
.. doxygenclass:: QMPotential
   :project: MRChem
   :members:
   :protected-members:
   :private-members:

XCOperator
----------
.. doxygenclass:: XCOperator
   :project: MRChem
   :members:
   :protected-members:
   :private-members:

XCPotential
-----------
.. doxygenclass:: XCPotential
   :project: MRChem
   :members:
   :protected-members:
   :private-members:

ReactionPotential
-----------------
.. doxygenclass:: mrchem::ReactionPotential
   :project: MRChem
   :members:
   :protected-members:
   :private-members:
Chemistry
=========

Classes for the chemistry overlay

:orphan:

----------------
User output file
----------------

.. raw:: html

    <style> .red {color:#aa0060; font-weight:bold; font-size:18px} </style>

.. role:: red

.. This documentation was autogenerated using parselglossy. Editing by hand is not recommended.

====================
User input reference
====================

- Keywords without a default value are **required**.
- Default values are either explicit or computed from the value of other keywords in the input.
- Sections where all keywords have a default value can be omitted.
- Predicates, if present, are the functions run to validate user input.

:red:`Keywords`
 :world_prec: Overall relative precision in the calculation. 

  **Type** ``float``

  **Predicates**
    - ``1.0e-10 < value < 1.0``

 :world_size: Total size of computational domain given as 2**(``world_size``). Always cubic and symmetric around the origin. Negative value means it will be computed from the molecular geometry. 

  **Type** ``int``

  **Default** ``-1``

  **Predicates**
    - ``value <= 10``

 :world_unit: Length unit for *all* coordinates given in user input. Everything will be converted to atomic units (bohr) before the main executable is launched, so the JSON input is *always* given in bohrs. 

  **Type** ``str``

  **Default** ``bohr``

  **Predicates**
    - ``value.lower() in ["bohr", "angstrom"]``

 :world_origin: Global gauge origin of the calculation. 

  **Type** ``List[float]``

  **Default** ``[0.0, 0.0, 0.0]``

  **Predicates**
    - ``len(value) == 3``

:red:`Sections`
 :Precisions: Define specific precision parameters. 

  :red:`Keywords`
   :exchange_prec: Precision parameter used in construction of Exchange operators. Negative value means it will follow the dynamic precision in SCF. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :helmholtz_prec: Precision parameter used in construction of Helmholtz operators. Negative value means it will follow the dynamic precision in SCF. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :poisson_prec: Precision parameter used in construction of Poisson operators. 
  
    **Type** ``float``
  
    **Default** ``user['world_prec']``
  
    **Predicates**
      - ``1.0e-10 < value < 1.0``
  
   :nuclear_prec: Precision parameter used in smoothing and projection of nuclear potential. 
  
    **Type** ``float``
  
    **Default** ``user['world_prec']``
  
    **Predicates**
      - ``1.0e-10 < value < 1.0``
  
 :Printer: Define variables for printed output. 

  :red:`Keywords`
   :print_level: Level of detail in the written output. Level 0 for production calculations, negative level for complete silence. 
  
    **Type** ``int``
  
    **Default** ``0``
  
   :print_mpi: Write separate output from each MPI to file called ``<file_name>-<mpi-rank>.out``. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :print_prec: Number of digits in property output (energies will get twice this number of digits). 
  
    **Type** ``int``
  
    **Default** ``6``
  
    **Predicates**
      - ``0 < value < 10``
  
   :print_width: Line width of printed output (in number of characters). 
  
    **Type** ``int``
  
    **Default** ``75``
  
    **Predicates**
      - ``50 < value < 100``
  
 :Plotter: Give details regarding the density and orbital plots. Three types of plots are available, line, surface and cube, and the plotting ranges are defined by three vectors (A, B and C) and an origin (O): ``line``: plots on line spanned by A, starting from O. ``surf``: plots on surface spanned by A and B, starting from O. ``cube``: plots on volume spanned by A, B and C, starting from O. 

  :red:`Keywords`
   :path: File path to plot directory. 
  
    **Type** ``str``
  
    **Default** ``plots``
  
    **Predicates**
      - ``value[-1] != '/'``
  
   :type: Type of plot: line (1D), surface (2D) or cube (3D). 
  
    **Type** ``str``
  
    **Default** ``cube``
  
    **Predicates**
      - ``value.lower() in ['line', 'surf', 'cube']``
  
   :points: Number of points in each direction on the cube grid. 
  
    **Type** ``List[int]``
  
    **Default** ``[20, 20, 20]``
  
    **Predicates**
      - ``all(p > 0 for p in value)``
      - ``not (user['Plotter']['type'] == 'line' and len(value) < 1)``
      - ``not (user['Plotter']['type'] == 'surf' and len(value) < 2)``
      - ``not (user['Plotter']['type'] == 'cube' and len(value) < 3)``
  
   :O: Origin of plotting ranges. 
  
    **Type** ``List[float]``
  
    **Default** ``[0.0, 0.0, 0.0]``
  
    **Predicates**
      - ``len(value) == 3``
  
   :A: First boundary vector for plot. 
  
    **Type** ``List[float]``
  
    **Default** ``[1.0, 0.0, 0.0]``
  
    **Predicates**
      - ``len(value) == 3``
  
   :B: Second boundary vector for plot. 
  
    **Type** ``List[float]``
  
    **Default** ``[0.0, 1.0, 0.0]``
  
    **Predicates**
      - ``len(value) == 3``
  
   :C: Third boundary vector for plot. 
  
    **Type** ``List[float]``
  
    **Default** ``[0.0, 0.0, 1.0]``
  
    **Predicates**
      - ``len(value) == 3``
  
 :MPI: Define MPI related parameters. 

  :red:`Keywords`
   :numerically_exact: This will use MPI algorithms that guarantees that the output is invariant wrt the number of MPI processes. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :shared_memory_size: Size (MB) of the MPI shared memory blocks of each shared function. 
  
    **Type** ``int``
  
    **Default** ``10000``
  
   :share_nuclear_potential: This will use MPI shared memory for the nuclear potential. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :share_coulomb_potential: This will use MPI shared memory for the Coulomb potential. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :share_xc_potential: This will use MPI shared memory for the exchange-correlation potential. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :bank_size: Number of MPI processes exclusively dedicated to manage orbital bank. 
  
    **Type** ``int``
  
    **Default** ``-1``
  
 :Basis: Define polynomial basis. 

  :red:`Keywords`
   :order: Polynomial order of multiwavelet basis. Negative value means it will be set automatically based on the world precision. 
  
    **Type** ``int``
  
    **Default** ``-1``
  
   :type: Polynomial type of multiwavelet basis. 
  
    **Type** ``str``
  
    **Default** ``interpolating``
  
    **Predicates**
      - ``value.lower() in ['interpolating', 'legendre']``
  
 :Derivatives: Define various derivative operators used in the code. 

  :red:`Keywords`
   :kinetic: Derivative used in kinetic operator. 
  
    **Type** ``str``
  
    **Default** ``abgv_55``
  
   :h_b_dip: Derivative used in magnetic dipole operator. 
  
    **Type** ``str``
  
    **Default** ``abgv_00``
  
   :h_m_pso: Derivative used in paramagnetic spin-orbit operator. 
  
    **Type** ``str``
  
    **Default** ``abgv_00``
  
 :Molecule: Define molecule. 

  :red:`Keywords`
   :charge: Total charge of molecule. 
  
    **Type** ``int``
  
    **Default** ``0``
  
   :multiplicity: Spin multiplicity of molecule. 
  
    **Type** ``int``
  
    **Default** ``1``
  
   :translate: Translate coordinates such that center of mass coincides with the global gauge origin. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :coords: Coordinates in xyz format. Atoms can be given either using atom symbol or atom number 
  
    **Type** ``str``
  
 :WaveFunction: Define the wavefunction method. 

  :red:`Keywords`
   :method: Wavefunction method. See predicates for valid methods. ``hf``, ``hartreefock`` and ``hartree-fock`` all mean the same thing, while ``lda`` is an alias for ``svwn5``. You can set a non-standard DFT functional (e.g. varying the amount of exact exchange) by choosing ``dft`` and specifing the functional(s) in the ``DFT`` section below. 
  
    **Type** ``str``
  
    **Predicates**
      - ``value.lower() in ['core', 'hartree', 'hf', 'hartreefock', 'hartree-fock', 'dft', 'lda', 'svwn3', 'svwn5', 'pbe', 'pbe0', 'bpw91', 'bp86', 'b3p86', 'b3p86-g', 'blyp', 'b3lyp', 'b3lyp-g', 'olyp', 'kt1', 'kt2', 'kt3']``
  
   :restricted: Use spin restricted wavefunction. 
  
    **Type** ``bool``
  
    **Default** ``True``
  
 :DFT: Define the exchange-correlation functional in case of DFT. 

  :red:`Keywords`
   :density_cutoff: Hard cutoff for passing density values to XCFun. 
  
    **Type** ``float``
  
    **Default** ``0.0``
  
   :functionals: List of density functionals with numerical coefficient. E.g. for PBE0 ``EXX 0.25``, ``PBEX 0.75``, ``PBEC 1.0``, see XCFun documentation <https://xcfun.readthedocs.io/>_. 
  
    **Type** ``str``
  
    **Default** `` ``
  
   :spin: Use spin separated density functionals. 
  
    **Type** ``bool``
  
    **Default** ``not(user['WaveFunction']['restricted'])``
  
 :Properties: Provide a list of properties to compute (total SCF energy and orbital energies are always computed). 

  :red:`Keywords`
   :dipole_moment: Compute dipole moment. 
  
    **Type** ``bool``
  
    **Default** ``True``
  
   :quadrupole_moment: Compute quadrupole moment. Note: Gauge origin dependent, should be used with ``translate = true`` in Molecule. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :polarizability: Compute polarizability tensor. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :magnetizability: Compute magnetizability tensor. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :nmr_shielding: Compute NMR shielding tensor. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :geometric_derivative: Compute geometric derivative. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :plot_density: Plot converged electron density. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :plot_orbitals: Plot converged molecular orbitals from list of indices, negative index plots all orbitals. 
  
    **Type** ``List[int]``
  
    **Default** ``[]``
  
 :ExternalFields: Define external electromagnetic fields. 

  :red:`Keywords`
   :electric_field: Strength of external electric field. 
  
    **Type** ``List[float]``
  
    **Default** ``[]``
  
    **Predicates**
      - ``len(value) == 0 or len(value) == 3``
  
 :Polarizability: Give details regarding the polarizability calculation. 

  :red:`Keywords`
   :frequency: List of external field frequencies. 
  
    **Type** ``List[float]``
  
    **Default** ``[0.0]``
  
 :NMRShielding: Give details regarding the NMR shileding calculation. 

  :red:`Keywords`
   :nuclear_specific: Use nuclear specific perturbation operator (h_m_pso). 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :nucleus_k: List of nuclei to compute. Negative value computes all nuclei. 
  
    **Type** ``List[int]``
  
    **Default** ``[-1]``
  
 :Files: Defines file paths used for program input/output. Note: all paths must be given in quotes if they contain slashes "path/to/file". 

  :red:`Keywords`
   :guess_basis: File name for GTO basis set, used with ``gto`` guess. 
  
    **Type** ``str``
  
    **Default** ``initial_guess/mrchem.bas``
  
   :guess_gto_p: File name for paired orbitals, used with ``gto`` guess. 
  
    **Type** ``str``
  
    **Default** ``initial_guess/mrchem.mop``
  
   :guess_gto_a: File name for alpha orbitals, used with ``gto`` guess. 
  
    **Type** ``str``
  
    **Default** ``initial_guess/mrchem.moa``
  
   :guess_gto_b: File name for beta orbitals, used with ``gto`` guess. 
  
    **Type** ``str``
  
    **Default** ``initial_guess/mrchem.mob``
  
   :guess_phi_p: File name for paired orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/phi_p_scf_idx_<0...Np>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_p``
  
   :guess_phi_a: File name for alpha orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/phi_a_scf_idx_<0...Na>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_a``
  
   :guess_phi_b: File name for beta orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/phi_b_scf_idx_<0...Nb>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_b``
  
   :guess_x_p: File name for paired response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/x_p_rsp_idx_<0...Np>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/X_p``
  
   :guess_x_a: File name for alpha response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/x_a_rsp_idx_<0...Na>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/X_a``
  
   :guess_x_b: File name for beta response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/x_b_rsp_idx_<0...Nb>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/X_b``
  
   :guess_y_p: File name for paired response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/y_p_rsp_idx_<0...Np>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/Y_p``
  
   :guess_y_a: File name for alpha response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/y_a_rsp_idx_<0...Na>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/Y_a``
  
   :guess_y_b: File name for beta response orbitals, used with ``mw`` guess. Expected path is ``<path_orbitals>/y_b_rsp_idx_<0...Nb>_<re/im>.mw 
  
    **Type** ``str``
  
    **Default** ``initial_guess/Y_b``
  
   :guess_cube_p: File name for paired orbitals, used with ``cube`` guess. Expected path is ``<path_orbitals>/phi_p_scf_idx_<0...Np>_<re/im>.cube 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_p``
  
   :guess_cube_a: File name for alpha orbitals, used with ``cube`` guess. Expected path is ``<path_orbitals>/phi_a>_scf_idx_<0...Na>_<re/im>.cube 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_a``
  
   :guess_cube_b: File name for beta orbitals, used with ``cube`` guess. Expected path is ``<path_orbitals>/phi_b_scf_idx_<0...Nb>_<re/im>.cube 
  
    **Type** ``str``
  
    **Default** ``initial_guess/phi_b``
  
   :cube_vectors: Directory where cube vectors are stored for mrchem calculation. 
  
    **Type** ``str``
  
    **Default** ``cube_vectors/``
  
 :SCF: Includes parameters related to the ground state SCF orbital optimization. 

  :red:`Keywords`
   :run: Run SCF solver. Otherwise properties are computed on the initial orbitals. 
  
    **Type** ``bool``
  
    **Default** ``True``
  
   :max_iter: Maximum number of SCF iterations. 
  
    **Type** ``int``
  
    **Default** ``100``
  
   :kain: Length of KAIN iterative history. 
  
    **Type** ``int``
  
    **Default** ``5``
  
   :rotation: Number of iterations between each diagonalization/localization. 
  
    **Type** ``int``
  
    **Default** ``0``
  
   :localize: Use canonical or localized orbitals. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :energy_thrs: Convergence threshold for SCF energy. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :guess_prec: Precision parameter used in construction of initial guess. 
  
    **Type** ``float``
  
    **Default** ``0.001``
  
    **Predicates**
      - ``1.0e-10 < value < 1.0``
  
   :guess_screen: Screening parameter used in GTO evaluations, in number of standard deviations. Every coordinate beyond N StdDev from the Gaussian center is evaluated to zero. Note that too aggressive screening is counter productive, because it leads to a sharp cutoff in the resulting function which requires higher grid refinement. Negative value means no screening. 
  
    **Type** ``float``
  
    **Default** ``12.0``
  
   :start_prec: Incremental precision in SCF iterations, initial value. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :final_prec: Incremental precision in SCF iterations, final value. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :guess_type: Type of initial guess for ground state orbitals. ``chk`` restarts a previous calculation which was dumped using the ``write_checkpoint`` keyword. This will load MRA and electron spin configuration directly from the checkpoint files, which are thus required to be identical in the two calculations. ``mw`` will start from final orbitals in a previous calculation written using the ``write_orbitals`` keyword. The orbitals will be re-projected into the new computational setup, which means that the electron spin configuration and MRA can be different in the two calculations. ``gto`` reads precomputed GTO orbitals (requires extra non-standard input files for basis set and MO coefficients). ``core`` and ``sad`` will diagonalize the Fock matrix in the given AO basis (SZ, DZ, TZ or QZ) using a Core or Superposition of Atomic Densities Hamiltonian, respectively. 
  
    **Type** ``str``
  
    **Default** ``sad_dz``
  
    **Predicates**
      - ``value.lower() in ['mw', 'chk', 'gto', 'core_sz', 'core_dz', 'core_tz', 'core_qz', 'sad_sz', 'sad_dz', 'sad_tz', 'sad_qz', 'sad_gto', 'cube']``
  
   :write_checkpoint: Write orbitals to disk in each iteration, file name ``<path_checkpoint>/phi_scf_idx_<0..N>``. Can be used as ``chk`` initial guess in subsequent calculations. Note: must be given in quotes if there are slashes in the path "path/to/checkpoint". 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :path_checkpoint: Path to checkpoint files during SCF, used with ``write_checkpoint`` and ``chk`` guess. 
  
    **Type** ``str``
  
    **Default** ``checkpoint``
  
    **Predicates**
      - ``value[-1] != '/'``
  
   :write_orbitals: Write final orbitals to disk, file name ``<path_orbitals>/phi_<p/a/b>_scf_idx_<0..Np/Na/Nb>``. Can be used as ``mw`` initial guess in subsequent calculations. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :path_orbitals: Path to where converged orbitals will be written in connection with the ``write_orbitals`` keyword. Note: must be given in quotes if there are slashes in the path "path/to/orbitals". 
  
    **Type** ``str``
  
    **Default** ``orbitals``
  
    **Predicates**
      - ``value[-1] != '/'``
  
   :orbital_thrs: Convergence threshold for orbital residuals. 
  
    **Type** ``float``
  
    **Default** ``10 * user['world_prec']``
  
 :Response: Includes parameters related to the response SCF optimization. 

  :red:`Keywords`
   :run: In which Cartesian directions to run response solver. 
  
    **Type** ``List[bool]``
  
    **Default** ``[True, True, True]``
  
   :max_iter: Maximum number of response iterations. 
  
    **Type** ``int``
  
    **Default** ``100``
  
   :kain: Length of KAIN iterative history. 
  
    **Type** ``int``
  
    **Default** ``5``
  
   :property_thrs: Convergence threshold for symmetric property. Symmetric meaning the property computed from the same operator as the response purturbation, e.g. for external magnetic field the symmetric property corresponds to the magnetizability (NMR shielding in non-symmetric, since one of the operators is external magnetic field, while the other is nuclear magnetic moment). 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :start_prec: Incremental precision in SCF iterations, initial value. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :final_prec: Incremental precision in SCF iterations, final value. 
  
    **Type** ``float``
  
    **Default** ``-1.0``
  
   :guess_prec: Precision parameter used in construction of initial guess. 
  
    **Type** ``float``
  
    **Default** ``0.001``
  
    **Predicates**
      - ``1.0e-10 < value < 1.0``
  
   :guess_type: Type of initial guess for response. ``none`` will start from a zero guess for the response functions. ``chk`` restarts a previous calculation which was dumped using the ``write_checkpoint`` keyword. ``mw`` will start from final orbitals in a previous calculation written using the ``write_orbitals`` keyword. The orbitals will be re-projected into the new computational setup. 
  
    **Type** ``str``
  
    **Default** ``none``
  
    **Predicates**
      - ``value.lower() in ['none', 'chk', 'mw']``
  
   :write_checkpoint: Write perturbed orbitals to disk in each iteration, file name ``<path_checkpoint>/<X/Y>_rsp_<direction>_idx_<0..N>``. Can be used as ``chk`` initial guess in subsequent calculations. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :path_checkpoint: Path to checkpoint files during SCF, used with ``write_checkpoint`` and ``chk`` guess. 
  
    **Type** ``str``
  
    **Default** ``checkpoint``
  
    **Predicates**
      - ``value[-1] != '/'``
  
   :write_orbitals: Write final perturbed orbitals to disk, file name ``<path_orbitals>/<X/Y>_<p/a/b>_rsp_<direction>_idx_<0..Np/Na/Nb>``. Can be used as ``mw`` initial guess in subsequent calculations. 
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :path_orbitals: Path to where converged orbitals will be written in connection with the ``write_orbitals`` keyword. 
  
    **Type** ``str``
  
    **Default** ``orbitals``
  
    **Predicates**
      - ``value[-1] != '/'``
  
   :orbital_thrs: Convergence threshold for orbital residuals. 
  
    **Type** ``float``
  
    **Default** ``10 * user['world_prec']``
  
   :localize: Use canonical or localized unperturbed orbitals. 
  
    **Type** ``bool``
  
    **Default** ``user['SCF']['localize']``
  
 :Environment: Includes parameters related to the computation of the reaction field energy of a system in an environment. 

  :red:`Keywords`
   :max_iter: Max number of iterations allowed in the nested procedure. 
  
    **Type** ``int``
  
    **Default** ``100``
  
   :run_environment: Perform the reaction field calculation of the reaction potential of the interaction between environment and molecule.  
  
    **Type** ``bool``
  
    **Default** ``False``
  
   :algorithm: What algorithm to use for the reaction field ``scrf`` runs a nested algorithm where the generalized Poisson equation is solved iterativelly until self consistency wrt. the convergence threshold. 
  
    **Type** ``str``
  
    **Default** ``scrf``
  
    **Predicates**
      - ``value.lower() in ['scrf']``
  
   :convergence_criterion: Adjust the convergence threshold for the nested procedure. ``dynamic`` Uses the absolute value of the latest orbital update as convergence threshold. When the orbitals are close to convergence (``mo_residual < world_prec*10``) the convergence threshold will be equal to ``world_prec``. ``static`` uses ``world_prec`` as convergence threshold. 
  
    **Type** ``str``
  
    **Default** ``dynamic``
  
    **Predicates**
      - ``value.lower() in ['dynamic', 'static']``
  
   :extrapolate_Vr: Extrapolate on the reaction potential if true, or on the surface charge distribution in the convergence acceleration. 
  
    **Type** ``bool``
  
    **Default** ``True``
  
   :density_type: What part of the total molecular charge density to use in the algorithm. ``total`` uses the total charge density. ``nuclear`` uses only the nuclear part of the total charge density. ``electronic`` uses only the electronic part of the total charge density. 
  
    **Type** ``str``
  
    **Default** ``total``
  
    **Predicates**
      - ``value.lower() in ['total', 'nuclear', 'electronic']``
  
   :kain: Number of previous reaction field iterates kept for convergence acceleration during the nested precedure. 
  
    **Type** ``int``
  
    **Default** ``user['SCF']['kain']``
  
  :red:`Sections`
   :Cavity: Define the interlocking spheres cavity. 
  
      :red:`Keywords`
       :spheres: Coordinates and radii  of the spheres written as $spheres x_0    y_0    z_0    R_0 ... x_N    y_N    z_N    R_N $end The units used are the same specified with the `world_unit` keyword. 
      
        **Type** ``str``
      
        **Default** ````
      
       :cavity_width: Width of cavity boundary 
      
        **Type** ``float``
      
        **Default** ``0.2``
      
   :Permittivity: Parameters for the permittivity function. 
  
      :red:`Keywords`
       :epsilon_in: Permittivity inside the cavity. 1.0 is the permittivity of free space, anything other than this is undefined behaviour. 
      
        **Type** ``float``
      
        **Default** ``1.0``
      
       :epsilon_out: Permittivity outside the cavity. This is characteristic of the solvent used. 
      
        **Type** ``float``
      
        **Default** ``2.0``
      
       :formulation: Formulation of the Permittivity function. Currently only the exponential is used. 
      
        **Type** ``str``
      
        **Default** ``exponential``
      
        **Predicates**
          - ``value.lower() in ['exponential']``
      ---------------
User input file
---------------

The input file is organized in sections and keywords that can be of different
type. Input keywords and sections are **case-sensitive**, while `values` are
**case-insensitive**.

.. code-block:: bash

    Section {
      keyword_1 = 1                         # int
      keyword_2 = 3.14                      # float
      keyword_3 = [1, 2, 3]                 # int array
      keyword_4 = foo                       # string
      keyword_5 = true                      # boolean
    }

Valid options for booleans are ``true/false``, ``on/off`` or ``yes/no``. Single
word strings can be given without quotes (be careful of special characters, like
slashes in file paths). A  complete list of available input keywords can be
found in the :ref:`User input reference`.

Top section
-----------

The main input section contain four keywords: the relative precision
:math:`\epsilon_{rel}` that will be guaranteed in the calculation and the size,
origin and unit of the computational domain. The top section is not specified
by name, just write the keywords directly, e.g

.. code-block:: bash

    world_prec = 1.0e-5                     # Overall relative precision
    world_size = 5                          # Size of domain 2^{world_size}
    world_unit = bohr                       # Global length unit
    world_origin = [0.0, 0.0, 0.0]          # Global gauge origin

The relative precision sets an upper limit for the number of correct digits
you are expected to get out of the computation (note that
:math:`\epsilon_{rel}=10^{-6}` yields :math:`\mu` Ha accuracy for the hydrogen
molecule, but only mHa accuracy for benzene).

The computational domain is always symmetric around the origin, with *total*
size given by the ``world_size`` parameter as :math:`[2^n]^3`, e.i.
``world_size = 5`` gives a domain of :math:`[-16,16]^3`.
Make sure that the world is large enough to allow the molecular density to
reach zero on the boundary. The ``world_size`` parameter can be left out,
in which case the size will be estimated based on the molecular geometry.
The ``world_unit`` relates to **all** coordinates given in the input file and
can be one of two options: ``angstrom`` or ``bohr``.

.. note::

    The ``world_size`` will be only approximately scaled by the angstrom unit,
    by adding an extra factor of 2 rather than the appropriate factor of ~1.89.
    This means that e.g. ``world_size = 5`` (:math:`[-16,16]^3`) with
    ``world_unit = angstrom`` will be translated into :math:`[-32,32]^3` bohrs.

Precisions
----------

MRChem uses a smoothed nuclear potential to avoid numerical problems in
connection with the :math:`Z/|r-R|` singularity. The smoothing is controlled by
a single parameter ``nuc_prec`` that is related to the expected error in the
energy due to the smoothing. There are also different precision parameters for
the `construction` of the Poisson and Helmholtz integral operators.

.. code-block:: bash

    Precisions {
      nuclear_prec = 1.0e-6                 # For construction of nuclear potential
      poisson_prec = 1.0e-6                 # For construction of Poisson operators
      helmholtz_prec = 1.0e-6               # For construction of Helmholtz operatos
    }

By default, all precision parameters follow ``world_prec`` and usually don't
need to be changed.

Printer
-------

This section controls the format of the printed output file (``.out``
extension). The most important option is the ``print_level``, but it also gives
options for number of digits in the printed output, as well as the line width
(defaults are shown):

.. code-block:: bash

    Printer {
      print_level = 0                       # Level of detail in the printed output
      print_width = 75                      # Line width (in characters) of printed output
      print_prec = 6                        # Number of digits in floating point output
    }

Note that energies will be printed with *twice* as many digits.
Available print levels are:

- ``print_level=-1`` no output is printed
- ``print_level=0`` prints mainly properties
- ``print_level=1`` adds timings for individual steps
- ``print_level=2`` adds memory and timing information on ``OrbitalVector`` level
- ``print_level=3`` adds memory and timing information on ``Orbital`` level
- ``print_level>10`` adds a *lot* more output from deep within MRCPP


MPI
---

This section defines some parameters that are used in MPI runs (defaults shown):

.. code-block:: bash

    MPI {
      bank_size = -1                        # Number of processes used as memory bank
      numerically_exact = false             # Guarantee MPI invariant results
      share_nuclear_potential = false       # Use MPI shared memory window
      share_coulomb_potential = false       # Use MPI shared memory window
      share_xc_potential = false            # Use MPI shared memory window
    }

The memory bank will allow larger molecules to get though if memory is the
limiting factor, but it will be slower, as the bank processes will not take
part in any computation. For calculations involving exact exchange (Hartree-Fock
or hybrid DFT functionals) a memory bank is **required** whenever there's more
than one MPI process. A negative bank size will set it automatically based on
the number of available processes. For pure DFT functionals on smaller molecules
it is likely more efficient to set `bank_size = 0`, otherwise it's recommended
to use the default. If a particular calculation runs out of memory, it might
help to increase the number of bank processes from the default value.

The ``numerically_exact`` keyword will trigger algorithms that guarantee that
the computed results are invariant (within double precision) with respect to
the number or MPI processes. These exact algorithms require more memory and are
thus not default. Even when the numbers are *not* MPI invariant they should be
correct and identical within the chosen ``world_prec``.

The ``share_potential`` keywords are used to share the memory space for the
particular functions between all processes located on the same physical machine.
This will save memory but it might slow the calculation down, since the shared
memory cannot be "fast" memory (NUMA) for all processes at once.


Basis
-----

This section defines the polynomial MultiWavelet basis

.. code-block:: bash

    Basis {
      type = Interpolating                  # Legendre or Interpolating
      order = 7                             # Polynomial order of MW basis
    }

The MW basis is defined by the polynomial order :math:`k`, and the type of
scaling functions: Legendre or Interpolating polynomials (in the current
implementation it doesn't really matter which type you choose). Note that
increased precision requires higher polynomial order (use e.g :math:`k = 5`
for :math:`\epsilon_{rel} = 10^{-3}`, and :math:`k = 13` for
:math:`\epsilon_{rel} = 10^{-9}`, and interpolate in between). If the ``order``
keyword is left out it will be set automatically according to

.. math:: k=-1.5*log_{10}(\epsilon_{rel})

The Basis section can usually safely be omitted in the input.

Molecule
--------

This input section specifies the geometry (given in ``world_unit`` units),
charge and spin multiplicity of the molecule, e.g. for water (coords must be
specified, otherwise defaults are shown):

.. code-block:: bash

    Molecule {
      charge = 0                            # Total charge of molecule
      multiplicity = 1                      # Spin multiplicity
      translate = false                     # Translate CoM to world_origin
    $coords
    O   0.0000     0.0000     0.0000        # Atomic symbol and coordinate
    H   0.0000     1.4375     1.1500        # Atomic symbol and coordinate
    H   0.0000    -1.4375     1.1500        # Atomic symbol and coordinate
    $end
    }

Since the computational domain is always cubic and symmetric around the origin
it is usually a good idea to ``translate`` the molecule to the origin (as long
as the ``world_origin`` is the true origin).

WaveFunction
------------

Here we give the wavefunction method and whether we run spin restricted (alpha
and beta spins are forced to occupy the same spatial orbitals) or not (method
must be specified, otherwise defaults are shown):

.. code-block:: bash

    WaveFunction {
      method = <wavefunction_method>        # Core, Hartree, HF or DFT
      restricted = true                     # Spin restricted/unrestricted
    }

There are currently four methods available: Core Hamiltonian, Hartree,
Hartree-Fock (HF) and Density Functional Theory (DFT). When running DFT you can
*either* set one of the default functionals in this section (e.g. ``method =
B3LYP``), *or* you can set ``method = DFT`` and specify a "non-standard"
functional in the separate DFT section (see below). See
:ref:`User input reference` for a list of available default functionals.

.. note::

    Restricted open-shell wavefunctions are not supported.

DFT
---

This section can be omitted if you are using a default functional, see above.
Here we specify the exchange-correlation functional used in DFT
(functional names must be specified, otherwise defaults are shown)

.. code-block:: bash

    DFT {
      spin = false                          # Use spin-polarized functionals
      density_cutoff = 0.0                  # Cutoff to set XC potential to zero
    $functionals
    <func1>     1.0                         # Functional name and coefficient
    <func2>     1.0                         # Functional name and coefficient
    $end
    }

You can specify as many functionals as you want, and they will be added on top
of each other with the given coefficient. Both exchange and correlation
functionals must be set explicitly, e.g. ``SLATERX`` and ``VWN5C`` for the
standard LDA functional. For hybrid functionals you must specify the amount
of exact Hartree-Fock exchange as a separate functional
``EXX`` (``EXX 0.2`` for B3LYP and ``EXX 0.25`` for PBE0 etc.). Option to use
spin-polarized functionals or not. Unrestricted calculations will use
spin-polarized functionals by default. The XC functionals are provided by the
`XCFun <https://github.com/dftlibs/xcfun>`_ library.

Properties
----------

Specify which properties to compute. By default, only the ground state SCF
energy as well as orbital energies will be computed. Currently the following
properties are available (all but the dipole moment are ``false`` by default)

.. code-block:: bash

    Properties {
      dipole_moment = true                  # Compute dipole moment
      quadrupole_moment = false             # Compute quadrupole moment
      polarizabiltity = false               # Compute polarizability
      magnetizability = false               # Compute magnetizability
      nmr_shielding = false                 # Compute NMR shieldings
      geometric_derivative = false          # Compute geometric derivative
      plot_density = false                  # Plot converged density
      plot_orbitals = []                    # Plot converged orbitals
    }

Some properties can be further specified in dedicated sections.

.. warning:: The computation of the molecular gradient suffers greatly from
   numerical noise.  The code replaces the nucleus-electron attraction with a
   smoothed potential. This can only partially recover the nuclear cusps, even
   with tight precision.  The molecular gradient is only suited for use in
   geometry optimization of small molecules and with tight precision thresholds.

Polarizability
++++++++++++++
The polarizability can be computed with several frequencies (by default only
static polarizability is computed):


.. code-block:: bash

    Polarizability {
      frequency = [0.0, 0.0656]             # List of frequencies to compute
    }

NMRShielding
++++++++++++

For the NMR shielding we can specify a list of nuclei to compute (by default
all nuclei are computed):

.. code-block:: bash

    NMRShielding {
      nuclear_specific = false              # Use nuclear specific perturbation operator
      nucleus_k = [0,1,2]                   # List of nuclei to compute (-1 computes all)
    }

The ``nuclear_specific`` keyword triggers response calculations using the
nuclear magnetic moment operator instead of the external magnetic field. For
small molecules this is not recommended since it requires a separate response
calculation for each nucleus, but it might be beneficial for larger systems if
you are interested only in a single shielding constant. Note that the components
of the *perturbing* operator defines the *row* index in the output tensor, so
``nuclear_specific = true`` will result in a shielding tensor which is
the transpose of the one obtained with ``nuclear_specific = false``.

Plotter
+++++++

The ``plot_density`` and ``plot_orbitals`` properties will use the Plotter
section to specify the parameters of the plots (by default you will get a
``cube`` plot on the unit cube):

.. code-block:: bash

    Plotter {
      path = plots                          # File path to store plots
      type = cube                           # Plot type (line, surf, cube)
      points = [20, 20, 20]                 # Number of grid points
      O = [-4.0,-4.0,-4.0]                  # Plot origin
      A = [8.0, 0.0, 0.0]                   # Boundary vector
      B = [0.0, 8.0, 0.0]                   # Boundary vector
      C = [0.0, 0.0, 8.0]                   # Boundary vector
    }


The plotting grid is computed from the vectors ``O``, ``A``, ``B`` and ``C`` in
the following way:

    1.  ``line`` plot: along the vector ``A`` starting from ``O``, using
        ``points[0]`` number of points.
    2.  ``surf`` plot: on the area spanned by the vectors ``A`` and ``B`` starting
        from ``O``, using ``points[0]`` and ``points[1]`` points in each direction.
    3.  ``cube`` plot: on the volume spanned by the vectors ``A``, ``B`` and ``C``
        starting from ``O``, using ``points[0]``, ``points[1]`` and ``points[2]``
        points in each direction.

The above example will plot on a 20x20x20 grid in the volume [-4,4]^3, and the
generated files (e.g. ``plots/phi_1_re.cube``) can be viewed directly in a
web browser by `blob <https://github.com/densities/blob/>`_ , like this benzene
orbital:

.. image:: ../gfx/blob.png


SCF
---

This section specifies the parameters for the SCF optimization of the ground
state wavefunction.

SCF solver
++++++++++

The optimization is controlled by the following keywords (defaults shown):

.. code-block:: bash

    SCF {
      run = true                            # Run SCF solver
      kain = 5                              # Length of KAIN iterative subspace
      max_iter = 100                        # Maximum number of SCF iterations
      rotation = 0                          # Iterations between diagonalize/localize
      localize = false                      # Use canonical or localized  orbitals
      start_prec = -1.0                     # Dynamic precision, start value
      final_prec = -1.0                     # Dynamic precision, final value
      orbital_thrs = 10 * world_prec        # Convergence threshold orbitals
      energy_thrs = -1.0                    # Convergence threshold energy
    }

If ``run = false`` no SCF is performed, and the properties are computed directly
on the initial guess wavefunction.

The ``kain`` (Krylov Accelerated Inexact Newton) keyword gives the length of
the iterative subspace accelerator (similar to DIIS). The ``rotation`` keyword
gives the number of iterations between every orbital rotation, which can be
either localization or diagonalization, depending on the ``localize`` keyword.
The first two iterations in the SCF are always rotated, otherwise it is
controlled by the ``rotation`` keyword (usually this is not very important, but
sometimes it fails to converge if the orbitals drift too far away from the
localized/canonical forms).

The dynamic precision keywords control how the numerical precision is changed
throughout the optimization. One can choose to use a lower ``start_prec`` in
the first iterations which is gradually increased to ``final_prec`` (both are
equal to ``world_prec`` by default). Note that lower initial precision might
affect the convergence rate.

In general, the important convergence threshold is that of the orbitals,
and by default this is set one order of magnitude higher than the overall
``world_prec``. For simple energy calculations, however, it is not necessary to
converge the orbitals this much due to the quadratic convergence of the energy.
This means that the number of correct digits in the total energy will be
saturated well before this point, and one should rather use the ``energy_thrs``
keyword in this case in order to save a few iterations.

.. note::

    It is usually not feasible to converge the orbitals *beyond* the overall
    precision ``world_prec`` due to numerical noise.

Initial guess
+++++++++++++

Several types of initial guess are available:

 - ``core`` and ``sad`` requires no further input and computes guesses from
   scratch.
 - ``chk`` and ``mw`` require input files from previous MW calculations.
 - ``gto`` requires non-standard Gaussian-type orbital input files,
   and is thus not fully supported.

The ``core`` and ``sad`` guesses are computed by diagonalizing the Hamiltonian
matrix using a Core or Superposition of Atomic Densities (SAD) Hamiltonian,
respectively. The matrix is constructed in a small AO basis with a given
"zeta quality", which should be added as a suffix in the keyword. Available AO
bases are hydrogenic orbitals of single ``sz``, double ``dz``, triple ``tz``
and quadruple ``qz`` zeta size.

The ``core`` and ``sad`` guesses are fully specified with the following keywords
(defaults shown):

.. code-block:: bash

    SCF {
      guess_prec = 1.0e-3                   # Numerical precision used in guess
      guess_type = sad_dz                   # Type of inital guess (chk, mw, gto, core, sad)
    }

Checkpointing
+++++++++++++

The program can dump checkpoint files at every iteration using the
``write_checkpoint`` keyword (defaults shown):

.. code-block:: bash

    SCF {
      path_checkpoint = checkpoint          # Path to checkpoint files
      write_checkpoint = false              # Save checkpoint files every iteration
    }

This allows the calculation to be restarted in case it crashes e.g. due to time
limit or hardware failure on a cluster. This is done by setting ``guess_type =
chk`` in the subsequent calculation:

.. code-block:: bash

    SCF {
      guess_type = chk                      # Type of inital guess (chk, mw, gto, core, sad)
    }

In this case the ``path_checkpoint`` must be the same as the previous
calculation, as well as all other parameters in the calculation (Molecule and
Basis in particular).

Write orbitals
++++++++++++++

The converged orbitals can be saved to file with the ``write_orbitals`` keyword
(defaults shown):

.. code-block:: bash

    SCF {
      path_orbitals = orbitals              # Path to orbital files
      write_orbitals = false                # Save converged orbitals to file
    }

This will make individual files for each orbital under the ``path_orbitals``
directory. These orbitals can be used as starting point for subsequent
calculations using the ``guess_type = mw`` initial guess:

.. code-block:: bash

    SCF {
      guess_prec = 1.0e-3                   # Numerical precision used in guess
      guess_type = mw                       # Type of inital guess (chk, mw, gto, core, sad)
    }

Here the orbitals will be re-projected onto the current MW basis with precision
``guess_prec``. We also need to specify the paths to the input files:

.. code-block:: bash

    Files {
      guess_phi_p = initial_guess/phi_p     # Path to paired MW orbitals
      guess_phi_a = initial_guess/phi_a     # Path to alpha MW orbitals
      guess_phi_b = initial_guess/phi_b     # Path to beta MW orbitals
    }

Note that by default orbitals are written to the directory called ``orbitals``
but the ``mw`` guess reads from the directory ``initial_guess`` (this is to
avoid overwriting the files by default). So, in order to use MW orbitals from a
previous calculation, you must either change one of the paths
(``SCF.path_orbitals`` or ``Files.guess_phi_p`` etc), or manually copy the files
between the default locations.

.. note::

    The ``mw`` guess must not be confused with the ``chk`` guess, although they
    are similar. The ``chk`` guess will blindly read in the orbitals that are
    present, regardless of the current molecular structure and computational
    setup (if you run with a different computational domain or MW basis
    type/order the calculation will crash). The ``mw`` guess will re-project
    the old orbitals onto the new computational setup and populate the orbitals
    based on the *new* molecule (here the computation domain and MW basis do
    *not* have to match).


Response
--------

This section specifies the parameters for the SCF optimization of the linear
response functions. There might be several independent response calculations
depending on the requested properties, e.g.

.. code-block:: bash

    Polarizability {
      frequency = [0.0, 0.0656]             # List of frequencies to compute
    }

will run one response for each frequency (each with three Cartesian components),
while

.. code-block:: bash

    Properties {
      magnetizability = true                # Compute magnetizability
      nmr_shielding = true                  # Compute NMR shieldings
    }

will combine both properties into a single response calculation, since the
perturbation operator is the same in both cases (unless you choose
``NMRShielding.nuclear_specific = true``, in which case there will be a
different response for each nucleus).

Response solver
+++++++++++++++

The optimization is controlled by the following keywords (defaults shown):

.. code-block:: bash

    Response {
      run = [true,true,true]                # Run response solver [x,y,z] direction
      kain = 5                              # Length of KAIN iterative subspace
      max_iter = 100                        # Maximum number of SCF iterations
      localize = false                      # Use canonical or localized  orbitals
      start_prec = -1.0                     # Dynamic precision, start value
      final_prec = -1.0                     # Dynamic precision, final value
      orbital_thrs = 10 * world_prec        # Convergence threshold orbitals
    }

Each linear response calculation involves the three Cartesian components of the
appropriate perturbation operator. If any of the components of ``run`` is
``false``, no response is performed in that particular direction, and the
properties are computed directly on the initial guess response functions
(usually zero guess).

The ``kain`` (Krylov Accelerated Inexact Newton) keyword gives the length of
the iterative subspace accelerator (similar to DIIS). The ``localize`` keyword
relates to the unperturbed orbitals, and can be set independently of the
``SCF.localize`` keyword.

The dynamic precision keywords control how the numerical precision is changed
throughout the optimization. One can choose to use a lower ``start_prec`` in
the first iterations which is gradually increased to ``final_prec`` (both are
equal to ``world_prec`` by default). Note that lower initial precision might
affect the convergence rate.

For response calculations, the important convergence threshold is that of the
orbitals, and by default this is set one order of magnitude higher than the
overall ``world_prec``. 

.. note::

    The quality of the response property depends on both the perturbed as well
    as the unperturbed orbitals, so they should be equally well converged.

Initial guess
+++++++++++++

The following initial guesses are available:

 - ``none`` start from a zero guess for the response functions.
 - ``chk`` and ``mw`` require input files from previous MW calculations.

By default, no initial guess is generated for the response functions, but the
``chk`` and ``mw`` guesses work similarly as for the SCF.

Checkpointing
+++++++++++++

The program can dump checkpoint files at every iteration using the
``write_checkpoint`` keyword (defaults shown):

.. code-block:: bash

    Response {
      path_checkpoint = checkpoint          # Path to checkpoint files
      write_checkpoint = false              # Save checkpoint files every iteration
    }

This allows the calculation to be restarted in case it crashes e.g. due to time
limit or hardware failure on a cluster. This is done by setting ``guess_type =
chk`` in the subsequent calculation:

.. code-block:: bash

    Response {
      guess_type = chk                      # Type of inital guess (none, chk, mw)
    }

In this case the ``path_checkpoint`` must be the same as the previous
calculation, as well as all other parameters in the calculation (Molecule and
Basis in particular).

Write orbitals
++++++++++++++

The converged response orbitals can be saved to file with the
``write_orbitals`` keyword (defaults shown):

.. code-block:: bash

    Response {
      path_orbitals = orbitals              # Path to orbital files
      write_orbitals = false                # Save converged orbitals to file
    }

This will make individual files for each orbital under the ``path_orbitals``
directory. These orbitals can be used as starting point for subsequent
calculations using the ``guess_type = mw`` initial guess:

.. code-block:: bash

    Response {
      guess_prec = 1.0e-3                   # Numerical precision used in guess
      guess_type = mw                       # Type of inital guess (chk, mw, gto, core, sad)
    }

Here the orbitals will be re-projected onto the current MW basis with precision
``guess_prec``. We also need to specify the paths to the input files (only X
for static perturbations, X and Y for dynamic perturbations):

.. code-block:: bash

    Files {
      guess_X_p = initial_guess/X_p         # Path to paired MW orbitals
      guess_X_a = initial_guess/X_a         # Path to alpha MW orbitals
      guess_X_b = initial_guess/X_b         # Path to beta MW orbitals
      guess_Y_p = initial_guess/Y_p         # Path to paired MW orbitals
      guess_Y_a = initial_guess/Y_a         # Path to alpha MW orbitals
      guess_Y_b = initial_guess/Y_b         # Path to beta MW orbitals
    }

Note that by default orbitals are written to the directory called ``orbitals``
but the ``mw`` guess reads from the directory ``initial_guess`` (this is to
avoid overwriting the files by default). So, in order to use MW orbitals from a
previous calculation, you must either change one of the paths
(``Response.path_orbitals`` or ``Files.guess_X_p`` etc), or manually copy the
files between the default locations.

=============
User's Manual
=============

The MRChem program comes as two executables::

    <install-path>/bin/mrchem                   # Python input parser and launcher
    <install-path>/bin/mrchem.x                 # MRChem main executable

where the former is a Python script that reads and validates the *user input
file* and produces a new *program input file* which is then passed as argument
to the latter, which is the actual C++ executable.

The input and output of the program is thus organized as *three* separate files:

+-------------------+---------------------------+---------------+
| File extension    | Description               | Format        |
+===================+===========================+===============+
| ``.inp``          | User input file           | GETKW/JSON    | 
+-------------------+---------------------------+---------------+
| ``.json``         | Program input/output      | JSON          |
+-------------------+---------------------------+---------------+
| ``.out``          | User output file          | Text          |
+-------------------+---------------------------+---------------+

The name of the user input file can be anything, as long as it has the ``.inp``
extension, and the corresponding ``.json`` and ``.out`` files will get the
same name prefix. The JSON program file will get both an ``"input"`` and an
``"output"`` section. This ``"input"`` section is rather detailed and contains
very implementation specific keywords, but it is automatically generated by the
``mrchem`` script, based on the more generic keywords of the user input file.
The ``mrchem`` script will further launch the ``mrchem.x`` main executable,
which will produce the text output file as well as the ``"output"`` section
of the JSON in/out file. The contents of all these files will be discussed
in more detail in the sections below.

.. toctree::
   :maxdepth: 1

   running
   user_inp
   user_ref
   program_json
-------------------------
Program input/output file
-------------------------


Input schema
------------

.. literalinclude:: schema_input.json
  :language: JSON

Output schema
-------------

.. literalinclude:: schema_output.json
  :language: JSON
-------------------
Running the program
-------------------

In the following we will assume to have a valid user input file for the water
molecule called ``h2o.inp``, e.g. like this

.. literalinclude:: h2o_getkw.inp

To run the calculation, pass the file name (without extension) as argument
to the ``mrchem`` script (make sure you understand the difference between the
``.inp``, ``.json`` and ``.out`` file, as described in the previous section)::

    $ mrchem h2o

This will under the hood actually do the following two steps::

    $ mrchem h2o.inp > h2o.json
    $ mrchem.x h2o.json > h2o.out

The first step includes input validation, which means that everything
that passes this step is a well-formed computation.


Dry-running the input parser
----------------------------

The execution of the two steps above can be done separately by dry-running the
parser script::

    $ mrchem --dryrun h2o

This will run only the input validation part and generate the ``h2o.json``
program input, but it will *not* launch the main executable ``mrchem.x``.
This can then be done manually in a subsequent step by calling::

    $ mrchem.x h2o.json

This separation can be useful for instance for developers or advanced users
who want to change some automatically generated input values before launching
the actual program, see :ref:`Input schema`.

Printing to standard output
---------------------------

By default the program will write to the text output file (``.out`` extension),
but if you rather would like it printed in the terminal you can add the
``--stdout`` option (then no text output file is created)::

    $ mrchem --stdout h2o

Reproducing old calculations
----------------------------

The JSON in/out file acts as a full record of the calculation, and can be
used to reproduce old results. Simply pass the JSON file once more
to ``mrchem.x``, and the ``"output"`` section will be overwritten::

    $ mrchem.x h2o.json

User input in JSON format
-------------------------

The user input file can be written in JSON format instead of the standard
syntax which is described in detail below. This is very convenient if you
have for instance a Python script to generate input files. The water
example above in JSON format reads (the ``coords`` string is not very elegant,
but unfortunately that's just how JSON works...):

.. literalinclude:: h2o_json.inp

which can be passed to the input parser with the ``--json`` option::

    $ mrchem --json h2o

.. note::

    A *user input file* in JSON format must **NOT** be confused with the JSON
    in/out file for the ``mrchem.x`` program. The file should still have a
    ``.inp`` extension, and contain all the same keywords which have to be
    validated and translated by the ``mrchem`` script into the ``.json``
    *program input file*.


Parallel execution
------------------

The MRChem program comes with support for both shared memory and distributed
memory parallelization, as well as a hybrid combination of the two. In order
to activate these capabilities, the code needs to be compiled with OpenMP
and/or MPI support (``--omp`` and/or ``--mpi`` options to the CMake ``setup``
script, see :ref:`Installation` instructions).


Shared memory OpenMP
++++++++++++++++++++

For the shared memory part, the program will automatically pick up the number
of threads from the environment variable ``OMP_NUM_THREADS``. If this variable
is *not* set it will usually default to the maximum available. So, to run
the code on 16 threads (all sharing the same physical memory space)::

    $ OMP_NUM_THREADS=16 mrchem h2o


Distributed memory MPI
++++++++++++++++++++++

In order to run a program in an MPI parallel fashion, it must be executed with an MPI
launcher like ``mpirun``, ``mpiexec``, ``srun``, etc. Note that it is only
the main executable ``mrchem.x`` that should be launched in parallel, **not**
the ``mrchem`` input parser script. This can be achieved *either* by running
these separately in a dry-run (here two MPI processes)::

    $ mrchem --dryrun h2o
    $ mpirun -np 2 mrchem.x h2o.json

*or* in a single command by passing the launcher string as argument to the
parser::

    $ mrchem --launcher="mpirun -np 2" h2o

This string can contain any argument you would normally pass to ``mpirun``
as it will be literally prepended to the ``mrchem.x`` command when the
``mrchem`` script executes the main program.


.. hint::

    For best performance, it is recommended to use shared memory *within*
    each `NUMA <https://en.wikipedia.org/wiki/Non-uniform_memory_access>`_
    domain (usually one per socket) of your CPU, and MPI across NUMA domains and
    ultimately machines. Ideally, the number of OpenMP threads should be
    between 8-20. E.g. on hardware with two sockets of 16 cores each, use
    OMP_NUM_THREADS=16 and scale the number of MPI processes by the size
    of the molecule, typically one process per ~5 orbitals or so (and
    definitely not *more* than one process per orbital).


Job example (Betzy)
+++++++++++++++++++

This job will use 4 compute nodes, with 12 MPI processes on each, and the MPI
process will use up to 15 OpenMP threads. 4 MPI process per node are used for
the "Bank". The Bank processes are using only one thread, therefore there is
in practice no overallocation. It is however important that bank_size is set
to be at least 4*4 = 16 (it is by default set, correctly, to one third of total
MPI size, i.e. 4*12/3=16).
It would also be possible to set 16 tasks per node, and set the bank size
parameter accordingly to 8*4=32.
The flags are optimized for OpenMPI (foss) library on Betzy.

.. literalinclude:: betzy_example.job

``--rank-by node``
  Tells the system to place the first MPI rank on the first node, the second MPI
  rank on the second node, until the last node, then start at the first node again.

``--map-by socket``
  Tells the system to map (group) MPI ranks according to socket before distribution
  between nodes. This will ensure that for example two bank cores will access
  different parts of memory and be placed as the 16th thread of a numa group.

``--bind-to numa``
  Tells the system to bind cores to one NUMA (Non Uniform Memory Access) group.
  On Betzy memory configuration groups cores by groups of 16, with cores in the same
  group having the same access to memory (other cores will have access to that part
  of the memory too, but slower).
  That means that a process will
  only be allowed to use one of the 16 cores of the group. (The operating system may
  change the core assigned to a thread/process and, without precautions, it may be
  assigned to any other core, which would result in much reduced performance). The 16
  cores of the group may then be used by the threads initiated by that MPI process.

``--oversubscribe``
  To tell MPI that it is should accept that the number of MPI processes times
  the number of threads is larger than the number of available cores.

More examples can be found in the `mrchem-examples <https://github.com/MRChemSoft/mrchem-examples>`_
repository on GitHub.

Parallel pitfalls
-----------------

.. warning::

    Parallel program execution is not a black box procedure, and the behavior and
    efficiency of the run depends on several factors, like hardware configuration,
    operating system, compiler type and flags, libraries for OpenMP and MPI, type
    of queing system on a shared cluster, etc. Please make sure that the program
    runs correctly on *your* system and is able to utilize the computational
    resources before commencing production calculations.

Typical pitfalls for OpenMP
+++++++++++++++++++++++++++

- Not compiling with correct OpenMP support.
- Not setting number of threads correctly.
- **Hyper-threads:** the round-robin thread distribution might fill all
  hyper-threads on each core before moving on to the next physical core.
  In general we discourage the use of hyper-threads, and recommend a single
  thread per physical core.
- **Thread binding:** all threads may be bound to the same core, which means you
  can have e.g. 16 threads competing for the limited resources available on
  this single core (typically two hyper-threads) while all other cores are
  left idle.

Typical pitfalls for MPI
++++++++++++++++++++++++

- Not compiling with the correct MPI support.
- Default launcher options might not give correct behavior.
- **Process binding:** if a process is bound to a core, then all its spawned
  threads will also be bound to the same core. In general we recommend binding
  to socket/NUMA.
- **Process distribution:** in a multinode setup, all MPI processes might land
  on the same machine, or the round-robin procedure might count each core as
  a separate machine.

How to verify a parallel MRChem run
+++++++++++++++++++++++++++++++++++

- In the printed output, verify that MRCPP has actually been compiled with
  correct support for MPI and/or OpenMP::

    ----------------------------------------------------------------------

    MRCPP version         : 1.2.0
    Git branch            : master
    Git commit hash       : 686037cb78be601ac58b
    Git commit author     : Stig Rune Jensen
    Git commit date       : Wed Apr 8 11:35:00 2020 +0200

    Linear algebra        : EIGEN v3.3.7
    Parallelization       : MPI/OpenMP

    ----------------------------------------------------------------------

- In the printed output, verify that the correct number of processes and
  threads has been detected::

    ----------------------------------------------------------------------

     MPI processes         :      (no bank)                             2
     OpenMP threads        :                                           16
     Total cores           :                                           32

    ----------------------------------------------------------------------

- Monitor your run with ``top`` to see that you got the expected number of
  ``mrchem.x`` processes (MPI), and that they actually run at the expected
  CPU percentage (OpenMP)::

    PID   USER      PR  NI    VIRT    RES    SHR S   %CPU  %MEM     TIME+ COMMAND
    9502  stig      25   5  489456 162064   6628 R 1595,3   2,0   0:14.50 mrchem.x
    9503  stig      25   5  489596 162456   6796 R 1591,7   2,0   0:14.33 mrchem.x

- Monitor your run with ``htop`` to see which core/hyper-thread is being used
  by each process. This is very useful to get the correct binding/pinning of
  processes and threads. In general you want one threads per core, which means
  that every other hyper-thread should remain idle. In a hybrid MPI/OpenMP
  setup it is rather common that each MPI process becomes bound to a single
  core, which means that all threads spawned by this process will occupy the
  same core (possibly two hyper-threads). This is then easily detected with
  ``htop``.

- Perform dummy executions of your parallel launcher (``mpirun``, ``srun``, etc)
  to check whether it picks up the correct parameters from the resource manager
  on your cluster (SLURM, Torque, etc). You can then for instance report
  bindings and host name for each process::

    $ mpirun --print-rank-map hostname

  Play with the launcher options until you get it right. Note that Intel and
  OpenMPI have slightly different options for their ``mpirun`` and usually
  different behavior. Beware that the behavior can also change when you move
  from single- to multinode execution, so it is in general not sufficient to
  verify you runs on a single machine.

- Perform a small scaling test on e.g. 1, 2, 4 processes and/or 1, 2, 4 threads
  and verify that the total computation time is reduced as expected (don't
  expect 100% efficiency at any step).
