# Contributing to chemfiles

:+1::tada: First off, thanks for taking the time to contribute to chemfiles! :tada::+1:

If you want to contribute but feel a bit lost, do not hesitate to
[contact][Gitter] us and ask your questions! We will happily mentor you through
your first contributions.

## Area of contributions

The first and best way to contribute to chemfiles is to use it and advertise it
to other potential users. Other than that, you can help with:

- the documentation: correcting typos, improving the sentences to make things
  more clear;
- bug fixes and improvement to existing C++ code;
- implementing new formats;
- improving the C interface;
- and many more …

All these contributions are very welcome. We accept contributions either via
Github pull request (have a look [here][PR] for Github model of pull request);
or you can send patches by email at luthaf@luthaf.fr.

If you want to work on the code and pick something easy to get started, have a
look at the [easy issues][easy-issues].


## Bug reports and feature requests

Bug and feature requests should be reported as [Github issue][issue]. For bugs,
you should provide information so that we can reproduce it: what did you try?
What did you expect? What happened instead? Please provide any useful code
snippet or input file with your bug report.

If you want to add a new feature to chemfiles, please create an [issue] so that
we can discuss it, and you have more chances to see your changes incorporated.

### Code contribution check-list

Every item in this list is explained in the next section

- [ ] Fork chemfiles;
- [ ] Create a local branch;
- [ ] Add code / correct typos / ...;
    - [ ] Add new tests with your code;
    - [ ] Add documentation for your code;
- [ ] Check that the tests still pass;
- [ ] Push to Github;
- [ ] Create a Pull-Request;
- [ ] Discuss your changes with the reviewers;
- [ ] Have your code merged in chemfiles
- [ ] Celebrate! :tada: :cake: :tada:

### Contribution tutorial

In this small tutorial, you should replace `<angle brackets>` as needed. If
anything is unclear, please [ask][Gitter] for clarifications! There are no dumb
questions.

---

Start by [forking chemfiles][fork], and then clone and build your fork locally
by running:

```bash
git clone https://github.com/<YOUR USERNAME>/chemfiles
cd chemfiles
mkdir build
cd build
cmake ..
make
```

Then create a new branch for your changes

```bash
git checkout -b <new-branch>
```

Do your changes, containing the documentation (in `doc`) and associated unit
tests (in `tests`) for new code.

Then, run all the test and check that they still pass:

```bash
cmake -DCHFL_BUILD_TESTS=ON .
make
ctest
```

Finally, you can push your code to Github, and create a [Pull-Request][PR] to
the `chemfiles/chemfiles` repository.

```bash
git commit  # ask for help if you don't know how to use git
git push -u origin <new-branch>
```

[Gitter]: https://gitter.im/chemfiles/chemfiles
[PR]: https://help.github.com/articles/using-pull-requests/
[easy-issues]: https://github.com/chemfiles/chemfiles/labels/Help%20wanted
[fork]: https://help.github.com/articles/fork-a-repo/
[issue]: https://github.com/chemfiles/chemfiles/issues/new
## Chemfiles: a library for reading and writing chemistry files

[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](http://chemfiles.org/chemfiles/)
[![Build Status](https://img.shields.io/travis/chemfiles/chemfiles/master.svg)](https://travis-ci.org/chemfiles/chemfiles)
[![Code Coverage](http://codecov.io/github/chemfiles/chemfiles/coverage.svg?branch=master)](http://codecov.io/github/chemfiles/chemfiles?branch=master)
[![Gitter](https://badges.gitter.im/chemfiles/chemfiles.svg)](https://gitter.im/chemfiles/chemfiles)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3653157.svg)](https://doi.org/10.5281/zenodo.3653157)

Chemfiles is a high-quality library for reading and writing trajectory files
created by computational chemistry simulations programs. To help you access
information (atomic positions, velocities, names, topology, etc.) about these
files, Chemfiles provides a **simple** and **unified** interface to a variety of
file formats.

- **unified**: the same code will work with all supported formats;
- **simple**: the interface is easy to use and extensively documented.

You can use Chemfiles to conduct post-processing analysis and extract physical
information about the systems you're simulating, to convert files from one
format to another, to write trajectories with your own simulation software, and
anything that requires reading or writing the file formats used in computational
chemistry.

Chemfiles is used in multiple scientific software
- [cfiles](https://github.com/chemfiles/cfiles) provides ready-to-use analysis
  algorithms simulations trajectories as a command line tool;
- [lemon](https://github.com/chopralab/lemon) is a framework for rapidly mining
  structural information from the Protein Data Bank;
- [lumol](https://github.com/lumol-org/lumol) is a prototype of universal
  extensible molecular simulation engine, supporting both molecular dynamics
  and Metropolis Monte Carlo simulations;
- [ANA](https://ana.run/) detects cavities, calculates their volume and their
  flexibility in macromolecular structures and molecular dynamics trajectories;

This repository contains the core of the chemfiles library — written in C++11,
with a C99 interface. You can also use chemfiles from other languages: [Python
2&3](https://github.com/chemfiles/chemfiles.py),
[Fortran](https://github.com/chemfiles/chemfiles.f03),
[Rust](https://github.com/chemfiles/chemfiles.rs), and
[Julia](https://github.com/chemfiles/chemfiles.jl).

## Quick Links

- [Is chemfiles for you?](#is-chemfiles-for-you)
- [Main features of chemfiles](#chemfiles-features)
- [Contact / Contribute / Cite](#contact-contribute-cite)
- [Getting Started](#getting-started)
- [Supported File Formats](http://chemfiles.org/chemfiles/latest/formats.html)
- [Full documentation](http://chemfiles.org/chemfiles/)
- Documentation for using Chemfiles from various languages:
    - [Python 2 and 3](http://chemfiles.org/chemfiles.py/)
    - [Fortran](http://chemfiles.org/chemfiles.f03/)
    - [C and C++](http://chemfiles.org/chemfiles/)
    - [Julia](http://chemfiles.org/Chemfiles.jl/)
    - [Rust](http://chemfiles.org/chemfiles.rs/)

## Is chemfiles for you?

You might want to use chemfiles if any of these points appeals to you:

- you don't want to spend time writing and debugging a file parser;
- you use binary formats because they are faster and take up less disk space;
- you write analysis algorithms and want to read more than one trajectory
  format;
- you write simulation software and want to use more than one format for input
  or output.

There are [other libraries](http://chemfiles.org/chemfiles/latest/others.html)
doing the roughly the same job as chemfiles, have a look at them if chemfiles is
not for you. Here we also say why we could not use them instead of creating a
new library.

- [OpenBabel](https://openbabel.org/wiki/Main_Page) is a C++ library providing
  convertions between more than 110 formats. It is more complex than chemfiles,
  and distributed under the GPL license.
- [VMD molfile plugins](http://www.ks.uiuc.edu/Research/vmd/) are a collection
  of plugins witten in C and C++ used by VMD to read/write trajectory files.
  They do not support a variable number of atoms in a trajectory.
- [MDTraj](http://mdtraj.org/latest/), [MDAnalyis](http://www.mdanalysis.org/),
  [cclib](https://cclib.github.io/) are Python libraries providing analysis and
  read capacities for trajectories. Unfortunely, they are only usable from
  Python.

## Chemfiles Features

- Reads both text (XYZ, PDB, ...) and binary (NetCDF, TNG, ...) file formats;
- Transparently read and write compressed files (`.gz`, `.xz` and `.bz2`);
- Filters atoms with a rich selection language, including constrains on
  multiple atoms;
- Supports non-constant numbers of atoms in trajectories;
- Easy-to-use programming interface in Python, C++, C, Fortran 95, Julia and
  Rust;
- Cross-platform and usable from Linux, OS X and Windows;
- Open source and freely available (3-clauses BSD license);

## Contact / Contribute / Cite

Chemfiles is free and open source. Your [contributions](Contributing.md) are
always welcome!

If you have questions or suggestions, or need help, please open an [issue] or
join us on our [Gitter] chat room.

If you are using Chemfiles in a published scientific study, please cite us using
the following DOI: https://doi.org/10.5281/zenodo.3653157.

## Getting Started

Here, we'll help you get started with the C++ and C interface. If you want to
use Chemfiles with another language, please refer to the corresponding
documentation.

### Installing Compiled Packages

We provide compiled packages of the latest Chemfiles release for Linux
distributions. You can use your package manager to download them
[here][OSB-download].

We also provide conda packages in the `conda-forge` community channel for Linux
and OS X. This package provides the C++, C and Python interfaces. Install the conda package by running:

```
conda install -c conda-forge chemfiles
```

Find more information about pre-compiled packages in the [documentation][install].

### Building from Source

You will need [cmake](http://cmake.org/) and a C++11 compiler.

```bash
git clone https://github.com/chemfiles/chemfiles
cd chemfiles
mkdir build
cd build
cmake ..
make
make install
```

### Usage Examples

This is what the interface looks like in C++:

```cpp
#include <iostream>
#include "chemfiles.hpp"

int main() {
    chemfiles::Trajectory trajectory("filename.xyz");

    auto frame = trajectory.read();
    std::cout << "There are " << frame.size() << " atoms in the frame" << std::endl;

    auto positions = frame.positions();
    // Do awesome science with the positions here !
}
```

## License

Guillaume Fraux created and maintains Chemfiles, which is distributed under the
[3 clauses BSD license](LICENSE). By contributing to Chemfiles, you agree to
distribute your contributions under the same license. 

Chemfiles depends on multiple external libraries, which are distributed under [their
respective licenses](external/README.md). All external libraries licenses should
be compatible with chemfiles's 3 clauses BSD. One notable execption depending on
your use case is [Gemmi](https://gemmi.readthedocs.io) which is distributed
under the Mozilla Public License version 2. You can use `CHFL_DISABLE_GEMMI=ON`
CMake flag to remove this dependency.

The [AUTHORS](AUTHORS) file lists all contributors to Chemfiles. Many thanks to
all of them!

[Gitter]: https://gitter.im/chemfiles/chemfiles
[issue]: https://github.com/chemfiles/chemfiles/issues/new
[install]: http://chemfiles.org/chemfiles/latest/installation.html
[OpenSuseBuild]: https://build.opensuse.org/package/show/home:Luthaf/chemfiles
[OSB-download]: https://software.opensuse.org/download.html?project=home%3ALuthaf&package=chemfiles
# Change Log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](https://semver.org/).

## Next Release (current master)

### Deprecation and removals

- Remove support for configuration files (`chemfiles.toml`) and associated
  functions.

### New features

- properties are now sorted, and iterating over properties will always yield
  them sorted by the associated key.
- added `chemfiles::guess_format` and `chfl_guess_format` to get the format
  chemfiles would use for a given file based on its filename

### Changes in supported formats

- Added read and write support for Amber Restart (.ncrst) files.
- Added native read and write support for LAMMPS trajectory (.lammpstrj) files,
  replacing the VMD molfile implementation.
- Added read support for PSF files using VMD molfile plugin.

### Changes to the C API

## 0.10.0 (14 Feb 2021)

### Changes in supported formats

- Added read and write support for CIF format.
- Added read and write support for extended XYZ. The XYZ format now default to
  extended output, and read extended files. Extended XYZ allow storing unit
  cell and arbitrary atomic properties.

### New features

- Added ability to read and write files directly in-memory. See
  `Trajectory::memory_reader`; `Trajectory::memory_writer`;
  `Trajectory::memory_buffer`; `chfl_trajectory_memory_reader`;
  `chfl_trajectory_memory_writer` and `chfl_trajectory_memory_buffer`.
- Added support for appending to gzip (.gz) compressed trajectories.
- Added support for sub-selection in numerical functions, for example
  `distance(#1, name O)`.
- Changed the UnitCell representation to the full cell matrix instead of
  a/b/c/alpha/beta/gamma
- Added `chemfiles::formats_list` function to get a list of formats
  and associated metadata.

### Changes to the C API

- Added missing `chfl_frame_clear_bonds` and `chfl_topology_clear_bonds`.
- Added `chfl_cell_from_matrix` and changed parameters to `chfl_cell`, removed
  `chfl_cell_triclinic`.
- Added `chfl_formats_list` function to get a list of formats and
  associated metadata.

## 0.9.3 (5 Feb 2020)

- Fix a bug in the PDB format where no atomic name/type was read from short
  ATOM/HETATM records.
- Fix a few bugs related to UnitCell Matrix and infinite UnitCell construction

## 0.9.2 (18 Dec 2019)

- When compiling chemfiles as a shared library, the dependencies symbols are
  now hidden. This should prevent clashes between say chemfiles's zlib and the
  system zlib.
- Cache sub-selection (the 'name O' in 'is_bonded(#1, name O)'), reducing
  selection evaluation time by a huge margin.

### Changes in supported formats

- Added read and write support for CML (Chemical Markup Language) files, a XML
  based format.
- Added a native implementation of XTC and TRR formats, replacing the VMD
  molfile version. The new code supports reading and writing files, should be
  faster and use less memory.
- Remove the ability to read frames from a trajectory that was opened in
  append mode. This mode is now write only.
- Added support for bzip2 (.bz2) compressed files when reading and writing

## 0.9.1 (13 Mar 2019)

- Fix a bug with memory allocation in the C API. The allocator did not remove
  pointers as soon as `chfl_free` was called, which leaded to an error when the
  system allocator re-used the pointers.

## 0.9.0 (18 Nov 2018)

- Direct reading and writing of compressed files. gzip and lzma (.xz) formats
  are supported.
- GROMACS .gro files now supported through custom implementation.
- Properties are now supported in the `Residue` class. They are accessed using
  `Residue::set` and `Residue::get`.
- The topology of residues is now automatically set using a lookup table for the
  PDB format.
- `Frame::guess_topology` was renamed to `Frame::guess_bonds`.
- The selection engine has been rewritten to add support for more complex
  selections:
  - it is now possible to use mathematical expressions in selections such as
    `x^2 - y ^2 < sqrt(z^2 + 25)`;
  - it is now possible to access geometrical properties in such mathematical
    expressions: `distance(#1, #2)`, `angle(#1, #2, #3)`, `dihedral(#1, #2, #3, #4)`, and `out_of_plane(#1, #2, #3, #4)` are supported;
  - it is now possible to add constrains on the topology of the system:
    `bonded(#1, #2)`, `is_angle(#1, #2, #3)`, `is_dihedral(#1, #2, #3, #4)`,
    and `improper(#1, #2, #3, #4)` are supported;
  - the topology constrains support sub-selections: instead of checking is
    `#1` and `#2` are bonded, one can check if `#1` is bonded to any atom
    matching a selection, for example `name O` with `bonded(#1, name O)`.
  - When using numbers as atomic names/types, they must now be inside double
    quotes (`name "45"`). This also allows for more exotic atomic names
    (`name "名"`).
  - Atomic properties can be checked, using the `[property] == Ow` syntax for
    string properties, `[property] == 2.3` for numeric properties and
    `[property]` for boolean properties.
- There is only one constructor for the `Frame` class: `Frame(UnitCell cell = UnitCell())`. The constructor taking a topology can be replaced with calls to
  `Frame::add_atom` and `Frame::add_bond`.
- Chemfiles will now read configuration from `.chemfiles.toml` or
  `chemfiles.toml` instead of `.chemfilesrc`
- Added `Trajectory::path` to get the file path used to create a trajectory
- Renamed `Property::get_kind` to `Property::kind`
- Added `Atom::properties`; `Frame::properties`; and `Residue::properties` to
  allow iteration over all the properties in an Atom/Frame/Residue.

### Changes in supported formats

- Added `MarcoMolecule Transmission Format (MMTF)` support, reading via mmtf-cpp.
- Added `Structure-Data File (SDF)` support, reading and writing.
- Added `Cambridge Structure Search and Retrieval (CSSR)` support, reading and writing.
- `LAMMPS Data` format now support triclinic unit cells.

### Changes to the C API

- Added `chfl_residue_get_property` and `chfl_residue_set_property` to provide
  access to residue properties.
- `chfl_frame_guess_topology` was renamed to `chfl_frame_guess_bonds`.
- Function accessing atoms/cell/residue/topology inside a frame/topology no
  longer make a copy. This allows for direct reading and writing inside the
  containing frame/topology.
- Added `chfl_trajectory_path` to get the file path used to create a trajectory
- Added `chfl_{atom,frame,residue}_properties_count` and
  `chfl_{atom,frame,residue}_list_properties` to list all properties in an
  Atom/Frame/Residue
- Replaced `chfl_*_free` by an unique `chfl_free` function

## 0.8 (14 Dec 2017)

### New features

- Change the license to the 3-clauses BSD license.
- Chemfiles will now read configuration files (by default in `.chemfilesrc`),
  and use the configuration data to rename atomic types to make sure they match
  element names. The `chemfiles::add_configuration` function can be used to add
  additional configuration files.
- Reading a `Frame` (with `Trajectory::read` or `Trajectory::read_step`) will
  now set the frame step.
- The `Atom` faillible methods (`atomic_number`, `vdw_radius`, `covalent_radius`
  and `full_name`) returns `optional<T>` instead of `T`.
- Functions taking an atomic index parameter can now throw `OutOfBounds` errors
  if the index is out of bounds.
- `Topology::append` is now called `Topology::add_atom`
- `Topology::natoms` and `Frame::natoms` are now called `Topology::size` and
  `Frame::size`
- `Topology::residue` is now called `Topology::residue_for_atom`
- Added `Frame::distance`, `Frame::angle`, `Frame::dihedral` and
  `Frame::out_of_plane` to get geometric information on the system, accounting
  for periodic boundary conditions.
- Added a `Property` class to store arbitrary properties in `Frame` and `Atom`.
- Added support for improper dihedral angles in `Topology`.
- `chemfiles::add_configuration` and `chemfiles::set_warning_callback` are now
  thread safe, and will block upon concurrent usage.
- `UnitCell::matricial` is renamed to `UnitCell::matrix`.
- The `UnitCell::shape` setter is renamed to `UnitCell::set_shape`.
- The `Trajectory::close` function can be used to close a trajectory and
  synchronize any buffered content with the storage.
- Some of the topology functions are now accsible directly on the frame:
  `Frame::add_bond`, `Frame::remove_bond`, `Frame::clear_bonds`,
  `Frame::add_residue` and `operator[]`. The non const version of
  `Frame::topology` is removed.

### Changes in supported formats

- Amber NetCDF format is now activated by default, by embedding the netcdf
  library in chemfiles.
- Added `LAMMPS Data` format, reading and writing [LAMMPS data files].
- Added `Tinker` format, reading and writing Tinker XYZ file format.
- Added `MOL2` format, reading mol2 files using VMD molfiles plugin.
- Added `Molden` format, reading molden files using VMD molfiles plugin.

[lammps data files]: https://lammps.sandia.gov/doc/read_data.html

### Changes to the C API

- Added `chfl_add_configuration` to add more configuration files.
- Renamed `chfl_vector_t` to `chfl_vector3d`, `chfl_match_t` to `cfl_match`; and
  `chfl_cell_shape_t` to `chfl_cellshape`.
- `chfl_atom_atomic_number`, `chfl_atom_vdw_radius` and
  `chfl_atom_covalent_radius` all returns 0 instead of -1 if the atom does not
  have a known value for this property. This allow `chfl_atom_atomic_number` to
  take a `uint64_t*` parameter instead of an `int64_t*`, following all the other
  functions in the C API.
- Added `CHFL_OUT_OF_BOUNDS` and `CHFL_PROPERTY_ERROR` variants to `chfl_status`
- Added `chfl_frame_distance`, `chfl_frame_angle`, `chfl_frame_dihedral`,
  `chfl_frame_out_of_plane` and `chfl_cell_wrap` to work with periodic boundary
  conditions.
- `chfl_residue` does not take the optional residue id as parameter, instead you
  should use `chfl_residue_with_id`.
- Added `chfl_residue_atoms` to get the list of atoms in a residue.
- Added `chfl_topology_impropers` and `chfl_topology_impropers_count` functions.
- Added `CHFL_PROPERTY` and related functions.
- `chfl_add_configuration` and `chfl_set_warning_callback` are now thread safe,
  and will block upon concurrent usage.
- Added `chfl_frame_add_bond`, `chfl_frame_remove_bond`, and
  `chfl_frame_add_residue`.

### Deprecation and removals

- `Topology::isbond`, `Topology::isangle`, `Topology::isdihedral`, and the
  corresponding C functions `chfl_topology_isbond`, `chfl_topology_isangle`
  `chfl_topology_isdihedral` are removed.

## 0.7 (25 Feb 2017)

### New features

- Add a public `Residue` class to C++ and C API to represent residue data.
  Residues are groups of atoms bonded together, which may or may not correspond
  to molecules; and are often used for bio-molecules.
- Add the `resname` and `resid` selector, to select atoms based on their
  residue.
- Account for the difference between the atom name ("H1") and atom type ("H")
  in some formats (PDB, TNG, ...). This introduces the `Atom::type` member
  function and the `chfl_atom_type` C API function.
- Add the `type` selector, to select atoms based on their type.
- Add "Frame::add_atom" function to add an atom and the corresponding position
  (and velocity) data to a frame, and the C API `chfl_frame_add_atom` function.
- Rename `UnitCell::type` to `UnitCell::shape`. This also affect
  `chfl_cell_shape_t`, `chfl_cell_shape`, and `chfl_cell_set_shape`.
- All the floating point data uses doubles instead of floats. This concerns
  atomic data, positions and velocities.
- Add "Selection::string" function and the corresponding `chfl_selection_string`
  to get the string used to build a selection.
- Selection variables uses the `#3` syntax instead of the `$3` syntax to allow
  passing selection string as shell arguments.
- Add `Frame::remove` and `chfl_frame_remove` to remove an atom in a frame.
- Allow to use chemfiles with Cmake `find_package`.

### Changes in supported formats

- Add read support for TNG files, an new portable and compressed binary format
  used by GROMACS.

### Changes to the C API

- All the integers at C boundary have a fixed size, most of the time using
  `uint64_t`.
- Add missing `chfl_topology_resize` function to C API.
- C API functions taking three lengths/angles now take a `double[3]` parameter
  instead.
- Rename `chfl_topology_are_linked` to `chfl_topology_residues_linked`.
- Rename `chfl_topology_append` to `chfl_topology_add_atom`.
- Remove `chfl_strerror`, as it is redundant with `chfl_last_error`.
- Merge `chfl_trajectory_set_topology_file` and `chfl_trajectory_set_topology_with_format`
  into `chfl_trajectory_topology_file`.
- The `chfl_frame` function no longer take the frame size as argument. It always
  creates an empty frame, that you can resize using `chfl_frame_resize`.
- `chfl_selection_evalutate` was a typo, it is renamed to `chfl_selection_evaluate`.

### Deprecation and removals

- Remove the `Atom::type` enum from C and C++ API.
- Remove the `Trajectory::sync` and the `chfl_trajectory_sync` functions.
  To ensure that all content of a file is written to the disk, the user need to
  close it.
- Remove the `Logger` and all the `chfl_log*` functions. Instead, the users can
  use `chemfiles::set_warning_callback` or `chfl_set_warning_callback` to set a
  global callback to call on warnings events.

## 0.6 (1 July 2016)

- Improve the selection language to allow selecting multiple atoms at once. For
  example, `"pairs: name($1) H and mass($2) > 5"` will select all pairs of atoms
  where the first atom name is `'H'` and the second atom mass is bigger than 5.
  - The implemented modes for selections are `one`, `atoms`, `two`, `pairs`,
    `three`, `four`, `bonds`, `angles` and `dihedrals`;
  - The `Selection` class is now directly exposed to the C API, as
    `CHFL_SELECTION*`. The `chfl_frame_selection` function is replaced by the
    `chfl_selection`, `chfl_selection_size`, `chfl_selection_evalutate`,
    `chfl_selection_matches` and `chfl_selection_free` functions, and the
    `chfl_match_t` helper struct.
- Add the `chfl_clear_errors` function, to cleanup the error state of the C API.
- Molfiles plugins are now incorporated in the Chemfiles library, and no longer
  distributed as shared libraries. The `CHEMFILES_PLUGINS` environment variable
  is a no-op.
- The caching of angles and dihedrals is now an implementation detail. That
  means that `Topology::recalculate` is gone, and that `Frame::guess_topology`
  and `chfl_frame_guess_topology` do not take a boolean parameter anymore.
- The opening mode is now a `char` instead of a string in `Trajectory`
  constructor, `chfl_trajectory_open`, and `chfl_trajectory_with_format`.
- Remove `operator<<` and `operator>>` for `Trajectory`. Users should use
  `Trajectory::read` and `Trajectory::write`
- Users can now specify the format when reading the topology associated with a
  trajectory from a file. The `chfl_trajectory_set_topology_with_format`
  function can be used to do so from the C API.
- The `chfl_atom_from_{frame,topology}` function now return `NULL` in case of
  out-of-bound access.

## 0.5 (19 Feb 2016)

- The C API now provide a direct view into the `positions` and `velocities`
  arrays. This remove the need for copy and separated getter
  (`chfl_frame_{position,velocities}`) and setter
  (`chfl_frame_{position,velocities}_set`) function. This also force usage of
  `chfl_frame_add_velocities` to add velocity data to a frame, and
  `chfl_frame_resize` to change the size of the frame.
- Add constants for error codes in C API. The following macro are defined:
  `CHFL_SUCCESS`, `CHFL_MEMORY_ERROR`, `CHFL_FILE_ERROR`, `CHFL_FORMAT_ERROR`,
  `CHFL_GENERIC_ERROR`, `CHFL_CXX_ERROR`.
- Add the `chfl_version` function in C API.
- Add a small selection language _a la_ VMD, allowing to select atoms matching
  a selection string like `"name H and x > 4"`. This is exposed to C++ with the
  public `Selection` class, and to C with the `chfl_frame_selection` function.
- Remove the periodicity handling from `UnitCell`. It was not implemented in
  boundaries conditions. The corresponding function where removed from the C
  API.
- Rename all setter function from `void xxx(const XXX& value)` to
  `void set_xxx(const XXX& value)` in C++ API.
- It is now possible to provide a callback for logging from the C API. The
  `chfl_log_stdout`, `chfl_log_silent` and `chfl_log_callback` function where
  added to the C API.

## 0.4 (30 Oct 2015)

- Chemharp can now be compiled as a static library! This should allow for easier
  embedding in external code, and easier distribution of binaries.
- Add a `chemfiles::Trajectory::sync` method to sync any buffered operation with
  the disk. The `chfl_trajectory_sync` function exposes it to the C API.
- Add a Rust binding
- Rewrite the Python binding to use ctypes. The same code can be used with
  Python 2 & 3, and with all numpy versions.
- Easier Python and Julia binding installation, using conda binary packaging.
- All the bindings now live on their own repository:
  - [Fortran](https://github.com/Luthaf/Chemharp.f03)
  - [Python](https://github.com/Luthaf/Chemharp.py)
  - [Julia](https://github.com/Luthaf/Chemharp.jl)
  - [Rust](https://github.com/Luthaf/Chemharp.rs)
- The library is now continuously tested on Visual Studio
- Various bug fixes and code improvements
- Renamed the library to Chemfiles.

## 0.3 (3 Aug 2015)

- Julia binding
- Initial Windows support, with both MSVC and mingw
- Add a binary frontend called `chrp`, implementing some analysis algorithms.
  For more information, see the [specific repository](https://github.com/Luthaf/chrp).

## 0.2 (31 May 2015)

- Add basic geometrical operations on vectors and implement basic periodic boundaries condition with the `UnitCell::wrap` function;
- Use VMD Molfiles plugins as a format provider to read trajectories. The following formats are added through Molfiles:
  _ PDB;
  _ GROMACS gro;
  _ GROMACS xtc;
  _ GROMACS trj;
  _ GROMACS trr;
  _ CHARMM dcd;

## 0.1 (16 May 2015)

Initial release. See the documentation for the full API.

Chemharp is usable from four languages:

- C++;
- C;
- Fortran;
- Python;

The following formats are supported:

- XYZ;
- AMBER NetCDF;
# Chemfiles documentation

This directory contains the documentation for the chemfiles library. This
documentation is divided in multiple sections:

- High level overview;
- C++ interface reference;
- C interface reference;
- Developers documentation.

The documentation is written using [sphinx-doc] and [RestructuredText]. It is
automatically compiled every time a modification is added on Github, and
deployed at http://chemfiles.org/chemfiles/.

Part of the documentation (the interface reference) is extracted from the source
code. The source code contains special comments following [Doxygen] conventions,
which are extracted and included in the documentation using [breathe].

To modify the documentation, you will need to modify either the `.rst` files in
the `doc/src` directory; or the [Doxygen] comments in header files from the
`include` directory.

## Building the documentation locally

First, get the source code of chemfiles (if you don't already have it):

```
git clone https://github.com/chemfiles/chemfiles
cd chemfiles
```

Then install [Doxygen] and Python using your favorite method. Finally, install
the required python packages by running:

```
pip install sphinx
pip install -r doc/requirements.txt
```

You can now build the documentation locally by running:

```
mkdir build
cd build
cmake -DCHFL_BUILD_DOCUMENTATION=ON ..
make doc_html
```

The documentation will be in the `build/doc/html` directory.

[sphinx-doc]: http://www.sphinx-doc.org/
[RestructuredText]: http://www.sphinx-doc.org/en/stable/rest.html
[Doxygen]: http://doxygen.org/
[breathe]: http://breathe.readthedocs.io/
# Example directory

This directory contains usage examples of the chemfiles library, for the C and C++
languages. Similar examples can be found in the corresponding repository for
all the bindings to chemfiles.

All of these example are compiled with the latest release to ensure that they are
up-to-date.

The files `indexes.*` present an example of reading a single frame from a
trajectory, and doing some analysis on it.

The files `rmsd.*` present an example of reading multiple frames from a
trajectory, and doing some analysis on the trajectory of a particle.

The files `convert.*` present an example of setting frame properties, and then
writing these frames to a trajectory.
# External libraries used in chemfiles

This directory contains external libraries used in chemfiles. Each of these
library is distributed under it's own license, which are all included in the
archive and reproduced at the end of this file.

## Updating

We maintain a set of patches to these libraries, in there respective forks, in
the `chemfiles` branch. These patches mainly concern integration of the code
with chemfiles build system.

- fmt: https://github.com/chemfiles/fmt
- zlib: https://github.com/chemfiles/zlib
- bzip2: https://github.com/chemfiles/bzip2
- liblzma: https://github.com/chemfiles/lzma
- pugixml: https://github.com/chemfiles/pugixml
- molfiles: https://github.com/chemfiles/molfiles
- TNG: https://github.com/chemfiles/tng
- NetCDF: https://github.com/chemfiles/netcdf-c
- mmtf-cpp: https://github.com/chemfiles/mmtf-c
- gemmi: https://github.com/chemfiles/gemmi

To update a library, update the corresponding repository, and then regenerate
the archive to be included in this directory with

```bash
# <lib> is the library name: tng, molfiles, netcdf, ...
git archive HEAD -9 --prefix=<lib>/ -o <lib>.tar.gz
```

For gemmi, since the library is header only you should use this command to
create the archive:

```bash
git archive HEAD -9 --prefix=gemmi/ -o gemmi.tar.gz include
```

## Licences

### fmt

Copyright (c) 2012 - 2016, Victor Zverovich

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### VMD molfiles

University of Illinois Open Source License
Copyright 2003 Theoretical and Computational Biophysics Group,
All rights reserved.

Developed by: Theoretical and Computational Biophysics Group
University of Illinois at Urbana-Champaign
http://www.ks.uiuc.edu/

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the Software), to deal with
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to
do so, subject to the following conditions:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimers.

Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimers in the documentation
and/or other materials provided with the distribution.

Neither the names of Theoretical and Computational Biophysics Group,
University of Illinois at Urbana-Champaign, nor the names of its contributors
may be used to endorse or promote products derived from this Software without
specific prior written permission.

THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS WITH THE SOFTWARE.

### TNG

Copyright (c) 2012-2013, The GROMACS development team,
check out http://www.gromacs.org for more information.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
- Neither the name of the GROMACS development team nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## NetCDF

The NetCDF Copyright.

Copyright 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014,
University Corporation for Atmospheric Research/Unidata.

Portions of this software were developed by the Unidata Program at the
University Corporation for Atmospheric Research.

Access and use of this software shall impose the following obligations
and understandings on the user. The user is granted the right, without
any fee or cost, to use, copy, modify, alter, enhance and distribute
this software, and any derivative works thereof, and its supporting
documentation for any purpose whatsoever, provided that this entire
notice appears in all copies of the software, derivative works and
supporting documentation. Further, UCAR requests that the user credit
UCAR/Unidata in any publications that result from the use of this
software or in any product that includes this software, although this
is not an obligation. The names UCAR and/or Unidata, however, may not
be used in any advertising or publicity to endorse or promote any
products or commercial entity unless specific written permission is
obtained from UCAR/Unidata. The user also understands that
UCAR/Unidata is not obligated to provide the user with any support,
consulting, training or assistance of any kind with regard to the use,
operation and performance of this software nor to provide the user
with any updates, revisions, new versions or "bug fixes."

THIS SOFTWARE IS PROVIDED BY UCAR/UNIDATA "AS IS" AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL UCAR/UNIDATA BE LIABLE FOR ANY SPECIAL,
INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.

## mmft-cpp

The MIT license applies to all files unless otherwise noted.

Copyright 2017, University of California San Diego

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

**mmtf-cpp depends on msgpack-c**

## msgpack-c

Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

## zlib

(C) 1995-2017 Jean-loup Gailly and Mark Adler

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1.  The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software
    in a product, an acknowledgment in the product documentation would be
    appreciated but is not required.
2.  Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.
3.  This notice may not be removed or altered from any source distribution.

Jean-loup Gailly Mark Adler
jloup@gzip.org madler@alumni.caltech.edu

## liblzma

liblzma is in the public domain.

## libbzip2

This program, "bzip2", the associated library "libbzip2", and all
documentation, are copyright (C) 1996-2010 Julian R Seward. All
rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. The origin of this software must not be misrepresented; you must
   not claim that you wrote the original software. If you use this
   software in a product, an acknowledgment in the product
   documentation would be appreciated but is not required.

3. Altered source versions must be plainly marked as such, and must
   not be misrepresented as being the original software.

4. The name of the author may not be used to endorse or promote
   products derived from this software without specific prior written
   permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Julian Seward, jseward@acm.org
bzip2/libbzip2 version 1.0.7 of 27 June 2019

## pugixml

MIT License

Copyright (c) 2006-2019 Arseny Kapoulkine

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## gemmi

# Mozilla Public License Version 2.0

1. Definitions

---

1.1. "Contributor"
means each individual or legal entity that creates, contributes to
the creation of, or owns Covered Software.

1.2. "Contributor Version"
means the combination of the Contributions of others (if any) used
by a Contributor and that particular Contributor's Contribution.

1.3. "Contribution"
means Covered Software of a particular Contributor.

1.4. "Covered Software"
means Source Code Form to which the initial Contributor has attached
the notice in Exhibit A, the Executable Form of such Source Code
Form, and Modifications of such Source Code Form, in each case
including portions thereof.

1.5. "Incompatible With Secondary Licenses"
means

    (a) that the initial Contributor has attached the notice described
        in Exhibit B to the Covered Software; or

    (b) that the Covered Software was made available under the terms of
        version 1.1 or earlier of the License, but not also under the
        terms of a Secondary License.

1.6. "Executable Form"
means any form of the work other than Source Code Form.

1.7. "Larger Work"
means a work that combines Covered Software with other material, in
a separate file or files, that is not Covered Software.

1.8. "License"
means this document.

1.9. "Licensable"
means having the right to grant, to the maximum extent possible,
whether at the time of the initial grant or subsequently, any and
all of the rights conveyed by this License.

1.10. "Modifications"
means any of the following:

    (a) any file in Source Code Form that results from an addition to,
        deletion from, or modification of the contents of Covered
        Software; or

    (b) any new file in Source Code Form that contains any Covered
        Software.

1.11. "Patent Claims" of a Contributor
means any patent claim(s), including without limitation, method,
process, and apparatus claims, in any patent Licensable by such
Contributor that would be infringed, but for the grant of the
License, by the making, using, selling, offering for sale, having
made, import, or transfer of either its Contributions or its
Contributor Version.

1.12. "Secondary License"
means either the GNU General Public License, Version 2.0, the GNU
Lesser General Public License, Version 2.1, the GNU Affero General
Public License, Version 3.0, or any later versions of those
licenses.

1.13. "Source Code Form"
means the form of the work preferred for making modifications.

1.14. "You" (or "Your")
means an individual or a legal entity exercising rights under this
License. For legal entities, "You" includes any entity that
controls, is controlled by, or is under common control with You. For
purposes of this definition, "control" means (a) the power, direct
or indirect, to cause the direction or management of such entity,
whether by contract or otherwise, or (b) ownership of more than
fifty percent (50%) of the outstanding shares or beneficial
ownership of such entity.

2. License Grants and Conditions

---

2.1. Grants

Each Contributor hereby grants You a world-wide, royalty-free,
non-exclusive license:

(a) under intellectual property rights (other than patent or trademark)
Licensable by such Contributor to use, reproduce, make available,
modify, display, perform, distribute, and otherwise exploit its
Contributions, either on an unmodified basis, with Modifications, or
as part of a Larger Work; and

(b) under Patent Claims of such Contributor to make, use, sell, offer
for sale, have made, import, and otherwise transfer either its
Contributions or its Contributor Version.

2.2. Effective Date

The licenses granted in Section 2.1 with respect to any Contribution
become effective for each Contribution on the date the Contributor first
distributes such Contribution.

2.3. Limitations on Grant Scope

The licenses granted in this Section 2 are the only rights granted under
this License. No additional rights or licenses will be implied from the
distribution or licensing of Covered Software under this License.
Notwithstanding Section 2.1(b) above, no patent license is granted by a
Contributor:

(a) for any code that a Contributor has removed from Covered Software;
or

(b) for infringements caused by: (i) Your and any other third party's
modifications of Covered Software, or (ii) the combination of its
Contributions with other software (except as part of its Contributor
Version); or

(c) under Patent Claims infringed by Covered Software in the absence of
its Contributions.

This License does not grant any rights in the trademarks, service marks,
or logos of any Contributor (except as may be necessary to comply with
the notice requirements in Section 3.4).

2.4. Subsequent Licenses

No Contributor makes additional grants as a result of Your choice to
distribute the Covered Software under a subsequent version of this
License (see Section 10.2) or under the terms of a Secondary License (if
permitted under the terms of Section 3.3).

2.5. Representation

Each Contributor represents that the Contributor believes its
Contributions are its original creation(s) or it has sufficient rights
to grant the rights to its Contributions conveyed by this License.

2.6. Fair Use

This License is not intended to limit any rights You have under
applicable copyright doctrines of fair use, fair dealing, or other
equivalents.

2.7. Conditions

Sections 3.1, 3.2, 3.3, and 3.4 are conditions of the licenses granted
in Section 2.1.

3. Responsibilities

---

3.1. Distribution of Source Form

All distribution of Covered Software in Source Code Form, including any
Modifications that You create or to which You contribute, must be under
the terms of this License. You must inform recipients that the Source
Code Form of the Covered Software is governed by the terms of this
License, and how they can obtain a copy of this License. You may not
attempt to alter or restrict the recipients' rights in the Source Code
Form.

3.2. Distribution of Executable Form

If You distribute Covered Software in Executable Form then:

(a) such Covered Software must also be made available in Source Code
Form, as described in Section 3.1, and You must inform recipients of
the Executable Form how they can obtain a copy of such Source Code
Form by reasonable means in a timely manner, at a charge no more
than the cost of distribution to the recipient; and

(b) You may distribute such Executable Form under the terms of this
License, or sublicense it under different terms, provided that the
license for the Executable Form does not attempt to limit or alter
the recipients' rights in the Source Code Form under this License.

3.3. Distribution of a Larger Work

You may create and distribute a Larger Work under terms of Your choice,
provided that You also comply with the requirements of this License for
the Covered Software. If the Larger Work is a combination of Covered
Software with a work governed by one or more Secondary Licenses, and the
Covered Software is not Incompatible With Secondary Licenses, this
License permits You to additionally distribute such Covered Software
under the terms of such Secondary License(s), so that the recipient of
the Larger Work may, at their option, further distribute the Covered
Software under the terms of either this License or such Secondary
License(s).

3.4. Notices

You may not remove or alter the substance of any license notices
(including copyright notices, patent notices, disclaimers of warranty,
or limitations of liability) contained within the Source Code Form of
the Covered Software, except that You may alter any license notices to
the extent required to remedy known factual inaccuracies.

3.5. Application of Additional Terms

You may choose to offer, and to charge a fee for, warranty, support,
indemnity or liability obligations to one or more recipients of Covered
Software. However, You may do so only on Your own behalf, and not on
behalf of any Contributor. You must make it absolutely clear that any
such warranty, support, indemnity, or liability obligation is offered by
You alone, and You hereby agree to indemnify every Contributor for any
liability incurred by such Contributor as a result of warranty, support,
indemnity or liability terms You offer. You may include additional
disclaimers of warranty and limitations of liability specific to any
jurisdiction.

4. Inability to Comply Due to Statute or Regulation

---

If it is impossible for You to comply with any of the terms of this
License with respect to some or all of the Covered Software due to
statute, judicial order, or regulation then You must: (a) comply with
the terms of this License to the maximum extent possible; and (b)
describe the limitations and the code they affect. Such description must
be placed in a text file included with all distributions of the Covered
Software under this License. Except to the extent prohibited by statute
or regulation, such description must be sufficiently detailed for a
recipient of ordinary skill to be able to understand it.

5. Termination

---

5.1. The rights granted under this License will terminate automatically
if You fail to comply with any of its terms. However, if You become
compliant, then the rights granted under this License from a particular
Contributor are reinstated (a) provisionally, unless and until such
Contributor explicitly and finally terminates Your grants, and (b) on an
ongoing basis, if such Contributor fails to notify You of the
non-compliance by some reasonable means prior to 60 days after You have
come back into compliance. Moreover, Your grants from a particular
Contributor are reinstated on an ongoing basis if such Contributor
notifies You of the non-compliance by some reasonable means, this is the
first time You have received notice of non-compliance with this License
from such Contributor, and You become compliant prior to 30 days after
Your receipt of the notice.

5.2. If You initiate litigation against any entity by asserting a patent
infringement claim (excluding declaratory judgment actions,
counter-claims, and cross-claims) alleging that a Contributor Version
directly or indirectly infringes any patent, then the rights granted to
You by any and all Contributors for the Covered Software under Section
2.1 of this License shall terminate.

5.3. In the event of termination under Sections 5.1 or 5.2 above, all
end user license agreements (excluding distributors and resellers) which
have been validly granted by You or Your distributors under this License
prior to termination shall survive termination.

---

-                                                                      *
- 6.  Disclaimer of Warranty \*
- ------------------------- \*
-                                                                      *
- Covered Software is provided under this License on an "as is" \*
- basis, without warranty of any kind, either expressed, implied, or \*
- statutory, including, without limitation, warranties that the \*
- Covered Software is free of defects, merchantable, fit for a \*
- particular purpose or non-infringing. The entire risk as to the \*
- quality and performance of the Covered Software is with You. \*
- Should any Covered Software prove defective in any respect, You \*
- (not any Contributor) assume the cost of any necessary servicing, \*
- repair, or correction. This disclaimer of warranty constitutes an \*
- essential part of this License. No use of any Covered Software is \*
- authorized under this License except under this disclaimer. \*
-                                                                      *

---

---

-                                                                      *
- 7.  Limitation of Liability \*
- -------------------------- \*
-                                                                      *
- Under no circumstances and under no legal theory, whether tort \*
- (including negligence), contract, or otherwise, shall any \*
- Contributor, or anyone who distributes Covered Software as \*
- permitted above, be liable to You for any direct, indirect, \*
- special, incidental, or consequential damages of any character \*
- including, without limitation, damages for lost profits, loss of \*
- goodwill, work stoppage, computer failure or malfunction, or any \*
- and all other commercial damages or losses, even if such party \*
- shall have been informed of the possibility of such damages. This \*
- limitation of liability shall not apply to liability for death or \*
- personal injury resulting from such party's negligence to the \*
- extent applicable law prohibits such limitation. Some \*
- jurisdictions do not allow the exclusion or limitation of \*
- incidental or consequential damages, so this exclusion and \*
- limitation may not apply to You. \*
-                                                                      *

---

8. Litigation

---

Any litigation relating to this License may be brought only in the
courts of a jurisdiction where the defendant maintains its principal
place of business and such litigation shall be governed by laws of that
jurisdiction, without reference to its conflict-of-law provisions.
Nothing in this Section shall prevent a party's ability to bring
cross-claims or counter-claims.

9. Miscellaneous

---

This License represents the complete agreement concerning the subject
matter hereof. If any provision of this License is held to be
unenforceable, such provision shall be reformed only to the extent
necessary to make it enforceable. Any law or regulation which provides
that the language of a contract shall be construed against the drafter
shall not be used to construe this License against a Contributor.

10. Versions of the License

---

10.1. New Versions

Mozilla Foundation is the license steward. Except as provided in Section
10.3, no one other than the license steward has the right to modify or
publish new versions of this License. Each version will be given a
distinguishing version number.

10.2. Effect of New Versions

You may distribute the Covered Software under the terms of the version
of the License under which You originally received the Covered Software,
or under the terms of any subsequent version published by the license
steward.

10.3. Modified Versions

If you create software not governed by this License, and you want to
create a new license for such software, you may create and use a
modified version of this License if you rename the license and remove
any references to the name of the license steward (except to note that
such modified license differs from this License).

10.4. Distributing Source Code Form that is Incompatible With Secondary
Licenses

If You choose to distribute Source Code Form that is Incompatible With
Secondary Licenses under the terms of this version of the License, the
notice described in Exhibit B of this License must be attached.

## Exhibit A - Source Code Form License Notice

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

If it is not possible or desirable to put the notice in a particular
file, then You may include the notice in a location (such as a LICENSE
file in a relevant directory) where a recipient would be likely to look
for such a notice.

You may add additional accurate notices of copyright ownership.

## Exhibit B - "Incompatible With Secondary Licenses" Notice

This Source Code Form is "Incompatible With Secondary Licenses", as
defined by the Mozilla Public License, v. 2.0.
.. _c-tutorials:

C Tutorials
===========

This section present some hand-on tutorials to the chemfiles C API. These
examples do exaclty the same thing as the :ref:`C++ ones <cxx-tutorials>`, using
the C API instead of the C++ API. As such, they are much more verbose. They also
do not correclty check for error conditions. When using the C API, you should
alway check that the return value of the functions is ``CHFL_SUCCESS``.

All the code here is under the `CC-0 Universal Licence`_ which means that you
are free to do whatever you want with it (*i.e.* it is Public Domain code)

.. _CC-0 Universal Licence: https://creativecommons.org/publicdomain/zero/1.0/

Read a single frame
-------------------

In this tutorials we will read a frame from a trajectory, and print the indexes
of all the atom in the half-space ``x < 5``.

We start by including the headers we need: chemfiles header, ```stdio.h`` for
printing and ``stdlib.h`` for ``malloc``. All the public functions of chemfiles
are accessible through the ``chemfiles.h`` header, which is the only header
needed.

.. literalinclude:: ../../examples/c/indexes.c
   :language: c
   :lines: 4-6

Then we open a :cpp:type:`CHFL_TRAJECTORY` in read (``'r'``) and read the first
frame.  We need to allocate memory for the :cpp:type:`CHFL_FRAME` before calling
:cpp:func:`chfl_trajectory_read`.

.. literalinclude:: ../../examples/c/indexes.c
   :language: c
   :lines: 9-11
   :dedent: 4

We can now and get the positions of the atoms and the number of atoms in the
frame with the :cpp:func:`chfl_frame_positions` function:

.. literalinclude:: ../../examples/c/indexes.c
   :language: c
   :lines: 13-15
   :dedent: 4

Knowning the total number of atoms in the frame, we can allocate memory to store
the indices of the atoms with ``x < 5``:

.. literalinclude:: ../../examples/c/indexes.c
   :language: c
   :lines: 17
   :dedent: 4

Iterating through the atoms in the frame, we get the ones matching our
condition. We need to track the number of ``matched`` atoms to know where to add
them in the ``less_than_five`` array.

.. literalinclude:: ../../examples/c/indexes.c
   :language: c
   :lines: 19-25
   :dedent: 4

At the end we can print our results

.. literalinclude:: ../../examples/c/indexes.c
   :language: c
   :lines: 27-30
   :dedent: 4

And free all allocated memory. We don't need to free ``positions``, as it points
into memory allocated and controlled by the frame.

.. literalinclude:: ../../examples/c/indexes.c
   :language: c
   :lines: 32-34
   :dedent: 4

.. htmlhidden::
    :toggle: Click here to see the whole program
    :before-not-html: The whole code looks like this

    .. literalinclude:: ../../examples/c/indexes.c
       :language: c
       :lines: 4-

For more information about reading frame in a trajectory, see the following
functions:

- :cpp:func:`chfl_trajectory_nsteps` to know when to stop reading
- :cpp:func:`chfl_trajectory_read_step` to directlty read a given step.
- :cpp:func:`chfl_trajectory_set_cell` and
  :cpp:func:`chfl_trajectory_set_topology` to specify an unit cell or a topology
  for all frames in a trajectory.

Generating a structure
----------------------

Now that we know how to read frames from a trajectory, let's try to create a new
structure and write it to a file. As previsouly, we start by including the
chemfiles header:

.. literalinclude:: ../../examples/c/generate.c
   :language: c
   :lines: 4

We can then create a :cpp:type:`CHFL_FRAME` and the two :cpp:type:`CHFL_ATOM` we
will need.

.. literalinclude:: ../../examples/c/generate.c
   :language: c
   :lines: 8-10
   :dedent: 4

We can now add the atoms to the frame, with their respective positions. The
third argument to :cpp:func:`chfl_frame_add_atom` is set to ``NULL`` as we don't
need the velocities.

.. literalinclude:: ../../examples/c/generate.c
   :language: c
   :lines: 12-14
   :dedent: 4

We can then add bonds between the atoms to fully define the topology

.. literalinclude:: ../../examples/c/generate.c
   :language: c
   :lines: 16-17
   :dedent: 4

Finally, we can set the :cpp:type:`CHFL_CELL` associated with this frame. We
also free the cell memory right away, as it is no longer needed.

.. literalinclude:: ../../examples/c/generate.c
   :language: c
   :lines: 19-21
   :dedent: 4

Now that our frame is constructed, it is time to write it to a file. For that,
we open a trajectory in write (``'w'``) mode, and write to it.

.. literalinclude:: ../../examples/c/generate.c
   :language: c
   :lines: 23-25
   :dedent: 4

And free all remaining memory with the right function.

.. literalinclude:: ../../examples/c/generate.c
   :language: c
   :lines: 27-29
   :dedent: 4

.. htmlhidden::
    :toggle: Click here to see the whole program
    :before-not-html: Wrapping everything up, the whole code looks like this:

    .. literalinclude:: ../../examples/c/generate.c
       :language: c
       :lines: 4-

Using selections
----------------

Now that we know how to read and write frame from trajectories, how about we do
a bit a filtering? In this tutorial, we will read all the frames from a file,
and use :ref:`selections <selection-language>` to filter which atoms we will
write back to another file. This example will also show how chemfiles can be
used to convert from a file format to another one.

We start by opening the two trajectories we will need

.. literalinclude:: ../../examples/c/select.c
   :language: c
   :lines: 10-11
   :dedent: 4

We create a :cpp:type:`CHFL_FRAME` and a :cpp:type:`CHFL_SELECTION` object to
filter the atoms we want to keep.

.. literalinclude:: ../../examples/c/select.c
   :language: c
   :lines: 12-13
   :dedent: 4

Then we get the number of steps in the trajectory

.. literalinclude:: ../../examples/c/select.c
   :language: c
   :lines: 15-16
   :dedent: 4

And iterate over the frames in the trajectory

.. literalinclude:: ../../examples/c/select.c
   :language: c
   :lines: 18-19
   :dedent: 4

From here, we need to use the selection to get the atoms we want to remove. This
is a two steps process: first we evaluate the selection and get the number of
matches

.. literalinclude:: ../../examples/c/select.c
   :language: c
   :lines: 21-22
   :dedent: 8

Second we allocate some memory and get all the matches (represented as
:cpp:class:`chfl_match`):

.. literalinclude:: ../../examples/c/select.c
   :language: c
   :lines: 23-24
   :dedent: 8

We can get the index of atoms in a `to_remove` array

.. literalinclude:: ../../examples/c/select.c
   :language: c
   :lines: 26-29
   :dedent: 8

In order to remove the atoms from the frame, we need to sort ``to_remove`` in
descending order: removing the atom at index i will shift the index of all the
atoms after i. So we need start from the end and work toward the start of the
frame. The ``compare_matches`` function is used to do this sorting in reverse
order with the standard `qsort function`_.

.. _qsort function: https://en.cppreference.com/w/c/algorithm/qsort

.. literalinclude:: ../../examples/c/select.c
   :language: c
   :lines: 31-34
   :dedent: 8

Finally, we can write the cleaned frame to the output file, and free the memory
we allocated:

.. literalinclude:: ../../examples/c/select.c
   :language: c
   :lines: 36-39
   :dedent: 8

The ``compare_matches`` function we used to sort the matches is defined as
follow:

.. literalinclude:: ../../examples/c/select.c
   :language: c
   :lines: 45-52

.. htmlhidden::
    :toggle: Click here to see the whole program
    :before-not-html: The whole program look like this:

    .. literalinclude:: ../../examples/c/select.c
       :language: c
       :lines: 4-
.. _selection-language:

Selection language
==================

Introduction
^^^^^^^^^^^^

Chemfiles implements a rich selection language that allows finding the atoms
matching a set of constraints in a :cpp:class:`Frame <chemfiles::Frame>`.  For
example, ``atom: name == H and x > 15`` would select all atoms with a name equal
to "H" and x cartesian coordinate bigger than 15. Here, ``name == H`` and ``x >
15`` are individual constraints, and they are combined with ``and``, meaning
both of them must be true for an atom to match to full selection.

Chemfiles atomic selection language differs from other atomic selection
languages (such as the ones in `VMD`_, `MDAnalysis`_, and many others) by the
fact that it is possible to formulate constraints not only on single atoms, but
also on pairs, triplets, and quadruplets of atoms. For example, ``angles:
name(#2) O and mass(#3) < 1.5`` will select all sets of three bonded atoms
forming an angle such that the name of the second atom is O and the mass of the
third atom is less than 1.5. Here, the first atom is left unconstrained.  Where
evaluating simple selections yields a list of matching atomic indexes,
evaluating triplet selection will return a list of triplets of atomic indexes
(and correspondingly for pairs and quadruplets).

.. _VMD: https://www.ks.uiuc.edu/Research/vmd/current/ug/node89.html
.. _MDAnalysis: https://docs.mdanalysis.org/stable/documentation_pages/selections.html

The number of atoms to select together is indicated in chemfiles by a context,
separated from the main selection by a colon. Seven contexts are available:
``atoms`` is the default context, matching single atoms. ``two``, ``three``, and
``four`` match arbitrary pairs, triplets and quadruplets of atoms respectively.
``bonds``, ``angles``, and ``dihedrals`` match pairs, triplets and quadruplets
of atoms bonded together to form the corresponding connectivity element.

Expressing the constraints
--------------------------

Selections are built by assembling simple constraints with the Boolean operators
``and``, ``or`` and ``not``. They follow the usual interpretation of logic: ``A
and B`` will be true only if both A and B are true, ``A or B`` will be true if
one of A or B is, and ``not A`` will be true is A is false. You can use three
types of constraints in your selections: Boolean constraints, string constraints,
and numeric constraints.

String constraints check text values, such as the atomic ``name`` or atomic
``type``, comparing it to a fixed value: ``name(#1) == Fe``. When using a
selection with more than one atom, constraints are applied to different atoms
using ``#1``, ``#2``, ``#3`` or ``#4`` to refer to a specific atom: ``angles:
name(#3) Co`` will check the name of the third atom of the angles, and so on.

Numeric constraints check numeric values such as the position of the atoms
(``x``, ``y``, and ``z``), the atomic ``mass``, the ``index`` of an atom in the
frame, *etc.* Numeric values can be combined with the usual mathematical
operations: ``x^2 + y^2 + z^2 < 10^2`` will check for atoms inside the sphere
with a radius of 10 Å centered on the origin. The usual ``+``, ``-``, ``*`` and
``/`` operators are supported, as well as ``^`` for exponentiation and ``%`` for
modulo (remainder of Euclidean division). These operations follow the standard
priority rules: ``1 + 2 * 3`` is 7, not 9. Numeric values are then compared with
one of ``==`` (equal), ``!=`` (not equal), ``<`` (less than),
``<=`` (less or equal), ``>`` (more than), or ``>=`` (more or equal).

Finally, Boolean constraints directly evaluate to either ``true`` or ``false``
for a given set of atoms. For example, ``is_bonded(i, j)`` will check if atoms i
and j are bonded together. Here, ``i`` and ``j`` can either be one of the atoms
currently being matched (``#1 / #2 / #3 / #4``) or another selection (called
sub-selection). In the latter case, all the atoms in the sub-selection are
checked to see if any of them verify the selection. This makes ``is_bonded(#1,
name O)`` select all atoms bonded to an oxygen; and ``is_angle(type C, #1, name
O)`` select all atoms in the middle of a C-X-O angle.

Constraints on atomic or residue properties
-------------------------------------------

It is possible to use atomic :cpp:class:`properties <chemfiles::Property>`
(unfortunately not frame properties) in constraints, with the ``[<property>]``
syntax. ``<property>`` should be replaced by the property name: ``[is_hetatm]``,
possibly using quotes around the property name if it contains spaces: ``["my own
property"]``. Depending on the context, a Boolean, string or numeric property
will be searched. If an atomic property with a given name cannot be found,
residue properties are searched instead, however, atomic properties take
precedence. If none can be found, or if the property type does not match, a
default value will be used instead: ``false`` for Boolean properties, ``""``
(empty string) for string properties, and ``NaN`` for numeric properties.

Selection language reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is the list of currently implemented constraints. Additional ideas are welcome!

Boolean constraints
-------------------

Boolean constraints directly evaluate to ``true`` or ``false``, and can be
combined with other constraints with ``and``, ``or`` and ``not``.

.. chemfiles-selection:: ``all``

    Matches any atom and always evaluate to ``true``


.. chemfiles-selection:: ``none``

    Never matches an atom and always evaluate to ``false``


.. chemfiles-selection:: ``is_bonded(i, j)``

    Check if atoms i and j are bonded together.

    i and j can either be one of the atoms currently being matched (i.e. ``#1 /
    #2 / #3 / #4``); or a sub-selection (e.g. ``is_bonded(#1, name O)``). In the
    latter case, the constrain will be checked for all atoms returned by the
    sub-selection, and evaluate to ``true`` if any of them matches.

    If i and j refers to the same atom, this evaluates to ``false``.


.. chemfiles-selection:: ``is_angle(i, j, k)``

    Check if atoms i, j and k are bonded together to form an angle, *i.e.* that
    i is bonded to j and j is bonded to k.

    i, j, and k can either be one of the atoms currently being matched (i.e.
    ``#1 / #2 / #3 / #4``); or a sub-selection (with the same behavior as
    ``is_bonded``).

    If any of i, j or k refer to the same atom, ``is_angle`` evaluate to
    ``false``.


.. chemfiles-selection::  ``is_dihedral(i, j, k, m)``

    Check if atoms i, j, k and m are bonded together to form a dihedral angle,
    *i.e.* that i is bonded to  j, j is bonded to k, and k is bonded to m.

    i, j, k, and m can either be one of the atoms currently being matched (i.e.
    ``#1 / #2 / #3 / #4``); or a sub-selection (with the same behavior as
    ``is_bonded``).

    If any of i, j, k or m refer to the same atom, ``is_dihedral`` evaluate to
    ``false``.


.. chemfiles-selection:: ``is_improper(i, j, k, m)``

    Check if atoms i, j, k and m are bonded together to form a dihedral angle,
    *i.e.* that all of i, k, and m are bonded to j.

    If any of i, j, k or m refer to the same atom, ``is_improper`` evaluate to
    ``false``.

.. |multiple-atoms-#i| replace::

    If multiple atoms are being matched simultaneously, the one to check can be
    specified by setting ``#i`` to one of ``#1 / #2 / #3 / #4``. If ``#i`` is
    not specified, ``#1`` is checked by default.


.. chemfiles-selection:: ``[<property>]`` / ``[<property>](#i)``

    Check if atoms or residues they belong to have a Boolean property named
    ``'property'``, and that this property is set to ``true``.

    This will evaluate to ``false`` if the property does not exist.

    |multiple-atoms-#i|

String constraints
------------------

String constraints are compared to a literal value for either equality with
``constrain == value``, or inequality with ``constrain != value``. Values
without spaces, starting with an ASCII letter and using only ASCII letters and
numbers can be written directly (``type == Fe``, ``resname != ALA``), other
values must be surrounded by quotation marks: ``name == "17"``, ``resname != "a
long name"``.

For string selectors, ``constrain value`` is equivalent to ``constrain ==
value``, and ``constrain value1 value2 value3`` is equivalent to ``constrain ==
value1 or constrain == value2 or constrain == value3``. This simplifies common
selections, such as ``name H O C N``.

.. chemfiles-selection:: ``type`` / ``type(#i)``

    Evaluates to the atomic type of the atom being checked.

    |multiple-atoms-#i|

.. chemfiles-selection:: ``name`` / ``name(#i)``

    Evaluates to the atomic name of the atom being checked. Some formats store
    both an atomic name (e.g. H3) and an atom type (e.g. H), this is why you can
    use two different selectors depending on the actual data;

    |multiple-atoms-#i|

.. chemfiles-selection:: ``resname`` / ``resname(#i)``

    Evaluates to the name of the :cpp:class:`chemfiles::Residue` containing the
    atom being checked. If an atom is not in a residue, this evaluates to the
    empty string.

    |multiple-atoms-#i|

.. chemfiles-selection:: ``[<property>]`` / ``[<property>](#i)``

    If the atom being matched or the residue it belongs to (if any) have a
    string property named ``'property'``, this evaluates to the value of the
    property.  Else this evaluates to an empty string (``""``).

    |multiple-atoms-#i|

Numeric constraints
-------------------

Numeric constraints can be combined using all of the usual numeric operations:
``+``, ``-``, ``*``, ``/``; as well as ``^`` for exponentiation and ``%`` for
modulo (remainder of Euclidean division). Parenthesis can be used for grouping
operations.

Additionally, mathematical functions are available to transform a number to
another value. Currently supported functions are: ``deg2rad`` and ``rad2deg``
functions for transforming radians to degrees and respectively; ``sin``,
``cos``, ``tan`` for the trigonometric functions (taking radians as input);
``asin`` and ``acos`` inverse trigonometric functions (giving output in
radians); ``sqrt`` for square root; ``exp`` for base *e* exponential; ``log``
for the natural or base *e* logarithm; ``log2`` for the base 2 logarithm and
``log10`` for the base 10 logarithm. If you need another of such mathematical
function, `open an issue`_ in the chemfiles repository!

.. _open an issue: https://github.com/chemfiles/chemfiles/issues/new/choose

Numeric values are finally compared using one of ``==`` (equal), ``!=`` (not
equal), ``>`` (greater), ``>=`` (greater or equal), ``<`` (lesser), or ``<=``
(lesser or equal). They can be compared with one another (``(x^2 + z) > y``) or
with literal values (``z <= 1.24e-1``).

.. note::

    Numeric selection operate on double precision floating point number, and as
    such are subject to the same limitations. In particular, while ``1 + 2 ==
    3`` will match all atoms, since this relation is always true, ``0.1 + 0.2 ==
    0.3`` will match no atoms, since ``0.1 + 0.2 == 0.30000000000000004`` when
    using floating point arithmetic.

Here are all the available numeric values

.. chemfiles-selection:: ``index`` / ``index(#i)``

    Evaluates to the atomic index of the atom being matched.

    |multiple-atoms-#i|

.. chemfiles-selection:: ``mass`` / ``mass(#i)``

    Evaluates to the mass of the atom being matched. For atoms with a type from
    the periodic table, the element mass is used as default. The atom mass can
    also be set manually with :cpp:func:`chemfiles::Atom::set_mass`.

    |multiple-atoms-#i|

.. chemfiles-selection:: ``x`` / ``y`` / ``z`` / ``x(#i)`` / ``y(#i)`` / ``z(#i)``

    Evaluates to the cartesian component of the position of the atom being
    matched, in Ångströms.

    |multiple-atoms-#i|

.. chemfiles-selection:: ``vx`` / ``vy`` / ``vz`` / ``vx(#i)`` / ``vy(#i)`` / ``vz(#i)``

    Evaluates to the cartesian component of the velocity of the atom being
    matched. If velocities are not defined in the :cpp:class:`chemfiles::Frame`
    being matched, this evaluates to `NaN`_.

    |multiple-atoms-#i|

.. chemfiles-selection:: ``resid`` / ``resid(#i)``

    Evaluates to the id of the :cpp:class:`chemfiles::Residue` containing the
    atom being checked. If an atom is not in a residue, or if the residue do not
    have an id, this evaluates to -1.

    |multiple-atoms-#i|

.. chemfiles-selection:: ``[<property>]`` / ``[<property>](#i)``

    If the atom being matched or the residue it belongs to (if any) have a
    numeric property named ``'property'``, this evaluates to the value of the
    property. Else this evaluates to `NaN`_.

    |multiple-atoms-#i|

.. _NaN: https://en.wikipedia.org/wiki/Na


.. chemfiles-selection:: ``distance(i, j)``

    Evaluates to the distance in Ångströms between atoms i and j, accounting
    for periodic boundary conditions.

    i and j can either be one of the atoms currently being matched (i.e. ``#1 /
    #2 / #3 / #4``); or a sub-selection (``distance(#1, name O)``). In the case
    of a sub-selection, the distance will be evaluated for all atoms matching
    the sub-selection, and then propagated through mathematical operations. In
    the end, the atoms will be considered a match if **any** of the values is a
    match. So ``distance(#1, name O) + 3 < 5`` will return all atoms at less
    than 2Å from **any** atom with name O.


.. chemfiles-selection:: ``angle(i, j, k)``

    Evaluates to the angle between atoms i, j and k in radians, accounting for
    periodic boundary conditions. The atoms do not need to be bonded together.

    i, j, and k can either be one of the atoms currently being matched (i.e.
    ``#1 / #2 / #3 / #4``); or a sub-selections. In the latter case,
    mathematical operations are combined as for ``distance(i, j)``.


.. chemfiles-selection:: ``dihedral(i, j, k, m)``

    Evaluates to the dihedral angle between atoms i, j, k and m in radians,
    accounting for periodic boundary conditions. The atoms do not need to be
    bonded together.

    i, j, k and m can either be one of the atoms currently being matched (i.e.
    ``#1 / #2 / #3 / #4``); or a sub-selections. In the latter case,
    mathematical operations are combined as for ``distance(i, j)``.


.. chemfiles-selection:: ``out_of_plane(i, j, k, m)``

    Evaluates to  the distance in Ångströms between the plane formed by the three
    atoms i, k, and m; and the atom j, accounting for periodic boundary
    conditions.

    i, j, k and m can either be one of the atoms currently being matched (i.e.
    ``#1 / #2 / #3 / #4``); or a sub-selections. In the latter case,
    mathematical operations are combined as for ``distance(i, j)``.


.. note::

    ``angle`` and ``dihedral`` constraints are different from the ``is_angle``
    and ``is_dihedral``. The firsts returns a number that can then be used in
    mathematical expressions, and can be used even if the atoms are not bonded
    together in the corresponding topology elements.

    The seconds directly evaluate to ``true`` or ``false``, and only check if
    the atoms are bonded together in the corresponding topology elements.
Overview
========

The figure below represents how the basic classes of chemfiles are organized and
how they interact together. Chemfiles is organized around a handful of public
classes: :cpp:class:`Trajectory <chemfiles::Trajectory>`, :cpp:class:`Frame
<chemfiles::Frame>`, :cpp:class:`Topology <chemfiles::Topology>`,
:cpp:class:`Residue <chemfiles::Residue>`, :cpp:class:`Atom <chemfiles::Atom>`,
:cpp:class:`UnitCell <chemfiles::UnitCell>` and :cpp:class:`Selection
<chemfiles::Selection>`, all of which are presented here.

.. image:: ../static/img/classes.*
    :align: center
    :width: 500px

A :cpp:class:`trajectory <chemfiles::Trajectory>` is the main entry point of
chemfiles. It reads one or many :cpp:class:`frames <chemfiles::Frame>` from a
file on the disk using a specific format. The file type and the format are
automatically determined from the extension.

A :cpp:class:`frame <chemfiles::Frame>` holds data for one step of a simulation,
consisting in the positions for all the atoms; optionally the velocities for all
the atoms; the :cpp:class:`topology <chemfiles::Topology>` and the
:cpp:class:`unit cell <chemfiles::UnitCell>` of the system.

The :cpp:class:`topology <chemfiles::Topology>` describes the organization of
the particles in the system. It contains a list of :cpp:class:`atoms
<chemfiles::Atom>` in the system, and information about which atoms are bonded
together. A :cpp:class:`residue <chemfiles::Residue>` is a group of atoms bonded
together, which may or may not corresponds to molecules. When working with
bio-molecules and specifically proteins from the PDB data bank, the residues
should correspond to amino-acids in the protein.

The :cpp:class:`Atom <chemfiles::Atom>` class contains basic information about
the atoms in the system: the name (if it is available), mass, kind of atom and
so on. Atoms are not limited to plain chemical elements.

The :cpp:class:`UnitCell <chemfiles::UnitCell>` class describes the boundary
conditions of the system: where are the boundaries, and what is the periodicity
of theses boundaries. An unit cell can be of three types: *Infinite*,
*Orthorhombic* or *Triclinic*. Infinite cells do not have any boundaries.
Orthorhombic cells are defined by three orthogonal vectors, and triclinic cells
are defined by three vectors without any constrain.

The :cpp:class:`Property <chemfiles::Property>` class store additional data or
metadata associated with :cpp:class:`frames <chemfiles::Frame>`,
:cpp:class:`residues <chemfiles::Residue>` or :cpp:class:`atoms
<chemfiles::Atom>`. Properties can store string values, numeric values, Boolean
values or vector values.

Chemfiles also provides a :ref:`selection language <selection-language>`,
implemented in the :cpp:class:`Selection <chemfiles::Selection>` class. This
selection language allows the users to select a group of atoms from a
:cpp:class:`frame <chemfiles::Frame>` using a selection string such as ``"(x <
45 and name O) or name C"``.
Installation
============

.. note::

    This page describe the different installation methods for the C and C++
    interfaces to chemfiles. If you want to use chemfiles from a different
    language, please see the corresponding pages:

    - Python: https://chemfiles.org/chemfiles.py/latest/
    - Julia: https://chemfiles.org/Chemfiles.jl/latest/
    - Rust: https://chemfiles.org/chemfiles.rs/latest/
    - Fortran: https://chemfiles.org/chemfiles.f03/latest/

Pre-compiled chemfiles
^^^^^^^^^^^^^^^^^^^^^^

We provide pre-compiled version of chemfiles (both a static library build and a
shared library build) on github. Please look for the version you need on the
`chemfiles/chemfiles-prebuilt
<https://github.com/chemfiles/chemfiles-prebuilt/releases>`_ repository.

--------------------------------------------------------------------------------

We also provide three `conda`_ packages in the `conda-forge` community channel for
Linux and Os X.

- `chemfiles <https://github.com/conda-forge/chemfiles-feedstock>`_ package
  contains both the C++ and Python libraries. This is the one you want most of
  the time.
- `chemfiles-lib <https://github.com/conda-forge/chemfiles-lib-feedstock>`_
  package contains the C++ library alone;
- `chemfiles-python <https://github.com/conda-forge/chemfiles-feedstock>`_
  package contains the Python binding to Chemfiles alone.

.. code-block:: bash

    # Get all the packages
    conda install -c conda-forge chemfiles

    # Get the C++ and C library
    conda install -c conda-forge chemfiles-lib

.. _conda: https://conda.io/
.. _OpenSUSE build service: https://software.opensuse.org/download.html?project=home%3ALuthaf&package=chemfiles

Building from sources
^^^^^^^^^^^^^^^^^^^^^

Provided you have all the required dependencies, you can also build chemfiles
from sources.

Core library dependencies
-------------------------

In order to build the core library, you will need a C++11 capable compiler.
Chemfiles is automatically tested against GCC (>= 4.8) on Linux, OSX and Windows
(mingw-w64); clang (>= 3.3) on Linux and OS X and MSCV 14 on Windows. It was
also compiled successfully with Intel C++ compilers. Please report any successful
compilation with other compilers! Chemfiles also uses the `CMake`_ build system,
which is available in all the package managers.

On UNIX-like systems (Linux, OS X, ...)
"""""""""""""""""""""""""""""""""""""""

All these dependencies can be installed in one command:

.. code-block:: bash

    # On apt-get based distributions
    apt-get update
    apt-get install build-essential cmake

    # On yum based distribution (CentOS/RHEL)
    yum install gcc gcc-c++ make cmake

    # On dnf based distribution (Fedora)
    dnf install @development-tools cmake

    # On OS X with Homebrew
    brew install cmake

.. _CMake: https://cmake.org/

On Windows
""""""""""

You can use either MSVC 2015 compiler, or gcc as provided by `mingw`_.
`MSYS2`_ offer a package manager to install all the needed libraries. After the
initial installation steps, you can run the following to install a recent C++
compiler:

.. code-block:: bash

    # On 64-bits windows
    pacman -S mingw64/mingw-w64-x86_64-gcc

    # On 32-bits windows
    pacman -S mingw32/mingw-w64-i686-gcc

You will also need to install cmake, which can be found `here <https://cmake.org/download/>`_.

.. _mingw: http://www.mingw.org/
.. _MSYS2: https://www.msys2.org/

Build steps
-----------

You can get the source code from either git, or from the `release`_ page of
Github. In the later case, just unpack the archive wherever you want the source
code to live. To get the latest development version, use git:

.. code-block:: bash

    cd where/you/whant/chemfiles/to/live
    git clone https://github.com/chemfiles/chemfiles
    cd chemfiles

.. _release: https://github.com/Luthaf/chemfiles/releases

The following command build and install chemfiles

.. code-block:: bash

    cd chemfiles
    mkdir build
    cd build
    cmake .. # various options are allowed here
    cmake --build .
    # if you whant to run the tests before installing:
    ctest
    cmake --build . --target install

The :command:`cmake` step can be further configured by using the curse-based GUI
(:command:`ccmake .`) or providing some command-line arguments. Here are the
most important options:

+---------------------------------------+---------------------+------------------------------+
| Option                                | Default value       | Effect/Informations          |
+=======================================+=====================+==============================+
| ``-DCMAKE_INSTALL_PREFIX=prefix``     | :file:`/usr/local`  | Set the installation prefix  |
|                                       |                     | to ``prefix``                |
+---------------------------------------+---------------------+------------------------------+
| ``-DCMAKE_BUILD_TYPE=type``           | ``release``         | Set to ``debug`` for debug   |
|                                       |                     | information                  |
+---------------------------------------+---------------------+------------------------------+
| ``-DBUILD_SHARED_LIBS=ON|OFF``        | ``OFF``             | Build shared library instead |
|                                       |                     | of static one.               |
+---------------------------------------+---------------------+------------------------------+
| ``-DCHFL_BUILD_DOCUMENTATION=ON|OFF`` | ``OFF``             | Build the documentation.     |
|                                       |                     | This needs `sphinx`_ and     |
|                                       |                     | `doxygen`_ to be installed   |
+---------------------------------------+---------------------+------------------------------+
| ``-DCHFL_BUILD_TESTS=ON|OFF``         | ``OFF``             | Build the test suite.        |
+---------------------------------------+---------------------+------------------------------+
| ``-DCHFL_SYSTEM_NETCDF=ON|OFF``       | ``OFF``             | Use the system-provided      |
|                                       |                     | netcdf library.              |
+---------------------------------------+---------------------+------------------------------+
| ``-DCHFL_SYSTEM_LZMA=ON|OFF``         | ``OFF``             | Use the system-provided      |
|                                       |                     | lzma library.                |
+---------------------------------------+---------------------+------------------------------+
| ``-DCHFL_SYSTEM_ZLIB=ON|OFF``         | ``OFF``             | Use the system-provided zlib |
+---------------------------------------+---------------------+------------------------------+

For instance, to install chemfiles to :file:`$HOME/local`, you should use:

.. code-block:: bash

    cmake -DCMAKE_INSTALL_PREFIX=$HOME/local ..

.. _doxygen: http://doxygen.org/
.. _sphinx: http://sphinx-doc.org/


Using chemfiles in your project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are multiple ways to use chemfiles in your own code:

- adding the include path and library manually (in a Makefile, or a Visual Studio project);
- using the CMake configuration file;
- including chemfiles inside a CMake based-project.

Manually setting include and library path
-----------------------------------------

After installing chemfiles on your computer, you can start using it with your
own C or C++ program by passing the corresponding include path and library path
to your compiler. For example, on \*nix (GNU/Linux or OS X) you can compile any
code depending on chemfiles with the following command

.. code-block:: bash

    # change <PREFIX> to the location where you installed chemfiles
    # (default is /usr/local)
    g++ my-code.cpp -o my-code -I<PREFIX>/include -lchemfiles -L<PREFIX>/lib

Here, ``-I <PREFIX>/include`` tells the compiler where to look for chemfiles
headers, ``-lchemfiles`` tells it to link the chemfiles library in the final
executable, and ``-L <PREFIX>/lib`` tells the compiler where to look for the
chemfiles library.

The same strategy should be possible to use with Visual Studio on Windows, or
any other IDE. Refer to your IDE documentation about how to add external
libraries.

Using cmake and ``find_package``
--------------------------------

If your project is already using CMake, and you installed chemfiles on your
computer, you can use the standard ``find_package`` to find the code and
automatically set the right include and library path.

.. code-block:: cmake

    add_executable(my-code my-code.cpp)

    find_package(chemfiles 0.8)
    # chemfiles_FOUND will be TRUE if the code was found

    target_link_library(my-code chemfiles)

Including chemfiles as a CMake subproject
-----------------------------------------

If your project is already using CMake, but you don't want to require your users
to install chemfiles separately, you can use cmake support for external
projects or subdirectories to include chemfiles sources directly inside your own
project. All CMake variables controlling chemfiles behavior are prefixed with
``CHFL_`` to prevent variable pollution.
Properties
==========

When reading files, chemfiles read multiple kind of data: positions, velocities,
atomic information (name, type, mass, charge, ...), bonds and other
connectivity elements. This model allow to read the most commonly available data
in various formats. But sometimes, a format defines additional information.
Instead of adding a new field/function for every kind of data there can be in a
file, chemfiles defines a generic interface to read and store this additional
data. These additional data are stored inside properties.

A property has a name and a value. The value can either be a real number, a
string, a Boolean value (true/false) or a 3 dimensional vector. A property is
either stored inside an atom, and associated with this atom (for example the
total atomic force), or stored in and associated with a frame. The later case is
used for general properties, such as the temperature of the system, or the
author of the file.

This section documents which format set and use properties.

Atomic properties
-----------------

.. the CSV files at properties/*.csv are generated from the corresponding TOML
.. files. Please edit the TOML files when adding new properties

.. csv-table::
   :file: properties/atom.csv
   :widths: 10, 10, 10, 70
   :header-rows: 1

Residue properties
------------------

.. csv-table::
   :file: properties/residue.csv
   :widths: 10, 10, 10, 70
   :header-rows: 1

Frame properties
-----------------

.. csv-table::
   :file: properties/frame.csv
   :widths: 10, 10, 10, 70
   :header-rows: 1

Additionally, the SDF format reads any property formatted as ``> <...>``, using
the value inside the angle brackets (``...`` here) as the property name.
.. _cxx-tutorials:

C++ Tutorials
=============

This section present some hand-on tutorials to the chemfiles C++ API. All the
code here is under the `CC-0 Universal Licence`_ which means that you are free
to do whatever you want with it (*i.e.* it is Public Domain code)

.. _CC-0 Universal Licence: https://creativecommons.org/publicdomain/zero/1.0/

Read a single frame
-------------------

In this tutorials we will read a frame from a trajectory, and print the indexes
of all the atom in the half-space ``x < 5``.

We start by including the headers we need: chemfiles header, iostream for
printing and vector for storing the indexes. All the public classes of chemfiles
are accessible through the ``chemfiles.hpp`` header, which is the only header
needed.

.. literalinclude:: ../../examples/cpp/indexes.cpp
   :language: cpp
   :lines: 4-6

Then we open a Trajectory and read the first frame:

.. literalinclude:: ../../examples/cpp/indexes.cpp
   :language: cpp
   :lines: 9-10
   :dedent: 4

We can now create a vector to store the indices of the atoms with ``x < 5``, and
get the positions of the atoms in the frame with the
:cpp:func:`chemfiles::Frame::positions` function

.. literalinclude:: ../../examples/cpp/indexes.cpp
   :language: cpp
   :lines: 12-13
   :dedent: 4

Iterating through the atoms in the frame, we get the ones matching our
condition. :cpp:func:`Frame::size` gives the number of atoms in the frame, which
is also the size of the ``positions`` array. This array contains
:cpp:class:`chemfiles::Vector3D` representing the positions of the atoms in
Angstroms.

.. literalinclude:: ../../examples/cpp/indexes.cpp
   :language: cpp
   :lines: 15-19
   :dedent: 4

And finally we can print our results

.. literalinclude:: ../../examples/cpp/indexes.cpp
   :language: cpp
   :lines: 21-24
   :dedent: 4

.. htmlhidden::
    :toggle: Click here to see the whole program
    :before-not-html: The whole code looks like this

    .. literalinclude:: ../../examples/cpp/indexes.cpp
       :language: cpp
       :lines: 4-

For more information about reading frame in a trajectory, see the following
functions:

- :cpp:func:`Trajectory::nsteps` and :cpp:func:`Trajectory::done` to know when
  to stop reading
- :cpp:func:`Trajectory::read_step` to directlty read a given step.
- :cpp:func:`Trajectory::set_cell` and :cpp:func:`Trajectory::set_topology` to
  specify an unit cell or a topology for all frames in a trajectory.

Generating a structure
----------------------

Now that we know how to read frames from a trajectory, let's try to create a new
structure and write it to a file. As previsouly, we start by including the
chemfiles header:

.. literalinclude:: ../../examples/cpp/generate.cpp
   :language: cpp
   :lines: 4

We start with a :cpp:type:`chemfiles::Frame` that we construct by giving it the
:cpp:type:`chemfiles::UnitCell` that defines the periodic boundary conditions.

.. literalinclude:: ../../examples/cpp/generate.cpp
   :language: cpp
   :lines: 7
   :dedent: 4

We can then add three :cpp:class:`chemfiles::Atom` to this frame with their
positions; and set the bonds between them to define a water molecule.

.. literalinclude:: ../../examples/cpp/generate.cpp
   :language: cpp
   :lines: 9-14
   :dedent: 4

Now that our frame is constructed, it is time to write it to a file. For that,
we open a trajectory in write (``'w'``) mode, and write to it.

.. literalinclude:: ../../examples/cpp/generate.cpp
   :language: cpp
   :lines: 16-17
   :dedent: 4

.. htmlhidden::
    :toggle: Click here to see the whole program
    :before-not-html: Wrapping everything up, the whole code looks like this:

    .. literalinclude:: ../../examples/cpp/generate.cpp
       :language: cpp
       :lines: 4-

Using selections
----------------

Now that we know how to read and write frame from trajectories, how about we do
a bit a filtering? In this tutorial, we will read all the frames from a file,
and use :ref:`selections <selection-language>` to filter which atoms we will
write back to another file. This example will also show how chemfiles can be
used to convert from a file format to another one.

We start by opening the two trajectories we will need

.. literalinclude:: ../../examples/cpp/select.cpp
   :language: cpp
   :lines: 8-9
   :dedent: 4

And we create a :cpp:class:`chemfiles::Selection` object to filter the atoms we
want to keep.

.. literalinclude:: ../../examples/cpp/select.cpp
   :language: cpp
   :lines: 11
   :dedent: 4

Then we can iterate over all the frames in the trajectory

.. literalinclude:: ../../examples/cpp/select.cpp
   :language: cpp
   :lines: 13-14
   :dedent: 4

And use the selection to get the list of atoms to remove. The result of
:cpp:func:`Selection::list` is a ``std::vector<size_t>`` containing the atoms
matching the selection.

.. literalinclude:: ../../examples/cpp/select.cpp
   :language: cpp
   :lines: 15
   :dedent: 8

In order to remove the atoms from the frame, we need to sort the ``to_remove``
vector in descending order: removing the atom at index i will shift the index of
all the atoms after i. So we start from the end and work toward the start of the
frame.

.. literalinclude:: ../../examples/cpp/select.cpp
   :language: cpp
   :lines: 16-19
   :dedent: 8

Finally, we can write the cleaned frame to the output file, and start the next
iteration.

.. literalinclude:: ../../examples/cpp/select.cpp
   :language: cpp
   :lines: 20
   :dedent: 8

.. htmlhidden::
    :toggle: Click here to see the whole program
    :before-not-html: The whole program look like this:

    .. literalinclude:: ../../examples/cpp/select.cpp
       :language: cpp
       :lines: 4-
Chemfiles, a modern library for chemistry file reading and writing
==================================================================

Chemfiles is a multi-language library written in modern C++ for reading and
writing from and to molecular trajectory files. These files are created by your
favorite theoretical chemistry program, and contains information about atomic
or residues names and positions. Some format also have additional information,
such as velocities, forces, energy, …

The main targeted audience of chemfiles are chemistry researchers working on
their own code to do some kind of awesome science, and who do not want to bother
about handling all the files format that may exist in the world.

Running simulation (either Quantum Dynamic, Monte Carlo, Molecular Dynamic, or
any other method) often produce enormous amounts of data, which had to be
post-processed in order to extract information. This post-processing step
involve reading and parsing the data, and computing physical values with the
help of statistical thermodynamic. Chemfiles tries to help you on the first
point, by providing the same interface to all the trajectory formats. If you
ever need to change your output format, your analysis tools will still work the
same way. Chemfiles is efficient because it allow you to write and debug your
code only once, and then to re-use it as needed.

Even if chemfiles is written in C++, it can be used from the most popular
scientific programming languages: C, Fortran, Python, … You can just pick up
your favorite language to use it. This part of the documentation presents the
data model used by chemfiles to store information about the trajectories, and
how to access that data in C and C++. The documentation for the other languages
interfaces to chemfiles are accessible at the following places:

* `Python interface`_, usable with Python 2 and 3;
* `Fortran interface`_, for Fortran 2003;
* `Julia interface`_, for `Julia`_ 0.6+;
* `Rust interface`_, for the `Rust`_ language;

.. note::

    Chemfiles follow `semantic versionning <https://semver.org/>`_. This means
    that all 0.x.y versions are compatible for all y, but 0.x and 0.(x+1) are
    not compatible. This also means that when chemfiles reaches 1.0.0, all code
    using 1.0.0 will be compatible with 1.x.y for any x and y.

.. _Python interface: https://chemfiles.github.io/chemfiles.py/latest/
.. _Fortran interface: https://chemfiles.github.io/chemfiles.f03/latest/
.. _Julia interface: https://chemfiles.github.io/Chemfiles.jl/latest/
.. _Rust interface: https://chemfiles.github.io/chemfiles.rs/latest/
.. _Julia: https://julialang.org/
.. _Rust: https://www.rust-lang.org/

User manual
-----------

This section of the documentation will cover how to install chemfiles, give an
high-level overview of the functionalities and the supported formats. In order
to use chemfiles in your code, you should refer to the :ref:`classes-reference`.
If you need help installing or using chemfiles, do not hesitate to ask questions
in our `Gitter chat room <https://gitter.im/chemfiles/chemfiles>`_.

.. toctree::
    :maxdepth: 2

    installation
    overview
    cxx-tutorials
    c-tutorials
    formats
    selections
    properties

.. _classes-reference:

API reference
-------------

All the public classes (for the C++ interface) and functions (for the C
interface) are extensively documented here.

.. toctree::
   :maxdepth: 2

   classes/index
   capi/index

Developer documentation
-----------------------

.. toctree::
   :maxdepth: 2

   devdoc/contributing
   devdoc/internals
Supported formats
=================

List of supported formats
-------------------------

This page list the supported formats in chemfiles, with a link to the format
specification. For each format, it is also indicated which kind of data can be
read or written with |yes| for 'yes' and |no| for 'no'. Even if chemfiles can
read some type of data from a file, this does not mean that all files of this
format will contain this type of data.

The *Extension* column gives the extension used when trying to detect the
format. If your file does not have this extension, you can also use the string
in the *Format* column as a parameter to the :cpp:class:`chemfiles::Trajectory`
constructor to manually specify which format to use.

.. role:: red

.. |yes| replace:: `✓`
.. |no| replace:: :red:`✗`

.. csv-table::
   :file: formats-overview.csv
   :header-rows: 1
   :class: formats-table
   :widths: 20, 11, 6, 6, 6, 7, 7, 7, 9, 9, 9


- **LAMMPS** format corresponds to trajectory files written by the LAMMPS `dump
  <https://lammps.sandia.gov/doc/dump.html>`_ command.
- **LAMMPS Data** format corresponds to LAMMPS data files, as read by the LAMMPS
  `read_data <https://lammps.sandia.gov/doc/read_data.html>`_ command.

.. note:: in-memory IO

    Additionally, some formats support reading and writing directly to memory,
    without going through a file. At this time, all text based files (excluding
    those backed by the Molfiles plugin) support both reading and writing directly
    to memory. The MMTF format supports reading from a memory buffer, but does not
    support writing. It is also possible to read a compressed GZ or XZ file directly
    to memory buffer, but writing compressed files is not supported.

Asking for a new format
-----------------------

If you want to use chemfiles with a format which is not yet implemented, you may
easily add it by yourself if you know some C++. See the ``src/formats/XYZ.cpp``
file for example. The list of planned formats can be found `here
<gh-new-format_>`_. If you can not find your favorite format in this list, you
can `open an issue <gh-new-issue_>`_ with a description of the format, or even
better a link to the format specification.

.. _gh-new-format: https://github.com/chemfiles/chemfiles/labels/A-formats
.. _gh-new-issue: https://github.com/chemfiles/chemfiles/issues/new
.. _capi-cell:

``CHFL_CELL``
-------------

.. doxygentypedef:: CHFL_CELL

.. only:: html

    Here is the full list of functions acting on :cpp:type:`CHFL_CELL`:

    - :cpp:func:`chfl_cell`
    - :cpp:func:`chfl_cell_from_matrix`
    - :cpp:func:`chfl_cell_from_frame`
    - :cpp:func:`chfl_cell_copy`
    - :cpp:func:`chfl_cell_volume`
    - :cpp:func:`chfl_cell_lengths`
    - :cpp:func:`chfl_cell_set_lengths`
    - :cpp:func:`chfl_cell_angles`
    - :cpp:func:`chfl_cell_set_angles`
    - :cpp:func:`chfl_cell_matrix`
    - :cpp:func:`chfl_cell_shape`
    - :cpp:func:`chfl_cell_set_shape`
    - :cpp:func:`chfl_cell_wrap`

    --------------------------------------------------------------------

.. doxygenfunction:: chfl_cell

.. doxygenfunction:: chfl_cell_from_matrix

.. doxygenfunction:: chfl_cell_from_frame

.. doxygenfunction:: chfl_cell_copy

.. doxygenfunction:: chfl_cell_volume

.. doxygenfunction:: chfl_cell_lengths

.. doxygenfunction:: chfl_cell_set_lengths

.. doxygenfunction:: chfl_cell_angles

.. doxygenfunction:: chfl_cell_set_angles

.. doxygenfunction:: chfl_cell_matrix

.. doxygenenum:: chfl_cellshape

.. doxygenfunction:: chfl_cell_shape

.. doxygenfunction:: chfl_cell_set_shape

.. doxygenfunction:: chfl_cell_wrap
.. _capi-topology:

``CHFL_TOPOLOGY``
-----------------

.. doxygentypedef:: CHFL_TOPOLOGY

.. only:: html

    Here is the full list of functions acting on :cpp:type:`CHFL_TOPOLOGY`:

    - :cpp:func:`chfl_topology`
    - :cpp:func:`chfl_topology_from_frame`
    - :cpp:func:`chfl_topology_copy`
    - :cpp:func:`chfl_topology_atoms_count`
    - :cpp:func:`chfl_topology_resize`
    - :cpp:func:`chfl_topology_add_atom`
    - :cpp:func:`chfl_topology_remove`
    - :cpp:func:`chfl_topology_add_bond`
    - :cpp:func:`chfl_topology_remove_bond`
    - :cpp:func:`chfl_topology_clear_bonds`
    - :cpp:func:`chfl_topology_bonds_count`
    - :cpp:func:`chfl_topology_angles_count`
    - :cpp:func:`chfl_topology_dihedrals_count`
    - :cpp:func:`chfl_topology_impropers_count`
    - :cpp:func:`chfl_topology_bonds`
    - :cpp:func:`chfl_topology_angles`
    - :cpp:func:`chfl_topology_dihedrals`
    - :cpp:func:`chfl_topology_impropers`
    - :cpp:func:`chfl_topology_add_residue`
    - :cpp:func:`chfl_topology_residues_count`
    - :cpp:func:`chfl_topology_residues_linked`
    - :cpp:func:`chfl_topology_bond_with_order`
    - :cpp:func:`chfl_topology_bond_orders`
    - :cpp:func:`chfl_topology_bond_order`

    --------------------------------------------------------------------

.. doxygenfunction:: chfl_topology

.. doxygenfunction:: chfl_topology_from_frame

.. doxygenfunction:: chfl_topology_copy

.. doxygenfunction:: chfl_topology_atoms_count

.. doxygenfunction:: chfl_topology_resize

.. doxygenfunction:: chfl_topology_add_atom

.. doxygenfunction:: chfl_topology_remove

.. doxygenfunction:: chfl_topology_add_bond

.. doxygenfunction:: chfl_topology_remove_bond

.. doxygenfunction:: chfl_topology_clear_bonds

.. doxygenfunction:: chfl_topology_bonds_count

.. doxygenfunction:: chfl_topology_angles_count

.. doxygenfunction:: chfl_topology_dihedrals_count

.. doxygenfunction:: chfl_topology_impropers_count

.. doxygenfunction:: chfl_topology_bonds

.. doxygenfunction:: chfl_topology_angles

.. doxygenfunction:: chfl_topology_dihedrals

.. doxygenfunction:: chfl_topology_impropers

.. doxygenfunction:: chfl_topology_add_residue

.. doxygenfunction:: chfl_topology_residues_count

.. doxygenfunction:: chfl_topology_residues_linked

.. doxygenfunction:: chfl_topology_bond_with_order

.. doxygenfunction:: chfl_topology_bond_orders

.. doxygenfunction:: chfl_topology_bond_order
.. _capi-property:

``CHFL_PROPERTY``
-----------------

.. doxygentypedef:: CHFL_PROPERTY

.. only:: html

    Here is the full list of functions acting on :cpp:type:`CHFL_PROPERTY`:

    - :cpp:func:`chfl_property_bool`
    - :cpp:func:`chfl_property_double`
    - :cpp:func:`chfl_property_string`
    - :cpp:func:`chfl_property_vector3d`
    - :cpp:func:`chfl_property_get_kind`
    - :cpp:func:`chfl_property_get_bool`
    - :cpp:func:`chfl_property_get_double`
    - :cpp:func:`chfl_property_get_string`
    - :cpp:func:`chfl_property_get_vector3d`

    --------------------------------------------------------------------

.. doxygenfunction:: chfl_property_bool

.. doxygenfunction:: chfl_property_double

.. doxygenfunction:: chfl_property_string

.. doxygenfunction:: chfl_property_vector3d

.. doxygenfunction:: chfl_property_get_kind

.. doxygenfunction:: chfl_property_get_bool

.. doxygenfunction:: chfl_property_get_double

.. doxygenfunction:: chfl_property_get_string

.. doxygenfunction:: chfl_property_get_vector3d
.. _c-api:

C interface reference
=====================

The C interface is define in the ``chemfiles.h``, in which all the functions and
enums have a ``chfl_`` prefix for namespacing. Main chemfiles types are the
following opaque pointer types:

* :ref:`CHFL_TRAJECTORY <capi-trajectory>` maps to the :ref:`Trajectory <class-Trajectory>` class;
* :ref:`CHFL_FRAME <capi-frame>` maps to the :ref:`Frame <class-Frame>` class;
* :ref:`CHFL_ATOM <capi-atom>` maps to the :ref:`Atom <class-Atom>` class;
* :ref:`CHFL_CELL <capi-cell>` maps to the :ref:`UnitCell <class-UnitCell>` class;
* :ref:`CHFL_TOPOLOGY <capi-topology>` maps to the :ref:`Topology <class-Topology>` class.
* :ref:`CHFL_RESIDUE <capi-residue>` maps to the :ref:`Residue <class-Residue>` class.
* :ref:`CHFL_SELECTION <capi-selection>` maps to the :ref:`Selection <class-Selection>` class.
* :ref:`CHFL_PROPERTY <capi-property>` maps to the :ref:`Property <class-Property>` class.

The user is reponsible for memory management when using these types.
Constructors functions (functions returning pointers to types defined above)
return freshly allocated memory, and calling the :cpp:func:`chfl_free` function
return the corresponding memory to the operating system.

.. doxygenfunction:: chfl_free

.. doxygentypedef:: chfl_vector3d

``chemfiles.h`` also exports the :ref:`same macro <exported-macro>` as in the
C++ interface, and the ``chfl_version`` function allow to access the version of
the Chemfiles library.

.. doxygenfunction:: chfl_version

.. toctree::
    :maxdepth: 2

    chfl_trajectory
    chfl_frame
    chfl_topology
    chfl_residue
    chfl_atom
    chfl_cell
    chfl_selection
    chfl_property
    misc
.. _capi-residue:

``CHFL_RESIDUE``
----------------

.. doxygentypedef:: CHFL_RESIDUE

.. only:: html

    Here is the full list of functions acting on :cpp:type:`CHFL_RESIDUE`:

    - :cpp:func:`chfl_residue`
    - :cpp:func:`chfl_residue_with_id`
    - :cpp:func:`chfl_residue_for_atom`
    - :cpp:func:`chfl_residue_from_topology`
    - :cpp:func:`chfl_residue_copy`
    - :cpp:func:`chfl_residue_name`
    - :cpp:func:`chfl_residue_id`
    - :cpp:func:`chfl_residue_add_atom`
    - :cpp:func:`chfl_residue_atoms_count`
    - :cpp:func:`chfl_residue_atoms`
    - :cpp:func:`chfl_residue_contains`
    - :cpp:func:`chfl_residue_set_property`
    - :cpp:func:`chfl_residue_get_property`
    - :cpp:func:`chfl_residue_properties_count`
    - :cpp:func:`chfl_residue_list_properties`

    --------------------------------------------------------------------

.. doxygenfunction:: chfl_residue

.. doxygenfunction:: chfl_residue_with_id

.. doxygenfunction:: chfl_residue_for_atom

.. doxygenfunction:: chfl_residue_from_topology

.. doxygenfunction:: chfl_residue_copy

.. doxygenfunction:: chfl_residue_name

.. doxygenfunction:: chfl_residue_id

.. doxygenfunction:: chfl_residue_add_atom

.. doxygenfunction:: chfl_residue_atoms_count

.. doxygenfunction:: chfl_residue_atoms

.. doxygenfunction:: chfl_residue_contains

.. doxygenfunction:: chfl_residue_set_property

.. doxygenfunction:: chfl_residue_get_property

.. doxygenfunction:: chfl_residue_properties_count

.. doxygenfunction:: chfl_residue_list_properties
.. _capi-trajectory:

``CHFL_TRAJECTORY``
-------------------

.. doxygentypedef:: CHFL_TRAJECTORY

.. only:: html

    Here is the full list of functions acting on :cpp:type:`CHFL_TRAJECTORY`:

    - :cpp:func:`chfl_trajectory_open`
    - :cpp:func:`chfl_trajectory_with_format`
    - :cpp:func:`chfl_trajectory_memory_reader`
    - :cpp:func:`chfl_trajectory_memory_writer`
    - :cpp:func:`chfl_trajectory_path`
    - :cpp:func:`chfl_trajectory_read`
    - :cpp:func:`chfl_trajectory_read_step`
    - :cpp:func:`chfl_trajectory_write`
    - :cpp:func:`chfl_trajectory_set_cell`
    - :cpp:func:`chfl_trajectory_set_topology`
    - :cpp:func:`chfl_trajectory_topology_file`
    - :cpp:func:`chfl_trajectory_nsteps`
    - :cpp:func:`chfl_trajectory_memory_buffer`
    - :cpp:func:`chfl_trajectory_close`

    --------------------------------------------------------------------

.. doxygenfunction:: chfl_trajectory_open

.. doxygenfunction:: chfl_trajectory_with_format

.. doxygenfunction:: chfl_trajectory_memory_reader

.. doxygenfunction:: chfl_trajectory_memory_writer

.. doxygenfunction:: chfl_trajectory_path

.. doxygenfunction:: chfl_trajectory_read

.. doxygenfunction:: chfl_trajectory_read_step

.. doxygenfunction:: chfl_trajectory_write

.. doxygenfunction:: chfl_trajectory_set_cell

.. doxygenfunction:: chfl_trajectory_set_topology

.. doxygenfunction:: chfl_trajectory_topology_file

.. doxygenfunction:: chfl_trajectory_nsteps

.. doxygenfunction:: chfl_trajectory_memory_buffer

.. doxygenfunction:: chfl_trajectory_close
.. _capi-selection:

``CHFL_SELECTION``
------------------

.. doxygentypedef:: CHFL_SELECTION

.. only:: html

    Here is the full list of functions acting on :cpp:type:`CHFL_SELECTION`:

    - :cpp:func:`chfl_selection`
    - :cpp:func:`chfl_selection_copy`
    - :cpp:func:`chfl_selection_size`
    - :cpp:func:`chfl_selection_string`
    - :cpp:func:`chfl_selection_evaluate`
    - :cpp:func:`chfl_selection_matches`

    --------------------------------------------------------------------

.. doxygenfunction:: chfl_selection

.. doxygenfunction:: chfl_selection_copy

.. doxygenfunction:: chfl_selection_size

.. doxygenfunction:: chfl_selection_string

.. doxygenfunction:: chfl_selection_evaluate

.. doxygenfunction:: chfl_selection_matches


.. doxygenstruct:: chfl_match
    :members:

.. doxygendefine:: CHFL_MAX_SELECTION_SIZE
.. _capi-atom:

``CHFL_ATOM``
-------------

.. doxygentypedef:: CHFL_ATOM

.. only:: html

    Here is the full list of functions acting on :cpp:type:`CHFL_ATOM`:

    - :cpp:func:`chfl_atom`
    - :cpp:func:`chfl_atom_from_frame`
    - :cpp:func:`chfl_atom_from_topology`
    - :cpp:func:`chfl_atom_copy`
    - :cpp:func:`chfl_atom_mass`
    - :cpp:func:`chfl_atom_set_mass`
    - :cpp:func:`chfl_atom_charge`
    - :cpp:func:`chfl_atom_set_charge`
    - :cpp:func:`chfl_atom_name`
    - :cpp:func:`chfl_atom_set_name`
    - :cpp:func:`chfl_atom_type`
    - :cpp:func:`chfl_atom_set_type`
    - :cpp:func:`chfl_atom_full_name`
    - :cpp:func:`chfl_atom_vdw_radius`
    - :cpp:func:`chfl_atom_covalent_radius`
    - :cpp:func:`chfl_atom_atomic_number`
    - :cpp:func:`chfl_atom_set_property`
    - :cpp:func:`chfl_atom_get_property`
    - :cpp:func:`chfl_atom_properties_count`
    - :cpp:func:`chfl_atom_list_properties`

    --------------------------------------------------------------------


.. doxygenfunction:: chfl_atom

.. doxygenfunction:: chfl_atom_from_frame

.. doxygenfunction:: chfl_atom_from_topology

.. doxygenfunction:: chfl_atom_copy

.. doxygenfunction:: chfl_atom_mass

.. doxygenfunction:: chfl_atom_set_mass

.. doxygenfunction:: chfl_atom_charge

.. doxygenfunction:: chfl_atom_set_charge

.. doxygenfunction:: chfl_atom_name

.. doxygenfunction:: chfl_atom_set_name

.. doxygenfunction:: chfl_atom_type

.. doxygenfunction:: chfl_atom_set_type

.. doxygenfunction:: chfl_atom_full_name

.. doxygenfunction:: chfl_atom_vdw_radius

.. doxygenfunction:: chfl_atom_covalent_radius

.. doxygenfunction:: chfl_atom_atomic_number

.. doxygenfunction:: chfl_atom_set_property

.. doxygenfunction:: chfl_atom_get_property

.. doxygenfunction:: chfl_atom_properties_count

.. doxygenfunction:: chfl_atom_list_properties
Miscelaneous functions
======================

Error handling
--------------

Apart from the *constructor* functions (the functions returning pointers to the
chemfiles types); all the functions return a status code from the
``chfl_status`` enum:

.. doxygenenum:: chfl_status

.. doxygenfunction:: chfl_last_error

.. doxygenfunction:: chfl_clear_errors

Format list and Metadata
------------------------

.. doxygenfunction:: chfl_guess_format

.. doxygenfunction:: chfl_formats_list

.. doxygenstruct:: chfl_format_metadata
    :members:

Warnings
--------

The chemfiles library send warnings when it encounters malformed files, or any
other condition that the user might want to know about. By default, these
warnings are printed to the standard error stream.
:cpp:func:`chfl_set_warning_callback` allow to redirect these warning by giving
it a callback function to be called on each warning event.

.. doxygentypedef:: chfl_warning_callback

.. doxygenfunction:: chfl_set_warning_callback
.. _capi-frame:

``CHFL_FRAME``
--------------

.. doxygentypedef:: CHFL_FRAME

.. only:: html

    Here is the full list of functions acting on :cpp:type:`CHFL_FRAME`:

    - :cpp:func:`chfl_frame`
    - :cpp:func:`chfl_frame_copy`
    - :cpp:func:`chfl_frame_atoms_count`
    - :cpp:func:`chfl_frame_resize`
    - :cpp:func:`chfl_frame_add_atom`
    - :cpp:func:`chfl_frame_remove`
    - :cpp:func:`chfl_frame_positions`
    - :cpp:func:`chfl_frame_velocities`
    - :cpp:func:`chfl_frame_has_velocities`
    - :cpp:func:`chfl_frame_add_velocities`
    - :cpp:func:`chfl_frame_set_cell`
    - :cpp:func:`chfl_frame_set_topology`
    - :cpp:func:`chfl_frame_add_bond`
    - :cpp:func:`chfl_frame_bond_with_order`
    - :cpp:func:`chfl_frame_remove_bond`
    - :cpp:func:`chfl_frame_clear_bonds`
    - :cpp:func:`chfl_frame_add_residue`
    - :cpp:func:`chfl_frame_step`
    - :cpp:func:`chfl_frame_set_step`
    - :cpp:func:`chfl_frame_guess_bonds`
    - :cpp:func:`chfl_frame_distance`
    - :cpp:func:`chfl_frame_angle`
    - :cpp:func:`chfl_frame_dihedral`
    - :cpp:func:`chfl_frame_out_of_plane`
    - :cpp:func:`chfl_frame_set_property`
    - :cpp:func:`chfl_frame_get_property`
    - :cpp:func:`chfl_frame_properties_count`
    - :cpp:func:`chfl_frame_list_properties`

    --------------------------------------------------------------------

.. doxygenfunction:: chfl_frame

.. doxygenfunction:: chfl_frame_copy

.. doxygenfunction:: chfl_frame_atoms_count

.. doxygenfunction:: chfl_frame_resize

.. doxygenfunction:: chfl_frame_add_atom

.. doxygenfunction:: chfl_frame_remove

.. doxygenfunction:: chfl_frame_positions

.. doxygenfunction:: chfl_frame_velocities

.. doxygenfunction:: chfl_frame_has_velocities

.. doxygenfunction:: chfl_frame_add_velocities

.. doxygenfunction:: chfl_frame_set_cell

.. doxygenfunction:: chfl_frame_set_topology

.. doxygenfunction:: chfl_frame_add_bond

.. doxygenfunction:: chfl_frame_bond_with_order

.. doxygenfunction:: chfl_frame_remove_bond

.. doxygenfunction:: chfl_frame_clear_bonds

.. doxygenfunction:: chfl_frame_add_residue

.. doxygenfunction:: chfl_frame_step

.. doxygenfunction:: chfl_frame_set_step

.. doxygenfunction:: chfl_frame_guess_bonds

.. doxygenfunction:: chfl_frame_distance

.. doxygenfunction:: chfl_frame_angle

.. doxygenfunction:: chfl_frame_dihedral

.. doxygenfunction:: chfl_frame_out_of_plane

.. doxygenfunction:: chfl_frame_set_property

.. doxygenfunction:: chfl_frame_get_property

.. doxygenfunction:: chfl_frame_properties_count

.. doxygenfunction:: chfl_frame_list_properties
Contributing
============

Chemfiles is an open-source project, that I make on my free time. Any
contribution, from documentation improvement to new features including bug fixes
are very welcome!

I would like to help, what can I do ?
-------------------------------------

Lot of things! Just pick one considering the time you can spend, and your
technical skills. And do not hesitate to come by out `gitter chat room`_ to say
hello and ask any question you can have!

.. _gitter chat room: https://gitter.im/chemfiles/chemfiles

Improving the code
^^^^^^^^^^^^^^^^^^

You can pick any `issue`_ in the list. `Help wanted`_ issues are specially
directed at first time contributors, and comes with step by step explanation of
how to solve the issue.

If you plan to add a new feature which is not in the issue list, please open a
new issue so that every one knows you are working on it, and so that the
implementation strategy can be discussed!

.. _issue: https://github.com/chemfiles/chemfiles/issues
.. _Help wanted: https://github.com/chemfiles/chemfiles/labels/Help%20wanted

Improve documentation
^^^^^^^^^^^^^^^^^^^^^

This documentation try to be easy to use, but there is always room for
improvements.  You can easily edit any :file:`.rst` file on `the github
repository <https://github.com/chemfiles/chemfiles/tree/master/doc>`_, and
propose your changes even with no git knowledge. All you need is a Github
account.

Share and spread the word
^^^^^^^^^^^^^^^^^^^^^^^^^

If you think that chemfiles is awesome, share it! The more users it will have,
the more features we can add to it, together!

External project used
---------------------

A few external projects are used in chemfiles developement, and you will need a
bit of knowledge of them to contribute. Depending on what you want to do, not
all these projects are needed.

- Source code versionning: `git`_, together with `github`_ web interface for
  issues and pull requests.
- documentation: `sphinx`_ for generating HTML and PDF documentation, `Doxygen`_
  for documenting the source code, and `breathe`_ to use Doxygen content in
  sphinx.
- Build system: `cmake`_ is used as a cross-plateform, cross-build system
  generator.
- Automatic testing: `Catch2`_ provide an nice unit test framework, and `travis`_
  run these tests each time the code is pushed to the repository.


.. _git: https://git-scm.com/
.. _github: https://github.com/about/
.. _sphinx: http://www.sphinx-doc.org/
.. _Doxygen: http://doxygen.org/
.. _breathe: https://breathe.readthedocs.io/
.. _cmake: https://cmake.org/
.. _Catch2: https://github.com/catchorg/Catch2
.. _travis: https://travis-ci.org/chemfiles/chemfiles/


Coding style, and other formating issues
----------------------------------------

When writing code for chemfiles, please respect the overall coding style. This
is not only a question of style, but make it easier to enter in the project if
all the files are formatted consistently. Do not use non standard features of
any programming language unless there is no other way to do it. In that case,
please wrap the code in ``#ifdef`` macros, and add cases for at least Linux,
macOS and Windows. You should add unit tests for each new piece of code, they
will be run on `travis`_ before your changes are merged.

Git messages should be informatives, and describe the *what*, not the *how*. If
your commit concern the documentation, please add the ``[doc]`` at the beggining
of the message.

Happy coding!
Chemfiles internals
===================

.. toctree::
   :maxdepth: 2

   file
   format

Sources organization
--------------------

You will find the following directory in chemfiles source code:

- ``cmake``: CMake modules used for build configuration;
- ``doc``: the source for this documentation;
- ``examples``: usage examples for C and C++ interfaces;
- ``external``: external libraries used by chemfiles;
- ``include``: the headers of chemfiles;
- ``scripts``: some python and bash scripts used in development.
- ``src``: the sources of the library;
- ``tests``: the sources of the unit tests;

Classes organization
--------------------

Chemfiles is written in C++11, in an object-oriented fashion. A ``Trajectory``
is built on the top of two other private classes: a :ref:`File <class-File>` and
a :ref:`Format <class-Format>`. These are pure abstract class defining the
interface for reading and writing data.

Adding new formats and tweaking behavior of existing formats should be done
either in the ``File`` implementation for everything related to interactions
with the actual file, or in the ``Format`` implementation for everything related
with parsing data from the file.

.. toctree::
   :maxdepth: 2

   file
   format

Every ``Format`` class can be associated to an extension and a format name, the
associations are managed by the ``FormatFactory`` class. New file and formats
should be registered with this class, by specializing the
:cpp:func:`chemfiles::format_metadata` template and calling
:cpp:func:`chemfiles::FormatFactory::add_format`.

.. doxygenclass:: chemfiles::FormatFactory
    :members:

.. doxygenfunction:: chemfiles::format_metadata
.. _class-File:

``File`` classes
================

The ``File`` classes provide abstraction of the IO operation, allowing for the
same format to be used with on-disk files, network files, memory-mapped files,
compressed files, *etc.*

Interface for files
-------------------

.. doxygenclass:: chemfiles::File
    :members:
    :protected-members:

.. doxygenclass:: chemfiles::TextFile
    :members:
    :protected-members:

Implemented classes
-------------------

These classes implement the ``File`` interface defined previously.

.. doxygenclass:: chemfiles::PlainFile
    :members:

.. doxygenclass:: chemfiles::GzFile
    :members:

.. doxygenclass:: chemfiles::XzFile
    :members:

.. doxygenclass:: chemfiles::NcFile
    :members:

.. doxygenclass:: chemfiles::TNGFile
    :members:

.. doxygenclass:: chemfiles::XDRFile
    :members:.. _class-Format:

``Format`` classes
==================

Interface for formats
---------------------

.. doxygenclass:: chemfiles::Format
    :members:

.. doxygenclass:: chemfiles::TextFormat
    :members:

Implemented formats
-------------------

These classes implement the format interface defined previously.

.. binary/non text based formats

.. doxygenclass:: chemfiles::Amber

.. doxygenenum:: chemfiles::AmberFormat

.. doxygenclass:: chemfiles::TNGFormat

.. doxygenclass:: chemfiles::TinkerFormat

.. doxygenclass:: chemfiles::Molfile

.. doxygenenum:: chemfiles::MolfileFormat

.. doxygenclass:: chemfiles::CIFFormat

.. doxygenclass:: chemfiles::CMLFormat

.. doxygenclass:: chemfiles::mmCIFFormat

.. doxygenclass:: chemfiles::MMTFFormat

.. doxygenclass:: chemfiles::TRRFormat

.. doxygenclass:: chemfiles::XTCFormat

.. text based formats

.. doxygenclass:: chemfiles::XYZFormat

.. doxygenclass:: chemfiles::CSSRFormat

.. doxygenclass:: chemfiles::GROFormat

.. doxygenclass:: chemfiles::LAMMPSDataFormat

.. doxygenclass:: chemfiles::MOL2Format

.. doxygenclass:: chemfiles::PDBFormat

.. doxygenclass:: chemfiles::SDFFormat

.. doxygenclass:: chemfiles::SMIFormat
.. _class-Selection:

Selection
=========

.. doxygenclass:: chemfiles::Selection
    :members:

.. doxygenclass:: chemfiles::Match
    :members:
.. _class-Trajectory:

Trajectory
==========

.. doxygenclass:: chemfiles::Trajectory
    :members:
.. _class-Atom:

Atom
====

.. doxygenclass:: chemfiles::Atom
    :members:
.. _cpp-api:

C++ interface reference
=======================


The full public C++ interface is contained in the ``chemfiles.hpp`` header, which
should be included in all the programs using chemfiles. All the classes are in the
``chemfiles`` namespace, so you may want to use ``using namespace chemfiles;`` at the
begining of your program in order to bring all the class in the main namespace.

.. toctree::
   :maxdepth: 2

   trajectory
   frame
   topology
   residue
   atom
   unitcell
   selection
   property
   misc
   helpers

.. _exported-macro:

Exported macro
--------------

In addition to these functions and classes, chemfiles header also export some macro, using
the ``CHEMFILES_`` prefix. The following macro are defined:

.. doxygendefine:: CHEMFILES_VERSION_MAJOR

.. doxygendefine:: CHEMFILES_VERSION_MINOR

.. doxygendefine:: CHEMFILES_VERSION_PATCH

.. doxygendefine:: CHEMFILES_VERSION
.. _class-Topology:

Topology
========

.. doxygenclass:: chemfiles::Topology
    :members:

Connectivity elements
=====================

.. doxygenclass:: chemfiles::Bond
    :members:

.. doxygenclass:: chemfiles::Angle
    :members:

.. doxygenclass:: chemfiles::Dihedral
    :members:

.. doxygenclass:: chemfiles::Improper
    :members:
.. _class-Frame:

Frame
=====

.. doxygenclass:: chemfiles::Frame
    :members:
Helper classes
==============

Chemfiles contains C++11 implementations of various goodies from more recent
versions of the language, or support libraries. Here is a small documentation on
how to use these types.

``optional<T>``
---------------

.. cpp:class:: template <class T> chemfiles::optional

    A class representing optional values. Basic usage example is given here, you
    can also refer to the documentation for ``std::optional`` on
    `cppreference`_.

    For most purposes, one can use an ``optional<T>`` value as if it was a
    pointer ``T*``, with ``nullopt`` indicating the absence of a value.

    .. code-block:: cpp

        // nullopt is used to indicate the optional data is missing
        chemfiles::optional<double> optional = nullopt;

        if (optional) {
            // chemfiles::optional is convertible to bool
        }

        // setting a value
        optional = 78.0;

        // extracting the value
        double a = *optional;
        double b = optional.value();

        // specifying a default value to be used if data is missing
        double c = optional.value_or(-1);


.. _cppreference: https://en.cppreference.com/w/cpp/utility/optional

``span<T>``
-----------

.. cpp:class:: template <class T> chemfiles::span

    A ``span<T>`` is a view inside a ``std::vector<T>`` providing all the
    operations of a vector except for memory allocation. The idea and
    implementation comes from the `GSL`_ library.

    .. code-block:: cpp

        // you should never need to create a span yourself
        chemfiles::span<double> span = /* ... */;

        // span supports range-based iteration
        for (auto a: span) {
            // ...
        }

        // but also the usual indexing
        for (size_t i=0; i<span.size(); i++) {
            double b = span[i];
        }

        // and direct indexing
        double b = span[5];
        span[6] = 78.3;

        // you can also use a span with standard algorithms
        auto sum = std::accumulate(span.begin(), span.end(), 0.0);


.. _GSL: https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#S-gsl
Miscelaneous
============

Basic types
-----------

.. doxygenclass:: chemfiles::Vector3D
    :members:

.. doxygenfunction:: chemfiles::dot

.. doxygenfunction:: chemfiles::cross

.. doxygenclass:: chemfiles::Matrix3D
    :members:

Format list and Metadata
------------------------

.. doxygenfunction:: chemfiles::guess_format

.. doxygenfunction:: chemfiles::formats_list

.. doxygenclass:: chemfiles::FormatMetadata
    :members:

Errors handling
---------------

In chemfiles, any error will throw an exception. Any program using chemfiles
should thus wrap any call in a ``try ... catch`` block, and handle the errors in
any pertinent way.

.. code-block:: cpp

    #include <iostream>
    #include <chemfiles.hpp>

    int main() {
        try {
            chemfiles::Trajectory file("filename.xyz");
            auto frame = file.read();
            auto positions = frame.positions();
            // Do something here
        } catch (const chemfiles::Error& e) {
            // Basic error handling logging error to stdout
            std::cout << "Error in chemfiles:" << e.what() << std::endl;
            return -1;
        }

        return 0;
    }

Any exceptions thrown by chemfiles will derive from ``chemfiles::Error``.
Catching  ``chemfiles::Error`` will then catch any exception thrown by
chemfiles. You also can catch any other error if you need finer grain control.
``chemfiles::Error`` derives from ``std::runtime_error``, so it should play
nicely with any existing C++ error handling.

.. doxygenstruct:: chemfiles::Error
    :members:

.. doxygenstruct:: chemfiles::FileError
    :members:

.. doxygenstruct:: chemfiles::MemoryError
    :members:

.. doxygenstruct:: chemfiles::FormatError
    :members:

.. doxygenstruct:: chemfiles::SelectionError
    :members:

.. doxygenstruct:: chemfiles::OutOfBounds
    :members:

.. doxygenstruct:: chemfiles::PropertyError
    :members:

Warnings
--------

Chemfiles send warnings when it encounters malformed files, or any other
condition that the user might want to know about. By default, these warnings are
printed to the standard error stream. :cpp:func:`chemfiles::set_warning_callback`
allow to redirect these warning by giving it a callback function to be called on
each warning event.

.. doxygenfunction:: chemfiles::set_warning_callback

.. doxygentypedef:: chemfiles::warning_callback_t
.. _class-Residue:

Residue
=======

.. doxygenclass:: chemfiles::Residue
    :members:
.. _class-UnitCell:

UnitCell
========

.. doxygenclass:: chemfiles::UnitCell
    :members:
.. _class-Property:

Property
========

.. doxygenclass:: chemfiles::Property
    :members:
