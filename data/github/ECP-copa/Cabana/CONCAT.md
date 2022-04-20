# Change Log

## (dev)

- Updated minimum Kokkos dependency to version 3.4

## 0.5.0

**New Features**

- Particle migration using Cajita grid added
- Random particle generation added
- Complete Cajita tutorial examples added
- Cajita performance benchmarks added

**Bug Fixes and Improvements**

- Remove all uses of `Kokkos::Impl`
- Redesign `SimdPolicy` to not modify the underlying `Kokkos::TeamPolicy`
- Rename `Cabana_REQUIRE_`{`PTHREAD` -> `THREADS`}
- Rename clang-format build rule `format` -> `cabana-format`
- Improved Doxygen coverage
- Improved wiki documentation

**Minimum dependency version updates**

- CMake minimum 3.16 required (previously 3.9)
- Optional dependency heFFTe minimum 2.1 (previously 2.0)
- Optional dependency HYPRE minimum 2.22.1 (previously 2.22.0)

**Experimental Features (subject to change in future releases)**

- Distributed particle output with SILO library interface
- Cajita load balancing added through ALL library interface

## 0.4.0

**New Features**

- C++14 required
- Updated minimum Kokkos dependency to version 3.2
- AMD HIP support and continuous integration testing
- Intel SYCL support and continuous integration testing
- OpenMP-Target support (with some disabled features)
- Hybrid particle-grid capability through the Cajita interfaces. Features include:
    - 2D/3D structured grid data structures
    - particle-grid interpolation
    - particle-grid communication
    - multidimensional distributed FFTs via heFFTe (including host, CUDA, and HIP)
    - linear solvers and preconditions via HYPRE (including host and CUDA)

**Bug Fixes and Improvements**

- Removed deprecated portability macros in favor of Kokkos macros (e.g. KOKKOS_INLINE_FUNCTION)
- General performance improvements including neighbor list and particle communication updates
- Improved Doxygen coverage, wiki documentation, and tutorials

**Experimental Features (subject to change in future releases)**

- Sparse grids support in Cajita
- Structured grid data I/O in Cajita

## 0.3.0

**New Features**

- Updated minimum Kokkos dependency to version 3.1
- CUDA and HIP support and testing in continuous integration
- Mirror view capability for AoSoA
- New performance benchmarks for sorting, communication, and neighbor lists
- Improving AoSoA memory managment with empty() and shrinkToFit()
- Second level neighbor parallel for and reduce algorithms for triplet operations
- Unmanaged AoSoA for wrapping user memory

**Bug Fixes and Improvements**
- Using new CMake target for linking Kokkos
- Removed numerous instances of default allocation of Kokkos Views
- Eliminated use of user-defined MPI tags in communication algorithms
- Cleaned usage of deprecated Kokkos code
- Update for compilation with C++14
- Significant performance enhancements to communication code

**Experimental Features (subject to change in future releases)**
- Tree-based neighbor lists using ArborX


## 0.2.0

**New Features**

- An optional MPI dependency has been added. Note that when CUDA is enabled the MPI implementation is expected to be CUDA-aware. [#45](https://github.com/ECP-copa/Cabana/pull/45)
- Particle redistribution via MPI. Implemented in the `Cabana::Distributor` [#43](https://github.com/ECP-copa/Cabana/pull/43)
- Particle halo exchange via MPI. Implemented in the `Cabana::Halo` [#43](https://github.com/ECP-copa/Cabana/pull/43)
- Parallel for concept for 2D indexing in AoSoA loops. Implemented via `Cabana::simd_parallel_for`. Includes a new execution space concept `Cabana::SimdPolicy`. [#49](https://github.com/ECP-copa/Cabana/pull/49)
- Parallel for concept for traversing neighbor lists. Implemented via `Cabana::neighbor_parallel_for` [#49](https://github.com/ECP-copa/Cabana/pull/49)
- Continuous integration for pull requests via GitHub [#9](https://github.com/ECP-copa/Cabana/pull/9)
- Support the ECP continuous integration infrastructure [#66](https://github.com/ECP-copa/Cabana/pull/66)
- New example using scafacos for long-range solvers [#46](https://github.com/ECP-copa/Cabana/pull/46)
- Additional tutorials documentation on the [Wiki](https://github.com/ECP-copa/Cabana/wiki)

**Bug Fixes and Improvements**
- Fixed a bug in the construction of slices on uninitialized `AoSoA` containers [#80](https://github.com/ECP-copa/Cabana/pull/80)
- Construct Verlet lists over the specified range of indices [#70](https://github.com/ECP-copa/Cabana/pull/70)
- Removed aliases of Kokkos macros and classes [#58](https://github.com/ECP-copa/Cabana/pull/58)

## 0.1.0

**New Features**

- Core portable data structures: `Tuple`, `SoA`, `AoSoA`, and `Slice`
- `DeepCopy`: deep copy of data between `AoSoA` data structures
- Sorting and binning `AoSoA` data structures
- `LinkedCellList`: linked cell list implementation
- Portable neighbor list interface
- `VerletList` - linked cell accelerated Verlet list implementation
- Basic tutorial

**Experimental Features (subject to change in future releases)**

- `parallel_for` - portable parallel loops over supported execution spaces
- `neighbor_parallel_for` - portable parallel loops over neighbor lists in supported execution spaces
- `RangePolicy` - defines an index range for parallel loops
# Contributing

Contributing to Cabana is easy: just open a [pull
request](https://help.github.com/articles/using-pull-requests/). Make
`master` the destination branch on the [Cabana
repository](https://github.com/ECP-copa/Cabana) and allow edits from
maintainers in the pull request.

Your pull request must pass Cabana's tests, which includes using the coding
style from `.clang-format` (tested with clang-10), and be reviewed by at least
one Cabana developer.

Other coding style includes:
* Camel case template parameters (`NewTemplateType`)
* Camel case class names (`NewClassName`)
* Lower camel case function names (`newFunctionName`)
  * Note: there are some exceptions to match Kokkos (e.g. `Cabana::deep_copy`
    and `Cabana::neighbor_parallel_for`)
* Lower case, underscore separated variables (`new_variable_name`)
* Class members which are `private` are preceeded by an underscore (`_private_class_variable`)
* Class/struct member type aliases use lower case, underscore separated names (`using integer_type = int;`)
# CoPA Cabana - The Exascale Co-Design Center for Particle Applications Toolkit

Cabana is a performance portable library for particle-based simulations.
Applications include, but are not limited to, Molecular Dynamics (MD) with
short- and/or long-range interactions; various flavors of Particle-in-Cell
(PIC) methods, including use within fluid and solid mechanics and plasma
physics; and N-body cosmology simulations.
Cabana provides particle data structures, algorithms, and utilities to enable
simulations on a variety of platforms including many-core architectures and
GPUs.

Cabana is developed as part of the Co-Design Center for Particle Applications
(CoPA) within the Exascale Computing Project (ECP) under the U.S. Department
of Energy. CoPA is a multi-institutional project with developers from ORNL,
LANL, SNL, LLNL, PPNL, and ANL.

## Documentation

Instructions for building Cabana on various platforms, an API reference with
tutorial links, and links to the Doxygen can be found in our
[wiki](https://github.com/ECP-copa/Cabana/wiki).

For Cabana-related questions you can open a GitHub issue to interact with the
developers.

## Contributing

We encourage you to contribute to Cabana! Please check the
[guidelines](CONTRIBUTING.md) on how to do so.

## License

Cabana is distributed under an [open source 3-clause BSD license](LICENSE).
