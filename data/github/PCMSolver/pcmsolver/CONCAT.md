## Individual Contributors

- Roberto Di Remigio (@robertodr)
- Luca Frediani (@ilfreddy)
- Monica Bugeanu (@mbugeanu)
- Arnfinn Hykkerud Steindal (@arnfinn)
- Radovan Bast (@bast)
- Lori A. Burns (@loriab)
- T. Daniel Crawford (@lothian)
- Krzysztof Mozgawa
- Ville Weijo (@vweijo)
- Ward Poelmans (@wpoely86)

This list was obtained 2020-11-29 by running `git shortlog -sn`
# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into PCMSolver.
Our contribution guide is based on [Psi4 contribution guide](https://github.com/psi4/psi4/blob/master/.github/CONTRIBUTING.md)

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free).
* [Fork](https://help.github.com/articles/fork-a-repo/) the
  [PCMSolver/pcmsolver](https://github.com/PCMSolver/pcmsolver) repository on GitHub.
* On your local machine,
  [clone](https://help.github.com/articles/cloning-a-repository/) your fork of
  the PCMSolver repository.

## Making Changes

* Add some really awesome code to your local fork.  It's usually a [good
  idea](http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/)
  to make changes on a
  [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/)
  with the branch name relating to the feature you are going to add.
* When you are ready for others to examine and comment on your new feature,
  navigate to your fork of PCMSolver on GitHub and open a
  [pull request](https://help.github.com/articles/using-pull-requests/) (PR)
  __towards the `master` branch__.
  Note that after you launch a PR from one of your fork's branches, all
  subsequent commits to that branch will be added to the open pull request
  automatically.
  Each commit added to the PR will be validated for mergeability, compilation
  and test suite compliance; the results of these tests will be visible on the
  PR page.
* The title of your pull request should be marked with `[WIP]` if it’s a work
  in progress and with `#trivial` if it is a set of trivial changes.
* If you're providing a new feature, you must add test cases, documentation and
  update the `CHANGELOG.md` file.
* When the code is ready to go, make sure you run the full or relevant portion
  of the test suite on your local machine to check that nothing is broken.
* When you're ready to be considered for merging, check the "Ready to go" box
  on the PR page to let the PCMSolver team know that the changes are complete.
  The code will not be merged until this box is checked, the continuous
  integration (Travis CI) returns checkmarks, and multiple core
  developers give "Approved" reviews.

## Pull Request Requirements

The project is integrated with [Danger.Systems](http://danger.systems/ruby/).
On each PR, one CI job will run the integration and a [bot](https://github.com/minazobot) will
report which requirements are **not met** in your PR.
These reports can be _warnings_ and _errors_. You will discuss and solve both
of them with the reviewers.
The automatic rules are laid out in the `Dangerfile` and are used to enforce an
adequate level of testing, documentation and code quality.

### Danger.Systems Warnings

* PRs classed as Work in Progress.
* Codebase was modified, but no tests were added.
* Nontrivial changes to the codebase, but no documentation added.
* Codebase was modified, but `CHANGELOG.md` was not updated.
* Source files were added or removed, but `.gitattributes` was not updated.

### Danger.Systems Errors

* Commit message linting, based on some of [these recommendations](https://chris.beams.io/posts/git-commit/):
  * Commit subject is more than one word.
  * Commit subject is no longer than 50 characters.
  * Commit subject and body are separated by an empty line.
* Clean commit history, without merge commits.
* Code style for `.hpp`, `.cpp`, `.h` files follows the conventions in
  `.clang-format`.

## Licensing

We do not require any formal copyright assignment or contributor license
agreement.
**Any contributions intentionally sent upstream are presumed to be offered under
terms of the OSI-approved LGPLv3 License.**

## Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)
[![DOI](https://zenodo.org/badge/23794148.svg)](https://zenodo.org/badge/latestdoi/23794148)
![Build and test PCMSolver](https://github.com/PCMSolver/pcmsolver/workflows/Build%20and%20test%20PCMSolver/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/pcmsolver/badge/?version=stable)](http://pcmsolver.readthedocs.org/en/latest/?badge=latest)
[![Coverage Status](https://codecov.io/gh/PCMSolver/pcmsolver/branch/release%2F1.2.Z/graph/badge.svg)](https://codecov.io/gh/PCMSolver/pcmsolver)

PCMSolver
=========

An API for the Polarizable Continuum Model. Copyright [Roberto Di Remigio](mailto:roberto.d.remigio@uit.no),
[Luca Frediani](mailto:luca.frediani@uit.no) and [contributors](https://github.com/PCMSolver/pcmsolver/blob/release/1.2.Z/AUTHORS.md)

- [Project website](https://github.com/PCMSolver/pcmsolver)
- [Changelog](CHANGELOG.md)
- [Documentation](http://pcmsolver.readthedocs.io)
- [Build and test history](https://travis-ci.org/PCMSolver/pcmsolver/builds)
- Licensed under [LGPLv3](LICENSE)
- CMake infrastructure managed *via* [Autocmake](http://autocmake.readthedocs.io/)
# Change Log

## [Version 1.3.0] - 2020-11-30

### Added

- An implementation of the fluctuating charge (FQ) molecular mechanics (MM) model.
- A new API initialization function `pcmsolver_new_read_host`. This allows
  creating an *uninitialized* `pcmsolver_context_t` object: the initialization is
  deferred to a later point and orchestrated by the host program.  The input
  options are then set *via* the other API functions:
  * `pcmsolver_set_bool_option`
  * `pcmsolver_set_int_option`
  * `pcmsolver_set_double_option`
  * `pcmsolver_set_string_option`
  A call to `pcmsolver_refresh` will trigger initialization of all internal
  objects managed through `pcmsolver_context_t`.  Use this function when you want
  to manage input parsing on the host side, rather than use a dedicated input file
 for PCMSolver.
- A `pcmsolver_fill_pcminput` function in the Fortran interface. This function
  simplifies the creation of the `PCMInput` object from Fortran.
- A `PCMSOLVER_WARNING` macro to emit warning to users.

### Changed

- Use GitHub Actions instead of Travis for continuous integration.
- Remove integration with Danger in CI.
- A C++11-compliant compilers is **required** to compile the library.
- The Fortran interface has been refactored to hide the transformation of
  strings from Fortran to C from the host program. 
- The check on the symmetric positive definiteness of the S matrix triggers a
  warning, not an error.
- We switched back to using LU decomposition for CPCM.
- For IEF, we switched to LU decomposition, instead of partially pivoted LU.
- UFF radii are scaled by 1.1 by default.

## Fixed

- The printout of sphere centers and radii always uses Angstrom.

## [Version 1.2.3] - 2019-02-14

### Added

- A `pcmsolver_citation` function was added which accepts a function pointer to
  printing facilities used in the host program. `pcmsolver_citation` will print
  versioning and citation information for the library. Please, include a call to
  this function and print the returned string somewhere in your code, so that
  users are aware of what literature to cite when using the library.

### Changed

- The list of finite elements return from the PEDRA Fortran code is now pruned
  to remove finite elements with area less than 1.0e-4 This avoids numerical
  artifacts in the formation of the PCM matrices.
  The number of pruned finite elements is reported when printing the cavity.
- The `pcmsolver_print` function now only prints out the set up for the PCM calculation.

## [Version 1.2.2] - 2019-01-27

### Changed

- Perform Gauss' theorem check on nuclear ASC upon initialization. In case a
  discrepancy larger than 10^-3 is detected in the absolute difference between the
  Gauss' theorem and computed values, calculation is aborted as this points to
  numerical problems, most likely in the construction of the cavity.
- Add check and test for S matrix positive-definiteness. Certain cavity
  discretizations lead to the S matrix not being symmetric positive-definite. We
  abort if that's the case.

## [Version 1.2.1] - 2018-05-02

### Fixed

- Import of local copy of the `pyparsing` module in Getkw. Thanks @loriab for
  the fix.
- The `conf.py` documentation build script is restored to normal operation.

## [Version 1.2.0] - 2018-04-27

### Added

- Green’s function for a spherical nanoparticle with **real** permittivity.
  The Green's function is known in analytical form from the work of
  Messina (J. Chem. Phys. 2002, 117 (24), 11062) and
  Delgado _et al._ (J. Chem. Phys. 2013, 139 (2), 024105)

### Deprecated

- C++03 support is effectively **deprecated** in favor of C++11. Support for GCC
  4.6 and earlier has been dropped. Please consider upgrading your C++ compiler
  to a [fully standard-compliant one](https://en.cppreference.com/w/cpp/compiler_support#cpp11)

### Changed

- The versioning machinery has been updated. The update was inspired by the
  `versioner.py` tool devised by @loriab.
- The documentation building scripts `conf.py` and `cloc_tools.py` have been
  thoroughly refactored and simplified:
  * `cloc_tools.py` is now called `cloc_wrapper.py`. It is no longer configured
  and can be found in the `cloc_tools` subfolder of `doc`.
  * A local build of documentation will only work if run from within the `doc` folder:
  ```
  sphinx-build . _build
  ```
  This choice was made to simplify the set up of the ReadTheDocs and local
  documentation building procedures and to minimize the chances of breaking
  either.
- The Fortran API bindings file `pcmsolver.f90` is now installed alongside the
  `pcmsolver.h` header file. The users will have to compile it explicitly to
  get the type checking from the API redeclaration in Fortran 90.
  The file is always installed.
- The `ENABLE_Fortran_API` configuration option has been renamed
  `TEST_Fortran_API`, since it now only triggers compilation of the
  `Fortran_host` test case.
- **BREAKING CHANGE** The minimum required version of CMake is now 3.3

### Fixed

- `std::string`-s are now used in the `Meddle` object functions manipulating
  the surface functions map. Explicit casts from `const char *` to `std::string`
  are handled in the API functions.
- Properly enforce `const`-correctness of the `Meddle` object and of its usage
  in the context API.
- [Cholesky decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html) is used
  in CPCMSolver to get the inverse of the S matrix. The robust Cholesky (LDLT)
  previously used is broken with the latest version of the Intel compilers.

## [Version 1.2.0-rc1] - 2018-03-02

### Added

- Double logarithmic scale for the integration of spherical diffuse
  interfaces: much more stable than the previous version, allowing for
  Runge-Kutta 4 integrator.
- A new CMake module `options_wrappers.cmake` that adds new wrapper macros for
  the CMake `option` command.

### Fixed

- Bug in the diffuse interface Green's function. Contrary to the sharp
  interface case, it is wrong to remove the monopole, which becomes
  identically zero when the corresponding differential equation is
  solved in extreme cases (e.g. charge far away from the sphere).
- Visibility of symbols in the shared library is _finally_ handled properly.
  The necessary flags to the C++ compiler were not set for the subtargets built
  as CMake `OBJECT` libraries. This results in a modest decrease in library
  size.
- `RPATH` handling for the standalone executable `run_pcm`.
  For the executable in `<build_dir>/bin` the `RPATH` will contain
  `<build_dir>/lib64` (or `<build_dir>/lib`) as path to `libpcm.so.1`.
  For the executable in `<install_prefix>/bin` the `RPATH` will contain
  `<install_prefix>/lib64` (or `<install_prefix>/lib`) as path to `libpcm.so.1`.
- Code coverage analysis was restored. We now use
  [Codecov](https://codecov.io). Thanks @arnfinn :tada:
- **BREAKING CHANGE** The layout of the installation for the Python scripts has
  been changed. This addresses issue #116. Given the user-defined installation
  prefix (`<prefix>`) one would obtain:
  ```
  <prefix>/pcmsolver
           ├── bin
           │   ├── go_pcm.py
           │   ├── plot_cavity.py
           │   └── run_pcm
           ├── include
           │   └── PCMSolver
           │       ├── bi_operators
           │       ├── cavity
           │       ├── Citation.hpp
           │       ├── Config.hpp
           │       ├── Cxx11Workarounds.hpp
           │       ├── ErrorHandling.hpp
           │       ├── external
           │       ├── GitInfo.hpp
           │       ├── green
           │       ├── interface
           │       ├── LoggerInterface.hpp
           │       ├── PCMInput.h
           │       ├── PCMSolverExport.h
           │       ├── pcmsolver.h
           │       ├── PhysicalConstants.hpp
           │       ├── solver
           │       ├── STLUtils.hpp
           │       ├── TimerInterface.hpp
           │       └── utils
           ├── lib64
           │   ├── libpcm.a
           │   ├── libpcm.so -> libpcm.so.1
           │   ├── libpcm.so.1
           │   └── python
           │       └── pcmsolver
           └── share
               └── cmake
                   └── PCMSolver
  ```

### Changed

- As a result of the visibility change, unit tests can only be linked against
  the static library, since all symbols are always visible in a static archive
  library.
- **BREAKING CHANGE** The `pcmsolver.py` script/module **was removed** to
  address issues #111 and #112. The `main` portion of the script has now been
  separated into a `go_pcm.py` script which can be used in the same way as
  `pcmsolver.py` was used before: to parse the PCMSolver input and to run the
  module standalone. This separation greatly simplifies the use as a module
  within the Python launcher scripts of host programs.

## [Version 1.1.12] - 2018-01-20

### Added

- The function `pcmsolver_fstring_to_carray` was added to the Fortran bindings
  for the library. As the name suggests, this function "translates" a Fortran
  string into a C `char` array. This supersedes and replaces the
  `pcmsolver_f2c_string` function.
- An updated version of the `pcmsolver_new` API function has been introduced.
  The full name of the function is `pcmsolver_new_v1112` (to avoid API breakage)
  and accepts as additional argument the name of the parsed input file
  generated by `pcmsolver.py`. This allows to decoupled parsing and writing out
  the parsed file, in case directory-specific operations have to be performed
  by the host program.

### Changed

- Dispatching of Green’s functions in the factory is purely label-based. The
  input reader generates a string that holds the type of the Green’s function,
  the strategy for calculating its normal derivatives and the dielectric
  profile (if needed): `TYPE_DERIVATIVE_PROFILE`. The Green’s functions are
  subscribed to the factory using this label, which is then also the one used
  for their creation. Iteration over _lists of types_ to dispatch the correct
  template parameters is thus completely avoided. This removes the dependency
  on the metafunctions in `ForId.hpp` and the dependency on Boost.MPL, the
  Boost metaprogramming library.
- Dispatching of quadrature rules for the numerical boundary integral operator
  has been rewritten, removing the dependency on Boost.MPL.
- The input wrapping data `struct`-s have been modified to include the string
  defining the type of object and to remove now unused information.
- Documentation building is fully handled _via_ `sphinx-build`: CMake will no longer generate a `doc` build target.
- Simplified `.travis.yml` and got rid of Conda to handle multiple Python versions.
- Moved python scripts from `bin/` to `bin/pcmsolver/` and introduced `bin/pcmsolver/__init__.py`.
- The `parse_pcm_input` function returns the parsed input as a string without
  saving it to file. The old behavior can be recovered by setting the
  `write_out` argument to `True`.
- The [Catch unit test framework](https://github.com/philsquared/Catch) has
  been updated to its latest version still supporting C++03
  [v1.11.0](https://github.com/philsquared/Catch/releases/tag/v1.11.0)
- Travis continuous integration is _no longer_ run on Mac infrastructure.

### Fixed

- Documentation building on ReadTheDocs is fully functional again, thanks @arnfinn :tada:
  The build had been failing for a while since docs were generated for all files, including
  documentation files from previous build. Besides, source code doxygen blocks were not
  extracted when inside namespaces.
- Compiler warning in the `Sphere` class due to a redeclaration of `operator<<`.
- The `plot_cavity.py` script is now Python 3 compatible.

### Removed

- `FortranCUtils` header and source files were removed. The `pcmsolver_f2c_string` and `pcmsolver_c2f_string`
  functions are thus gone. For the former, use the replacement
  `pcmsolver_fstring_to_carray` function provided in the Fortran bindings to the library.
- Some unused files have been removed:
  * `Interpolation.hpp`
  * `Interpolation.cpp`
  * `Vector2.hpp`
  * `Vector3.hpp`
  * `ForId.hpp`
  * `DerivativeUtils.hpp`
- The `cavityType()`, `solverType()`, `greenInsideType()`, `greenOutsideType()`
  and `integratorType()` functions in the `Input` object have been removed.
  This information is now wrapped into the corresponding input-wrapping
  `struct`-s: `CavityData`, `SolverData`, `GreenData` and `BIOperatorData`.

### Deprecated

- We are in the process of removing the dependency on the Boost libraries. The
  only dependency on Boost allowed in the 1.Y.Z release series will be for
  replacing C++11 functionality that is missing from the C++98 standard. For
  such code, we will require a pure C++11 counterpart.

## [Version 1.1.11] - 2017-10-25

### Added

- A Python script, using Matplotlib, to plot the cavity.
  The script can also color-map the finite elements according to the values of
  a surface function.
- The input learnt to parse the additional `ChargeDistribution` section.
  It is possible to specify a classical charge distribution of point multipoles.
  This can be an additional source of electrostatic potential for the calculation
  of the ASC.
- Restored compilation for g++ < v5.1.
- [Ninja](https://ninja-build.org/) can be used as a generator.
  Notice that at least [CMake 3.7.2](https://cmake.org/cmake/help/v3.7/generator/Ninja.html#fortran-support)
  **and** the [Kitware-maintained](https://github.com/Kitware/ninja) version of
  Ninja are required to successfully compile.

### Changed

- Use [`#pragma once`](https://en.wikipedia.org/wiki/Pragma_once) instead of
  `#ifndef, #define, #endif` to guard against multiple inclusion of header files.
- The uppercased contents of the `.pcm` input file are written to a temporary
  file, instead of overwriting the user provided file. The temporary file is
  removed after it has been parsed. Fixes #91 as noted by @ilfreddy.
- Use Runge-Kutta-Fehlberg 7(8) ODE solver to integrate the radial equation
  in the spherical diffuse Green's function class.

### Fixed

- A bug in the initialization of a restart cavity from old `.npz` files.
  Currently, the `.npz` file saves sphere center, arcs and vertices of each
  finite element. This information is in fact needed to plot the cavity using
  the Python script in `tools`. Older `.npz` files did not contain this
  information and could not be read in. The additional information is read in
  as arrays of zeros in case it is not present on file.

## [Version 1.1.10] - 2017-03-27

### Changed

- Updated the `cloc.pl` script to version 1.72
- Simplified the internal structure of the `Meddle` and `Input` objects.
- Export dependency on Zlib for the static libraries. Thanks @loriab for the pull request
  fixing [a build problem within Psi4](http://forum.psicode.org/t/crc32-undefined-symbol-at-runtime-when-built-with-pcmsolver-gcc-4-9-4/449/7)

## [Version 1.1.9] - 2017-02-16

### Changed

- PCMSolver is now exported as a proper [CMake target](https://cmake.org/cmake/help/v3.0/manual/cmake-buildsystem.7.html)
  See PR #38 for details. Thanks @loriab for the work.
- The Python scripts shipped with the library are now Python 2 and Python 3 compatible.
- `Factory` is no longer implemented as a Singleton.
- The [Catch unit test framework](https://github.com/philsquared/Catch) has
  been updated to its latest version
  [v1.7.2](https://github.com/philsquared/Catch/releases/tag/v1.7.2)
- Updated the version of Eigen bundled with the code.
  The minimum required version of Eigen is still 3.3.0, but we ship
  [Eigen 3.3.2](http://eigen.tuxfamily.org/index.php?title=ChangeLog#Eigen_3.3.2)

### Fixed

- Revert to use [Robust Cholesky decomposition](https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html)
  to compute the inverse of the S matrix in `CPCMSolver`.

## [Version 1.1.8] - 2017-02-06

### Added

- Namespaces for all of the internal code have been introduced.
  The top-level namespace is `pcm`. At finer levels the namespaces have the same
  names as the respective subdirectories. Read the programmers' documentation
  for further details on the use of namespaces in the project.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The top-level, convenience header `BoundaryIntegralOperator.hpp` includes all
  subclasses and utility headers in the `bi_operators` subdirectory.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The top-level, convenience header `Cavity.hpp` includes all subclasses and
  utility headers in the `cavity` subdirectory.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The top-level, convenience header `Green.hpp` includes all subclasses and
  utility headers in the `green` subdirectory.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The top-level, convenience header `Solver.hpp` includes all subclasses and
  utility headers in the `solver` subdirectory.
  Related to issue #34 on [GitHub] and #60 on [GitLab].

### Changed

- The abstract base class for the boundary integral operator integrators has
  been renamed `IBoundaryIntegralOperator`.
  The relevant factory is bootstrapped upon creation of the `Meddle` object,
  i.e. at library initialization.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The abstract base class for the cavities has been renamed `ICavity`.
  The relevant factory is bootstrapped upon creation of the `Meddle` object,
  i.e. at library initialization.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The abstract base class for the solvers has been renamed `ISolver`.
  The relevant factory is bootstrapped upon creation of the `Meddle` object,
  i.e. at library initialization.
  Related to issue #34 on [GitHub] and #60 on [GitLab].
- The `typedef` for numerical differentiation in the Green's function classes
  has been renamed `Stencil` to avoid name clashes with the `Numerical`
  boundary integral operator type.

### Fixed

- A bug in the selection of the extended diagnostics flags for the GNU C++
  compiler. These flags are now enabled only for versions >= 5.1.0 and when the
  C++11 standard is enable. Fixes issue #36 on [GitHub] and #62 on [GitLab].
- A bug in the initialization of the factory for the cavity classes was fixed.
  The bug manifested only in the static library `libpcm.a`
  Fixes issue #34 on [GitHub] and #60 on [GitLab].

## [Version 1.1.7] - 2016-12-01

### Added

- A pre-commit hook in `.githooks/pre-commit-clang-format` checking that the
  format of C++ header and source files conforms to the project style. The hook
  uses `clang-format` and the style mandated by the `.clang-format` file to
  check files in the tree. Commit is rejected if one or more files are non
  compliant. The hook generates a patch and shows the command needed to apply
  it.
  To enable the hooks you need to have a `.git/hooks/pre-commit` file
  containing this line `.githooks/pre-commit`
  _NOT recommended_ The hook can be skipped by passing the `--no-verify` option to `git commit`
- A pre-commit hook in `.githooks/pre-commit-license-maintainer` checking the
  license headers. The hook is configured based on the `.gitattributes` file.
  The hook will check the license headers and amend them, either by updating
  the year and authors information or by adding the whole header.
  To enable the hooks you need to have a `.git/hooks/pre-commit` file
  containing this line `.githooks/pre-commit`
  _NOT recommended_ The hook can be skipped by passing the `--no-verify` option to `git commit`
- An `UNUSED` preprocessor macro to mark arguments as unused.
- An `UNUSED_FUNCTION` preprocessor macro to mark functions as unused.
- A set of preprocessor macros with Git information (`GIT_COMMIT_HASH`,
  `GIT_COMMIT_AUTHOR`, `GIT_COMMIT_DATE` and `GIT_BRANCH`) is automatically
  generated and saved in the `git_info.h` header file.
- An API function to print the contents of a surface function to the host
  program output.
- An API function to get the dipole moment, relative to the origin, due to the ASC
  on the cavity. Both the norm and the components can be obtained.
- `.gitattributes` now instructs Git to ignore binary files in diff operations.
  PNG files are diff-ed using EXIF information. To set this up properly,
  install an EXIF tool on your machine and run `git config diff.exif.textconv
  exiftool` in your local copy of the repository.

### Changed

- The Fortran bindings file has been renamed `pcmsolver.f90`.
- The Green's function, solver and boundary integral operator classes have been
  radically redesigned. This avoids coupling between integrators and Green's
  function that existed in the previous design.
  See the [Green's function code
  reference](http://pcmsolver.readthedocs.io/en/latest/code-reference/greens-functions.html)
  for a more detailed explanation.
- **BREAKING CHANGE** The minimum required version of Eigen is now 3.3.0
  The version bundled with the code has been accordingly updated.
- The `PCMSOLVER_ERROR` macro now takes only one argument and prints out a more
  informative error message.
- Switched to the latest version of
  [Autocmake](http://autocmake.readthedocs.io/) The configuration file is now YAML-based.
  The PyYAML module is thus required.
- The extended diagnostic flags `-Wsuggest-attribute=pure
  -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wsuggest-final-types
  -Wsuggest-final-methods -Wsuggest-override -Wuseless-cast
  -Wunsafe-loop-optimizations` are always set when using the GNU C++ compiler
  in a debug configuration.
- The C++11 compatibility CMake macros now check for the availability of the
  `noreturn` attribute. A workaround macro, accessible _via_ ``, has
  been added to the `Cxx11Workarounds.hpp` header file.
- **BREAKING CHANGE** The ouput flushing function must be passed explicitly as
  a function pointer to the `pcmsolver_new` function during library
  initialization.
  The function pointer has the signature
  `typedef void (*HostWriter)(const char * message)`
  thus accepting a single argument instead of the previous two.
- [GNU standard installation
  directories](http://www.gnu.org/prep/standards/html_node/Directory-Variables.html)
  have been imposed, thanks to work by @loriab.
  Given a prefix, header files are now installed to `include/pcmsolver`,
  executables to `bin`, libraries to `lib` and scripting tools to `share`.
  The install prefix and the installation directories can be specified by the
  `--prefix`, `--bindir`, `--libdir`, `--includedir` and `--datadir` options to
  the `setup.py` script (or the corresponding CMake variables)

## [Version 1.1.6] - 2016-09-20

### Added

- A function returning a molecule object for the water molecule.

### Changed

- [Cholesky decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html) is used
  whenever the inverse of the S matrix has to be calculated.
  The S matrix is self-adjoint, positive-definite and the LLT decomposition is
  faster than LDLT.

### Fixed

- Some inconsistencies in input reading from host and a related memory leak in the radii
  initialization.

## [Version 1.1.5] - 2016-07-19

### Added

- A radii set derived from [Allinger's MM3 model](http://dx.doi.org/10.1016/S0166-1280(09)80008-0)
  can now be chosen to build the van der Waals cavity surface.
  Notice that the values reported in the original paper are **divided by** 1.2, to match the
  default radii set used in [ADF](https://www.scm.com/doc/ADF/Input/COSMO.html)
  The closest match to ADF can be obtained by using CPCM as solver, Allinger's radii and setting
  the scaling of radii to false.

## [Version 1.1.4] - 2016-07-05

### Changed

- The `CPCMSolver` object now stores the scaled, Hermitian, symmetry-adapted S matrix.
  Polarization weights are then directly computed from the incoming MEP.
- The `IEFSolver` object now stores the non-Hermitian, symmetry-adapted T and R matrices.
  The explicit calculation of the inverse of T is thus avoided.
  However, two square matrices of size equal to the cavity size are stored instead
  of just one. To obtain the polarization weights _two_ linear systems of equations are solved.
  A [partially pivoted LU decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU.html)
  is used to solve the linear system(s).
  The strategy used in v1.1.3 suffered from a reduced numerical accuracy, due to the fact that
  the polarization weights were not correctly defined.

### Removed

- The `hermitivitize` function will only work correctly on matrices. This
  reverts modifications in the previous release.

## [Version 1.1.3] - 2016-07-03

### Changed

- The `PEDRA.OUT` cavity generator log now reports the initial _and_ final
  lists of spheres. The final list only contains those spheres that were
  actually tesselated and, possibly, the added spheres.
- For all solvers, the symmetrization needed to obtain the polarization weights
  happens directly on the computed charges, instead of symmetrizing the system
  matrices.
- The `IEFSolver` object stores the unsymmetrized T^-1R matrices.
  A [partially pivoted LU decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU.html)
  is used to compute T^-1R.
  A [robust Cholesky decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html) is
  used to form the R matrix in the anisotropic IEF case.
- The `CPCMSolver` object now stores the scaled, unsymmetrized S matrix. The
  explicit calculation and storage of its inverse is thus avoided.
  A [robust Cholesky decomposition](http://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html) is
  used to solve the linear equation system.
- The `hermitivitize` function can now correctly symmetrize vectors.

### Fixed

- A fix for the initialization of the explicit list of spheres when running the
  standalone executable. The bug prevented the generation of the correct
  `Molecule` object, with subsequent failure in the cavity generator.
- A memory leak occuring in the cavity generator PEDRA was fixed. This was uncovered by @shoefener
  and manifested only with considerably large cavities (> 200 input spheres)

### Removed

- The function `CPCMMatrix` in the `SolverImpl.hpp` header file is no longer available.

## [Version 1.1.2] - 2016-05-31

### Fixed

- Signatures for strings in Fortran90 bindings. They have now the proper
  C interoperable type `character(kind=c_char, len=1) :: label(lbl_size)`.
  For the host this means that surface function labels will have to be declared
  as character arrays, for example: `character :: label(7) = (/'T', 'o', 't', 'M', 'E', 'P', char(0)/)`

### Changed

- More informative error messages for runtime crashes caused by access to
  surface functions.
- The signatures for the interface functions now accept and/or return `int` (`c_int`)
  instead of `size_t` (`c_size_t`). This simplifies interfacing with Fortran hosts.

## [Version 1.1.1] - 2016-03-10

### Added

- A runtime check to ensure that all atoms have a nonzero radius. API kills
  program execution if this is the case.
- An API function to retrieve the areas/weights of the cavity finite elements.
  The values in the returned array are in Bohr^2. Addresses a feature request
  from @shoefener (Issue #13)
- The standalone executable `run_pcm` is now tested within the unit tests
  suite. The tests cover the cases where the cavity is given implicitly,
  explicitly or by substitution of radii on chosen atoms.

### Changed

- Boundary integral operators classes learnt to accept a scaling factor for the
  diagonal elements of the approximate collocation matrices. The change is
  reflected in the Green's funtion classes and in the input parsing. Addresses
  a feature request from @shoefener (Issue #16)
- `GePolCavity` learnt to print also the list of spheres used to generate the
  cavity.
- Different internal handling of conversion factors from Bohr to Angstrom.
- CMake minimum required version is 2.8.10
- `Atom`, `Solvent` and `Sphere` are now PODs. The radii and solvent lists are free
  functions.
- `PCMSOLVER_ERROR` kills program execution when an error arises but does not
  use C++ exceptions.
- `include`-s are now specified on a per-directory basis (see programmers'
  manual for a more detailed explanation)
- Default types for template paramters `DerivativeTraits`,  `IntegratorPolicy`
  and `ProfilePolicy` are now given for the Green's functions classes. This
  reduced the verbosity in instatiating these objects significantly.

### Known Issues

- The new printer in `GePolCavity` might not work properly when an explicit list
  of spheres is provided in the input.
- On Ubuntu 12.10, 32 bit the Intel compiler version 2013.1 produces a faulty
  library. It is possibly a bug in the implementation of `iso_c_binding`, see
  Issue #25

### Removed

- `SurfaceFunction` as a class is no longer available. We keep track of surface
  functions at the interface level _via_ a label-vector map.

## [Version 1.1.0] - 2016-02-07

### Added

- Green's function for diffuse interfaces in spherical symmetry

### Changed

- CMake minimum required version is 2.8.8 (2016-01-08)
- Documentation is now served [here](http://pcmsolver.readthedocs.org/)

## v1.0.4 - 2015-07-22 [YANKED]

## v1.0.3 - 2015-03-29 [YANKED]

## v1.0.2 - 2015-03-28 [YANKED]

## v1.0.1 - 2015-01-06 [YANKED]

## v1.0.0 - 2014-09-30 [YANKED]

[Unreleased]: https://github.com/PCMSolver/pcmsolver/compare/v1.3.0...HEAD
[Version 1.3.0]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.3...v1.3.0
[Version 1.2.3]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.2...v1.2.3
[Version 1.2.2]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.1...v1.2.2
[Version 1.2.1]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.0...v1.2.1
[Version 1.2.0]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.0-rc1...v1.2.0
[Version 1.2.0-rc1]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.12...v1.2.0-rc1
[Version 1.1.12]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.11...v1.1.12
[Version 1.1.11]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.10...v1.1.11
[Version 1.1.10]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.9...v1.1.10
[Version 1.1.9]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.8...v1.1.9
[Version 1.1.8]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.7...v1.1.8
[Version 1.1.7]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.6...v1.1.7
[Version 1.1.6]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.5...v1.1.6
[Version 1.1.5]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.4...v1.1.5
[Version 1.1.4]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.3...v1.1.4
[Version 1.1.3]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.2...v1.1.3
[Version 1.1.2]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.1...v1.1.3
[Version 1.1.1]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.0...v1.1.1
[Version 1.1.0]: https://github.com/PCMSolver/pcmsolver/releases/tag/v1.1.0

[GitHub]: https://github.com/PCMSolver/pcmsolver
[GitLab]: https://gitlab.com/PCMSolver/pcmsolver
<!--- Provide a general summary of your changes in the Title above -->
<!--- Use labels to help the developers -->
<!--- Remove the unnecessary sections when filing -->

## Description
<!--- Describe your changes in detail -->

## Motivation and Context
<!--- Why is this change required? What problem does it solve? -->
<!--- If it fixes an open issue, please link to the issue here. -->

## How Has This Been Tested?
<!--- Please describe in detail how you tested your changes. -->
<!--- Include details of your testing environment, and the tests you ran to -->
<!--- see how your change affects other areas of the code, etc. -->

## Screenshots (if appropriate):

## Todos
<!--- Notable points that this PR has either accomplished or will accomplish. -->
* **Developer Interest**
<!--- Changes affecting developers -->
  - [ ] Feature1
* **User-Facing for Release Notes**
<!--- Changes affecting users -->
  - [ ] Feature2

## Types of changes
<!--- What types of changes does your code introduce? Put an `x` in all the boxes that apply: -->
- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to change)

## Questions
<!--- Questions to the developers -->
- [ ]  Question1

## Status
<!--- Check this box when ready to be merged -->
- [ ] Ready to go
- [ ] Cherry-pick to latest release branch
<!--- Provide a general summary of the issue in the Title above -->
<!--- Use labels to help the developers -->
<!--- Remove the unnecessary sections when filing -->

## Expected Behavior
<!--- If you're describing a bug, tell us what should happen -->
<!--- If you're suggesting a change/improvement, tell us how it should work -->

## Current Behavior
<!--- If describing a bug, tell us what happens instead of the expected behavior -->
<!--- If suggesting a change/improvement, explain the difference from current behavior -->

## Possible Solution
<!--- Not obligatory, but suggest a fix/reason for the bug, -->
<!--- or ideas how to implement the addition or change -->

## Steps to Reproduce (for bugs)
<!--- Provide a link to a live example, or an unambiguous set of steps to -->
<!--- reproduce this bug. Include code to reproduce, if relevant -->
1.
2.
3.
4.

## Context
<!--- How has this issue affected you? What are you trying to accomplish? -->
<!--- Providing context helps us come up with a solution that is most useful in the real world -->

## Your Environment
<!--- Include as many relevant details about the environment you experienced the bug in -->
* Version used:
* Environment name and version (e.g. PHP 5.4 on nginx 1.9.1):
* Server type and version:
* Operating System and version:
* Link to your project:
The CMake infrastructure for this project is generated using [Autocmake]
by Radovan Bast, Roberto Di Remigio, Jonas Juselius and contributors.
The `update.py` Python script and the contents of the directories `autocmake` and `downloaded` are licensed
under the terms of the [BSD-3-Clause license], unless otherwise stated.

[Autocmake]: http://autocmake.org
[BSD-3-Clause license]: https://tldrlegal.com/license/bsd-3-clause-license-(revised)# Namespaces


We use namespaces to delimit the visibility of functions and classes defined in
the various subdirectories of the project.
Namespaces provide a convenient layered structure to the project and we use
them as a convention to signal which functions and classes are supposed to be
used in any given layer.
The top-level namespace is called `pcm` and includes all functions and classes
that can be called from the outside world, i.e. a C++ API.
Each subdirectory introduces a new namespace of the same name, nested into `pcm`.
Code that can be used _outside_ of a given subdirectory is put directly in the
`pcm` namespace, i.e. the outermost layer.
Finally, the namespace `detail`, at the third level of nesting, is used for
functions and classes that are used exclusively within the code in a given
subdirectory.
Publications
============

Peer-reviewed journal articles
------------------------------

2015
~~~~

+ `Four-Component Relativistic Calculations in Solution with the Polarizable Continuum Model of Solvation: Theory, Implementation, and Application to the Group 16 Dihydrides H2X (X = O, S, Se, Te, Po) <http://pubs.acs.org/doi/abs/10.1021/jp507279y>`_
+ `Wavelet Formulation of the Polarizable Continuum Model. II. Use of Piecewise Bilinear Boundary Elements <http://pubs.rsc.org/en/content/articlelanding/2015/cp/c5cp03410h>`_

2016
~~~~

+ `A Polarizable Continuum Model for Molecules at Spherical Diffuse Interfaces <http://dx.doi.org/10.1063/1.4943782>`_

2017
~~~~

+ `Four-Component Relativistic Density Functional Theory with the Polarizable Continuum Model: Application to EPR Parameters and Paramagnetic NMR Shifts <http://dx.doi.org/10.1080/00268976.2016.1239846>`_
+ `Open-ended formulation of self-consistent field response theory with the polarizable continuum model for solvation <https://doi.org/10.1039/C6CP06814F>`_
+ `Psi4 1.1: An Open-Source Electronic Structure Program Emphasizing Automation, Advanced Libraries, and Interoperability <https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00174>`_
+ `Combining frozen-density embedding with the conductor-like screening model using Lagrangian techniques for response properties <http://onlinelibrary.wiley.com/doi/10.1002/jcc.24813/abstract>`_

Theses
------

+ `The Polarizable Continuum Model Goes Viral! Extensible, Modular and Sustainable Development of Quantum Mechanical Continuum Solvation Models <https://munin.uit.no/handle/10037/10786>`_ Doctoral thesis, Roberto Di Remigio, January 2017.

Presentations
-------------

+ `A modular implementation of the Polarizable Continuum Model for Solvation <https://www.dropbox.com/s/uzzv8c0wx8eswbc/talk_pisa.pdf?dl=0>`_ Presentation given by Roberto Di Remigio at the workshop in honour of professor Jacopo Tomasi's 80th birthday. Pisa, August 31 - September 1 2014.
+ `The Polarizable Continuum Model Goes Viral! <http://tinyurl.com/phd-forsvaring>`_ PhD defense, Roberto Di Remigio, January 16 2017.
+ PCMSolver: a modern, modular approach to include solvation in any quantum chemistry code. Presentation given by Luca Frediani at WATOC 2017. Munich, August 27 - September 1 2017.

Posters
-------

+ `Plug the solvent in your favorite QM program <https://www.dropbox.com/s/gmj6l54mdj6r9z7/posterICQC.pdf?dl=0>`_ Presented by Luca Frediani at the 14th International Congress of Quantum Chemistry. Boulder, Colorado, June 25-30 2012.
+ `4-Component Relativistic Calculations in Solution with the Polarizable Continuum Model of Solvation <https://www.dropbox.com/s/edvrimiwh5rlg9y/posterFemEx.pdf?dl=0>`_ Presented by Roberto Di Remigio at the FemEx-Oslo conference. Oslo, June 13-16 2014.


References
==========

.. bibliography:: pcmsolver.bib
.. PCMSolver documentation master file, created by
   sphinx-quickstart on Mon Oct 26 15:18:26 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PCMSolver's documentation!
=====================================

This is the documentation for the PCMSolver application programming interface.
PCMSolver is an API for solving the Polarizable Continuum Model electrostatic problem :cite:`Tomasi2005`

.. image:: gfx/pcmsolver-scheme.png
   :scale: 70 %
   :align: center

With PCMSolver we aim to:

1. provide a plug-and-play library for adding the PCM functionality to *any* quantum chemistry program;

2. create a playground for easily extending the implementation of the model.

PCMSolver is distributed under the terms of the GNU Lesser General Public License.
An archive with the currently released source can be found on `GitHub <https://github.com/PCMSolver/pcmsolver/releases>`_.

.. literalinclude:: snippets/citation.bib
   :language: tex

PCMSolver has been added to the following quantum chemistry programs

+ `Psi4 <http://www.psicode.org/>`_
+ `DALTON <http://daltonprogram.org/>`_
+ `LSDALTON <http://daltonprogram.org/>`_
+ `DIRAC <http://www.diracprogram.org/>`_
+ `ReSpect <http://www.respectprogram.org/>`_
+ `KOALA <https://dx.doi.org/10.1002/jcc.23679>`_

Don't see you code listed here? `Please contact us <roberto.d.remigio@uit.no>`_.

.. toctree::
   :caption: Table of Contents
   :name: mastertoc
   :hidden:

   users/users-manual
   publications
   programmers/programmers-manual
   code-reference/classes-and-functions
   zreferences

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Documentation
=============

This documentation is generated using `Sphinx <http://sphinx-doc.org/>`_ and
`Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_ The two softwares are
bridged by means of the `Breathe extension <https://breathe.readthedocs.org/>`_
The online version of this documentation is built and served by `Read The Docs
<https://readthedocs.org/>`_.  The webpage http://pcmsolver.readthedocs.org/ is
updated on each push to the public GitHub repository.


How and what to document
------------------------

Doxygen enables documenting the code in the source code files thus removing a
"barrier" for developers.  To avoid that the code degenerates into a Big Ball
of Mud, it is mandatory to document directly within the source code classes and
functions. To document general programming principles, design choices,
maintenance etc. you can create a .rst file in the ``doc`` directory. Remember
to refer the new file inside the ``index.rst`` file (it won't be parsed
otherwise).  Sphing uses `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ and `Markdown
<https://daringfireball.net/projects/markdown/>`_. Support for Markdown is not
as extensive as for reStructuredText, see `these comments
<https://blog.readthedocs.com/adding-markdown-support/>`_. Follow the guidelines
in :cite:`Wilson2014` regarding what to document.

Write the documentation in the header file. To document a class, put
``/*! \class <myclass>`` inside the namespace but before the class.
Add the following to a ``.rst`` file:

.. code-block:: rst

  .. doxygenclass:: <namespace>::<myclass>
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Do similar when documenting ``struct``-s and complete files.

.. note::

   Use ``/*! */`` to open and close a Doxygen comment.

Documenting methods in derived classes
--------------------------------------

Virtual methods should only be documented in the base classes.
This avoids unnecessary verbosity and conforms to the principle: "Document
_what_, not _how_" :cite:`Wilson2014`
If you feel the _how_ needs to be explicitly documented, add some notes in the
appropriate ``.rst`` file.

How does this work?
-------------------

To have an offline version of the documentation just issue
in the ``doc`` folder:

.. code-block:: bash

   sphinx-build . _build

The HTML will be stored in ``_build/``. Open the ``_build/index.html`` file with
your browser to see and browse the documentation.

.. warning::

   It is only possible to build documentation locally from within the ``doc``
   folder.
   This choice was made to simplify the set up of the ReadTheDocs and local
   documentation building procedures and to minimize the chances of breaking
   either.

.. note::

   Building the documentation requires Python, Doxygen, Sphinx, Perl and the
   Python modules breathe, matplotlib, sphinx-rtd-theme, sphinxcontrib-bibtex
   and recommonmark.
   The required python modules can be installed by running ``pip install -r
   requirements.txt``.
   There is also a ``Pipfile`` in case people prefer to use ``pipenv``.
CMake usage
===========

This is a brief guide to our CMake infrastructure which is managed
*via* `Autocmake <http://autocmake.readthedocs.org/en/latest/>`_

.. warning::

   The minimum required CMake version is 2.8.10

Adding new source subdirectories and/or files
---------------------------------------------

Developers **HAVE TO** manually list the sources in a given subdirectory
of the main source directory ``src/``. In our previous infrastructure this
was not necessary, but the developers needed to trigger CMake to regenerate the
Makefiles manually.

New subdirectory
................

First of all, you will have to let CMake know that a new source-containing
subdirectory has been added to the source tree. Due to the hierarchical
approach CMake is based upon you will need to modify the ``CMakeLists.txt`` in
the ``src`` directory and create a new one in your new subdirectory.  For the
first step:

   1. if your new subdirectory contains header files, add a line like
   the following to the ``CMakeLists.txt`` file contained in the ``src`` directory:

   .. code-block:: bash

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
Code contributions
==================

We have adopted a fully public *fork and pull request* workflow, where every proposed changeset
has to go through a code review and approval process.

The code changes are developed on a *branch* of the
*fork*. When completed, the developer submits the changes for review
through the web interface: a *pull request* (PR) is opened, requesting that the changes from
the *source branch* on the fork be merged into a *target branch* in
the canonical repository.
Once the PR is open, the new code is automatically tested.
Core developers of PCMSolver will then review the contribution and
discuss additional changes to be made.
Eventually, if all the tests are passing and a developer approves the suggested
contribution, the changes are merged into the target branch.
The target branch is (usually) the `master` branch, that is, the main development branch.

.. note::

  All PRs goes to the master branch

The creator of the PR is responsible for keeping the code up to date with master, 
so the code in the PR reflects what will be the code in the master branch after merging.

Branching Model
---------------

We are using the `stable mainline branching model for Git <http://www.bitsnbites.eu/a-stable-mainline-branching-model-for-git/>`_.
In the main repository on github there are two types of branches:

- *one* main developing branch, called ``master``
- release branches

A new release branch is created from the master branch for a new release, with the format ``release/vMAJOR.MINOR``.
A release branch will never be merged back to the master branch and will only receive bug fixes, thus no new features.
These bug fixes would be cherry picked from the master branch, to ensure that the master branch always contains all bug fixes.
In case a bug fix is only relevant for a given release, the bug should be fixed with a PR directly to the corresponding release branch.
In case a bug fix is easy to perform on a release branch but challenging to perform on the master branch, the fix can be directed to a release branch.
Then an issue *have* to be created to make sure it will also be fixed on the master branch.

Feature branches are not created on the main repository, but on forks. 
These are based on the master branch from the main repository and merged into the master branch through pull requests.

Changelog
=========

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
===========================

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
====================

`Git pre-commit hooks <https://git-scm.com/book/gr/v2/Customizing-Git-Git-Hooks>`_ are used to
keep track of code style and license header in source files.
Code style is checked using ``clang-format`` for C/C++ and ``yapf`` for Python.

.. warning::
   **You need to install ``clang-format`` (v3.9 recommended) and ``yapf``
   (v0.20 recommended) to run the code style validation hook!**

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

Do not include header files unnecessarily. Even if PCMSolver is not a big
project, unnecessary include directives and/or forward declarations introduce
nasty interdependencies among different parts of the code.  This reflects
mainly in longer compilation times, but also in uglier looking code (see also
the discussion in :cite:`Sutter1999`).

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

We perform unit testing of our API. The unit testing framework used is
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
   to the ``CMakeLists.txt``

   .. code-block:: cmake

      add_subdirectory(new_subdir)

#. create a ``CMakeLists.txt`` inside your new subdirectory.
   This ``CMakeLists.txt`` adds the source for a given unit test to the global ``UnitTestsSources``
   property and notifies CTest that a test with given name is part of the test suite.
   The generation of the ``CMakeLists.txt`` can be managed by ``make_cmake_files.py`` Python script.
   This will take care of also setting up CTest labels. This helps in further grouping
   the tests for our convenience.
   Catch uses tags to index tests and tags are surrounded by square brackets. The Python script
   inspects the sources and extracts labels from Catch tags.
   The ``add_Catch_test`` CMake macro takes care of the rest:

   .. code-block:: cmake

      add_Catch_test(
        NAME
          <test-name> # Mandatory!
        LABELS
          <test-labels> # Mandatory! One per line, for readability
        DEPENDS
          <test-dependencies> # Optional. One per line, for readability
        REFERENCE_FILES
          <test-refs> # Optional. One per line, for readability
        COST
          <test-cost> # Optional. Roughly the seconds it takes to run the test
      )

   We require that each source file containing tests follows the naming convention
   new_subdir_testname and that testname gives some clue to what is being tested.
   Depending on the execution of tests in a different subdirectory is bad practice.
   A possible workaround is to add some kind of input file and create a text fixture
   that sets up the test environment. Have a look in the ``tests/input`` directory
   for an example
#. create the ``.cpp`` files containing the tests. Use the following template:

   .. literalinclude:: ../snippets/test_example.cpp
      :language: cpp
      :linenos:

   In this example we are creating a test fixture. The fixture will instatiate
   a ``GePolCavity`` with fixed parameters. The result is then tested against reference values
   in the various ``SECTION`` s.
   It is **important** to add the documentation lines on top of the tests, to help other
   developers understand which class is being tested and what parameters are being tested.
   Within Catch fixtures are created behind the curtains, you do not need to worry about
   those details. This results in somewhat terser test source files.
Versioning and minting a new release
====================================

Our versioning machinery is based on a modified version of the ``versioner.py``
script devised by Lori A. Burns (Georgia Tech) for the `Psi4
<http://www.psicode.org>`_ quantum chemistry code.
The documentation that follows is also adapted from the corresponding Psi4
documentation, available at `this link <http://www.psicode.org/psi4manual/1.1/manage_git.html>`_

This guide will walk you through the actions to perform to mint a new release
of the code. Version numbering follows the guidelines of `semantic versioning
<http://semver.org/>`_. The allowed format is ``MAJOR.MINOR.PATCH-DESCRIBE``,
where ``DESCRIBE`` can be a string describing a prerelease state, such as
``rc2``, ``alpha1``, ``beta3`` and so forth.

Minting a new release
---------------------

The ``tools/metadata.py`` file records the versioning information for the current
release. The information in this file is used by the ``versioner.py`` script to
compute a *unique version number* for development snapshots.

.. note::

   To correctly mint a new release, you will have to be on the latest release
   branch of (i) a direct clone or (ii) clone-of-fork with release branch
   up-to-date with upstream (including tags!!!) and with upstream as remote.

This is the step-by-step guide to releasing a new version of PCMSolver:

#. **DECIDE** an upcoming version number, say ``1.2.0``.
#. **TIDY UP** ``CHANGELOG.md``:

   * **SET** the topmost header to the upcoming version number and release date.

     ::

       ## [Version 1.2.0] - 2018-03-31

   * **CHECK** that the links at the bottom of the document are correct.

     ::

       [Unreleased]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.0...HEAD
       [Version 1.2.0]: https://github.com/PCMSolver/pcmsolver/compare/v1.2.0-rc1...v1.2.0
       [Version 1.2.0-rc1]: https://github.com/PCMSolver/pcmsolver/compare/v1.1.12...v1.2.0-rc1

#. **UPDATE** the ``AUTHORS.md`` file:

   * Run ``git shortlog -sn`` and cross-check with the current contents of ``AUTHORS.md``.
     Edit where necessary and don't forget to include, where
     available, the GitHub handle. Authors are ordered by the number of commits.
   * Update the revision date at the bottom of this file.

     ::

       >>> cat AUTHORS.md
       ## Individual Contributors

       - Roberto Di Remigio (@robertodr)
       - Luca Frediani (@ilfreddy)
       - Monica Bugeanu (@mbugeanu)
       - Arnfinn Hykkerud Steindal (@arnfinn)
       - Radovan Bast (@bast)
       - T. Daniel Crawford (@lothian)
       - Krzysztof Mozgawa
       - Lori A. Burns (@loriab)
       - Ville Weijo (@vweijo)
       - Ward Poelmans (@wpoely86)

       This list was obtained 2018-03-02 by running `git shortlog -sn`

#. **CHECK** that the ``.mailmap`` file is up-to-date.
#. **CHECK** that the documentation builds locally.
#. **ACT** to check all the changed files in.
#. **OBSERVE** current versioning state

   * https://github.com/PCMSolver/pcmsolver/releases says ``v1.2.0-rc1`` & ``9a8c391``

    ::

      >>> git tag
      v1.1.0
      v1.1.1
      v1.1.10
      v1.1.11
      v1.1.12
      v1.1.2
      v1.1.3
      v1.1.4
      v1.1.5
      v1.1.6
      v1.1.7
      v1.1.8
      v1.1.9
      v1.2.0-rc1

      >>> cat tools/metadata.py
      __version__ = '1.2.0-rc1'
      __version_long = '1.2.0-rc1+9a8c391'
      __version_upcoming_annotated_v_tag = '1.2.0'
      __version_most_recent_release = '1.1.12'


      def version_formatter(dummy):
          return '(inplace)'

      >>> git describe --abbrev=7 --long --always HEAD
      v1.2.0-rc1-14-gfc02d9d

      >>> git describe --abbrev=7 --long --dirty
      v1.2.0-rc1-14-gfc02d9d-dirty

      >>> python tools/versioner.py
      Defining development snapshot version: 1.2.0.dev14+fc02d9d (computed)
      1.2.0.dev14 {versioning-script} fc02d9d 1.1.12.999 dirty  1.1.12 <-- 1.2.0.dev14+fc02d9d

      >>> git diff

   * Observe that current latest tag matches metadata script and git
     describe, that GH releases matches metadata script, that upcoming in
     metadata script matches current ``versioner.py`` version.

#. **ACT** to bump tag in code. The current tag is ``v1.2.0-rc1``, the imminent tag is ``v1.2.0``.

   * Edit current & prospective tag in ``tools/metadata.py``. Use your
     decided-upon tag ``v1.2.0`` and a speculative next tag, say ``v1.3.0``,
     and use 7 "z"s for the part you can't predict.

     ::

       >>> vim tools/metadata.py

       >>> git diff
       diff --git a/tools/metadata.py b/tools/metadata.py
       index 5d87b55..6cbc05e 100644
       --- a/tools/metadata.py
       +++ b/tools/metadata.py
       @@ -1,6 +1,6 @@
       -__version__ = '1.2.0-rc1'
       -__version_long = '1.2.0-rc1+9a8c391'
       -__version_upcoming_annotated_v_tag = '1.2.0'
       -__version_most_recent_release = '1.1.12'
       +__version__ = '1.2.0'
       +__version_long = '1.2.0+zzzzzzz'
       +__version_upcoming_annotated_v_tag = '1.3.0'
       +__version_most_recent_release = '1.2.0'

   * **COMMIT** changes to ``tools/metadata.py``.

     ::

       >>> git add tools/metadata.py
       >>> git commit -m "Bump version to v1.2.0"

#. **OBSERVE** undefined version state. Note the 7-character git hash for the new commit, here ``fc02d9d``.

   ::

     >>> git describe --abbrev=7 --long --always HEAD
     v1.2.0-rc1-14-gfc02d9d

     >>> git describe --abbrev=7 --long --dirty
     v1.2.0-rc1-14-gfc02d9d-dirty

     >>> python tools/versioner.py
     Undefining version for irreconcilable tags: 1.2.0-rc1 (computed) vs 1.2.0 (recorded)
     undefined {versioning-script} fc02d9d 1.2.0.999 dirty  1.2 <-- undefined+fc02d9d

#. **ACT** to bump tag in git, then bump git tag in code.

   * Use the decided-upon tag ``v1.2.0`` and the observed hash ``fc02d9d`` to
     mint a new *annotated* tag, minding that "v"s are present here.

   * Use the observed hash to edit ``tools/metadata.py`` and commit immediately.

   ::

     >>> git tag -a v1.2.0 fc02d9d -m "Version 1.2.0 released"

     >>> vim tools/metadata.py

     >>> git diff
     diff --git a/tools/metadata.py b/tools/metadata.py
     index 6cbc05e..fdc202e 100644
     --- a/tools/metadata.py
     +++ b/tools/metadata.py
     @@ -1,5 +1,5 @@
      __version__ = '1.2.0'
     -__version_long = '1.2.0+zzzzzzz'
     +__version_long = '1.2.0+fc02d9d'
      __version_upcoming_annotated_v_tag = '1.3.0'
      __version_most_recent_release = '1.2.0'

     >>> python tools/versioner.py
     Amazing, this can't actually happen that git hash stored at git commit.

     >>> git add tools/metadata.py

     >>> git commit -m "Records tag for v1.2.0"

#. **OBSERVE** current versioning state. There is nothing to take note of. This
   is just a snapshot to ensure that you did not mess up.

    ::

      >>> python tools/versioner.py
      Defining development snapshot version: 1.2.0.dev1+4e0596e (computed)
      1.2.0.dev1 {master} 4e0596e 1.2.0.999   1.2 <-- 1.2.0.dev1+4e0596e

      >>> git describe --abbrev=7 --long --always HEAD
      v1.2.0-1-g4e0596e

      >>> git describe --abbrev=7 --long --dirty
      v1.2.0-1-g4e0596e

      >>> git tag
      v1.1.0
      v1.1.1
      v1.1.10
      v1.1.11
      v1.1.12
      v1.1.2
      v1.1.3
      v1.1.4
      v1.1.5
      v1.1.6
      v1.1.7
      v1.1.8
      v1.1.9
      v1.2.0-rc1
      v1.2.0

      >>> cat tools/metadata.py
      __version__ = '1.2.0'
      __version_long = '1.2.0+fc02d9d'
      __version_upcoming_annotated_v_tag = '1.3.0'
      __version_most_recent_release = '1.2.0'

      >>> cat metadata.out.py | head -8
      __version__ = '1.2.0.dev1'
      __version_branch_name = 'master'
      __version_cmake = '1.2.0.999'
      __version_is_clean = 'True'
      __version_last_release = '1.2.0'
      __version_long = '1.2.0.dev1+4e0596e'
      __version_prerelease = 'False'
      __version_release = 'False'

      >>> git log --oneline
      4e0596e Records tag for v1.2.0
      fc02d9d Bump version to v1.2.0

#. **ACT** to inform remote of bump

   * Temporarily disengage "Include administrators" on protected release branch.

    ::

      >>> git push origin release/1.2

      >>> git push origin v1.2.0

   * Now https://github.com/PCMSolver/pcmsolver/releases says ``v1.2.0`` & ``fc023d9d``

#. **EDIT** release description in the `GitHub web UI <https://github.com/PCMSolver/pcmsolver/releases>`_.

`Zenodo <https://zenodo.org/>`_ will automatically generate a new, versioned
DOI for the new release. It is no longer necessary to update the badge
in the ``README.md`` since it will always resolve to the latest released by
Zenodo.


How to create and remove an annotated Git tag on a remote
---------------------------------------------------------

PCMSolver versioning only works with *annotated* tags, not *lightweight*
tags as are created with the `GitHub interface
<https://github.com/PCMSolver/pcmsolver/releases/new>`_

* Create *annotated* tag::

    >>> git tag -a v1.1.12 <git hash if not current> -m "Version 1.1.12 released"
    >>> git push upstream --tags

* Delete tag::

    >>> git tag -d v1.1.12
    >>> git push origin :refs/tags/v1.1.12

* Pull tags::

    >>> git fetch <remote> 'refs/tags/*:refs/tags/*'
Timer class
-----------

The ``Timer`` class enables timing of execution throughout the module.
Timer support is enabled by passing ``-DENABLE_TIMER=ON`` to the ``setup.py``
script.
Timing macros are available by inclusion of the ``Config.hpp`` header file.

The class is basically a wrapper around an ordered map of strings and cpu timers.
To time a code snippet:

.. code-block:: cpp

   TIMER_ON("code-snippet");
   // code-snippet
   TIMER_OFF("code-snippet");

The timings are printed out to the ``pcmsolver.timer.dat`` by a call
to the ``TIMER_DONE`` macro. This should obviously happen at the very end
of the execution!

.. doxygenfile:: TimerInterface.hpp
   :project: PCMSolver
=============================
PCMSolver Programmers' Manual
=============================

.. toctree::

   general-structure
   coding-standards
   documentation
   cmake-usage
   versioning
   maintenance
   profiling
   testing
   timer-class
General Structure
=================

.. image:: structure.png
   :scale: 90 %
   :align: center

External libraries:

+ parts of the C++ `Boost <http://www.boost.org/>`_ libraries are used to provide
  various functionality, like ordinary differential equations integrators.
  The source for the 1.54.0 release is shipped with the
  module's source code. Some of the libraries used
  need to be compiled. Boost is released under the terms
  of the `Boost Software License, v1.0 <http://opensource.org/licenses/BSL-1.0>`_ (see also
  http://www.boost.org/users/license.html)

  .. warning::

     As of v1.1.11 we have started removing the dependency from Boost.
     The use of Boost is thus deprecated.

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
+ the `libtaylor <https://github.com/uekstrom/libtaylor>`_ library implementing automatic differentiation and available
  under the terms of the `MIT License <(http://opensource.org/licenses/MIT>`_.

Third-party code snippets:

+ Fortran subroutines `dsyevv3`, `dsyevh3`, `dsyevj3` for the diagonalization
  of 3x3 Hermitian matrices.  These subroutines were copied verbatim from the
  source code provided by `Joachim Kopp <http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/>`_
  and described in :cite:`KOPP2008` (also available on the `arXiv <http://arxiv.org/abs/physics/0610206>`_) The diagonalization
  subroutines are made available under the terms of the `GNU Lesser General
  Public License, v2.1 <http://opensource.org/licenses/LGPL-2.1>`_
+ C++ cnpy library for saving arrays in C++ into Numpy arrays. The library is
  from `Carl Rogers <https://github.com/rogersce/cnpy>`_ under the terms of the
  `MIT License <(http://opensource.org/licenses/MIT>`_.
  The version in PCMSolver is slightly different.
Profiling
---------

You should obtain profiling information before attempting any optimization of
the code. There are many ways of obtaining this information, but we have only
experimented with the following:

#. Using Linux ``perf`` and related `tools <http://www.brendangregg.com/perf.html>`_.
#. Using ``gperftools``.
#. Using Intel VTune.

Profiling should be done using the standalone executable ``run_pcm`` and any of
the input files gathered under the ``tests/benchmark`` directory. These files
are copied to the build directory. If you are lazy, you can run the profiling
from the build directory:

.. code-block:: bash

   >>> cd tests/benchmark

   >>> env PYTHONPATH=<build_dir>/lib64/python:$PYTHONPATH
          python <build_dir>/bin/go_pcm.py --inp=standalone.pcm --exe=<build_dir>/bin

Using ``perf``
==============

``perf`` is a tool available on Linux. Though part of the kernel tools, it is
not usually preinstalled on most Linux distributions. For visualization
purposes we also need `additional tools <https://github.com/brendangregg/perf-tools>`_,
in particular the `flame graph generation scripts <https://github.com/brendangregg/FlameGraph>`_
Probably your distribution has them prepackaged already.
``perf`` will trace all CPU events on your system, hence you might need to
fiddle with some kernel set up files to get permissions to trace events.

.. note::
   ``perf`` **is NOT** available on ``stallo``. Even if it were, you would
   probably not have permissions to record kernel traces.

These are the instructions I used:

1. Trace execution. This will save CPU stack traces to a ``perf.data`` file.
   Successive runs do not overwrite this file.

   .. code-block:: bash

      >>> cd tests/benchmark

      >>> perf record -F 99 -g -- env PYTHONPATH=<build_dir>/lib64/python:$PYTHONPATH python
                    <build_dir>/bin/go_pcm.py --inp=standalone.pcm --exe=<build_dir>/bin

2. Get reports. There are different ways of getting a report from the
   ``perf.data`` file. The following will generate a call tree.

   .. code-block:: bash

      >>> perf report --stdio

3. Generate an interactive flame graph.

   .. code-block:: bash

      >>> perf script | stackcollapse-perf.pl > out.perf-folded

      >>> cat out.perf-folded | flamegraph.pl > perf-run_pcm.svg

Using ``gperftools``
====================

This set of tools was previously known as Google Performance Tools. The
executable needs to be linked against the ``profiler``, ``tcmalloc``
and ``unwind`` libraries.
CMake will attempt to find them. If this fails, you will have to install them,
you should either check if they are available for your distribution or compile
from source.
In principle, one could use the ``LD_PRELOAD`` mechanism to skip the *ad hoc*
compilation of the executable.

.. note::
   ``gperftools`` **is** available on ``stallo``, but it's an ancient version.

1. Configure the code with the ``--gperf`` option enabled. CPU and heap
   profiling, together with heap-checking will be available.

2. CPU profiling can be done with the following command:

   .. code-block:: bash

      >>> env CPUPROFILE=run_pcm.cpu.prof PYTHONPATH=<build_dir>/lib64/python:$PYTHONPATH
              python <build_dir>/bin/go_pcm.py --inp=standalone.pcm --exe=<build_dir>/bin

  This will save the data to the ``run_pcm.cpu.prof`` file. To analyze the gathered
  data we can use the ``pprof`` script:

  .. code-block:: bash

     >>> pprof --text <build_dir>/bin/run_pcm run_pcm.cpu.prof

  This will print a table. Any row will look like the following:

  .. code-block:: bash

     2228   7.2%  24.8%    28872  93.4% pcm::utils::splineInterpolation

  where the columns respectively report:

  #. Number of profiling samples in this function.
  #. Percentage of profiling samples in this function.
  #. Percentage of profiling samples in the functions printed so far.
  #. Number of profiling samples in this function and its callees.
  #. Percentage of profiling samples in this function and its callees.
  #. Function name.

  For more details look `here <https://gperftools.github.io/gperftools/cpuprofile.html>`_

3. Heap profiling can be done with the following command:

   .. code-block:: bash

      >>> env HEAPPROFILE=run_pcm.hprof PYTHONPATH=<build_dir>/lib64/python:$PYTHONPATH
              python <build_dir>/bin/go_pcm.py --inp=standalone.pcm --exe=<build_dir>/bin

  This will output a series of datafiles ``run_pcm.hprof.0000.heap``,
  ``run_pcm.hprof.0001.heap`` and so forth. You will have to kill execution
  when enough samples have been collected.
  Analysis of the heap profiling data can be done using ``pprof``. `Read more
  here <https://gperftools.github.io/gperftools/heapprofile.html>`_


Using Intel VTune
=================

This is probably the easiest way to profile the code.
`VTune <https://software.intel.com/en-us/intel-vtune-amplifier-xe>`_ is Intel software, it might be possible to get a personal, free license.
The instructions will hold on any machine where VTune is installed and you can
look for more details on the `online documentation <https://software.intel.com/en-us/vtune-amplifier-help>`_
You can, in principle, use the GUI. I haven't managed to do that though.

On ``stallo``, start an interactive job and load the following modules:

.. code-block:: bash

   >>> module load intel/2018a

   >>> module load CMake

   >>> module load VTune

   >>> export BOOST_INCLUDEDIR=/home/roberto/Software/boost/include

   >>> export BOOST_LIBRARYDIR=/home/roberto/Software/boost/lib

You will need to compile with optimizations activated, *i.e.* release mode.
It is better to first parse the input file and then call ``run_pcm``:

.. code-block:: bash

   >>> cd <build_dir>/tests/benchmark

   >>> env PYTHONPATH=../../lib64/python:$PYTHONPATH
       python ../../bin/go_pcm.py --inp=standalone_bubble.pcm

To start collecting hotspots:

.. code-block:: bash

   >>> amplxe-cl -collect hotspots ../../bin/run_pcm @standalone_bubble.pcm

VTune will generate a folder ``r000hs`` with the collected results. A report
for the hotspots can be generated with:

.. code-block:: bash

   >>> amplxe-cl -report hotspots -r r000hs > report
.. _fortran-example:

Interfacing with a Fortran host
===============================

.. literalinclude:: ../../tests/Fortran_host/Fortran_host.f90
   :language: Fortran
   :linenos:
Input description
=================

PCMSolver needs a number of input parameters at runtime. The API provides two
ways of providing them:

1. by means of an additional input file, parsed by the ``go_pcm.py`` script;
2. by means of a special section in the host program input.

Method 1 is more flexible: all parameters that can be modified by the user
are available.  The host program needs only copy the additional input file to
the scratch directory before execution.  Method 2 just gives access to the
core parameters.

In this page, input style and input parameters available in Method 1 will be
documented.

Note that it is also possible to run the module standalone and use a classical charge
distribution.
The classical charge distribution can be specified by giving a molecular geometry
in the molecule section and an additional point multipoles distribution
in the charge distribution section.
The ``run_pcm`` executable has to be compiled for a standalone run with:

.. code-block:: bash

   python <build-path/bin>/go_pcm.py --exe <build-path/bin> --inp molecule.inp 

where the ``molecule.inp`` input file looks like:

.. literalinclude:: ../snippets/molecule.inp

The script and the executable do not need to be in the same directory.

Input style
-----------

The input for PCMSolver is parsed through the `Getkw
<https://github.com/juselius/libgetkw>`_ library written by Jonas Juselius
and is organized in **sections** and **keywords**.  Input reading is
case-insensitive. An example input structure is shown below, there are also
some working examples in the directory ``examples``.  A general input
parameter has the following form (Keyword = [Data type]):

.. literalinclude:: ../snippets/example_input-all.inp

Array-valued keywords will expect the array to be given in comma-separated
format and enclosed in square brackets. The purpose of tags is to distinguish
between cases in which multiple instances of the same kind of object can be
managed by the program.  There exist only certain legal tagnames and these
are determined in the C++ code.  Be aware that the input parsing script does
not check the correctness of tags.

Input parameters
================

Available sections:

+ top section: sets up parameters affecting the module globally;
+ Cavity: sets up all information needed to form the cavity and discretize
  its surface;
+ Medium: sets up the solver to be used and the properties of the medium,
  i.e. the Green's functions inside and outside the cavity;
+ Green, subsection of medium. Sets up the Green's function inside and
  outside the cavity.
+ Molecule: molecular geometry to be used in a standalone run.
+ ChargeDistribution: sets up a classical multipolar (currently up to dipoles)
  charge distribution to use as additional source of electrostatic potential.

.. note::

   The Molecule and ChargeDistribution sections only make sense in a standalone run,
   i.e. when using the ``run_pcm`` executable.

.. warning::

   Exactly matching results obtained from implementations of IEFPCM and/or
   CPCM (COSMO) given in other program packages requires careful selection of
   all the parameters involved.
   A partial checklist of parameters you should always keep in mind:

   * solvent permittivities (static and optical)
   * atomic radii set
   * scaling of the atomic radii
   * cavity surface
   * cavity partition (tesselation)
   * PCM matrix formation algorithm
   * strategy used to solve the PCM linear equations system.

Top section keywords
--------------------

.. glossary::

   Units
      Units of measure used in the input file. If Angstrom is given, all
      relevant input parameters are first converted in au and subsequently
      parsed.

      * **Type**: string
      * **Valid values**: AU | Angstrom
      * **Default**: No Default

   CODATA
      Set of fundamental physical constants to be used in the module.

      * **Type**: integer
      * **Valid values**: 2010 | 2006 | 2002 | 1998
      * **Default**: 2010

Cavity section keywords
-----------------------

.. glossary::

   Type
     The type of the cavity. Completely specifies type of molecular surface
     and its discretization.  Only one type is allowed.  Restart cavity will
     read the file specified by NpzFile keyword and create a GePol cavity
     from that.

     * **Type**: string
     * **Valid values**: GePol | Restart
     * **Default**: none

   NpzFile
     The name of the ``.npz`` file to be used for the GePol cavity restart.

     * **Type**: string
     * **Default**: empty string

   Area
     Average area (weight) of the surface partition for the GePol cavity.

     * **Type**: double
     * **Valid values**: :math:`d \geq 0.01\,\text{a.u.}^2`
     * **Valid for**: GePol cavity
     * **Default value**: :math:`0.3\,\text{a.u.}^2`

   Scaling
     If true, the radii for the spheres will be scaled by 1.2. For finer
     control on the scaling factor for each sphere, select explicit creation
     mode.

     * **Type**: bool
     * **Valid for**: all cavities except Restart
     * **Default value**: True

   RadiiSet
     Select set of atomic radii to be used. Currently Bondi-Mantina
     :cite:`Bondi1964,Mantina2009`, UFF :cite:`Rappe1992` and Allinger's MM3
     :cite:`Allinger1994-tx` sets available, see :ref:`available-radii`.

     * **Type**: string
     * **Valid values**: Bondi | UFF | Allinger
     * **Valid for**: all cavities except Restart
     * **Default value**: Bondi

     .. note::

        Radii in Allinger's MM3 set are obtained by **dividing** the value in the original
        paper by 1.2, as done in the `ADF COSMO implementation <https://www.scm.com/doc/ADF/Input/COSMO.html>`_
        We advise to turn off scaling of the radii by 1.2 when using this set.

   MinRadius
     Minimal radius for additional spheres not centered on atoms. An
     arbitrarily big value is equivalent to switching off the use of added
     spheres, which is the default.

     * **Type**: double
     * **Valid values**: :math:`d \geq 0.4\,\text{a.u.}`
     * **Valid for**: GePol cavity
     * **Default value**: :math:`100.0\,\text{a.u.}`

   Mode
     How to create the list of spheres for the generation of the molecular
     surface:

     + in Implicit mode, the atomic coordinates and charges will be obtained
       from the QM host program.  Spheres will be centered on the atoms and
       the atomic radii, as specified in one the built-in sets, will be used.
       Scaling by 1.2 will be applied according to the keyword Scaling;
     + in Atoms mode, the atomic coordinates and charges will be obtained
       from the QM host program.  For the atoms specified by the array given
       in keyword Atoms, the built-in radii will be substituted by the radii
       provided in the keyword Radii.  Scaling by 1.2 will be applied
       according to the keyword Scaling;
     + in Explicit mode, both centers and radii of the spheres are to be
       specified in the keyword Spheres.  The user has full control over the
       generation of the list of spheres. Scaling by 1.2 is **not** applied,
       regardless of the value of the Scaling keyword.

     * **Type**: string
     * **Valid values**: Implicit | Atoms | Explicit
     * **Valid for**: all cavities except Restart
     * **Default value**: Implicit

   Atoms
     Array of atoms whose radius has to be substituted by a custom value.

     * **Type**: array of integers
     * **Valid for**: all cavities except Restart

   Radii
     Array of radii replacing the built-in values for the selected atoms.

     * **Type**: array of doubles
     * **Valid for**: all cavities except Restart

   Spheres
     Array of coordinates and centers for construction of the list of spheres
     in explicit mode.  Format is :math:`[\ldots, x_i, y_i, z_i, R_i, \ldots]`

     * **Type**: array of doubles
     * **Valid for**: all cavities except Restart

Medium section keywords
-----------------------

.. glossary::

   SolverType
     Type of solver to be used. All solvers are based on the Integral Equation Formulation of
     the Polarizable Continuum Model :cite:`Cances1998`

     + IEFPCM. Collocation solver for a general dielectric medium
     + CPCM. Collocation solver for a conductor-like approximation to the dielectric medium

     * **Type**: string
     * **Valid values**: IEFPCM | CPCM
     * **Default value**: IEFPCM

   Nonequilibrium
     Initializes an additional solver using the dynamic permittivity.
     To be used in response calculations.

     * **Type**: bool
     * **Valid for**: all solvers
     * **Default value**: False

   Solvent
     Specification of the dielectric medium outside the cavity. This keyword
     **must always** be given a value.
     If the solvent name given is different from Explicit any other settings in
     the Green's function section will be overridden by the built-in values for
     the solvent specified. See Table :ref:`available-solvents` for details.
     ``Solvent = Explicit``, triggers parsing of the Green's function sections.

     * **Type**: string
     * **Valid values**:

        + Water                , H2O;
        + Propylene Carbonate  , C4H6O3;
        + Dimethylsulfoxide    , DMSO;
        + Nitromethane         , CH3NO2;
        + Acetonitrile         , CH3CN;
        + Methanol             , CH3OH;
        + Ethanol              , CH3CH2OH;
        + Acetone              , C2H6CO;
        + 1,2-Dichloroethane   , C2H4CL2;
        + Methylenechloride    , CH2CL2;
        + Tetrahydrofurane     , THF;
        + Aniline              , C6H5NH2;
        + Chlorobenzene        , C6H5CL;
        + Chloroform           , CHCL3;
        + Toluene              , C6H5CH3;
        + 1,4-Dioxane          , C4H8O2;
        + Benzene              , C6H6;
        + Carbon Tetrachloride , CCL4;
        + Cyclohexane          , C6H12;
        + N-heptane            , C7H16;
        + Explicit.

   MatrixSymm
     If True, the PCM matrix obtained by the IEFPCM collocation solver is
     symmetrized :math:`\mathbf{K} := \frac{\mathbf{K} + \mathbf{K}^\dagger}{2}`

     * **Type**: bool
     * **Valid for**: IEFPCM solver
     * **Default**: True

   Correction
     Correction, :math:`k` for the apparent surface charge scaling factor in the
     CPCM solver :math:`f(\varepsilon) = \frac{\varepsilon - 1}{\varepsilon +
     k}`

     * **Type**: double
     * **Valid values**: :math:`k > 0.0`
     * **Valid for**: CPCM solver
     * **Default**: 0.0

   DiagonalIntegrator
     Type of integrator for the diagonal of the boundary integral operators

     * **Type**: string
     * **Valid values**: COLLOCATION
     * **Valid for**: IEFPCM, CPCM
     * **Default**: COLLOCATION
     * **Notes**: in future releases we will add PURISIMA and NUMERICAL as options

   DiagonalScaling
     Scaling factor for diagonal of collocation matrices

     * **Type**: double
     * **Valid values**: :math:`f > 0.0`
     * **Valid for**: IEFPCM, CPCM
     * **Default**: 1.07
     * **Notes**: values commonly used in the literature are 1.07 and 1.0694

   ProbeRadius
     Radius of the spherical probe approximating a solvent molecule. Used for
     generating the solvent-excluded surface (SES) or an approximation of it.
     Overridden by the built-in value for the chosen solvent.

     * **Type**: double
     * **Valid values**: :math:`d \in [0.1, 100.0]\,\text{a.u.}`
     * **Valid for**: all solvers
     * **Default**: 1.0

Green section keywords
----------------------

If ``Solvent = Explicit``, **two** Green's functions sections must be specified
with tags ``inside`` and ``outside``, i.e. ``Green<inside>`` and
``Green<outside>``.  The Green's function inside will always be the vacuum,
while the Green's function outside might vary.

.. glossary::

   Type
     Which Green's function characterizes the medium.

     * **Type**: string
     * **Valid values**: Vacuum | UniformDielectric | SphericalDiffuse | SphericalSharp
     * **Default**: Vacuum

   Der
     How to calculate the directional derivatives of the Green's function:

       + Numerical, perform numerical differentiation **debug option**;
       + Derivative, use automatic differentiation to get the directional derivative;
       + Gradient, use automatic differentiation to get the full gradient **debug option**;
       + Hessian, use automatic differentiation to get the full hessian **debug option**;

     * **Type**: string
     * **Valid values**: Numerical | Derivative | Gradient | Hessian
     * **Default**: Derivative

     .. note::

        The spherical diffuse Green's function **always** uses numerical differentiation.


   Eps
     Static dielectric permittivity of the medium

     * **Type**: double
     * **Valid values**: :math:`\varepsilon \geq 1.0`
     * **Default**: 1.0

   EpsDyn
     Dynamic dielectric permittivity of the medium

     * **Type**: double
     * **Valid values**: :math:`\varepsilon \geq 1.0`
     * **Default**: 1.0

   Profile
      Functional form of the dielectric profile

      * **Type**: string
      * **Valid values**: Tanh | Erf | Log
      * **Valid for**: SphericalDiffuse
      * **Default**: Log

   Eps1
     Static dielectric permittivity inside the interface

     * **Type**: double
     * **Valid values**: :math:`\varepsilon \geq 1.0`
     * **Valid for**: SphericalDiffuse, SphericalSharp
     * **Default**: 1.0

   EpsDyn1
     Dynamic dielectric permittivity inside the interface

     * **Type**: double
     * **Valid values**: :math:`\varepsilon \geq 1.0`
     * **Valid for**: SphericalDiffuse, SphericalSharp
     * **Default**: 1.0

   Eps2
     Static dielectric permittivity outside the interface

     * **Type**: double
     * **Valid values**: :math:`\varepsilon \geq 1.0`
     * **Valid for**: SphericalDiffuse, SphericalSharp
     * **Default**: 1.0

   EpsDyn2
     Dynamic dielectric permittivity outside the interface

     * **Type**: double
     * **Valid values**: :math:`\varepsilon \geq 1.0`
     * **Valid for**: SphericalDiffuse, SphericalSharp
     * **Default**: 1.0

   Center
     Center of the interface layer. This corresponds to the radius
     of the spherical droplet.

     * **Type**: double
     * **Valid for**: SphericalDiffuse, SphericalSharp
     * **Default**: 100.0 a.u.

   Width
     Physical width of the interface layer. This value is divided by 6.0
     internally.

     * **Type**: double
     * **Valid for**: SphericalDiffuse
     * **Default**: 5.0 a.u.

     .. warning::

        Numerical instabilities may arise if a too small value is selected.

   InterfaceOrigin
     Center of the spherical droplet

     * **Type**: array of doubles
     * **Valid for**: SphericalDiffuse, SphericalSharp
     * **Default**: :math:`[0.0, 0.0, 0.0]`

   MaxL
     Maximum value of the angular momentum in the expansion of the
     Green's function for the spherical diffuse Green's function

     * **Type**: integer
     * **Valid for**: SphericalDiffuse, SphericalSharp
     * **Default**: 30

Molecule section keywords
-------------------------

It is possible to run the module standalone and use a classical charge
distribution as specified in this section of the input.
The ``run_pcm`` executable has to be compiled for a standalone run with:

.. code-block:: bash

   python go_pcm.py -x molecule.inp

where the ``molecule.inp`` input file looks like:

.. literalinclude:: ../snippets/molecule.inp

.. glossary::

   Geometry
     Coordinates and charges of the molecular aggregate.
     Format is :math:`[\ldots, x_i, y_i, z_i, Q_i, \ldots]`
     Charges are always assumed to be in atomic units

     * **Type**: array of doubles

ChargeDistribution section keywords
-----------------------------------

Set a classical charge distribution, inside or outside the cavity
No additional spheres will be generated.

.. glossary::

   Monopoles
     Array of point charges
     Format is :math:`[\ldots, x_i, y_i, z_i, Q_i, \ldots]`

     * **Type**: array of doubles

   Dipoles
     Array of point dipoles.
     Format is :math:`[\ldots, x_i, y_i, z_i, \mu_{x_i}, \mu_{y_i}, \mu_{z_i} \ldots]`
     The dipole moment components are always read in atomic units.

     * **Type**: array of doubles

MMFQ section keywords
---------------------

Set a classical fluctuating charge force field. This is incompatible with any
options specifying a continuum model.  No additional spheres will be generated.

.. glossary::

   SitesPerFragment
     Number of sites per MM fragment. For water this is 3.

     * **Type**: integer
     * **Default**: 3

   Sites
     Array of MM sites for the FQ model
     Format is :math:`[\ldots, x_i, y_i, z_i, chi_i, eta_i \ldots]`

     * **Type**: array of doubles

   NonPolarizable
     Whether to make this force field nonpolarizable.

     * **Type**: bool
     * **Default**: false


.. _available-radii:

Available radii
---------------

.. image::  ../gfx/bondi_mantina.png
   :scale: 70 %
   :align: center

.. image::  ../gfx/uff.png
   :scale: 70 %
   :align: center

.. image::  ../gfx/allingerMM3.png
   :scale: 70 %
   :align: center

.. _available-solvents:

Available solvents
------------------

The macroscopic properties for the built-in list of solvents are:

  + static permittivity, :math:`\varepsilon_s`
  + optical permittivity, :math:`\varepsilon_\infty`
  + probe radius, :math:`r_\mathrm{probe}` in Angstrom.

The following table summarizes the built-in solvents and their properties.
Solvents are ordered by decreasing static permittivity.

 ==================== ======== ===================== ========================== ========================
 Name                 Formula  :math:`\varepsilon_s` :math:`\varepsilon_\infty` :math:`r_\mathrm{probe}`
 ==================== ======== ===================== ========================== ========================
 Water                H2O              78.39                 1.776                     1.385
 Propylene Carbonate  C4H6O3           64.96                 2.019                     1.385
 Dimethylsulfoxide    DMSO             46.7                  2.179                     2.455
 Nitromethane         CH3NO2           38.20                 1.904                     2.155
 Acetonitrile         CH3CN            36.64                 1.806                     2.155
 Methanol             CH3OH            32.63                 1.758                     1.855
 Ethanol              CH3CH2OH         24.55                 1.847                     2.180
 Acetone              C2H6CO           20.7                  1.841                     2.38
 1,2-Dichloroethane   C2H4Cl2          10.36                 2.085                     2.505
 Methylenechloride    CH2Cl2            8.93                 2.020                     2.27
 Tetrahydrofurane     THF               7.58                 1.971                     2.9
 Aniline              C6H5NH2           6.89                 2.506                     2.80
 Chlorobenzene        C6H5Cl            5.621                2.320                     2.805
 Chloroform           CHCl3             4.90                 2.085                     2.48
 Toluene              C6H5CH3           2.379                2.232                     2.82
 1,4-Dioxane          C4H8O2            2.250                2.023                     2.630
 Benzene              C6H6              2.247                2.244                     2.630
 Carbon tetrachloride CCl4              2.228                2.129                     2.685
 Cyclohexane          C6H12             2.023                2.028                     2.815
 N-heptane            C7H16             1.92                 1.918                     3.125
 ==================== ======== ===================== ========================== ========================
.. _C-example:

Interfacing with a C host
=========================

.. warning::

   Multidimensional arrays are handled in *column-major ordering*
   (i.e. Fortran ordering) by the module.

.. literalinclude:: ../../tests/C_host/C_host.c
   :language: C
   :linenos:
Interfacing a QM program and PCMSolver
======================================

For the impatients: tl;dr
-------------------------

In these examples, we want to show how *every function* in the API works.
If your program is written in Fortran, head over to :ref:`fortran-example`
If your program is written in C/C++, head over to :ref:`C-example`

How PCMSolver handles potentials and charges: surface functions
---------------------------------------------------------------

Electrostatic potential vectors and the corresponding apparent surface
charge vectors are handled internally as `surface functions`.
The actual values are stored into Eigen vectors and saved into a
map. The mapping is between the name of the surface function, given by
the programmer writing the interface to the library, and the vector holding
the values.

What you should care about: API functions
-----------------------------------------

These are the contents of the ``pcmsolver.h`` file defining
the public API of the PCMSolver library. The Fortran bindings
for the API are in the ``pcmsolver.f90`` file.
The indexing of symmetry operations and their mapping to a bitstring
is explained in the following Table. This is important when passing
symmetry information to the :cpp:func:`pcmsolver_new` function.

.. _symmetry-ops:
.. table:: Symmetry operations indexing within the module

   ===== === ========= ======
   Index zyx Generator Parity
   ===== === ========= ======
     0   000     E       1.0
     1   001    Oyz     -1.0
     2   010    Oxz     -1.0
     3   011    C2z      1.0
     4   100    Oxy     -1.0
     5   101    C2y      1.0
     6   110    C2x      1.0
     7   111     i      -1.0
   ===== === ========= ======


.. doxygenfile:: mock_pcmsolver.h
   :project: PCMSolver

Host input forwarding
---------------------

.. doxygenstruct:: PCMInput
   :project: PCMSolver

Internal details of the API
---------------------------

.. doxygenclass:: pcm::Meddle
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

.. doxygenclass:: pcm::Input
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:
Building the module
===================

PCMSolver configuration and build process is managed through CMake.

Prerequisites and dependencies
------------------------------

A number of prerequisites and dependencies are to be satisfied to successfully
build the module. It will be here assumed that you want to perform a "full"
build, i.e. you want to build the static libraries to be linked to your QM
program, the unit test suite and an offline copy of this documentation.

Compilers
~~~~~~~~~

+ a C++ compiler, compliant with the 2011 ISO C++ standard. The build system
  will downgrade to using the 1998 ISO C++ standard plus the 2003 technical
  corrigendum and some additional defect reports, if no suitable support if
  found.

  .. warning::

     Backwards compatibility support for the C++03 standard is **deprecated**
     and will be removed in upcoming releases of the library.

+ a C compiler, compliant with the ISO C99 standard.
+ a Fortran compiler, compliant with the Fortran 2003 standard.

The list of primary test environments can be found in the `README.md
<https://github.com/PCMSolver/pcmsolver/blob/master/README.md>`_ file. It is
entirely possible that using other compiler versions you might be able to build
the module. In order to ensure that you have a sane build, you will have to run
the unit test suite.

Libraries and toolchain programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+ CMake version 3.3 and higher;
+ Git version 1.7.1 and higher;
+ Python interpreter 2.7 and higher;
+ Boost libraries version 1.54.0 and higher;

.. note::

   Version 1.54.0 of Boost libraries is shipped with the module and resides in the ``cmake/downloaded`` subdirectory.
   Unless you want to use another version of Boost, you should not worry about satisfying this dependency.

+ `zlib <http://www.zlib.net/>`_ version 1.2 and higher (unit test suite only);
+ Doxygen version 1.7.6 and higher (documentation only)
+ Perl (documentation only)
+ Sphinx (documentation only)

PCMSolver relies on the Eigen template libraries version 3.3.0 and higher.
Version 3.3.0 of Eigen libraries is shipped with the module and resides in the ``external`` subdirectory.

Configuration
-------------

Configuration is managed through the front-end script ``setup.py`` residing in the
repository main directory. Issuing:

.. code-block:: bash

   ./setup [options] [build path]

will create the build directory in build path and run CMake with the given
options. By default, files are configured in the ``build`` directory. The ``-h`` or
``--help`` option will list the available options and their effect. Options can
be forwarded directly to CMake by using the ``--cmake-options`` flag and listing
the ``-D...`` options. Usually the following command is sufficient to get the
configuration done for a debug build, including compilation of the unit test
suite:

.. code-block:: bash

   ./setup --type=debug

The unit tests suite is **always** compiled in standalone mode, unless the
``-DENABLE_TESTS=OFF`` option is forwarded to CMake.

Getting Boost
~~~~~~~~~~~~~

You can get Boost libraries in two ways:

 + already packaged by your Linux distribution or through MacPorts/Brew;
 + by downloading the archive from http://www.boost.org/ and building it yourself.

In case your distribution packages a version older than 1.54.0 you might chose
to either build Boost on your own or to rely on the automated build of the
necessary Boost libraries when compiling the module (recommended).  Full
documentation on how to build Boost on Unix variants is available
`here <http://www.boost.org/doc/libs/1_56_0/more/getting_started/unix-variants.html>`_.
It is here assumed that the user **does not** have root access to the machine
and will install the libraries to a local prefix, a subdirectory of
``/home/user-name`` tipically.
Once you've downloaded and unpacked the archive, run the bootstrap script to configure:

.. code-block:: bash

   cd path/to/boost
   ./bootstrap.sh --prefix=/home/user-name/boost

Running ``./bootstrap.sh --help`` will list the available options for the script. To build run:

.. code-block:: bash

   ./b2 install

This might take a while. After a successful build you will find the headers in
``/home/user-name/boost/include`` and libraries in ``/home/user-name/boost/lib``
Now, you will have Boost in a nonstandard location. Without hints CMake will
not be able to find it and configuration of `PCMSolver` will fail.  To avoid
this, you will have to pass the location of the headers and libraries to the
setup script, either with:

.. code-block:: bash

   ./setup --boost-headers=/home/user-name/boost/include --boost-libs=/home/user-name/boost/lib

or with:

.. code-block:: bash

   ./setup -DBOOST_INCLUDEDIR=/home/user-name/boost/include -DBOOST_LIBRARYDIR=/home/user-name/boost/lib

Advanced configuration options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These options are marked as advanced as it is highly unlikely they will
be useful when not programming the library:

* ``--exdiag`` Enable C++ extended diagnostics flags. Disabled by default.
* ``--ccache`` Enable use of ccache for C/C++ compilation caching.
  Enabled by default, unless ccache is not available.
* ``--build-boost`` Deactivate Boost detection and build on-the-fly. Disabled by default.
* ``--eigen`` Root directory for Eigen3. Search for Eigen3 in the location provided by the
  user. If search fails, fall back to the version bundled with the library.
* ``--static`` Create only static library. Disabled by default.

Some options can only be tweaked `via` ``--cmake-options`` to the setup script:

* ``ENABLE_DOCS`` Enable build of documentation. This requires a number of additional dependencies.
  If any of these are not met, documentation is not built. Enabled by default.
* ``ENABLE_LOGGER`` Enable compilation of logger sources. Disabled by default.

  .. warning::

     The logger is not currently in use in any part of the code.

* ``ENABLE_TIMER`` Enable compilation of timer sources. Enabled by default.
* ``BUILD_STANDALONE`` Enable compilation of standalone ``run_pcm`` executable. Enabled by default.
* ``TEST_Fortran_API`` Test the Fortran 90 bindings for the API. Enabled by default.
* ``ENABLE_GENERIC`` Enable mostly static linking in shared library. Disabled by default.
* ``ENABLE_TESTS`` Enable compilation of unit tests suite. Enabled by default.
* ``SHARED_LIBRARY_ONLY`` Create only shared library. Opposite of ``--static``.
* ``PYMOD_INSTALL_LIBDIR`` *If set*, installs python scripts/modules to
  ``${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pcmsolver`` rather than the
  default ``${CMAKE_INSTALL_BINDIR}`` (i.e., ``bin``).
* ``CMAKE_INSTALL_BINDIR`` Where to install executables, if not to ``bin``.
* ``CMAKE_INSTALL_LIBDIR`` Where to install executables, if not to ``bin``.
* ``CMAKE_INSTALL_INCLUDESDIR`` Where to install executables, if not to ``bin``.

* ``CMAKE_INSTALL_BINDIR`` Location within ``CMAKE_INSTALL_PREFIX`` (``--prefix``) to
  which executables are installed (default: ``bin``).
* ``CMAKE_INSTALL_LIBDIR`` Location within ``CMAKE_INSTALL_PREFIX`` (``--prefix``) to
  which libraries are installed (default: ``lib``).
* ``CMAKE_INSTALL_INCLUDEDIR`` Location within ``CMAKE_INSTALL_PREFIX`` (``--prefix```)
  to which headers are installed (default: ``include``).
* ``PYMOD_INSTALL_LIBDIR`` *If set*, location within ``CMAKE_INSTALL_LIBDIR`` to which
  python modules are installed,
  ``${CMAKE_INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}/pcmsolver``. *If not set*,
  python modules installed to default ``${CMAKE_INSTALL_LIBDIR}/python/pcmsolver``.

Build and test
--------------

To compile and link, just go to the build directory and run:

.. code-block:: bash

   make -j N

where ``N`` is the number of cores you want to use when building.

.. note::

   Building on more than one core can sometimes result in a "race condition"
   and a crash. If that happens, please report the problem as an issue on our
   issue tracker on GitHub. Running ``make`` on a single core might get you through
   compilation.

To run the whole test suite:

.. code-block:: bash

   ctest -j N

You can also use CTest to run a specific test or a set of tests. For example:

.. code-block:: bash

   ctest -R gepol

will run all the test containing the string "gepol" in their name.

=======================
PCMSolver Users' Manual
=======================

.. toctree::
   building
   input
   interfacing
   fortran-example
   C-example
Helper classes and functions
============================

.. image:: ../gfx/bar_charts/utils.svg
   :scale: 70 %
   :align: center

Sphere
------
.. doxygenstruct:: pcm::utils::Sphere
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Atom
----
.. doxygenstruct:: pcm::utils::Atom
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:


ChargeDistribution
------------------
.. doxygenstruct:: pcm::utils::ChargeDistribution
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Molecule
--------
.. doxygenclass:: pcm::Molecule
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Solvent
-------
.. doxygenstruct:: pcm::utils::Solvent
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Symmetry
--------
.. doxygenclass:: Symmetry
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Mathematical utilities
----------------------
.. doxygenfile:: MathUtils.hpp
   :project: PCMSolver
Boundary integral operators
===========================

.. image:: ../gfx/bar_charts/bi_operators.svg
   :scale: 70 %
   :align: center

IBoundaryIntegralOperator
-------------------------
.. doxygenclass:: pcm::IBoundaryIntegralOperator
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Collocation
-----------

.. doxygenclass:: pcm::bi_operators::Collocation
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Purisima
--------

.. doxygenclass:: pcm::bi_operators::Purisima
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Numerical
---------

.. doxygenclass:: pcm::bi_operators::Numerical
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:
Solvers
=======

We will here describe the inheritance hierarchy for generating solvers, in
order to use and extend it properly.  The runtime creation of solver objects
relies on the Factory Method pattern :cite:`Gamma1994,Alexandrescu2001`,
implemented through the generic Factory class.

.. image:: ../gfx/bar_charts/solver.svg
   :scale: 70 %
   :align: center

ISolver
-------
.. doxygenclass:: pcm::ISolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

IEFSolver
---------
.. doxygenclass:: pcm::solver::IEFSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

CPCMSolver
----------
.. doxygenclass:: pcm::solver::CPCMSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:
Cavities
========

We will here describe the inheritance hierarchy for generating cavities, in
order to use and extend it properly.  The runtime creation of cavity objects
relies on the Factory Method pattern :cite:`Gamma1994,Alexandrescu2001`,
implemented through the generic Factory class.

.. image:: ../gfx/bar_charts/cavity.svg
   :scale: 70 %
   :align: center

ICavity
-------
.. doxygenclass:: pcm::ICavity
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

GePolCavity
-----------

.. doxygenclass:: pcm::cavity::GePolCavity
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

RestartCavity
-------------

.. doxygenclass:: pcm::cavity::RestartCavity
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:
===============================
Classes and functions reference
===============================

.. image:: ../gfx/bar_charts/total.svg
   :scale: 70 %
   :align: center

.. toctree::

   cavities
   greens-functions
   dielectric-profiles
   solvers
   bi-operators
   helper-classes
   namespaces
Dielectric profiles
===================

.. image:: ../gfx/bar_charts/dielectric_profile.svg
   :scale: 70 %
   :align: center

Uniform
-------
.. doxygenstruct:: pcm::dielectric_profile::Uniform
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Anisotropic
-----------
.. doxygenclass:: pcm::dielectric_profile::Anisotropic
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Yukawa
-------
.. doxygenstruct:: pcm::dielectric_profile::Yukawa
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

OneLayerLog
-----------
.. doxygenclass:: pcm::dielectric_profile::OneLayerLog
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

OneLayerTanh
------------
.. doxygenclass:: pcm::dielectric_profile::OneLayerTanh
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

OneLayerErf
-----------
.. doxygenclass:: pcm::dielectric_profile::OneLayerErf
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Sharp
-----
.. doxygenstruct:: pcm::dielectric_profile::Sharp
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:
Green's Functions
=================

We will here describe the inheritance hierarchy for generating Green's
functions, in order to use and extend it properly.  The runtime creation of
Green's functions objects relies on the Factory Method pattern
:cite:`Gamma1994,Alexandrescu2001`, implemented through the
generic Factory class.

The top-level header, _i.e._ to be included in client code, is ``Green.hpp``.
The common interface to all Green's function classes is specified by the ``IGreensFunction`` class,
this is non-templated.
All other classes are templated.
The Green's functions are registered to the factory based on a label encoding: type, derivative, and dielectric profile.
The only allowed labels must be listed in ``src/green/Green.hpp``. If they are not, they can not be selected at run time.

.. image:: ../gfx/bar_charts/green.svg
   :scale: 70 %
   :align: center

IGreensFunction
---------------

.. doxygenclass:: pcm::IGreensFunction
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

GreensFunction
--------------
.. doxygenclass:: pcm::green::GreensFunction
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

Vacuum
------
.. doxygenclass:: pcm::green::Vacuum
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

UniformDielectric
-----------------
.. doxygenclass:: pcm::green::UniformDielectric
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

IonicLiquid
-----------
.. doxygenclass:: pcm::green::IonicLiquid
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

AnisotropicLiquid
-----------------
.. doxygenclass:: pcm::green::AnisotropicLiquid
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

SphericalDiffuse
----------------
.. doxygenclass:: pcm::green::SphericalDiffuse
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:
