# freud Contributor Agreement

These terms apply to your contribution to the freud Open Source Project ("Project") owned or managed by the Regents of the University of Michigan ("Michigan"), and set out the intellectual property rights you grant to Michigan in the contributed materials. If this contribution is on behalf of a company, the term "you" will also mean the company you identify below. If you agree to be bound by these terms, fill in the information requested below and provide your signature.

1. The term "contribution" means any source code, object code, patch, tool, sample, graphic, specification, manual, documentation, or any other material posted or submitted by you to a project.
2. With respect to any worldwide copyrights, or copyright applications and registrations, in your contribution:
    * you hereby assign to Michigan joint ownership, and to the extent that such assignment is or becomes invalid, ineffective or unenforceable, you hereby grant to Michigan a perpetual, irrevocable, non-exclusive, worldwide, no-charge, royalty-free, unrestricted license to exercise all rights under those copyrights. This includes, at Michigan's option, the right to sublicense these same rights to third parties through multiple levels of sublicensees or other licensing arrangements;
    * you agree that both Michigan and you can do all things in relation to your contribution as if each of us were the sole owners, and if one of us makes a derivative work of your contribution, the one who makes the derivative work (or has it made) will be the sole owner of that derivative work;
    * you agree that you will not assert any moral rights in your contribution against us, our licensees or transferees;
    * you agree that we may register a copyright in your contribution and exercise all ownership rights associated with it; and
    * you agree that neither of us has any duty to consult with, obtain the consent of, pay or render an accounting to the other for any use or distribution of your contribution.
3. With respect to any patents you own, or that you can license without payment to any third party, you hereby grant to Michigan a perpetual, irrevocable, non-exclusive, worldwide, no-charge, royalty-free license to:
    * make, have made, use, sell, offer to sell, import, and otherwise transfer your contribution in whole or in part, alone or in combination with or included in any product, work or materials arising out of the project to which your contribution was submitted; and
    * at Michigan's option, to sublicense these same rights to third parties through multiple levels of sublicensees or other licensing arrangements.
4. Except as set out above, you keep all right, title, and interest in your contribution. The rights that you grant to Michigan under these terms are effective on the date you first submitted a contribution to Michigan, even if your submission took place before the date you sign these terms. Any contribution Michigan makes available under any license will also be made available under a suitable Free Software Foundation or Open Source Initiative approved license.
5. With respect to your contribution, you represent that:
    * it is an original work and that you can legally grant the rights set out in these terms;
    * it does not to the best of your knowledge violate any third party's copyrights, trademarks, patents, or other intellectual property rights; and
you are authorized to sign this contract on behalf of your company (if identified below).
6. The terms will be governed by the laws of the State of Michigan and applicable U.S. Federal Law. Any choice of law rules will not apply.

**By making contribution, you electronically sign and agree to the terms of the freud Contributor Agreement.**

![by-sa.png](https://licensebuttons.net/l/by-sa/3.0/88x31.png)

Based on the Sun Contributor Agreement - version 1.5.
This document is licensed under a Creative Commons Attribution-Share Alike 3.0 Unported License
https://creativecommons.org/licenses/by-sa/3.0/
Contributions are welcomed via pull requests. First, contact the _freud_ developers prior to beginning
your work to ensure that your plans mesh well with the planned development direction and standards set for the project.
Then implement your code.

Submit a pull request. Multiple developers and/or users will review requested changes and make comments.
The lead developer(s) will merge into the `master` branch after the review is complete and approved.

# Features

## Implement functionality in a general and flexible fashion

The _freud_ library provides a lot of flexibility to the user. Your pull request should provide something that is
applicable to a variety of use-cases and not just the one thing you might need it to do for your research. Speak to the
lead developers before writing your code, and they will help you make design choices that allow flexibility.

## Do not degrade performance of existing code paths

New functionalities should only activate expensive code paths when they are requested by the user.
Do not slow down existing code.

## Add dependencies only if absolutely necessary

In order to make _freud_ as widely available as possible, we try to keep the number of dependencies to a minimum.
If you need a feature present in an external library, follow the following steps:

1. Add to _freud_ itself if it's simple or if other modules would benefit:
    * Example: Added simple tensor math for cubatic order parameter.
2. Add via submodule if the code exists externally:
    * Example: _fsph_ for spherical harmonics.
3. Contact _freud_ developers to inquire if the library you'd like as a dependency fits in with the overall design/goals
of _freud_.

# Version control

## Base your work off the correct branch

Use the [OneFlow](https://www.endoflineblog.com/oneflow-a-git-branching-model-and-workflow) model of development:
  * Both new features and bug fixes should be developed in branches based on `master`.
  * Hotfixes (critical bugs that need to be released *fast*) should be developed in a branch based on the latest tagged release.

## Propose a single set of related changes

Changes proposed in a single topic branch / pull request should all be related to each other. Don't propose too
many changes at once, review becomes challenging. Multiple new features that are loosely coupled should be completed
in separate topic branches. It is OK if the branch for `feature2` is based on `feature1` - as long as it is made clear
that `feature1` should be merged before the review of `feature2`. It is better to merge both `feature1` and `feature2`
into a temporary integration branch during testing.

## Keep changes to a minimum

Avoid whitespace changes and use separate pull requests to propose fixes in spelling, grammar, etc. :)

## Agree to the contributor agreement

All contributors must agree to the Contributor Agreement ([ContributorAgreement.md](ContributorAgreement.md))
before their pull request can be merged.

# Source code

## Use a consistent style

It is important to have a consistent style throughout the source code.
Follow the source conventions defined in the documentation for all code.

## Document code with comments

Add comments to code so that other developers can understand it.

## Compiles without warnings

Your changes should compile without warnings.

# Tests

## Write unit tests

All new functionality in _freud_ should be tested with automatic unit tests that execute in a few seconds (if your
specific test requires a long amount of time, please alert the _freud_ developers as to why this is required so that
your test can be opted-out of for "regular" unit-testing). High level features should be tested from Python, and the
Python tests should attempt to cover all options that the user can select.

## Validity tests

In addition to the unit tests, the developer should run research-scale analysis using the new functionality and
ensure that it behaves as intended.

# User documentation

## Write user documentation

User documentation for the user facing script commands should be documented with docstrings in [Google format](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).
Include examples on using new functionality.

## Add developer to the credits

Developers need to be credited for their work. Update the [credits documentation](doc/source/reference/credits.rst)
to reference what each developer has contributed to the code.

## Update Change Log

Add a short concise entry describing the change to the [ChangeLog.md](ChangeLog.md).
# Change Log
The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.8.0 -- 2022-01-25

### Added
* `freud.diffraction.StaticStructureFactorDirect` class (unstable) can be used to compute the static structure factor S(k) by sampling reciprocal space vectors.
* Python 3.10 is supported.
* Documentation examples are tested with pytest.
* Use clang-format as pre-commit hook.
* Add related tools section to the documentation.

### Fixed
* `freud.diffraction.DiffractionPattern` normalization changed such that `S(k=0) = N`.
* Added error checking for `r_min`, `r_max` arguments in `freud.density.RDF`, `freud.locality.NeighborList`, `freud.locality.NeighborQuery`, and `freud.density.LocalDensity` classes.
* CMake build system only uses references to TBB target.

### Changed
* Re-organized tests for the static structure factor classes.
* Move `util::Histogram<T>::Axes` to `util::Axes`.
* Use new `flake8` plugin `flake8-force` for linting Cython code.

## v2.7.0 -- 2021-10-01

### Added
* `freud.diffraction.StaticStructureFactorDebye` class (unstable) can be used to compute the static structure factor S(k) using the Debye formula.

### Fixed
* Updated lambda functions to capture `this` by reference, to ensure compatibility with C++20 and above.
* Fixed ``Box.contains`` to run in linear time, ``O(num_points)``.
* Fixed compilation to pass compiler optimization flags when build type is ReleaseWithDocs (major perf regression since 2.4.1).

## v2.6.2 -- 2021-06-26

### Fixed
* Upgrade to auditwheel 4.0.0 in cibuildwheel to ensure RPATH is patched properly for `libfreud.so` in Linux wheels.

## v2.6.1 -- 2021-06-23

### Fixed
* Added missing git submodules to source distribution.

## v2.6.0 - 2021-06-22

### Added
* Added `out` option for the `wrap`, `unwrap`, `make_absolute`, and `make_fractional` methods of `Box`.
* The `Steinhardt` and `SolidLiquid` classes expose the raw `qlmi` arrays.
* The `Steinhardt` class supports computing order parameters for multiple `l`.

### Changed
* Improvements to plotting for the `DiffractionPattern`.
* Wheels are now built with cibuildwheel.

### Fixed
* Fixed/Improved the `k` values and vectors in the `DiffractionPattern` (more improvement needed).
* Fixed incorrect computation of `Steinhardt` averaged quantities. Affects all previous versions of freud 2.
* Fixed documented formulas for `Steinhardt` class.
* Fixed broken arXiv links in bibliography.

## v2.5.1 - 2021-04-06

### Added
* The `compute` method of `DiffractionPattern` class has a `reset` argument.

### Fixed
* Documentation on ReadTheDocs builds and renders.

## v2.5.0 - 2021-03-16

### Changed
* NeighborList `filter` method has been optimized.
* TBB 2021 is now supported (removed use of deprecated TBB features).
* Added new pre-commit hooks for `black`, `isort`, and `pyupgrade`.
* Testing framework now uses `pytest`.

## v2.4.1 - 2020-11-16

### Fixed
* Python 3.8 builds with Windows MSVC were broken due to an unrecognized CMake compiler option.
* Fixed broken documentation by overriding scikit-build options.
* RPATH on Linux is now set correctly to find TBB libraries not on the global search path.
* 2D box image calculations now return zero for the image z value.
* Fixed wrong attribute name in `EnvironmentCluster.plot`.

## v2.4.0 - 2020-11-09

### Added
* The Box class has a method `contains` to determine particle membership in a box.
* NeighborList class exposes `num_points` and `num_query_points` attributes.
* `compute` method of `GaussianDensity` class has a `values` argument.
* Support for pre-commit hooks.
* Python 3.9 is supported.

### Changed
* NeighborList raises a `ValueError` instead of a `RuntimeError` if provided invalid constructor arguments.
* freud now builds using scikit-build (requires CMake).

### Deprecated
* `freud.order.Translational`

### Fixed
* Source distributions now include Cython source files.
* Hexatic order parameter (unweighted) normalizes by number of neighbors instead of the symmetry order k.
* Particles with an i-j normal vector of [0, 0, 0] are excluded from 2D Voronoi NeighborList computations for numerical stability reasons.
* Memory leak in `makeDefaultNlist` function where a NeighborList was being allocated and not freed.

## v2.3.0 - 2020-08-03

### Added
* Support for garnett 0.7.
* Custom NeighborLists can be created from a set of points using `from_points`. Distances will be calculated automatically.
* The Box class has methods `compute_distances` and `compute_all_distances` to calculate distances between arrays of points and query points.
* Hexatic can now compute 2D Minkowski Structure Metrics, using `weighted=True` along with a Voronoi NeighborList.
* Examples have been added to the Cluster, Density, Environment, and Order Modules.
* Module examples have been integrated with doctests to ensure they are up to date with API.
* SphereVoxelization class in the `density` module computes a grid of voxels occupied by spheres.
* `freud.diffraction.DiffractionPattern` class (unstable) can be used to compute 2D diffraction patterns.

### Changed
* Cython is now a required dependency (not optional). Cythonized `.cpp` files have been removed.
* An instance of GaussianDensity cannot compute 3D systems if it has been previously computed 2D systems.

### Fixed
* Histogram bin locations are computed in a more numerically stable way.
* Improved error handling of Cubatic input parameters.
* PMFTs are now properly normalized such that the pair correlation function tends to unity for an ideal gas.
* PMFTXYT uses the correct orientations when `points` and `query_points` differ.
* GaussianDensity Gaussian normalization in 2D systems has been corrected.

### Removed
* Python 3.5 is no longer supported. Python 3.6+ is required.

## v2.2.0 - 2020-02-24

### Added
* NeighborQuery objects can now create NeighborLists with neighbors sorted by bond distance.
* LocalDescriptors `compute` takes an optional maximum number of neighbors to compute for each particle.

### Fixed
* Corrected calculation of neighbor distances in the Voronoi NeighborList.
* Added finite tolerance to ensure stability of 2D Voronoi NeighborList computations.

## v2.1.0 - 2019-12-19

### Added
* The Box class has methods `center_of_mass` and `center` for periodic-aware center of mass and shifting points to center on the origin.

### Changed
* The make\_random\_box system method no longer overwrites the NumPy global random number generator state.
* The face\_orientations argument of PMFTXYZ has been renamed to equiv\_orientations and must be provided as an Mx4 array, where M is the number of symmetrically equivalent particle orientations.
* Improved documentation about query modes.
* The Voronoi class uses smarter heuristics for its voro++ block sizes, resulting in significant performance gains for large systems.

### Fixed
* The from\_box method correctly passes user provided dimensions to from\_matrix it if is called.
* Correctly recognize Ovito DataCollection objects in from\_system.
* Corrected `ClusterProperties` calculation of centers of mass in specific systems.
* Set z positions to 0 for 2D GSD systems in from\_system.
* PMFTXY and PMFTXYZ index into query orientations using the query point index instead of the point index.

## v2.0.1 - 2019-11-08

### Added
* Rewrote development documentation to match the conventions and logic in version 2.0 of the code.

### Fixed
* Automatic conversion of 2D systems from various data sources.
* Mybinder deployment works with freud v2.0.
* Minor errors in freud-examples have been corrected.

## v2.0.0 - 2019-10-31

### Added
* Ability to specify "system-like" objects that contain a box and set of points for most computes.
* NeighborLists and query arguments are now accepted on equal footing by compute methods that involve neighbor finding via the `neighbors=...` argument.
* Extensive new documentation including tutorial for new users and reference sections on crucial topics.
* Standard method for preprocessing arguments of pair computations.
* New internal ManagedArray object that allows data persistence and improves indexing in C++.
* Internal threaded storage uses the standard ManagedArray object.
* C++ Histogram class to standardize n-dimensional binning and simplify writing new methods.
* Upper bound r\_max option for number of neighbors queries.
* Lower bound r\_min option for all queries.
* Steinhardt now supports l = 0, 1.
* C++ BondHistogramCompute class encapsulates logic of histogram-based methods.
* 2D PMFTs accept quaternions as well as angles for their orientations.
* ClusterProperties computes radius of gyration from the gyration tensor for each cluster.
* `freud.data` module for generating example particle systems.
* Optional normalization for RDF, useful for small systems.
* `plot()` methods for `NeighborQuery` and `Box` objects.
* Added support for reading system data directly from MDAnalysis, garnett, gsd, HOOMD-blue, and OVITO.
* Various validation tests.

### Changed
* All compute objects that perform neighbor computations now use NeighborQuery internally.
* Neighbor-based compute methods now accept NeighborQuery (or "system-like") objects as the first argument.
* All compute objects that perform neighbor computations now loop over NeighborBond objects.
* Renamed (ref\_points, points) to (points, query\_points) to clarify their usage.
* Bond vector directionality is standardized for all computes that use it (always from query\_point to point).
* Standardized naming of various common parameters across freud such as the search distance r\_max.
* Accumulation is now performed with `compute(..., reset=False)`.
* Arrays returned to Python persist even after the compute object is destroyed or resizes its arrays.
* All class attributes are stored in the C++ members and accessed via getters wrapped as Python properties.
* Code in the freud.common has been moved to freud.util.
* NeighborQuery objects require z == 0 for all points if the box is 2D.
* Renamed several Box methods, box.ParticleBuffer is now locality.PeriodicBuffer.
* Cluster now finds connected components of the neighbor graph (the cluster cutoff distance is given through query arguments).
* Refactored and renamed attributes of Cluster and ClusterProperties modules.
* CorrelationFunction of complex inputs performs the necessary conjugation of the values before computing.
* Updated GaussianDensity constructor to accept tuples as width instead of having 2 distinct signatures.
* RDF bin centers are now strictly at the center of bins.
* RDF no longer performs parallel accumulation of cumulative counts (provided no performance gains and was substantially more complex code).
* MatchEnv has been split into separate classes for the different types of computations it is capable of performing, and these classes all use v2.0-style APIs.
* The Voronoi class was rewritten to use voro++ for vastly improved performance and correctness in edge cases.
* Improved Voronoi plotting code.
* Cubatic uses standard library random functions instead of Saru (which has been removed from the repo).
* APIs for several order parameters have been standardized.
* SolidLiquid order parameter has been completely rewritten, fixing several bugs and simplifying its C++ code.
* Steinhardt uses query arguments.
* PMFTXY2D has been renamed to PMFTXY.
* Removed unused orientations from PMFTXYZ and PMFTXY.
* PMFTXY and PMFTXYZ include the phase space volume of coordinates that are implicitly integrated out (one angle in PMFTXY, and three angles in PMFTXYZ).
* Documentation uses automodule instead of autoclass.
* Citations are now included using bibtex and sphinxcontrib-bibtex.

### Fixed
* Removed all neighbor exclusion logic from all classes, depends entirely on locality module now.
* Compute classes requiring 2D systems check the dimensionality of their input boxes.
* LinkCell nearest neighbor queries properly check the largest distance found before proceeding to next shell.
* LocalDensity uses the correct number of points/query points.
* RDF no longer forces the first bin of the PCF and first two bins of the cumulative counts to be 0.
* Steinhardt uses the ThreadStorage class and properly resets memory where needed.

### Removed
* The freud.util module.
* Python 2 is no longer supported. Python 3.5+ is required.
* LinkCell no longer exposes the internals of the cell list data structure.
* Cubatic no longer returns the per-particle tensor or the constant r4 tensor.

## v1.2.2 - 2019-08-15

### Changed
* LocalWl return values are real instead of complex.

### Fixed
* Fixed missing Condon-Shortley phase affecting LocalWl and Steinhardt Wl
  computations. This missing factor of -1 caused results for third-order (Wl)
  Steinhardt order parameters to be incorrect, shown by their lack of
  rotational invariance. This problem was introduced in v0.5.0.
* Reduced various compiler warnings.
* Possible out of bounds LinkCell access.
* RDF plots now use the provided `ax` object.

## v1.2.1 - 2019-07-26

### Changed
* Optimized performance for `RotationalAutocorrelation`.
* Added new tests for cases with two different sets of points.

### Fixed
* Fixed bug resulting in the `LocalQlNear` and `LocalWlNear` class wrongly
  using a hard instead of a soft cut-off, which may have resulted in an
  incorrect number of neighbors. This would cause incorrect results especially
  for systems with an average n-th nearest-neighbor distance smaller than
  `rmax`. This problem was introduced in v0.6.4.
* Fixed duplicate neighbors found by `LinkCell` `NeighborQuery` methods
* Corrected data in `LocalQl`, `LocalWl` documentation example
* Repeated Cubatic Order Parameter computations use the correct number of
  replicates.
* Repeated calls to `LocalQl.computeNorm` properly reset the underlying data.
* Clarified documentation for `LocalBondProjection` and `MSD`

## v1.2.0 - 2019-06-27

### Added
* Added `.plot()` method and IPython/Jupyter PNG representations for many
  classes.
* `AttributeError` is raised when one tries to access an attribute that has not
  yet been computed.
* Added `freud.parallel.getNumThreads()` method.
* New examples for integration with simulation and visualization workflows.

### Changed
* Removed extra C++ includes to speed up builds.
* The C++ style is now based on clang-format.
* Refactored C++ handling of thread-local storage.
* SolLiq order parameter computations are parallelized with TBB.
* Optimized performance of Voronoi.
* Several Box properties are now given as NumPy arrays instead of tuples.
* Box methods handling multiple vectors are parallelized with TBB.
* Eigen is now used for all matrix diagonalizations.

### Fixed
* Calling setNumThreads works correctly even if a parallel compute method has
  already been called.
* Fixed segfault with chained calls to NeighborQuery API.
* Correct `exclude_ii` logic.

### Removed
* Removed outdated `computeNList` function from `LocalDescriptors`.

## v1.1.0 - 2019-05-23

### Added
* New neighbor querying API to enable reuse of query data structures (see NeighborQuery class).
* AABBQuery (AABB tree-based neighbor finding) added to public API.
* Ability to dynamically select query method based on struct of arguments.
* All compute objects have `__repr__` and `__str__` methods defined.
* NeighborLists can be accessed as arrays of particle indices via
  `__getitem__`.
* ParticleBuffer supports different buffer sizes in x, y, z.
* Box makeCoordinates, makeFraction, getImage now support 2D arrays with
  multiple points.

### Changed
* Use constant memoryviews to prevent errors with read-only inputs.
* LocalQl is now parallelized with TBB.
* Optimized performance of RotationalAutocorrelation.
* NematicOrderParameter uses SelfAdjointEigenSolver for improved stability.
* Added build flags for Cython debugging.
* LinkCell computes cell neighbors on-demand and caches the results for
  significant speedup.

### Fixed
* Corrected type of `y_max` argument to PMFTXY2D from int to float.
* Reduce logging verbosity about array conversion.
* Fixed number of threads set upon exiting the NumThreads context manager.
* Corrected quaternion array sizes and added missing defaults in the
  documentation.
* Empty ParticleBuffers return valid array shapes for concatenation.
* Wheels are built against NumPy 1.10 for improved backwards compatibility.

## v1.0.0 - 2019-02-08

### Added
* Freshly updated README and documentation homepage.
* Moved to [GitHub](https://github.com/glotzerlab/freud).
* New msd.MSD class for computing mean-squared displacements.
* New order.RotationalAutocorrelation class.
* Cython memoryviews are now used to convert between C++ and Cython.
* New and improved freud logo.
* Internal-only AABB tree (faster for many large systems, but API is unstable).

### Changed
* Improved module documentation, especially for PMFT.
* Refactored internals of LocalQl and related classes.
* Upgraded ReadTheDocs configuration.

### Fixed
* Improved CubaticOrderParameter handling of unusable seeds.
* Fixed box error in NearestNeighbors.

### Removed
* All long-deprecated methods and classes were removed.
* Bond module removed.

## v0.11.4 - 2018-11-09

### Added
* Builds are now tested on Windows via Appveyor, though officially unsupported.

### Fixed
* Multiple user-reported issues in setup.py were resolved.
* C++ errors are handled more cleanly as Python exceptions.
* Fixed bug in SolLiq box parameters.
* Documentation corrected for NeighborList.
* Various minor compiler errors on Windows were resolved.

## v0.11.3 - 2018-10-18

### Fixed
* Linux wheels are now pushed to the real PyPI server instead of the test
  server.
* macOS deployment pyenv requires patch versions to be specified.

## v0.11.2 - 2018-10-18

### Fixed
* Error in Python versions in macOS automatic deployment.

## v0.11.1 - 2018-10-18

### Added
* PyPI builds automatically deploy for Mac and Linux.

### Changed
* macOS deployment target is now 10.12 instead of 10.9 to ensure TBB
  compatibility.
* Unwrapping positions with images is now vectorized.
* Minor documentation fixes.

### Fixed
* TBB includes were not always detected correctly by setup.py.

## v0.11.0 - 2018-09-27

### Added
* Example notebooks are now shown in the documentation.
* Many unit tests were added.
* New class: `freud.environment.LocalBondProjection`.
* `freud` is now available on the Python Package Index (PyPI) as
  `freud-analysis`.

### Changed
* Documentation was revised for several modules.
* New class `freud.box.ParticleBuffer` was adapted from the previous
  `VoronoiBuffer` to include support for triclinic boxes.
* The `bond` and `pmft` modules verify system dimensionality matches the
  coordinate system used.
* Minor optimization: arrays are reduced across threads only when necessary.

### Fixed
* NumPy arrays of lengths 2, 3, 6 are now correctly ducktyped into boxes.
* Removed internal use of deprecated code.
* C++ code using `uint` has been changed to `unsigned int`, to improve compiler
  compatibility.

### Deprecated
* In `freud.locality.LinkCell`, `computeCellList()` has been replaced by
  `compute()`.

### Removed
* The `kspace` module has been removed.

## v0.10.0 - 2018-08-27

### Added
* Codecov to track test coverage.
* Properties were added to MatchEnv, AngularSeparation, Cubatic/Nematic order
  parameters, Voronoi.

### Changed
* freud uses Cython and setup.py instead of CMake for installation.
* Properties (not get functions) are the official way to access computed
  results.
* Interface module has been improved significantly.
* density.FloatCF, density.ComplexCF, order parameter documentation is
  improved.
* Many compute methods now use points, orientations from ref\_points,
  ref\_orientations if not provided.
* Reset methods have been renamed to `reset`.

### Fixed
* `kspace` module had a missing factor of pi in the volume calculation of
  `FTsphere`.

### Deprecated
* Get functions have been deprecated.
* Setter methods have been deprecated.
* Reduce methods are called internally, so the user-facing methods have been
  deprecated.

### Removed
* GaussianDensity.resetDensity() is called internally.

## v0.9.0 - 2018-07-30

### Added
* Allow specification of rmin for LocalWl (previously was only possible for
  LocalQl).
* New environment module. Contains classes split from the order module.
* Box duck-typing: methods accepting a box argument will convert box-like
  objects into freud.box.Box objects.
* All Python/Cython code is now validated with flake8 during continuous
  integration.

### Changed
* Refactoring of LocalQl and LocalWl Steinhardt order parameters.
* MatchEnv uses BiMap instead of boost::bimap.
* All boost shared\_arrays have been replaced with std::shared\_ptr.
* Replaced boost geometry with standard containers in brute force registration
  code.
* NearestNeighbors automatically uses ref\_points as the points if points are
  not provided.
* Box::unwrap and Box::wrap return the vectors after updating.
* Everything other than true order parameters moved from Order module to
  Environment module.
* Use lambda function in parallel\_for in CorrelationFunction.
* Tests no longer depend on nose. Python's unittest is used instead.
* Vastly improved documentation clarity and correctness across all modules.
* Docstrings are now in Google format. The developer guide offers guidance for
  module authors.

### Fixed
* Fixed LocalDescriptors producing NaN's in some cases.
* Fixed cython passing C++ the default argument force\_resize to
  NeighborList::resize.
* Standardize freud.common.convert\_array error message.

### Removed
* Boost is no longer needed to build or run freud.
* Removed undocumented shapesplit module.
* Removed extra argument from TransOrderParam in C++.

## v0.8.2 - 2018-06-07

### Added
* Allow specification of maximum number of neighbors to use when computing
  LocalDescriptors

### Changed
* Using the default neighbor list with LocalDescriptors requires specifying the
  precompute argument
* Updated and improved tests
* Cleaned AngularSeparation module and documentation

## v0.8.1 - 2018-05-09

### Fixed
* Memory issue in nlist resolved

## v0.8.0 - 2018-04-06

### Added
* Voronoi neighborlist now includes periodic neighbors
* Voronoi neighborlist computes weight according to the facet area in 3D
* Box module exposes `getImage(vec)`
* Voronoi module can compute and return cell volumes/areas

### Changed
* Cluster module supports box argument in compute methods.
* Refactored C++ code to reduce extraneous #includes
* Refactored PMFT code
* Refactored box module to remove unused methods
* Resolved bug in `kspace.AnalyzeSFactor3D`

### Deprecated
* Box module `getCoordinates()` in favor of duplicate `box.makeCoordinates()`

### Removed
* Removed deprecated API for ComplexWRDF and FloatWRDF

## v0.7.0 - 2018-03-02

### Added
* Added nematic order parameter
* Added optional rmin argument to density.RDF
* Added credits file
* Wrote development guide
* Added Python interface for box periodicity

### Changed
* Various bug fixes and code cleaning
* Fixed all compile-time warnings
* Ensured PEP 8 compliance everywhere
* Minimized boost dependence
* Many documentation rewrites
* Wrote development guide
* Made tests deterministic (seeded RNGs)
* Removed deprecated Box API warnings
* Standardized numpy usage

## v0.6.4 - 2018-02-05

* Added a generic neighbor list interface
* Set up CircleCI for continuous integration
* Set up documentation on ReadTheDocs
* Added bumpversion support
* Various bug fixes
* Added python-style properties for accessing data
* Fixed issues with voronoi neighbor list

## v0.6.0

* trajectory module removed
* box constructor API updated
* PMFTXYZ API updated to take in more quaternions for `face_orientations`, or
  have a sensible default value
* NearestNeighbors:
    - over-expanding box fixed
    - strict rmax mode added
    - ability to get wrapped vectors added
    - minor updates to C-API to return full lists from C
* Addition of Bonding modules
* Addition of local environment matching

## v0.5.0

* Replace boost::shared\_array with std::shared\_ptr (C++ 11)
* Moved all tbb template classes to lambda expressions
* Moved trajectory.Box to box.Box
* trajectory is deprecated
* Fixed Bond Order Diagram and allow for global, local, or orientation
  correlation
* Added python-level voronoi calculation
* Fixed issues with compiling on OS X, including against conda python installs
* Added code to compute bonds between particles in various coordinate systems

## v0.4.1
* PMFT: Fixed issue involving binning of angles correctly
* PMFT: Fixed issue in R12 which prevented compute/accumulate from being called
  with non-flattened arrays
* PMFT: Updated xyz api to allow simpler symmetric orientations to be supplied
* PMFT: Updated pmftXY2D api
* PMFT: Histograms are properly normalized, allowing for comparison between
  systems without needing to "zero" the system
* fsph: Added library to calculate spherical harmonics via cython
* Local Descriptors: Uses fsph, updates to API
* Parallel: Added default behavior to setNumThreads and added context manager


## v0.4.0

* Add compiler flags for C++11 features
* Added Saru RNG (specifically for Cubatic Order Parameter, available to all)
* Cubatic Order Parameter
* Rank 4 tensor struct
* Environment Matching/Cluster Environment
* Shape aware fourier transform for structure factor calculation
* Added deprecation warnings; use python -W once to check warnings
* Added Change Log
* Moved all documentation in Sphinx; documentation improvements
* C++ wrapping moved from Boost to Cython
* Itercell only available in python 3.x
* PMFT:
    - XY
    - XYZ with symmetry, API updated
    - R, T1, T2
    - X, Y, T2
* viz removed (is not compatible with cython)

## No change logs prior to v0.4.0
<!-- Provide a general summary of your changes in the Title above. -->

## Description
<!-- Describe your changes in detail. -->

## Motivation and Context
<!-- Why is this change required? What problem does it solve? -->
<!-- Replace ??? with the issue number that this pull request resolves, if applicable. -->
Resolves: #???

## How Has This Been Tested?
<!-- Please describe in detail how you tested your changes. -->
<!-- Include details of your testing environment, and the tests you ran to
     see how your changes affect other areas of the code, etc. -->

## Screenshots
<!-- (if appropriate) -->

## Types of changes
<!-- What types of changes does your code introduce? Put an `x` in all the boxes that apply: -->
- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds or improves functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to change)
- [ ] Documentation improvement (updates to user guides, docstrings, or developer docs)

## Checklist:
<!-- Put an `x` in all the boxes that apply. If you're unsure about any of
     these, don't hesitate to ask. We're here to help! -->
- [ ] I have read the [**CONTRIBUTING**](https://github.com/glotzerlab/freud/blob/master/CONTRIBUTING.md) document.
- [ ] My code follows the code style of this project.
- [ ] I have updated the documentation (if relevant).
- [ ] I have added tests that cover my changes (if relevant).
- [ ] All new and existing tests passed.
- [ ] I have updated the [credits](https://github.com/glotzerlab/freud/blob/master/doc/source/reference/credits.rst).
- [ ] I have updated the [Changelog](https://github.com/glotzerlab/freud/blob/master/ChangeLog.md).
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is and what you expected to happen.

**To Reproduce**
Steps to reproduce the behavior:
1. Start with a system with a box and particles like ...
2. Create an analysis method for ...
3. Call the function ...
4. See error

**Error output**
If applicable, add terminal output or screenshots to help explain your problem.

**System configuration (please complete the following information):**
 - OS: [e.g. macOS]
 - Version of Python [e.g. 3.7]
 - Version of freud [e.g. 1.0]

**Additional context**
Add any other context about the problem here.
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
# Lennard-Jones Sample Data

The script `lj.py` generates a sample GSD and DCD file using HOOMD-blue.
This data was generated using HOOMD 2.7.
# Minkowski Structure Metrics Validations

Data files in this folder come from [this repository](https://github.com/ArvDee/Minkowski-structure-metrics-calculator/tree/master/testLattices) with work by Robin van Damme.

## Data Description

> Both the averaged and the non-averaged output Minkowski structure metrics should match the Steinhardt bond order parameters for the simple cubic (SC), face-centered cubic (FCC) and hexagonal close-packed (HCP) lattices, whose values are reported in Table 1 of Ref. [1].
> The BCC will not match the Steinhardt case, as the definition of the neighbours differs between the two.
> Two reference values for the Minkowski structure metrics of the BCC lattice are reported in Ref 2, section IIc.
>
> For each lattice there are 4 files that can be used: either a unit cell or a large lattice, and either in .dat or in .gsd format.
> All these should provide the same results.
>
> [1] Mickel et al. 2013, DOI 10.1063/1.4774084
> [2] Dietz et al. 2016, DOI 10.1103/PhysRevE.94.033207

## License

```
MIT License

Copyright (c) 2019 van Damme, R. (Robin)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
# Steinhardt Average Order Parameters Validation

Data files in this folder come from Robin van Damme and collaborators, with permission. Two reference codes ([this repository](https://github.com/ArvDee/Minkowski-structure-metrics-calculator/) and another implementation labeled "GC") were used to compute Steinhardt order parameters.
# Lock-free parallel disjoint set data structure

This is a small self-contained C++11 implementation of the UNION-FIND data
structure with path compression and union by rank and a few extras It supports
concurrent `find()`, `same()` and `unite()` calls as described in the paper

*Wait-free Parallel Algorithms for the Union-Find Problem*
by Richard J. Anderson and Heather Woll

In addition, this class supports optimistic locking (`try_lock()`/`unlock()`)
of disjoint sets and a *combined* unite+unlock operation for pairs of sets.
=====
freud
=====

|Citing freud|
|PyPI|
|conda-forge|
|ReadTheDocs|
|Binder|
|GitHub-Stars|

.. |Citing freud| image:: https://img.shields.io/badge/cite-freud-informational.svg
   :target: https://freud.readthedocs.io/en/stable/reference/citing.html
.. |PyPI| image:: https://img.shields.io/pypi/v/freud-analysis.svg
   :target: https://pypi.org/project/freud-analysis/
.. |conda-forge| image:: https://img.shields.io/conda/vn/conda-forge/freud.svg
   :target: https://anaconda.org/conda-forge/freud
.. |ReadTheDocs| image:: https://readthedocs.org/projects/freud/badge/?version=latest
   :target: https://freud.readthedocs.io/en/latest/?badge=latest
.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/glotzerlab/freud-examples/master?filepath=index.ipynb
.. |GitHub-Stars| image:: https://img.shields.io/github/stars/glotzerlab/freud.svg?style=social
   :target: https://github.com/glotzerlab/freud

Overview
========

The **freud** Python library provides a simple, flexible, powerful set of tools
for analyzing trajectories obtained from molecular dynamics or Monte Carlo
simulations. High performance, parallelized C++ is used to compute standard
tools such as radial distribution functions, correlation functions, order
parameters, and clusters, as well as original analysis methods including
potentials of mean force and torque (PMFTs) and local environment matching. The
**freud** library supports
`many input formats <https://freud.readthedocs.io/en/stable/topics/datainputs.html>`__
and outputs `NumPy arrays <https://numpy.org/>`__, enabling integration
with the scientific Python ecosystem for many typical materials science
workflows.

Resources
=========

- `Reference Documentation <https://freud.readthedocs.io/>`__: Examples, tutorials, topic guides, and package Python APIs.
- `Installation Guide <https://freud.readthedocs.io/en/stable/gettingstarted/installation.html>`__: Instructions for installing and compiling **freud**.
- `freud-users Google Group <https://groups.google.com/d/forum/freud-users>`__: Ask questions to the **freud** user community.
- `GitHub repository <https://github.com/glotzerlab/freud>`__: Download the **freud** source code.
- `Issue tracker <https://github.com/glotzerlab/freud/issues>`__: Report issues or request features.

Related Tools
=============

- `HOOMD-blue <https://hoomd-blue.readthedocs.io/>`__: Perform MD / MC simulations that can be analyzed with **freud**.
- `signac <https://signac.io/>`__: Manage your workflow with **signac**.

Citation
========

When using **freud** to process data for publication, please `use this citation
<https://freud.readthedocs.io/en/stable/reference/citing.html>`__.


Installation
============

The easiest ways to install **freud** are using pip:

.. code:: bash

   pip install freud-analysis

or conda:

.. code:: bash

   conda install -c conda-forge freud

**freud** is also available via containers for `Docker
<https://hub.docker.com/r/glotzerlab/software>`__ or `Singularity
<https://glotzerlab.engin.umich.edu/downloads/glotzerlab>`__.  If you need more detailed
information or wish to install **freud** from source, please refer to the
`Installation Guide
<https://freud.readthedocs.io/en/stable/gettingstarted/installation.html>`__ to compile
**freud** from source.


Examples
========

The **freud** library is called using Python scripts. Many core features are
`demonstrated in the freud documentation
<https://freud.readthedocs.io/en/stable/examples.html>`_. The examples come in
the form of Jupyter notebooks, which can also be downloaded from the `freud
examples repository <https://github.com/glotzerlab/freud-examples>`_ or
`launched interactively on Binder
<https://mybinder.org/v2/gh/glotzerlab/freud-examples/master?filepath=index.ipynb>`_.
Below is a sample script that computes the radial distribution function for a
simulation run with `HOOMD-blue <https://hoomd-blue.readthedocs.io/>`__ and
saved into a `GSD file <https://gsd.readthedocs.io/>`_.

.. code:: python

   import freud
   import gsd.hoomd

   # Create a freud compute object (RDF is the canonical example)
   rdf = freud.density.RDF(bins=50, r_max=5)

   # Load a GSD trajectory (see docs for other formats)
   traj = gsd.hoomd.open('trajectory.gsd', 'rb')
   for frame in traj:
       rdf.compute(system=frame, reset=False)

   # Get bin centers, RDF data from attributes
   r = rdf.bin_centers
   y = rdf.rdf


Support and Contribution
========================

Please visit our repository on `GitHub <https://github.com/glotzerlab/freud>`__ for the library source code.
Any issues or bugs may be reported at our `issue tracker <https://github.com/glotzerlab/freud/issues>`__, while questions and discussion can be directed to our `user forum <https://groups.google.com/forum/#!forum/freud-users>`__.
All contributions to **freud** are welcomed via pull requests!
.. include:: ../../README.rst

Table of Contents
=================

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   gettingstarted/introduction
   gettingstarted/installation
   gettingstarted/quickstart
   gettingstarted/tutorial
   gettingstarted/examples

.. toctree::
   :maxdepth: 2
   :caption: Topic Guides

   topics/querying
   topics/optimizing
   topics/datainputs

.. toctree::
   :maxdepth: 2
   :caption: API

   modules/box
   modules/cluster
   modules/data
   modules/density
   modules/diffraction
   modules/environment
   modules/interface
   modules/locality
   modules/msd
   modules/order
   modules/parallel
   modules/pmft

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/development
   reference/citing
   reference/zreferences
   reference/license
   reference/credits


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. _optimizing:

===========================
Using **freud** Efficiently
===========================

The **freud** library is designed to be both fast and easy-to-use.
In many cases, the library's performance is good enough that users don't need to worry about their usage patterns.
However, in highly performance-critical applications (such as real-time visualization or on-the-fly calculations mid-simulation), uses can benefit from knowing the best ways to make use of **freud**.
This page provides some guidance on this topic.

Reusing Locality Information
============================

Perhaps the most powerful method users have at their disposal for speeding up calculations is proper reuse of the data structures in :mod:`freud.locality`.
As one example, consider using **freud** to calculate multiple neighbor-based quantities for the same set of data points.
It is important to recognize that internally, each time such a calculation is performed using a ``(box, points)`` :class:`tuple`, the compute class is internally rebuilding a neighbor-finding accelerator such a :class:`freud.locality.AABBQuery` object and then using it to find neighbors:

.. code-block:: python

    # Behind the scenes, freud is essentially running
    # freud.locality.AABBQuery(box, points).query(points, dict(r_max=5, exclude_ii=True))
    # and feeding the result to the RDF calculation.
    rdf = freud.density.RDF(bins=50, r_max=5)
    rdf.compute(system=(box, points))


If users anticipate performing many such calculations on the same system of points, they can amortize the cost of rebuilding the :class:`AABBQuery <freud.locality.AABBQuery>` object by constructing it once and then passing it into multiple computations:

.. code-block:: python

    # Now, let's instead reuse the object for a pair of calculations:
    nq = freud.locality.AABBQuery(box=box, points=points)
    rdf = freud.density.RDF(bins=50, r_max=5)
    rdf.compute(system=nq)

    r_max = 4
    orientations = np.array([[1, 0, 0, 0]] * num_points)
    pmft = freud.pmft.PMFTXYZ(r_max, r_max, r_max, bins=100)
    pmft.compute(system=nq, orientations=orientations)

This reuse can significantly improve performance in e.g. visualization contexts where users may wish to calculate a :class:`bond order diagram <freud.environment.BondOrder>` and an :class:`RDF <freud.density.RDF>` at each frame, perhaps for integration with a visualization toolkit like `OVITO <https://www.ovito.org/>`_.

A slightly different use-case would be the calculation of multiple quantities based on *exactly the same set of neighbors*.
If the user in fact expects to perform computations with the exact same pairs of neighbors (for example, to compute :class:`freud.order.Steinhardt` for multiple :math:`l` values), then the user can further speed up the calculation by precomputing the entire :class:`freud.NeighborList` and storing it for future use.

.. code-block:: python

    r_max = 3
    nq = freud.locality.AABBQuery(box=box, points=points)
    nlist = nq.query(points, dict(r_max=r_max)).toNeighborList()
    q6_arrays = []
    for l in range(3, 6):
        ql = freud.order.Steinhardt(l=l)
        q6_arrays.append(ql.compute((box, points), neighbors=nlist).particle_order)


Notably, if the user calls a compute method with ``compute(system=(box, points))``, unlike in the examples above **freud** **will not construct** a :class:`freud.locality.NeighborQuery` internally because the full set of neighbors is completely specified by the :class:`NeighborList <freud.NeighborList>`.
In all these cases, **freud** does the minimal work possible to find neighbors, so judicious use of these data structures can substantially accelerate your code.

Proper Data Inputs
==================

Minor speedups may also be gained from passing properly structured data to **freud**.
The package was originally designed for analyzing particle simulation trajectories, which are typically stored in single-precision binary formats.
As a result, the **freud** library also operates in single precision and therefore converts all inputs to single-precision.
However, NumPy will typically work in double precision by default, so depending on how data is streamed to **freud**, the package may be performing numerous data copies in order to ensure that all its data is in single-precision.
To avoid this problem, make sure to specify the appropriate data types (:attr:`numpy.float32`) when constructing your NumPy arrays.
.. _querying:

=========
Query API
=========

This page provides a thorough review of how neighbor finding is structured in **freud**.
It assumes knowledge at the level of the :ref:`neighbors` level of the tutorial; if you're not familiar with using the :meth:`query <freud.locality.NeighborQuery.query>` method with query arguments to find neighbors of points, please familiarize yourself with that section of the tutorial.

The central interface for neighbor finding is the :py:class:`freud.locality.NeighborQuery` family of classes, which provide methods for dynamically finding neighbors given a :py:class:`freud.box.Box`.
The :py:class:`freud.locality.NeighborQuery` class defines an abstract interface for neighbor finding that is implemented by its subclasses, namely the :py:class:`freud.locality.LinkCell` and :py:class:`freud.locality.AABBQuery` classes.
These classes represent specific data structures used to accelerate neighbor finding.
These two different methods have different performance characteristics, but in most cases :class:`freud.locality.AABBQuery` performs at least as well as, if not better than, :class:`freud.locality.LinkCell` and is entirely parameter free, so it is the default method of choice used internally in **freud**'s ``PairCompute`` classes.

In general, these data structures operate by constructing them using one set of points, after which they can be queried to efficiently find the neighbors of arbitrary other points using :py:meth:`freud.locality.NeighborQuery.query`.


Query Arguments
===============

The table below describes the set of valid query arguments.

+----------------+-----------------------------------------------------------------------+-----------+---------------------------+---------------------------------------------------------------------+
| Query Argument | Definition                                                            | Data type | Legal Values              | Valid for                                                           |
+================+=======================================================================+===========+===========================+=====================================================================+
| mode           | The type of query to perform (distance cutoff or number of neighbors) | str       | 'none', 'ball', 'nearest' | :class:`freud.locality.AABBQuery`, :class:`freud.locality.LinkCell` |
+----------------+-----------------------------------------------------------------------+-----------+---------------------------+---------------------------------------------------------------------+
| r_max          | Maximum distance to find neighbors                                    | float     | r_max > 0                 | :class:`freud.locality.AABBQuery`, :class:`freud.locality.LinkCell` |
+----------------+-----------------------------------------------------------------------+-----------+---------------------------+---------------------------------------------------------------------+
| r_min          | Minimum distance to find neighbors                                    | float     | 0 <= r_min < r_max        | :class:`freud.locality.AABBQuery`, :class:`freud.locality.LinkCell` |
+----------------+-----------------------------------------------------------------------+-----------+---------------------------+---------------------------------------------------------------------+
| num_neighbors  | Number of neighbors                                                   | int       | num_neighbors > 0         | :class:`freud.locality.AABBQuery`, :class:`freud.locality.LinkCell` |
+----------------+-----------------------------------------------------------------------+-----------+---------------------------+---------------------------------------------------------------------+
| exclude_ii     | Whether or not to include neighbors with the same index in the array  | bool      | True/False                | :class:`freud.locality.AABBQuery`, :class:`freud.locality.LinkCell` |
+----------------+-----------------------------------------------------------------------+-----------+---------------------------+---------------------------------------------------------------------+
| r_guess        | Initial search distance for sequence of ball queries                  | float     | r_guess > 0               | :class:`freud.locality.AABBQuery`                                   |
+----------------+-----------------------------------------------------------------------+-----------+---------------------------+---------------------------------------------------------------------+
| scale          | Scale factor for r_guess when not enough neighbors are found          | float     | scale > 1                 | :class:`freud.locality.AABBQuery`                                   |
+----------------+-----------------------------------------------------------------------+-----------+---------------------------+---------------------------------------------------------------------+

Query Modes
===========

Ball Query (Distance Cutoff)
----------------------------

A ball query finds all particles within a specified radial distance of the provided query points.
This query is executed when ``mode='ball'``.
As described in the table above, this mode can be coupled with filters for a minimum distance (``r_min``) and/or self-exclusion (``exclude_ii``).

Nearest Neighbors Query (Fixed Number of Neighbors)
---------------------------------------------------

A nearest neighbor query (sometimes called :math:`k`-nearest neighbors) finds a desired number of neighbor points for each query point, ordered by distance to the query point.
This query is executed when ``mode='nearest'``.
As described in the table above, this mode can be coupled with filters for a maximum distance (``r_max``), minimum distance (``r_min``), and/or self-exclusion (``exclude_ii``).

Mode Deduction
--------------

The ``mode`` query argument specifies the type of query that is being performed, and it therefore governs how other arguments are interpreted.
In most cases, however, the query mode can be deduced from the set of query arguments.
Specifically, any query with the ``num_neighbors`` key set is assumed to be a query with ``mode='nearest'``.
One of ``num_neighbors`` or ``r_max`` must always be specified to form a valid set of query arguments.
Specifying the ``mode`` key explicitly will ensure that querying behavior is consistent if additional query modes are added to **freud**.


Query Results
=============

Although they don't typically need to be operated on directly, it can be useful to know a little about the objects returned by queries.
The :class:`freud.locality.NeighborQueryResult` stores the ``query_points`` passed to a ``query`` and returns neighbors for them one at a time (like any Python :class:`iterator`).
The primary goal of the result class is to support easy iteration and conversion to more persistent formats.
Since it is an iterator, you can use any typical Python approach to consuming it, including passing it to :class:`list` to build a list of the neighbors.
For a more **freud**-friendly approach, you can use the :meth:`toNeighborList <freud.locality.NeighborQueryResult.toNeighborList>` method to convert the object into a **freud** :class:`freud.locality.NeighborList`.
Under the hood, the underlying C++ classes loop through candidate points and identifying neighbors for each ``query_point``; this is the same process that occurs when ``Compute classes`` employ :class:`NeighborQuery <freud.locality.NeighborQuery>` objects for finding neighbors on-the-fly, but in that case it all happens on the C++ side.


Custom NeighborLists
====================

Thus far, we've mostly discussed :class:`NeighborLists <freud.locality.NeighborList>` as a way to persist neighbor information beyond a single query.
In :ref:`optimizing`, more guidance is provided on how you can use these objects to speed up certain uses of **freud**.
However, these objects are also extremely useful because they provide a *completely customizable* way to specify neighbors to **freud**.
Of particular note here is the :meth:`freud.locality.NeighborList.from_arrays` factory function that allows you to make :class:`freud.locality.NeighborList` objects by directly specifying the ``(i, j)`` pairs that should be in the list.
This kind of explicit construction of the list enables custom analyses that would otherwise be impossible.
For example, consider a molecular dynamics simulation in which particles only interact via extremely short-ranged patches on their surface, and that particles should only be considered bonded if their patches are actually interacting, irrespective of how close together the particles themselves are.
This type of neighbor interaction cannot be captured by any normal querying mode, but could be constructed by the user and then fed to **freud** for downstream analysis.

Nearest Neighbor Asymmetry
==========================

There is one important but easily overlooked detail associated with using query arguments with mode ``'nearest'``.
Consider a simple example of three points on the x-axis located at -1, 0, and 2 (and assume the box is of dimensions :math:`(100, 100, 100)`, i.e. sufficiently large that periodicity plays no role):

.. code-block:: python

    box = [100, 100, 100]
    points = [[-1, 0, 0], [0, 0, 0], [2, 0, 0]]
    query_args = dict(mode='nearest', num_neighbors=1, exclude_ii=True)
    list(freud.locality.AABBQuery(box, points).query(points, query_args))
    # Output: [(0, 1, 1), (1, 0, 1), (2, 1, 2)]

Evidently, the calculation is not symmetric.
This feature of nearest neighbor queries can have unexpected side effects if a ``PairCompute`` is performed using distinct ``points`` and ``query_points`` and the two are interchanged.
In such cases, users should always keep in mind that **freud** promises that every ``query_point`` will end up with ``num_neighbors`` points (assuming no hard cutoff ``r_max`` is imposed and enough points are present in the system).
However, it is possible (and indeed likely) that any given ``point`` will have more or fewer than that many neighbors.
This distinction can be particularly tricky for calculations that depend on vector directionality: **freud** imposes the convention that bond vectors always point from ``query_point`` to ``point``, so users of calculations like PMFTs where directionality is important should keep this in mind.
.. _datainputs:

=====================================
Reading Simulation Data for **freud**
=====================================

The **freud** package is designed for maximum flexibility by making minimal assumptions about its data.
However, users accustomed to the more restrictive patterns of most other tools may find this flexibility confusing.
In particular, knowing how to provide data from specific simulation sources can be a significant source of confusion.
This page is intended to describe how various types of data may be converted into a form suitable for **freud**.

To simplify the examples below, we will assume in all cases that the user wishes to compute a :class:`radial distribution function <freud.density.RDF>` over all frames in the trajectory and that the following code has already been run:

.. code-block:: python

    import freud
    rdf = freud.density.RDF(bins=50, r_max=5)

Native Integrations
===================

The **freud** library offers interoperability with several popular tools for particle simulations, analysis, and visualization.
Below is a list of file formats and tools that are directly supported as "system-like" objects (see :class:`freud.locality.NeighborQuery.from_system`).
Such system-like objects are data containers that store information about a periodic box and particle positions.
Other attributes, such as particle orientations, are not included automatically in the system representation and must be loaded as separate NumPy arrays.

GSD Trajectories
----------------

Using the GSD Python API, GSD files can be easily integrated with **freud** as shown in the :ref:`quickstart`.
This format is natively supported by `HOOMD-blue <https://hoomd-blue.readthedocs.io/>`__.
Note: the GSD format can also be read by :ref:`MDAnalysis <mdanalysisreaders>` and :ref:`garnett <garnetttrajectories>`.
Here, we provide an example that reads data from a GSD file.

.. code-block:: python

    import gsd.hoomd
    traj = gsd.hoomd.open('trajectory.gsd', 'rb')

    for frame in traj:
        rdf.compute(system=frame, reset=False)

.. _mdanalysisreaders:

MDAnalysis Readers
------------------

The `MDAnalysis <https://www.mdanalysis.org/>`__ package can read `many popular trajectory formats <https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#supported-coordinate-formats>`__, including common output formats from CHARMM, NAMD, LAMMPS, Gromacs, Tinker, AMBER, GAMESS, HOOMD-blue, and more.

DCD files are among the most familiar simulation outputs due to their longevity.
Here, we provide an example that reads data from a DCD file.

.. code-block:: python

    import MDAnalysis
    reader = MDAnalysis.coordinates.DCD.DCDReader('trajectory.dcd')

    for frame in reader:
        rdf.compute(system=frame, reset=False)

.. _mdtrajreaders:

MDTraj Readers
--------------

The `MDTraj <http://mdtraj.org/>`__ package can read `many popular trajectory formats <http://mdtraj.org/latest/load_functions.html#format-specific-loading-functions>`__, including common output formats from AMBER, MSMBuilder2, Protein Data Bank files, OpenMM, Tinker, Gromacs, LAMMPS, HOOMD-blue, and more.

To use data read with MDTraj in freud, a system-like object must be manually constructed because it does not have a "frame-like" object containing information about the periodic box and particle positions (both quantities are provided as arrays over the whole trajectory).
Here, we provide an example of how to construct a system:

.. code-block:: python

    import mdtraj
    traj = mdtraj.load_xtc('output/prd.xtc', top='output/prd.gro')

    for system in zip(np.asarray(traj.unitcell_vectors), traj.xyz):
        rdf.compute(system=system, reset=False)

.. _garnetttrajectories:

garnett Trajectories
--------------------

The `garnett <https://garnett.readthedocs.io/>`__ package can read `several trajectory formats <https://garnett.readthedocs.io/en/stable/readerswriters.html#file-formats>`__ that have historically been supported by the HOOMD-blue simulation engine, as well as other common types such as DCD and CIF.
The **garnett** package will auto-detect supported file formats by the file extension.
Here, we provide an example that reads data from a POS file.

.. code-block:: python

    import garnett

    with garnett.read('trajectory.pos') as traj:
        for frame in traj:
            rdf.compute(system=frame, reset=False)

OVITO Modifiers
---------------

The `OVITO Open Visualization Tool <https://www.ovito.org/>`__ supports user-written Python modifiers.
The **freud** package can be installed alongside OVITO to enable user-written `Python script modifiers <https://www.ovito.org/docs/current/particles.modifiers.python_script.php>`_ that leverage analyses from **freud**.
Below is an example modifier that creates a user particle property in the OVITO pipeline for Steinhardt bond order parameters.

.. code-block:: python

    import freud

    def modify(frame, data):
        ql = freud.order.Steinhardt(l=6)
        ql.compute(system=data, neighbors={'num_neighbors': 6})
        data.create_user_particle_property(
            name='ql', data_type=float, data=ql.ql)
        print('Created ql property for {} particles.'.format(data.particles.count))

HOOMD-blue Snapshots
--------------------

`HOOMD-blue <https://hoomd-blue.readthedocs.io/>`__ supports analyzers, callback functions that can perform analysis.
Below is an example demonstrating how to use an anlyzer to log the Steinhardt bond order parameter :math:`q_6` during the simulation run.

.. code-block:: python

    import hoomd
    from hoomd import md
    import freud

    hoomd.context.initialize()

    # Create a 10x10x10 simple cubic lattice of particles with type name A
    system = hoomd.init.create_lattice(
        unitcell=hoomd.lattice.sc(a=2.0, type_name='A'), n=10)

    # Specify Lennard-Jones interactions between particle pairs
    nl = md.nlist.cell()
    lj = md.pair.lj(r_cut=3.0, nlist=nl)
    lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)

    # Integrate at constant temperature
    md.integrate.mode_standard(dt=0.005)
    hoomd.md.integrate.langevin(group=hoomd.group.all(), kT=1.2, seed=4)

    # Create a Steinhardt object to analyze snapshots
    ql = freud.order.Steinhardt(l=6)

    def compute_q6(timestep):
        snap = system.take_snapshot()
        ql.compute(system=snap, neighbors={'num_neighbors': 6})
        return ql.order

    # Register a logger that computes q6 and saves to a file
    ql_logger = hoomd.analyze.log(filename='ql.dat', quantities=['q6'], period=100)
    ql_logger.register_callback('q6', compute_q6)

    # Run for 10,000 time steps
    hoomd.run(10e3)

Reading Text Files
==================

Typically, it is best to use one of the natively supported data readers described above; however it is sometimes necessary to parse trajectory information directly from a text file.
One example of a plain text format is the XYZ file format, which can be generated and used by many tools for particle simulation and analysis, including LAMMPS and VMD.
Note that various readers do exist for XYZ files, including MDAnalysis, but in this example we read the file manually to demonstrate how to read these inputs as plain text.
Though they are easy to parse, XYZ files usually contain no information about the system box, so this must already be known by the user.
Assuming knowledge of the box used in the simulation, a LAMMPS XYZ file could be used as follows:

.. code-block:: python

    N = int(np.genfromtxt('trajectory.xyz', max_rows=1))
    traj = np.genfromtxt(
        'trajectory.xyz', skip_header=2,
        invalid_raise=False)[:, 1:4].reshape(-1, N, 3)
    box = freud.box.Box.cube(L=20)

    for frame_positions in traj:
        rdf.compute(system=(box, frame_positions), reset=False)

The first line is the number of particles, so we read this line and use it to determine how to reshape the contents of the rest of the file into a NumPy array.

Other External Readers
======================

For many trajectory formats, high-quality readers already exist.
However sometimes these readers' data structures must be converted into a format understood by **freud**.
Below, we show an example that converts the MDAnalysis box dimensions from a matrix into a :class:`freud.box.Box`.
Note that :ref:`MDAnalysis inputs <mdanalysisreaders>` are natively supported by **freud** without this extra step.
For other formats not supported by a reader listed above, a similar process can be followed.

.. code-block:: python

    import MDAnalysis
    reader = MDAnalysis.coordinates.DCD.DCDReader('trajectory.dcd')

    for frame in reader:
        box = freud.box.Box.from_matrix(frame.triclinic_dimensions)
        rdf.compute(system=(box, frame.positions), reset=False)
.. _environment:

==================
Environment Module
==================

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.environment.BondOrder
    freud.environment.LocalDescriptors
    freud.environment.EnvironmentCluster
    freud.environment.EnvironmentMotifMatch
    freud.environment.AngularSeparationGlobal
    freud.environment.AngularSeparationNeighbor
    freud.environment.LocalBondProjection

.. rubric:: Details

.. automodule:: freud.environment
    :synopsis: Analyze local particle environments.
    :members:
==========
Box Module
==========

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.box.Box

.. rubric:: Details

.. automodule:: freud.box
    :synopsis: Represents periodic boxes.
    :members:
==============
Cluster Module
==============

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.cluster.Cluster
    freud.cluster.ClusterProperties

.. rubric:: Details

.. automodule:: freud.cluster
    :synopsis: Find clusters of points.
    :members:
==============
Density Module
==============

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.density.CorrelationFunction
    freud.density.GaussianDensity
    freud.density.LocalDensity
    freud.density.RDF
    freud.density.SphereVoxelization

.. rubric:: Details

.. automodule:: freud.density
    :synopsis: Analyze system density.
    :members:
============
Order Module
============

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.order.Cubatic
    freud.order.Nematic
    freud.order.Hexatic
    freud.order.Translational
    freud.order.Steinhardt
    freud.order.SolidLiquid
    freud.order.RotationalAutocorrelation

.. rubric:: Details

.. automodule:: freud.order
    :synopsis: Compute order parameters.
    :members:
==================
Diffraction Module
==================

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.diffraction.DiffractionPattern
    freud.diffraction.StaticStructureFactorDebye
    freud.diffraction.StaticStructureFactorDirect

.. rubric:: Details

.. automodule:: freud.diffraction
    :synopsis: Analyze diffraction patterns and structure factors.
    :members:
==========
MSD Module
==========

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.msd.MSD

.. rubric:: Details

.. automodule:: freud.msd
    :synopsis: Compute mean squared displacement.
    :members:
===============
Locality Module
===============

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.locality.AABBQuery
    freud.locality.LinkCell
    freud.locality.NeighborList
    freud.locality.NeighborQuery
    freud.locality.NeighborQueryResult
    freud.locality.PeriodicBuffer
    freud.locality.Voronoi

.. rubric:: Details

.. automodule:: freud.locality
    :synopsis: Data structures for finding neighboring points.
    :members:
===========
PMFT Module
===========

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.pmft.PMFTR12
    freud.pmft.PMFTXYT
    freud.pmft.PMFTXY
    freud.pmft.PMFTXYZ

.. rubric:: Details

.. automodule:: freud.pmft
    :synopsis: Compute potentials of mean force and torque.
    :members:
===========
Data Module
===========

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.data.UnitCell
    freud.data.make_random_system

.. rubric:: Details

.. automodule:: freud.data
    :synopsis: Utility functions for generating standard types of systems of points.
    :members:
================
Interface Module
================

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.interface.Interface

.. rubric:: Details

.. automodule:: freud.interface
    :synopsis: Measure interfaces.
    :members:
===============
Parallel Module
===============

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    freud.parallel.NumThreads
    freud.parallel.get_num_threads
    freud.parallel.set_num_threads

.. rubric:: Details

.. automodule:: freud.parallel
    :synopsis: Manage TBB thread usage.
    :members:
=====================
How to cite **freud**
=====================

Please acknowledge the use of this software within the body of your publication for example by copying or adapting the following formulation:

*Data analysis for this publication utilized the freud library[1].*

  [1] V. Ramasubramani, B. D. Dice, E. S. Harper, M. P. Spellings, J. A. Anderson, and S. C. Glotzer. freud: A Software Suite for High Throughput Analysis of Particle Simulation Data. Computer Physics Communications Volume 254, September 2020, 107275. doi:10.1016/j.cpc.2020.107275.

The paper is available online from `Computer Physics Communications <https://www.sciencedirect.com/science/article/pii/S0010465520300916>`_ and a pre-print is freely available on `arXiv <https://arxiv.org/abs/1906.06317>`_.

To cite this reference, you can use the following BibTeX entry:

.. code-block:: bibtex

    @article{freud2020,
        title = {freud: A Software Suite for High Throughput
                 Analysis of Particle Simulation Data},
        author = {Vyas Ramasubramani and
                  Bradley D. Dice and
                  Eric S. Harper and
                  Matthew P. Spellings and
                  Joshua A. Anderson and
                  Sharon C. Glotzer},
        journal = {Computer Physics Communications},
        volume = {254},
        pages = {107275},
        year = {2020},
        issn = {0010-4655},
        doi = {https://doi.org/10.1016/j.cpc.2020.107275},
        url = {http://www.sciencedirect.com/science/article/pii/S0010465520300916},
        keywords = {Simulation analysis, Molecular dynamics,
                    Monte Carlo, Computational materials science},
    }

Optionally, publications using **freud** in the context of machine learning or data visualization may also wish to cite this reference.

  [2] B. D. Dice, V. Ramasubramani, E. S. Harper, M. P. Spellings, J. A. Anderson, and S. C. Glotzer. Analyzing Particle Systems for Machine Learning and Data Visualization with freud. Proceedings of the 18th Python in Science Conference, 2019, 27-33. doi:10.25080/Majora-7ddc1dd1-004.

The paper is freely available from the `SciPy Conference website <http://conference.scipy.org/proceedings/scipy2019/bradley_dice.html>`_.

To cite this reference, you can use the following BibTeX entry:

.. code-block:: bibtex

    @InProceedings{freud2019,
        title = {Analyzing Particle Systems for Machine Learning
                 and Data Visualization with freud},
        author = {Bradley D. Dice and
                  Vyas Ramasubramani and
                  Eric S. Harper and
                  Matthew P. Spellings and
                  Joshua A. Anderson and
                  Sharon C. Glotzer },
        booktitle = {Proceedings of the 18th Python in Science Conference},
        pages = {27-33},
        year = {2019},
        editor = {Chris Calloway and David Lippa and Dillon Niederhut and David Shupe},
        doi = {https://doi.org/10.25080/Majora-7ddc1dd1-004},
        url = {http://conference.scipy.org/proceedings/scipy2019/bradley_dice.html}
    }
Credits
=======

freud Developers
----------------

The following people contributed to the development of freud.

Vyas Ramasubramani - **Lead developer**

* Ensured PEP8 compliance.
* Added CircleCI continuous integration support.
* Create environment module and refactored order module.
* Rewrote most of freud docs, including order, density, and environment modules.
* Fixed nematic order parameter.
* Add properties for accessing class members.
* Various minor bug fixes.
* Refactored PMFT code.
* Refactored Steinhardt order parameter code.
* Wrote numerous examples of freud usage.
* Rewrote most of freud tests.
* Replaced CMake-based installation with setup.py using Cython.
* Add code coverage metrics.
* Added support for installing from PyPI, including ensuring that NumPy is installed.
* Converted all docstrings to Google format, fixed various incorrect docs.
* Debugged and added rotational autocorrelation code.
* Added MSD module.
* Wrote NeighborQuery, _QueryArgs, NeighborQueryResult classes.
* Wrote neighbor iterator infrastructure.
* Wrote PairCompute and SpatialHistogram parent classes.
* Wrote ManagedArray class.
* Wrote C++ histogram-related classes.
* Initial design of freud 2.0 API (NeighborQuery objects, neighbor computations, histograms).
* Standardized neighbor API in Python to use dictionaries of arguments or NeighborList objects for all pair computations.
* Standardized all attribute access into C++ with Python properties.
* Standardized variable naming of points/query\_points across all of freud.
* Standardized vector directionality in computes.
* Enabled usage of quaternions in place of angles for orientations in 2D PMFT calculations.
* Wrote new freud 2.0 compute APIs based on neighbor\_query objects and neighbors as either dictionaries or NeighborLists.
* Rewrote MatchEnv code to fit freud 2.0 API, splitting it into 3 separate calculations and rewriting internals using NeighborQuery objects.
* Wrote tutorial and reference sections of documentation.
* Unified util and common packages.
* Rewrote all docstrings in the package for freud 2.0.
* Changed Cubatic to use Mersenne Twisters for rng.
* Moved all citations into Bibtex format.
* Created data module.
* Standardized PMFT normalization.
* Enabled optional normalization of RDF.
* Changed correlation function to properly take the complex conjugate of inputs.
* Wrote developer documentation for version 2.0.
* Fixed handling of 2D systems from various data sources.
* Fixed usage of query orientations in PMFTXY, PMFTXYT and PMFTXYZ when points and query points are not identical.
* Refactored and standardized PMFT tests.
* Rewrote build system to use scikit-build.
* Added support for pre-commit hooks.
* Added the `out` option for the `unwrap`, `make_fractional`, and `make_absolute` methods of `Box`.
* Enabled access to the qlmi arrays in Steinhardt and SolidLiquid and added rigorous tests of correctness.

Bradley Dice - **Lead developer**

* Cleaned up various docstrings.
* Fixed bugs in HexOrderParameter.
* Cleaned up testing code.
* Added bumpversion support.
* Reduced all compile warnings.
* Added Python interface for box periodicity.
* Added Voronoi support for neighbor lists across periodic boundaries.
* Added Voronoi weights for 3D.
* Added Voronoi cell volume computation.
* Incorporated internal BiMap class for Boost removal.
* Wrote numerous examples of freud usage.
* Added some freud tests.
* Added ReadTheDocs support.
* Rewrote interface module into pure Cython.
* Added box duck-typing.
* Removed nose from unit testing.
* Use lambda function for parallelizing CorrelationFunction with TBB.
* Finalized boost removal.
* Wrote AABBQuery class.
* Consolidated cluster module functionality.
* Rewrote SolidLiquid order parameter class.
* Updated AngularSeparation class.
* Rewrote Voronoi implementation to leverage voro++.
* Implemented Voronoi bond weighting to enable Minkowski structure metrics.
* Refactored methods in Box and PeriodicBuffer for v2.0.
* Added checks to C++ for 2D boxes where required.
* Refactored cluster module.
* Standardized vector directionality in computes.
* NeighborQuery support to ClusterProperties, GaussianDensity, Voronoi, PeriodicBuffer, Interface.
* Standardized APIs for order parameters.
* Added radius of gyration to ClusterProperties.
* Improved Voronoi plotting code.
* Corrected number of points/query points in LocalDensity.
* Made PeriodicBuffer inherit from _Compute.
* Removed cudacpu and HOOMDMath includes.
* Added plotting functionality for Box and NeighborQuery objects.
* Added support for reading system data directly from MDAnalysis, garnett, gsd, HOOMD-blue, and OVITO.
* Revised tutorials and documentation on data inputs.
* Updated MSD to perform accumulation with ``compute(..., reset=False)``.
* Added test PyPI support to continuous integration.
* Added continuous integration to freud-examples.
* Implemented periodic center of mass computations in C++.
* Revised docs about query modes.
* Implemented smarter heuristics in Voronoi for voro++ block sizes, resulting in significant performance gains for large systems.
* Corrected calculation of neighbor distances in the Voronoi NeighborList.
* Added finite tolerance to ensure stability of 2D Voronoi NeighborList computations.
* Improved stability of Histogram bin calculations.
* Improved error handling of Cubatic input parameters.
* Added 2D Minkowski Structure Metrics to Hexatic, enabled by using ``weighted=True`` along with a Voronoi NeighborList.
* Worked with Tommy Waltmann to add the SphereVoxelization feature.
* Fixed GaussianDensity normalization in 2D systems.
* Prevented GaussianDensity from computing 3D systems after it has computed 2D systems.
* Contributed code, design, and testing for ``DiffractionPattern`` class.
* Fixed ``Hexatic`` order parameter (unweighted) to normalize by number of neighbors instead of the symmetry order k.
* Added ``num_query_points`` and ``num_points`` attributes to NeighborList class.
* Added scikit-build support for Windows.
* Fixed 2D image calculations.
* Optimized NeighborList ``filter`` method.
* Fixed documented formulas for ``Steinhardt`` class.
* Fixed incorrect computation of ``Steinhardt`` averaged quantities.
* Fixed RPATH problems affecting ``libfreud.so`` in Linux wheels.
* Updated lambda functions to capture ``this`` by reference, to ensure compatibility with C++20 and above.
* Contributed code, design, documentation, and testing for ``StaticStructureFactorDebye`` class.
* Fixed ``Box.contains`` to run in linear time, ``O(num_points)``.
* Contributed code, design, documentation, and testing for ``StaticStructureFactorDirect`` class.
* Fixed doctests to run with pytest.

Eric Harper, University of Michigan - **Former lead developer**

* Added TBB parallelism.
* Wrote PMFT module.
* Added NearestNeighbors (since removed).
* Wrote RDF.
* Added bonding module (since removed).
* Added cubatic order parameter.
* Added hexatic order parameter.
* Added Pairing2D (since removed).
* Created common array conversion logic.

Joshua A. Anderson, University of Michigan - **Creator and former lead developer**

* Initial design and implementation.
* Wrote LinkCell and IteratorLinkCell.
* Wrote GaussianDensity, LocalDensity.
* Added parallel module.
* Added indexing modules (since removed).
* Wrote Cluster and ClusterProperties modules.

Matthew Spellings - **Former lead developer**

* Added generic neighbor list.
* Enabled neighbor list usage across freud modules.
* Added correlation functions.
* Added LocalDescriptors class.
* Added interface module.

Erin Teich

* Wrote environment matching (MatchEnv) class.
* Wrote BondOrder class (with Julia Dshemuchadse).
* Wrote AngularSeparation class (with Andrew Karas).
* Contributed to LocalQl development.
* Wrote LocalBondProjection class.

M. Eric Irrgang

* Authored kspace module (since removed).
* Fixed numerous bugs.
* Contributed to freud.shape (since removed).

Chrisy Du

* Authored Steinhardt order parameter classes.
* Fixed support for triclinic boxes.

Antonio Osorio

* Developed TrajectoryXML class.
* Various bug fixes.
* OpenMP support.

Richmond Newman

* Developed the freud box.
* Solid liquid order parameter.

Carl Simon Adorf

* Developed the Python box module.

Jens Glaser

* Wrote kspace front-end (since removed).
* Modified kspace module (since removed).
* Wrote Nematic order parameter class.

Benjamin Schultz

* Wrote Voronoi class.
* Fix normalization in GaussianDensity.
* Bug fixes in shape module (since removed).

Bryan VanSaders

* Make Cython catch C++ exceptions.
* Add shiftvec option to PMFT.

Ryan Marson

* Various GaussianDensity bugfixes.

Yina Geng

* Co-wrote Voronoi neighbor list module.
* Add properties for accessing class members.

Carolyn Phillips

* Initial design and implementation.
* Package name.

Ben Swerdlow

* Documentation and installation improvements.

James Antonaglia

* Added number of neighbors as an argument to HexOrderParameter.
* Bugfixes.
* Analysis of deprecated kspace module.

Mayank Agrawal

* Co-wrote Voronoi neighbor list module.

William Zygmunt

* Helped with Boost removal.

Greg van Anders

* Bugfixes for CMake and SSE2 installation instructions.

James Proctor

* Cythonization of the cluster module.

Rose Cersonsky

* Enabled TBB-parallelism in density module.
* Fixed how C++ arrays were pulled into Cython.

Wenbo Shen

* Translational order parameter.

Andrew Karas

* Angular separation.
* Wrote reference implementation for rotational autocorrelation.

Paul Dodd

* Fixed CorrelationFunction namespace, added ComputeOCF class for TBB parallelization.

Tim Moore

* Added optional rmin argument to density.RDF.
* Enabled NeighborList indexing.
* Documentation fixes.

Alex Dutton

* BiMap class for MatchEnv.

Matthew Palathingal

* Replaced use of boost shared arrays with shared ptr in Cython.
* Helped incorporate BiMap class into MatchEnv.

Kelly Wang

* Enabled NeighborList indexing.
* Added methods ``compute_distances`` and ``compute_all_distances`` to Box.
* Added method ``crop`` to Box.
* Added 2D Box tests for ``get_image`` and ``contains``.
* Added the ``reset`` argument to the ``compute`` method of ``DiffractionPattern`` class.

Yezhi Jin

* Added support for 2D arrays in the Python interface to Box functions.
* Rewrote Voronoi implementation to leverage voro++.
* Implemented Voronoi bond weighting to enable Minkowski structure metrics.
* Contributed code, design, and testing for ``DiffractionPattern`` class.

Brandon Butler

* Added support for multiple ``l`` in ``Steinhardt`` along with performance improvements.
* Rewrote Steinhardt order parameter.

Jin Soo Ihm

* Added benchmarks.
* Contributed to NeighborQuery classes.
* Refactored C++ to perform neighbor queries on-the-fly.
* Added plotting functions to analysis classes.
* Wrote RawPoints class.
* Created Compute parent class with decorators to ensure properties have been computed.
* Updated common array conversion logic.
* Added many validation tests.

Mike Henry

* Fixed syntax in freud-examples notebooks for v2.0.
* Updated documentation links

Michael Stryk

* Added short examples into Cluster, Density, Environment, and Order Modules.

Tommy Waltmann

* Worked with Bradley Dice to add the SphereVoxelization feature.
* Contributed code, design, and testing for ``DiffractionPattern`` class.
* Contributed code, design, and testing for ``StaticStructureFactorDebye`` class.
* Contributed code, design, and testing for ``StaticStructureFactorDirect`` class.
* Refactor tests for ``StaticStructureFactor`` classes.
* Improve CMake build system to use more modern style.
* Remove CI build configurations from CircleCI which were already covered by CIBuildWheel

Maya Martirossyan

* Added test for Steinhardt for particles without neighbors.

Pavel Buslaev

* Added ``values`` argument to compute method of ``GaussianDensity`` class.

Charlotte Zhao

* Worked with Vyas Ramasubramani and Bradley Dice to add the ``out`` option for ``box.Box.wrap``.

Domagoj Fijan

* Contributed code, design, documentation, and testing for ``StaticStructureFactorDebye`` class.
* Contributed code, design, documentation, and testing for ``StaticStructureFactorDirect`` class.

Andrew Kerr

* Contributed documentation for ``StaticStructureFactorDebye`` class.

Emily Siew

* Contributed documentation for ``StaticStructureFactorDebye`` class.

Dylan Marx

* Contributed documentation for ``NeighborQuery`` class.

Source code
-----------

.. highlight:: none

Eigen (http://eigen.tuxfamily.org) is included as a git submodule in freud.
Eigen is made available under the Mozilla Public License v2.0
(http://mozilla.org/MPL/2.0/). Its linear algebra routines are used for
various tasks including the computation of eigenvalues and eigenvectors.

fsph (https://github.com/glotzerlab/fsph) is included as a git submodule in
freud. It is used for the calculation of spherical harmonics. fsph is made
available under the MIT license::

    Copyright (c) 2016 The Regents of the University of Michigan

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

HOOMD-blue (https://github.com/glotzerlab/hoomd-blue) is the original source of
some algorithms and tools for vector math implemented in freud. HOOMD-blue is
made available under the BSD 3-Clause license::

	BSD 3-Clause License for HOOMD-blue

	Copyright (c) 2009-2019 The Regents of the University of Michigan All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	1. Redistributions of source code must retain the above copyright notice,
	   this list of conditions and the following disclaimer.

	2. Redistributions in binary form must reproduce the above copyright notice,
	   this list of conditions and the following disclaimer in the documentation
	   and/or other materials provided with the distribution.

	3. Neither the name of the copyright holder nor the names of its contributors
	   may be used to endorse or promote products derived from this software without
	   specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
	ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

voro++ (https://github.com/chr1shr/voro) is included as a git submodule in
freud. It is used for computing Voronoi diagrams. voro++ is made available
under the following license::

    Voro++ Copyright (c) 2008, The Regents of the University of California, through
    Lawrence Berkeley National Laboratory (subject to receipt of any required
    approvals from the U.S. Dept. of Energy). All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    (1) Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.

    (2) Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

    (3) Neither the name of the University of California, Lawrence Berkeley
    National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    You are under no obligation whatsoever to provide any bug fixes, patches, or
    upgrades to the features, functionality or performance of the source code
    ("Enhancements") to anyone; however, if you choose to make your Enhancements
    available either publicly, or directly to Lawrence Berkeley National
    Laboratory, without imposing a separate written license agreement for such
    Enhancements, then you hereby grant the following license: a non-exclusive,
    royalty-free perpetual license to install, use, modify, prepare derivative
    works, incorporate into other computer software, distribute, and sublicense
    such enhancements or derivative works thereof, in binary and source code form.
==========
References
==========

.. It's necessary to have a file that creates the bibliography that is
   processed by Sphinx AFTER all of the :cite: commands have been parsed. To
   fix this, we create this dummy file that comes as late as possible and
   import the bibliography here.

.. bibliography:: freud.bib
.. _development:

=================
Development Guide
=================

Contributions to **freud** are highly encouraged.
The pages below offer information about how the project is structured, the goals of the **freud** library, and how to contribute new modules.

.. toctree::
   :maxdepth: 2

   development/design
   development/howtoadd
   development/specialtopics
   development/releases
=========================
Making **freud** Releases
=========================

Release Process
===============

The steps to make a new release of **freud** are documented on the GitHub wiki:

https://github.com/glotzerlab/freud/wiki/Creating-a-Release
=================
Design Principles
=================

Vision
======

The **freud** library is designed to be a powerful and flexible library
for the analysis of simulation output. To support a variety of
analysis routines, **freud** places few restrictions on its components.
The primary requirement for an analysis routine in **freud** is that it
should be substantially computationally intensive so as to require
coding up in C++: **all freud code should be composed of fast C++
routines operating on systems of particles in periodic boxes.** To
remain easy-to-use, all C++ modules should be wrapped in Python
code so they can be easily accessed from Python scripts or through
a Python interpreter.

In order to achieve this goal, **freud** takes the following viewpoints:

* **freud** works directly with `NumPy <https://numpy.org/>`__
  arrays to retain maximum flexibility. Integrations with other tools
  should be performed via the common data representations of NumPy arrays.
* For ease of maintenance, **freud** uses Git for version control;
  GitHub for code hosting and issue tracking; and the PEP 8
  standard for code, stressing explicitly written code which is easy
  to read.
* To ensure correctness, **freud** employs unit testing using the
  Python :mod:`pytest` framework. In addition, **freud** utilizes
  `CircleCI <https://circleci.com>`__ for continuous integration to
  ensure that all of its code works correctly and that any changes or
  new features do not break existing functionality.

Language choices
================

The **freud** library is written in two languages: Python and C++. C++ allows for
powerful, fast code execution while Python allows for easy, flexible
use. Intel Threading Building Blocks parallelism provides further power to
C++ code. The C++ code is wrapped with Cython, allowing for user
interaction in Python. NumPy provides the basic data structures in
**freud**, which are commonly used in other Python plotting libraries and
packages.

Unit Tests
==========

All modules should include a set of unit tests which test the correct
behavior of the module. These tests should be simple and short, testing
a single function each, and completing as quickly as possible
(ideally < 10 sec, but times up to a minute are acceptable if justified).

Benchmarks
==========

Modules can be benchmarked in the following way.
The following code is an example benchmark for the :code:`freud.density.RDF` module.

.. literalinclude:: ../../../../benchmarks/benchmark_density_RDF.py
   :language: python
   :linenos:

in a file :code:`benchmark_density_RDF.py` in the :code:`benchmarks` directory.
More examples can be found in the :code:`benchmarks` directory.
The runtime of :code:`BenchmarkDensityRDF.bench_run` will be timed for :code:`number`
of times on the input sizes of :code:`Ns`. Its runtime with respect to the number of threads
will also be measured. Benchmarks are run as a part of continuous integration,
with performance comparisons between the current commit and the master branch.

Make Execution Explicit
=======================

While it is tempting to make your code do things "automatically", such
as have a calculate method find all :code:`_calc` methods in a class, call
them, and add their returns to a dictionary to return to the user, it is
preferred in **freud** to execute code explicitly. This helps avoid issues
with debugging and undocumented behavior:

.. code-block:: python

    # this is bad
    class SomeFreudClass(object):
        def __init__(self, **kwargs):
            for key in kwargs.keys:
                setattr(self, key, kwargs[key])

    # this is good
    class SomeOtherFreudClass(object):
        def __init__(self, x=None, y=None):
            self.x = x
            self.y = y

Code Duplication
================

When possible, code should not be duplicated. However, being explicit is
more important. In **freud** this translates to many of the inner loops of
functions being very similar:

.. code-block:: c++

    // somewhere deep in function_a
    for (int i = 0; i < n; i++)
        {
        vec3[float] pos_i = position[i];
        for (int j = 0; j < n; j++)
            {
            pos_j = = position[j];
            // more calls here
            }
        }

    // somewhere deep in function_b
    for (int i = 0; i < n; i++)
        {
        vec3[float] pos_i = position[i];
        for (int j = 0; j < n; j++)
            {
            pos_j = = position[j];
            // more calls here
            }
        }

While it *might* be possible to figure out a way to create a base C++
class all such classes inherit from, run through positions, call a
calculation, and return, this would be rather complicated. Additionally,
any changes to the internals of the code may result in performance
penalties, difficulty in debugging, etc. As before, being explicit is
better.

However, if you have a class which has a number of methods, each of
which requires the calling of a function, this function should be
written as its own method (instead of being copy-pasted into each
method) as is typical in object-oriented programming.

Python vs. Cython vs. C++
=========================

The **freud** library is meant to leverage the power of C++ code imbued with
parallel processing power from TBB with the ease of writing Python code.
The bulk of your calculations should take place in C++, as shown in the
snippet below:

.. code-block:: python

    # this is bad
    def badHeavyLiftingInPython(positions):
        # check that positions are fine
        for i, pos_i in enumerate(positions):
            for j, pos_j in enumerate(positions):
                if i != j:
                    r_ij = pos_j - pos_i
                    # ...
                    computed_array[i] += some_val
        return computed_array

    # this is good
    def goodHeavyLiftingInCPlusPlus(positions):
        # check that positions are fine
        cplusplus_heavy_function(computed_array, positions, len(pos))
        return computed_array

In the C++ code, implement the heavy lifting function called above from Python:

.. code-block:: c++

    void cplusplus_heavy_function(float* computed_array,
                                  float* positions,
                                  int n)
        {
        for (int i = 0; i < n; i++)
            {
            for (int j = 0; j < n; j++)
                {
                if (i != j)
                    {
                    r_ij = pos_j - pos_i;
                    // ...
                    computed_array[i] += some_val;
                    }
                }
            }
        }

Some functions may be necessary to write at the Python level due to a Python
library not having an equivalent C++ library, complexity of coding, etc. In
this case, the code should be written in Cython and a *reasonable* attempt
to optimize the code should be made.
.. _specialtopics:

==============
Special Topics
==============

Some of the most central components have a high level of abstraction.
This abstraction has multiple advantages: it dramatically simplifies the process of implementing new code, it reduces duplication throughout the code base, and ensures that bug fixes and optimization can occur along a single path for the entire module.
However, this abstraction comes at the cost of significant complexity.
This documentation should help orient developers familiarizing themselves with these topics by providing high-level overviews of how these parts of the code are structured and how the pieces fit together.

.. toctree::
   :maxdepth: 2

   specialtopics/memory
   specialtopics/neighbors
   specialtopics/histograms
=========================
Contributing to **freud**
=========================

Code Conventions
================

Pre-commit
----------

It is strongly recommended to `set up a pre-commit hook <https://pre-commit.com/>`__ to ensure code is compliant with all automated linters and style checks before pushing to the repository:

.. code-block:: bash

    pip install pre-commit
    pre-commit install

To manually run `pre-commit <https://pre-commit.com/>`__ for all the files present in the repository, run the following command:

.. code-block:: bash

    pre-commit run --all-files --show-diff-on-failure


Python
------

Python (and Cython) code in **freud** should follow `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_.

During continuous integration (CI), all Python and Cython code in **freud** is analyzed using automated linters and formatters including :code:`flake8`, :code:`black`, :code:`isort`, and :code:`pyupgrade`.
Documentation is written in reStructuredText and generated using `Sphinx <http://www.sphinx-doc.org/en/stable/index.html>`_.
It should be written according to the `Google Python Style Guide <https://github.com/google/styleguide/blob/gh-pages/pyguide.md#38-comments-and-docstrings>`_.
A few specific notes:

- The shapes of NumPy arrays should be documented as part of the type in the following manner::

    points ((:math:`N_{points}`, 3) :class:`numpy.ndarray`):

- Optional arguments should be documented as such within the type after the actual type, and the default value should be included within the description::

    box (:class:`freud.box.Box`, optional): Simulation box (Default value = None).


C++
---

C++ code should follow the result of running :code:`clang-format` with the style specified in the file :code:`.clang-format`.
Please refer to the `clang-format documentation <https://clang.llvm.org/docs/ClangFormat.html>`__ for details.
The :code:`clang-format` style will be automatically enforced by pre-commit via CI.
When in doubt, run :code:`clang-format -style=file FILE_WITH_YOUR_CODE` in the top directory of the **freud** repository.

The :code:`check-style` step of continuous integration (CI) runs :code:`clang-tidy` and :code:`cppcheck`.
If the :code:`check-style` CI fails, please read the output log for information on what to fix.

Additionally, all CMake code is tested using `cmakelang's cmake-format <https://cmake-format.readthedocs.io/en/latest/index.html>`__.

Doxygen docstrings should be used for classes, functions, etc.


Code Organization
=================

The code in **freud** is a mix of Python, Cython, and C++.
From a user's perspective, methods in **freud** correspond to ``Compute`` classes, which are contained in Python modules that group methods by topic.
To keep modules well-organized, **freud** implements the following structure:

- All C++ code is stored in the ``cpp`` folder at the root of the repository, with subdirectories corresponding to each module (e.g. ``cpp/locality``).
- Python code is stored in the ``freud`` folder at the root of the repository.
- C++ code is exposed to Python using Cython code contained in pxd files with the following convention: ``freud/_MODULENAME.pxd`` (note the preceding underscore).
- The core Cython code for modules is contained in ``freud/MODULENAME.pyx`` (no underscore).
- Generated Cython C++ code (e.g. ``freud/MODULENAME.cxx``) should not be committed during development. These files are generated using Cython when building from source, and are unnecessary when installing compiled binaries.
- If a Cython module contains code that must be imported into other Cython modules (such as the :class:`freud.box.Box` class), the ``pyx`` file must be accompanied by a ``pxd`` file with the same name: ``freud/MODULENAME.pxd`` (distinguished from ``pxd`` files used to expose C++ code by the lack of a preceding underscore). For more information on how ``pxd`` files work, see the `Cython documentation <https://cython.readthedocs.io/en/latest/src/tutorial/pxd_files.html>`_.
- All tests in **freud** are based on the Python :mod:`pytest` library and are contained in the ``tests`` folder. Test files are named by the convention ``tests/test_MODULENAME_CLASSNAME.py``.
- Benchmarks for **freud** are contained in the ``benchmarks`` directory and are named analogously to tests: ``benchmarks/benchmark_MODULENAME_CLASSNAME.py``.

Benchmarks
----------

Benchmarking in **freud** is performed by running the ``benchmarks/benchmarker.py`` script.
This script finds all benchmarks (using the above naming convention) and attempts to run them.
Each benchmark is defined by extending the ``Benchmark`` class defined in ``benchmarks/benchmark.py``, which provides the standard benchmarking utilities used in **freud**.
Subclasses just need to define a few methods to parameterize the benchmark, construct the **freud** object being benchmarked, and then call the relevant compute method.
Rather than describing this process in detail, we consider the benchmark for the :code:`freud.density.RDF` module as an example.

.. literalinclude:: ../../../../benchmarks/benchmark_density_RDF.py
   :language: python
   :linenos:

The ``__init__`` method defines basic parameters of the run, the ``bench_setup`` method is called to build up the :class:`~freud.density.RDF` object, and the ``bench_run`` is used to time and call ``compute``.
More examples can be found in the :code:`benchmarks` directory.
The runtime of :code:`BenchmarkDensityRDF.bench_run` will be timed for :code:`number` of times on the input sizes of :code:`Ns`.
Its runtime with respect to the number of threads will also be measured.
Benchmarks are run as a part of continuous integration, with performance comparisons between the current commit and the master branch.

Steps for Adding New Code
=========================

Once you've determined to add new code to **freud**, the first step is to create a new branch off of :code:`master`.
The process of adding code differs based on whether or not you are editing an existing module in **freud**.
Adding new methods to an existing module in **freud** requires creating the new C++ files in the ``cpp`` directory, modifying the corresponding ``_MODULENAME.pxd`` file in the ``freud`` directory, and creating a wrapper class in ``freud/MODULENAME.pyx``.
If the new methods belong in a new module, you must create the corresponding ``cpp`` directory and the ``pxd`` and ``pyx`` files accordingly.

In order for code to compile, it must be added to the relevant ``CMakeLists.txt`` file.
New C++ files for existing modules must be added to the corresponding ``cpp/MODULENAME/CMakeLists.txt`` file.
For new modules, a ``cpp/NEWMODULENAME/CMakeLists.txt`` file must be created, and in addition the new module must be added to the ``cpp/CMakeLists.txt`` file in the form of both an ``add_subdirectory`` command and addition to the ``libfreud`` library in the form of an additional source in the ``add_library`` command.
Similarly, new Cython modules must be added to the appropriate list in the ``freud/CMakeLists.txt`` file depending on whether or not there is C++ code associated with the module.
Finally, you will need to import the new module in ``freud/__init__.py`` by adding :code:`from . import MODULENAME` so that your module is usable as ``freud.MODULENAME``.

Once the code is added, appropriate tests should be added to the ``tests`` folder.
Test files are named by the convention ``tests/test_MODULENAME_CLASSNAME.py``.
The final step is updating documentation, which is contained in ``rst`` files named with the convention ``doc/source/modules/MODULENAME.rst``.
If you have added a class to an existing module, all you have to do is add that same class to the ``autosummary`` section of the corresponding ``rst`` file.
If you have created a new module, you will have to create the corresponding ``rst`` file with the summary section listing classes and functions in the module followed by a more detailed description of all classes.
All classes and functions should be documented inline in the code, which allows automatic generation of the detailed section using the ``automodule`` directive (see any of the module ``rst`` files for an example).
Finally, the new file needs to be added to ``doc/source/index.rst`` in the ``API`` section.
=================
Memory Management
=================

Memory handling in **freud** is a somewhat intricate topic.
Most **freud** developers do not need to be aware of such details; however, certain practices must be followed to ensure that the expected behavior is achieved.
This page provides an overview of how data should be handled in **freud** and how module developers should use **freud**'s core classes to ensure proper memory management.
A thorough description of the process is also provided for developers who need to understand the internal logic for further development.

.. note::

    This page specifically deals with modules primarily written in C++. These
    concepts do not apply to pure Python/Cython modules.

Problem Statement
=================

The standard structure for **freud** modules involves a core implementation in a C++ class wrapped in a Cython class that owns a pointer to the C++ object.
Python ``compute`` methods call through to C++ ``compute`` methods, which perform the calculation and populate class member arrays that are then accessed via properties of the owning Cython class.
These classes are designed to be reusable, i.e. ``compute`` may be called many times on the same object with different data, and the accessed properties will return the most current data.
Users have a reasonable expectation that if the accessed property is saved to another variable it will remain unchanged by future calls to ``compute`` or if the originating C++ object is destructed, but a naive implementation that ensures this invariant would involve reallocating memory on every call to compute, an unnecessarily expensive operation.
Ultimately, what we want is a method that performs the minimal number of memory allocations while allowing users to operate transparently on outputs without worrying about whether the data will be invalidated by future operations.

ManagedArray
============

The **freud** ``ManagedArray`` template class provides a solution to this problem for arbitrary types of numerical data.
Proper usage of the class can be summarized by the following steps:

#. Declaring ``ManagedArray`` class members in C++.
#. Calling the ``prepare`` method in every ``compute``.
#. Making the array accessible via a getter method that **returns a const reference**.
#. Calling ``make_managed_numpy_array`` in Cython and returning the output as a property.

Plenty of examples of following this pattern can be found throughout the codebase, but for clarity we provide a complete description with examples below.
If you are interested in more details on the internals of ``ManagedArray`` and how it actually works, you can skip to :ref:`managedarray_explained`.


Using ManagedArrays
-------------------

We'll use :mod:`freud.cluster.Cluster` to illustrate how the four steps above may be implemented.
This class takes in a set of points and assigns each of them to clusters, which are store in the C++ array ``m_cluster_idx``.

Step 1 is simple: we note that ``m_cluster_idx`` is a member variable of type ``ManagedArray<unsigned int>``.
For step 2, we look at the first few lines of ``Cluster::compute``, where we see a call to ``m_cluster_idx.prepare``.
This method encapsulates the core logic of ``ManagedArray``, namely the intelligent reallocation of memory whenever other code is still accessing it.
This means that, if a user saves the corresponding Python property :attr:`freud.cluster.Cluster.cluster_idx` to a local variable in a script and then calls :attr:`freud.cluster.Cluster.compute`, the saved variable will still reference the original data, and the new data may be accessed again using :attr:`freud.cluster.Cluster.cluster_idx`.

Step 3 for the cluster indices is accomplished in the following code block:

.. code-block:: C++

    //! Get a reference to the cluster ids.
    const util::ManagedArray<unsigned int> &getClusterIdx() const
    {
        return m_cluster_idx;
    }

**The return type of this method is crucial: all such methods must return const references to the members**.

The final step is accomplished on the Cython side.
Here is how the cluster indices are exposed in :class:`freud.cluster.Cluster`:

.. code-block:: python

    @_Compute._computed_property
    def cluster_idx(self):
        """:math:`N_{points}` :class:`numpy.ndarray`: The cluster index for
        each point."""
        return freud.util.make_managed_numpy_array(
            &self.thisptr.getClusterIdx(),
            freud.util.arr_type_t.UNSIGNED_INT)

Essentially all the core logic is abstracted away from the user through the :func:`freud.data.make_managed_numpy_array`, which creates a NumPy array that is a view on an existing ``ManagedArray``.
This NumPy array will, in effect, take ownership of the data in the event that the user keeps a reference to it and requests a recomputation.
Note the signature of this function: the first argument must be **a pointer to the ManagedArray** (which is why we had to return it by reference), and the second argument indicates the type of the data (the possible types can be found in ``freud/util.pxd``).
There is one other point to note that is not covered by the above example; if the template type of the ``ManagedArray`` is not a scalar, you also need to provide a third argument indicating the size of this vector.
The most common use-case is for methods that return an object of type ``ManagedArray<vec3<float>>``: in this case, we would call ``make_managed_numpy_array(&GETTER_FUNC, freud.util.arr_type_t.FLOAT, 3)``.


Indexing ManagedArrays
----------------------

With respect to indexing, the ``ManagedArray`` class behaves like any standard array-like container and can be accessed using e.g. ``m_cluster_idx[index]``.
In addition, because many calculations in **freud** output multidimensional information, ``ManagedArray`` also supports multidimensional indexing using ``operator()``.
For example, setting the element at second row and third column of a 2D ``ManagedArray`` ``array`` to one can be done using ``array(1, 2) = 1`` (indices beginning from 0).
Therefore, ``ManagedArray`` objects can be used easily inside the core C++ calculations in **freud**.


.. _managedarray_explained:

Explaining ManagedArrays
------------------------

We now provide a more detailed accounting of how the ``ManagedArray`` class actually works.
Consider the following block of code:

.. code-block:: python

    rdf = freud.density.RDF(bins=100, r_max=3)

    rdf.compute(system=(box1, points1))
    rdf1 = rdf.rdf

    rdf2.compute(system=(box2, points2))
    rdf2 = rdf.rdf

We require that ``rdf1`` and ``rdf2`` be distinct arrays that are only equal if the results of computing the RDF are actually equivalent for the two systems, and we want to achieve this with the minimal number of memory allocations.
In this case, that means there are two required allocations; returning copies would double that.

To achieve this goal, ``ManagedArray`` objects store a pointer to a pointer.
Multiple ``ManagedArray`` objects can point to the same data array, and the pointers are all shared pointers to automate deletion of arrays when no references remain.
The key using the class properly is the ``prepare`` method, which checks the reference count to determine whether it's safe to simply zero out the existing memory or if it needs to allocate a new array.
In the above example, when ``compute`` is called a second time the ``rdf1`` object in Python still refers to the computed data, so ``prepare`` will detect that there are multiple (two) shared pointers pointing to the data and choose to reallocate the class's ``ManagedArray`` storing the RDF.
By calling ``prepare`` at the top of every ``compute`` method, developers ensure that the array used for the rest of the method has been properly zeroed out, and they do not need to worry about whether reallocation is needed (including cases where array sizes change).

To ensure that all references to data are properly handled, some additional logic is required on the Python side as well.
The Cython ``make_managed_numpy_array`` instantiates a ``_ManagedArrayContainer`` class, which is essentially just a container for a ``ManagedArray`` that points to the same data as the ``ManagedArray`` provided as an argument to the function.
This link is what increments the underlying shared pointer reference counter.
The ``make_managed_numpy_array`` uses the fact that a ``_ManagedArrayContainer`` can be transparently converted to a NumPy array that points to the container; as a result, no data copies are made, but all NumPy arrays effectively share ownership of the data along with the originating C++ class.
If any such arrays remain in scope for future calls to ``compute``, ``prepare`` will recognize this and reallocate memory as needed.
================
Neighbor Finding
================

Neighbor finding is central to many methods in **freud**.
The purpose of this page is to introduce the various C++ classes and functions involved in the neighbor finding process.
This page focuses on aspects of the neighbor finding utilities in **freud** that are important for developers.
As such, it assumes knowledge on the level of :ref:`neighbors` and :ref:`paircompute`, so please familiarize yourself with the contents of those pages before proceeding.

There are two primary use-cases for the neighbor finding code in **freud**.
One is to directly expose this functionality to the user, via the :class:`~freud.locality.NeighborQuery` abstract class and its subclasses.
The second it to enable looping over nearest neighbors (as defined by arbitrary query arguments or a precomputed :class:`~freud.locality.NeighborList`) inside of compute methods defined in C++.
To support both of these use-cases, **freud** defines how to find neighbors inside iterator classes, which can be naturally looped over in either case.
In this page, we first discuss these iterators and how they are structured with respect to the ``locality::NeighborQuery`` C++ class.
We then discuss the utility functions built around this class to enable easier C++ computations, at which point we also discuss how :class:`~freud.locality.NeighborList` objects fit into this framework.

Per-Point Iterators and the NeighborQuery class
===============================================

The lowest-level unit of the neighbor finding infrastructure in **freud** is the ``locality::NeighborPerPointIterator``.
This class defines an abstract interface for all neighbor iteration in **freud**, an interface essentially composed of the ``next`` and ``end`` methods.
Given an instance of this class, these two methods provide the means for client code to get the next neighbor in the iteration until there are no further neighbors.
Calls to ``next`` produces instances of ``locality::NeighborBond``, a simple data class that contains the core information defining a bond (a pair of points, a distance, and any useful ancillary information).

The rationale for making per-point iteration the central element is twofold.
The first is conceptual: all logic for finding neighbors is naturally reducible to a set of conditions on the neighbors of each query point.
The second is more practical: since finding neighbors for each point must be sequential in many cases (such as nearest neighbor queries), per-point iteration is the smallest logical unit that can be parallelized.

Instances of ``locality::NeighborPerPointIterator`` should not be constructed directly; instead, they should be constructed via the ``locality::NeighborQuery::querySingle`` method.
The ``locality::NeighborQuery`` class is an abstract data type that defines an interface for any method implemented for finding neighbors.
All subclasses must implement the ``querySingle`` method, which should return a subclass of ``locality::NeighborPerPointIterator``.
For instance, the ``locality::AABBQuery`` class implements ``querySingle``, which returns an instance of the ``locality::AABBIterator`` subclass.
In general, different ``NeighborQuery`` subclasses will need to implement separate per-point iterators for each query mode; for ``locality::AABBQuery``, these are the ``locality::AABBQueryIterator`` and the ``locality::AABBQueryBallIterator``, which encode the logic for nearest neighbor and ball queries, respectively.

Although the ``querySingle`` method is what subclasses should implement, the primary interface to ``NeighborQuery`` subclasses is the ``query`` method, which accepts an arbitrary set of query points and query arguments and simply generates a ``locality::NeighborQueryIterator`` object.
The ``locality::NeighborQueryIterator`` class is an intermediary that allows lazy generation of neighbors.
It essentially functions as the container for a set of points, query points, and query arguments; once iteration of this object begins, it produces ``NeighborPerPointIterator`` objects on demand.
This mode of operation also enables the generator approach to looping over neighbors in Python, since iterating in Python corresponds directly to calling ``next`` on the underlying ``NeighborQueryIterator``.

There is one conceptual complexity associated with this class that is important to understand.
Since all of the logic for finding neighbors is contained in the per-point iterator classes, the ``NeighborQueryIterator`` actually retains a reference to the constructing ``NeighborQuery`` object so that it can call ``querySingle`` for each point.
This bidirectionally linked structure enables encapsulation of the neighbor finding logic while also supporting easily parallelization (via parallel calls to ``querySingle``).
Additionally, this structure makes it natural to generate ``NeighborList`` objects.

NeighborLists
=============

The ``NeighborList`` class represents a static list of neighbor pairs.
The results of any query can be converted to a ``NeighborList`` by calling the ``toNeighborList`` method of the ``NeighborQueryIterator``, another reason why the iterator class logic is separated from the ``NeighborQuery`` object: it allows generation of a ``NeighborList`` from -- and more generally, independent operation on -- the result of a query.
The ``NeighborList`` is simply implemented as a collection of raw arrays, one of which holds pairs of neighbors.
The others hold any additional information associated with each bond, such as distances or weights.

By definition, the bonds in a ``NeighborList`` are stored in the form ``(query_point, point)`` (i.e. this is how the underlying array is indexed) and ordered by ``query_point``.
This ordering makes the structure amenable to a fast binary search algorithm.
Looping over neighbors of a given ``query_point`` is then simply a matter of finding the first index in the list where that ``query_point`` appears and then iterating until the ``query_point`` index in the list no longer matches the one under consideration.

Computing With Neighbors
========================

One of the most common operations in **freud** is performing some computation over all neighbor-bonds in the system.
Users have multiple ways of specifying neighbors (using query arguments or by a :class:`~freud.locality.NeighborList`), so **freud** provides some utility functions to abstract the process of looping over neighbors.
These functions are defined in ``locality/NeighborComputeFunctional.h``; the two most important ones are ``loopOverNeighbors`` and ``loopOverNeighborsIterator``.
Compute functions that perform neighbor computations typically accept a ``NeighborQuery``, a ``QueryArgs``, and a ``NeighborList`` object.
These objects can then be passed to either of the utility functions, which loop over the ``NeighborList`` if it was provided (if no :class:`~freud.locality.NeighborList` is provided by the Python user, a ``NULL`` pointer is passed through), and if not, perform a query on the ``NeighborQuery`` object using the provided ``QueryArgs`` to generate the required neighbors.
The actual computation should be encapsulated as a lambda function that is passed as an argument to these utilities.

The distinction between the two utility functions lies in the signature of the accepted lambda functions, which enables a slightly different form of computation.
The default ``loopOverNeighbors`` function does exactly what is described above, namely it calls the provided compute function for every single bond.
However, some computations require some additional code to be executed for each ``query_point``, such as some sort of normalization.
To enable this mode of operation, the ``loopOverNeighborsIterator`` method instead requires a lambda function that accepts two arguments, the ``query_point`` index and a ``NeighborPerPointIterator``.
This way, the client code can loop over the neighbors of a given ``query_point`` and perform the needed computation, then execute additional code (which may optionally depend on the index of the ``query_point``, e.g. to update a specific array index).

Default Systems
---------------

There is one important implementation detail to note.
The user is permitted to simply provide a set of points rather than a ``NeighborQuery`` object on the Python side (i.e. any valid argument to :meth:`~freud.locality.NeighborQuery.from_system`), but we need a natural way to mirror this in C++, ideally without too many method overloads.
To implement this, we provide the ``RawPoints`` C++ class and its Python :class:`~freud.locality._RawPoints` mirror, which is essentially a plain container for a box and a set of query points.
This object inherits from ``NeighborQuery``, allowing it to be passed directly into the C++ compute methods.

However, neighbor computations still need to know how to find neighbors.
In this case, they must construct a ``NeighborQuery`` object capable of neighbor finding and then use the provided query arguments to find neighbors.
To enable this calculation, the ``RawPoints`` class implements a query method that simply constructs an ``AABBQuery`` internally and queries it for neighbors.

Default NeighborLists
---------------------

Some compute methods are actually computations that produce quantities per bond.
One example is the ``SolidLiquid`` order parameter, which computes an order parameter value for each bond.
The ``NeighborComputeFunctional.h`` file implements a ``makeDefaultNList`` function that supports this calculation by creating a ``NeighborList`` object from whatever inputs are provided on demand.
==========
Histograms
==========

Histograms are a common type of calculation implemented in **freud** because custom histograms are hard to compute efficiently in pure Python.
The C++ ``Histogram`` class support weighted N-dimensional histograms with different spacings in each dimension.
The key to this flexibility is the ``Axis`` class, which defines the range spacing along a single axis; an N-dimensional ``Histogram`` is composed of a sequence of N ``Axis`` objects.
Binning values into the histogram is performed by binning along each axis.
The standard ``RegularAxis`` subclass of ``Axis`` defines an evenly spaced axis with bin centers defined as the center of each bin; additional subclasses may be defined to add different spacing if desired.

Multithreading is achieved through the ``ThreadLocalHistogram`` class, which is a simple wrapper around the ``Histogram`` that creates an equivalent histogram on each thread.
The standard pattern for parallel histogramming is to generate a ``ThreadLocalHistogram`` and add data into it, then call the ``Histogram::reduceOverThreads`` method to accumulate these data into a single histogram.
In case any additional post-processing is required per bin, it can also be executed in parallel by providing it as a lambda function to ``Histogram::reduceOverThreadsPerBin``.


Computing with Histograms
=========================

The ``Histogram`` class is designed as a data structure for the histogram.
Most histogram computations in **freud** involve standard neighbor finding to get bonds, followed by binning some function of these bonds into a histogram.
Examples include RDFs (binning bond distances), PMFTs (binning bonds by the different vector components of the bond), and bond order diagrams (binning bond angles).
An important distinction between histogram computations and most others is that histograms naturally support an accumulation of information over multiple frames of data, an operation that is ill-defined for many other computations.
As a result, histogram computations also need to implement some boilerplate for handling accumulating and averaging data over multiple frames.

The details of these computations are encapsulated by the ``BondComputeHistogram`` class, which contains a histogram, provides accessors to standard histogram properties like bin counts and axis sizes, and has a generic accumulation method that accepts a lambda compute function.
This signature is very similar to the utility functions for looping over neighbors, and in fact the function is transparently forwarded to ``locality::loopOverNeighbors``.
Any compute that matches this pattern should inherit from the ``BondComputeHistogram`` class and must implement an ``accumulate`` method to perform the computation and a ``reduce`` to reduce thread local histograms into a single histogram..
.. _installation:

============
Installation
============

Installing freud
================

The **freud** library can be installed via `conda <https://conda.io/projects/conda/>`_ or pip, or compiled from source.

Install via conda
-----------------

The code below will install **freud** from `conda-forge <https://anaconda.org/conda-forge/freud>`_.

.. code-block:: bash

    conda install -c conda-forge freud

Install via pip
-----------------

The code below will install **freud** from `PyPI <https://pypi.org/project/freud-analysis/>`_.

.. code-block:: bash

    pip install freud-analysis

Compile from source
-------------------

The following are **required** for installing **freud**:

- A C++14-compliant compiler
- `Python <https://www.python.org/>`__ (>=3.6)
- `NumPy <https://www.numpy.org/>`__ (>=1.14)
- `Intel Threading Building Blocks <https://www.threadingbuildingblocks.org/>`__ (>=2017.03.R2)
- `Cython <https://cython.org/>`__ (>=0.29.14)
- `scikit-build <https://scikit-build.readthedocs.io/>`__ (>=0.10.0)
- `CMake <https://cmake.org/>`__ (>=3.12.0)

.. note::

    Depending on the generator you are using, you may require a newer version of CMake.
    In particular, Visual Studio 2019 requires CMake >= 3.14.0.
    For more information on specific generators, see the `CMake generator documentation <https://cmake.org/cmake/help/git-stage/manual/cmake-generators.7.html>`__.

The **freud** library uses scikit-build and CMake to handle the build process itself, while the other requirements are required for compiling code in **freud**.
These requirements can be met by installing the following packages from the `conda-forge channel <https://conda-forge.org/>`__:

.. code-block:: bash

    conda install -c conda-forge tbb tbb-devel numpy cython scikit-build cmake

All requirements other than TBB can also be installed via the `Python Package Index <https://pypi.org/>`__

.. code-block:: bash

    pip install numpy cython scikit-build cmake

Wheels for tbb and tbb-devel exist on PyPI, but only for certain operating systems, so your mileage may vary.
For non-conda users, we recommend using OS-specific package managers (e.g. `Homebrew <https://brew.sh/>`__ for Mac) to install TBB.
As in the snippets above, it may be necessary to install both both a TBB and a "devel" package in order to get both the headers and the shared libraries.

The code that follows builds **freud** and installs it for all users (append ``--user`` if you wish to install it to your user site directory):

.. code-block:: bash

    git clone --recurse-submodules https://github.com/glotzerlab/freud.git
    cd freud
    python setup.py install

You can also build **freud** in place so that you can run from within the folder:

.. code-block:: bash

    # Run tests from the tests directory
    python setup.py build_ext --inplace

Building **freud** in place has certain advantages, since it does not affect your Python behavior except within the **freud** directory itself (where **freud** can be imported after building).

CMake Options
+++++++++++++

The scikit-build tool allows setup.py to accept three different sets of options separated by ``--``, where each set is provided directly to scikit-build, to CMake, or to the code generator of choice, respectively.
For example, the command ``python setup.py build_ext --inplace -- -DCOVERAGE=ON -G Ninja -- -j 4`` tell scikit-build to perform an in-place build, it tells CMake to turn on the ``COVERAGE`` option and use Ninja for compilation, and it tells Ninja to compile with 4 parallel threads.
For more information on these options, see the `scikit-build docs <scikit-build.readthedocs.io/>`__.

.. note::

    The default CMake build configuration for freud is ``ReleaseWithDocs`` (not a standard build configuration like ``Release`` or ``RelWithDebInfo``).
    On installation, ``setup.py`` assumes ``--build-type=ReleaseWithDocs`` by default if no build type is specified.
    Using this build configuration is a workaround for `this issue <https://github.com/scikit-build/scikit-build/issues/518>`__ with scikit-build and Cython embedding docstrings.

In addition to standard CMake flags, the following CMake options are available for **freud**:

.. glossary::

    \--COVERAGE
      Build the Cython files with coverage support to check unit test coverage.


The **freud** CMake configuration also respects the following environment variables (in addition to standards like ``LD_LIBRARY_PATH``).

.. glossary::

    TBBROOT
      The root directory where TBB is installed.
      Useful if TBB is installed in a non-standard location or cannot be located for some other reason.
      This variable is set by the ``tbbvars.sh`` script included with TBB when building from source.

    TBB_INCLUDE_DIR
      The directory where the TBB headers (e.g. ``tbb.h``) are located.
      Useful if TBB is installed in a non-standard location or cannot be located for some other reason.

.. note::

    **freud** makes use of git submodules. To manually update git submodules, execute:

    .. code-block:: bash

        git submodule update --init --recursive

Unit Tests
==========

The unit tests for **freud** are included in the repository and are configured to be run using the Python :mod:`pytest` library:

.. code-block:: bash

    # Run tests from the tests directory
    cd tests
    python -m pytest .

Note that because **freud** is designed to require installation to run (i.e. it cannot be run directly out of the build directory), importing **freud** from the root of the repository will fail because it will try and import the package folder.
As a result, unit tests must be run from outside the root directory if you wish to test the installed version of **freud**.
If you want to run tests within the root directory, you can instead build **freud** in place:

.. code-block:: bash

    # Run tests from the tests directory
    python setup.py build_ext --inplace

This build will place the necessary files alongside the **freud** source files so that **freud** can be imported from the root of the repository.

Documentation
=============

The documentation for **freud** is `hosted online at ReadTheDocs <https://freud.readthedocs.io/>`_.
You may also build the documentation yourself.

Building the documentation
--------------------------

The following are **required** for building **freud** documentation:

- `Sphinx <http://www.sphinx-doc.org/>`_
- `Read the Docs Sphinx Theme <https://sphinx-rtd-theme.readthedocs.io/>`_
- `nbsphinx <https://nbsphinx.readthedocs.io/>`_
- `jupyter_sphinx <https://jupyter-sphinx.readthedocs.io/>`_
- `sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/>`_

You can install these dependencies using conda:

.. code-block:: bash

    conda install -c conda-forge sphinx sphinx_rtd_theme nbsphinx jupyter_sphinx sphinxcontrib-bibtex

or pip:

.. code-block:: bash

    pip install sphinx sphinx-rtd-theme nbsphinx jupyter-sphinx sphinxcontrib-bibtex

To build the documentation, run the following commands in the source directory:

.. code-block:: bash

    cd doc
    make html
    # Then open build/html/index.html

To build a PDF of the documentation (requires LaTeX and/or PDFLaTeX):

.. code-block:: bash

    cd doc
    make latexpdf
    # Then open build/latex/freud.pdf
.. _tutorial:

========
Tutorial
========

This tutorial provides a complete introduction to **freud**.
Rather than attempting to touch on all features in **freud**, it focuses on common core concepts that will help understand how **freud** works with data and exposes computations to the user.
The tutorial begins by introducing the fundamental concepts of periodic systems as implemented in **freud** and the concept of ``Compute classes``, which consitute the primary API for performing calculations with **freud**.
The tutorial then discusses the most common calculation performed in **freud**, finding neighboring points in periodic systems.
The package's neighbor finding tools are tuned for high performance neighbor finding, which is what enables most of other calculations in **freud**, which typically involve characterizing local environments of points in some way.
The next part of the tutorial discusses the role of histograms in **freud**, focusing on the common features and properties that all histograms share.
Finally, the tutorial includes a few more complete demonstrations of using **freud** that should provide reasonable templates for use with almost any other features in **freud**.

.. toctree::
   :maxdepth: 2

   tutorial/periodic.rst
   tutorial/computeclass.rst
   tutorial/neighborfinding
   tutorial/paircompute.rst
.. _quickstart:

================
Quickstart Guide
================

Once you have `installed freud <installation.rst>`_, you can start using **freud** with any simulation data that you have on hand.
As an example, we'll assume that you have run a simulation using the `HOOMD-blue <https://glotzerlab.engin.umich.edu/hoomd-blue/>`_ and used the :class:`hoomd.dump.gsd` command to output the trajectory into a file ``trajectory.gsd``.
The `GSD file format <https://gsd.readthedocs.io/en/stable/>`_ provides its own convenient Python file reader that offers access to data in the form of NumPy arrays, making it immediately suitable for calculation with **freud**. Many other file readers and data formats are supported, see :ref:`datainputs` for a full list and more examples.

We start by reading the data into a NumPy array:

.. code-block:: python

    import gsd.hoomd
    traj = gsd.hoomd.open('trajectory.gsd', 'rb')


We can now immediately calculate important quantities.
Here, we will compute the radial distribution function :math:`g(r)` using the :class:`freud.density.RDF` compute class.
Since the radial distribution function is in practice computed as a histogram, we must specify the histogram bin widths and the largest interparticle distance to include in our calculation.
To do so, we simply instantiate the class with the appropriate parameters and then perform a computation on the given data:

.. code-block:: python

    import freud
    rdf = freud.density.RDF(bins=50, r_max=5)
    rdf.compute(system=traj[-1])

We can now access the data through properties of the ``rdf`` object.

.. code-block:: python

    r = rdf.bin_centers
    y = rdf.rdf

Many classes in **freud** natively support plotting their data using `Matplotlib <https://matplotlib.org/>`:

.. code-block:: python

    import matplotlib as plt
    fig, ax = plt.subplots()
    rdf.plot(ax=ax)

You will note that in the above example, we computed :math:`g(r)` only using the final frame of the simulation trajectory, ``traj[-1]``.
However, in many cases, radial distributions and other similar quantities may be noisy in simulations due to the natural fluctuations present.
In general, what we are interested in are *time-averaged* quantities once a system has equilibrated.
To perform such a calculation, we can easily modify our original calculation to take advantage of **freud**'s *accumulation* features.
To accumulate, just add the argument ``reset=False`` with a supported compute object (such as histogram-like computations).
Assuming that you have some method for identifying the frames you wish to include in your sample, our original code snippet would be modified as follows:

.. code-block:: python

    import freud
    rdf = freud.density.RDF(bins=50, r_max=5)
    for frame in traj:
        rdf.compute(frame, reset=False)

You can then access the data exactly as we previously did.
And that's it!

Now that you've seen a brief example of reading data and computing a radial distribution function, you're ready to learn more.
If you'd like a complete walkthrough please see the :ref:`tutorial`.
The tutorial walks through many of the core concepts in **freud** in greater detail, starting with the basics of the simulation systems we analyze and describing the details of the neighbor finding logic in **freud**.
To see specific features of **freud** in action, look through the :ref:`examples`.
More detailed documentation on specific classes and functions can be found in the `API documentation <modules>`_.
============
Introduction
============

The **freud** library is a Python package for analyzing particle simulations.
The package is designed to directly use numerical arrays of data, making it easy to use for a wide range of use-cases.
The most common use-case of **freud** is for computing quantities from molecular dynamics simulation trajectories, but it can be used for analyzing any type of particle simulation.
By operating directly on numerical arrays of data, **freud** allows users to parse custom simulation outputs into a suitable structure for input, rather than relying specific file types or data structures.

The core of **freud** is analysis of periodic systems, which are represented through the :class:`freud.box.Box` class.
The :class:`freud.box.Box` supports arbitrary triclinic systems for maximum flexibility, and is used throughout the package to ensure consistent treatment of these systems.
The package's many methods are encapsulated in various *compute classes*, which perform computations and populate class attributes for access.
Of particular note are the various computations based on nearest neighbor finding in order to characterize particle environments.
Such methods are simplified and accelerated through a centralized neighbor finding interface defined in the :class:`freud.locality.NeighborQuery` family of classes in the :mod:`freud.locality` module of freud.
.. _examples:

========
Examples
========

Examples are provided as `Jupyter <https://jupyter.org/>`_ notebooks in a separate
`freud-examples <https://github.com/glotzerlab/freud-examples>`_ repository.
These notebooks may be launched `interactively on Binder <https://mybinder.org/v2/gh/glotzerlab/freud-examples/master?filepath=index.ipynb>`_
or downloaded and run on your own system.
Visualization of data is done via `Matplotlib <https://matplotlib.org/>`_ and `Bokeh <https://bokeh.pydata.org/>`_, unless otherwise noted.


Key concepts
============

There are a few critical concepts, algorithms, and data structures that are central to all of **freud**.
The :class:`freud.box.Box` class defines the concept of a periodic simulation box, and the :mod:`freud.locality` module defines methods for finding nearest neighbors of particles.
Since both of these are used throughout **freud**, we recommend reading the :ref:`tutorial` first, before delving into the workings of specific **freud** analysis modules.

.. toctree::
    :maxdepth: 1
    :glob:

    examples/module_intros/box*
    examples/module_intros/locality*

Analysis Modules
================

These introductory examples showcase the functionality of specific modules in **freud**, showing how they can be used to perform specific types of analyses of simulations.

.. toctree::
    :maxdepth: 1
    :glob:

    examples/module_intros/cluster*
    examples/module_intros/density*
    examples/module_intros/diffraction*
    examples/module_intros/environment*
    examples/module_intros/interface*
    examples/module_intros/order*
    examples/module_intros/pmft*

Example Analyses
================

The examples below go into greater detail about specific applications of **freud** and use cases that its analysis methods enable, such as user-defined analyses, machine learning, and data visualization.

.. toctree::
    :maxdepth: 1
    :glob:

    examples/examples/NetworkX-CNA
    examples/examples/HOOMD-MC-W6/HOOMD-MC-W6
    examples/examples/GROMACS-MDTRAJ-WATER-RDF/Compute_RDF
    examples/examples/LAMMPS-LJ-MSD/LAMMPS-LJ-MSD
    examples/examples/Using Machine Learning for Structural Identification
    examples/examples/Handling Multiple Particle Types (A-B Bonds)
    examples/examples/Calculating RDF from GSD files
    examples/examples/Calculating Strain via Voxelization
    examples/examples/Visualization with fresnel
    examples/examples/Visualization with plato
    examples/examples/Visualizing 3D Voronoi and Voxelization


Benchmarks
==========

Performance is a central consideration for **freud**. Below are some benchmarks comparing **freud** to other tools offering similar analysis methods.

.. toctree::
    :maxdepth: 1
    :glob:

    examples/examples/Benchmarking*
.. _paircompute:

=================
Pair Computations
=================

Some computations in **freud** do not depend on bonds at all.
For example, :class:`freud.density.GaussianDensity` creates a "continuous equivalent" of a system of points by placing normal distributions at each point's location to smear out its position, then summing the value of these distributions at a set of fixed grid points.
This calculation can be quite useful because it allows the application of various analysis tools like fast Fourier transforms, which require regular grids.
For the purposes of this tutorial, however, the importance of this class is that it is an example of a calculation where neighbors are unimportant: the calculation is performed on a per-point basis only.

The much more common pattern in **freud**, though, is that calculations involve the local neighborhoods of points.
To support efficient, flexible computations of such quantities, various ``Compute classes`` essentially expose the same API as the query interface demonstrated in the previous section.
These ``PairCompute classes`` are designed to mirror the querying functionality of **freud** as closely as possible.

As an example, let's consider :class:`freud.density.LocalDensity`, which calculates the density of points in the local neighborhood of each point.
Adapting our code from the `previous section <neighborfinding>`_, the simplest usage of this class would be as follows:

.. code-block:: python

    import numpy as np
    import freud

    L = 10
    num_points = 100

    points = np.random.rand(num_points)*L - L/2
    box = freud.box.Box.cube(L)

    # r_max specifies how far to search around each point for neighbors
    r_max = 2

    # For systems where the points represent, for instance, particles with a
    # specific size, the diameter is used to add fractional volumes for
    # neighbors that would be overlapping the sphere of radius r_max around
    # each point.
    diameter = 0.001

    ld = freud.density.LocalDensity(r_max, diameter)
    ld.compute(system=(box, points))

    # Access the density.
    ld.density

Using the same example system we've been working with so far, we've now calculated an estimate for the number of points in the neighborhood of each point.
Since we already told the computation how far to search for neighbors based on ``r_max``, all we had to do was pass a :class:`tuple` ``(box, points)`` to compute indicate where the points were located.

Binary Systems
==============

Imagine that instead of a single set of points, we actually had two different types of points and we were interested in finding the density of one set of points in the vicinity of the other.
In that case, we could modify the above calculation as follows:


.. code-block:: python

    import numpy as np
    import freud
    L = 10
    num_points = 100
    points = np.random.rand(num_points)*L - L/2
    query_points = np.random.rand(num_points/10)*L - L/2

    r_max = 2
    diameter = 0.001

    ld = freud.density.LocalDensity(r_max, diameter)
    ld.compute(system=(box, points), query_points=query_points)

    # Access the density.
    ld.density

The choice of names names here is suggestive of exactly what this calculation is now doing.
Internally, :class:`freud.density.LocalDensity` will search for all ``points`` that are within the cutoff distance ``r_max`` of every ``query_point`` (essentially using the query interface we introduced previously) and use that to calculate ``ld.density``.
Note that this means that ``ld.density`` now contains densities for every ``query_point``, i.e. it is of length 10, not 100.
Moreover, recall that one of the features of the querying API is the specification of whether or not to count particles as their own neighbors.
``PairCompute classes`` will attempt to make an intelligent determination of this for you; if you do not pass in a second set of ``query_points``, they will assume that you are computing with a single set of points and automatically exclude self-neighbors, but otherwise all neighbors will be included.

So far, we have included all points within a fixed radius; however, one might instead wish to consider the density in some shell, such as the density between 1 and 2 distance units away.
To address this need, you could simply adapt the call to ``compute`` above as follows:

.. code-block:: python

    ld.compute(system=(box, points), query_points=query_points,
               neighbors=dict(r_max=2, r_min=1))

The ``neighbors`` argument to ``PairCompute`` classes allows users to specify arbitary query arguments, making it possible to easily modify **freud** calculations on-the-fly.
The ``neighbors`` argument is actually more general than query arguments you've seen so far: if query arguments are not precise enough to specify the exact set of neighbors you want to compute with, you can instead provide a :class:`NeighborList <freud.locality.NeighborList>` directly

.. code-block:: python

    ld.compute(system=(box, points), query_points=query_points,
               neighbors=nlist)

This feature allows users essentially arbitrary flexibility to specify the bonds that should be included in any bond-based computation.
A common use-case for this is constructing a :class:`NeighborList <freud.locality.NeighborList>` using :class:`freud.locality.Voronoi`; Voronoi constructions provide a powerful alternative method of defining neighbor relationships that can improve the accuracy and robustness of certain calculations in **freud**.

You may have noticed in the last example that all the arguments are specified using keyword arguments.
As the previous examples have attempted to show, the ``query_points`` argument defines a second set of points to be used when performing calculations on binary systems, while the ``neighbors`` argument is how users can specify which neighbors to consider in the calculation.

The ``system`` argument is what, to this point, we have been specifying as a :class:`tuple` ``(box, points)``.
However, we don't have to use this tuple.
Instead, we can pass in any :class:`freud.locality.NeighborQuery`, the central class in **freud**'s querying infrastructure.
In fact, you've already seen examples of :class:`freud.locality.NeighborQuery`: the :class:`freud.locality.AABBQuery` object that we originally used to find neighbors.
There are also a number of other input types that can be converted via :meth:`freud.locality.NeighborQuery.from_system`, see also :ref:`datainputs`.
Since these objects all contain a :class:`freud.box.Box` and a set of points, they can be directly passed to computations:

.. code-block:: python

    aq = freud.locality.AABBQuery(box, points)
    ld.compute(system=aq, query_points=query_points, neighbors=nlist)

For more information on why you might want to use :class:`freud.locality.NeighborQuery` objects instead of the tuples, see :ref:`optimizing`.
For now, just consider this to be a way in which you can simplify your calls to many **freud** computes in one script by storing ``(box, points)`` into another objects.

You've now covered the most important information needed to use **freud**!
To recap, we've discussed how **freud** handles periodic boundary conditions, the structure and usage of ``Compute classes``, and methods for finding and performing calculations with pairs of neighbors.
For more detailed information on specific methods in **freud**, see the :ref:`examples` page or look at the API documentation for specific modules.
.. _pbcs:

============================
Periodic Boundary Conditions
============================

The central goal of **freud** is the analysis of simulations performed in periodic boxes.
Periodic boundary conditions are ubiquitous in simulations because they permit the simulation of quasi-infinite systems with minimal computational effort.
As long as simulation systems are sufficiently large, i.e. assuming that points in the system experience correlations over length scales substantially smaller than the system length scale, periodic boundary conditions ensure that the system appears effectively infinite to all points.

In order to consistently define the geometry of a simulation system with periodic boundaries, **freud** defines the :class:`freud.box.Box` class.
The class encapsulates the concept of a triclinic simulation box in a right-handed coordinate system.
Triclinic boxes are defined as parallelepipeds: three-dimensional polyhedra where every face is a parallelogram.
In general, any such box can be represented by three distinct, linearly independent box vectors.
Enforcing the requirement of right-handedness guarantees that the box can be represented by a matrix of the form

.. math::
    :nowrap:

    \begin{eqnarray*}
    \mathbf{h}& =& \left(\begin{array}{ccc} L_x & xy L_y & xz L_z \\
                                            0   & L_y    & yz L_z \\
                                            0   & 0      & L_z    \\
                         \end{array}\right)
    \end{eqnarray*}

where each column is one of the box vectors.

.. note::
    All **freud** boxes are centered at the origin, so for a given box the
    range of possible positions is :math:`[-L/2, L/2)`.

As such, the box is characterized by six parameters: the box vector lengths :math:`L_x`, :math:`L_y`, and :math:`L_z`, and the tilt factors :math:`xy`, :math:`xz`, and :math:`yz`.
The tilt factors are directly related to the angles between the box vectors.
All computations in **freud** are built around this class, ensuring that they naturally handle data from simulations conducted in non-cubic systems.
There is also native support for two-dimensional (2D) systems when setting :math:`L_z = 0`.

Boxes can be constructed in a variety of ways.
For simple use-cases, one of the factory functions of the :py:class:`freud.box.Box` provides the easiest possible interface:

.. code-block:: python

    # Make a 10x10 square box (for 2-dimensional systems).
    freud.box.Box.square(10)

    # Make a 10x10x10 cubic box.
    freud.box.Box.cube(10)

For more complex use-cases, the :py:meth:`freud.box.Box.from_box` method provides an interface to create boxes from any object that can easily be interpreted as a box.

.. code-block:: python

    # Create a 10x10 square box from a list of two items.
    freud.box.Box.from_box([10, 10])

    # Create a 10x10x10 cubic box from a list of three items.
    freud.box.Box.from_box([10, 10, 10])

    # Create a triclinic box from a list of six items (including tilt factors).
    freud.box.Box.from_box([10, 5, 2, 0.1, 0.5, 0.7])

    # Create a triclinic box from a dictionary.
    freud.box.Box.from_box(dict(Lx=8, Ly=7, Lz=10, xy=0.5, xz=0.7, yz=0.2))

    # Directly call the constructor.
    freud.box.Box(Lx=8, Ly=7, Lz=10, xy=0.5, xz=0.7, yz=0.2, dimensions=3)

More examples on how boxes can be created may be found in the API documentation of the :class:`Box <freud.box.Box>` class.
.. _computeclass:

===============
Compute Classes
===============

Calculations in **freud** are built around the concept of ``Compute classes``, Python objects that encode a given method and expose it through a ``compute`` method.
In general, these methods operate on a system composed of a triclinic box and a NumPy array of particle positions.
The box can be provided as any object that can be interpreted as a **freud** box (as demonstrated in the examples above).
We can look at the :class:`freud.order.Hexatic` order parameter calculator as an example:

.. code-block:: python

    import freud
    positions = ...  # Read positions from trajectory file.
    op = freud.order.Hexatic(k=6)
    op.compute(
        system=({'Lx': 5, 'Ly': 5, 'dimensions': 2}, positions),
        neighbors=dict(r_max=3)
    )

    # Plot the value of the order parameter.
    from matplotlib import pyplot as plt
    plt.hist(np.absolute(op.particle_order))

Here, we are calculating the hexatic order parameter, then using Matplotlib to plot.
The :class:`freud.order.Hexatic` class constructor accepts a single argument :code:`k`, which represents the periodicity of the calculation.
If you're unfamiliar with this order parameter, the most important piece of information here is that many compute methods in **freud** require parameters that are provided when the ``Compute class`` is constructed.

To calculate the order parameter we call :meth:`compute <freud.order.Hexatic.compute>`, which takes two arguments, a :class:`tuple` `(box, points)` and a :class:`dict`.
We first focus on the first argument.
The ``box`` is any object that can be coerced into a :class:`freud.box.Box` as described in the previous section; in this case, we use a dictionary to specify a square (2-dimensional) box.
The ``points`` must be anything that can be coerced into a 2-dimensional NumPy array of shape :code:`(N, 3)`
In general, the points may be provided as anything that can be interpreted as an :math:`N\times 3` list of positions; for more details on valid inputs here, see :func:`numpy.asarray`.
Note that because the hexatic order parameter is designed for two-dimensional systems, the points must be provided of the form `[x, y, 0]` (i.e. the z-component must be 0).
We'll go into more detail about the ``(box, points)`` tuple `soon <paircompute>`_, but for now, it's sufficient to just think of it as specifying the system of points we want to work with.

Now let's return to the ``neighbors`` argument to ``compute``, which is a dictionary is used to determine which particle neighbors to use.
Many computations in **freud** (such as the hexatic order parameter) involve the bonds in the system (for example, the average length of bonds or the average number of bonds a given point has).
However, the concept of a bond is sufficiently variable between different calculations; for instance, should points be considered bonded if they are within a certain distance of each other?
Should every point be considered bonded to a fixed number of other points?

To accommodate this variability, **freud** offers a very general framework by which bonds can be specified, and we'll go into more details in the `next section <neighborfinding>`_.
In the example above, we've simply informed the :class:`Hexatic <freud.order.Hexatic>` class that we want it to define bonds as pairs of particles that are less than 3 distance units apart.
We then access the computed order parameter as ``op.particle_order`` (we use :func:`np.absolute` because the output is a complex number and we just want to see its magnitude).


Accessing Computed Properties
=============================

In general, ``Compute classes`` expose their calculations using `properties <https://docs.python.org/3/library/functions.html#property>`_.
Any parameters to the ``Compute`` object (e.g. :code:`k` in the above example) can typically be accessed as soon as the object is constructed:

.. code-block:: python

    op = freud.order.Hexatic(k=6)
    op.k

Computed quantities can also be accessed in a similar manner, but only after the ``compute`` method is called.
For example:

.. code-block:: python

    op = freud.order.Hexatic(k=6)

    # This will raise an exception.
    op.particle_order

    op.compute(
        system=({'Lx': 5, 'Ly': 5, 'dimensions': 2}, positions),
        neighbors=dict(r_max=3)
    )

    # Now you can access this.
    op.particle_order

.. note::
    Most (but not all) of **freud**'s ``Compute classes`` are Python wrappers
    around high-performance implementations in C++. As a result, none of the
    data or the computations is actually stored in the Python object. Instead,
    the Python object just stores an instance of the C++ object that actually
    owns all its data, performs calculations, and returns computed quantities
    to the user. Python properties provide a nice way to hide this logic so
    that the Python code involves just a few lines.

Compute objects is that they can be used many times to calculate quantities, and the most recently calculated output can be accessed through the property.
If you need to perform a series of calculations and save all the data, you can also easily do that:

.. code-block:: python

    # Recall that lists of length 2 automatically convert to 2D freud boxes.
    box = [5, 5]

    op = freud.order.Hexatic(k=6)

    # Assuming that we have a list of Nx3 NumPy arrays that represents a
    # simulation trajectory, we can loop over it and calculate the order
    # parameter values in sequence.
    trajectory  = ...  # Read trajectory file into a list of positions by frame.
    hexatic_values = []
    for points in trajectory:
        op.compute(system=(box, points), neighbors=dict(r_max=3))
        hexatic_values.append(op.particle_order)


To make using **freud** as simple as possible, all ``Compute classes`` are designed to return ``self`` when compute is called.
This feature enables a very concise *method-chaining* idiom in **freud** where computed properties are accessed immediately:

.. code-block:: python

    particle_order = freud.order.Hexatic(k=6).compute(
        system=(box, points)).particle_order
.. _neighbors:

=================
Finding Neighbors
=================

Now that you've been introduced to the basics of interacting with **freud**, let's dive into the central feature of **freud**: efficiently and flexibly finding neighbors in periodic systems.

Problem Statement
=================

Neighbor-Based Calculations
---------------------------
As discussed in :ref:`the previous section <computeclass>`, a central task in many of the computations in **freud** is finding particles' neighbors.
These calculations typically only involve a limited subset of a particle's neighbors that are defined as characterizing its local environment.
This requirement is analogous to the force calculations typically performed in molecular dynamics simulations, where a cutoff radius is specified beyond which pair forces are assumed to be small enough to neglect.
Unlike in simulation, though, many analyses call for different specifications than simply selecting all points within a certain distance.

An important example is the calculation of order parameters, which can help characterize phase transitions.
Such parameters can be highly sensitive to the precise way in which neighbors are selected.
For instance, if a hard distance cutoff is imposed in finding neighbors for the hexatic order parameter, a particle may only be found to have five neighbors when it actually has six neighbors except the last particle is slightly outside the cutoff radius.
To accomodate such differences in a flexible manner, **freud** allows users to specify neighbors in a variety of ways.

Finding Periodic Neighbors
--------------------------

.. image:: ../images/PeriodicBoundaryConditions.png
    :width: 200px

Finding neighbors in periodic systems is significantly more challenging than in aperiodic systems.
To illustrate the difference, consider the figure above, where the black dashed line indicates the boundaries of the system.
If this system were aperiodic, the three nearest neighbors for point 1 would be points 5, 6, and 7.
However, due to periodicity, point 2 is actually closer to point 1 than any of the others if you consider moving straight through the top (or equivalently, the bottom) boundary.
Although many tools provide efficient implementations of algorithms for finding neighbors in aperiodic systems, they seldom generalize to periodic systems.
Even more rare is the ability to work not just in *cubic* periodic systems, which are relatively tractable, but in arbitrary triclinic geometries as described in :ref:`pbcs`.
This is precisely the type of calculation **freud** is designed for.


Neighbor Querying
=================

To understand how ``Compute classes`` find neighbors in **freud**, it helps to start by learning about **freud**'s neighbor finding classes directly.
Note that much more detail on this topic is available in the :ref:`querying` topic guide; in this section we will restrict ourselves to a higher-level overview.
For our demonstration, we will make use of the :class:`freud.locality.AABBQuery` class, which implements one fast method for periodic neighbor finding.
The primary mode of interfacing with this class (and other neighbor finding classes) is through the :meth:`query <freud.locality.AABBQuery>` interface.

.. code-block:: python

    import numpy as np
    import freud

    # As an example, we randomly generate 100 points in a 10x10x10 cubic box.
    L = 10
    num_points = 100

    # We shift all points into the expected range for freud.
    points = np.random.rand(num_points, 3)*L - L/2
    box = freud.box.Box.cube(L)
    aq = freud.locality.AABBQuery(box, points)

    # Now we generate a smaller sample of points for which we want to find
    # neighbors based on the original set.
    query_points = np.random.rand(num_points//10, 3)*L - L/2
    distances = []

    # Here, we ask for the 4 nearest neighbors of each point in query_points.
    for bond in aq.query(query_points, dict(num_neighbors=4)):
        # The returned bonds are tuples of the form
        # (query_point_index, point_index, distance). For instance, a bond
        # (1, 3, 0.2) would indicate that points[3] was one of the 4 nearest
        # neighbors for query_points[1], and that they are separated by a
        # distance of 0.2
        # (i.e. np.linalg.norm(query_points[1] - points[3]) == 2).
        distances.append(bond[2])

    avg_distance = np.mean(distances)

Let's dig into this script a little bit.
Our first step is creating a set of 100 points in a cubic box.
Note that the shifting done in the code above could also be accomplished using the :meth:`Box.wrap <freud.box.Box.wrap>` method like so: ``box.wrap(np.random.rand(num_points)*L)``.
The result would appear different, because if plotted without considering periodicity, the points would range from :code:`-L/2` to :code:`L/2` rather than from 0 to :code:`L`.
However, these two sets of points would be equivalent in a periodic system.

We then generate an additional set of ``query_points`` and ask for neighbors using the :meth:`query <freud.locality.AABBQuery>` method.
This function accepts two arguments: a set of points, and a :class:`dict` of **query arguments**.
Query arguments are a central concept in **freud** and represent a complete specification of the set of neighbors to be found.
In general, the most common forms of queries are those requesting either a fixed number of neighbors, as in the example above, or those requesting all neighbors within a specific distance.
For example, if we wanted to rerun the above example but instead find all bonds of length less than or equal to 2, we would simply replace the for loop above with:

.. code-block:: python

    for bond in aq.query(query_points, dict(r_max=2)):
        distances.append(bond[2])

Query arguments constitute a powerful method for specifying a query request.
Many query arguments may be combined for more specific purposes.
A common use-case is finding all neighbors within a single set of points (i.e. setting ``query_points = points`` in the above example).
In this situation, however, it is typically not useful for a point to find itself as a neighbor since it is trivially the closest point to itself and falls within any cutoff radius.
To avoid this, we can use the ``exclude_ii`` query argument:

.. code-block:: python

    query_points = points
    for bond in aq.query(query_points, dict(num_neighbors=4, exclude_ii=True)):
        pass

The above example will find the 4 nearest neighbors to each point, excepting the point itself.
A complete description of valid query arguments can be found in :ref:`querying`.

Neighbor Lists
==============

Query arguments provide a simple but powerful language with which to express neighbor finding logic.
Used in the manner shown above, :meth:`query <freud.locality.AABBQuery>` can be used to express many calculations in a very natural, Pythonic way.
By itself, though, the API shown above is somewhat restrictive because the output of :meth:`query <freud.locality.AABBQuery>` is a `generator <https://docs.python.org/3/glossary.html#term-generator>`_.
If you aren't familiar with generators, the important thing to know is that they can be looped over, *but only once*.
Unlike objects like lists, which you can loop over as many times as you like, once you've looped over a generator once, you can't start again from the beginning.

In the examples above, this wasn't a problem because we simply iterated over the bonds once for a single calculation.
However, in many practical cases we may need to reuse the set of neighbors multiple times.
A simple solution would be to simply to store the bonds into a list as we loop over them.
However, because this is such a common use-case, **freud** provides its own containers for bonds: the :class:`freud.locality.NeighborList`.

Queries can easily be used to generate :class:`NeighborList <freud.locality.NeighborList>` objects using their :meth:`toNeighborList <freud.locality.NeighborQueryResult.toNeighborList>` method:

.. code-block:: python

    query_result = aq.query(query_points, dict(num_neighbors=4, exclude_ii=True))
    nlist = query_result.toNeighborList()

The resulting object provides a persistent container for bond data.
Using :class:`NeighborLists <freud.locality.NeighborList>`, our original example might instead look like this:

.. code-block:: python

    import numpy as np
    import freud

    L = 10
    num_points = 100

    points = np.random.rand(num_points, 3)*L - L/2
    box = freud.box.Box.cube(L)
    aq = freud.locality.AABBQuery(box, points)

    query_points = np.random.rand(num_points//10, 3)*L - L/2
    distances = []

    # Here, we ask for the 4 nearest neighbors of each point in query_points.
    query_result = aq.query(query_points, dict(num_neighbors=4))
    nlist = query_result.toNeighborList()
    for (i, j) in nlist:
        # Note that we have to wrap the bond vector before taking the norm;
        # this is the simplest way to compute distances in a periodic system.
        distances.append(np.linalg.norm(box.wrap(query_points[i] - points[j])))

    avg_distance = np.mean(distances)

Note that in the above example we looped directly over the ``nlist`` and recomputed distances.
However, the ``query_result`` contained information about distances: here's how we access that through the ``nlist``:

.. code-block:: python

    assert np.allclose(nlist.distances, distances)

The indices are also accessible through properties, or through a NumPy-like slicing interface:


.. code-block:: python

    assert np.all(nlist.query_point_indices == nlist[:, 0])
    assert np.all(nlist.point_indices == nlist[:, 1])

Note that the ``query_points`` are always in the first column, while the ``points`` are in the second column.
:class:`freud.locality.NeighborList` objects also store other properties; for instance, they may assign different weights to different bonds.
This feature can be used by, for example, :class:`freud.order.Steinhardt`, which is typically used for calculating `Steinhardt order parameters <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784>`_, a standard tool for characterizing crystalline order.
When provided appropriately weighted neighbors, however, the class instead computes `Minkowski structure metrics <https://iopscience.iop.org/article/10.1088/1367-2630/15/8/083028/meta>`_, which are much more sensitive measures that can differentiate a wider array of crystal structures.
