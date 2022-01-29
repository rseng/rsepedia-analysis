# <img src="https://github.com/gridap/Gridap.jl/blob/master/images/color-text.png" width="250" title="Gridap logo">


| **Documentation** |
|:------------ |
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/Gridap.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/Gridap.jl/dev) |
|**Build Status** |
| [![Build Status](https://github.com/gridap/Gridap.jl/workflows/CI/badge.svg?branch=master)](https://github.com/gridap/Gridap.jl/actions?query=workflow%3ACI) [![Codecov](https://codecov.io/gh/gridap/Gridap.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/Gridap.jl) |
| **Community** |
| [![Join the chat at https://gitter.im/Gridap-jl/community](https://badges.gitter.im/Gridap-jl/community.svg)](https://gitter.im/Gridap-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) |
| **Citation** |
| [![DOI](https://joss.theoj.org/papers/10.21105/joss.02520/status.svg)](https://doi.org/10.21105/joss.02520) |

## What

Gridap provides a set of tools for the grid-based approximation of partial differential equations (PDEs) written in the
[Julia programming language](https://julialang.org/). The library currently supports linear and nonlinear PDE systems for scalar and vector fields, single and multi-field problems, conforming and nonconforming finite element (FE) discretizations, on structured and unstructured meshes of simplices and n-cubes. Gridap is extensible and modular. One can implement new FE spaces, new reference elements, use external mesh generators, linear solvers, post-processing tools, etc. See, e.g., the list of available [Gridap plugins](https://github.com/gridap/Gridap.jl#plugins).

Gridap has a very expressive API allowing one to solve complex PDEs with very few lines of code. The user can write the underlying weak form with a syntax almost 1:1 to the mathematical notation, and Gridap generates an efficient FE assembly loop automatically by leveraging the Julia JIT compiler. For instance, the weak form for an interior penalty DG method for the Poisson equation can be simply specified as: 
```julia
a(u,v) =
  ∫( ∇(v)⋅∇(u) )*dΩ +
  ∫( (γ/h)*v*u - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )*dΓ +
  ∫(
    (γ/h)*jump(v*n_Λ)⋅jump(u*n_Λ) -
    jump(v*n_Λ)⋅mean(∇(u)) -
    mean(∇(v))⋅jump(u*n_Λ)
    )*dΛ

l(v) =
  ∫( v*f )*dΩ +
  ∫( (γ/h)*v*u - (n_Γ⋅∇(v))*u )*dΓ
```
See the complete code [here](https://github.com/gridap/Gridap.jl/blob/master/test/GridapTests/PoissonDGTests.jl). As an example for multi-field PDEs, this is how the weak form for the Stokes equation with Neumann boundary conditions can be specified:
```julia
a((u,p),(v,q)) =
  ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )*dΩ

l((v,q)) =
  ∫( v⋅f + q*g )*dΩ +
  ∫( v⋅(n_Γ⋅∇u) - (n_Γ⋅v)*p )*dΓ
```
See the complete code [here](https://github.com/gridap/Gridap.jl/blob/master/test/GridapTests/StokesTaylorHoodTests.jl).


## Documentation

- [**STABLE**](https://gridap.github.io/Gridap.jl/stable) &mdash; **Documentation for the most recently tagged version of Gridap.jl.**
- [**DEVEL**](https://gridap.github.io/Gridap.jl/dev) &mdash; *Documentation for the in-development version of Gridap.*

## Tutorials

A hands-on user-guide to the library is available as a set of [tutorials](https://github.com/gridap/Tutorials). They are available as Jupyter notebooks and html pages.

## Installation

Gridap is a registered package in the official [Julia package registry](https://github.com/JuliaRegistries/General).  Thus, the installation of Gridap is straight forward using the [Julia's package manager](https://julialang.github.io/Pkg.jl/v1/). Open the Julia REPL, type `]` to enter package mode, and install as follows
```julia
pkg> add Gridap
```

## Plugins
- [GridapDistributed](https://github.com/gridap/GridapDistributed.jl) Distributed-memory extension of Gridap.
- [GridapEmbedded](https://github.com/gridap/GridapEmbedded.jl) Embedded finite elements in Julia.
- [GridapGmsh](https://github.com/gridap/GridapGmsh.jl) Generate a FE mesh with [GMSH](www.gmsh.info) and use it in Gridap.
- [GridapMakie](https://github.com/gridap/GridapMakie.jl) Makie plotting recipies for Gridap.
- [GridapPardiso](https://github.com/gridap/GridapPardiso.jl) Use the [Intel Pardiso MKL direct sparse solver](https://software.intel.com/en-us/mkl-developer-reference-fortran-intel-mkl-pardiso-parallel-direct-sparse-solver-interface) in Gridap.
- [GridapPETSc](https://github.com/gridap/GridapPETSc.jl) Use [PETSc](https://petsc.org/) linear and nonlinear solvers in Gridap.
- [GridapODEs](https://github.com/gridap/GridapODEs.jl) Gridap support for time-dependent PDEs.


## Examples

These are some popular PDEs solved with the Gridap library. Examples taken from the [Gridap Tutorials](https://github.com/gridap/Tutorials).

| ![](https://gridap.github.io/Tutorials/dev/assets/poisson/fig_uh.png)   |  ![](https://gridap.github.io/Tutorials/dev/assets/elasticity/disp_ux_40.png) | ![](https://gridap.github.io/Tutorials/dev/assets/hyperelasticity/neo_hook_3d.png)  | ![](https://gridap.github.io/Tutorials/dev/assets/p_laplacian/sol-plap.png)  |
|:-------------:|:-------------:|:-----:|:----:|
| [Poisson equation](https://gridap.github.io/Tutorials/dev/pages/t001_poisson/) |  [Linear elasticity](https://gridap.github.io/Tutorials/dev/pages/t003_elasticity/) |  [Hyper-elasticity](https://gridap.github.io/Tutorials/dev/pages/t005_hyperelasticity/)  | [p-Laplacian](https://gridap.github.io/Tutorials/dev/pages/t004_p_laplacian/)   |
| ![](https://gridap.github.io/Tutorials/dev/assets/dg_discretization/jump_u.png) | ![](https://gridap.github.io/Tutorials/dev/assets/darcy/darcy_results.png) |![](https://gridap.github.io/Tutorials/dev/assets/inc_navier_stokes/ins_solution.png) | ![](https://gridap.github.io/Tutorials/dev/assets/isotropic_damage/damage_end.png) |
| [Poisson eq. with DG](https://gridap.github.io/Tutorials/dev/pages/t006_dg_discretization/)  |  [Darcy eq. with RT](https://gridap.github.io/Tutorials/dev/pages/t007_darcy/)  |  [Incompressible Navier-Stokes](https://gridap.github.io/Tutorials/dev/pages/t008_inc_navier_stokes/)  | [Isotropic damage](https://gridap.github.io/Tutorials/dev/pages/t010_isotropic_damage/)  |

## Known issues

Since Julia 1.6 ownwards we have noticed large first call latencies of Gridap.jl codes with the default compiler optimization level (i.e., `-O2`). 
In general, while developing code, but specially if you are noting high first call latencies, we recommend to run `julia` with the `-O1` flag. For production runs use `-O2` or `-O3`.  

 ## Gridap community

You can ask questions and interact with the Gridap community on the Julia Slack channel #gridap (see [here](https://julialang.org/slack/) how to join). or our [gitter](https://gitter.im/Gridap-jl/community).

## Contributing to Gridap

Gridap is a collaborative project open to contributions. If you want to contribute, please take into account:

  - Before opening a PR with a significant contribution, contact the project administrators, e.g., by writing a message in [our gitter chat](https://gitter.im/Gridap-jl/community) or by opening an issue describing what you are willing to implement. Wait for feed-back.
  - Carefully read and follow the instructions in the [CONTRIBUTING.md](https://github.com/gridap/Gridap.jl/blob/master/CONTRIBUTING.md) file.
  - Carefully read and follow the instructions in the [CODE_OF_CONDUCT.md](https://github.com/gridap/Gridap.jl/blob/master/CODE_OF_CONDUCT.md) file.
  - Open a PR with your contribution.

Want to help? We have a number of [issues waiting for help](https://github.com/gridap/Gridap.jl/labels/help%20wanted). You can start contributing to the Gridap project by solving some of those issues.


## How to cite Gridap

In order to give credit to the `Gridap` contributors, we simply ask you to cite the reference below in any publication in which you have made use of the `Gridap` project. If you are using other `Gridap` sub-packages, please cite them as indicated in their repositories.

```
@article{Badia2020,
  doi = {10.21105/joss.02520},
  url = {https://doi.org/10.21105/joss.02520},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2520},
  author = {Santiago Badia and Francesc Verdugo},
  title = {Gridap: An extensible Finite Element toolbox in Julia},
  journal = {Journal of Open Source Software}
}
```

## Contact


Please, contact the project administrators, [Santiago Badia](mailto:santiago.badia@monash.edu) and [Francesc Verdugo](mailto:fverdugo@cimne.upc.edu), for further questions about licenses and terms of use.

# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Support for periodic conditions in `CartesianDiscreteModel`s built from `cmin`, `cmax`. Since PR [#738](https://github.com/gridap/Gridap.jl/pull/738).

## [0.17.7]- 2021-12-03

### Fixed
- Improving compile times by removing unnecessary `@inline` statements. Since PR [#726](https://github.com/gridap/Gridap.jl/pull/726).

### Added
- CellFE constructor now gets optional arguments and pass them down. Since PR [#728](https://github.com/gridap/Gridap.jl/pull/728).


## [0.17.6] - 2021-12-01

### Added
- Implemented `Base.unaliascopy(A::SubVector)`. Since PR [#715](https://github.com/gridap/Gridap.jl/pull/715).

### Fixed
- Bugfix in  `Base.view(glue::FaceToFaceGlue,ids::AbstractArray)`. Since PR [#724](https://github.com/gridap/Gridap.jl/pull/724).


## [0.17.5] - 2021-11-08

### Added
- Hiding the creation of `FESolver` and `LinearFESolver` from user code. Since PR [#705](https://github.com/gridap/Gridap.jl/pull/705).

## [0.17.4] - 2021-10-27

### Fixed
 - Using implementation of `pvtk_grid` provided in WriteVTK.  Since PR [#699](https://github.com/gridap/Gridap.jl/pull/699).

## [0.17.3] - 2021-10-27

### Fixed
 - Adding a newer version of WriteVTK in the [compat] section.  Since PR [#698](https://github.com/gridap/Gridap.jl/pull/698).

## [0.17.2] - 2021-10-26

### Fixed
- AD with multi-field residuals with different num dofs per field. Since PR [#687](https://github.com/gridap/Gridap.jl/pull/687).

## [0.17.1] - 2021-10-26

### Fixed

- Laplacian `Δ` operator on unstructured linear grids for quantities defined in the reference space (i.e. shape funcitons in standard FEM). Since PR [#691](https://github.com/gridap/Gridap.jl/pull/691).
- Laplacian `Δ` operator on triangulations using `GridView` (e.g., when interpolating functions in a sub-domain or on the boundary). Since PR [#691](https://github.com/gridap/Gridap.jl/pull/691).
- Fixed typo in , function `solve! of `LinearSolvers.jl`. Since PR [#692](https://github.com/gridap/Gridap.jl/pull/692).

## [0.17.0] - 2021-10-22

### Added

- Aliases `Interior`, `Boundary`, `Skeleton`, and `Interface` for the `Triangulation`, `BoundaryTriangulation`, `SkeletonTriangulation`, and `InterfaceTriangulation` constructors. Since PR [#662](https://github.com/gridap/Gridap.jl/pull/662).
- Function `create_pvtk_file` for exporting results in `pvtu` format. Since PR [#685](https://github.com/gridap/Gridap.jl/pull/685).

### Changed

- Major refactoring in the `Triangulation` interface to properly support the solution of PDEs defined on domains of different dimension. The major change from the user perspective is that `Triangulation` objects can be used both to integrate the weak form (as before) but also to define FE spaces (except for unfitted triangulations obviously). It is still possible to define FE spaces from `DiscreteModels`, but it is safer and more idiomatic (closer to the math notation) to use `Triangulation` objects from now on. Since PR [#662](https://github.com/gridap/Gridap.jl/pull/662).
- Changes in assembly interface to allow optimization when assembling matrices and vectors simultaneously. Since PR [#685](https://github.com/gridap/Gridap.jl/pull/685).

### Removed

- `BoundaryDiscreteModel`, `RestrictedDiscreteMdeol`, `RestrictedTriangulation`, `TriangulationStyle`, `BackgroundTriangulation`, `SubTriangulation`, `get_cell_to_bgcell`, `get_cell_ref_map`, `get_background_triangulation`, and `have_compatible_domains`. Since PR [#662](https://github.com/gridap/Gridap.jl/pull/662). 
- Functions `scale_entries!` and `fill_entries!`. Replaced by Julia functions `LinearAlgebra.rmul!` and `LinearAlgebra.fillstored!`. Since PR [#680](https://github.com/gridap/Gridap.jl/pull/680).

## [0.16.5] - 2021-09-08

### Added
- Implemented DIV operator for FE functions in RT space. Since PR [#650](https://github.com/gridap/Gridap.jl/pull/650).
- `GenericAssemblyStrategy`. Since PR [#655](https://github.com/gridap/Gridap.jl/pull/655).
- Additional high level API assembly functions. Since PR [#652](https://github.com/gridap/Gridap.jl/pull/652).

### Fixed
- Bug related with the release of ChainRulesCore version 1.3.1. Since [#654](https://github.com/gridap/Gridap.jl/pull/654).
- Inheritance relationship for DiscreteModelPortion. Since PR [#645](https://github.com/gridap/Gridap.jl/pull/645).
- Optimization to RT FEs. Since PR [#638](https://github.com/gridap/Gridap.jl/pull/638).
- Bug in boundary discrete model. Since PR [#651](https://github.com/gridap/Gridap.jl/pull/651).

## [0.16.4] - 2021-08-17

### Added
- Trait to CellQuadrature to support the evaluation of integrals on the reference domain. Since PR [#636](https://github.com/gridap/Gridap.jl/pull/636).
- Type `Interpolable` allowing to interpolate data from two independent meshes. Since PR [#632](https://github.com/gridap/Gridap.jl/pull/632).

## [0.16.3] - 2021-06-28

### Fixed
- Deactivating optimizations related with `MemoArray` since they are not reliable. Since PR [#624](https://github.com/gridap/Gridap.jl/pull/624).
- Bug related with `ArrayBlock`. Since PR [#623](https://github.com/gridap/Gridap.jl/pull/623).

## [0.16.2] - 2021-06-21

### Fixed
- Bug related with boundary integration caused by some optimization introduced in v0.16. Fixed via PR [#616](https://github.com/gridap/Gridap.jl/pull/616).

## [0.16.1] - 2021-06-04

### Fixed
- Bug for 1st order FE spaces in combination of 1st order models with periodic BCs. Since PR [#611](https://github.com/gridap/Gridap.jl/pull/611).

## [0.16.0] - 2021-06-04

### Added
- User API to select specific quadrature rules. Since PR [#578](https://github.com/gridap/Gridap.jl/pull/578).
- Experimental support for mixed dimensional PDEs. Since PR [#567](https://github.com/gridap/Gridap.jl/pull/567).
- Added `get_cell_dof_basis(model,cell_reffes,::Conformity)` and `get_cell_shapefuns(model,cell_reffes,::Conformity)`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Implemented `get_cell_dof_basis` and `get_cell_shapefuns` for global RT FE spaces in a new file `DivConformingFESpaces.jl`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Added support to allow evaluation of FE functions at arbitrary points. Since PR [#523](https://github.com/gridap/Gridap.jl/pull/523).
- Implemented `compute_cell_points_from_vector_of_points` to build `CellPoint` from a vector of points. Since PR [#523](https://github.com/gridap/Gridap.jl/pull/523).

### Changed
- Major refactoring in the handling of blocks (e.g. in multi-field and skeleton terms). The new code follows a much more simple approach based in the new type `ArrayBlock`. Since PR [#583](https://github.com/gridap/Gridap.jl/pull/583).
- The default quadrature rule for tets has changed. Since PR [#578](https://github.com/gridap/Gridap.jl/pull/578).
- Refactoring in `SparseMatrixAssembler` to make it more extensible and efficient. Since PR [#568](https://github.com/gridap/Gridap.jl/pull/568).
- Renamed `get_free_values` -> `get_free_dof_values`. Since PR [#567](https://github.com/gridap/Gridap.jl/pull/567).
- Renamed `get_dirichlet_values` -> `get_dirichlet_dof_values`. Since PR [#606](https://github.com/gridap/Gridap.jl/pull/606).
- Renamed `object` -> `value` the variable in `ConstantField`. Since PR [#606](https://github.com/gridap/Gridap.jl/pull/606).
- Miscellaneous changes in the FE assembly to allow the solution of mixed dimensional problems. Since PR [#567](https://github.com/gridap/Gridap.jl/pull/567).
- Renamed `get_cell_shapefuns` by `get_fe_basis`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Renamed `get_cell_shapefuns_trial` by `get_trial_fe_basis`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Renamed `get_cell_dof_basis` by `get_fe_dof_basis`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Removed `conformity` optional keyword argument from `FESpace(::DiscreteModel,::CellFE; kwargs...)` constructor. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).
- Replaced `CellFE(::AbstractArray{<:Field},::AbstractArray{<:ReferenceFE})` by `CellFE(::DiscreteModel,::AbstractArray{<:ReferenceFE},::Conformity)`. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).

### Removed
- All code associated with with `BlockArrayCoo`. Since PR [#583](https://github.com/gridap/Gridap.jl/pull/583).
- Module `Gridap.Integration` has been deleted and its contents have been merged into `Gridap.ReferenceFEs` module.
- Types `SparseMatrixCSR` and `SymSparseMatrixCSR` have been moved to the registered package [`SparseMatricesCSR`](https://github.com/gridap/SparseMatricesCSR.jl). To use them simply add `SparseMatricesCSR` into your environment and type `using SparseMatricesCSR`. Since  Since PR [#568](https://github.com/gridap/Gridap.jl/pull/568).
- Removed `PushForwardMap` and all code depending upon it. Since PR [#579](https://github.com/gridap/Gridap.jl/pull/579).

## [0.15.5] - 2021-05-25

### Added
- Differential operators `(∇+k)(u)`, `(∇+k)⋅u`, `(∇+k)×u`, `(∇+k)⊗u`, and `u⊗(∇+k)` for some `u::CellField` and `k::VectorValue`. Since PR [#597](https://github.com/gridap/Gridap.jl/pull/597).
- Definition of `u.*v` between instances of vector-valued `CellField` objects `u` and `v`. Also differential operators `∇.*u`  and `(∇+k).*u`. Since PR [#597](https://github.com/gridap/Gridap.jl/pull/597).

## [0.15.4] - 2021-03-29

### Fixed
- Bug in `CartesianDiscreteModel` with periodic boundary conditions that shows up in Julia 1.6 but not in Julia 1.5. Since commit [da005cf](https://github.com/gridap/Gridap.jl/commit/da005cf4cde68617f92d76744e307798ef7e8340).

## [0.15.3] - 2021-03-16

### Added
- `get_cell_map` now returns array of `AffineMap` for linear grids of simplices. Needed to compute Laplacian operator, inverse maps etc. Since PR [#553](https://github.com/gridap/Gridap.jl/pull/553).

### Fixed
- Bug in `print_op_tree`. Since PR [#563](https://github.com/gridap/Gridap.jl/pull/563)

## [0.15.2] - 2021-03-08

### Added
- Method `inverse_map` for `AffineMap`. Since PR [#552](https://github.com/gridap/Gridap.jl/pull/552).
- Method `get_cell_points` for `CellDof`. Since PR [#551](https://github.com/gridap/Gridap.jl/pull/551).
- Evaluation of `MonomialBasis` objects at a single point. Since PR [#550](https://github.com/gridap/Gridap.jl/pull/550).
- `rand` function for `MultiValue` objects. Since PR [#530](https://github.com/gridap/Gridap.jl/pull/530).

### Fixed
- Bug in `return_value` for `Broadcasting(∇∇)`. Since PR [#554](https://github.com/gridap/Gridap.jl/pull/554).
- Bug in `dot` for third order tensors. Since PR [#544](https://github.com/gridap/Gridap.jl/pull/544).

## [0.15.1] - 2021-01-22

### Added
- Added support for Hessian and Laplacian operators. Only implemented for Finite Elements with an `AffineMap`. Since  PR [#514](https://github.com/gridap/Gridap.jl/pull/514).

### Fixed
- Bug in `RestrictedDiscreteModel` for periodic boundary conditions. Since  PR [#517](https://github.com/gridap/Gridap.jl/pull/517).
- Bug in `sum(a::LazyArray)` when  `eltype(a) <: AbstractArray`. Since PR [#513](https://github.com/gridap/Gridap.jl/pull/513).

## [0.15.0] - 2020-12-14

This version is a major (backwards-incompatible) refactoring of the project which is not summarized here for the sake of brevity. Most of the functionality of v0.14.0 is available in v0.15.0, but possibly with a significantly different API. See [here](https://github.com/gridap/Tutorials/compare/v0.14.0...v0.15.0) the changes in the sources of the Gridap Tutorials between versions 0.14.0 and 0.15.0 to effectively see the major changes in the API.

## [0.14.2] - 2020-11-24

### Added
- Added additional tensor operations and new double contraction notation `⋅²`. Implemented a `zero` constructor for `ThirdOrderTensorValues` to allow integration of 3-tensors. Since PR [#415](https://github.com/gridap/Gridap.jl/pull/415/).

### Fixed
 - Bug-fix for 32-bit Julia: Replace all occurences of Int64 by Int. Since PR [#445](https://github.com/gridap/Gridap.jl/pull/445).
 - Bug-fix for 32-bit Julia. Using inttype=Int keyword argument for JSON parsing. Since PR [#456](https://github.com/gridap/Gridap.jl/pull/456).

## [0.14.1] - 2020-09-17

### Added
 - Added VectorWithEntryInserted and VectorWithEntryRemoved. Since PR [#401](https://github.com/gridap/Gridap.jl/pull/401/).
 - Added missing get_constant_approach() getter to FESpaceWithConstantFixed. Since PR [#409](https://github.com/gridap/Gridap.jl/pull/409).

### Deprecated
 - The name FESpaceWithLastDofRemoved has been deprecated in favor of its generalization FESpaceWithConstantFixed. Since PR [#396](https://github.com/gridap/Gridap.jl/pull/396) and PR [#404](https://github.com/gridap/Gridap.jl/pull/404).

## [0.14.0] - 2020-08-27

### Removed
 - Support for Julia v1.0. Now, the minimum supported is Julia v1.3. Since PR [#376](https://github.com/gridap/Gridap.jl/pull/376/).

### Changed
 - Major refactoring associated with the handling of elemental matrices and vectors in multi-field computations and also on the skeleton. Since PR [#376](https://github.com/gridap/Gridap.jl/pull/376/).
 - First and second argument switch in `update_state_variables!` in order to have function-first style. Since PR [#376](https://github.com/gridap/Gridap.jl/pull/376/).
 - Table struct has been generalized such that data and ptrs arrays can be of an arbitrary type extending AbstractArray. Since PR [#310](https://github.com/gridap/Gridap.jl/pull/310/)
 - `interpolate, interpolate!, interpolate_dirichlet...` switched argument order to function first style. For instance `interpolate(u, V)` instead of `interpolate(V, u)`

### Added
 - Allowing the construction of an `HomogeneousTrialFESpace` from a `TrialFESpace`. Since PR [#384](https://github.com/gridap/Gridap.jl/pull/384).
 - Support for automatic differentiation of residuals and Jacobians in multi-field computations since PR [#383](https://github.com/gridap/Gridap.jl/pull/383/).
 - New `FilterKernel` since PR [#379](https://github.com/gridap/Gridap.jl/pull/379/).

### Fixed
 - Bug associated with boundary triangulation in 1D discrete models. Since PR [#393](https://github.com/gridap/Gridap.jl/pull/393).

## [0.13.4] - 2020-08-23

### Added
  - New `FilteredCellArray` since PR [#372](https://github.com/gridap/Gridap.jl/pull/372/).

## [0.13.3] - 2020-08-12

### Added
  - `Visualization.visualization_data` function that makes it easier to bring fields into
  visualization library friendly formats. Since PR [#354](https://github.com/gridap/Gridap.jl/pull/354).
  - Gradient of a product binary operation (`*`) between a scalar and a field. Since PR [#340](https://github.com/gridap/Gridap.jl/pull/340).

## [0.13.2] - 2020-07-31

### Added
  - Automatic differentiation of the Jacobian from a given residual and the Jacobian and the residual from a given energy. Not working at this moment on the Skeleton nor for multi-field (WIP), but yes for other cases.
  Now, the user can omit `jac` from `FETerm(res,jac,trian,quad)`, i.e. `FETerm(res,trian,quad)` and the Jacobian will be automatically generated. In addition, the user can write `FEEnergy(ener,trian,quad)` for a given `ener(uh)` function
  and the residual and the Jacobian will be automatically generated. Since PR [#338](https://github.com/gridap/Gridap.jl/pull/338/).

## [0.13.1] - 2020-07-24

### Fixed
  - Bugs associated with the degenerated case of 0-length arrays. Since PR [#331](https://github.com/gridap/Gridap.jl/pull/331/) and [#332](https://github.com/gridap/Gridap.jl/pull/332/).

## [0.13.0] - 2020-07-23

### Added
  - Automatic differentiation for symmetric gradient, i.e. `ε(u)` for a given vector-valued function `u`. Since PR [#327](https://github.com/gridap/Gridap.jl/pull/327/).
  - Added missing SparseMatrixAssembler constructor for MultiFieldFESpaces. Since PR [#320](https://github.com/gridap/Gridap.jl/pull/320/).
  - kw-argument `space` to `LagrangianRefFE` constructor in order to select the type of underlying polynomial space, i.e., `:Q`, `:S`, or `:P`. Since PR [#321](https://github.com/gridap/Gridap.jl/pull/321).

### Changed
  - The meaning of `inward/outward` has slightly changed for `SkeletonCellBasis` objects. Now, by accessing to these properties a `ReducedSkeletonCellBasis` is returned, which allows to use the result in a more flexible way (in particular, the result can be used in a similar way than the result of `jump` or `mean`). Since PR [#317](https://github.com/gridap/Gridap.jl/pull/317).
  - Major refactoring in `ReferenceFEs` module. Since PR [#319](https://github.com/gridap/Gridap.jl/pull/319) and [#321](https://github.com/gridap/Gridap.jl/pull/321). In particular:
    - `NodalReferenceFE` has been replaced by a new abstract type `LagrangianRefFE`.
    - `GenericNodalCartesianRefFE` has been replaced by `GenericLagrangianRefFE`.

### Removed
  - Removals associated with the `ReferenceFEs` refactoring in PR [#319](https://github.com/gridap/Gridap.jl/pull/319):
    - Removed `QDiscRefFE` constructor. Use a standard `LagrangianRefFE` and `L2Conformity` instead.
    - Removed `PDiscRefFE` constructor. Use `LagrangianRefFE` constructor with the kw-argument `space=:P`.
    - Removed `CDLagrangianRefFE` constructor. Use a standard `LagrangianRefFE` and `CDConformity` instead.
    - Removed fields `face_own_dofs` and `face_own_dof_permutations` from `GenericRefFE`.
    - Removed struct `DiscRefFE`.

### Fixed
  - Better handling of FE terms defined on empty triangulations. Since PR [#329](https://github.com/gridap/Gridap.jl/pull/329).
  - Replaced `+=` by `add_entry!`. Since PR [#316](https://github.com/gridap/Gridap.jl/pull/316).
  - Minor fix to let Vtk.jl support changes in Vtk 1.7.X versus 1.6.X. Since PR [#324](https://github.com/gridap/Gridap.jl/pull/324).

## [0.12.0] - 2020-07-07

### Added

  - Added `SkeletonTriangulation` constructor in order to integrate, where a given interpolation is discontinuous. Since PR [#304](https://github.com/gridap/Gridap.jl/pull/304).
  - New `ConformingFESpace` constructor. Since PR [#293](https://github.com/gridap/Gridap.jl/pull/293).
  - Added `QDiscRefFE` constructor for `DiscRefFE`. Since PR [#293](https://github.com/gridap/Gridap.jl/pull/293).
  - New `FESpace` constructor that takes an instance of `ReferenceFE`. Since PR [#294](https://github.com/gridap/Gridap.jl/pull/294).
  - New `FESpace` constructor that takes an instance of `Conformity`. Since PR [#311](https://github.com/gridap/Gridap.jl/pull/311).
  - New `CDLagrangianRefFE` struct, that provides a Lagrangian reference FE with different conformity per direction. Since PR [#299](https://github.com/gridap/Gridap.jl/pull/299).
  - New `FESpace` method that takes a model and a `RefFE`. Since PR [#299](https://github.com/gridap/Gridap.jl/pull/299).
  - Possibility to have 0 order in `DISC` directions of a `CDLagrangianRefFE`. Since PR [#308](https://github.com/gridap/Gridap.jl/pull/308).
  - Added setindex! method for Reindexed. Since PR [#309](https://github.com/gridap/Gridap.jl/pull/309).


### Changed

  - Changed the interfaces of `ReferenceFE` and `NodalReferenceFE` in relation of DOF ownership. Now function `get_face_own_dofs` and related ones are parametrized by a `Conformity` object. Since PR [#311](https://github.com/gridap/Gridap.jl/pull/311).
  - The constructors `GenericRefFE`, `GenericNodalCartesianRefFE`, and `compute_conforming_cell_dofs` take an extra argument of type `Conformity`. Since PR [#311](https://github.com/gridap/Gridap.jl/pull/311).
  - Renamed `PDiscRefFE` -> `DiscRefFE` struct keeping the name for constructor. Since PR [#293](https://github.com/gridap/Gridap.jl/pull/293).
  - One of the `GradConformingFESpace` methods now more general `ConformingFESpace`. Since PR [#293](https://github.com/gridap/Gridap.jl/pull/293).
  - `DivConformingFESpace` and `CurlConformingFESpace` constructors eliminated. Since PR [#293](https://github.com/gridap/Gridap.jl/pull/293).
  - Extend table to support arbitrary vector types. Since PR [#310](https://github.com/gridap/Gridap.jl/pull/310).

### Fixed

  - Construction of `VectorValue`, `TensorValue`, et al. objects from non-homogeneous arguments.  This solves some problems associated with automatic differentiation. Since PR [#298](https://github.com/gridap/Gridap.jl/pull/298).
  - `CDLagrangianRefFE` node ordering. Since PR [#305](https://github.com/gridap/Gridap.jl/pull/305).

## [0.11.2] - 2020-06-22

### Added

  - Method `solve!(x,ls,op::AffineOperator,cache::Nothing,newmatrix)`. Since PR [#288](https://github.com/gridap/Gridap.jl/pull/288).

### Fixed

  - Bug related with `WriteVTK` version 1.7. Fixed via PR [#287](https://github.com/gridap/Gridap.jl/pull/287).
  - Bug in outer constructor of Table{...} for input arrays of abstract type. Fixed via PR [#285](https://github.com/gridap/Gridap.jl/pull/285).

## [0.11.1] - 2020-06-19

### Fixed

  - Bug in the handling of caches in `NLSolver`. Fixed via PR [#283](https://github.com/gridap/Gridap.jl/pull/283).
  - Bug that showed up when interpolating a FE function defined on an
  `ExtendedFESpace` onto a non-extended `FESpace`. Fixed via PR [#282](https://github.com/gridap/Gridap.jl/pull/282).

## [0.11.0] - 2020-06-16

### Added

  - Operator `⊙` (\odot) as an alias of `inner`. Since PR [#239](https://github.com/gridap/Gridap.jl/pull/239).
  - Operator `⊗` (\otimes) as an alias of `outer`. Since PR [#239](https://github.com/gridap/Gridap.jl/pull/239).
  - Support for (symmetric) 4th order tensors. Since PR [#239](https://github.com/gridap/Gridap.jl/pull/239).
  - Optimizations for symmetric 2nd order tensors. Since PR [#239](https://github.com/gridap/Gridap.jl/pull/239).
  - Methods for `cross` function (aka `×` (\times)) to operate with `VectorValues`. Since PR [#280](https://github.com/gridap/Gridap.jl/pull/280).
  - Interpolation is now supported also for multifield spaces. Since PR [#279](https://github.com/gridap/Gridap.jl/pull/279).

### Changed

  - Major refactoring in the module `Gridap.TensorValues`.
  Since PR [#239](https://github.com/gridap/Gridap.jl/pull/239).
   **The following changes are likely to affect all users:**
    - The operator `*` is not allowed for expressing the dot product anymore. Use `LinearAlgebra.dot`
  function aka `⋅` (\cdot).
    - The syntax `∇*u` is not allowed anymore.  Use `∇⋅u` instead.
    - Gridap re-exports `dot`, `⋅`, and other names from LinearAlbegra that are used
  often in Gridap code.
    - Function `n_components` is renamed to `num_components`.
  - The `SingleFieldFESpace` interface has changed. The function `gather_free_and_dirichlet_values!`
  has been added as mandatory for all FE space implementations and the old function `gather_free_and_dirichlet_values`
  is now optional. Since PR [#279](https://github.com/gridap/Gridap.jl/pull/279).

## [0.10.4] - 2020-06-8

### Added

- Functions `create_vtk_file` and `createvtk`. Since PR [#273](https://github.com/gridap/Gridap.jl/pull/273).

## [0.10.3] - 2020-05-29

### Added

 - Function `print_op_tree` to visualize lazy operation trees. Since PR [#270](https://github.com/gridap/Gridap.jl/pull/270).
 - Exported `apply` and `reindex` from `Gridap` top level. Since PR [#270](https://github.com/gridap/Gridap.jl/pull/270).
 - Extended support of `CartesianDiscreteModel` to models with periodic boundary conditions.
 PR [#266](https://github.com/gridap/Gridap.jl/pull/266).

### Deprecated
 - Optional argument `map` for CartesianDescriptor converted to a key-word argument. Since PR [#266](https://github.com/gridap/Gridap.jl/pull/266).

### Fixed
 - Fixed some methods of the `sparsecsr` generic function. Since PR [#262](https://github.com/gridap/Gridap.jl/pull/262).
 - Fixed BUG in `findnz` function for `SparseMatrixCSR`. Since PR [#264](https://github.com/gridap/Gridap.jl/pull/264).
 - Fixed `restrict(::AbstractArray,::TriangulationPortion)` for portions of triangulations extending `BoundaryTriangulation`. Since PR [#267](https://github.com/gridap/Gridap.jl/pull/267).

## [0.10.2] - 2020-05-21

### Added

 - New key-word arguments `zeromean_trian` and `zeromean_quad` in the `FESpace` constructor. Since
 PR [#257](https://github.com/gridap/Gridap.jl/pull/257).
 - New method `reindex(::Triangulation,indices)`. Since
 PR [#257](https://github.com/gridap/Gridap.jl/pull/257).
 - New functions `get_face_to_face(::BoundaryTriangulation)` and `get_cell_around(::BoundaryTriangulation)`. Since
 PR [#256](https://github.com/gridap/Gridap.jl/pull/256).

## [0.10.1] - 2020-05-19

### Fixed

  - Added missing implementation of `simplexify(SEGMENT)` and `simplexify(VERTEX)`. Since PR [#252](https://github.com/gridap/Gridap.jl/pull/252).

## [0.10.0] - 2020-05-14

### Added

  - Extended support of `TriangulationPortion` to boundary and skeleton triangulations.  Since PR [#249](https://github.com/gridap/Gridap.jl/pull/249).
  - Added `FESpaceWithLinearConstraints`. Since PR [#247](https://github.com/gridap/Gridap.jl/pull/247).
  - Added inner constructor to `CartesianDiscreteModel` allowing to build a model that represents a subgrid of
    a larger grid. Since PR [#245](https://github.com/gridap/Gridap.jl/pull/245).

### Changed

  - The part associated with the imposition of constraints in the `FESpace` interface has changed slightly. Since PR [#247](https://github.com/gridap/Gridap.jl/pull/247).
  - Simplified the signature of `zero_free_values(::FESpace)`. Since PR [#249](https://github.com/gridap/Gridap.jl/pull/249).
  - Simplified the signature of `zero_initial_guess(op::NonlinearOperator)`. Since PR [#249](https://github.com/gridap/Gridap.jl/pull/249).
  - Major refactoring in the `Assembler` interface.
    **Important change:** Now, assembly-related functions take the data returned by functions like
    `collect_cell_matrix` as it is. Example: the old user code `assemble_matrix(assembler,collect_cell_matrix(du,dv,terms)...)`
    now is written simply as `assemble_matrix(assembler,collect_cell_matrix(du,dv,terms))`, i.e., the unpack of the last argument is not
    used anymore.  In addition, with the new assembler interface, it is possible to customize the assembly process
    via a so-called `AssemblerStrategy` object. Since PR [#249](https://github.com/gridap/Gridap.jl/pull/249).
  - Change the types of the sizes and partition fields of CartesianDescriptor to tuples instead of points.
    Since PR [#246](https://github.com/gridap/Gridap.jl/pull/246).

## [0.9.2] - 2020-04-26

### Added

  - Automatic differentiation of manufactured solutions. Since PR [#236](https://github.com/gridap/Gridap.jl/pull/236).

## [0.9.1] - 2020-04-20

### Added

  - Function `cell_measure`. Since PR [#234](https://github.com/gridap/Gridap.jl/pull/234).

### Fixed

  - Several bugs associated with `ExtendedFESpace`. In particular, we have fixed a bug that showed up when combining `ZeroMeanFESpace` and `ExtendedFESpace`. Since PR [#234](https://github.com/gridap/Gridap.jl/pull/234).

## [0.9.0] - 2020-04-18

### Added

  - Function `HomogeneousTrialFESpace`. Since PR [#226](https://github.com/gridap/Gridap.jl/pull/226).
  - Function `lazy_append` in order to lazily append two objects (implemented for `AbstractVector`, `Triangulation`, and `CellQuadrature`).  Since PR [#220](https://github.com/gridap/Gridap.jl/pull/220).
  - Support for FE spaces with DOFs defined in the physical space. Since PR [#216](https://github.com/gridap/Gridap.jl/pull/216) and [#218](https://github.com/gridap/Gridap.jl/pull/218).

### Changed

  - Replaced `non_linear` -> `nonlinear` and `NonLinear` -> `Nonlinear`. Since PR [#223](https://github.com/gridap/Gridap.jl/pull/223).
  - The `FESpace` interface has slightly changed, mainly the return type of functions `get_cell_basis` and `get_cell_dof_basis.`. Since PR [#216](https://github.com/gridap/Gridap.jl/pull/216) and [#218](https://github.com/gridap/Gridap.jl/pull/218).

### Fixed

- Bug that showed up in multi-field computations when some field had no contribution to the rhs vector. Since [#229](https://github.com/gridap/Gridap.jl/pull/229).
- Bug in gradient operator in the void part of `ExtendedFESpace` objects. Since PR [#219](https://github.com/gridap/Gridap.jl/pull/219).
- Bug in jumps of quantities restricted to `InterfaceTriangulation` objects.  Since PR [#215](https://github.com/gridap/Gridap.jl/pull/215).

## [0.8.0] - 2020-03-17

### Added

- Support for surface-coupled multi-physics. See [`SurfaceCouplingTests.jl`](https://github.com/gridap/Gridap.jl/blob/master/test/GridapTests/SurfaceCouplingTests.jl) for further details. Since PR [#209](https://github.com/gridap/Gridap.jl/pull/209).
- Support for constitutive laws with state / historical variables. See [`IsotropicDamageTests.jl`](https://github.com/gridap/Gridap.jl/blob/master/test/GridapTests/IsotropicDamageTests.jl) for further details. Since PR [#208](https://github.com/gridap/Gridap.jl/pull/208).
- Curl-conforming reference FE `NedelecRefFE` and corresponding FE space constructor since PR [#199](https://github.com/gridap/Gridap.jl/pull/199).
- New constructors `AffineFETermFromCellMatVec` and `FETermFromCellJacRes` that provides full control in the definition of cell matrices and vectors. Since PR [#191](https://github.com/gridap/Gridap.jl/pull/191).
- Support for simultaneous integration of matrices and vectors. Since PR [#191](https://github.com/gridap/Gridap.jl/pull/191).

### Changed

- Renaming NonLinear to Nonlinear since it is one word and it is not consistent with style
- The definition of interpolation order in Raviart-Thomas and Nédélec reference FEs has changed. Now, the divergence of functions in the Raviart-Thomas space of order `k` belongs to `P_k` or `Q_k` depending on the underlying polytope. Idem for Nédelec, but using the curl instead of the divergence. Since PR [#212](https://github.com/gridap/Gridap.jl/pull/212).

- The order in which test and trial spaces are written in the code has changed and also the other in the arguments of functions defining bi-linear and linear forms, and weak residuals and Jacobians. **This affects everybody that is using Gridap, even the most basic users**. Now, we write the trial space before the test one in all methods taking two spaces in their arguments.  E.g., we have changed `AffineFEOperator(V,U,terms...)` to `AffineFEOperator(U,V,terms...)`, where `U` is the trial and `V` is the test space. For functions defining weak forms, now we have: The new signatures for bi-linear and a linear forms are `a(u,v)`, `l(v)`, where `u` is a trial function and `v` is a test one. For weak Jacobians and residuals `jac(u,du,v)` and `res(u,v)`, where `u` is the (trial) function in which we evaluate these quantities, `du` is the direction in which we evaluate the Jacobian and `v` is a test function. Since PR [#195](https://github.com/gridap/Gridap.jl/pull/195) and PR [#197](https://github.com/gridap/Gridap.jl/pull/197).

- The part related with the application of constraints in the `FESpace` interface has changed. Since PR [#191](https://github.com/gridap/Gridap.jl/pull/191).

### Fixed

- Bug in 1d Cartesian grids. Since PR [#192](https://github.com/gridap/Gridap.jl/pull/192).

## [0.7.1] - 2020-02-18

### Added

- New `DirichletFESpace` that can be used to compute matrices and vectors associated with the Dirichlet DOFs. Since commit [972afcc](https://github.com/gridap/Gridap.jl/commit/972afcc6dd8e024a7daeebd160a9dabe44ff5921)

## [0.7.0] - 2020-02-13

This version is a major refactoring of the project which is not summarized here for the sake of brevity. Most of the functionality of v0.6.0 is available in v0.7.0, but with a possibly slightly different API. See [here](https://github.com/gridap/Tutorials/compare/v0.6.0...v0.7.0) the changes in the sources of the Gridap Tutorials between versions 0.6.0 and 0.7.0 to effectively see the major changes in the API.

## [0.6.0] - 2020-01-24
### Added
- New `GenericRefFE`. Since commit [876ef1e](https://github.com/gridap/Gridap.jl/commit/c3c9010177432b8f07aaecf4a0baa4b93876ef1e)
- New `NedelecRefFE` constructor that generates Nedelec FEs of arbitrary order in 2D and 3D on hex. Since commit [876ef1e](https://github.com/gridap/Gridap.jl/commit/c3c9010177432b8f07aaecf4a0baa4b93876ef1e)
- New keyword argument `map` in the constructor of `CartesianModel`, which allows one to transform the original domain, by defaut [0,1]^d to a new domain through a homeomorphic map. Since commit [30cc4bc](https://github.com/gridap/Gridap.jl/commit/1c51b18f7e21c4915c0b379585dc5d98130cc4bc)
- New keyword argument `map` in the constructor of `CartesianGrid` and a new `map` attribute in this structure, since commit [30cc4bc](https://github.com/gridap/Gridap.jl/commit/1c51b18f7e21c4915c0b379585dc5d98130cc4bc)
- `CartesianGridPoints` has new attribute `map` since commit [30cc4bc](https://github.com/gridap/Gridap.jl/commit/1c51b18f7e21c4915c0b379585dc5d98130cc4bc)
- Added [`SparseMatricesCSR`](https://github.com/gridap/SparseMatricesCSR.jl) support to `SparseMatrixAssembler` and `MultiSparseMatrixAssembler` in [PR #118](https://github.com/gridap/Gridap.jl/pull/118#).

### Changed
- The `RaviartThomasRefFE` has now been replaced by `GenericRefFE`, and the constructor for Raviart-Thomas FEs is called `RTRefFE`. Since commit [876ef1e](https://github.com/gridap/Gridap.jl/commit/c3c9010177432b8f07aaecf4a0baa4b93876ef1e)
- The default map in the `CartesianModel` constructor is [0,1]^d instead of [-1,1]^d, since commit [30cc4bc](https://github.com/gridap/Gridap.jl/commit/1c51b18f7e21c4915c0b379585dc5d98130cc4bc)
- `CartesianGrid` has attribute `map` since commit [30cc4bc](https://github.com/gridap/Gridap.jl/commit/1c51b18f7e21c4915c0b379585dc5d98130cc4bc)

## [0.5.2] - 2019-10-22
### Fixed
- Incompatibility problem with `TensorValues` version 0.3.5. Via commit [3c0682a](https://github.com/gridap/Gridap.jl/commit/3c0682a84250e17086457ca3a90a49d9bce133d0).


## [0.5.1] - 2019-10-03
### Added
- Pretty printing for the types most exposed to users. Since PR [#109](https://github.com/gridap/Gridap.jl/pull/109).
### Fixed
- Bug related to `ZeroMeanFESpace`. Via PR [#111](https://github.com/gridap/Gridap.jl/pull/111).

## [0.5.0] - 2019-09-27

### Added
- Added a high level constructor, namely `FESpace`, to create different types of FE spaces. See issue [#100](https://github.com/gridap/Gridap.jl/issues/100) for more details. Since PR [#102](https://github.com/gridap/Gridap.jl/pull/102).
- Added `ZeroMeanFESpace` to construct FE spaces whose functions have zero mean value. Since PR [#102](https://github.com/gridap/Gridap.jl/pull/102).
- Added Hdiv FE space using Raviart-Thomas reference FEs in [34bfa34](https://github.com/gridap/Gridap.jl/commit/34bfa344efd1bc6a5d3c5993d9639259ed21671a)
- Added the corresponding DOF basis for Raviart-Thomas reference FEs for interpolation of fields [60b9021](https://github.com/gridap/Gridap.jl/commit/60b9021b6d4b5e66a9ec4fe2067aa8278f8ccb52)
- Added an arbitrary order div-conforming Raviart-Thomas reference FE of arbitrary order on quads in commit
[60b9021](https://github.com/gridap/Gridap.jl/commit/60b9021b6d4b5e66a9ec4fe2067aa8278f8ccb52)
- Now, the `tags` argument is optional when constucting `SkeletonTriangulation` and `BoundaryTriangulation` objects from a `DiscreteModel`. Since commit [e6424a3](https://github.com/gridap/Gridap.jl/commit/e6424a304feb38547241e86de07a821e26344a7e).
- Added `mean` operator for quantities restricted to a `SkeletonTriangulation`. Since commit [83798b4](https://github.com/gridap/Gridap.jl/commit/83798b4f38aaf482b968ffd0359eb75c79a21385).
- Extended `NormalVector` to `SkeletonTriangulations`. Since commit [5fb8487](https://github.com/gridap/Gridap.jl/commit/5fb84871128c4388559cc5052d9ff00f0be19462).
- Now, `TrialFESpaces` can be constructed from values instead of functions if the corresponding Dirichlet conditions are constant. Since commit [bae237e](https://github.com/gridap/Gridap.jl/commit/bae237e881db6569622f3559f82bcc3999560526).
- Added the possibility of adding new tags to a `FaceLabels` object via the function `add_tag_from_tags!` and using it to construct FE spaces. Since commit [e9dfac4](https://github.com/gridap/Gridap.jl/commit/e9dfac4489047c0b7e1c62507f4335e9fc76dfd8).
- Added `BackslashSolver` to facilitate the usage in Gridap of the build-in Julia backslash linear solver. Since commit [8e3a9b7](https://github.com/gridap/Gridap.jl/commit/8e3a9b71c64b032c5a572a7ef696f4cbf875190b).
- Added `NLSolver` to facilitat the usage in Gridap of the non-linear solvers available in the official Julia package `NLsolve`. Introduced in commit [e5a933f](https://github.com/gridap/Gridap.jl/commit/e5a933f3093faea221a50bdd796d7f02113ed52c) as `JuliaNLSolver`. Renamed to `NLSolver` in  PR [#108](https://github.com/gridap/Gridap.jl/pull/108).

### Changed
- The Signature of `solve!` for `NumericalSetup` objects. The argument for the system matrix has been removed. The information about the matrix is already in the `NumericalSetup` object. Since commit  [ac212d3](https://github.com/gridap/Gridap.jl/commit/ac212d30205700a919a37f9abf9dac6cbde03e38).
- The signature of `solve!(::FEFunction,::FESolver,::FEOperator)`. Before it was used as `cache = solve!(uh,solver,op)`, now it is used as `uh, cache = solve!(uh,solver,op)`. Since PR [#102](https://github.com/gridap/Gridap.jl/pull/102).
- Previous ConformingFESpace constructor is H1ConformingFESpace since [34bfa34](https://github.com/gridap/Gridap.jl/commit/34bfa344efd1bc6a5d3c5993d9639259ed21671a)

### Deprecated
- `JuliaNLSolver`. Renamed to  `NLSolver`. Since PR [#108](https://github.com/gridap/Gridap.jl/pull/108).
- Key-word argument `order` in `CellQuadrature` constructor.  Renamed to `degree`. Since PR [#108](https://github.com/gridap/Gridap.jl/pull/108).

### Fixed
- Bug in `@law` macro for more than one `FEBasis` arguments. Solved via PR [#104](https://github.com/gridap/Gridap.jl/pull/104).
- Bug in `NonlinearFEOperator` constructor with default assembler in multi-field computations. Solved via PR [#104](https://github.com/gridap/Gridap.jl/pull/104).
- Bug in `NormalVector` for non-Cartesian grids. Solved via PR [#98](https://github.com/gridap/Gridap.jl/pull/98).

## [0.4.0] - 2019-09-07

### Added

- Added support to high order simplicial Lagrangian finite elements. Since commit [cbefe9b](https://github.com/gridap/Gridap.jl/commit/cbefe9bbea83d00e7f6ccbef50396ddc7dc49b80).
- Now the built-in simplicial grids are oriented. Since commit [cbefe9b](https://github.com/gridap/Gridap.jl/commit/cbefe9bbea83d00e7f6ccbef50396ddc7dc49b80).
- Added binary operations between `FEFuntion` and `Number`, and `FEBasis` and `Number`.  Since PR [#88](https://github.com/gridap/Gridap.jl/pull/88).
- Added `PDiscRefFE`, `DiscFESpace`, and `ConstrainedFESpace`. Since PR [#88](https://github.com/gridap/Gridap.jl/pull/88).
- Now its possible to pass a `CellNumer` or an `Array` of numbers into a constitutive law. Usefull to identify which is the material of the current Gauss point in multi-material problems. Since commit [62cb2c3](https://github.com/gridap/Gridap.jl/commit/62cb2c354e2a09c556324a4fe9861329989299f4).
- `LinearFESolver` is now optional for solving a `LinearFEOperator`. Since commit [5c1caa8](https://github.com/gridap/Gridap.jl/commit/5c1caa8c92b260db72f5902e778ec5c0eb88728b).
- `Assembler` is now optional to build `FEOperator` objects. Since commit [b1bf517](https://github.com/gridap/Gridap.jl/commit/b1bf5172955b940f6b3c9d027bd4a839c6486199).
- Binary operations between `Function` and `FEFunction`. Since commit [a7f22f5](https://github.com/gridap/Gridap.jl/commit/a7f22f5ac1f45d9e8f53906472257aa582726e87).
- Extended constructions of `CLagrangianFESpace` and `DLagrangianFESpace`. `diritags` and `dirimasks` are now optional. `diritags` can now be also a vector of `String`. Since commit [776b402](https://github.com/gridap/Gridap.jl/commit/776b40238365f145037fc5e490600bf5b45434ef).
- Added `div`, `curl`, and `trace` operators. Since commit [5a0f322](https://github.com/gridap/Gridap.jl/commit/5a0f322c5b938f12e26e9c0a7c9361aa649e014f).
- Macro `@law` to facilitate the definition of constitutive laws. Since commit [30b67f2](https://github.com/gridap/Gridap.jl/commit/30b67f29009b872944be94486dc4a1b0134a0a60).
- Definition of linear forms `b(v) = inner(v, f)` directly from a function `f`. Since commit [bb42847](https://github.com/gridap/Gridap.jl/commit/bb42847c702a99b9b5f2c2d922fbe4c95b23f646)
- Serialization and de-serialization of `DiscreteModel` objects into and from `json` format. Since PR [#76](https://github.com/gridap/Gridap.jl/pull/76).
- Support for boundary integration (e.g., Neumann BCs) for multi-field computations. Since PR [#75](https://github.com/gridap/Gridap.jl/pull/75).

### Changed

- Signature of `LagrangianRefFE` constructor. Since commit [529c764](https://github.com/gridap/Gridap.jl/commit/529c7646a531db6910a00f04a925dadec3a50b7c).

### Fixed

- Bug in `LinearFETerm` for multi-field computations. Fixed via commit [2b957d1](https://github.com/gridap/Gridap.jl/commit/2b957d1b3a9a9a4396075801d8c837f6aff921c8).
- Bug in `MultiCellArray` constructor. Fixed via commit [bbc3b1c](https://github.com/gridap/Gridap.jl/commit/bbc3b1c91752f8efa978731cb90c6198dc0e5227).
- Bug in binary operations between FEFunction and FEBasis. Fixed via commit [aa49689](https://github.com/gridap/Gridap.jl/commit/aa49689be2a8dc14e052a6409c8348f492b52b3e).

## [0.3.0] - 2019-08-06
### Added
- `CurlGradMonomialBasis` spanning the polynomial space needed for RT elements on n-cubes.
- `CLagrangianFESpace` and `DLagrangianFESpace` types providing an efficient implementation for continuous and discontinuous Lagrangian FE spaces respectivelly. In contrast to `ConfirmingFESpace`, the new types allow to select which are the components that are actually prescribed on the Dirichlet boundary. Since PR [#64](https://github.com/gridap/Gridap.jl/pull/64).
- `simplexify` funciton to convert `Grid` and `DiscreteModel` objects made of n-cubes to the corresponding counterparts made of n-simplices. Since PR [#62](https://github.com/gridap/Gridap.jl/pull/62).
- Duffy transformation based integration for n-simplices of arbitrary dimension. Since PR [#61](https://github.com/gridap/Gridap.jl/pull/61).
- `NormalVector` to construct the outward normal vector to a given `BoundaryTriangulation`. Since PR [#60](https://github.com/gridap/Gridap.jl/pull/60).
- Support for tensor-valued FE computations. Since PR [#57](https://github.com/gridap/Gridap.jl/pull/57).
- Support for integration on the skeleton of the mesh. This includes `SkeletonTriangulation`, an integration mesh for the skeleton, `restrict` function is extended to restrict to the skeleton, `jump` function to compute jumps of `CellFields` and `CellBasis` restricted to the skeleton, extension of `FETerms` to allow integration on the skeleton. See PR [#47](https://github.com/gridap/Gridap.jl/pull/47)
- Support for Robin boundary conditions. Since commit [946054a](https://github.com/gridap/Gridap.jl/commit/946054a028e658afa87c7a7c71e973957a2c4877)
- Support for Neumann boundary conditions. Since commit [4dcd16f](https://github.com/gridap/Gridap.jl/commit/4dcd16fbff9edf66fb66efb748ef01901c20a4aa)
- `FETerm` and `AffineFETerm` abstract types and several concrete implementations. They allow to deal with problems whose weak form has terms integrated over different geometrical entities. `NonlinearFEOperator` and `LinearFEOperator` can be constructed using several terms. Since commit [0f99234](https://github.com/gridap/Gridap.jl/commit/0f99234156dd0174485ca83431de76aa3825584a)
- Extended `Assembler` and `MultiAssembler` to deal with several terms. See issue [#42](https://github.com/gridap/Gridap.jl/issues/42) and PR [#43](https://github.com/gridap/Gridap.jl/pull/43).
- `IdentityCellNumber`, an indexable cell number that simply returns the given index. Also efficient implementation of `reindex` for this type (i.e. do nothing). Available since commit [b6b4c32](https://github.com/gridap/Gridap.jl/commit/b6b4c32c8c4b826a41ba64c770ac8a1c394e16f0)
- Function `restrict` for restricting `CellField` and `CellBasis` objects to surfaces. Available since commit [e981f3c](https://github.com/gridap/Gridap.jl/commit/e981f3c221f3624cfc6764efa47f22652fc22b4f)
- `BoundaryTriangulation` an integration mesh used to integrate `CellField` and `CellBasis` objects restricted on a surface. Available since commit [e981f3c](https://github.com/gridap/Gridap.jl/commit/e981f3c221f3624cfc6764efa47f22652fc22b4f)
- `NonIterableCellMap`, a cell map that has iteration intentionally disabled. Available since commit [956a537](https://github.com/gridap/Gridap.jl/commit/956a5374db6c3b9546e85e0d4d49ae0560057565).
- `CompressedCellValue`, `CompressedCellArray`, and `CompressedCellMap`, as well as efficient versions of `apply`, `evaluate`, and `reindex` for these types. See PR [#41](https://github.com/gridap/Gridap.jl/pull/41) for more details.
- `NEWS.md` file (a changelog file)

### Changed
- Domains are now in [0,1] instead of [-1,1]. Quadratures and nodes arrays modified accordingly. Since commit [268dfe1](https://github.com/gridap/Gridap.jl/commit/268dfe12ef7d736fcd9ad0b9b256740aaf15b2e7).
- Changed the signature of `assemble`, `apply_constraints`, `apply_constraints_rows`, and `apply_constraints_cols` to support FE assembly of several terms, which  are integrated in different domains. The old API of `asseble` is still functional, but not for the `apply_constraints` et al. Since PR [#43](https://github.com/gridap/Gridap.jl/pull/43). Further changed in commit   [a335aed](https://github.com/gridap/Gridap.jl/commit/a335aede65c92a1f61f0ff0dbb0fb44cc20cf906).

### Fixed

- Bug in generation of the cellwise local to global DOF map for high order interpolations. Fixed via PR [#56](https://github.com/gridap/Gridap.jl/pull/56).
- Bug in numerical integration. There was a bug for computations where the number of cell DOFs was different from the number of integration points. Fixed via commit [0b3d4bf](https://github.com/gridap/Gridap.jl/commit/0b3d4bfadea48707c748fca0de65a51a598b6ca6)

## [0.2.0] - 2019-06-29

A changelog is not maintained for this version.

This version introduces the core finite element machinery for linear and non-linear problems,
single field and multi-field problems with terms integrated over the interior of the computational domain.

## [0.1.0] - 2019-05-20

A changelog is not maintained for this version.

This version is non functional. It is just a tag for registering the package.
# Gridap Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project administrators, either 
[Santiago Badia](mailto:santiago.badia@monash.edu) or 
[Francesc Verdugo](mailto:fverdugo@cimne.upc.edu). All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org/), version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
Contributing to `Gridap`
==

By contributing to `Gridap`, you accept and agree to the following *Developer Certificate of Origin Version 1.1* (see below) for all your contributions.

```
Developer Certificate of Origin
Version 1.1

Copyright (C) 2004, 2006 The Linux Foundation and its contributors.
1 Letterman Drive
Suite D4700
San Francisco, CA, 94129

Everyone is permitted to copy and distribute verbatim copies of this
license document, but changing it is not allowed.


Developer's Certificate of Origin 1.1

By making a contribution to this project, I certify that:

(a) The contribution was created in whole or in part by me and I
    have the right to submit it under the open source license
    indicated in the file; or

(b) The contribution is based upon previous work that, to the best
    of my knowledge, is covered under an appropriate open source
    license and I have the right under that license to submit that
    work with modifications, whether created in whole or in part
    by me, under the same open source license (unless I am
    permitted to submit under a different license), as indicated
    in the file; or

(c) The contribution was provided directly to me by some other
    person who certified (a), (b) or (c) and I have not modified
    it.

(d) I understand and agree that this project and the contribution
    are public and that a record of the contribution (including all
    personal information I submit with it, including my sign-off) is
    maintained indefinitely and may be redistributed consistent with
    this project or the open source license(s) involved.
```


Gridap Style Guides
===

* 2 spaces for indentation level
* No trailing white spaces
* CamelCase for typenames
* Pluralized CamelCase for files that implement a type
* CamelCasesTests for CamelCase type test file
* Use lowercase for methods, with underscores only when necessary
* Use whitespace for readability
* 80 characterl line length limit
* Use method! for muting methods
* Wrap multiline expressions in parentheses to avoid errors

See [the Julia CONTRIBUTING.md](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md) for further information.

```@meta
CurrentModule = Gridap.Polynomials
```

# Gridap.Polynomials

```@autodocs
Modules = [Polynomials,]
```
# Gridap.jl

Documentation of the Gridap library.

!!! note

     These documentation pages are under construction.

## Introduction

Gridap provides a set of tools for the grid-based approximation of
partial differential equations (PDEs) written in the
[Julia programming language](https://julialang.org/).
The main motivation behind the development of this library is to provide an easy-to-use framework for the development of complex PDE solvers in a dynamically typed style without sacrificing the performance of statically typed languages.
The library currently supports linear and nonlinear PDE systems for scalar and vector fields, single and multi-field problems, conforming and nonconforming finite element discretizations, on structured and unstructured meshes of simplices and hexahedra.

## How to use this documentation

* The first step for new users is to visit the [Getting Started](@ref) page.

* A set of tutorials written as Jupyter notebooks and html pages are available [here](https://github.com/gridap/Tutorials).

* The detailed documentation is in the [Manual](@ref) section.

* Guidelines for developers of the Gridap project is found in the [Gridap wiki](https://github.com/gridap/Gridap.jl/wiki) page.

## Julia educational resources

A basic knowledge of the Julia programming language is needed to use the Gridap package.
Here, one can find a list of resources to get started with this programming language.

* First steps to learn Julia form the [Gridap wiki](https://github.com/gridap/Gridap.jl/wiki/Start-learning-Julia) page.
* Official webpage [docs.julialang.org](https://docs.julialang.org/)
* Official list of learning resources [julialang.org/learning](https://julialang.org/learning/)

## Manual

```@contents
Pages = [
  "Gridap.md",
  "Helpers.md",
  "Io.md",
  "Algebra.md",
  "Arrays.md",
  "TensorValues.md",
  "Fields.md",
  "Polynomials.md",
  "Integration.md",
  "ReferenceFEs.md",
  "Geometry.md",
  "CellData.md",
  "Visualization.md",
  "FESpaces.md",
  "MultiField.md",
  ]
```


```@meta
CurrentModule = Gridap.Geometry
```

# Gridap.Geometry

```@autodocs
Modules = [Geometry,]
```
```@meta
CurrentModule = Gridap.TensorValues
```

# Gridap.TensorValues

```@autodocs
Modules = [TensorValues,]
```

```@meta
CurrentModule = Gridap.CellData
```

# Gridap.CellData

```@autodocs
Modules = [CellData,]
```

```@meta
CurrentModule = Gridap.Fields
```

# Gridap.Fields

```@autodocs
Modules = [Fields,]
```

```@meta
CurrentModule = Gridap.Visualization
```

# Gridap.Visualization

```@autodocs
Modules = [Visualization,]
```
# Gridap

```@docs
Gridap
```

```@meta
CurrentModule = Gridap.Helpers
```

# Gridap.Helpers

```@autodocs
Modules = [Helpers,]
```


```@meta
CurrentModule = Gridap.Algebra
```

# Gridap.Algebra

```@autodocs
Modules = [Algebra,]
```
# Getting Started

## Installation requirements

Gridap is tested on Linux, but it should be also possible to use it on Mac OS and Windows since it is written exclusively in Julia and it only depends on registered Julia packages.

## Installation

Gridap is a registered package. Thus, the installation should be straight forward using the Julia's package manager [Pkg](https://julialang.github.io/Pkg.jl/v1/). To this end, open the Julia REPL (i.e., execute the `julia` binary), type `]` to enter package mode, and install Gridap as follows

```julia
pkg> add Gridap
```

That's all.

For further information about how to install and manage Julia packages, see the
[Pkg documentation](https://julialang.github.io/Pkg.jl/v1/).

## Further steps

We recommend to follow the [Gridap Tutorials](https://gridap.github.io/Tutorials/dev/) in order to get familiar with the library.


```@meta
CurrentModule = Gridap.ReferenceFEs
```

# Gridap.ReferenceFEs

```@autodocs
Modules = [ReferenceFEs,]
```


```@meta
CurrentModule = Gridap.Io
```

# Gridap.Io

```@autodocs
Modules = [Io,]
```

```@meta
CurrentModule = Gridap.MultiField
```

# Gridap.MultiField

```@autodocs
Modules = [MultiField,]
```

```@meta
CurrentModule = Gridap.Arrays
```
# Gridap.Arrays

```@autodocs
Modules = [Arrays,]
```


```@meta
CurrentModule = Gridap.FESpaces
```

# Gridap.FESpaces

```@autodocs
Modules = [FESpaces,]
```
