# SummationByPartsOperators.jl: A Julia library of provably stable discretization techniques with mimetic properties

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ranocha.github.io/SummationByPartsOperators.jl/stable)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ranocha.github.io/SummationByPartsOperators.jl/dev)
[![Build Status](https://github.com/ranocha/SummationByPartsOperators.jl/workflows/CI/badge.svg)](https://github.com/ranocha/SummationByPartsOperators.jl/actions?query=workflow%3ACI)
[![Codecov](http://codecov.io/github/ranocha/SummationByPartsOperators.jl/coverage.svg?branch=main)](http://codecov.io/github/ranocha/SummationByPartsOperators.jl?branch=main)
[![Coveralls](https://coveralls.io/repos/github/ranocha/SummationByPartsOperators.jl/badge.svg?branch=main)](https://coveralls.io/github/ranocha/SummationByPartsOperators.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![JOSS](https://joss.theoj.org/papers/c1bc6f211c4cce38bfdd0d312816bc69/status.svg)](https://joss.theoj.org/papers/c1bc6f211c4cce38bfdd0d312816bc69)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4773575.svg)](https://doi.org/10.5281/zenodo.4773575)
<!-- [![GitHub commits since tagged version](https://img.shields.io/github/commits-since/ranocha/SummationByPartsOperators.jl/v0.5.5.svg?style=social&logo=github)](https://github.com/ranocha/SummationByPartsOperators.jl) -->
<!-- [![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/SummationByPartsOperators.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html) -->

The Julia library
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
provides a unified interface of different discretization approaches including
finite difference, Fourier pseudospectral, continuous Galerkin, and discontinuous
Galerkin methods.
This unified interface is based on the notion of summation-by-parts (SBP)
operators. Originally developed for finite difference methods, SBP operators
are discrete derivative operators designed specifically to get provably stable
(semi-) discretizations, mimicking energy/entropy estimates from the continuous
level discretely and paying special attention to boundary conditions.

SummationByPartsOperators.jl is mainly written to be useful for both students
learning the basic concepts and researchers developing new numerical algorithms
based on SBP operators. Thus, this package uses Julia's multiple dispatch and
strong type system to provide a unified framework of all of these seemingly
different discretizations while being reasonably optimized at the same time,
achieving good performance without sacrificing flexibility.


## Installation

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is a registered Julia package. Thus, you can install it from the Julia REPL via
```julia
julia> using Pkg; Pkg.add("SummationByPartsOperators")
```

If you want to update SummationByPartsOperators.jl, you can use
```julia
julia> using Pkg; Pkg.update("SummationByPartsOperators")
```
As usual, if you want to update SummationByPartsOperators.jl and all other
packages in your current project, you can execute
```julia
julia> using Pkg; Pkg.update()
```
A brief list of notable changes is available in [`NEWS.md`](NEWS.md).


## Basic examples

Compute the derivative on a periodic domain using a central finite difference operator.
```julia
julia> using SummationByPartsOperators

julia> using Plots: plot, plot!

julia> D = periodic_derivative_operator(derivative_order=1, accuracy_order=2,
                                        xmin=0.0, xmax=2.0, N=20)
Periodic first-derivative operator of order 2 on a grid in [0.0, 2.0] using 20 nodes,
stencils with 1 nodes to the left, 1 nodes to the right, and coefficients of Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.

julia> x = grid(D); u = sinpi.(x);

julia> plot(x, D * u, label="numerical")

julia> plot!(x, π .* cospi.(x), label="analytical")
```
You should see a plot like the following.

<p align="center">
  <img width="300px" src="https://user-images.githubusercontent.com/12693098/118977199-2ef4b280-b976-11eb-8e02-aec722d75bfa.png">
</p>


Compute the derivative on a bounded domain using an SBP finite difference operator.
```julia
julia> using SummationByPartsOperators

julia> using Plots: plot, plot!

julia> D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=2,
                               xmin=0.0, xmax=1.0, N=21)
SBP first-derivative operator of order 2 on a grid in [0.0, 1.0] using 21 nodes
and coefficients of Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
    derivatives.
  Journal of Computational Physics 199, pp. 503-540.

julia> x = grid(D); u = exp.(x);

julia> plot(x, D * u, label="numerical")

julia> plot!(x, exp.(x), label="analytical")
```
You should see a plot like the following.

<p align="center">
  <img width="300px" src="https://user-images.githubusercontent.com/12693098/118978404-93fcd800-b977-11eb-80b3-3dbfce5ecfd6.png">
</p>



## Brief overview

The following derivative operators are implemented as "lazy"/matrix-free
operators, i.e. no large (size of the computational grid) matrix is formed
explicitly. They are linear operators and implement the same interface as
matrices in Julia (at least partially). In particular, `*` and `mul!` are
supported.


### Periodic domains

- `periodic_derivative_operator(; derivative_order, accuracy_order, xmin, xmax, N)`

  These are classical central finite difference operators using `N` nodes on the
  interval `[xmin, xmax]`.

- `periodic_derivative_operator(Holoborodko2008(); derivative_order, accuracy_order, xmin, xmax, N)`

  These are central finite difference operators using `N` nodes on the
  interval `[xmin, xmax]` and the coefficients of
  [Pavel Holoborodko](http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/).

- `fourier_derivative_operator(; xmin, xmax, N)`

  Fourier derivative operators are implemented using the fast Fourier transform of
  [FFTW.jl](https://github.com/JuliaMath/FFTW.jl).

All of these periodic derivative operators support multiplication and addition
such that polynomials and rational functions of them can be represented efficiently,
e.g. to solve elliptic problems of the form `u = (D^2 + I) \ f`.


### Finite (nonperiodic) domains

- `derivative_operator(source_of_coefficients; derivative_order, accuracy_order, xmin, xmax, N)`

  Finite difference SBP operators for first and second derivatives can be obtained
  by using `MattssonNordström2004()` as `source_of_coefficients`.
  Other sources of coefficients are implemented as well. To obtain a full list
  of all operators, use `subtypes(SourceOfCoefficients)`.

- `legendre_derivative_operator(; xmin, xmax, N)`

  Use Lobatto Legendre polynomial collocation schemes on `N`, i.e.
  polynomials of degree `N-1`, implemented via
  [PolynomialBases.jl](https://github.com/ranocha/PolynomialBases.jl).


### Dissipation operators

Additionally, some artificial dissipation/viscosity operators are implemented.
The most basic usage is `Di = dissipation_operator(D)`,
where `D` can be a (periodic, Fourier, Legendre, SBP FD) derivative
operator. Use `?dissipation_operator` for more details.


### Continuous and discontinuous Galerkin methods

SBP operators on bounded domains can be coupled continuously or discontinuously
to obtain CG//DG-type methods. You need to create an appropriate `mesh` and
a basic operator `D` that should be used on each element.
Then, global CG/DG operators are constructed lazily/matrix-free by calling
`couple_continuously(D, mesh)` or
`couple_discontinuously(D, mesh, coupling::Union{Val{:plus}, Val{:central}, Val{:minus}}=Val(:central))`.
Choosing `coupling=Val(:central)` yields a classical SBP operator; the other two
`coupling` types result in upwind SBP operators. Currently, only uniform meshes

- `UniformMesh1D(xmin::Real, xmax::Real, Nx::Integer)`
- `UniformPeriodicMesh1D(xmin::Real, xmax::Real, Nx::Integer)`

are implemented.


### Conversion to other forms

Sometimes, it can be convenient to obtain an explicit (sparse, banded) matrix form
of the operators. Therefore, some conversion functions are supplied, e.g.
```julia
julia> using SummationByPartsOperators

julia> D = derivative_operator(MattssonNordström2004(),
                               derivative_order=1, accuracy_order=2,
                               xmin=0.0, xmax=1.0, N=5)
SBP first-derivative operator of order 2 on a grid in [0.0, 1.0] using 5 nodes
and coefficients of Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
    derivatives.
  Journal of Computational Physics 199, pp. 503-540.

julia> Matrix(D)
5×5 Array{Float64,2}:
 -4.0   4.0   0.0   0.0  0.0
 -2.0   0.0   2.0   0.0  0.0
  0.0  -2.0   0.0   2.0  0.0
  0.0   0.0  -2.0   0.0  2.0
  0.0   0.0   0.0  -4.0  4.0

julia> using SparseArrays

julia> sparse(D)
5×5 SparseMatrixCSC{Float64, Int64} with 10 stored entries:
 -4.0   4.0    ⋅     ⋅    ⋅
 -2.0    ⋅    2.0    ⋅    ⋅
   ⋅   -2.0    ⋅    2.0   ⋅
   ⋅     ⋅   -2.0    ⋅   2.0
   ⋅     ⋅     ⋅   -4.0  4.0

julia> using BandedMatrices

julia> BandedMatrix(D)
5×5 BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 -4.0   4.0    ⋅     ⋅    ⋅
 -2.0   0.0   2.0    ⋅    ⋅
   ⋅   -2.0   0.0   2.0   ⋅
   ⋅     ⋅   -2.0   0.0  2.0
   ⋅     ⋅     ⋅   -4.0  4.0
```


## Documentation

The latest documentation is available
[online](https://ranocha.github.io/SummationByPartsOperators.jl/stable)
and under [`docs/src`](docs/src).
Some additional examples can be found in the directory
[`notebooks`](https://github.com/ranocha/SummationByPartsOperators.jl/tree/main/notebooks).
In particular, examples of complete discretizations of
[the linear advection equation](https://github.com/ranocha/SummationByPartsOperators.jl/blob/main/notebooks/Advection_equation.ipynb),
[the heat equation](https://github.com/ranocha/SummationByPartsOperators.jl/blob/main/notebooks/Heat_equation.ipynb),
and the [wave equation](https://github.com/ranocha/SummationByPartsOperators.jl/blob/main/notebooks/Wave_equation.ipynb) are available.
Further examples are supplied as
[tests](https://github.com/ranocha/SummationByPartsOperators.jl/tree/main/test).


## Referencing

If you use
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
for your research, please cite it using the bibtex entry
```bibtex
@article{ranocha2021sbp,
  title={{SummationByPartsOperators.jl}: {A} {J}ulia library of provably stable
         semidiscretization techniques with mimetic properties},
  author={Ranocha, Hendrik},
  journal={Journal of Open Source Software},
  year={2021},
  month={08},
  doi={10.21105/joss.03454},
  volume={6},
  number={64},
  pages={3454},
  publisher={The Open Journal},
  url={https://github.com/ranocha/SummationByPartsOperators.jl}
}
```


## License and contributing

This project is licensed under the MIT license (see [LICENSE.md](LICENSE.md)).
Since it is an open-source project, we are very happy to accept contributions
from the community. Please refer to [CONTRIBUTING.md](CONTRIBUTING.md) for more
details.
# Changelog

SummationByPartsOperators.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.


## Changes in the v0.5 lifecycle

#### Deprecated

- The (keyword) argument `parallel::Union{Val{:serial}, Val{:threads}}`
  is deprecated in favor of `mode` with possible values
  `FastMode()` (default), `SafeMode()`, and `ThreadedMode()`


## Breaking changes from v0.4.x to v0.5

- Switch from British English to American English consistently, e.g.,
  `semidiscretise` → `semidiscretize`
- `add_transpose_derivative_left!` and `add_transpose_derivative_right!`
  were replaced by the more general functions
  `mul_transpose_derivative_left!` and `mul_transpose_derivative_right!`,
  which use the same interface as `mul!`
- The number of nodes passed to `periodic_central_derivative_operator`, and
  `periodic_derivative_operator` changed from the number of visualization nodes
  to the number of compute nodes (= number of visualization nodes minus one),
  in accordance with `fourier_derivative_operator`
# Contributing

SummationByPartsOperators.jl is an open-source project and we are very happy
to accept contributions from the community. Please feel free to open issues
or submit patches (preferably as pull requests) any time. For planned larger
contributions, it is often beneficial to get in contact first, for example via
issues.

SummationByPartsOperators.jl and its contributions are licensed under the
MIT license (see [LICENSE.md](LICENSE.md)). As a contributor, you certify
that all your contributions are in conformance with the *Developer Certificate
of Origin (Version 1.1)*, which is reproduced below.

## Developer Certificate of Origin (Version 1.1)
The following text was taken from
[https://developercertificate.org](https://developercertificate.org):

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
# Notebooks for SummationByPartsOperators.jl

These [Jupyter](https://jupyter.org) notebooks contain some examples using the Julia package
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl).
# SummationByPartsOperators.jl

The Julia library
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
provides a unified interface of different discretization approaches including
finite difference, Fourier pseudospectral, continuous Galerkin, and discontinuous
Galerkin methods.
This unified interface is based on the notion of summation-by-parts (SBP)
operators. Originally developed for finite difference methods, SBP operators
are discrete derivative operators designed specifically to get provably stable
(semi-) discretizations, mimicking energy/entropy estimates from the continuous
level discretely and paying special attention to boundary conditions.

SummationByPartsOperators.jl is mainly written to be useful for both students
learning the basic concepts and researchers developing new numerical algorithms
based on SBP operators. Thus, this package uses Julia's multiple dispatch and
strong type system to provide a unified framework of all of these seemingly
different discretizations while being reasonably optimized at the same time,
achieving good performance without sacrificing flexibility.


## Installation

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
is a registered Julia package. Thus, you can install it from the Julia REPL via
```julia
julia> using Pkg; Pkg.add("SummationByPartsOperators")
```

If you want to update SummationByPartsOperators.jl, you can use
```julia
julia> using Pkg; Pkg.update("SummationByPartsOperators")
```
As usual, if you want to update SummationByPartsOperators.jl and all other
packages in your current project, you can execute
```julia
julia> using Pkg; Pkg.update()
```


## Basic examples

Compute the derivative on a periodic domain using a central finite difference operator.
```julia
julia> using SummationByPartsOperators

julia> using Plots: plot, plot!

julia> D = periodic_derivative_operator(derivative_order=1, accuracy_order=2,
                                        xmin=0.0, xmax=2.0, N=20)
Periodic first-derivative operator of order 2 on a grid in [0.0, 2.0] using 20 nodes,
stencils with 1 nodes to the left, 1 nodes to the right, and coefficients of Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.

julia> x = grid(D); u = sinpi.(x);

julia> plot(x, D * u, label="numerical")

julia> plot!(x, π .* cospi.(x), label="analytical")
```
You should see a plot like the following.

![](https://user-images.githubusercontent.com/12693098/118977199-2ef4b280-b976-11eb-8e02-aec722d75bfa.png)


Compute the derivative on a bounded domain using an SBP finite difference operator.
```julia
julia> using SummationByPartsOperators

julia> using Plots: plot, plot!

julia> D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=2,
                               xmin=0.0, xmax=1.0, N=21)
SBP first-derivative operator of order 2 on a grid in [0.0, 1.0] using 21 nodes
and coefficients of Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
    derivatives.
  Journal of Computational Physics 199, pp. 503-540.

julia> x = grid(D); u = exp.(x);

julia> plot(x, D * u, label="numerical")

julia> plot!(x, exp.(x), label="analytical")
```
You should see a plot like the following.

![](https://user-images.githubusercontent.com/12693098/118978404-93fcd800-b977-11eb-80b3-3dbfce5ecfd6.png)


## Referencing

If you use
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
for your research, please cite it using the bibtex entry
```bibtex
@article{ranocha2021sbp,
  title={{SummationByPartsOperators.jl}: {A} {J}ulia library of provably stable
         semidiscretization techniques with mimetic properties},
  author={Ranocha, Hendrik},
  journal={Journal of Open Source Software},
  year={2021},
  month={08},
  doi={10.21105/joss.03454},
  volume={6},
  number={64},
  pages={3454},
  publisher={The Open Journal},
  url={https://github.com/ranocha/SummationByPartsOperators.jl}
}
```
Please also cite the appropriate references for specific SBP operators
you use, which can be obtained via [`source_of_coefficients`](@ref).


## License and contributing

This project is licensed under the MIT license (see [License](@ref)).
Since it is an open-source project, we are very happy to accept contributions
from the community. Please refer to the section [Contributing](@ref) for more
details.
# Introduction

Summation-by-parts (SBP) operators are discrete derivative operators designed to
enable (semi-) discrete stability proofs mimicking the energy method from the
continuous level. To do so, SBP operators mimic integration-by-parts discretely.
Here, we will briefly explain the basic concepts. If you want to learn more about
this subject, the classical review articles of [^SvärdNordström2014] and
[^FernándezHickenZingg2014] are good starting points. More recent references and
applications of SBP operators from many classes implemented in
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
are given by [^RanochaMitsotakisKetcheson2021].

Since SBP operators are designed to mimic integration-by-parts, they need a notion
of derivatives and integrals. Here, derivatives are interpreted as linear operators `D`
(derivative matrices) and integrals are interpreted as discrete inner products,
represented by the associated mass/norm matrices `M`. Thus, the discrete derivative
of a grid function `u` is `D * u` and the discrete inner product of two grid functions
`u` and `v` is `dot(u, M, v)`, where `M = mass_matrix(D)`. Here, we have already
introduced some basic interfaces provided by
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl):
- Derivative operators act as linear operators implementing `*` (and `mul!` for
  more efficient in-place updates avoiding allocations).
- The mass matrix associated to an SBP derivative operator can be retrieved via
  [`mass_matrix`](@ref).


## Periodic domains

Periodic (central) SBP operators mimic the properties of differential operators
on periodic domains. Hence, they are

- skew-symmetric if they approximate odd derivatives
- symmetric and semi-definite if they approximate even derivatives;
  second-derivative operators are negative semi-definite,
  fourth-derivative operators are positive semi-definite etc.

Classical central finite difference operators on periodic domains are periodic
SBP operators. They can be constructed via [`periodic_derivative_operator`](@ref).
Similarly, Fourier collocation methods can be interpreted as periodic SBP operators,
which can be constructed via [`fourier_derivative_operator`](@ref).

```jldoctest
julia> using SummationByPartsOperators, LinearAlgebra

julia> D = periodic_derivative_operator(derivative_order=1, accuracy_order=2,
                                        xmin=0.0, xmax=2.0, N=20)
Periodic first-derivative operator of order 2 on a grid in [0.0, 2.0] using 20 nodes,
stencils with 1 nodes to the left, 1 nodes to the right, and coefficients of Fornberg (1998)
  Calculation of Weights in Finite Difference Formulas.
  SIAM Rev. 40.3, pp. 685-691.

julia> M = mass_matrix(D)
UniformScaling{Float64}
0.1*I

julia> M * Matrix(D) + Matrix(D)' * M |> norm
0.0

julia> D = fourier_derivative_operator(xmin=0.0, xmax=2.0, N=20)
Periodic 1st derivative Fourier operator {T=Float64}
on a grid in [0.0, 2.0] using 20 nodes and 11 modes

julia> M = mass_matrix(D)
UniformScaling{Float64}
0.1*I

julia> norm(M * Matrix(D) + Matrix(D)' * M) < 10 * eps(eltype(D))
true
```

As you have seen above, conversion methods to other common types such as `Matrix`,
`sparse` from the standard library SparseArrays, and `BandedMatrix` from
[BandedMatrices.jl](https://github.com/JuliaMatrices/BandedMatrices.jl) are
available.


## Non-periodic domains

On non-periodic domains, additional boundary terms appear. Thus, the basic
symmetry properties of SBP operators are the same as the ones of periodic SBP
operators modulo boundary terms. Note that the correct handling of boundary terms
is the basic reason of the success of SBP operators. In particular for hyperbolic
problems, other boundary treatments that might appear senseful can result in
catastrophic failure.

### First-derivative operators

First-derivative SBP operators need to mimic
```math
  \int_{x_\mathrm{min}}^{x_\mathrm{max}} u(x) \bigl( \partial_x v(x) \bigr) \mathrm{d}x
+ \int_{x_\mathrm{min}}^{x_\mathrm{max}} \bigl( \partial_x u(x) \bigr) v(x) \mathrm{d}x
= u(x_\mathrm{max}) v(x_\mathrm{max}) - u(x_\mathrm{min}) v(x_\mathrm{min}).
```
Thus, a discrete evaluation at the boundary of the domain is necessary. For
SBP operators with a grid including the boundary nodes, this can be achieved
by simply picking the first/last nodal coefficient of a grid function `u`.
If boundary nodes are not included, some interpolation is necessary in general.
Nevertheless, getting a boundary value is a linear functional that is often
represented in the literature using (transposed) vectors `tL, tR`. Then,
an SBP operator has to satisfy `M * D + D' * M == tR * tR' - tL * tL'`.
The boundary operators are represented matrix-free via
[`derivative_left`](@ref) and [`derivative_right`](@ref) for zeroth-order
derivatives.

```jldoctest; filter = r"((┌.*[\n\r]+)(│ `LoopVectorization.check_args` on your inputs failed; running fallback `@inbounds @fastmath` loop instead\.*[\n\r]+)(│.*[\n\r]+)(└.*[\n\r]+))"
julia> using SummationByPartsOperators, LinearAlgebra

julia> D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=2,
                               xmin=0//1, xmax=1//1, N=9)
SBP first-derivative operator of order 2 on a grid in [0//1, 1//1] using 9 nodes
and coefficients of Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
    derivatives.
  Journal of Computational Physics 199, pp. 503-540.

julia> tL = zeros(eltype(D), size(D, 1)); tL[1] = 1; tL'
1×9 adjoint(::Vector{Rational{Int64}}) with eltype Rational{Int64}:
 1//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1

julia> tR = zeros(eltype(D), size(D, 1)); tR[end] = 1; tR'
1×9 adjoint(::Vector{Rational{Int64}}) with eltype Rational{Int64}:
 0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  1//1

julia> M = mass_matrix(D)
9×9 Diagonal{Rational{Int64}, Vector{Rational{Int64}}}:
 1//16   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅     1//8   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅    1//8   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅    1//8   ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅    1//8   ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅    1//8   ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅     ⋅    1//8   ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    1//8   ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    1//16

julia> M * Matrix(D) + Matrix(D)' * M == tR * tR' - tL * tL'
true

julia> u = randn(size(grid(D))); derivative_left(D, u, Val(0)) == u[begin]
true

julia> u = randn(size(grid(D))); derivative_right(D, u, Val(0)) == u[end]
true
```

Here, we have introduced some additional features. Firstly, exact rational
coefficients are provided, based on the type of `xmin` and `xmax` (if available).
Secondly, a [`source_of_coefficients`](@ref) has to be provided when constructing
the SBP operator. You can list them using
```@example
using InteractiveUtils, SummationByPartsOperators
subtypes(SourceOfCoefficients)
```
Here and in the following, the order of accuracy of (finite difference) SBP
operators refers to the local order of accuracy in the interior, cf.
[`accuracy_order`](@ref).

A special case of first-derivative SBP operators are polynomial derivative operators
on Lobatto-Legendre nodes, implemented in [`legendre_derivative_operator`](@ref).

### Second-derivative operators

To mimic integration-by-parts of second derivatives,
```math
  \int_{x_\mathrm{min}}^{x_\mathrm{max}} u(x) \bigl( \partial_x^2 v(x) \bigr) \mathrm{d}x
= - \int_{x_\mathrm{min}}^{x_\mathrm{max}} \bigl( \partial_x u(x) \bigr) \bigl( \partial_x v(x) \bigr) \mathrm{d}x
  + u(x_\mathrm{max}) \bigl( \partial_x v(x_\mathrm{max}) \bigr)
  - \bigl( \partial_x u(x_\mathrm{min})) v(x_\mathrm{min}),
```
the evaluation of the first derivative at the boundaries is necessary. These
linear functionals are available as [`derivative_left`](@ref) and
[`derivative_right`](@ref). In the literature, they are often called `dL` and
`dR`. Then, a second-derivative SBP operator has to be of the form
`M * D == -A + tR * dR' - tL * dL'`, where `A` is symmetric and positive
semidefinite.

```jldoctest
julia> using SummationByPartsOperators, LinearAlgebra

julia> D = derivative_operator(MattssonNordström2004(), derivative_order=2, accuracy_order=2,
                               xmin=0//1, xmax=1//1, N=9)
SBP second-derivative operator of order 2 on a grid in [0//1, 1//1] using 9 nodes
and coefficients of Mattsson, Nordström (2004)
  Summation by parts operators for finite difference approximations of second
    derivatives.
  Journal of Computational Physics 199, pp. 503-540.

julia> M = mass_matrix(D)
9×9 Diagonal{Rational{Int64}, Vector{Rational{Int64}}}:
 1//16   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅     1//8   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅    1//8   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅    1//8   ⋅     ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅    1//8   ⋅     ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅    1//8   ⋅     ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅     ⋅    1//8   ⋅     ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    1//8   ⋅
  ⋅      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    1//16

julia> tL = derivative_left(D, Val(0)); tL'
1×9 adjoint(::Vector{Bool}) with eltype Bool:
 1  0  0  0  0  0  0  0  0

julia> tR = derivative_right(D, Val(0)); tR'
1×9 adjoint(::Vector{Bool}) with eltype Bool:
 0  0  0  0  0  0  0  0  1

julia> dL = derivative_left(D, Val(1)); dL'
1×9 adjoint(::Vector{Rational{Int64}}) with eltype Rational{Int64}:
 -12//1  16//1  -4//1  0//1  0//1  0//1  0//1  0//1  0//1

julia> dR = derivative_right(D, Val(1)); dR'
1×9 adjoint(::Vector{Rational{Int64}}) with eltype Rational{Int64}:
 0//1  0//1  0//1  0//1  0//1  0//1  4//1  -16//1  12//1

julia> A = -M * Matrix(D) + tR * dR' - tL * dL'
9×9 Matrix{Rational{Int64}}:
  8//1  -8//1   0//1   0//1   0//1   0//1   0//1   0//1   0//1
 -8//1  16//1  -8//1   0//1   0//1   0//1   0//1   0//1   0//1
  0//1  -8//1  16//1  -8//1   0//1   0//1   0//1   0//1   0//1
  0//1   0//1  -8//1  16//1  -8//1   0//1   0//1   0//1   0//1
  0//1   0//1   0//1  -8//1  16//1  -8//1   0//1   0//1   0//1
  0//1   0//1   0//1   0//1  -8//1  16//1  -8//1   0//1   0//1
  0//1   0//1   0//1   0//1   0//1  -8//1  16//1  -8//1   0//1
  0//1   0//1   0//1   0//1   0//1   0//1  -8//1  16//1  -8//1
  0//1   0//1   0//1   0//1   0//1   0//1   0//1  -8//1   8//1

julia> isposdef(A)
true
```
Usually, there is no need to form `dL, dR` explicitly. Instead, you can use the
matrix-free variants [`derivative_left`](@ref) and [`derivative_right`](@ref).
Some procedures imposing boundary conditions weakly require adding the transposed
boundary derivatives to a grid function, which can be achieved by
[`mul_transpose_derivative_left!`](@ref) and [`mul_transpose_derivative_right!`](@ref).
You can find applications of these operators in the source code of
[`WaveEquationNonperiodicSemidiscretization`](@ref).

A special case of second-derivative SBP operators are polynomial derivative operators
on Lobatto-Legendre nodes, implemented in [`legendre_second_derivative_operator`](@ref).


## Upwind operators

Upwind SBP operators were introduced by [`Mattsson2017`](@ref). They combine
two derivative operators `Dp` (`:plus`) and `Dm` (`:minus`) such that
`M * Dp + Dm' * M == tR * tR' - tL * tL'` and `M * (Dp - Dm)` is negative
semidefinite.

```jldoctest
julia> using SummationByPartsOperators, LinearAlgebra

julia> Dp = derivative_operator(Mattsson2017(:plus), derivative_order=1, accuracy_order=2,
                                xmin=0//1, xmax=1//1, N=9)
SBP first-derivative operator of order 2 on a grid in [0//1, 1//1] using 9 nodes
and coefficients of Mattsson (2017)
  Diagonal-norm upwind SBP operators.
  Journal of Computational Physics 335, pp. 283-310.
  (upwind coefficients plus)

julia> Dm = derivative_operator(Mattsson2017(:minus), derivative_order=1, accuracy_order=2,
                                xmin=0//1, xmax=1//1, N=9)
SBP first-derivative operator of order 2 on a grid in [0//1, 1//1] using 9 nodes
and coefficients of Mattsson (2017)
  Diagonal-norm upwind SBP operators.
  Journal of Computational Physics 335, pp. 283-310.
  (upwind coefficients minus)

julia> M = mass_matrix(Dp)
9×9 Diagonal{Rational{Int64}, Vector{Rational{Int64}}}:
 1//32   ⋅      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅
  ⋅     5//32   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅
  ⋅      ⋅     1//8   ⋅     ⋅     ⋅     ⋅     ⋅      ⋅
  ⋅      ⋅      ⋅    1//8   ⋅     ⋅     ⋅     ⋅      ⋅
  ⋅      ⋅      ⋅     ⋅    1//8   ⋅     ⋅     ⋅      ⋅
  ⋅      ⋅      ⋅     ⋅     ⋅    1//8   ⋅     ⋅      ⋅
  ⋅      ⋅      ⋅     ⋅     ⋅     ⋅    1//8   ⋅      ⋅
  ⋅      ⋅      ⋅     ⋅     ⋅     ⋅     ⋅    5//32   ⋅
  ⋅      ⋅      ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     1//32

julia> M * Matrix(Dp) + Matrix(Dm)' * M
9×9 Matrix{Rational{Int64}}:
 -1//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1
  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1
  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1
  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1
  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1
  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1
  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1
  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1
  0//1  0//1  0//1  0//1  0//1  0//1  0//1  0//1  1//1

julia> minimum(eigvals(-M * (Matrix(Dp) - Matrix(Dm)))) > -100 * eps() # tolerance for zero eigenvalues
true
```


## Continuous and discontinuous Galerkin methods

SBP operators can be coupled to obtain (nodal) continuous Galerkin (CG) methods.
If the underlying SBP operators are [`LegendreDerivativeOperator`](@ref)s,
these are CG spectral element methods (CGSEM). However, a continuous coupling
of arbitrary SBP operators is supported.

```jldoctest
julia> using SummationByPartsOperators, LinearAlgebra

julia> D = couple_continuously(
               legendre_derivative_operator(xmin=-1.0, xmax=1.0, N=3),
               UniformMesh1D(xmin=0.0, xmax=1.0, Nx=3))
First derivative operator {T=Float64} on 3 Lobatto Legendre nodes in [-1.0, 1.0]
coupled continuously on UniformMesh1D{Float64} with 3 cells in (0.0, 1.0)

julia> Matrix(D)
7×7 Matrix{Float64}:
 -9.0  12.0  -3.0   0.0   0.0    0.0   0.0
 -3.0   0.0   3.0   0.0   0.0    0.0   0.0
  1.5  -6.0   0.0   6.0  -1.5    0.0   0.0
  0.0   0.0  -3.0   0.0   3.0    0.0   0.0
  0.0   0.0   1.5  -6.0   0.0    6.0  -1.5
  0.0   0.0   0.0   0.0  -3.0    0.0   3.0
  0.0   0.0   0.0   0.0   3.0  -12.0   9.0

julia> mass_matrix(D)
7×7 Diagonal{Float64, Vector{Float64}}:
 0.0555556   ⋅         ⋅         ⋅         ⋅         ⋅         ⋅
  ⋅         0.222222   ⋅         ⋅         ⋅         ⋅         ⋅
  ⋅          ⋅        0.111111   ⋅         ⋅         ⋅         ⋅
  ⋅          ⋅         ⋅        0.222222   ⋅         ⋅         ⋅
  ⋅          ⋅         ⋅         ⋅        0.111111   ⋅         ⋅
  ⋅          ⋅         ⋅         ⋅         ⋅        0.222222   ⋅
  ⋅          ⋅         ⋅         ⋅         ⋅         ⋅        0.0555556
```

SBP operators can also be coupled as in discontinuous Galerkin (DG) methods.
Using a central numerical flux results in central SBP operators; upwind fluxes
yield upwind SBP operators. If [`LegendreDerivativeOperator`](@ref)s are used,
the discontinuous coupling yields DG spectral element methods (DGSEM).

```jldoctest
julia> using SummationByPartsOperators, LinearAlgebra

julia> D = couple_discontinuously(
               legendre_derivative_operator(xmin=-1.0, xmax=1.0, N=3),
               UniformPeriodicMesh1D(xmin=0.0, xmax=1.0, Nx=3),
               Val(:central))
First derivative operator {T=Float64} on 3 Lobatto Legendre nodes in [-1.0, 1.0]
coupled discontinuously (upwind: Val{:central}()) on UniformPeriodicMesh1D{Float64} with 3 cells in (0.0, 1.0)

julia> M = mass_matrix(D);

julia> M * Matrix(D) + Matrix(D)' * M |> iszero
true
```

Right now, only uniform meshes [`UniformMesh1D`](@ref) and [`UniformPeriodicMesh1D`](@ref)
are implemented.


## Basic interfaces and additional features

To actually compute and plot the discrete grid functions, a few additional ingredients
are necessary.
- The discrete coefficients of a function on the [`grid`](@ref) of an SBP
  operator can usually be computed as `x = grid(D); u = u_function.(x)`,
  at least for nodal bases. In general, [`compute_coefficients`](@ref)
  (or the in-place version [`compute_coefficients!`](@ref))
  can also be used for this task.
- To get a grid and discrete values suitable for plotting, you can use
  [`evaluate_coefficients`](@ref) (or the in-place version
  [`evaluate_coefficients!`](@ref)).
  The plot nodes returned from [`evaluate_coefficients`](@ref) can be different
  from the nodes of the [`grid`](@ref) associated to an SBP operator.
- To implement boundary procedures, the weights of the mass matrix at the boundary
  are often needed. These can be obtained without forming `M = mass_matrix(D)`
  explicitly via [`left_boundary_weight`](@ref) and [`right_boundary_weight`](@ref).
- Instead of forming a mass matrix explicitly, discrete integrals can be evaluated
  efficiently using [`integrate`](@ref).
- Dissipation operators based on the same discrete inner product as SBP derivative
  operators can be obtained via [`dissipation_operator`](@ref).


## Next steps

If you are familiar with SBP operators in general, this introduction might already
be enough for you to apply
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
to your problems. Otherwise, you might want to have a look at the references,
the tutorials coming next,
or some ready-to-use semidiscretizations of the following partial differential
equations (PDEs). These are shipped with this package and you are encouraged to
look at their source code to learn more about it.

- Linear scalar advection with variable coefficient:
  [`VariableLinearAdvectionNonperiodicSemidiscretization`](@ref)
- Burgers' equation (inviscid):
  [`BurgersPeriodicSemidiscretization`](@ref), [`BurgersNonperiodicSemidiscretization`](@ref)
- Scalar conservation law with cubic flux:
  [`CubicPeriodicSemidiscretization`](@ref), [`CubicNonperiodicSemidiscretization`](@ref)
- A scalar conservation law with quartic, non-convex flux:
  [`QuarticNonconvexPeriodicSemidiscretization`](@ref)
- The second-order wave equation:
  [`WaveEquationNonperiodicSemidiscretization`](@ref)

Some additional examples are included as [Jupyter](https://jupyter.org) notebooks
in the directory [`notebooks`](https://github.com/ranocha/SummationByPartsOperators.jl/tree/main/notebooks).
Even more examples and research articles making use of
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
are listed in the section [Applications](@ref).
If you want to know even more, you can have a look at the
[test](https://github.com/ranocha/SummationByPartsOperators.jl/tree/main/test).


## References

[^SvärdNordström2014]:
    Svärd, Nordström (2014).
    Review of summation-by-parts schemes for initial–boundary-value problems.
    [DOI: 10.1016/j.jcp.2014.02.031](https://doi.org/10.1016/j.jcp.2014.02.031)

[^FernándezHickenZingg2014]:
    Fernández, Hicken, Zingg (2014).
    Review of summation-by-parts operators with simultaneous approximation terms
    for the numerical solution of partial differential equations.
    [DOI: 10.1016/j.compfluid.2014.02.016](https://doi.org/10.1016/j.compfluid.2014.02.016)

[^RanochaMitsotakisKetcheson2021]:
    Ranocha, Mitsotakis, Ketcheson (2021).
    A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations.
    [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)

# Applications

Here is a (non-exhaustive) list of research using
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl).

- Hendrik Ranocha, Manuel Quezada de Luna, and David I Ketcheson (2021).
  On the Rate of Error Growth in Time for Numerical Solutions of
  Nonlinear Dispersive Wave Equations.
  [arXiv: 2102.07376 [math.NA]](https://arxiv.org/abs/2102.07376)
  A reproducibility repository containing source code for all numerical
  experiments is available at
  [DOI: 10.5281/zenodo.4540467](https://doi.org/10.5281/zenodo.4540467)
- Hendrik Ranocha, Dimitrios Mitsotakis, and David I Ketcheson (2021).
  A Broad Class of Conservative Numerical Methods for Dispersive
  Wave Equations.
  [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
  A reproducibility repository containing source code for all numerical
  experiments is available at
  [DOI: 10.5281/zenodo.3908803](https://doi.org/10.5281/zenodo.3908803)
- Philippe G LeFloch and Hendrik Ranocha (2021).
  Kinetic functions for nonclassical shocks, entropy stability, and
  discrete summation by parts.
  [DOI: 10.1007/s10915-021-01463-6](https://doi.org/10.1007/s10915-021-01463-6)
- Jan Nordström and Hendrik Ranocha (2021).
  A New Class of *A* Stable Summation by Parts Time Integration Schemes
  with Strong Initial Conditions.
  [DOI: 10.1007/s10915-021-01454-7](https://doi.org/10.1007/s10915-021-01454-7)
  A reproducibility repository containing source code for all numerical
  experiments is available at
  [DOI: 10.5281/zenodo.3699173](https://doi.org/10.5281/zenodo.3699173)
- Hendrik Ranocha (2021).
  On Strong Stability of Explicit Runge-Kutta Methods for
  Nonlinear Semibounded Operators.
  [DOI: 10.1093/imanum/drz070](https://doi.org/10.1093/imanum/drz070)
- Hendrik Ranocha and David I Ketcheson (2020).
  Energy Stability of Explicit Runge-Kutta Methods for Nonautonomous
  or Nonlinear Problems.
  [DOI: 10.1137/19M1290346](https://doi.org/10.1137/19M1290346)
  A reproducibility repository containing source code for all numerical
  experiments is available at
  [DOI: 10.5281/zenodo.3464243](https://doi.org/10.5281/zenodo.3464243)
- Hendrik Ranocha, Katharina Ostaszewski, and Philip Heinisch (2020).
  Discrete Vector Calculus and Helmholtz Hodge Decomposition
  for Classical Finite Difference Summation by Parts Operators.
  [DOI: 10.1007/s42967-019-00057-2](https://doi.org/10.1007/s42967-019-00057-2)
  A reproducibility repository containing source code for all numerical
  experiments is available at
  [DOI: 10.5281/zenodo.3375170](https://doi.org/10.5281/zenodo.3375170)
- Hendrik Ranocha and Gregor J Gassner (2020).
  Preventing pressure oscillations does not fix local linear stability
  issues of entropy-based split-form high-order schemes.
  [arXiv: 2009.13139 [math.NA]](https://arxiv.org/abs/2009.13139)
  A reproducibility repository containing source code for all numerical
  experiments is available at
  [DOI: 10.5281/zenodo.4054366](https://doi.org/10.5281/zenodo.4054366)
- Philipp Öffner and Hendrik Ranocha (2019).
  Error Boundedness of Discontinuous Galerkin Methods with Variable
  Coefficients.
  [DOI: 10.1007/s10915-018-00902-1](https://doi.org/10.1007/s10915-018-00902-1)

If you use this package for your own research, please cite it as described
[in the documentation](https://ranocha.de/SummationByPartsOperators.jl/stable/#Referencing)
and make a PR to add your work to the list above.
# SummationByPartsOperators.jl API

```@meta
CurrentModule = SummationByPartsOperators
```

```@autodocs
Modules = [SummationByPartsOperators]
```
# Benchmarks

Here are some simple benchmarks. Take them with a grain of salt since they run
on virtual machines in the cloud to generate the documentation automatically.


## First-derivative operators

#### Periodic domains

Let's set up some benchmark code.

```@example first-derivative-periodic
using BenchmarkTools
using LinearAlgebra, SparseArrays
using SummationByPartsOperators, DiffEqOperators

BLAS.set_num_threads(1) # make sure that BLAS is serial to be fair

T = Float64
xmin, xmax = T(0), T(1)

D_SBP = periodic_derivative_operator(derivative_order=1, accuracy_order=2,
                                     xmin=xmin, xmax=xmax, N=100)
x = grid(D_SBP)
D_DEO = CenteredDifference(derivative_order(D_SBP), accuracy_order(D_SBP),
                           step(x), length(x)) * PeriodicBC(eltype(D_SBP))

D_sparse = sparse(D_SBP)

u = randn(eltype(D_SBP), length(x)); du = similar(u);
@show D_SBP * u ≈ D_DEO * u ≈ D_sparse * u

function doit(D, text, du, u)
  println(text)
  sleep(0.1)
  show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D, $u))
  println()
end
```

First, we benchmark the implementation from SummationByPartsOperators.jl.
```@example first-derivative-periodic
doit(D_SBP, "D_SBP:", du, u)
```

Next, we compare this to the runtime obtained using a sparse matrix representation
of the derivative operator. Depending on the hardware etc., this can be an order
of magnitude slower than the optimized implementation from SummationByPartsOperators.jl.
```@example first-derivative-periodic
doit(D_sparse, "D_sparse:", du, u)
```

Finally, we benchmark the implementation of the same derivative operator in
DiffEqOperators.jl.
```@example first-derivative-periodic
doit(D_DEO, "D_DEO:", du, u)
```


#### Bounded domains

We start again by setting up some benchmark code.

```@example first-derivative-bounded
using BenchmarkTools
using LinearAlgebra, SparseArrays
using SummationByPartsOperators, BandedMatrices

BLAS.set_num_threads(1) # make sure that BLAS is serial to be fair

T = Float64
xmin, xmax = T(0), T(1)

D_SBP = derivative_operator(MattssonNordström2004(), derivative_order=1,
                            accuracy_order=6, xmin=xmin, xmax=xmax, N=10^3)
D_sparse = sparse(D_SBP)
D_banded = BandedMatrix(D_SBP)

u = randn(eltype(D_SBP), size(D_SBP, 1)); du = similar(u);
@show D_SBP * u ≈ D_sparse * u ≈ D_banded * u

function doit(D, text, du, u)
  println(text)
  sleep(0.1)
  show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D, $u))
  println()
end
```

First, we benchmark the implementation from SummationByPartsOperators.jl.
```@example first-derivative-bounded
doit(D_SBP, "D_SBP:", du, u)
```

Again, we compare this to a representation of the derivative operator as a
sparse matrix. No surprise - it is again much slower, as in periodic domains.
```@example first-derivative-bounded
doit(D_sparse, "D_sparse:", du, u)
```

FInally, we compare it to a representation as banded matrix. Disappointingly,
this is still much slower than the optimized implementation from
SummationByPartsOperators.jl.
```@example first-derivative-bounded
doit(D_banded, "D_banded:", du, u)
```


## Dissipation operators

We follow the same structure as before. At first, we set up some benchmark code.

```@example dissipation
using BenchmarkTools
using LinearAlgebra, SparseArrays
using SummationByPartsOperators, BandedMatrices

BLAS.set_num_threads(1) # make sure that BLAS is serial to be fair

T = Float64
xmin, xmax = T(0), T(1)

D_SBP = derivative_operator(MattssonNordström2004(), derivative_order=1,
                            accuracy_order=6, xmin=xmin, xmax=xmax, N=10^3)
Di_SBP  = dissipation_operator(MattssonSvärdNordström2004(), D_SBP)
Di_sparse = sparse(Di_SBP)
Di_banded = BandedMatrix(Di_SBP)
Di_full   = Matrix(Di_SBP)

u = randn(eltype(D_SBP), size(D_SBP, 1)); du = similar(u);
@show Di_SBP * u ≈ Di_sparse * u ≈ Di_banded * u ≈ Di_full * u

function doit(D, text, du, u)
  println(text)
  sleep(0.1)
  show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D, $u))
  println()
end
```

At first, let us benchmark the derivative and dissipation operators implemented
in SummationByPartsOperators.jl.
```@example dissipation
doit(D_SBP, "D_SBP:", du, u)
doit(Di_SBP, "Di_SBP:", du, u)
```

Next, we compare the results to sparse matrix representations. It will not
come as a surprise that these are again much (around an order of magnitude)
slower.
```@example dissipation
doit(Di_sparse, "Di_sparse:", du, u)
doit(Di_banded, "Di_banded:", du, u)
```

Finally, let's benchmark the same computation if a full (dense) matrix is used
to represent the derivative operator. This is obviously a bad idea but 🤷
```@example dissipation
doit(Di_full, "Di_full:", du, u)
```


## Structure-of-Arrays (SoA) and Array-of-Structures (AoS)

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
tries to provide efficient support of

- `StaticVector`s from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)
- [StructArrays.jl](https://github.com/JuliaArrays/StructArrays.jl)

To demonstrate this, let us set up some benchmark code.

```@example soa-aos
using BenchmarkTools
using StaticArrays, StructArrays
using LinearAlgebra, SparseArrays
using SummationByPartsOperators, BandedMatrices

BLAS.set_num_threads(1) # make sure that BLAS is serial to be fair

struct Vec5{T} <: FieldVector{5,T}
  x1::T
  x2::T
  x3::T
  x4::T
  x5::T
end

# Apply `mul!` to each component of a plain array of structures one after another
function mul_aos!(du, D, u, args...)
  for i in 1:size(du, 1)
    mul!(view(du, i, :), D, view(u, i, :), args...)
  end
end

T = Float64
xmin, xmax = T(0), T(1)

D_SBP = derivative_operator(MattssonNordström2004(), derivative_order=1,
                            accuracy_order=4, xmin=xmin, xmax=xmax, N=101)
D_sparse = sparse(D_SBP)
D_full   = Matrix(D_SBP)
```

At first, we benchmark the application of the operators implemented in
SummationByPartsOperators.jl and their representations as sparse and dense
matrices in the scalar case. As before, the sparse matrix representation
is around an order of magnitude slower and the dense matrix representation
is far off.
```@example soa-aos
println("Scalar case")
u = randn(T, size(D_SBP, 1)); du = similar(u)
println("D_SBP")
show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D_SBP, $u))
println("\nD_sparse")
show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D_sparse, $u))
println("\nD_full")
show(stdout, MIME"text/plain"(), @benchmark mul!($du, $D_full, $u))
```

Next, we use a plain array of structures (AoS) in the form of a two-dimensional
array and our custom `mul_aos!` implementation that loops over each component,
using `mul!` on `view`s.
Here, the differences between the timings are less pronounced.
```@example soa-aos
println("Plain Array of Structures")
u_aos_plain = randn(T, 5, size(D_SBP, 1)); du_aos_plain = similar(u_aos_plain)
println("D_SBP")
show(stdout, MIME"text/plain"(), @benchmark mul_aos!($du_aos_plain, $D_SBP, $u_aos_plain))
println("\nD_sparse")
show(stdout, MIME"text/plain"(), @benchmark mul_aos!($du_aos_plain, $D_sparse, $u_aos_plain))
println("\nD_full")
show(stdout, MIME"text/plain"(), @benchmark mul_aos!($du_aos_plain, $D_full, $u_aos_plain))
```

Now, we use an array of structures (AoS) based on `reinterpret` and standard
`mul!`. This is much more efficient for the implementation in SummationByPartsOperators.jl.
In Julia v1.6, this is also more efficient for sparse matrices but less efficient
for dense matrices (compared to the plain AoS approach with `mul_aos!` above).
```@example soa-aos
println("Array of Structures (reinterpreted array)")
u_aos_r = reinterpret(reshape, Vec5{T}, u_aos_plain); du_aos_r = similar(u_aos_r)
@show D_SBP * u_aos_r ≈ D_sparse * u_aos_r ≈ D_full * u_aos_r
mul!(du_aos_r, D_SBP, u_aos_r)
@show reinterpret(reshape, T, du_aos_r) ≈ du_aos_plain
println("D_SBP")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos_r, $D_SBP, $u_aos_r))
println("\nD_sparse")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos_r, $D_sparse, $u_aos_r))
println("\nD_full")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos_r, $D_full, $u_aos_r))
```

Next, we still use an array of structures (AoS), but copy the data into a plain
`Array` instead of using the `reinterpret`ed versions. There is no significant
difference to the previous version in this case.
```@example soa-aos
println("Array of Structures")
u_aos = Array(u_aos_r); du_aos = similar(u_aos)
@show D_SBP * u_aos ≈ D_sparse * u_aos ≈ D_full * u_aos
mul!(du_aos, D_SBP, u_aos)
@show du_aos ≈ du_aos_r
println("D_SBP")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos, $D_SBP, $u_aos))
println("\nD_sparse")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos, $D_sparse, $u_aos))
println("\nD_full")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_aos, $D_full, $u_aos))
```

Finally, let's look at a structure of arrays (SoA). Interestingly, this is
slower than the array of structures we used above. On Julia v1.6, the sparse
matrix representation performs particularly bad in this case.
```@example soa-aos
println("Structure of Arrays")
u_soa = StructArray(u_aos); du_soa = similar(u_soa)
@show D_SBP * u_soa ≈ D_sparse * u_soa ≈ D_full * u_soa
mul!(du_soa, D_SBP, u_soa)
@show du_soa ≈ du_aos
println("D_SBP")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_soa, $D_SBP, $u_soa))
println("\nD_sparse")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_soa, $D_sparse, $u_soa))
println("\nD_full")
show(stdout, MIME"text/plain"(), @benchmark mul!($du_soa, $D_full, $u_soa))
```
# Wave equation

Consider the linear wave equation

```math
\begin{aligned}
    \partial_t^2 u(t,x) &= \partial_x^2 u(t,x), && t \in (0,T), x \in (x_{min}, x_{max}), \\
    u(0,x) &= u_0(x), && x \in (x_{min}, x_{max}), \\
    \partial_t u(0,x) &= v_0(x), && x \in (x_{min}, x_{max}), \\
    \text{boundary conditions}, &&& x \in \partial (x_{min}, x_{max}).
\end{aligned}
```

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
includes a pre-built semidiscretization of this equation:
[`WaveEquationNonperiodicSemidiscretization`](@ref).
Have a look at the source code if you want to dig deeper.
In particular, you can find applications of
[`derivative_left`](@ref), [`derivative_right`](@ref)
[`mul_transpose_derivative_left!`](@ref), and [`mul_transpose_derivative_right!`](@ref).
Below is an example demonstrating how to use this semidiscretization.


```@example wave_equation
using SummationByPartsOperators, OrdinaryDiffEq
using LaTeXStrings; using Plots: Plots, plot, plot!, savefig

# general parameters
xmin = -1.
xmax = +1.
tspan = (0., 8.0)
u0_func(x) = exp(-20x^2)
v0_func(x) = zero(x)
# HomogeneousNeumann, HomogeneousDirichlet, and NonReflecting BCs are available
left_bc  = Val(:HomogeneousNeumann)
right_bc = Val(:HomogeneousDirichlet)

# setup spatial semidiscretization
D2 = derivative_operator(MattssonSvärdShoeybi2008(), derivative_order=2,
                         accuracy_order=4, xmin=xmin, xmax=xmax, N=101)
semi = WaveEquationNonperiodicSemidiscretization(D2, left_bc, right_bc)
ode = semidiscretize(v0_func, u0_func, semi, tspan)

# solve second-order ODE using a Runge-Kutta-Nyström method
sol = solve(ode, DPRKN6(), saveat=range(first(tspan), stop=last(tspan), length=200))

# visualize the result
plot(xguide=L"x")
plot!(evaluate_coefficients(sol[end].x[2], semi), label=L"u")
plot!(evaluate_coefficients(sol[end].x[1], semi), label=L"\partial_t u")
savefig("example_wave_equation.png");
```

![](example_wave_equation.png)


## Advanced visualization of different boundary conditions

Let's create animations of the numerical solutions for different
boundary conditions.

```julia
using Printf; using Plots: Animation, frame, gif

function create_gif(left_bc::Val{LEFT_BC}, right_bc::Val{RIGHT_BC}) where {LEFT_BC, RIGHT_BC}
    xmin = -1.
    xmax = +1.
    tspan = (0., 8.0)
    u0_func(x) = exp(-20x^2)
    v0_func(x) = zero(x)

    D2 = derivative_operator(MattssonSvärdShoeybi2008(), derivative_order=2,
                            accuracy_order=4, xmin=xmin, xmax=xmax, N=101)
    semi = WaveEquationNonperiodicSemidiscretization(D2, left_bc, right_bc)
    ode = semidiscretize(v0_func, u0_func, semi, tspan)

    sol = solve(ode, DPRKN6(), saveat=range(first(tspan), stop=last(tspan), length=200))

    anim = Animation()
    idx = 1
    x, u = evaluate_coefficients(sol[idx].x[2], D2)
    fig = plot(x, u, xguide=L"x", yguide=L"u", xlim=extrema(x), ylim=(-1.05, 1.05),
              label="", title=@sprintf("\$t = %6.2f \$", sol.t[idx]))
    for idx in 1:length(sol.t)
        fig[1] = x, sol.u[idx].x[2]
        plot!(title=@sprintf("\$t = %6.2f \$", sol.t[idx]))
        frame(anim)
    end
    gif(anim, "wave_equation_$(LEFT_BC)_$(RIGHT_BC).gif")
end

create_gif(Val(:HomogeneousNeumann), Val(:HomogeneousNeumann))
```

![wave_equation_HomogeneousNeumann_HomogeneousNeumann](https://user-images.githubusercontent.com/12693098/119228021-3603f800-bb11-11eb-9703-157503308ec3.gif)

```julia
create_gif(Val(:HomogeneousNeumann), Val(:HomogeneousDirichlet))
```

![wave_equation_HomogeneousNeumann_HomogeneousDirichlet](https://user-images.githubusercontent.com/12693098/119228026-3a301580-bb11-11eb-8354-de23104fe285.gif)

```julia
create_gif(Val(:HomogeneousNeumann), Val(:NonReflecting))
```

![wave_equation_HomogeneousNeumann_NonReflecting](https://user-images.githubusercontent.com/12693098/119228041-5633b700-bb11-11eb-9c17-bc56c906dae3.gif)
# Linear advection equation with constant coefficients

This tutorial is concerned with the basic linear advection equation

```math
\begin{aligned}
    \partial_t u(t,x) + \partial_x u(t,x) &= 0, && t \in (0,T), x \in (x_{min}, x_{max}), \\
    u(0,x) &= u_0(x), && x \in (x_{min}, x_{max}), \\
    u(t,x_{min}) &= u_L(t).
\end{aligned}
```

Note that the advection velocity is positive (unity). Thus, a boundary condition
needs to be specified exactly at the left boundary. Otherwise, the problem will
not be well-posed (under-specified or over-specified).

## Basic example using finite difference SBP operators

Let's create an appropriate discretization of this equation step by step. At first,
we load packages that we will use in this example.

```@example linear_advection
using SummationByPartsOperators, OrdinaryDiffEq
using LaTeXStrings; using Plots: Plots, plot, plot!, savefig
```

Next, we specify the initial and boundary data as Julia functions as well as the
spatial domain and the time span.

```@example linear_advection
xmin, xmax = -1.0, 1.0
u0_func(x) = sinpi(x)
uL_func(t) = t >= 3 ? sinpi(t) : zero(t)
tspan = (0., 8.0)
```

This choice of the domain and boundary condition ensures that the initial profile
is transported out of the domain before non-homogeneous boundary data influences
the solution.

Next, we implement the semidiscretization using the interface of
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)
which is part of [DifferentialEquations.jl](https://diffeq.sciml.ai/latest/).

```@example linear_advection
function rhs!(du, u, params, t)
  D = params.D

  # Set `du = - D * u` using in-place multiplication avoiding allocations
  # for efficiency
  mul!(du, D, u, -one(eltype(D)))

  # Next, we impose the boundary conditions weakly using an SAT at the left
  # boundary. Since we use the strong form of the equation, we do not need to
  # do anything at the right boundary.
  # Assuming that boundary nodes are included in the grid, adding this SAT
  # can be achieved by
  du[begin] += (uL_func(t) - u[begin]) / left_boundary_weight(D)

  return nothing
end
```

Here, we have used a simultaneous approximation term (SAT) to impose the boundary
condition weakly. In general, this approach is related to the weak imposition
of boundary conditions using numerical fluxes in finite volume and discontinuous
Galerkin methods; they are even equivalent for the linear advection equation
considered here.

Next, we choose an SBP operator `D`, evaluate the initial data on the grid, and
set up the semidiscretization as an ODE problem.

```@example linear_advection
D = derivative_operator(MattssonNordström2004(), derivative_order=1, accuracy_order=4,
                        xmin=xmin, xmax=xmax, N=101)
u0 = compute_coefficients(u0_func, D)
params = (D=D, )
ode = ODEProblem(rhs!, u0, tspan, params);
```

Finally, we can solve the ODE using an explicit Runge-Kutta method with adaptive
time stepping.

```@example linear_advection
sol = solve(ode, Tsit5(), saveat=range(first(tspan), stop=last(tspan), length=200));

plot(xguide=L"x", yguide=L"u")
plot!(evaluate_coefficients(sol[1], D), label=L"u_0")
plot!(evaluate_coefficients(sol[end], D), label=L"u_\mathrm{numerical}")
savefig("example_linear_advection.png");
```

![](example_linear_advection.png)

## Advanced visualization

Let's create an animation of the numerical solution.

```julia
using Printf; using Plots: Animation, frame, gif

let anim = Animation()
    idx = 1
    x, u = evaluate_coefficients(sol[idx], D)
    fig = plot(x, u, xguide=L"x", yguide=L"u", xlim=extrema(x), ylim=(-1.05, 1.05),
              label="", title=@sprintf("\$t = %6.2f \$", sol.t[idx]))
    for idx in 1:length(sol.t)
        fig[1] = x, sol.u[idx]
        plot!(title=@sprintf("\$t = %6.2f \$", sol.t[idx]))
        frame(anim)
    end
    gif(anim, "example_linear_advection.gif")
end
```

![example_linear_advection_animation](https://user-images.githubusercontent.com/12693098/119224994-7bb8c480-bb01-11eb-9c3e-c4fea709da71.gif)

## Continuous and discontinuous Galerkin methods

You can use a CG or DG method by swapping out the derivative operator `D`.

```@example linear_advection
plot(xguide=L"x", yguide=L"u")
plot!(evaluate_coefficients(sol[1], D), label=L"u_0")
plot!(evaluate_coefficients(sol[end], D), label=L"u_\mathrm{FD}")

# CGSEM using polynomials of degree 3, i.e. 4 nodes per element, and 30 elements
D_CGSEM = couple_continuously(
            legendre_derivative_operator(xmin=-1.0, xmax=1.0, N=4),
            UniformMesh1D(xmin=xmin, xmax=xmax, Nx=30))
ode_CGSEM = ODEProblem(rhs!, compute_coefficients(u0_func, D_CGSEM), tspan, (D=D_CGSEM,))
sol_CGSEM = solve(ode_CGSEM, Tsit5(), save_everystep=false)
plot!(evaluate_coefficients(sol_CGSEM[end], D_CGSEM), label=L"u_\mathrm{CG}")

# DGSEM using polynomials of degree 3, i.e. 4 nodes per element, and 30 elements
# which are coupled using upwind fluxes
D_DGSEM = couple_discontinuously(
            legendre_derivative_operator(xmin=-1.0, xmax=1.0, N=4),
            UniformMesh1D(xmin=xmin, xmax=xmax, Nx=30),
            Val(:minus))
ode_DGSEM = ODEProblem(rhs!, compute_coefficients(u0_func, D_DGSEM), tspan, (D=D_DGSEM,))
sol_DGSEM = solve(ode_DGSEM, Tsit5(), save_everystep=false)
plot!(evaluate_coefficients(sol_DGSEM[end], D_DGSEM), label=L"u_\mathrm{DG}")

savefig("example_linear_advection_Galerkin.png");
```

![](example_linear_advection_Galerkin.png)
# Linear advection equation with variable coefficients

This tutorial is concerned with the linear advection equation

```math
\begin{aligned}
    \partial_t u(t,x) + \partial_x (a(x) u(t,x)) &= 0, && t \in (0,T), x \in (x_{min}, x_{max}), \\
    u(0,x) &= u_0(x), && x \in (x_{min}, x_{max}), \\
    \text{boundary conditions}, &&& x \in \partial (x_{min}, x_{max})
\end{aligned}
```

with variable coefficient ``a``.

The boundary conditions depend on the sign of the transport velocity ``a``
at the boundary. In particular, specifying a Dirichlet type boundary condition
is only allowed for inflow boundaries, e.g. ``a(x_{min}) > 0`` at ``x = x_{min}``.

[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
includes a pre-built semidiscretization of this equation:
[`VariableLinearAdvectionNonperiodicSemidiscretization`](@ref).
Have a look at the source code if you want to dig deeper. Below is an example
demonstrating how to use this semidiscretization.

```@example variable_linear_advection
using SummationByPartsOperators, OrdinaryDiffEq
using LaTeXStrings; using Plots: Plots, plot, plot!, savefig

# general parameters
xmin = -1.
xmax = +1.
tspan = (0., 8.0)
afunc(x) = one(x)
u0func(x) = sinpi(x)
# Dirichlet type boundary conditions; they are used only at inflow boundaries
left_bc(t) = t >= 3 ? sinpi(t) : zero(t)
right_bc(t) = zero(t)

# discretization parameters
interior_order = 4
N = 101
# whether a split form should be applied or not
split_form = Val(false)

# setup spatial semidiscretization
D = derivative_operator(MattssonSvärdShoeybi2008(), 1, interior_order, xmin, xmax, N)
# whether or not artificial dissipation should be applied: nothing, dissipation_operator(D)
Di = nothing
semi = VariableLinearAdvectionNonperiodicSemidiscretization(D, Di, afunc, split_form, left_bc, right_bc)
ode = semidiscretize(u0func, semi, tspan)

# solve ODE
sol = solve(ode, SSPRK104(), dt=D.Δx, adaptive=false,
            save_everystep=false)

# visualise the result
plot(xguide=L"x", yguide=L"u")
plot!(evaluate_coefficients(sol[1], semi), label=L"u_0")
plot!(evaluate_coefficients(sol[end], semi), label=L"u_\mathrm{numerical}")
savefig("example_linear_advection.png");
```

![](example_linear_advection.png)
# Linear advection diffusion equation with periodic boundary conditions

Let's consider the linear advection diffusion equation

```math
\begin{aligned}
    \partial_t u(t,x) + a \partial_x u(t,x) &= \varepsilon \partial_x^2 u(t,x), && t \in (0,T), x \in (x_{min}, x_{max}), \\
    u(0,x) &= u_0(x), && x \in (x_{min}, x_{max}), \\
\end{aligned}
```

with periodic boundary conditions. Here, `a` is the constant advection velocity
and `ε > 0` is the constant diffusion coefficient.

## Basic example using finite difference SBP operators

Let's create an appropriate discretization of this equation step by step. At first,
we load packages that we will use in this example.

```@example advection_diffusion
using SummationByPartsOperators, OrdinaryDiffEq
using LaTeXStrings; using Plots: Plots, plot, plot!, savefig
```

Next, we specify the initial data as Julia function as well as the
spatial domain and the time span.

```@example advection_diffusion
xmin, xmax = -1.0, 1.0
u0_func(x) = sinpi(x)
tspan = (0., 10.0)
```

Next, we implement the semidiscretization using the interface of
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)
which is part of [DifferentialEquations.jl](https://diffeq.sciml.ai/latest/).

```@example advection_diffusion
function advection_diffusion!(du, u, params, t)
    # In-place version of du = -a * D1 * u
    mul!(du, params.D1, u, -params.a)
    # In-place version of du = du + ε * D2 * u
    mul!(du, params.D2, u, params.ε, true)
end
```

Next, we choose first- and second-derivative SBP operators `D1, D2`, evaluate
the initial data on the grid, and set up the semidiscretization as an ODE problem.

```@example advection_diffusion
N = 100 # number of grid points
D1 = periodic_derivative_operator(derivative_order=1, accuracy_order=4,
                                  xmin=xmin, xmax=xmax, N=N)
D2 = periodic_derivative_operator(derivative_order=2, accuracy_order=4,
                                  xmin=xmin, xmax=xmax, N=N)
u0 = u0_func.(grid(D1))
params = (D1=D1, D2=D2, a=1.0, ε=0.03)
ode = ODEProblem(advection_diffusion!, u0, tspan, params);
```

Finally, we can solve the ODE using an explicit Runge-Kutta method with adaptive
time stepping.

```@example advection_diffusion
sol = solve(ode, Tsit5(), saveat=range(first(tspan), stop=last(tspan), length=200));

plot(xguide=L"x", yguide=L"u")
plot!(evaluate_coefficients(sol[1], D1), label=L"u_0")
plot!(evaluate_coefficients(sol[end], D1), label=L"u_\mathrm{numerical}")
savefig("example_advection_diffusion.png");
```

![](example_advection_diffusion.png)


## Advanced visualization

Let's create an animation of the numerical solution.

```julia
using Printf; using Plots: Animation, frame, gif

let anim = Animation()
    idx = 1
    x, u = evaluate_coefficients(sol[idx], D1)
    fig = plot(x, u, xguide=L"x", yguide=L"u", xlim=extrema(x), ylim=(-1.05, 1.05),
              label="", title=@sprintf("\$t = %6.2f \$", sol.t[idx]))
    for idx in 1:length(sol.t)
        fig[1] = x, sol.u[idx]
        plot!(title=@sprintf("\$t = %6.2f \$", sol.t[idx]))
        frame(anim)
    end
    gif(anim, "example_advection_diffusion.gif")
end
```

![example_advection_diffusion_animation](https://user-images.githubusercontent.com/12693098/119226459-7b242c00-bb09-11eb-848b-d09590aa1c31.gif)

