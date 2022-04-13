# Contributing to `Manopt.jl`

First, thanks for taking the time to contribute.
Any contribution is appreciated and welcome.

The following is a set of guidelines to [`Manopt.jl`](https://juliamanifolds.github.io/Manopt.jl/).

#### Table of Contents

- [Contributing to `Manopt.jl`](#Contributing-to-manoptjl)
      - [Table of Contents](#Table-of-Contents)
  - [I just have a question](#I-just-have-a-question)
  - [How can I file an issue?](#How-can-I-file-an-issue)
  - [How can I contribute?](#How-can-I-contribute)
    - [Add a missing method](#Add-a-missing-method)
    - [Provide a new algorithm](#Provide-a-new-algorithm)
    - [Provide a new example](#Provide-a-new-example)
    - [Code style](#Code-style)

## I just have a question

The developer can most easily be reached in the Julia Slack channel [#manifolds](https://julialang.slack.com/archives/CP4QF0K5Z).
You can apply for the Julia Slack workspace [here](https://julialang.org/slack/) if you haven't joined yet.
You can also ask your question on [discourse.julialang.org](https://discourse.julialang.org).

## How can I file an issue?

If you found a bug or want to propose a feature, we track our issues within the [GitHub repository](https://github.com/JuliaManifolds/Manopt.jl/issues).

## How can I contribute?

### Add a missing method

There is still a lot of methods for within the optimisation framework of  `Manopt.jl`, may it be functions, gradients, differentials, proximal maps, step size rules or stopping criteria.
If you notice a method missing and can contribute an implementation, please do so!
Even providing a single new method is a good contribution.

### Provide a new algorithm

A main contribution you can provide is another algorithm that is not yet included in the
package.
An alorithm is always based on a is a concrete type of a [`Problem`](https://manoptjl.org/stable/plans/index.html#Problems-1) storing the main information of the task and a concrete type of an [`Option`](https://manoptjl.org/stable/plans/index.html#Options-1from) storing all information that needs to be known to the solver in general. The actual algorithm is split into an initialization phase, see [`initialize_solver!`](https://manoptjl.org/stable/solvers/index.html#Manopt.initialize_solver!), and the implementation of the `i`th step of the sovler itself, see  before the iterative procedure, see [`step_solver!`](https://manoptjl.org/stable/solvers/index.html#Manopt.step_solver!).
For these two functions it would be great if a new algorithm uses functions from the [`ManifoldsBase.jl`](https://juliamanifolds.github.io/Manifolds.jl/latest/interface.html) interface as generic as possible. For example, if possible use [`retract!(M,q,p,X)`](https://juliamanifolds.github.io/Manifolds.jl/latest/interface.html#ManifoldsBase.retract!-Tuple{AbstractManifold,Any,Any,Any}) in favour of [`exp!(M,q,p,X)`](https://juliamanifolds.github.io/Manifolds.jl/latest/interface.html#ManifoldsBase.exp!-Tuple{AbstractManifold,Any,Any,Any}) to perform a step starting in `p` in direction `X` (in place of `q`), since the exponential map might be too expensive to evaluate or might not be available on a certain manifold. See [Retractions and inverse retractions](https://juliamanifolds.github.io/Manifolds.jl/latest/interface.html#Retractions-and-inverse-Retractions) for more details.
Further, if possible, prefer [`retract!(M,q,p,X)`](https://juliamanifolds.github.io/Manifolds.jl/latest/interface.html#ManifoldsBase.retract!-Tuple{AbstractManifold,Any,Any,Any}) in favour of [`retract(M,p,X)`](https://juliamanifolds.github.io/Manifolds.jl/latest/interface.html#ManifoldsBase.retract-Tuple{AbstractManifold,Any,Any}), since a computation in place of a suitable variable `q` reduces memory allocations.

Usually, the methods implemented in `Manopt.jl` also have a high-level interface, that is easier to call, creates the necessary problem and options structure and calls the solver.

The two technical functions `initialize_solver!` and `step_solver!` should be documented with technical details, while the high level interface should usually provide a general description and some literature references to the algorithm at hand.

### Provide a new example

The `examples/` folder features several examples covering all solvers. Still, if you have a new example that you implemented yourself for fun or for a paper, feel free to add it to the repository as well. Also if you have a [Pluto](https://github.com/fonsp/Pluto.jl) notebook of your example, feel free to contribute that.

### Code style

We try to follow the [documentation guidelines](https://docs.julialang.org/en/v1/manual/documentation/) from the Julia documentation as well as [Blue Style](https://github.com/invenia/BlueStyle).
We run [`JuliaFormatter.jl`](https://github.com/domluna/JuliaFormatter.jl) on the repo in the way set in the `.JuliaFormatter.toml` file, which enforces a number of conventions consistent with the Blue Style.

We also follow a few internal conventions:

- It is preferred that the `Problem`'s struct contains information about the general structure of the problem
- Any implemented function should be accompanied by its mathematical formulae if a closed form exists.
- Problem and option structures are stored within the `plan/` folder and sorted by properties of the problem and/or solver at hand
- Within the source code of one algorithm, the high level interface should be first, then the initialisation, then the step.
- Otherwise an alphabetical order is preferrable.
- The above implies that the mutating variant of a function follows the non-mutating variant.
- There should be no dangling `=` signs.
- Always add a newline between things of different types (struct/method/const).
- Always add a newline between methods for different functions (including mutating/nonmutating variants).
- Prefer to have no newline between methods for the same function; when reasonable, merge the docstrings.
- All `import`/`using`/`include` should be in the main module file.
# Manopt.jl

Optimization on Manifolds.

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://manoptjl.org/stable)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![CI](https://github.com/JuliaManifolds/Manopt.jl/workflows/CI/badge.svg)](https://github.com/JuliaManifolds/Manopt.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![codecov](https://codecov.io/gh/JuliaManifolds/Manopt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaManifolds/Manopt.jl)
[![DOI](https://zenodo.org/badge/74746729.svg)](https://zenodo.org/badge/latestdoi/74746729)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03866/status.svg)](https://doi.org/10.21105/joss.03866)

For a function f that maps from a [Riemannian manifold](https://en.wikipedia.org/wiki/Riemannian_manifold)
ℳ onto the real line, we aim to solve

> Find the minimizer x on ℳ, i.e. the (or a) point where f attains its minimum.

`Manopt.jl` provides a framework for optimization on manifolds.
Based on [Manopt](https://manopt.org) and
[MVIRT](https://ronnybergmann.net/mvirt/), both implemented in Matlab,
this toolbox aims to provide an easy access to optimization methods on manifolds
for [Julia](https://julialang.org), including example data and visualization methods.

## Getting started

In Julia you can get started by just typing

```julia
] add Manopt
```

then checkout the [Get Started: Optimize!](https://manoptjl.org/stable/tutorials/MeanAndMedian.html) tutorial or the
[examples](https://github.com/JuliaManifolds/Manopt.jl/tree/master/examples)
in this repository.

## Citation

If you use `Manopt.jl`in your work, please cite the following

```biblatex
@article{Bergmann2022,
    Author    = {Ronny Bergmann},
    Doi       = {10.21105/joss.03866},
    Journal   = {Journal of Open Source Software},
    Number    = {70},
    Pages     = {3866},
    Publisher = {The Open Journal},
    Title     = {Manopt.jl: Optimization on Manifolds in {J}ulia},
    Volume    = {7},
    Year      = {2022},
}
```

To refer to a certain version or the source code in general we recommend to cite for example

```biblatex
@software{manoptjl-zenodo-mostrecent,
    Author = {Ronny Bergmann},
    Copyright = {MIT License},
    Doi = {10.5281/zenodo.4290905},
    Publisher = {Zenodo},
    Title = {Manopt.jl},
    Year = {2022},
}
```

for the most recent version or a corresponding version specific DOI, see [the list of all versions](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%224290905%22&sort=-version&all_versions=True).
Note that both citations are in [BibLaTeX](https://ctan.org/pkg/biblatex) format.
# Notation

In this package, we follow the notation introduced in [Manifolds.jl – Notation](https://juliamanifolds.github.io/Manifolds.jl/latest/misc/notation.html)

with the following additional or slightly changed notation

| Symbol | Description | Also used | Comment |
|:--:|:--------------- |:--:|:-- |
| ``∇`` | The [Levi-Cevita connection](https://en.wikipedia.org/wiki/Levi-Civita_connection) | | |
| ``\operatorname{grad}f`` | The Riemannian gradient | ``∇f``| due to possible confusion with the connection, we try to avoid ``∇f`` |
| ``\operatorname{Hess}f``| The Riemannian Hessian | |
| ``p,q``| points on a manifold | | when definiting functions |
| ``x_k, y_k`` | points on a manifold | | iterates of an algorithm |
# Table of Contents, Types and Functions

This page lists all pages of this documentations, all available types and functions.

## Complete List of Contents

```@contents
Depth = 3
```

## Available Types

```@index
Modules = [Manopt]
Order   = [:type]
```

## Solver Functions

```@index
Modules = [Manopt]
Pages = ["plans/index.md", "solvers/index.md"]
```

## Functions

```@index
Modules = [Manopt]
Pages = ["functions/adjointDifferentials.md", "functions/costs.md", "functions/differentials.md", "functions/gradients.md", "functions/jacobiFields.md", "functions/proximalMaps.md"]
```
# About

Manopt.jl inherited its name from [Manopt](https://manopt.org), a Matlab toolbox.
It is currently Maintained by [Ronny Bergmann](https://ronnybergmann.net/about.html) (manopt@ronnybergmann.net) with contributions from Tom Christian Riemer, who implemented the [trust regions](@ref trust_regions) solver.

If you want to contribute a manifold or algorithm or have any questions, visit
the [GitHub repository](https://github.com/JuliaManifolds/Manopt.jl/)
to clone/fork the repository or open an issue.
# Welcome to Manopt.jl

```@meta
CurrentModule = Manopt
```

```@docs
Manopt.Manopt
```

For a function $f:\mathcal M → ℝ$ defined on a [Riemannian manifold](https://en.wikipedia.org/wiki/Riemannian_manifold) $\mathcal M$ we aim to solve

$\operatorname*{argmin}_{x ∈ \mathcal M} f(x),$

or in other words: find the point $x$ on the manifold, where $f$ reaches its minimal function value.

`Manopt.jl` provides a framework for optimization on manifolds.
Based on [Manopt](https://manopt.org) and
[MVIRT](https://ronnybergmann.net/mvirt/), both implemented in Matlab,
this toolbox provide an easy access to optimization methods on manifolds
for [Julia](https://julialang.org), including example data and visualization methods.

If you want to delve right into `Manopt.jl` check out the
[Get started: Optimize!](@ref Optimize) tutorial.

`Manopt.jl` makes it easy to use an algorithm for your favorite
manifold as well as a manifold for your favorite algorithm. It already provides
many manifolds and algorithms, which can easily be enhanced, for example to
[record](@ref RecordOptions) certain data or
[display information](@ref DebugOptions) throughout iterations.

If you use `Manopt.jl`in your work, please cite the following

```biblatex
@article{Bergmann2022,
    Author    = {Ronny Bergmann},
    Doi       = {10.21105/joss.03866},
    Journal   = {Journal of Open Source Software},
    Number    = {70},
    Pages     = {3866},
    Publisher = {The Open Journal},
    Title     = {Manopt.jl: Optimization on Manifolds in {J}ulia},
    Volume    = {7},
    Year      = {2022},
}
```

To refer to a certain version or the source code in general we recommend to cite for example

```biblatex
@software{manoptjl-zenodo-mostrecent,
    Author = {Ronny Bergmann},
    Copyright = {MIT License},
    Doi = {10.5281/zenodo.4290905},
    Publisher = {Zenodo},
    Title = {Manopt.jl},
    Year = {2022},
}
```

for the most recent version or a corresponding version specific DOI, see [the list of all versions](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%224290905%22&sort=-version&all_versions=True).
Note that both citations are in [BibLaTeX](https://ctan.org/pkg/biblatex) format.

## Main Features

### Functions on Manifolds

Several functions are available, implemented on an arbitrary manifold, [cost functions](@ref CostFunctions), [differentials](@ref DifferentialFunctions), and [gradients](@ref GradientFunctions) as well as [proximal maps](@ref proximalMapFunctions), but also several [jacobi Fields](@ref JacobiFieldFunctions) and their [adjoints](@ref adjointDifferentialFunctions).

### Optimization Algorithms (Solvers)

For every optimization algorithm, a [solver](@ref Solvers) is implemented based on a [`Problem`](@ref) that describes the problem to solve and its [`Options`](@ref) that set up the solver, store interims values. Together they
form a [plan](@ref planSection).

### Visualization

To visualize and interpret results, `Manopt.jl` aims to provide both easy plot functions as well as [exports](@ref Exports). Furthermore a system to get [debug](@ref DebugOptions) during the iterations of an algorithms as well as [record](@ref RecordOptions) capabilities, i.e. to record a specified tuple of values per iteration, most prominently [`RecordCost`](@ref) and
[`RecordIterate`](@ref). Take a look at the [Get started: Optimize!](@ref Optimize) tutorial how to easily activate this.

## Manifolds

This project is build upon [ManifoldsBase.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/interface.html), a generic interface to implement manifolds. Certain functions are extended for specific manifolds from [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/), but all other manifolds from that package can be used here, too.

The notation in the documentation aims to follow the same [notation](https://juliamanifolds.github.io/Manifolds.jl/stable/notation.html) from these packages.

## Literature

If you want to get started with manifolds, one book is [[do Carmo, 1992](#doCarmo1992)],
and if you want do directly dive into optimization on manifolds, my favourite reference is
[[Absil, Mahony, Sepulchre, 2008](#AbsilMahonySepulchre2008)], which is also available
online for free.

```@raw html
<ul>
<li id="AbsilMahonySepulchre2008">
    [<a>Absil, Mahony, Sepulchre, 2008</a>]
    P.-A. Absil, R. Mahony and R. Sepulchre,
    <emph>Optimization Algorithms on Matrix Manifolds</emph>,
    Princeton University Press, 2008,
    doi: <a href="https://doi.org/10.1515/9781400830244">10.1515/9781400830244</a>,
    <a href="http://press.princeton.edu/chapters/absil/">open access</a>.
</li>
<li id="doCarmo1992">
    [<a>doCarmo, 1992</a>]
    M. P. do Carmo,
    <emph>Riemannian Geometry</emph>,
    Birkhäuser Boston, 1992,
    ISBN: 0-8176-3490-8.
</li>
</ul>
```
# [Stochastic Gradient Descent](@id StochasticGradientDescentSolver)

```@meta
CurrentModule = Manopt
```

```@docs
stochastic_gradient_descent
stochastic_gradient_descent!
```

## Options

```@docs
StochasticGradientDescentOptions
```

Additionally, the options share a [`DirectionUpdateRule`](@ref),
so you can also apply [`MomentumGradient`](@ref) and [`AverageGradient`](@ref) here.
The most inner one should always be.

```@docs
AbstractStochasticGradientProcessor
StochasticGradient
```
# [Nelder Mead Method](@id NelderMeadSolver)

```@meta
CurrentModule = Manopt
```

```@docs
    NelderMead
    NelderMead!
```

## Options

```@docs
    NelderMeadOptions
```
# [Steihaug-Toint Truncated Conjugate-Gradient Method](@id tCG)

The aim is to solve the trust-region subproblem

```math
\operatorname*{arg\,min}_{η  ∈  T_{x}\mathcal{M}} m_{x}(η) = F(x) +
⟨\operatorname{grad}F(x), η⟩_{x} + \frac{1}{2} ⟨
\operatorname{Hess}F(x)[η], η⟩_{x}
```

```math
\text{s.t.} \; ⟨η, η⟩_{x} \leq {Δ}^2
```

on a manifold by using the Steihaug-Toint truncated conjugate-gradient method.
All terms involving the trust-region radius use an inner product w.r.t. the
preconditioner; this is because the iterates grow in length w.r.t. the
preconditioner, guaranteeing that we do not re-enter the trust-region.

## Initialization

Initialize ``η_0 = η`` if using randomized approach and
``η`` the zero tangent vector otherwise, ``r_0 = \operatorname{grad}F(x)``,
``z_0 = \operatorname{P}(r_0)``, ``δ_0 = z_0`` and ``k=0``

## Iteration

Repeat until a convergence criterion is reached

1. Set ``κ = ⟨δ_k, \operatorname{Hess}F(x)[δ_k]⟩_x``,
    ``α =\frac{⟨r_k, z_k⟩_x}{κ}`` and
    ``⟨η_k, η_k⟩_{x}^* = ⟨η_k, \operatorname{P}(η_k)⟩_x +
    2α ⟨η_k, \operatorname{P}(δ_k)⟩_{x} +  {α}^2
    ⟨δ_k, \operatorname{P}(δ_k)⟩_{x}``.
2. If ``κ ≤ 0`` or ``⟨η_k, η_k⟩_x^* ≥ Δ^2``
    return ``η_{k+1} = η_k + τ δ_k`` and stop.
3. Set ``η_{k}^*= η_k + α δ_k``, if
    ``⟨η_k, η_k⟩_{x} + \frac{1}{2} ⟨η_k,
    \operatorname{Hess}[F] (η_k)_{x}⟩_{x} ≤ ⟨η_k^*,
    η_k^*⟩_{x} + \frac{1}{2} ⟨η_k^*,
    \operatorname{Hess}[F] (η_k)_ {x}⟩_{x}``
    set ``η_{k+1} = η_k`` else set ``η_{k+1} = η_{k}^*``.
4. Set ``r_{k+1} = r_k + α \operatorname{Hess}[F](δ_k)_x``,
     ``z_{k+1} = \operatorname{P}(r_{k+1})``,
     ``β = \frac{⟨r_{k+1},  z_{k+1}⟩_{x}}{⟨r_k, z_k
   ⟩_{x}}`` and ``δ_{k+1} = -z_{k+1} + β δ_k``.
5. Set ``k=k+1``.

## Result

The result is given by the last computed ``η_k``.

## Remarks

The ``\operatorname{P}(⋅)`` denotes the symmetric, positive deﬁnite
preconditioner. It is required if a randomized approach is used i.e. using
a random tangent vector ``η_0`` as initial
vector. The idea behind it is to avoid saddle points. Preconditioning is
simply a rescaling of the variables and thus a redeﬁnition of the shape of
the trust region. Ideally ``\operatorname{P}(⋅)`` is a cheap, positive
approximation of the inverse of the Hessian of ``F`` at ``x``. On
default, the preconditioner is just the identity.

To step number 2: Obtain ``τ`` from the positive root of
``\left\lVert η_k + τ δ_k \right\rVert_{\operatorname{P}, x} = Δ``
what becomes after the conversion of the equation to

````math
 τ = \frac{-⟨η_k, \operatorname{P}(δ_k)⟩_{x} +
 \sqrt{⟨η_k, \operatorname{P}(δ_k)⟩_{x}^{2} +
 ⟨δ_k, \operatorname{P}(δ_k)⟩_{x} ( Δ^2 -
 ⟨η_k, \operatorname{P}(η_k)⟩_{x})}}
 {⟨δ_k, \operatorname{P}(δ_k)⟩_{x}}.
````

It can occur that ``⟨δ_k, \operatorname{Hess}[F] (δ_k)_{x}⟩_{x}
= κ ≤ 0`` at iteration ``k``. In this case, the model is not strictly
convex, and the stepsize ``α =\frac{⟨r_k, z_k⟩_{x}}
{κ}`` computed in step 1. does not give a reduction in the modelfunction
``m_x(⋅)``. Indeed, ``m_x(⋅)`` is unbounded from below along the
line ``η_k + α δ_k``. If our aim is to minimize the model within
the trust-region, it makes far more sense to reduce ``m_x(⋅)`` along
``η_k + α δ_k`` as much as we can while staying within the
trust-region, and this means moving to the trust-region boundary along this
line. Thus when ``κ ≤ 0`` at iteration k, we replace ``α =
\frac{⟨r_k, z_k⟩_{x}}{κ}`` with ``τ`` described as above.
The other possibility is that ``η_{k+1}`` would lie outside the trust-region at
iteration k (i.e. ``⟨η_k, η_k⟩_{x}^{* }
≥ {Δ}^2`` what can be identified with the norm of ``η_{k+1}``). In
particular, when ``\operatorname{Hess}[F] (⋅)_{x}`` is positive deﬁnite
and ``η_{k+1}`` lies outside the trust region, the solution to the
trust-region problem must lie on the trust-region boundary. Thus, there
is no reason to continue with the conjugate gradient iteration, as it
stands, as subsequent iterates will move further outside the trust-region
boundary. A sensible strategy, just as in the case considered above, is to
move to the trust-region boundary by ﬁnding ``τ``.

## Interface

```@docs
  truncated_conjugate_gradient_descent
  truncated_conjugate_gradient_descent!
```

## Options

```@docs
TruncatedConjugateGradientOptions
```

## Additional Stopping Criteria

```@docs
StopIfResidualIsReducedByPower
StopIfResidualIsReducedByFactor
StopWhenTrustRegionIsExceeded
StopWhenCurvatureIsNegative
StopWhenModelIncreased
```
# [The Riemannian Chambolle-Pock Algorithm](@id ChambollePockSolver)

The Riemannian Chambolle–Pock is a generalization of the Chambolle–Pock algorithm[^ChambollePock2011].
It is also known as primal dual hybrig gradient (PDHG) or primal dual proximal splitting (PDPS) algorithm.

In order to minimize over $p∈\mathcal M§ the cost function consisting of

```math
F(p) + G(Λ(p)),
```

where $F:\mathcal M → \overline{ℝ}$, $G:\mathcal N → \overline{ℝ}$, and
$Λ:\mathcal M →\mathcal N$.
If the manifolds $\mathcal M$ or $\mathcal N$ are not Hadamard, it has to be considered locally,
i.e. on geodesically convex sets $\mathcal C \subset \mathcal M$ and $\mathcal D \subset\mathcal N$
such that $Λ(\mathcal C) \subset \mathcal D$.

The algorithm is available in four variants: exact versus linearized (see `variant`)
as well as with primal versus dual relaxation (see `relax`). For more details, see
[^BergmannHerzogSilvaLouzeiroTenbrinckVidalNunez2020].
In the following we note the case of the exact, primal relaxed Riemannian Chambolle–Pock algorithm.

Given base points $m∈\mathcal C$, $n=Λ(m)∈\mathcal D$,
initial primal and dual values $p^{(0)} ∈\mathcal C$, $ξ_n^{(0)} ∈T_n^*\mathcal N$,
and primal and dual step sizes $\sigma_0$, $\tau_0$, relaxation $\theta_0$,
as well as acceleration $\gamma$.

As an initialization, perform $\bar p^{(0)} \gets p^{(0)}$.

The algorithms performs the steps $k=1,…,$ (until a [`StoppingCriterion`](@ref) is fulfilled with)

1. ```math
   ξ^{(k+1)}_n = \operatorname{prox}_{\tau_k G_n^*}\Bigl(ξ_n^{(k)} + \tau_k \bigl(\log_n Λ (\bar p^{(k)})\bigr)^\flat\Bigr)
   ```
2. ```math
   p^{(k+1)} = \operatorname{prox}_{\sigma_k F}\biggl(\exp_{p^{(k)}}\Bigl( \operatorname{PT}_{p^{(k)}\gets m}\bigl(-\sigma_k DΛ(m)^*[ξ_n^{(k+1)}]\bigr)^\sharp\Bigr)\biggr)
   ```
3. Update
   * ``\theta_k = (1+2\gamma\sigma_k)^{-\frac{1}{2}}``
   * ``\sigma_{k+1} = \sigma_k\theta_k``
   * ``\tau_{k+1} =  \frac{\tau_k}{\theta_k}``
4. ```math
   \bar p^{(k+1)}  = \exp_{p^{(k+1)}}\bigl(-\theta_k \log_{p^{(k+1)}} p^{(k)}\bigr)
   ```

Furthermore you can exchange the exponential map, the logarithmic map, and the parallel transport
by a retraction, an in verse retraction and a vector transport.

Finally you can also update the base points $m$ and $n$ during the iterations.
This introduces a few additional vector transports. The same holds for the case that
$Λ(m^{(k)})\neq n^{(k)}$ at some point. All these cases are covered in the algorithm.

```@meta
CurrentModule = Manopt
```

```@docs
ChambollePock
ChambollePock!
```

## Problem & Options

```@docs
PrimalDualOptions
ChambollePockOptions
```

## Useful Terms

```@docs
primal_residual
dual_residual
```

## Debug

```@docs
DebugDualBaseIterate
DebugDualBaseChange
DebugPrimalBaseIterate
DebugPrimalBaseChange
DebugDualChange
DebugDualIterate
DebugDualResidual
DebugPrimalChange
DebugPrimalIterate
DebugPrimalResidual
DebugPrimalDualResidual
```

## Record

```@docs
RecordDualBaseIterate
RecordDualBaseChange
RecordDualChange
RecordDualIterate
RecordPrimalBaseIterate
RecordPrimalBaseChange
RecordPrimalChange
RecordPrimalIterate
```

## Internals

```@docs
Manopt.update_prox_parameters!
```

[^ChambollePock2011]:
    > A. Chambolle, T. Pock:
    > _A first-order primal-dual algorithm for convex problems with applications to imaging_,
    > Journal of Mathematical Imaging and Vision 40(1), 120–145, 2011.
    > doi: [10.1007/s10851-010-0251-1](https://dx.doi.org/10.1007/s10851-010-0251-1)# [Riemannian quasi-Newton methods](@id quasiNewton)

```@meta
    CurrentModule = Manopt
```

```@docs
    quasi_Newton
    quasi_Newton!
```

## Background

The aim is to minimize a real-valued function on a Riemannian manifold, i.e.

```math
\min f(x), \quad x ∈ \mathcal{M}.
```

Riemannian quasi-Newtonian methods are as generalizations of their Euclidean counterparts Riemannian line search methods. These methods determine a search direction ``η_k ∈ T_{x_k} \mathcal{M}`` at the current iterate ``x_k`` and a suitable stepsize ``α_k`` along ``\gamma(α) = R_{x_k}(α η_k)``, where ``R: T \mathcal{M} →\mathcal{M}`` is a retraction. The next iterate is obtained by

```math
x_{k+1} = R_{x_k}(α_k η_k).
```

In quasi-Newton methods, the search direction is given by

```math
η_k = -{\mathcal{H}_k}^{-1}[\operatorname{grad}f (x_k)] = -\mathcal{B}_k [\operatorname{grad} (x_k)],
```

where ``\mathcal{H}_k : T_{x_k} \mathcal{M} →T_{x_k} \mathcal{M}`` is a positive definite self-adjoint operator, which approximates the action of the Hessian ``\operatorname{Hess} f (x_k)[⋅]`` and ``\mathcal{B}_k = {\mathcal{H}_k}^{-1}``. The idea of quasi-Newton methods is instead of creating a complete new approximation of the Hessian operator ``\operatorname{Hess} f(x_{k+1})`` or its inverse at every iteration, the previous operator ``\mathcal{H}_k`` or ``\mathcal{B}_k`` is updated by a convenient formula using the obtained information about the curvature of the objective function during the iteration. The resulting operator ``\mathcal{H}_{k+1}`` or ``\mathcal{B}_{k+1}`` acts on the tangent space ``T_{x_{k+1}} \mathcal{M}`` of the freshly computed iterate ``x_{k+1}``.
In order to get a well-defined method, the following requirements are placed on the new operator ``\mathcal{H}_{k+1}`` or ``\mathcal{B}_{k+1}`` that is created by an update. Since the Hessian ``\operatorname{Hess} f(x_{k+1})`` is a self-adjoint operator on the tangent space ``T_{x_{k+1}} \mathcal{M}``, and ``\mathcal{H}_{k+1}`` approximates it, we require that ``\mathcal{H}_{k+1}`` or ``\mathcal{B}_{k+1}`` is also self-adjoint on ``T_{x_{k+1}} \mathcal{M}``. In order to achieve a steady descent, we want ``η_k`` to be a descent direction in each iteration. Therefore we require, that ``\mathcal{H}_{k+1}`` or ``\mathcal{B}_{k+1}`` is a positive definite operator on ``T_{x_{k+1}} \mathcal{M}``. In order to get information about the curvature of the objective function into the new operator ``\mathcal{H}_{k+1}`` or ``\mathcal{B}_{k+1}``, we require that it satisfies a form of a Riemannian quasi-Newton equation:

```math
\mathcal{H}_{k+1} [T_{x_k \rightarrow x_{k+1}}({R_{x_k}}^{-1}(x_{k+1}))] = \operatorname{grad}(x_{k+1}) - T_{x_k \rightarrow x_{k+1}}(\operatorname{grad}f(x_k))
```

or

```math
\mathcal{B}_{k+1} [\operatorname{grad}f(x_{k+1}) - T_{x_k \rightarrow x_{k+1}}(\operatorname{grad}f(x_k))] = T_{x_k \rightarrow x_{k+1}}({R_{x_k}}^{-1}(x_{k+1}))
```

where ``T_{x_k \rightarrow x_{k+1}} : T_{x_k} \mathcal{M} →T_{x_{k+1}} \mathcal{M}`` and the chosen retraction ``R`` is the associated retraction of ``T``. We note that, of course, not all updates in all situations will meet these conditions in every iteration.
For specific quasi-Newton updates, the fulfilment of the Riemannian curvature condition, which requires that

```math
g_{x_{k+1}}(s_k, y_k) > 0
```

holds, is a requirement for the inheritance of the self-adjointness and positive definiteness of the ``\mathcal{H}_k`` or ``\mathcal{B}_k`` to the operator ``\mathcal{H}_{k+1}`` or ``\mathcal{B}_{k+1}``. Unfortunately, the fulfillment of the Riemannian curvature condition is not given by a step size ``\alpha_k > 0`` that satisfies the generalised Wolfe conditions. However, in order to create a positive definite operator ``\mathcal{H}_{k+1}`` or ``\mathcal{B}_{k+1}`` in each iteration, in [^HuangGallivanAbsil2015] the so-called locking condition was introduced, which requires that the isometric vector transport ``T^S``, which is used in the update formula, and its associate retraction ``R`` fulfill

```math
T^{S}{x, ξ_x}(ξ_x) = β T^{R}{x, ξ_x}(ξ_x), \quad β = \frac{\lVert ξ_x \rVert_x}{\lVert T^{R}{x, ξ_x}(ξ_x) \rVert_{R_{x}(ξ_x)}},
```

where ``T^R`` is the vector transport by differentiated retraction. With the requirement that the isometric vector transport ``T^S`` and its associated retraction ``R`` satisfies the locking condition and using the tangent vector

```math
y_k = {β_k}^{-1} \operatorname{grad}f(x_{k+1}) - T^{S}{x_k, α_k η_k}(\operatorname{grad}f(x_k)),
```

where

```math
β_k = \frac{\lVert α_k η_k \rVert_{x_k}}{\lVert T^{R}{x_k, α_k η_k}(α_k η_k) \rVert_{x_{k+1}}},
```

in the update, it can be shown that choosing a stepsize ``α_k > 0`` that satisfies the Riemannian Wolfe conditions leads to the fulfilment of the Riemannian curvature condition, which in turn implies that the operator generated by the updates is positive definite.
In the following we denote the specific operators in matrix notation and hence use ``H_k`` and ``B_k``, respectively.

## Direction Updates

In general there are different ways to compute a fixed [`AbstractQuasiNewtonUpdateRule`](@ref).
In general these are represented by

```@docs
AbstractQuasiNewtonDirectionUpdate
QuasiNewtonMatrixDirectionUpdate
QuasiNewtonLimitedMemoryDirectionUpdate
QuasiNewtonCautiousDirectionUpdate
```

## Hessian Update Rules

Using

```@docs
update_hessian!
```

the following update formulae for either ``H_{k+1}`` or `` B_{k+1}`` are available.

```@docs
AbstractQuasiNewtonUpdateRule
BFGS
DFP
Broyden
SR1
InverseBFGS
InverseDFP
InverseBroyden
InverseSR1
```

## Options

The quasi Newton algorithm is based on a [`GradientProblem`](@ref).

```@docs
QuasiNewtonOptions
```

## Literature
# [Alternating Gradient Descent](@id AlternatingGradientDescentSolver)

```@meta
CurrentModule = Manopt
```

```@docs
alternating_gradient_descent
alternating_gradient_descent!
```

## Problem

```@docs
AlternatingGradientProblem
```

## Options

```@docs
AlternatingGradientDescentOptions
```

Additionally, the options share a [`DirectionUpdateRule`](@ref),
which chooses the current component, so they can be decorated further;
The most inner one should always be the following one though.

```@docs
AlternatingGradient
```
# [Subgradient Method](@id SubgradientSolver)

```@docs
subgradient_method
subgradient_method!
```

## Options

```@docs
SubGradientMethodOptions
```

For [`DebugAction`](@ref)s and [`RecordAction`](@ref)s to record (sub)gradient,
its norm and the step sizes, see the [steepest Descent](@ref GradientDescentSolver)
actions.

# [Conjugate Gradient Descent](@id CGSolver)

```@meta
CurrentModule = Manopt
```

```@docs
conjugate_gradient_descent
conjugate_gradient_descent!
```

## Options

```@docs
ConjugateGradientDescentOptions
```

## Available Coefficients

The update rules act as [`DirectionUpdateRule`](@ref), which internally always first evaluate the gradient itself.

```@docs
ConjugateDescentCoefficient
DaiYuanCoefficient
FletcherReevesCoefficient
HagerZhangCoefficient
HeestenesStiefelCoefficient
LiuStoreyCoefficient
PolakRibiereCoefficient
SteepestDirectionUpdateRule
```

# Literature
# [The Riemannian Trust-Regions Solver](@id trust_regions)

The aim is to solve an optimization problem on a manifold

```math
\operatorname*{min}_{x  ∈  \mathcal{M}} F(x)
```

by using the Riemannian trust-regions solver. It is number one choice for smooth
optimization. This trust-region method uses the Steihaug-Toint truncated
conjugate-gradient method [`truncated_conjugate_gradient_descent`](@ref)
to solve the inner minimization problem called the
trust-regions subproblem. This inner solve can be preconditioned by providing
a preconditioner (symmetric and positive deﬁnite, an approximation of the
inverse of the Hessian of ``F``). If no Hessian of the cost function ``F`` is
provided, a standard approximation of the Hessian based on the gradient
``\operatorname{grad}F`` with [`ApproxHessianFiniteDifference`](@ref) will be computed.

## Initialization

Initialize ``x_0 = x`` with an initial point ``x`` on the manifold. It can be
given by the caller or set randomly. Set the initial trust-region radius
``\Delta =\frac{1}{8} \bar{\Delta}`` where ``\bar{\Delta}`` is the maximum radius
the trust-region can have. Usually one uses
the root of the manifold dimension ``\operatorname{dim}(\mathcal{M})``.
For accepting the next iterate and evaluating the new trust-region radius one
needs an accept/reject threshold ``\rho'  ∈  [0,\frac{1}{4})``, which is
``\rho' = 0.1`` on default. Set ``k=0``.

## Iteration

Repeat until a convergence criterion is reached

1. Set ``η`` as a random tangent vector if using randomized approach. Else
    set ``η`` as the zero vector in the tangential space ``T_{x_k}\mathcal{M}``.
2. Set ``η^*`` as the solution of the trust-region subproblem, computed by
    the tcg-method with ``η`` as initial vector.
3. If using randomized approach compare ``η^*`` with the Cauchy point
    ``η_{c}^* = -\tau_{c} \frac{\Delta}{\lVert \operatorname{Grad}[F] (x_k) \rVert_{x_k}} \operatorname{Grad}[F] (x_k)`` by the model function ``m_{x_k}(⋅)``. If the
    model decrease is larger by using the Cauchy point, set
    ``η^* = η_{c}^*``.
4. Set ``{x}^* = \operatorname{retr}_{x_k}(η^*)``.
5. Set ``\rho = \frac{F(x_k)-F({x}^*)}{m_{x_k}(η)-m_{x_k}(η^*)}``, where
    ``m_{x_k}(⋅)`` describes the quadratic model function.
6. Update the trust-region radius:``\Delta = \begin{cases}\frac{1}{4} \Delta &\text{ if } \rho < \frac{1}{4} \, \text{or} \, m_{x_k}(η)-m_{x_k}(η^*) \leq 0 \, \text{or}  \, \rho = \pm  ∈ fty , \\\operatorname{min}(2 \Delta, \bar{\Delta}) &\text{ if } \rho > \frac{3}{4} \, \text{and the tcg-method stopped because of negative curvature or exceeding the trust-region},\\\Delta & \, \text{otherwise.}\end{cases}``
7. If ``m_{x_k}(η)-m_{x_k}(η^*) \geq 0`` and ``\rho > \rho'`` set
    ``x_k = {x}^*``.
8. Set ``k = k+1``.

## Result

The result is given by the last computed ``x_k``.

## Remarks

To the Initialization: A random point on the manifold.

To step number 1: Using randomized approach means using a random tangent
vector as initial vector for the approximal solve of the trust-regions
subproblem. If this is the case, keep in mind that the vector must be in the
trust-region radius. This is achieved by multiplying
`η` by `sqrt(4,eps(Float64))` as long as
its norm is greater than the current trust-region radius ``\Delta``.
For not using randomized approach, one can get the zero tangent vector.

To step number 2: Obtain ``η^*`` by (approximately) solving the
trust-regions subproblem

```math
\operatorname*{arg\,min}_{η  ∈  T_{x_k}\mathcal{M}} m_{x_k}(η) = F(x_k) +
\langle \operatorname{grad}F(x_k), η \rangle_{x_k} + \frac{1}{2} \langle
\operatorname{Hess}[F](η)_ {x_k}, η \rangle_{x_k}
```

```math
\text{s.t.} \; \langle η, η \rangle_{x_k} \leq {\Delta}^2
```

with the Steihaug-Toint truncated conjugate-gradient (tcg) method. The problem
as well as the solution method is described in the
[`truncated_conjugate_gradient_descent`](@ref).

To step number 3: If using a random tangent vector as an initial vector, compare
the result of the tcg-method with the Cauchy point. Convergence proofs assume
that one achieves at least (a fraction of) the reduction of the Cauchy point.
The idea is to go in the direction of the gradient to an optimal point. This
can be on the edge, but also before.
The parameter ``\tau_{c}`` for the optimal length is defined by

```math
\tau_{c} = \begin{cases} 1 & \langle \operatorname{Grad}[F] (x_k), \,
\operatorname{Hess}[F] (η_k)_ {x_k}\rangle_{x_k} \leq 0 , \\
\operatorname{min}(\frac{{\operatorname{norm}(\operatorname{Grad}[F] (x_k))}^3}
{\Delta \langle \operatorname{Grad}[F] (x_k), \,
\operatorname{Hess}[F] (η_k)_ {x_k}\rangle_{x_k}}, 1) & \, \text{otherwise.}
\end{cases}
```

To check the model decrease one compares

```math
m_{x_k}(η_{c}^*) = F(x_k) + \langle η_{c}^*,
\operatorname{Grad}[F] (x_k)\rangle_{x_k} + \frac{1}{2}\langle η_{c}^*,
\operatorname{Hess}[F] (η_{c}^*)_ {x_k}\rangle_{x_k}
```

with

```math
m_{x_k}(η^*) = F(x_k) + \langle η^*,
\operatorname{Grad}[F] (x_k)\rangle_{x_k} + \frac{1}{2}\langle η^*,
\operatorname{Hess}[F] (η^*)_ {x_k}\rangle_{x_k}.
```

If ``m_{x_k}(η_{c}^*) < m_{x_k}(η^*)`` then ``m_{x_k}(η_{c}^*)`` is the better choice.

To step number 4: ``\operatorname{retr}_{x_k}(⋅)`` denotes the retraction, a
mapping ``\operatorname{retr}_{x_k}:T_{x_k}\mathcal{M} \rightarrow \mathcal{M}``
wich approximates the exponential map. In some cases it is cheaper to use this
instead of the exponential.

To step number 6: One knows that the [`truncated_conjugate_gradient_descent`](@ref) algorithm stopped for
these reasons when the stopping criteria [`StopWhenCurvatureIsNegative`](@ref),
[`StopWhenTrustRegionIsExceeded`](@ref) are activated.

To step number 7: The last step is to decide if the new point ``{x}^*`` is
accepted.

## Interface

```@docs
trust_regions
trust_regions!
```

## Options

```@docs
AbstractHessianOptions
TrustRegionsOptions
```

## Approximation of the Hessian

```@docs
ApproxHessianFiniteDifference
```
# [Douglas–Rachford Algorithm](@id DRSolver)

The (Parallel) Douglas–Rachford ((P)DR) Algorithm was generalized to Hadamard
manifolds in [[Bergmann, Persch, Steidl, 2016](#BergmannPerschSteidl2016)].

The aim is to minimize the sum

$F(x) = f(x) + g(x)$

on a manifold, where the two summands have proximal maps
$\operatorname{prox}_{λ f}, \operatorname{prox}_{λ g}$ that are easy
to evaluate (maybe in closed form or not too costly to approximate).
Further define the Reflection operator at the proximal map as

$\operatorname{refl}_{λ f}(x) = \exp_{\operatorname{prox}_{λ f}(x)} \bigl( -\log_{\operatorname{prox}_{λ f}(x)} x \bigr)$.

Let $\alpha_k ∈  [0,1]$ with $\sum_{k ∈ \mathbb N} \alpha_k(1-\alpha_k) =  ∈ fty$
and $λ > 0$ which might depend on iteration $k$ as well) be given.

Then the (P)DRA algorithm for initial data $x_0 ∈ \mathcal H$ as

## Initialization

Initialize $t_0 = x_0$ and $k=0$

## Iteration

Repeat  until a convergence criterion is reached

1. Compute $s_k = \operatorname{refl}_{λ f}\operatorname{refl}_{λ g}(t_k)$
2. within that operation store $x_{k+1} = \operatorname{prox}_{λ g}(t_k)$ which is the prox the inner reflection reflects at.
3. Compute $t_{k+1} = g(\alpha_k; t_k, s_k)$
4. Set $k = k+1$

## Result

The result is given by the last computed $x_K$.

For the parallel version, the first proximal map is a vectorial version, where
in each component one prox is applied to the corresponding copy of $t_k$ and
the second proximal map corresponds to the indicator function of the set,
where all copies are equal (in $\mathcal H^n$, where $n$ is the number of copies),
leading to the second prox being the Riemannian mean.

## Interface

```@docs
  DouglasRachford
  DouglasRachford!
```

## Options

```@docs
DouglasRachfordOptions
```

For specific [`DebugAction`](@ref)s and [`RecordAction`](@ref)s see also
[Cyclic Proximal Point](@ref CPPSolver).

## Literature

```@raw html
<ul>
<li id="BergmannPerschSteidl2016">[<a>Bergmann, Persch, Steidl, 2016</a>]
  Bergmann, R; Persch, J.; Steidl, G.: <emph>A Parallel Douglas–Rachford
  Algorithm for Minimizing ROF-like Functionals on Images with Values in
  Symmetric Hadamard Manifolds.</emph>
  SIAM Journal on Imaging Sciences, Volume 9, Number 3, pp. 901–937, 2016.
  doi: <a href="https://doi.org/10.1137/15M1052858">10.1137/15M1052858</a>,
  arXiv: <a href="https://arxiv.org/abs/1512.02814">1512.02814</a>.
</li>
</ul>
```

# Solvers

```@meta
CurrentModule = Manopt
```

Solvers can be applied to [`Problem`](@ref)s with solver
specific [`Options`](@ref).

# List of Algorithms

The following algorithms are currently available

| Solver  | File   | Problem & Option  |
----------|--------|-------------------|
[Alternating Gradient Descent](@ref AlternatingGradientDescentSolver) | `alterating_gradient_descent.jl` | [`AlternatingGradientProblem`](@ref), [`AlternatingGradientDescentOptions`](@ref)
[Chambolle-Pock](@ref ChambollePockSolver) | `Chambolle-Pock.jl` | [`PrimalDualProblem`](@ref), [`ChambollePockOptions`](@ref)
[Cyclic Proximal Point](@ref CPPSolver) | `cyclic_proximal_point.jl` | [`ProximalProblem`](@ref), [`CyclicProximalPointOptions`](@ref)
[Douglas–Rachford](@ref DRSolver) | `DouglasRachford.jl` | [`ProximalProblem`](@ref), [`DouglasRachfordOptions`](@ref)
[Gradient Descent](@ref GradientDescentSolver) | `gradient_descent.jl` |  [`GradientProblem`](@ref), [`GradientDescentOptions`](@ref)
[Nelder-Mead](@ref NelderMeadSolver) | `NelderMead.jl` | [`CostProblem`](@ref), [`NelderMeadOptions`](@ref)
[Particle Swarm](@ref ParticleSwarmSolver) | `particle_swarm.jl` | [`CostProblem`](@ref), [`ParticleSwarmOptions`](@ref)
[Quasi-Newton Method](@ref quasiNewton) | `quasi_newton.jl`| [`GradientProblem`](@ref), [`QuasiNewtonOptions`](@ref)
[Subgradient Method](@ref SubgradientSolver) | `subgradient_method.jl` | [`SubGradientProblem`](@ref), [`SubGradientMethodOptions`](@ref)
[Steihaug-Toint Truncated Conjugate-Gradient Method](@ref tCG) | `truncated_conjugate_gradient_descent.jl` | [`HessianProblem`](@ref), [`TruncatedConjugateGradientOptions`](@ref)
[The Riemannian Trust-Regions Solver](@ref trust_regions) | `trust_regions.jl` | [`HessianProblem`](@ref), [`TrustRegionsOptions`](@ref)

Note that the solvers (or their [`Options`](@ref) to be precise) can also be decorated to enhance your algorithm by general additional properties, see [Decorated Solvers](@ref DecoratedSolvers).

## [StoppingCriteria](@id StoppingCriteria)

Stopping criteria are implemented as a `functor`, i.e. inherit from the base type

```@docs
StoppingCriterion
StoppingCriterionSet
```

```@autodocs
Modules = [Manopt]
Pages = ["plans/stopping_criterion.jl"]
Order = [:type]
```

as well as the functions

```@docs
Base.:&(::StoppingCriterion, ::StoppingCriterion)
Base.:|(::StoppingCriterion, ::StoppingCriterion)
get_reason
get_stopping_criteria
get_active_stopping_criteria
are_these_stopping_critera_active
```

further stopping criteria might be available for individual Solvers.

## [Decorated Solvers](@id DecoratedSolvers)

The following decorators are available.

### [Debug Solver](@id DebugSolver)

The decorator to print debug during the iterations can be activated by
decorating the [`Options`](@ref) with [`DebugOptions`](@ref) and implementing
your own [`DebugAction`](@ref)s.
For example printing a gradient from the [`GradientDescentOptions`](@ref) is
automatically available, as explained in the [`gradient_descent`](@ref) solver.

```@autodocs
Modules = [Manopt]
Pages   = ["debug_solver.jl"]
```

### [Record Solver](@id RecordSolver)

The decorator to record certain values during the iterations can be activated by
decorating the [`Options`](@ref) with [`RecordOptions`](@ref) and implementing
your own [`RecordAction`](@ref)s.
For example recording the gradient from the [`GradientDescentOptions`](@ref) is
automatically available, as explained in the [`gradient_descent`](@ref) solver.

```@autodocs
Modules = [Manopt]
Pages   = ["record_solver.jl"]
```

## Technical Details

 The main function a solver calls is

```@docs
solve(p::Problem, o::Options)
```

which is a framework, that you in general should not change or redefine.
It uses the following methods, which also need to be implemented on your own
algorithm, if you want to provide one.

```@docs
initialize_solver!
step_solver!
get_solver_result
stop_solver!(p::Problem, o::Options, i::Int)
```
# [Particle Swarm Optimization](@id ParticleSwarmSolver)

```@meta
CurrentModule = Manopt
```

```@docs
  particle_swarm
  particle_swarm!
```

## Options

```@docs
ParticleSwarmOptions
```

## Literature# [Cyclic Proximal Point](@id CPPSolver)

The Cyclic Proximal Point (CPP) algorithm is a [Proximal Problem](@ref ProximalProblem).

It aims to minimize

```math
F(x) = \sum_{i=1}^c f_i(x)
```

assuming that the [proximal maps](@ref proximalMapFunctions) $\operatorname{prox}_{λ f_i}(x)$
are given in closed form or can be computed efficiently (at least approximately).

The algorithm then cycles through these proximal maps, where the type of cycle
might differ and the proximal parameter $λ_k$ changes after each cycle $k$.

For a convergence result on
[Hadamard manifolds](https://en.wikipedia.org/wiki/Hadamard_manifold)
see [[Bačák, 2014](#Bačák2014)].

```@docs
cyclic_proximal_point
cyclic_proximal_point!
```

## Options

```@docs
CyclicProximalPointOptions
```

## Debug Functions

```@docs
DebugProximalParameter
```

## Record Functions

```@docs
RecordProximalParameter
```

## Literature

```@raw html
<ul>
<li id="Bačák2014">[<a>Bačák, 2014</a>]
  Bačák, M: <emph>Computing Medians and Means in Hadamard Spaces.</emph>,
  SIAM Journal on Optimization, Volume 24, Number 3, pp. 1542–1566,
  doi: <a href="https://doi.org/10.1137/140953393">10.1137/140953393</a>,
  arxiv: <a href="https://arxiv.org/abs/1210.2145">1210.2145</a>.
  </li>
</ul>
```
# [Gradient Descent](@id GradientDescentSolver)

```@meta
CurrentModule = Manopt
```

```@docs
  gradient_descent
  gradient_descent!
```

## Options

```@docs
AbstractGradientOptions
GradientDescentOptions
```

## Direction Update Rules

A field of the options is the `direction`, a [`DirectionUpdateRule`](@ref), which by default [`IdentityUpdateRule`](@ref) just evaluates the gradient but can be enhanced for example to

```@docs
DirectionUpdateRule
IdentityUpdateRule
MomentumGradient
AverageGradient
Nesterov
```

## Debug Actions

```@docs
DebugGradient
DebugGradientNorm
DebugStepsize
```

## Record Actions

```@docs
RecordGradient
RecordGradientNorm
RecordStepsize
```
# [Bézier curves](@id BezierCurves)

```@autodocs
Modules = [Manopt]
Pages   = ["bezier_curves.jl"]
```

## Literature# Specific manifold functions

This small section extends the functions available from [ManifoldsBase.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/interface.html) and [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/), espcially a few random generators, that are simpler than the functions available.

```@autodocs
Modules = [Manopt]
Pages   = ["manifold_functions.jl"]
```

## Simplified random functions

While statistics are available in [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/index.html),
the following functions provide default random points and vectors on manifolds.

```@autodocs
Modules = [Manopt]
Pages   = ["random.jl"]
```

## Initialize data

```@autodocs
Modules = [Manopt]
Pages   = ["initialize_data.jl"]
```
# [Jacobi Fields](@id JacobiFieldFunctions)

A smooth tangent vector field $J: [0,1] → T\mathcal M$
along a geodesic $g(⋅;x,y)$ is called _Jacobi field_
if it fulfills the ODE

$\displaystyle 0 = \frac{D}{dt}J + R(J,\dot g)\dot g,$

where $R$ is the Riemannian curvature tensor.
Such Jacobi fields can be used to derive closed forms for the exponential map,
the logarithmic map and the geodesic, all of them with respect to both arguments:
Let $F:\mathcal N → \mathcal M$ be given (for the $\exp_x⋅$
  we have $\mathcal N = T_x\mathcal M$, otherwise $\mathcal N=\mathcal M$) and denote by
$ξ_1,…,ξ_d$ an orthonormal frame along $g(⋅;x,y)$ that diagonalizes
the curvature tensor with corresponding eigenvalues $κ_1,…,κ_d$.
Note that on symmetric manifolds such a frame always exists.

Then $DF(x)[η] = \sum_{k=1}^d \langle η,ξ_k(0)\rangle_xβ(κ_k)ξ_k(T)$ holds,
where $T$ also depends on the function $F$ as the weights $β$. The values
stem from solving the corresponding system of (decoupled) ODEs.

Note that in different references some factors might be a little different,
for example when using unit speed geodesics.

The following weights functions are available

```@autodocs
Modules = [Manopt]
Pages   = ["Jacobi_fields.jl"]
```
# [Cost Functions](@id CostFunctions)

The following cost functions are available

```@autodocs
Modules = [Manopt]
Pages   = ["costs.jl"]
```
# [Adjoint Differentials](@id adjointDifferentialFunctions)

```@autodocs
Modules = [Manopt]
Pages   = ["adjoint_differentials.jl"]
```
# [Differentials](@id DifferentialFunctions)

```@autodocs
Modules = [Manopt]
Pages   = ["functions/differentials.jl"]
```
# Functions

There are several functions required within optimization, most prominently
[costFunctions](@ref CostFunctions) and [gradients](@ref GradientFunctions). This package includes
several cost functions and corresponding gradients, but also corresponding
[proximal maps](@ref proximalMapFunctions) for variational methods
manifold-valued data. Most of these functions require the evaluation of
[Differential](@ref DifferentialFunctions)s or their `AdjointDifferential`s as well
as [JacobiFields](@ref JacobiFieldFunctions) (e.g. easily to evaluate for symmetric manifolds).# [Gradients](@id GradientFunctions)

For a function $f:\mathcal M→ℝ$
the Riemannian gradient $\operatorname{grad}f(x)$ at $x∈\mathcal M$
is given by the unique tangent vector fulfilling

$\langle \operatorname{grad}f(x), ξ\rangle_x = D_xf[ξ],\quad
\forall ξ ∈ T_x\mathcal M,$
where $D_xf[ξ]$ denotes the differential of $f$ at $x$ with respect to
the tangent direction (vector) $ξ$ or in other words the directional
derivative.

This page collects the available gradients.

```@autodocs
Modules = [Manopt]
Pages   = ["gradients.jl"]
```
# [Proximal Maps](@id proximalMapFunctions)

For a function $\varphi:\mathcal M →ℝ$ the proximal map is defined
as

$\displaystyle\operatorname{prox}_{λ\varphi}(x)
= \operatorname*{argmin}_{y ∈ \mathcal M} d_{\mathcal M}^2(x,y) + \varphi(y),
\quad λ > 0,$

where $d_{\mathcal M}: \mathcal M \times \mathcal M → ℝ$ denotes
the geodesic distance on \(\mathcal M\). While it might still be difficult to
compute the minimizer, there are several proximal maps known (locally) in closed
form. Furthermore if $x^{\star} ∈ \mathcal M$ is a minimizer of $\varphi$, then

$\displaystyle\operatorname{prox}_{λ\varphi}(x^\star) = x^\star,$

i.e. a minimizer is a fixed point of the proximal map.

This page lists all proximal maps available within Manopt. To add you own, just
extend the `functions/proximal_maps.jl` file.

```@autodocs
Modules = [Manopt]
Pages   = ["proximal_maps.jl"]
```
# [Plans for solvers](@id planSection)

```@meta
CurrentModule = Manopt
```

In order to start a solver, both a [`Problem`](@ref) and [`Options`](@ref) are required.
Together they form a __plan__ and these are stored in this folder. For
sub-problems there are maybe also only [`Options`](@ref), since they than refer to the
same problem.

## Options

For most algorithms a certain set of options can either be
generated beforehand of the function with keywords can be used.
Generally the type

```@docs
Options
get_options
```

Since the `Options` directly relate to a solver, they are documented with the
corresponding [Solvers](@ref).
You can always access the options (since they
might be decorated) by calling [`get_options`](@ref).

### Decorators for Options

Options can be decorated using the following trait and function to initialize

```@docs
dispatch_options_decorator
is_options_decorator
decorate_options
```

In general decorators often perform actions so we introduce

```@docs
AbstractOptionsAction
```

as well as a helper for storing values using keys, i.e.

```@docs
StoreOptionsAction
get_storage
has_storage
update_storage!
```

#### [Debug Options](@id DebugOptions)

```@autodocs
Modules = [Manopt]
Pages = ["plans/debug_options.jl"]
Order = [:type, :function]
```

see [DebugSolver](@ref DebugSolver) for details on the decorated solver.

Further specific [`DebugAction`](@ref)s can be found at the specific Options.

#### [Record Options](@id RecordOptions)

```@autodocs
Modules = [Manopt]
Pages = ["plans/record_options.jl"]
Order = [:type, :function]
Private = false
```

```@docs
getindex(ro::RecordOptions, s::Symbol)
getindex(::RecordGroup,::Any...)
```

see [RecordSolver](@ref RecordSolver) for details on the decorated solver.

Further specific [`RecordAction`](@ref)s can be found at the specific Options.

there's one internal helper that might be useful for you own actions, namely

```@docs
record_or_reset!
```

### [Stepsize and Linesearch](@id Stepsize)

The step size determination is implemented as a `Functor` based on

```@docs
Stepsize
```

in general there are

```@autodocs
Modules = [Manopt]
Pages = ["plans/stepsize.jl"]
Order = [:type,:function]
```

## Problems

A problem usually contains its cost function and provides and
implementation to access the cost

```@docs
Problem
get_cost
```

A problem can be of different type, more specifically, whether its containing functions,
for example to compute the gradient work with allocation or without. To be precise, an
allocation function `X = gradF(x)` allocates memory for its result `X`, while `gradF!(X,x)` does not.

```@docs
AbstractEvaluationType
AllocatingEvaluation
MutatingEvaluation
```

### Cost based problem

```@docs
CostProblem
```

### Gradient based problem

```@docs
AbstractGradientProblem
GradientProblem
StochasticGradientProblem
get_gradient
get_gradients
```

### Subgradient based problem

```@docs
SubGradientProblem
get_subgradient
```

### [Proximal Map(s) based problem](@id ProximalProblem)

```@docs
ProximalProblem
get_proximal_map
```

### [Hessian based problem](@id HessianProblem)

```@docs
HessianProblem
get_hessian
get_preconditioner
```

### [Primal dual based problem](@id PrimalDualProblem)

```@docs
PrimalDualProblem
get_primal_prox
get_dual_prox
forward_operator
linearized_forward_operator
adjoint_linearized_operator
```
# [Error Measures](@id ErrorMeasures)

```@docs
meanSquaredError
meanAverageError
```
# [Exports](@id Exports)

Exports aim to provide a consistent generation of images of your results. For example if you [record](@ref RecordOptions) the trace your algorithm walks on the [Sphere](https://juliamanifolds.github.io/Manifolds.jl/stable/manifolds/sphere.html), you yan easily export this trace to a rendered image using [`asymptote_export_S2_signals`](@ref) and render the result with [Asymptote](https://sourceforge.net/projects/asymptote/).
Despite these, you can always [record](@ref RecordOptions) values during your iterations,
and export these, for example to `csv`.

## Asymptote

The following functions provide exports both in graphics and/or raw data using [Asymptote](https://sourceforge.net/projects/asymptote/).

```@autodocs
Modules = [Manopt]
Pages   = ["Asymptote.jl"]
```

# [Data](@id Data)

For some manifolds there are artificial or real application data available
that can be loaded using the following data functions

```@autodocs
Modules = [Manopt]
Pages   = ["artificialDataFunctions.jl"]
```---
title: 'Manopt.jl: Optimization on Manifolds in Julia'
tags:
  - Julia
  - Riemannian manifolds
  - optimization
  - numerical analysis
authors:
  - name: Ronny Bergmann
    orcid: 0000-0001-8342-7218
    affiliation: 1
affiliations:
 - name: Norwegian University of Science and Technology, Department of Mathematical Sciences, Trondheim, Norway
   index: 1
date: 22 July 2021
bibliography: bibliography.bib

---

# Summary

[`Manopt.jl`](https://manoptjl.org) provides a set of optimization algorithms for optimization problems given on a Riemannian manifold $\mathcal M$.
Based on a generic optimization framework, together with the interface [`ManifoldsBase.jl`](https://github.com/JuliaManifolds/ManifoldsBase.jl) for Riemannian manifolds, classical and recently developed methods are provided in an efficient implementation. Algorithms include the derivative-free Particle Swarm and Nelder–Mead algorithms, as well as classical gradient, conjugate gradient and stochastic gradient descent. Furthermore, quasi-Newton methods like a Riemannian L-BFGS [@HuangGallivanAbsil:2015:1] and nonsmooth optimization algorithms like a Cyclic Proximal Point Algorithm [@Bacak:2014:1], a (parallel) Douglas-Rachford algorithm [@BergmannPerschSteidl:2016:1] and a Chambolle-Pock algorithm [@BergmannHerzogSilvaLouzeiroTenbrinckVidalNunez:2021:1] are provided, together with several basic cost functions, gradients and proximal maps as well as debug and record capabilities.

# Statement of Need

In many applications and optimization tasks, non-linear data appears naturally.
For example, when data on the sphere is measured [@GousenbourgerMassartMusolasAbsilJaquesHendrickxMarzouk:2017], diffusion data can be captured as a signal or even multivariate data of symmetric positive definite matrices [@ValkonenBrediesKnoll2013], and orientations like they appear for electron backscattered diffraction (EBSD) data [@BachmannHielscherSchaeben2011]. Another example are fixed rank matrices, appearing in matrix completion [@Vandereycken:2013:1].
Working on these data, for example doing data interpolation and approximation [@BergmannGousenbourger:2018:2], denoising [@LellmannStrekalovskiyKoetterCremers:2013:1; @BergmannFitschenPerschSteidl:2018], inpainting [@BergmannChanHielscherPerschSteidl:2016], or performing matrix completion [@GaoAbsil:2021], can usually be phrased as an optimization problem

$$ \text{Minimize}\quad f(x) \quad \text{where } x\in\mathcal M, $$

where the optimization problem is phrased on a Riemannian manifold $\mathcal M$.

A main challenge of these algorithms is that, compared to the (classical) Euclidean case, there is no addition available. For example, on the unit sphere $\mathbb S^2$ of unit vectors in $\mathbb R^3$, adding two vectors of unit length yields a vector that is not of unit norm.
The solution is to generalize the notion of a shortest path from the straight line to what is called a (shortest) geodesic, or acceleration-free curve.
Similarly, other features and properties also have to be rephrased and generalized when performing optimization on a Riemannian manifold.
Algorithms to perform the optimization can still often be stated in a generic way, i.e. on an arbitrary Riemannian manifold $\mathcal M$.
Further examples and a thorough introduction can be found in @AbsilMahonySepulchre:2008:1; @Boumal:2020:1.

For a user facing an optimization problem on a manifold, there are two obstacles to the actual numerical optimization: firstly, a suitable implementation of the manifold at hand is required, for example how to evaluate the above-mentioned geodesics; and secondly, an implementation of the optimization algorithm that employs said methods from the manifold, such that the algorithm can be applied to the cost function $f$ a user already has.

Using the interface for manifolds from the `ManifoldsBase.jl` package, the algorithms are implemented in the optimization framework. They can then be used with any manifold from [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/) [@AxenBaranBergmannRzecki:2021:1], a library of efficiently-implemented Riemannian manifolds.
`Manopt.jl` provides a low-bar entry to optimization on manifolds, while also providing efficient implementations, that can easily be extended to cover  manifolds specified by the user.

# Functionality

`Manopt.jl` provides a comprehensive framework for optimization on Riemannian manifolds and a variety of algorithms using this framework.
The framework includes a generic way to specify a step size and a stopping criterion, as well as enhance the algorithm with debug and recording capabilities.
Each of the algorithms has a high-level interface to make it easy to use the algorithms directly.

An optimization task in `Manopt.jl` consists of a `Problem p` and `Options o`.
The `Problem` consists of all static information, like the cost function and a potential gradient of the optimization task. The `Options` specify the type of algorithm and the settings and data required to run the algorithm. For example, by default most options specify that the exponential map, which generalizes the notion of addition to the manifold, should be used and the algorithm steps are performed following an acceleration-free curve on the manifold. This might not be known in closed form for some manifolds, e.g. the [`Spectrahedron`](https://juliamanifolds.github.io/Manifolds.jl/v0.7/) does not have -- to the best of the author's knowledge -- a closed-form expression for the exponential map; hence more general arbitrary *retractions* can be specified for this instead.
Retractions are first-order approximations for the exponential map. They provide an alternative to the acceleration-free form, if no closed form solution is known. Otherwise, a retraction might also be chosen, when their evaluation is computationally cheaper than to use the exponential map, especially if their approximation error can be stated; see e.g. @BendokatZimmermann:2021.

Similarly, tangent vectors at different points are identified by a vector transport, which by default is the parallel transport.
By always providing a default, a user can start immediately, without thinking about these details. They can then modify these settings to improve speed or accuracy by specifying other retractions or vector transport to their needs.

The main methods to implement for a user-defined solver are `initialize_solver!(p,o)`, which fills the data in the options with an initial state, and `step_solver!(p,o,i)`, which performs the $i$th iteration.

Using a decorator pattern, `Options` can be encapsulated in `DebugOptions` and `RecordOptions`, which print and record arbitrary data stored within `Options`, respectively. This enables to investigate how the optimization is performed in detail and use the algorithms from within this package also for numerical analysis.

In the current version 0.3.17 of `Manopt.jl` the following algorithms are available:

* Alternating Gradient Descent ([`alternating_gradient_descent`](https://manoptjl.org/v0.3/solvers/alternating_gradient_descent.html))
* Chambolle-Pock ([`ChambollePock`](https://manoptjl.org/v0.3/solvers/ChambollePock.html)) [@BergmannHerzogSilvaLouzeiroTenbrinckVidalNunez:2021:1]
* Conjugate Gradient Descent ([`conjugate_gradient_descent`](https://manoptjl.org/v0.3/solvers/conjugate_gradient_descent.html)), which includes eight direction update rules using the `coefficient` keyword:
  [`SteepestDirectionUpdateRule`](https://manoptjl.org/v0.3/solvers/conjugate_gradient_descent.html#Manopt.SteepestDirectionUpdateRule),   [`ConjugateDescentCoefficient`](https://manoptjl.org/v0.3/solvers/conjugate_gradient_descent.html#Manopt.ConjugateDescentCoefficient). [`DaiYuanCoefficient`](https://manoptjl.org/v0.3/solvers/conjugate_gradient_descent.html#Manopt.DaiYuanCoefficient), [`FletcherReevesCoefficient`](https://manoptjl.org/v0.3/solvers/conjugate_gradient_descent.html#Manopt.FletcherReevesCoefficient), [`HagerZhangCoefficient`](https://manoptjl.org/v0.3/solvers/conjugate_gradient_descent.html#Manopt.HagerZhangCoefficient), [`HeestenesStiefelCoefficient`](https://manoptjl.org/v0.3/solvers/conjugate_gradient_descent.html#Manopt.HeestenesStiefelCoefficient), [`LiuStoreyCoefficient`](https://manoptjl.org/v0.3/solvers/conjugate_gradient_descent.html#Manopt.LiuStoreyCoefficient), and [`PolakRibiereCoefficient`](https://manoptjl.org/v0.3/solvers/conjugate_gradient_descent.html#Manopt.PolakRibiereCoefficient)
* Cyclic Proximal Point ([`cyclic_proximal_point`](https://manoptjl.org/v0.3/solvers/cyclic_proximal_point.html)) [@Bacak:2014:1]
* (parallel) Douglas–Rachford ([`DouglasRachford`](https://manoptjl.org/v0.3/solvers/DouglasRachford.html)) [@BergmannPerschSteidl:2016:1]
* Gradient Descent ([`gradient_descent`](https://manoptjl.org/v0.3/solvers/gradient_descent.html)), including direction update rules ([`IdentityUpdateRule`](https://manoptjl.org/v0.3/solvers/gradient_descent.html#Manopt.IdentityUpdateRule) for the classical gradient descent) to perform [`MomentumGradient`](https://manoptjl.org/v0.3/solvers/gradient_descent.html#Manopt.MomentumGradient), [`AverageGradient`](https://manoptjl.org/v0.3/solvers/gradient_descent.html#Manopt.AverageGradient), and [`Nesterov`](https://manoptjl.org/v0.3/solvers/gradient_descent.html#Manopt.Nesterov) types
* Nelder-Mead ([`NelderMead`](https://manoptjl.org/v0.3/solvers/NelderMead.html))
* Particle-Swarm Optimization ([`particle_swarm`](https://manoptjl.org/v0.3/solvers/particle_swarm.html)) [@BorckmansIshtevaAbsil2010]
* Quasi-Newton ([`quasi_Newton`](https://manoptjl.org/v0.3/solvers/quasi_Newton.html)), with [`BFGS`](https://manoptjl.org/v0.3/solvers/quasi_Newton.html#Manopt.BFGS), [`DFP`](https://manoptjl.org/v0.3/solvers/quasi_Newton.html#Manopt.DFP), [`Broyden`](https://manoptjl.org/v0.3/solvers/quasi_Newton.html#Manopt.Broyden) and a symmetric rank 1 ([`SR1`](https://manoptjl.org/v0.3/solvers/quasi_Newton.html#Manopt.SR1)) update, their inverse updates as well as a limited memory variant of (inverse) BFGS  (using the `memory` keyword) [@HuangGallivanAbsil:2015:1]
* Stochastic Gradient Descent ([`stochastic_gradient_descent`](https://manoptjl.org/v0.3/solvers/stochastic_gradient_descent.html))
* Subgradient Method ([`subgradient_method`](https://manoptjl.org/v0.3/solvers/subgradient.html))
* Trust Regions ([`trust_regions`](https://manoptjl.org/v0.3/solvers/trust_regions.html)), with inner Steihaug-Toint ([`truncated_conjugate_gradient_descent`](https://manoptjl.org/v0.3/solvers/truncated_conjugate_gradient_descent.html)) solver [@AbsilBakerGallivan2006]

# Example

`Manopt.jl` is registered in the general Julia registry and can hence be installed typing `]add Manopt` in the Julia REPL.
Given the [`Sphere`](https://juliamanifolds.github.io/Manifolds.jl/v0.7/manifolds/sphere.html) from `Manifolds.jl` and a set of unit vectors $p_1,...,p_N\in\mathbb R^3$, where $N$ is the number of data points,
we can compute the generalization of the mean, called the Riemannian Center of Mass [@Karcher:1977:1], defined as the minimizer of the squared distances to the given data – a property that the mean in vector spaces fulfills:

$$ \operatorname*{arg\,min}_{x\in\mathcal M}\quad \displaystyle\sum_{k=1}^Nd_{\mathcal M}(x, p_k)^2, $$

where $d_{\mathcal M}$ denotes the length of a shortest geodesic connecting the points specified by its two arguments;
this is called the Riemannian distance. For the sphere this [`distance`](https://juliamanifolds.github.io/Manifolds.jl/v0.7/manifolds/sphere.html#ManifoldsBase.distance-Tuple{AbstractSphere,%20Any,%20Any}) is given by the length of the shorter great arc connecting the two points.

```julia
using Manopt, Manifolds, LinearAlgebra, Random
Random.seed!(42)
M = Sphere(2)
n = 40
p = 1/sqrt(3) .* ones(3)
B = DefaultOrthonormalBasis()
pts = [ exp(M, p, get_vector(M, p, 0.425*randn(2), B)) for _ in 1:n ]

F(M, y) = sum(1/(2*n) * distance.(Ref(M), pts, Ref(y)).^2)
gradF(M, y) = sum(1/n * grad_distance.(Ref(M), pts, Ref(y)))

x_mean = gradient_descent(M, F, gradF, pts[1])
```

The resulting `x_mean` minimizes the (Riemannian) distances squared, but is especially a point of unit norm.
This should be compared to `mean(pts)`, which computes the mean in the embedding of the sphere, $\mathbb R^3$, and yields a point “inside” the sphere,
since its norm is approximately `0.858`. But even projecting this back onto the sphere yields a point that does not fulfill the property of minimizing the squared distances.

In the following figure the data `pts` (teal) and the resulting mean (orange) as well as the projected Euclidean mean (small, cyan) are shown.

![40 random points `pts` and the result from the gradient descent to compute the `x_mean` (orange) compared to a projection of their (Eucliean) mean onto the sphere (cyan).](src/img/MeanIllustr.png)

In order to print the current iteration number, change and cost every iteration as well as the stopping reason, you can provide a `debug` keyword with the corresponding symbols interleaved with strings. The Symbol `:Stop` indicates that the reason for stopping reason should be printed at the end. The last integer in this array specifies that debugging information should be printed only every $i$th iteration.
While `:x` could be used to also print the current iterate, this usually takes up too much space.

It might be more reasonable to *record* these data instead.
The `record` keyword can be used for this, for example to record the current iterate `:x`,  the `:Change` from one iterate to the next and the current function value or `:Cost`.
To access the recorded values, set `return_options` to `true`, to obtain not only the resulting value as in the example before, but the whole `Options` structure.
Then the values can be accessed using the `get_record` function.
Just calling `get_record` returns an array of tuples, where each tuple stores the values of one iteration.
To obtain an array of values for one recorded value,
use the access per symbol, i.e. from the `Iteration`s we want to access the recorded iterates `:x` as follows:

```julia
o = gradient_descent(M, F, gradF, pts[1],
    debug=[:Iteration, " | ", :Change, " | ", :Cost, "\n", :Stop],
    record=[:x, :Change, :Cost],
    return_options=true
)
x_mean_2 = get_solver_result(o) # the solver result
all_values = get_record(o) # a tuple of recorded data per iteration
iterates = get_record(o, :Iteration, :x) # iterates recorded per iteration
```

The debugging output of this example looks as follows:

```julia-repl
Initial |  | F(x): 0.20638171781316278
# 1 | Last Change: 0.22025631624261213 | F(x): 0.18071614247165613
# 2 | Last Change: 0.014654955252636971 | F(x): 0.1805990319857418
# 3 | Last Change: 0.0013696682667046617 | F(x): 0.18059800144857607
# 4 | Last Change: 0.00013562945413135856 | F(x): 0.1805979913344784
# 5 | Last Change: 1.3519139571830234e-5 | F(x): 0.1805979912339798
# 6 | Last Change: 1.348534506171897e-6 | F(x): 0.18059799123297982
# 7 | Last Change: 1.3493575361575816e-7 | F(x): 0.1805979912329699
# 8 | Last Change: 2.580956827951785e-8 | F(x): 0.18059799123296988
# 9 | Last Change: 2.9802322387695312e-8 | F(x): 0.18059799123296993
The algorithm reached approximately critical point after 9 iterations;
    the gradient norm (1.3387605239861564e-9) is less than 1.0e-8.
```

For more details on more algorithms to compute the mean and other statistical functions on manifolds like the median
see [https://juliamanifolds.github.io/Manifolds.jl/v0.7/features/statistics.html](https://juliamanifolds.github.io/Manifolds.jl/v0.7/features/statistics.html).

# Related research and software

The two projects that are most similar to `Manopt.jl` are [`Manopt`](https://manopt.org) [@manopt] in Matlab and [`pymanopt`](https://pymanopt.org) [@pymanopt] in Python.
Similarly [`ROPTLIB`](https://www.math.fsu.edu/~whuang2/Indices/index_ROPTLIB.html) [@HuangAbsilGallivanHand:2018:1] is a package for optimization on Manifolds in C++.
While all three packages cover some algorithms, most are less flexible, for example in stating the stopping criterion, which is fixed to mainly the maximal number of iterations or a small gradient. Most prominently, `Manopt.jl` is the first package that also covers methods for high-performance and high-dimensional nonsmooth optimization on manifolds.

The Riemannian Chambolle-Pock algorithm presented in @BergmannHerzogSilvaLouzeiroTenbrinckVidalNunez:2021:1 was developed using `Manopt.jl`. Based on this theory and algorithm, a higher-order algorithm was introduced in @DiepeveenLellmann:2021:1. Optimized examples from @BergmannGousenbourger:2018:2 performing data interpolation and approximation with manifold-valued Bézier curves are also included in `Manopt.jl`.

# References
