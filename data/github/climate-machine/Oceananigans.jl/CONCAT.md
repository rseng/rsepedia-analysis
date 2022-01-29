<!-- Title -->
<h1 align="center">
  Oceananigans.jl
</h1>

<!-- description -->
<p align="center">
  <strong>🌊 Fast and friendly ocean-flavored Julia software for simulating incompressible fluid dynamics in Cartesian and spherical shell domains on CPUs and GPUs. https://clima.github.io/OceananigansDocumentation/stable</strong>
</p>

<!-- Information badges -->
<p align="center">
  <a href="https://www.repostatus.org/#active">
    <img alt="Repo status" src="https://www.repostatus.org/badges/latest/active.svg?style=flat-square" />
  </a>
  <a href="https://mit-license.org">
    <img alt="MIT license" src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square">
  </a>
  <a href="https://github.com/CliMA/Oceananigans.jl/discussions">
    <img alt="Ask us anything" src="https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg?style=flat-square">
  </a>
  <a href="https://github.com/SciML/ColPrac">
    <img alt="ColPrac: Contributor's Guide on Collaborative Practices for Community Packages" src="https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet?style=flat-square">
  </a>
  <a href="https://doi.org/10.21105/joss.02018">
    <img alt="JOSS" src="https://joss.theoj.org/papers/10.21105/joss.02018/status.svg">
  </a>
</p>

<!-- Version and documentation badges -->
<p align="center">
  <a href="https://github.com/CliMA/Oceananigans.jl/releases">
    <img alt="GitHub tag (latest SemVer pre-release)" src="https://img.shields.io/github/v/tag/CliMA/Oceananigans.jl?include_prereleases&label=latest%20version&logo=github&sort=semver&style=flat-square">
  </a>
  <a href="https://clima.github.io/OceananigansDocumentation/stable">
    <img alt="Stable documentation" src="https://img.shields.io/badge/documentation-stable%20release-blue?style=flat-square">
  </a>
  <a href="https://clima.github.io/OceananigansDocumentation/dev">
    <img alt="Development documentation" src="https://img.shields.io/badge/documentation-in%20development-orange?style=flat-square">
  </a>
</p>

<!-- CI/CD badges -->
<p align="center">
  <a href="https://buildkite.com/clima/oceananigans">
    <img alt="Buildkite CPU+GPU build status" src="https://img.shields.io/buildkite/4d921fc17b95341ea5477fb62df0e6d9364b61b154e050a123/main?logo=buildkite&label=Buildkite%20CPU%2BGPU&style=flat-square">
  </a>
  <a href="https://hub.docker.com/r/aliramadhan/oceananigans">
    <img alt="Docker build status" src="https://img.shields.io/docker/cloud/build/aliramadhan/oceananigans?label=Docker&logo=docker&logoColor=white&style=flat-square">
  </a>
</p>

Oceananigans.jl is a fast and friendly fluid flow solver written in Julia that can be run in 1-3 dimensions on CPUs and GPUs. It can simulate the incompressible Boussinesq equations, the shallow water equations, or the hydrostatic Boussinesq equations with a free surface. Oceananigans.jl comes with user-friendly features for simulating rotating stratified fluids including user-defined boundary conditions and forcing functions, arbitrary tracers, large eddy simulation turbulence closures, high-order advection schemes, immersed boundaries, Lagrangian particle tracking, and more!

We strive for a user interface that makes Oceananigans.jl as friendly and intuitive to use as possible, allowing users to focus on the science. Internally, we have attempted to write the underlying algorithm so that the code runs as fast as possible for the configuration chosen by the user --- from simple two-dimensional setups to complex three-dimensional simulations --- and so that as much code as possible is shared between the different architectures, models, and grids.

## Contents

* [Installation instructions](#installation-instructions)
* [Running your first model](#running-your-first-model)
* [Getting help](#getting-help)
* [Contributing](#contributing)
* [Movies](#movies)
* [Performance benchmarks](#performance-benchmarks)

## Installation instructions

You can install the latest version of Oceananigans using the built-in package manager (accessed by pressing `]` in the Julia command prompt) to add the package and instantiate/build all depdendencies

```julia
julia>]
(v1.6) pkg> add Oceananigans
(v1.6) pkg> instantiate
```

We recommend installing Oceananigans with the built-in Julia package manager, because this installs a stable, tagged release. Oceananigans.jl can be updated to the latest tagged release from the package manager by typing

```julia
(v1.6) pkg> update Oceananigans
```

At this time, updating should be done with care, as Oceananigans is under rapid development and breaking changes to the user API occur often. But if anything does happen, please open an issue!

**Note**: The latest version of Oceananigans requires at least Julia v1.6 to run. Installing Oceananigans with an older version of Julia will install an older version of Oceananigans (the latest version compatible with your version of Julia).

## Running your first model

Let's initialize a 3D horizontally periodic model with 100×100×50 grid points on a 2×2×1 km domain and simulate it for 1 hour using a constant time step of 60 seconds.

```julia
using Oceananigans
grid = RectilinearGrid(size=(100, 100, 50), extent=(2π, 2π, 1))
model = NonhydrostaticModel(grid=grid)
simulation = Simulation(model, Δt=60, stop_time=3600)
run!(simulation)
```

You just simulated what might have been a 3D patch of ocean, it's that easy! It was a still lifeless ocean so nothing interesting happened but now you can add interesting dynamics and visualize the output.

### More interesting example

Let's add something to make the dynamics a bit more interesting. We can add a hot bubble in the middle of the domain and watch it rise to the surface. This example also shows how to set an initial condition.

```julia
using Oceananigans

N = Nx = Ny = Nz = 128   # Number of grid points in each dimension.
L = Lx = Ly = Lz = 2000  # Length of each dimension.
topology = (Periodic, Periodic, Bounded)

model = NonhydrostaticModel(
            grid = RectilinearGrid(CPU(); topology=topology, size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz)),
         tracers = (:T, :S),
        buoyancy = SeawaterBuoyancy(),
         closure = IsotropicDiffusivity(ν=4e-2, κ=4e-2)
)

# Set a temperature perturbation with a Gaussian profile located at the center.
# This will create a buoyant thermal bubble that will rise with time.
x₀, z₀ = Lx/2, Lz/2
T₀(x, y, z) = 20 + 0.01 * exp(-100 * ((x - x₀)^2 + (z - z₀)^2) / (Lx^2 + Lz^2))
set!(model, T=T₀)

simulation = Simulation(model, Δt=10, stop_iteration=5000)
run!(simulation)
```

By changing `architecture = CPU()` to `architecture = GPU()` the example will run on a CUDA-enabled Nvidia GPU!

You can see some movies from GPU simulations below along with CPU and GPU [performance benchmarks](https://github.com/clima/Oceananigans.jl#performance-benchmarks).

## Getting help

If you are interested in using Oceananigans.jl or are trying to figure out how to use it, please feel free to ask us questions and get in touch! If you're trying to set up a model then check out the examples and model setup documentation. Check out the [examples](https://github.com/clima/Oceananigans.jl/tree/main/examples) and please feel free to [start a discussion](https://github.com/CliMA/Oceananigans.jl/discussions) if you have any questions, comments, suggestions, etc! There is also an #oceananigans channel on the [Julia Slack](https://julialang.org/slack/).

## Citing

If you use Oceananigans.jl as part of your research, teaching, or other activities, we would be grateful if you could cite our work and mention Oceananigans.jl by name.

```bibtex
@article{OceananigansJOSS,
  doi = {10.21105/joss.02018},
  url = {https://doi.org/10.21105/joss.02018},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {53},
  pages = {2018},
  author = {Ali Ramadhan and Gregory LeClaire Wagner and Chris Hill and Jean-Michel Campin and Valentin Churavy and Tim Besard and Andre Souza and Alan Edelman and Raffaele Ferrari and John Marshall},
  title = {Oceananigans.jl: Fast and friendly geophysical fluid dynamics on GPUs},
  journal = {Journal of Open Source Software}
}
```

We also maintain a [list of publication using Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/#Papers-and-preprints-using-Oceananigans.jl). If you have work using Oceananigans.jl that you would like to have listed there, please open a pull request to add it or let us know!

## Contributing

If you're interested in contributing to the development of Oceananigans we want your help no matter how big or small a contribution you make! It's always great to have new people look at the code with fresh eyes: you will see errors that other developers have missed.

Let us know by [opening an issue](https://github.com/clima/Oceananigans.jl/issues/new) if you'd like to work on a new feature or if you're new to open-source and want to find a cool little project or issue to work on that fits your interests! We're more than happy to help along the way.

For more information, check out our [contributor's guide](https://github.com/clima/Oceananigans.jl/blob/main/CONTRIBUTING.md).

## Movies

### [Deep convection](https://www.youtube.com/watch?v=kpUrxnKKMjI)

[![Watch deep convection in action](https://raw.githubusercontent.com/ali-ramadhan/ali-ramadhan.Github.io/master/img/surface_temp_3d_00130_halfsize.png)](https://www.youtube.com/watch?v=kpUrxnKKMjI)

### [Free convection](https://www.youtube.com/watch?v=yq4op9h3xcU)

[![Watch free convection in action](https://raw.githubusercontent.com/ali-ramadhan/ali-ramadhan.Github.io/master/img/free_convection_0956.png)](https://www.youtube.com/watch?v=yq4op9h3xcU)

### [Winds blowing over the ocean](https://www.youtube.com/watch?v=IRncfbvuiy8)

[![Watch winds blowing over the ocean](https://raw.githubusercontent.com/ali-ramadhan/ali-ramadhan.Github.io/master/img/wind_stress_0400.png)](https://www.youtube.com/watch?v=IRncfbvuiy8)

### [Free convection with wind stress](https://www.youtube.com/watch?v=ob6OMQgPfI4)

[![Watch free convection with wind stress in action](https://raw.githubusercontent.com/ali-ramadhan/ali-ramadhan.Github.io/master/img/wind_stress_unstable_7500.png)](https://www.youtube.com/watch?v=ob6OMQgPfI4)

## Performance benchmarks

We've performed some preliminary performance benchmarks (see the [performance benchmarks](https://clima.github.io/OceananigansDocumentation/stable/appendix/benchmarks/) section of the documentation) by initializing models of various sizes and measuring the wall clock time taken per model iteration (or time step).

This is not really a fair comparison as we haven't parallelized across all the CPU's cores so we will revisit these benchmarks once Oceananigans.jl can run on multiple CPUs and GPUs.

To make full use of or fully saturate the computing power of a GPU such as an Nvidia Tesla V100 or
a Titan V, the model should have around ~10 million grid points or more.

Sometimes counter-intuitively running with `Float32` is slower than `Float64`. This is likely due
to type mismatches causing slowdowns as floats have to be converted between 32-bit and 64-bit, an
issue that needs to be addressed meticulously. Due to other bottlenecks such as memory accesses and
GPU register pressure, `Float32` models may not provide much of a speedup so the main benefit becomes
lower memory costs (by around a factor of 2).

![Performance benchmark plots](https://user-images.githubusercontent.com/20099589/89906791-d2c85b00-dbb9-11ea-969a-4b8db2c31680.png)
# Contributors Guide

Thank you for considering contributing to Oceananigans! 

Feel free to ask us questions and chat with us at any time about any topic at all
by 

* [Opening a GitHub issue](https://github.com/CliMA/Oceananigans.jl/issues/new/choose)
 
* [Creating a GitHub discussion](https://github.com/CliMA/Oceananigans.jl/discussions/new)

* Sending a message to the [#oceananigans channel](https://julialang.slack.com/archives/C01D24C0CAH) on [Julia Slack](https://julialang.org/slack/).

We follow the [ColPrac guide](https://github.com/SciML/ColPrac) for collaborative
practices. We ask that new contributors read that guide before submitting a pull request.

## Creating issues

The simplest way to contribute to Oceananigans is to create or comment on issues and discussions.

The most useful bug reports:

* Provide an explicit code snippet --- not just a link --- that reproduces the bug in the latest tagged version of Oceananigans. This is sometimes called the ["minimal working example"](https://en.wikipedia.org/wiki/Minimal_working_example). Reducing bug-producing code to a minimal example can dramatically decrease the time it takes to resolve an issue.

* Paste the _entire_ error received when running the code snippet, even if it's unbelievably long.

* Use triple backticks (```` ``` ````) to enclose code snippets, and other [markdown formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) to make your issue easy and quick to read.

* Report the Oceananigans version, Julia version, machine (especially if using a GPU) and any other possibly useful details of the computational environment in which the bug was created.

Discussions are recommended for asking questions about (for example) the user interface, implementation details, science, and life in general.

## But I want to _code_!

* New users help write Oceananigans code and documentation by [forking the Oceananigans repository](https://docs.github.com/en/github/collaborating-with-pull-requests/working-with-forks), [using git](https://guides.github.com/introduction/git-handbook/) to edit code and docs, and then creating a [pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork). Pull requests are reviewed by Oceananigans collaborators.

* A pull request can be merged once it is reviewed and approved by collaborators. If the pull request author has write access, they have the reponsibility of merging their pull request. Otherwise, Oceananigans.jl collabators will execute the merge with permission from the pull request author.

* Note: for small or minor changes (such as fixing a typo in documentation), the [GitHub editor](https://docs.github.com/en/github/managing-files-in-a-repository/managing-files-on-github/editing-files-in-your-repository) is super useful for forking and opening a pull request with a single click.

* Write your code with love and care. In particular, conform to existing Oceananigans style and formatting conventions. For example, we love verbose and explicit variable names, use `TitleCase` for types, `snake_case` for objects, and always.put.spaces.after.commas. For formatting decisions we loosely follow the [YASGuide](https://github.com/jrevels/YASGuide). It's worth few extra minutes of our time to leave future generations with well-written, readable code.

## What is a "collaborator" and how can I become one?

* Collaborators have permissions to review pull requests and  status allows a contributor to review pull requests in addition to opening them. Collaborators can also create branches in the main Oceananigans repository.

* We ask that new contributors try their hand at forking Oceananigans, and opening and merging a pull request before requesting collaborator status.

## What's a good way to start developing Oceananigans?

* Tackle an existing issue. We keep a list of [good first issues](https://github.com/CLiMA/Oceananigans.jl/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22)
  that are self-contained and suitable for a newcomer to try and work on.

* Try to run Oceananigans and play around with it to simulate your favorite
  fluids and ocean physics. If you run into any bugs/problems or find it difficult
  to use or understand, please open an issue!

* Write up an example or tutorial on how to do something useful with
  Oceananigans, like how to set up a new physical configuration.

* Improve documentation or comments if you found something hard to use.

* Implement a new feature if you need it to use Oceananigans.

If you're interested in working on something, let us know by commenting on
existing issues or by opening a new issue. This is to make sure no one else
is working on the same issue and so we can help and guide you in case there
is anything you need to know beforehand.

We also hang out on the #oceananigans channel on Julia Slack, which is a great
place to discuss anything Oceananigans-related, especially contributions! To
join the Julia Slack, go to [https://julialang.org/slack/](https://julialang.org/slack/).
---
title: 'Oceananigans.jl: Fast and friendly geophysical fluid dynamics on GPUs'
tags:
  - fluid
  - ocean
  - climate
  - Julia
  - gpu
authors:
  - name: Ali Ramadhan
    orcid: 0000-0003-1102-1520
    affiliation: 1
  - name: Gregory LeClaire Wagner
    orcid: 0000-0001-5317-2445
    affiliation: 1
  - name: Chris Hill
    affiliation: 1
  - name: Jean-Michel Campin
    affiliation: 1
  - name: Valentin Churavy
    affiliation: 1
  - name: Tim Besard
    affiliation: 2
  - name: Andre Souza
    affiliation: 1
  - name: Alan Edelman
    affiliation: 1
  - name: Raffaele Ferrari
    affiliation: 1
  - name: John Marshall
    affiliation: 1
affiliations:
 - name: Massachusetts Institute of Technology
   index: 1
 - name: Julia Computing, Inc.
   index: 2
date: 11 August 2020
bibliography: paper.bib
---

# Summary

``Oceananigans.jl`` is a fast and friendly software package for the numerical
simulation of incompressible, stratified, rotating fluid flows on CPUs and GPUs.
``Oceananigans.jl`` is fast and flexible enough for research yet simple enough
for students and first-time programmers. ``Oceananigans.jl`` is being developed
as part of the Climate Modeling Alliance project for the simulation of
small-scale ocean physics at high-resolution that affect the evolution of
Earth’s climate.

``Oceananigans.jl`` is designed for high-resolution simulations in idealized
geometries and supports direct numerical simulation, large eddy simulation,
arbitrary numbers of active and passive tracers, and linear and nonlinear
equations of state for seawater. Under the hood, ``Oceananigans.jl`` employs a
finite volume algorithm similar to that used by the Massachusetts Institute of
Technology general circulation model [@Marshall1997].

![Fig. 1](free_convection_and_baroclinic_instability.png)
Fig. 1: (Left) Large eddy simulation of small-scale oceanic boundary layer
turbulence forced by a surface cooling in a horizontally periodic domain using
$256^3$ grid points. The upper layer is well-mixed by turbulent convection and
bounded below by a strong buoyancy interface. (Right) Simulation of
instability of a horizontal density gradient in a rotating channel using
$256\times512\times128$ grid points. A similar process called baroclinic
instability acting on basin-scale temperature gradients fills the oceans with
eddies that stir carbon and heat. Plots made with `matplotlib` [@Hunter2007]
and `cmocean` [@Thyng2016].

``Oceananigans.jl`` leverages the Julia programming language [@Bezanson2017] to
implement high-level, low-cost abstractions, a friendly user interface, and a
high-performance model in one language and a common code base for execution on
the CPU or GPU with Julia’s native GPU compiler [@Besard2019]. Because Julia is
a high-level language, development is streamlined and users can flexibly specify
model configurations, set up arbitrary diagnostics and output, extend the code
base, and implement new features. Configuring a model with `architecture=CPU()`
or `architecture=GPU()` will execute the model on the CPU or GPU. By pinning a
simulation script against a specific version of Oceananigans, simulation results
are reproducible up to hardware differences.

Performance benchmarks show significant speedups when running on a GPU. Large
simulations on an Nvidia Tesla V100 GPU require ~1 nanosecond per grid point per
iteration. GPU simulations are therefore roughly 3x more cost-effective
than CPU simulations on cloud computing platforms such as Google Cloud. A GPU
with 32 GB of memory can time-step models with ~150 million grid points assuming
five fields are being evolved; for example, three velocity components and
tracers for temperature and salinity. These performance gains permit the
long-time integration of realistic simulations, such as large eddy simulation of
oceanic boundary layer turbulence over a seasonal cycle or the generation of
training data for turbulence parameterizations in Earth system models.

``Oceananigans.jl`` is continuously tested on CPUs and GPUs with unit tests,
integration tests, analytic solutions to the incompressible Navier-Stokes
equations, convergence tests, and verification experiments against published
scientific results. Future development plans include support for distributed
parallelism with CUDA-aware MPI as well as topography.

Ocean models that are similar to ``Oceananigans.jl`` include MITgcm
[@Marshall1997] and MOM6 [@Adcroft2019], both written in Fortran. However,
``Oceananigans.jl`` features a more efficient non-hydrostatic pressure solver
than MITgcm (and MOM6 is strictly hydrostatic). PALM [@Maronga2020] is Fortran
software for large eddy simulation of atmospheric and oceanic boundary layers
with complex boundaries on parallel CPU and GPU architectures. ``Oceananigans.jl``
is distinguished by its use of Julia which allows for a script-based interface as
opposed to a configuration-file-based interface used by MITgcm, MOM6, and PALM.
Dedalus [@Burns2020] is Python software with an intuitive script-based interface
that solves general partial differential equations, including the incompressible
Navier-Stokes equations, with spectral methods.

# Acknowledgements

Our work is supported by the generosity of Eric and Wendy Schmidt by
recommendation of the Schmidt Futures program, and by the National Science
Foundation under grant AGS-6939393.

# References
# Convergence Tests

This directory contains scripts and modules for testing the numerical
convergence of `Oceananigans` time stepping algorithms and spatial discretiation.

To instantiate the convergence test environment, run

```
julia -e 'using Pkg; Pkg.activate(pwd()); Pkg.instantiate(); Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..", "..")))'
```

## Time stepping convergence tests

```
julia --project point_exponential_decay.jl
```

produces `figs/point_exponential_decay_time_stepper_convergence.png`.

## One-dimensional advection-diffusion tests

### Advection and diffusion of a cosine

```
julia --project one_dimensional_cosine_advection_diffusion.jl
```

produces

* `figs/cosine_advection_diffusion_solutions.png`
* `figs/cosine_advection_diffusion_error_convergence.png`

### Advection and diffusion of a Gaussian

```
julia --project one_dimensional_gaussian_advection_diffusion.jl
```

produces

* `figs/gaussian_advection_diffusion_solutions.png`
* `figs/gaussian_advection_diffusion_error_convergence.png`

## Two-dimensional diffusion

```
julia --project two_dimensional_diffusion.jl
```

produces `figs/two_dimensional_diffusion_convergence.png`.

## Two-dimensional Taylor-Green vortex

```
julia --project run_taylor_green.jl
```

and then

```
julia --project analyze_taylor_green.jl
```

produces `figs/taylor_green_convergence.png`.

## Two-dimensional forced flow with free-slip boundary conditions

```
julia --project run_forced_free_slip.jl
```

followed by

```
julia --project analyze_forced_free_slip.jl
```

produces `figs/forced_free_slip_convergence.png`.

## Two-dimensional forced flow with fixed-slip boundary conditions

```
julia --project run_forced_fixed_slip.jl
```

followed by

```
julia --project analyze_forced_fixed_slip.jl
```

produces `figs/forced_fixed_slip_convergence.png`.
## Differentiation and interpolation operators

The geometry of the staggerd grid used by Oceananigans.jl (sometimes called the "C grid")
is (in one dimension) shown below
```
face   cell   face   cell   face

        i-1            i
         ↓             ↓
  |      ×      |      ×      |
  ↑             ↑             ↑
 i-1            i            i+1
```
Difference operators are denoted by a `δ` (`\delta`). Calculating the difference
of a cell-centered quantity `c` at cell `i` returns the difference at face `i`
```
δcᵢ = cᵢ - cᵢ₋₁
```
and so this operation, if applied along the x-dimension, is denoted by `δxᶠᵃᵃ`.

The difference of a face-centered quantity `u` at face `i` returns the difference at cell `i`
```
δuᵢ = uᵢ₊₁ - uᵢ
```
and is thus denoted `δxᶜᵃᵃ` when applied along the x-dimension.

The three characters at the end of the function name, `faa` for example, indicates that the
output lies on the cell faces in the x-dimension but remains at their original positions in 
the y- and z-dimensions. Thus we further identify this operator by the superscript `ᶠᵃᵃ`, where
the `a` stands for "any" as the location is unchanged by the operator and is determined by
the input.

As a result, the interpolation of a quantity `c` from a cell `i` to face `i` (which is denoted
"`ℑxᶠᵃᵃ`" in the code below) is
```
ℑxᶠᵃᵃ(c)ᵢ = (cᵢ + cᵢ₋₁) / 2
```
Conversely, the interpolation of a quantity `u` from a face `i` to cell `i` is given by
```
ℑxᶠᵃᵃ(u)ᵢ = (uᵢ₊₁ + uᵢ) / 2
```
The `ℑ` (`\Im`) symbol indicates that an interpolation is being performed. For example, `ℑx`
indicates that the interpolation is performed along the x-dimension.
# Oceananigans.jl performance benchmarks

This directory contains scripts and modules for benchmarking various features of Oceananigans.

To instantiate the benchmarks environment, run

```
julia -e 'using Pkg; Pkg.activate(pwd()); Pkg.instantiate(); Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))'
```

Once the environment has been instantiated, benchmarks can be run via, e.g.

```
julia --project benchmark_nonhydrostatic_model.jl
```

Most scripts benchmark one feature (e.g. advection schemes, arbitrary tracers). If your machine contains a CUDA-compatible GPU, benchmarks will also run on the GPU. Tables with benchmark results will be printed (and each table will also be saved to an HTML file).

## Distributed benchmarks

Run distributed benchmarks by running the launcher scripts for either the shallow water model: `distributed_shallow_water_model.jl` or the nonhydrostatic model: `distributed_nonhydrostatic_model.jl`. Change settings within the scripts to toggle between strong or weak scaling and threaded or MPI architecture. The single and serial scripts executed by the launcher scripts can also be executed manually from the command line with the appropriate arguments.

## Measuring performance regression

Running the `benchmark_regression.jl` script will run the nonhydrostatic model tests on the current branch and on the main branch for comparison. This is useful to test whether the current branch slows down the code or introduces any performance regression.

# Oceananigans.jl

*🌊 Fast and friendly fluid dynamics on CPUs and GPUs.*

Oceananigans.jl is a fast and friendly fluid flow solver written in Julia that can be run in 1-3 dimensions on CPUs
and GPUs. It can simulate the incompressible Boussinesq equations, the shallow water equations, or the hydrostatic
Boussinesq equations with a free surface. Oceananigans.jl comes with user-friendly features for simulating rotating
stratified fluids including user-defined boundary conditions and forcing functions, arbitrary tracers, large eddy
simulation turbulence closures, high-order advection schemes, immersed boundaries, Lagrangian particle tracking, and
more!

We strive for a user interface that makes Oceananigans.jl as friendly and intuitive to use as possible,
allowing users to focus on the science. Internally, we have attempted to write the underlying algorithm
so that the code runs as fast as possible for the configuration chosen by the user --- from simple
two-dimensional setups to complex three-dimensional simulations --- and so that as much code
as possible is shared between the different architectures, models, and grids.

## Getting help

If you are interested in using Oceananigans.jl or are trying to figure out how to use it, please feel free to ask us
questions and get in touch! If you're trying to set up a model then check out the examples and model setup
documentation. Please feel free to [start a discussion](https://github.com/CliMA/Oceananigans.jl/discussions)
if you have any questions, comments, suggestions, etc! There is also an #oceananigans channel on the
[Julia Slack](https://julialang.org/slack/).

## Citing

If you use Oceananigans.jl as part of your research, teaching, or other activities, we would be grateful if you could
cite our work and mention Oceananigans.jl by name.

```bibtex
@article{OceananigansJOSS,
  doi = {10.21105/joss.02018},
  url = {https://doi.org/10.21105/joss.02018},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {53},
  pages = {2018},
  author = {Ali Ramadhan and Gregory LeClaire Wagner and Chris Hill and Jean-Michel Campin and Valentin Churavy and Tim Besard and Andre Souza and Alan Edelman and Raffaele Ferrari and John Marshall},
  title = {Oceananigans.jl: Fast and friendly geophysical fluid dynamics on GPUs},
  journal = {Journal of Open Source Software}
}
```

## Papers and preprints using Oceananigans.jl

If you have work using Oceananigans.jl that you would like to have listed here, please open a pull request to add it or let us know!

1. Wagner, G. L., Chini, G. P., Ramadhan, A., Gallet, B., & Ferrari, R. (2021). [Near-inertial waves and turbulence driven by the growth of swell](https://doi.org/10.1175/JPO-D-20-0178.1), _Journal of Physical Oceanography_, **51(5)**, 1337-1351. DOI: [10.1175/JPO-D-20-0178.1](https://doi.org/10.1175/JPO-D-20-0178.1)

1. Buffett, B. A. (2021). [Conditions for turbulent Ekman layers in precessionally driven flow](https://doi.org/10.1093/gji/ggab088), _Geophysical Journal International_, **226(1)**, 56–65. DOI: [10.1093/gji/ggab088](https://doi.org/10.1093/gji/ggab088)

1. Bhamidipati, N., Souza, A.N. & Flierl, G.R., (2020). [Turbulent mixing of a passive scalar in the ocean mixed layer](https://doi.org/10.1016/j.ocemod.2020.101615). _Ocean Modelling_, **149**, 101615. DOI: [10.1016/j.ocemod.2020.101615](https://doi.org/10.1016/j.ocemod.2020.101615)

1. Souza, A. N., Wagner, G. L., Ramadhan, A., Allen, B., Churavy, V., Schloss, J., Campin, J. M., Hill, C., Edelman, A., Marshall, J., Flierl, G., & Ferrari, R. (2020). [Uncertainty quantification of ocean parameterizations: Application to the K‐Profile‐Parameterization for penetrative convection](https://doi.org/10.1029/2020MS002108). _Journal of Advances in Modeling Earth Systems_, **12**, e2020MS002108. DOI: [10.1029/2020MS002108](https://doi.org/10.1029/2020MS002108)
# References

```@bibliography
```
# Using GPUs

A big feature of Oceananigans is being able to run on graphical processing units (GPUs)
for increased performance. Depending on your CPU and GPU combination, speedups of >150x
are possible, for example on Google Cloud where running on GPUs is more cost-effective.
See the [performance benchmarks](@ref performance_benchmarks) for more details.

See [Architecture](@ref) for instructions on setting up a model on a GPU.

Oceananigans does not yet support distributed parallelism (multi-CPU or multi-GPU).

!!! tip "Running on GPUs"
    If you are having issues with running Oceananigans on a GPU or setting things up,
    please [open an issue](https://github.com/CLiMA/Oceananigans.jl/issues/new)
    and we'll do our best to help out!

## When to use a GPU

GPUs are very useful for running large simulations. If your simulation uses over
1,000,000 grid points, you will probably benefit significantly from running your
simulation on a GPU.

GPU simulations tend to be memory-limited. That is, you'll probably fill the GPU's
memory long before the model becomes unbearably slow. High-end GPUs such as the
Nvidia Tesla V100 only come with up to 32 GB of memory. On a GPU with 16 GB of memory,
you can run a simulation (with 2 tracers) with up to ~50 million grid points.

## Getting access to GPUs

If you don't have a GPU there are a few resources you can try to acquire a GPU from.

In general, to get good performance you'll want a GPU with good 64-bit floating point
performance although Oceananigans can be used with 32-bit floats. Most recent gaming GPUs
should work but might have poor 64-bit float performance.

If you have access to any supercomputer clusters, check to see if they have any GPUs.
See also this Stack Overflow post:
[Where can I get access to GPU cluster for educational purpose?](https://scicomp.stackexchange.com/questions/8508/where-can-i-get-access-to-gpu-cluster-for-educational-purpose)

Cloud computing providers such as Google Cloud and Amazon EC2 allow you to rent GPUs per
hour. Sometimes they offer free trials or credits that can be used towards GPUs although
they seem to be getting less common.

See the [Julia on Google Colab: Free GPU-Accelerated Shareable Notebooks](https://discourse.julialang.org/t/julia-on-google-colab-free-gpu-accelerated-shareable-notebooks/15319)
post on the Julia Discourse.

[Code Ocean](https://codeocean.com/) also has
[GPU support](https://help.codeocean.com/en/articles/1053107-gpu-support) and allows you
to spin up capsules with pretty decent Tesla K80 GPUs for free (for now) if you want to
play around with them. They may not be powerful enough for huge simulations though. You'll
want to use their "Ubuntu Linux with GPU support (18.04.3)" with the ability to compile
CUDA code. Then you'll have to install Julia manually.

## I have a GPU. Now what?

Make sure you have an Nvidia GPU that is CUDA compatible:
[https://developer.nvidia.com/cuda-gpus](https://developer.nvidia.com/cuda-gpus). Most
recent GPUs should be but older GPUs and many laptop GPUs may not be.

Then download and install the CUDA Toolkit:
[https://developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads)

Once the CUDA Toolkit is installed, you might have to build Oceananigans again
```
julia>]
(v1.6) pkg> build Oceananigans
```
The [ocean wind mixing and convection](@ref gpu_example) example is a good one to test out on the GPU.
# Contributors Guide

Thank you for considering contributions to Oceananigans! We hope this guide
helps you make a contribution.

Feel free to ask us questions and chat with us at any time about any topic at all
by:

* [Opening a GitHub issue](https://github.com/CliMA/Oceananigans.jl/issues/new)

* [Creating a GitHub discussion](https://github.com/CliMA/Oceananigans.jl/discussions/new)

* Sending a message to the [#oceananigans channel](https://julialang.slack.com/archives/C01D24C0CAH) on [Julia Slack](https://julialang.org/slack/).

## Creating issues

The simplest way to contribute to Oceananigans is to create or comment on issues and discussions.

The most useful bug reports:

* Provide an explicit code snippet --- not just a link --- that reproduces the bug in the latest tagged version of Oceananigans. This is sometimes called the ["minimal working example"](https://en.wikipedia.org/wiki/Minimal_working_example). Reducing bug-producing code to a minimal example can dramatically decrease the time it takes to resolve an issue.

* Paste the _entire_ error received when running the code snippet, even if it's unbelievably long.

* Use triple backticks (e.g., ````` ```some_code; and_some_more_code;``` `````) to enclose code snippets, and other [markdown formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) to make your issue easy and quick to read.

* Report the Oceananigans version, Julia version, machine (especially if using a GPU) and any other possibly useful details of the computational environment in which the bug was created.

Discussions are recommended for asking questions about (for example) the user interface, implementation details, science, and life in general.

## But I want to _code_!

* New users help write Oceananigans code and documentation by [forking the Oceananigans repository](https://docs.github.com/en/github/collaborating-with-pull-requests/working-with-forks), [using git](https://guides.github.com/introduction/git-handbook/) to edit code and docs, and then creating a [pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork). Pull requests are reviewed by Oceananigans collaborators.

* A pull request can be merged once it is reviewed and approved by collaborators. If the pull request author has write access, they have the reponsibility of merging their pull request. Otherwise, Oceananigans.jl collabators will execute the merge with permission from the pull request author.

* Note: for small or minor changes (such as fixing a typo in documentation), the [GitHub editor](https://docs.github.com/en/github/managing-files-in-a-repository/managing-files-on-github/editing-files-in-your-repository) is super useful for forking and opening a pull request with a single click.

* Write your code with love and care. In particular, conform to existing Oceananigans style and formatting conventions. For example, we love verbose and explicit variable names, use `TitleCase` for types, `snake_case` for objects, and always.put.spaces.after.commas. For formatting decisions we loosely follow the [YASGuide](https://github.com/jrevels/YASGuide). It's worth few extra minutes of our time to leave future generations with well-written, readable code.

## What is a "collaborator" and how can I become one?

* Collaborators have permissions to review pull requests and  status allows a contributor to review pull requests in addition to opening them. Collaborators can also create branches in the main Oceananigans repository.

* We ask that new contributors try their hand at forking Oceananigans, and opening and merging a pull request before requesting collaborator status.

## What's a good way to start developing Oceananigans?

* Tackle an existing issue. We keep a list of [good first issues](https://github.com/CLiMA/Oceananigans.jl/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22)
  that are self-contained and suitable for a newcomer to try and work on.

* Try to run Oceananigans and play around with it to simulate your favorite
  fluids and ocean physics. If you run into any problems or find it difficult
  to use or understand, please open an issue!

* Write up an example or tutorial on how to do something useful with
  Oceananigans, like how to set up a new physical configuration.

* Improve documentation or comments if you found something hard to use.

* Implement a new feature if you need it to use Oceananigans.

If you're interested in working on something, let us know by commenting on existing issues or 
by opening a new issue. This is to make sure no one else is working on the same issue and so 
we can help and guide you in case there is anything you need to know beforehand.

## Ground Rules

* Each pull request should consist of a logical collection of changes. You can
  include multiple bug fixes in a single pull request, but they should be related.
  For unrelated changes, please submit multiple pull requests.

* Do not commit changes to files that are irrelevant to your feature or bugfix
  (eg: `.gitignore`).

* Be willing to accept criticism and work on improving your code; we don't want
  to break other users' code, so care must be taken not to introduce bugs. We
  discuss pull requests and keep working on them until we believe we've done a
  good job.

* Be aware that the pull request review process is not immediate, and is
  generally proportional to the size of the pull request.

## Reporting a bug

The easiest way to get involved is to report issues you encounter when using
Oceananigans or by requesting something you think is missing.

* Head over to the [issues](https://github.com/CLiMA/Oceananigans.jl/issues) page.

* Search to see if your issue already exists or has even been solved previously.

* If you indeed have a new issue or request, click the "New Issue" button.

* Please be as specific as possible. Include the version of the code you were using, as
  well as what operating system you are running. The output of Julia's `versioninfo()`
  and `] status` is helpful to include. Try your best to include a complete, ["minimal working example"](https://en.wikipedia.org/wiki/Minimal_working_example) that reproduces the issue.

## Setting up your development environment

* Install [Julia](https://julialang.org/) on your system.

* Install `git` on your system if it is not already there (install XCode command line tools on
  a Mac or `git bash` on Windows).

* Login to your GitHub account and make a fork of the
  [Oceananigans repository](https://github.com/CLiMA/Oceananigans.jl) by
  clicking the "Fork" button.

* Clone your fork of the Oceananigans repository (in terminal on Mac/Linux or git shell/
  GUI on Windows) in the location you'd like to keep it.
  ```
  git clone https://github.com/your-user-name/Oceananigans.jl.git
  ```

* Navigate to that folder in the terminal or in Anaconda Prompt if you're on Windows.

* Connect your repository to the upstream (main project).
  ```
  git remote add oceananigans https://github.com/CLiMA/Oceananigans.jl.git
  ```

* Create the development environment by opening Julia via `julia --project` then
  typing in `] instantiate`. This will install all the dependencies in the Project.toml
  file.

* You can test to make sure Oceananigans works by typing in `] test`. Doing so will run all
  the tests (and this can take a while).

Your development environment is now ready!

## Pull Requests

We follow the [ColPrac guide](https://github.com/SciML/ColPrac) for collaborative practices.
We ask that new contributors read that guide before submitting a pull request.

Changes and contributions should be made via GitHub pull requests against the ``master`` branch.

When you're done making changes, commit the changes you made. Chris Beams has written a 
[guide](https://chris.beams.io/posts/git-commit/) on how to write good commit messages.

When you think your changes are ready to be merged into the main repository, push to your fork
and [submit a pull request](https://github.com/CLiMA/Oceananigans.jl/compare/).

**Working on your first Pull Request?** You can learn how from this _free_ video series
[How to Contribute to an Open Source Project on GitHub](https://egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github), Aaron Meurer's [tutorial on the git workflow](https://www.asmeurer.com/git-workflow/), or the guide [“How to Contribute to Open Source"](https://opensource.guide/how-to-contribute/).

## Documentation

Now that you've made your awesome contribution, it's time to tell the world how to use it.
Writing documentation strings is really important to make sure others use your functionality
properly. Didn't write new functions? That's fine, but be sure that the documentation for
the code you touched is still in great shape. It is not uncommon to find some strange wording
or clarification that you can take care of while you are here.

You can preview how the Documentation will look like after merging by building the documentation 
locally. From the main directory of your local repository call

```
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project=docs/ -e 'using Pkg; Pkg.instantiate()'
JULIA_DEBUG=Documenter julia --project=docs/ docs/make.jl
```

and then open `docs/build/index.html` in your favorite browser. Providing the environment variable 
`JULIA_DEBUG=Documenter` will provide with more information in the documentation build process and
thus help figuring out a potential bug.

## Credits

This contributor's guide is heavily based on the excellent [MetPy contributor's guide](https://github.com/Unidata/MetPy/blob/master/CONTRIBUTING.md).
# Installation instructions

You can install the latest version of Oceananigans using the built-in package manager (accessed by pressing `]` in the
Julia command prompt) to add the package and instantiate/build all dependencies

```julia
julia>]
(v1.6) pkg> add Oceananigans
(v1.6) pkg> instantiate
```

We recommend installing Oceananigans with the built-in Julia package manager, because this installs a stable, tagged
release. Oceananigans.jl can be updated to the latest tagged release from the package manager by typing

```julia
(v1.6) pkg> update Oceananigans
```

At this time, updating should be done with care, as Oceananigans is under rapid development and breaking changes to the user API occur often. But if anything does happen, please open an issue!

But if anything does happen or your code stops working, please open an issue and ask! We're more than happy to help with getting your simulations up and running.

!!! compat "Julia 1.6 or newer"
    The latest version of Oceananigans requires at least Julia v1.6 to run.
    Installing Oceananigans with an older version of Julia will install an older version of Oceananigans (the latest version compatible with your version of Julia).
# Gallery

Collection of cool movies!

## [Deep convection](https://www.youtube.com/watch?v=kpUrxnKKMjI)

An idealized simulation of deep convection in the ocean. The simulation employs a resolution of 256x256x128 volumes in
a 2x2x1 km horizontally periodic domain. Heat is sucked out of the ocean surface within a cooling disk of radius 600 m
at a rate of 800 W/m² which cools the surface water and making it denser. This cold dense water then sinks into the
ocean interior, initiating a convective process that penetrates deep into the ocean.

This deep convection process can happen when a cold storm passes through warmer waters, which happens for example in
the Labrador Sea.

The video shows the temperature field and the domain is sliced in half so the convection happening under the cooling
disk is clear.

[![Watch deep convection in action](https://raw.githubusercontent.com/ali-ramadhan/ali-ramadhan.Github.io/master/img/surface_temp_3d_00130_halfsize.png)](https://www.youtube.com/watch?v=kpUrxnKKMjI)


## [Free convection](https://www.youtube.com/watch?v=yq4op9h3xcU)

An idealized simulation of free convection in the ocean. The simulation employs a resolution of 256x256x256 volumes in
a 100x100x100 m horizontally periodic domain. Heat is sucked out of the ocean surface at a rate of 75 W/m² which cools
the surface water and making it denser. This cold dense water then sinks into the ocean interior, initiating a
convective process that keeps mixing the upper layer of the ocean. This "mixed layer" has a relatively constant
temperature and keeps deepening as the surface is cooled.

The video shows the temperature field and the domain is sliced in half.

[![Watch free convection in action](https://raw.githubusercontent.com/ali-ramadhan/ali-ramadhan.Github.io/master/img/free_convection_0956.png)](https://www.youtube.com/watch?v=yq4op9h3xcU)


## [Winds blowing over the ocean](https://www.youtube.com/watch?v=IRncfbvuiy8)

An idealized simulation of a strong wind stress acting on the surface of a stratified ocean. The simulation employs a
resolution of 256x256x256 volumes in a 100x100x100 m horizontally periodic domain. A pretty strong wind stress of
0.1 N/m² is applied in the x direction which mechanically mixes the upper layer of the ocean. This leads to a "mixed
layer" of constant temperature near the surface of the ocean. You can also see the onset of Kelvin-Helmholtz
instabilities as the mechanical mixing sets in.

The video shows the temperature field in the top 25 meters and the domain is sliced in half for visualization. The line
plots show the horizontally averaged temperature profile (top right), horizontally averaged turbulent kinetic energy
(middle right), and the horizontally averaged buoyancy flux (or temperature flux).

[![Watch winds blowing over the ocean](https://raw.githubusercontent.com/ali-ramadhan/ali-ramadhan.Github.io/master/img/wind_stress_0400.png)](https://www.youtube.com/watch?v=IRncfbvuiy8)


## [Free convection with wind stress](https://www.youtube.com/watch?v=ob6OMQgPfI4)

An idealized simulation of a strong wind stress acting on the surface of a stratified ocean along with a cooling flux
that sucks heat out of the surface. The simulation employs a resolution of 256x256x256 volumes in a 100x100x100 m
horizontally periodic domain. A pretty strong wind stress of 0.1 N/m² is applied in the x direction which mechanically
mixes the upper layer of the ocean. Also, heat is sucked out of the ocean surface at a rate of 75 W/m² which cools the
surface water and making it denser. This cold dense water then sinks into the ocean interior, initiating a convective
process that keeps mixing the upper layer of the ocean. This leads to a "mixed layer" of constant temperature near the
surface of the ocean. You can also see the onset of Kelvin-Helmholtz instabilities as the mechanical mixing sets in.

The video shows the temperature field and the domain is sliced in half for visualization. The line plots show the
horizontally averaged temperature profile (top right), horizontally averaged turbulent kinetic energy (middle right),
and the horizontally averaged buoyancy flux (or temperature flux). The unusual periodic prism colormap is used to show
the fine details at the surface as it cools and the layers of different temperatures (the isopycnals) being perturbed
by internal waves.

[![Watch free convection with wind stress in action](https://raw.githubusercontent.com/ali-ramadhan/ali-ramadhan.Github.io/master/img/wind_stress_unstable_7500.png)](https://www.youtube.com/watch?v=ob6OMQgPfI4)

# Simulation tips

Oceananigans attemps to optimize as much of a computation as possible "behind the scenes".
Yet Oceananigans' flexibility places some responsibility on users to ensure high performance simulations,
especially for complex setups with user-defined forcing functions, boundary condition functions, and diagnostics.
Furthermore, in case of more complex GPU runs, some details could
sometimes prevent your simulation from running altogether. While Julia knowledge is obviously
desirable here, a user that is unfamiliar with Julia can get away with efficient simulations by
learning a few rules of thumb. It is nonetheless recommended that users go through Julia's
[performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/), which contains more
in-depth explanations of some of the aspects discussed here.

## General (CPU/GPU) simulation tips

### Avoid global variables whenever possible

In general using a [global
variable](https://docs.julialang.org/en/v1/manual/variables-and-scoping/#Global-Scope) (which can be
loosely defined as a variable defined in the main script) inside functions slows down the code. One
way to circumvent this is to always [use local variables or pass them as arguments to
functions](https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-global-variables). This
helps the compiler optimize the code.

Another way around this is to [define global variables as constants whenever
possible](https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-global-variables). One
thing to keep in mind when doing this is that when a `const` is defined, its value can't be changed
until you restart the Julia session. So this latter approach is good for production-ready code, but
may be undesirable in the early stages of development while you still have to change the parameters
of the simulation for exploration.

It is especially important to avoid global variables in functions that are meant to be executed in
GPU kernels (such as functions defining boundary conditions and forcings). Otherwise the Julia GPU
compiler can fail with obscure errors. This is explained in more detail in the GPU simulation tips
section below.

### Consider inlining small functions

Inlining is when the compiler [replaces a function call with the body of the function that is being
called before compiling](https://en.wikipedia.org/wiki/Inline_expansion). The advantage of inlining
(which in julia can be done with the [`@inline`
macro](https://docs.julialang.org/en/v1/devdocs/meta/)) is that gets rid of the time spent calling
the function. The Julia compiler automatically makes some calls as to what functions it should or
shouldn't inline, but you can force a function to be inlined by including the macro `@inline` before
its definition. This is more suited for small functions that are called often. Here's an example of
an implementation of the Heaviside function that forces it to be inlined:

```julia
@inline heaviside(X) = ifelse(X < 0, zero(X), one(X))
```

In practice it's hard to say whether inlining a function will bring runtime benefits _with
certainty_, since Julia and KernelAbstractions.jl (needed for GPU runs) already inline some
functions automatically. However, it is generally a good idea to at least investigate this aspect in
your code as the benefits can potentially be significant.

## GPU simulation tips

Running on GPUs can be very different from running on CPUs. Oceananigans makes most of the necessary
changes in the background, so that for very simple simulations changing between CPUs and GPUs is
just a matter of changing the `architecture` argument in the model from `CPU()` to `GPU()`. However,
for more complex simulations some care needs to be taken on the part of the user. While knowledge of
GPU computing (and Julia) is again desirable, an inexperienced user can also achieve high efficiency
in GPU simulations by following a few simple principles.

### Global variables that need to be used in GPU computations need to be defined as constants or passed as parameters

Any global variable that needs to be accessed by the GPU needs to be a constant or the simulation
will crash. This includes any variables that are referenced as global variables in functions
used for forcing of boundary conditions. For example,

```julia
T₀ = 20 # ᵒC
surface_temperature(x, y, t) = T₀ * sin(2π / 86400 * t)
T_bcs = FieldBoundaryConditions(bottom = GradientBoundaryCondition(surface_temperature))
```

will throw an error if run on the GPU (and will run more slowly than it should on the CPU).
Replacing the first line above with

```julia
const T₀ = 20 # ᵒC
```

fixes the issue by indicating to the compiler that `T₀` will not change.

Note that the _literal_ `2π / 86400` is not an issue -- it's only the
_variable_ `T₀` that must be declared `const`.

Alternatively, passing the variable as a parameter to `GradientBoundaryCondition` also works:

```julia
T₀ = 20 # ᵒC
surface_temperature(x, y, t, p) = p.T₀ * sin(2π / 86400 * t)
T_bcs = FieldBoundaryConditions(bottom = GradientBoundaryCondition(surface_temperature, parameters=(T₀=T₀,)))
```

### Complex diagnostics using computed `Field`s may not work on GPUs

`Field`s are the most convenient way to calculate diagnostics for your simulation. They will
always work on CPUs, but when their complexity is high (in terms of number of abstract operations)
the compiler can't translate them into GPU code and they fail for GPU runs. (This limitation is summarized 
in [this Github issue](https://github.com/CliMA/Oceananigans.jl/issues/1886) and contributions are welcome.)
For example, in the example below, calculating `u²` works in both CPUs and GPUs, but calculating 
`ε` will not compile on GPUs when we call the command `compute!`:

```julia
using Oceananigans
grid = RectilinearGrid(size=(4, 4, 4), extent=(1, 1, 1))
model = NonhydrostaticModel(grid=grid, closure=IsotropicDiffusivity(ν=1e-6))
u, v, w = model.velocities
ν = model.closure.ν
u² = Field(u^2)
ε = Field(ν*(∂x(u)^2 + ∂x(v)^2 + ∂x(w)^2 + ∂y(u)^2 + ∂y(v)^2 + ∂y(w)^2 + ∂z(u)^2 + ∂z(v)^2 + ∂z(w)^2))
compute!(u²)
compute!(ε)
```

There are a few ways to work around this issue.
One is to compute `ε` in steps by nesting computed `Field`s,
```julia
ddx² = Field(∂x(u)^2 + ∂x(v)^2 + ∂x(w)^2)
ddy² = Field(∂y(u)^2 + ∂y(v)^2 + ∂y(w)^2)
ddz² = Field(∂z(u)^2 + ∂z(v)^2 + ∂z(w)^2)
ε = Field(ν * (ddx² + ddy² + ddz²))
compute!(ε)
```

This method is expensive because it requires computing and storing 3 intermediate terms.
`ε` may also be calculated via `KernelFunctionOperations`s, which
requires explicitly building a "kernel function" from low-level Oceananigans
operators.

```julia
using Oceananigans.Operators
using Oceananigans.AbstractOperations: KernelFunctionOperation

@inline fψ_plus_gφ²(i, j, k, grid, f, ψ, g, φ) = @inbounds (f(i, j, k, grid, ψ) + g(i, j, k, grid, φ))^2

function isotropic_viscous_dissipation_rate_ccc(i, j, k, grid, u, v, w, ν)
    Σˣˣ² = ∂xᶜᵃᵃ(i, j, k, grid, u)^2
    Σʸʸ² = ∂yᵃᶜᵃ(i, j, k, grid, v)^2
    Σᶻᶻ² = ∂zᵃᵃᶜ(i, j, k, grid, w)^2

    Σˣʸ² = ℑxyᶜᶜᵃ(i, j, k, grid, fψ_plus_gφ², ∂yᵃᶠᵃ, u, ∂xᶠᵃᵃ, v) / 4
    Σˣᶻ² = ℑxzᶜᵃᶜ(i, j, k, grid, fψ_plus_gφ², ∂zᵃᵃᶠ, u, ∂xᶠᵃᵃ, w) / 4
    Σʸᶻ² = ℑyzᵃᶜᶜ(i, j, k, grid, fψ_plus_gφ², ∂zᵃᵃᶠ, v, ∂yᵃᶠᵃ, w) / 4

    return ν * 2 * (Σˣˣ² + Σʸʸ² + Σᶻᶻ² + 2 * (Σˣʸ² + Σˣᶻ² + Σʸᶻ²))
end

ε_op = KernelFunctionOperation{Center, Center, Center}(isotropic_viscous_dissipation_rate_ccc,
                                                       grid;
                                                       computed_dependencies=(u, v, w, ν))

ε = Field(ε_op)

compute!(ε)
```

Writing kernel functions like `isotropic_viscous_dissipation_rate_ccc`
requires understanding the C-grid, but incurs only one iteration over the domain.

`KernelFunctionOperation`s for some diagnostics common to large eddy simulation are defined in
[Oceanostics.jl](https://github.com/tomchor/Oceanostics.jl/blob/3b8f67338656557877ef8ef5ebe3af9e7b2974e2/src/TurbulentKineticEnergyTerms.jl#L35-L57),

```julia
using Oceanostics: IsotropicPseudoViscousDissipationRate
ε = IsotropicViscousDissipationRate(model, u, v, w, ν)
compute!(ε)
```
[Start an issue on Github](https://github.com/CliMA/Oceananigans.jl/issues/new) if more help is needed.


### Try to decrease the memory-use of your runs

GPU runs are sometimes memory-limited. A state-of-the-art Tesla V100 GPU has 32GB of
memory --- enough memory for simulations with about 100 million points, or grids a bit smaller
than 512 × 512 × 512. (The maximum grid size depends on some user-specified factors,
like the number of passive tracers or computed diagnostics.)
For large simulations on the GPU, careful management of memory allocation may be required:

- Use the [`nvidia-smi`](https://developer.nvidia.com/nvidia-system-management-interface) command
  line utility to monitor the memory usage of the GPU. It should tell you how much memory there is
  on your GPU and how much of it you're using and you can run it from Julia with the command ``run(`nvidia-smi`)``.

- Try to use higher-order advection schemes. In general when you use a higher-order scheme you need
  fewer grid points to achieve the same accuracy that you would with a lower-order one. Oceananigans
  provides two high-order advection schemes: 5th-order WENO method (WENO5) and 3rd-order upwind.

- Manually define scratch space to be reused in diagnostics. By default, every time a user-defined
  diagnostic is calculated the compiler reserves a new chunk of memory for that calculation, usually
  called scratch space. In general, the more diagnostics, the more scratch space needed and the bigger
  the memory requirements. However, if you explicitly create a scratch space and pass that same
  scratch space for as many diagnostics as you can, you minimize the memory requirements of your
  calculations by reusing the same chunk of memory. As an example, you can see scratch space being
  created
  [here](https://github.com/CliMA/LESbrary.jl/blob/cf31b0ec20219d5ad698af334811d448c27213b0/examples/three_layer_ constant_fluxes.jl#L380-L383)
  and then being used in calculations
  [here](https://github.com/CliMA/LESbrary.jl/blob/cf31b0ec20219d5ad698af334811d448c27213b0/src/TurbulenceStatistics/first_through_third_order.jl#L109-L112).

### Arrays in GPUs are usually different from arrays in CPUs

Oceananigans.jl uses [`CUDA.CuArray`](https://cuda.juliagpu.org/stable/usage/array/) to store 
data for GPU computations. One limitation of `CuArray`s compared to the `Array`s used for 
CPU computations is that `CuArray` elements in general cannot be accessed outside kernels
launched through CUDA.jl or KernelAbstractions.jl. (You can learn more about GPU kernels 
[here](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#kernels) and 
[here](https://cuda.juliagpu.org/stable/usage/overview/#Kernel-programming-with-@cuda).)
Doing so requires individual elements to be copied from or to the GPU for processing,
which is very slow and can result in huge slowdowns. For this reason, Oceananigans.jl disables CUDA
scalar indexing by default. See the [scalar indexing](https://juliagpu.github.io/CUDA.jl/dev/usage/workflow/#UsageWorkflowScalar)
section of the CUDA.jl documentation for more information on scalar indexing.

For example if can be difficult to just view a `CuArray` since Julia needs to access 
its elements to do that. Consider the example below:

```julia
julia> using Oceananigans; using Adapt

julia> grid = RectilinearGrid(size=(1,1,1), extent=(1,1,1))
RectilinearGrid{Float64, Periodic, Periodic, Bounded}
                   domain: x ∈ [0.0, 1.0], y ∈ [0.0, 1.0], z ∈ [-1.0, 0.0]
                 topology: (Periodic, Periodic, Bounded)
        size (Nx, Ny, Nz): (1, 1, 1)
        halo (Hx, Hy, Hz): (1, 1, 1)
grid spacing (Δx, Δy, Δz): (1.0, 1.0, 1.0)

julia> model = NonhydrostaticModel(grid=grid, architecture=GPU())
NonhydrostaticModel{GPU, Float64}(time = 0 seconds, iteration = 0) 
├── grid: RectilinearGrid{Float64, Periodic, Periodic, Bounded}(Nx=1, Ny=1, Nz=1)
├── tracers: (:T, :S)
├── closure: IsotropicDiffusivity{Float64,NamedTuple{(:T, :S),Tuple{Float64,Float64}}}
├── buoyancy: SeawaterBuoyancy{Float64,LinearEquationOfState{Float64},Nothing,Nothing}
└── coriolis: Nothing

julia> typeof(model.velocities.u.data)
OffsetArrays.OffsetArray{Float64,3,CUDA.CuArray{Float64,3}}

julia> adapt(Array, model.velocities.u.data)
3×3×3 OffsetArray(::Array{Float64,3}, 0:2, 0:2, 0:2) with eltype Float64 with indices 0:2×0:2×0:2:
[:, :, 0] =
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0

[:, :, 1] =
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0

[:, :, 2] =
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
```

Notice that in order to view the `CuArray` that stores values for `u` we needed to transform
it into a regular `Array` first using `Adapt.adapt`. If we naively try to view the `CuArray`
without that step we get an error:

```julia
julia> model.velocities.u.data
3×3×3 OffsetArray(::CUDA.CuArray{Float64,3}, 0:2, 0:2, 0:2) with eltype Float64 with indices 0:2×0:2×0:2:
[:, :, 0] =
Error showing value of type OffsetArrays.OffsetArray{Float64,3,CUDA.CuArray{Float64,3}}:
ERROR: scalar getindex is disallowed
```

Here `CUDA.jl` throws an error because scalar `getindex` is not `allowed`. Another way around 
this limitation is to allow scalar operations on `CuArray`s. We can temporarily
do that with the `CUDA.@allowscalar` macro or by calling `CUDA.allowscalar(true)`.


```julia
julia> using CUDA; CUDA.allowscalar(true)

julia> model.velocities.u.data
3×3×3 OffsetArray(::CuArray{Float64,3}, 0:2, 0:2, 0:2) with eltype Float64 with indices 0:2×0:2×0:2:
[:, :, 0] =
┌ Warning: Performing scalar operations on GPU arrays: This is very slow, consider disallowing these operations with `allowscalar(false)`
└ @ GPUArrays ~/.julia/packages/GPUArrays/WV76E/src/host/indexing.jl:43
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0

[:, :, 1] =
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0

[:, :, 2] =
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
```

Notice the warning we get when we do this. Scalar operations on GPUs can be very slow, so it is
advised to only use this last method when using the REPL or prototyping --- never in
production-ready scripts.


You might also need to keep these differences in mind when using arrays
to define initial conditions, boundary conditions or
forcing functions on a GPU. To learn more about working with `CuArray`s, see the
[array programming](https://juliagpu.github.io/CUDA.jl/dev/usage/array/) section
of the CUDA.jl documentation.
# Time-stepping and the fractional step method

The time-integral of the momentum equation with the pressure decomposition from time step ``n`` at ``t = t_n`` 
to time step ``n+1`` at ``t_{n+1}`` is
```math
    \begin{equation}
    \label{eq:momentum-time-integral}
    \boldsymbol{v}^{n+1} - \boldsymbol{v}^n = 
        \int_{t_n}^{t_{n+1}} \Big [ - \boldsymbol{\nabla} p_{\rm{non}} 
                                    - \boldsymbol{\nabla}_{h} p_{\rm{hyd}} 
                                    - \left ( \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} \right ) \boldsymbol{v} 
                                    - \boldsymbol{f} \times \boldsymbol{v} 
                                    + \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{\tau} 
                                    + \boldsymbol{F}_{\boldsymbol{v}} \Big ] \, \mathrm{d} t \, ,
    \end{equation}
```
where the superscript ``n`` and ``n+1`` imply evaluation at ``t_n`` and ``t_{n+1}``, such that 
``\boldsymbol{v}^n \equiv \boldsymbol{v}(t=t_n)``. The crux of the fractional step method is 
to treat the pressure term ``\boldsymbol{\nabla} p_{\rm{non}}`` implicitly using the approximation
```math
\int_{t_n}^{t_{n+1}} \boldsymbol{\nabla} p_{\rm{non}} \, \mathrm{d} t \approx 
    \Delta t \boldsymbol{\nabla} p_{\rm{non}}^{n+1} \, ,
```
while treating the rest of the terms on the right hand side of \eqref{eq:momentum-time-integral} 
explicitly. The implicit treatment of pressure ensures that the velocity field obtained at 
time step ``n+1`` is divergence-free.

To effect such a fractional step method, we define an intermediate velocity field ``\boldsymbol{v}^\star`` such that
```math
    \begin{equation}
    \label{eq:intermediate-velocity-field}
    \boldsymbol{v}^\star - \boldsymbol{v}^n = \int_{t_n}^{t_{n+1}} \boldsymbol{G}_{\boldsymbol{v}} \, \mathrm{d} t \, ,
    \end{equation}
```
where, e.g., for the non-hydrostatic model, 
```math
\boldsymbol{G}_{\boldsymbol{v}} \equiv - \boldsymbol{\nabla}_h p_{\rm{hyd}} 
                       - \left ( \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} \right ) \boldsymbol{v} 
                       - \boldsymbol{f} \times \boldsymbol{v} 
                       + \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{\tau} 
                       + \boldsymbol{F}_{\boldsymbol{v}}
```
collects all terms on the right side of the time-integral of the momentum equation except the 
contribution of non-hydrostatic pressure ``\boldsymbol{\nabla} p_n``. The integral on the right 
of the equation for ``\boldsymbol{v}^\star`` may be approximated by a variety of  explicit methods: 
for example, a forward Euler method uses
```math
    \begin{equation}
    \int_{t_n}^{t_{n+1}} G \, \mathrm{d} t \approx \Delta t G^n \, ,
    \label{eq:forward-euler}
    \end{equation}
```
for any time-dependent function ``G(t)``, while a second-order Adams-Bashforth method uses the approximation
```math
    \begin{equation}
    \label{eq:adams-bashforth}
    \int_{t_n}^{t_{n+1}} G \, \mathrm{d} t \approx 
        \Delta t \left [ \left ( \tfrac{3}{2} + \chi \right ) G^n 
        - \left ( \tfrac{1}{2} + \chi \right ) G^{n-1} \right ] \, ,
    \end{equation}
```
where ``\chi`` is a parameter. [Ascher95](@cite) claim that ``\chi = \tfrac{1}{8}`` is optimal; 
``\chi = -\tfrac{1}{2}`` yields the forward Euler scheme.

Combining the equations for ``\boldsymbol{v}^\star`` and the time integral of the momentum equation yields
```math
    \begin{equation}
    \label{eq:fractional-step}
    \boldsymbol{v}^{n+1} - \boldsymbol{v}^\star = - \Delta t \boldsymbol{\nabla} p_{\rm{non}}^{n+1} \, .
    \end{equation}
```
Taking the divergence of fractional step equation and requiring that 
``\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{v}^{n+1} = 0`` yields a Poisson equation 
for the kinematic pressure ``p_{\rm{non}}`` at time-step ``n+1``:
```math
    \nabla^2 p_{\rm{non}}^{n+1} = \frac{\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{v}^{\star}}{\Delta t} \, .
```
With ``\boldsymbol{v}^\star`` and ``p_{\rm{non}}`` in hand, ``\boldsymbol{v}^{n+1}`` is then 
computed via the fractional step equation.

Tracers are stepped forward explicitly via
```math
    \begin{equation}
    \label{eq:tracer-timestep}
    c^{n+1} - c^n = \int_{t_n}^{t_{n+1}} G_c \, \mathrm{d} t \, ,
    \end{equation}
```
where 
```math
    G_c \equiv - \boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{v} c \right ) - \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{q}_c + F_c \, ,
```
and the same forward Euler or Adams-Bashforth scheme as for the explicit evaluation of the time-integral of
``\boldsymbol{G}_u`` is used to evaluate the integral of ``G_c``.
# [Large eddy simulation](@id numerical_les)

The idea behind large eddy simulation (LES) is to resolve the "large eddies" while modeling the effect of unresolved
sub-grid scale motions. This is done usually be assuming eddy viscosity and eddy diffusivity models and providing an
estimate for the eddy viscosity ``\nu_e`` and diffusivity ``\kappa_e``.

Much of the early work on LES was motivated by the study of atmospheric boundary layer turbulence, being developed
by [Smagorinsky63](@cite) and [Lilly66](@cite), then first implemented by [Deardorff70](@cite) and [Deardorff74](@cite).

In the LES framework, the Navier-Stokes equations are averaged in the same way as [Reynolds1895](@cite) except that the
mean field ``\overline{\boldsymbol{v}}`` is obtained via convolution with a filter convolution kernel ``G``
```math
\overline{\boldsymbol{v}(\boldsymbol{x}, t)} = G \star \boldsymbol{v} =
  \int_{-\infty}^\infty \int_{-\infty}^\infty
  \boldsymbol{v}(\boldsymbol{x}^\prime, t) G(\boldsymbol{x} - \boldsymbol{x}^\prime, t - \tau) \, d\boldsymbol{x}^\prime \, \mathrm{d} \tau \, ,
```
as described by [Leonard75](@cite) who introduced the general filtering formalism.

The ``\overline{v_i^\prime v_j^\prime}`` terms are now components of what is called the sub-grid scale (SGS) stress
tensor ``\tau^\text{SGS}_{ij}``, which looks the same as the Reynolds stress tensor so we will drop the SGS superscript.

It is probably important to note that the large eddy simulation filtering operation does not satisfy the properties
of a Reynolds operator (§2.1)[sagaut06](@cite) and that in general, the filtered residual is not zero:
``\overline{\boldsymbol{v}^\prime(\boldsymbol{x}, t)} \ne 0``.

§13.2 of [Pope00](@cite) lists a number of popular choices for the filter function ``G``. For practical reasons we
simply employ the box kernel
```math
  \begin{equation}
  \label{eq:box-kernel}
  G_\Delta = G(\boldsymbol{x}, t) = \frac{1}{\Delta} H \left( \frac{1}{2}\Delta - |\boldsymbol{x}| \right) \delta(t - t_n) \, ,
  \end{equation}
```
where ``H`` is the Heaviside function, ``\Delta`` is the grid spacing, and ``t_n`` is the current time step. With
\eqref{eq:box-kernel} we get back the averaging operator originally used by [Deardorff70](@cite)
```math
\overline{\boldsymbol{v}(x, y, z, t)} =
  \frac{1}{\Delta x \Delta y \Delta z}
  \int_{x - \frac{1}{2}\Delta x}^{x + \frac{1}{2}\Delta x}
  \int_{y - \frac{1}{2}\Delta y}^{y + \frac{1}{2}\Delta y}
  \int_{z - \frac{1}{2}\Delta z}^{z + \frac{1}{2}\Delta z}
  \boldsymbol{v}(\xi, \eta, \zeta, t) \, \mathrm{d} \xi \, \mathrm{d} \eta \, \mathrm{d} \zeta \, ,
```
which if evaluated at the cell centers just returns the cell averages we already compute in the finite volume method.


## Smagorinsky-Lilly model

[Smagorinsky63](@cite) estimated the eddy viscosity ``\nu_e`` via a characteristic length scale ``\Delta`` times a velocity
scale given by ``\Delta |\overline{S}|`` where ``|\overline{S}| = \sqrt{2\overline{S}_{ij}\overline{S}_{ij}}``. Thus the
SGS stress tensor is given by
```math
\tau_{ij} = -2 \nu_e \overline{S}_{ij} = -2 (C_s \Delta)^2 |\overline{S}| \overline{S}_{ij} \, ,
```
where ``C_s`` is a dimensionless constant. The grid spacing is usually used for the characteristic length scale ``\Delta``.
The eddy diffusivities are calculated via ``\kappa_e = \nu_e / \text{Pr}_t`` where the turbulent Prandtl number
``\text{Pr}_t`` is usually chosen to be ``\mathcal{O}(1)`` from experimental observations.

Assuming that the SGS energy cascade is equal to the overall dissipation rate ``\varepsilon`` from the
[Kolmogorov41](@cite) theory, [Lilly66](@cite) was able to derive a value of
```math
C_s = \left( \frac{3}{2}C_K\pi^\frac{4}{3} \right)^{-\frac{3}{4}} \approx 0.16 \, ,
```
using an empirical value of ``C_K \approx 1.6`` for the Kolmogorov constant. This seems reasonable for isotropic
turbulence if the grid spacing ``\Delta`` falls in the inertial range. In practice, ``C_s`` is a tunable parameter.

Due to the presence of the constant ``C_s``, the model is sometimes referred to as the *constant Smagorinsky* model
in contrast to *dynamic Smagorinsky* models that dynamically compute ``C_s`` to account for effects such as buoyant
convection.

## Anisotropic minimum dissipation models

Minimum-dissipation eddy-viscosity models are a class of LES closures that use the minimum eddy dissipation required to
dissipate the energy of sub-grid scale motion. [Rozema15](@cite) proposed the first minimum-dissipation model
appropriate for use on anisotropic grids, termed the *anisotropic minimum dissipation* (AMD) model.

It has a number of desirable properties over Smagorinsky-type closures: it is more cost-effective than dynamic
Smagorinsky, it appropriately switches off in laminar and transitional flows, and it is consistent with the exact SGS
stress tensor on both isotropic and anisotropic grids. [Abkar16](@cite) extended the AMD model to model SGS scalar
fluxes for tracer transport. [Abkar17](@cite) further extended the model to include a buoyancy term that accounts for
the contribution of buoyant forces to the production and suppression of turbulence.

[Vreugdenhil18](@cite) derive a modified AMD model by following the requirement suggested by [Verstappen18](@cite),
which entail normalising the displacement, the velocity, and the velocity gradient by the filter width to ensure that
the resulting eddy dissipation properly counteracts the spurious kinetic energy transferred by convective nonlinearity,
to derive a modified AMD model.

The eddy viscosity and diffusivity are defined in terms of eddy viscosity and diffusivity *predictors*
``\nu_e^\dagger`` and ``\kappa_e^\dagger``, such that
```math
\nu_e = \max \lbrace 0, \nu_e^\dagger \rbrace
\quad \text{and} \quad
\kappa_e = \max \lbrace 0, \kappa_e^\dagger \rbrace \, ,
```
to ensure that ``\nu_e \ge 0`` and ``\kappa_e \ge 0``. Leaving out the overlines and understanding that all variables
represent the resolved/filtered variables, the eddy viscosity predictor is given by
```math
    \begin{equation}
    \label{eq:nu-dagger}
    \nu_e^\dagger = -(C\Delta)^2
      \frac
        {\left( \hat{\partial}_k \hat{v}_i \right) \left( \hat{\partial}_k \hat{v}_j \right) \hat{S}_{ij}
        + C_b\hat{\delta}_{i3} \alpha g \left( \hat{\partial}_k \hat{v_i} \right) \hat{\partial}_k \theta}
        {\left( \hat{\partial}_l \hat{v}_m \right) \left( \hat{\partial}_l \hat{v}_m \right)} \, ,
    \end{equation}
```
and the eddy diffusivity predictor by
```math
    \begin{equation}
    \kappa_e^\dagger = -(C\Delta)^2
    \frac
        {\left( \hat{\partial}_k \hat{v}_i \right) \left( \hat{\partial}_k \hat{\theta} \right) \hat{\partial}_i \theta}
        {\left( \hat{\partial}_l \hat{\theta} \right) \left( \hat{\partial}_l \hat{\theta} \right)} \, ,
    \end{equation}
```
where
```math
  \begin{equation}
  \hat{x}_i = \frac{x_i}{\Delta_i}, \quad
  \hat{v}_i(\hat{x}, t) = \frac{v_i(x, t)}{\Delta_i}, \quad
  \hat{\partial}_i \hat{v}_j(\hat{x}, t) = \frac{\Delta_i}{\Delta_j} \partial_i v_j(x, t), \quad
  \hat{\delta}_{i3} = \frac{\delta_{i3}}{\Delta_3} \, ,
  \end{equation}
```
so that the normalized rate of strain tensor is
```math
    \begin{equation}
    \label{eq:S-hat}
    \hat{S}_{ij} =
      \frac{1}{2} \left[ \hat{\partial}_i \hat{v}_j(\hat{x}, t) + \hat{\partial}_j \hat{v}_i(\hat{x}, t) \right] \, .
    \end{equation}
```

In equations \eqref{eq:nu-dagger}--\eqref{eq:S-hat}, ``C`` is a modified Poincaré "constant" that is independent from
the filter width ``\Delta`` but does depend on the accuracy of the discretization method used. [Abkar16](@cite) cite
``C^2 = \frac{1}{12}`` for a spectral method and ``C^2 = \frac{1}{3}`` for a second-order accurate scheme. ``\Delta_i`` is
the filter width in the ``x_i``-direction, and ``\Delta`` is given by the square root of the harmonic mean of the squares
of the filter widths in each direction
```math
    \frac{1}{\Delta^2} = \frac{1}{3} \left( \frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2} \right) \, .
```
The term multiplying ``C_b`` is the buoyancy modification introduced by [Abkar17](@cite) and is small for weakly
stratified flows. We have introduced the ``C_b`` constant so that the buoyancy modification term may be turned on and off.
# Poisson solvers

## The elliptic problem for the pressure

The 3D non-hydrostatic pressure field is obtained by taking the divergence of the horizontal 
component of the momentum equation and invoking the vertical component to yield an elliptic 
Poisson equation for the non-hydrostatic kinematic pressure
```math
   \begin{equation}
   \label{eq:poisson-pressure}
   \nabla^2 p_{NH} = \frac{\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{v}^n}{\Delta t} + \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{G}_{\boldsymbol{v}} \equiv \mathscr{F} \, ,
   \end{equation}
```
along with homogenous Neumann boundary conditions ``\boldsymbol{v} \cdot \boldsymbol{\hat{n}} = 0`` 
(Neumann on ``p`` for wall-bounded directions and periodic otherwise) and where ``\mathscr{F}`` 
denotes the source term for the Poisson equation.

For hydrostatic problems the Poisson equation above only needs to be solved for the vertically 
integrated flow and the pressure field is a two dimensional term ``p_{S}``. In this case a fully 
three-dimensional solve is not needed.

## Direct method

Discretizing elliptic problems that can be solved via a classical separation-of-variables approach, such as Poisson's
equation, results in a linear system of equations ``M \boldsymbol{x} = \boldsymbol{y}`` where ``M`` is a real symmetric matrix of block
tridiagonal form. This allows for the matrix to be decomposed and solved efficiently, provided that the eigenvalues and
eigenvectors of the blocks are known (§2) [Buzbee70](@cite). In the case of Poisson's equation on a rectangle,
[Hockney65](@cite) has taken advantage of the fact that the fast Fourier transform can be used to perform the matrix
multiplication steps resulting in an even more efficient method. [Schumann88](@cite) describe the implementation of such
an algorithm for Poisson's equation on a staggered grid with Dirichlet, Neumann, and periodic boundary conditions.

The method can be explained easily by taking the Fourier transform of both sides of \eqref{eq:poisson-pressure} to yield
```math
    \begin{equation}
    \label{eq:poisson-spectral}
    -(k_x^2 + k_y^2 + k_z^2) \widehat{p}_{NH} = \widehat{\mathscr{F}}
    \quad \implies \quad
    \widehat{p}_{NH} = - \frac{\widehat{\mathscr{F}}}{k_x^2 + k_y^2 + k_z^2} \, ,
    \end{equation}
```
where ``\widehat{\cdot}`` denotes the Fourier component. Here ``k_x``, ``k_y``, and ``k_z`` are the wavenumbers. However, when
solving the equation on a staggered grid we require a solution for ``p_{NH}`` that is second-order accurate such that
when when its Laplacian is computed, ``\nabla^2 p_{NH}`` matches ``\mathscr{F}`` to machine precision. This is crucial to
ensure that the projection step in §\ref{sec:fractional-step} works. To do this, the wavenumbers are replaced by
eigenvalues ``\lambda_x``, ``\lambda_y``, and ``\lambda_z`` satisfying the discrete form of Poisson's equation with
appropriate boundary conditions. Thus, Poisson's equation is diagonalized in Fourier space and the Fourier
coefficients of the solution are easily solved for
```math
\widehat{p}_{NH}(i, j, k) = - \frac{\widehat{\mathscr{F}}(i, j, k)}{\lambda^x_i + \lambda^y_j + \lambda^z_k} \, .
```

The eigenvalues are given by [Schumann88](@cite) and can also be tediously derived by plugging in the definition of the
discrete Fourier transform into \eqref{eq:poisson-spectral}:
```math
\begin{align}
    \lambda^x_i &= 4\frac{N_x^2}{L_x^2} \sin^2 \left [ \frac{(i-1)\pi}{N_x}  \right ], \quad i=0,1, \dots,N_x-1 \, , \\
    \lambda^x_j &= 4\frac{N_y^2}{L_y^2} \sin^2 \left [ \frac{(j-1)\pi}{N_y}  \right ], \quad j=0,1, \dots,N_y-1 \, , \\
    \lambda^x_k &= 4\frac{N_z^2}{L_z^2} \sin^2 \left [ \frac{(k-1)\pi}{2N_z} \right ], \quad k=0,1, \dots,N_z-1 \, ,
\end{align}
```
where ``\lambda_x`` and ``\lambda_y`` correspond to periodic boundary conditions in the horizontal and ``\lambda_z`` to
Neumann boundary conditions in the vertical.

There is also an ambiguity in the solution to Poisson's equation as it's only defined up to a constant. To resolve this
we choose the solution with zero mean by setting the zeroth Fourier coefficient ``p_{000}`` (corresponding to
``k_x = k_y = k_z = 0``) to zero. This also has the added benefit of discarding the zero eigenvalue so we don't divide by
it.

The Fast Fourier transforms are computed using FFTW.jl [[Frigo98](@cite) and [Frigo05](@cite)] on the CPU and using the
cuFFT library on the GPU. Along wall-bouded dimensions, the cosine transform is used. In particular, as the transforms
are performed on a staggered grid, DCT-II (`REDFT10`) is used to perform the forward cosine transform and DCT-III
(`REDFT01`) is used to perform the inverse cosine transform.

## Direct method with a vertically stretched grid

Using Fourier transforms for all three dimensions results in a method requiring ``\mathcal{O}(N \log_2 N)`` operations
where ``N`` is the total number of grid points. This algorithm can be made even more efficient by solving a tridiagonal
system along one of the dimensions and utilizing cyclic reduction. This results in the *Fourier analysis cyclic
reduction* or ``\text{FACR}(\ell)`` algorithm (with ``\ell`` cyclic reduction steps) which requires only
``\mathcal{O}(N \log_2\log_2 N)`` operations provided the optimal number of cyclic reduction steps is taken, which is
``\ell = \log_2 \log_2 n`` where ``n`` is the number of grid points in the cyclic reduction dimension. The FACR algorithm
was first developed by [Hockney69](@cite) and is well reviewed by [Swarztrauber77](@cite) then further benchmarked and
extended by [Temperton79](@cite) and [Temperton80](@cite).

Furthermore, the FACR algorithm removes the restriction that the grid is uniform in one of the dimensions so it can
be utilized to implement a fast Poisson solver for vertically stretched grids if the cyclic reduction is applied in the
along the vertical dimension.

Expanding ``p_{NH}`` and ``\mathscr{F}`` into Fourier modes along the ``x`` and ``y`` directions
```math
p_{ijk} = \sum_{m=1}^{N_x} \sum_{n=1}^{N_y} \tilde{p}_{mnk} \; e^{-\mathrm{i} 2\pi i m / N_x} \;  e^{-\mathrm{i} 2\pi j n / N_y} \, ,
```
and recalling that Fourier transforms do ``\partial_x \rightarrow \mathrm{i} k_x`` and ``\partial_y \rightarrow \mathrm{i} k_y`` we can write
\eqref{eq:poisson-pressure} as
```math
\sum_{m=1}^{N_x} \sum_{n=1}^{N_y}
\left\lbrace
    \partial_z^2 \tilde{p}_{mnk} - (k_x^2 + k_y^2) \tilde{p}_{mnk} - \tilde{\mathscr{F}}_{mnk}
\right\rbrace e^{-\mathrm{i} 2 \pi i m / N_x}  e^{-\mathrm{i} 2 \pi j n / N_y} = 0 \, .
```
Discretizing the ``\partial_z^2`` derivative and equating the term inside the brackets to zero we arrive at
``N_x\times N_y`` symmetric tridiagonal systems of ``N_z`` linear equations for the Fourier modes:
```math
\frac{\tilde{p}_{mn,k-1}}{\Delta z^C_k}
- \left\lbrace \frac{1}{\Delta z^C_k} + \frac{1}{\Delta z^C_{k+1}} + \Delta z^F_k (k_x^2 + k_y^2) \right\rbrace
  \tilde{p}_{mnk}
+ \frac{\tilde{p}_{mn,k+1}}{\Delta z^C_{k+1}}
= \Delta z^F_k \tilde{\mathscr{F}}_{mnk} \, .
```

## Cosine transforms on the GPU

Unfortunately cuFFT does not provide cosine transforms and so we must write our own fast cosine 
transforms for the GPU. We implemented the fast 1D and 2D cosine transforms described by [Makhoul80](@cite) 
which compute it by applying the regular Fourier transform to a permuted version of the array.

In this section we will be using the DCT-II as the definition of the forward cosine transform 
for a real signal of length ``N``
```math
    \begin{equation}
    \label{eq:FCT}
    \text{DCT}(X): \quad Y_k = 2 \sum_{j=0}^{N-1} \cos \left[ \frac{\pi(j + \frac{1}{2})k}{N} \right] X_j \, ,
    \end{equation}
```
and the DCT-III as the definition of the inverse cosine transform
```math
    \begin{equation}
    \label{eq:IFCT}
    \text{IDCT}(X): \quad Y_k = X_0 + 2 \sum_{j=1}^{N-1} \cos \left[ \frac{\pi j (k + \frac{1}{2})}{N} \right] X_j \, ,
    \end{equation}  
```
and will use ``\omega_M = e^{-2 \pi \mathrm{i} / M}`` to denote the ``M^\text{th}`` root of unity, sometimes called the twiddle factors
in the context of FFT algorithms.

### 1D fast cosine transform
To calculate \eqref{eq:FCT} using the fast Fourier transform, we first permute the input signal along the appropriate
dimension by ordering the odd elements first followed by the even elements to produce a permuted signal
```math
    X^\prime_n =
    \begin{cases}
        \displaystyle X_{2N}, \quad 0 \le n \le \left[ \frac{N-1}{2} \right] \, , \\
        \displaystyle X_{2N - 2n - 1}, \quad \left[ \frac{N+1}{2} \right] \le n \le N-1 \, ,
    \end{cases}
```
where ``[a]`` indicates the integer part of ``a``. This should produce, for example,
```math
    \begin{equation}
    \label{eq:permutation}
    (a, b, c, d, e, f, g, h) \quad \rightarrow \quad (a, c, e, g, h, f, d, b) \, ,
    \end{equation}
```
after which \eqref{eq:FCT} is computed using
```math
  Y = \text{DCT}(X) = 2 \text{Re} \left\lbrace \omega_{4N}^k \text{FFT} \lbrace X^\prime \rbrace \right\rbrace \, .
```

### 1D fast inverse cosine transform
The inverse \eqref{eq:IFCT} can be computed using
```math
  Y = \text{IDCT}(X) = \text{Re} \left\lbrace \omega_{4N}^{-k} \text{IFFT} \lbrace X \rbrace \right\rbrace \, ,
```
after which the inverse permutation of \eqref{eq:permutation} must be applied.

### 2D fast cosine transform
Unfortunately, the 1D algorithm cannot be applied dimension-wise so the 2D algorithm is more 
complicated. Thankfully, the permutation \eqref{eq:permutation} can be applied dimension-wise. 
The forward cosine transform for a real signal of length ``N_1 \times N_2`` is then given by
```math
Y_{k_1, k_2} = \text{DCT}(X_{n_1, n_2}) =
2 \text{Re} \left\lbrace
    \omega_{4N_1}^k \left( \omega_{4N_2}^k \tilde{X} + \omega_{4N_2}^{-k} \tilde{X}^- \right)
\right\rbrace \, ,
```
where ``\tilde{X} = \text{FFT}(X^\prime)`` and ``\tilde{X}^-`` indicates that ``\tilde{X}`` is indexed in reverse.

### 2D fast inverse cosine transform
The inverse can be computed using
```math
Y_{k_1, k_2} = \text{IDCT}(X_{n_1, n_2}) =
\frac{1}{4} \text{Re} \left\lbrace
    \omega_{4N_1}^{-k} \omega_{4N_2}^{-k}
    \left( \tilde{X} - M_1 M_2 \tilde{X}^{--} \right)
    - \mathrm{i} \left( M_1 \tilde{X}^{-+} + M_2 \tilde{X}^{+-} \right)
\right\rbrace \, ,
```
where ``\tilde{X} = \text{IFFT}(X)`` here, ``\tilde{X}^{-+}`` is indexed in reverse along the first dimension,
``\tilde{X}^{-+}`` along the second dimension, and ``\tilde{X}^{--}`` along both. ``M_1`` and ``M_2`` are masks of lengths
``N_1`` and ``N_2`` respectively, both containing ones except at the first element where ``M_0 = 0``. Afterwards, the inverse
permutation of \eqref{eq:permutation} must be applied.

Due to the extra steps involved in calculating the cosine transform in 2D, running with two 
wall-bounded dimensions typically slows the model down by a factor of 2. Switching to the FACR 
algorithm may help here as a 2D cosine transform won't be necessary anymore.

## Iterative Solvers

For problems with irregular grids the eigenvectors of the discrete Poisson operator are no longer simple Fourier
series sines and cosines. This means discrete Fast Fourier Transforms can't be used to generate the projection 
of the equation right hand side onto eigenvectors. So an eigenvector based approach to solving
the Poisson equation is not computationally efficient.

An pre-conditioned conjugate gradient iterative solver is used instead for problems with grids
that are non uniform in multiple directions. This includes curvilinear grids on the sphere and
also telescoping cartesian grids that stretch along more than one dimension. There are two forms 
of the pressure operator in this approach. One is rigid lid form and one is an implicit 
free-surface form.

### Rigid lid pressure operator

The rigid lid operator is based on the same continuous form as is used in the Direct Method
solver.

### Implicit free surface pressure operator

The implicit free surface solver solves for the free-surface, ``\eta``, in the vertically
integrated continuity equation

```math
    \begin{equation}
    \label{eq:vertically-integrated-continuity}
    \partial_t \eta + \partial_x (H \hat{u}) + \partial_y (H \hat{v}) = M \, ,
    \end{equation}
```

where ``M`` is some surface volume flux (e.g., terms such as precipitation, evaporation and runoff); 
currently ``M=0`` is assumed. To form a linear system that can be solved implicitly we recast
the continuity equation into a discrete integral form

```math
    \begin{equation}
    \label{eq:semi-discrete-integral-continuity}
    A_z \partial_t \eta + \delta_{x}^{caa} \sum_{k} A_{x} u + \delta_y^{caa} \sum_k A_y v = A_z M \, ,
    \end{equation}
```

and apply the discrete form to the hydrostatic form of the velocity fractional step equation

```math
    \begin{equation}
    \label{eq:hydrostatic-fractional-step}
    \boldsymbol{v}^{n+1} = \boldsymbol{v}^{\star} - g \Delta t \boldsymbol{\nabla} \eta^{n+1} \, .
    \end{equation}
```

as follows.

Assuming ``M=0`` (for now), for the ``n+1`` timestep velocity we want the following to hold

```math
    A_z \frac{\eta^{n+1} - \eta^{n}}{\Delta t} = -\delta_x^{caa} \sum_k A_x u^{n+1} - \delta_y^{caa} \sum_k A_y v^{n+1} \, .
```

Substituting ``u^{n+1}`` and ``v^{n+1}`` from the discrete form of the  right-hand-side of
\eqref{eq:hydrostatic-fractional-step} then gives an implicit equation for ``\eta^{n+1}``,

```math
\begin{align}
   \delta_x^{caa}\sum_k A_x \partial_x^{faa} \eta^{n+1} & + \delta_y^{aca} \sum_k A_y \partial_y^{afa} \eta^{n+1} - \frac{1}{g \, \Delta t^2} A_z \eta^{n+1} = \nonumber \\
   & = \frac{1}{g \, \Delta t} \left( \delta_x^{caa} \sum_k A_x u^{\star} + \delta_y^{aca} \sum_k A_y v^{\star} \right ) - \frac{1}{g \, \Delta t^{2}} A_z \eta^{n} \, .
\end{align}
```

Formulated in this way, the linear operator will be symmetric and so can be solved using a
preconditioned conjugate gradient algorithmn.
# Spatial operators

To calculate the various terms and perform the time-stepping, discrete difference and interpolation 
operators must be designed from which all the terms, such as momentum advection and Laplacian 
diffusion, may be constructed. Much of the material in this section is derived from [Marshall97FV](@cite).

## Differences

Difference operators act as the discrete form of the derivative operator. Care must be taken 
when calculating differences on a staggered grid. For example, the the difference of a cell-centered 
variable such as temperature ``T`` lies on the faces  in the direction of the difference, and 
vice versa. In principle, there are three difference operators, one for each  direction
```math
  \delta_x f = f_E - f_W , \quad
  \delta_y f = f_N - f_S , \quad
  \delta_z f = f_T - f_B ,
```
where the ``E`` and ``W`` subscripts indicate that the value is evaluated the eastern or western 
wall of the cell, ``N`` and ``S`` indicate the northern and southern walls, and ``T`` and ``B`` 
indicate the top and bottom walls.

Additionally, two ``\delta`` operators must be defined for each direction to account for the 
staggered nature of the grid. One for taking the difference of a cell-centered variable and 
projecting it onto the cell faces
```math
\begin{align}
    \delta_x^{faa} f_{i, j, k} &= f_{i, j, k} - f_{i-1, j, k} \, , \\
    \delta_y^{afa} f_{i, j, k} &= f_{i, j, k} - f_{i, j-1, k} \, , \\
    \delta_z^{aaf} f_{i, j, k} &= f_{i, j, k} - f_{i, j, k-1} \, , 
\end{align}
```
and another for taking the difference of a face-centered variable and projecting it onto the
cell centers
```math
\begin{align}
    \delta_x^{caa} f_{i, j, k} &= f_{i+1, j, k} - f_{i, j, k} \, , \\
    \delta_y^{aca} f_{i, j, k} &= f_{i, j+1, k} - f_{i, j, k} \, , \\
    \delta_z^{aac} f_{i, j, k} &= f_{i, j, k+1} - f_{i, j, k} \, .
\end{align}
```

## Interpolation

In order to add or multiply variables that are defined at different points they are interpolated. 
In our case, linear interpolation or averaging is employed. Once again, there are two averaging 
operators, one for each direction,
```math
\begin{equation}
  \overline{f}^x = \frac{f_E + f_W}{2} \, , \quad
  \overline{f}^y = \frac{f_N + f_S}{2} \, , \quad
  \overline{f}^z = \frac{f_T + f_B}{2} \, .
\end{equation}
```

Additionally, three averaging operators must be defined for each direction. One for taking the 
average of a cell-centered  variable and projecting it onto the cell faces
```math
\begin{align}
    \overline{f_{i, j, k}}^{faa} = \frac{f_{i, j, k} + f_{i-1, j, k}}{2} \, , \\
    \overline{f_{i, j, k}}^{afa} = \frac{f_{i, j, k} + f_{i, j-1, k}}{2} \, , \\
    \overline{f_{i, j, k}}^{aaf} = \frac{f_{i, j, k} + f_{i, j, k-1}}{2} \, ,
\end{align}
```
and another for taking the average of a face-centered variable and projecting it onto the cell centers
```math
\begin{align}
    \overline{f_{i, j, k}}^{caa} = \frac{f_{i+1, j, k} + f_{i, j, k}}{2} \, , \\
    \overline{f_{i, j, k}}^{aca} = \frac{f_{i, j+1, k} + f_{i, j, k}}{2} \, , \\
    \overline{f_{i, j, k}}^{aac} = \frac{f_{i, j, k+1} + f_{i, j, k}}{2} \, .
\end{align}
```

## Divergence and flux divergence

The divergence of the flux of a cell-centered quantity over the cell can be calculated as
```math
\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{f}
= \frac{1}{V} \left[ \delta_x^{faa} (A_x f_x)
                   + \delta_y^{afa} (A_y f_y)
                   + \delta_z^{aaf} (A_z f_z) \right] \, ,
```
where ``\boldsymbol{f} = (f_x, f_y, f_z)`` is the flux with components defined normal to the 
faces, and ``V`` is the volume of the cell. The presence of a solid boundary is indicated by 
setting the appropriate flux normal to the boundary to zero.

A similar divergence operator can be defined for a face-centered quantity. The divergence of,
e.g., the flux of ``T`` over a cell, ``\boldsymbol{\nabla} \boldsymbol{\cdot} (\boldsymbol{v} T)``, 
is then
```math
\renewcommand{\div}[1] {\boldsymbol{\nabla} \boldsymbol{\cdot} \left ( #1 \right )}
\div{\boldsymbol{v} T}
= \frac{1}{V} \left[ \delta_x^{caa} (A_x u \overline{T}^{faa})
                   + \delta_y^{aca} (A_y v \overline{T}^{afa})
                   + \delta_z^{aac} (A_z w \overline{T}^{aaf}) \right] \, ,
```
where ``T`` is interpolated onto the cell faces where it can be multiplied by the velocities, 
which are then differenced and projected onto the cell centers where they added together.

## Momentum advection

The advection terms that appear in model equations can be rewritten using the incompressibility 
(``\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{v} = 0``) as, e.g,
```math
\renewcommand{\div}[1] {\boldsymbol{\nabla} \boldsymbol{\cdot} \left ( #1 \right )}
\begin{align}
\boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} u & = \div{u \boldsymbol{v}} - u ( \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{v} ) \nonumber \\
    & = \div{u \boldsymbol{v}} \, ,
\end{align}
```
which can then be discretized similarly to the flux divergence operator, however, they must 
be discretized differently for each direction.

For example, the ``x``-momentum advection operator is discretized as
```math
\boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} u
= \frac{1}{\overline{V}^x} \left[
    \delta_x^{faa} \left( \overline{A_x u}^{caa} \overline{u}^{caa} \right)
  + \delta_y^{afa} \left( \overline{A_y v}^{aca} \overline{u}^{aca} \right)
  + \delta_z^{aaf} \left( \overline{A_z w}^{aac} \overline{u}^{aac} \right)
\right] \, ,
```
where ``\overline{V}^x`` is the average of the volumes of the cells on either side of the face 
in question. Calculating ``\partial_x (uu)`` can be performed by interpolating ``A_x u`` and 
``u`` onto the cell centers then multiplying them and differencing them back onto the faces. 
However, in the case of the the two other terms, ``\partial_y (vu)`` and ``\partial_z (wu)``, 
the two variables must be interpolated onto the cell edges to be multiplied then differenced 
back onto the cell faces.

## Discretization of isotropic diffusion operators

An isotropic viscosity operator acting on vertical momentum is discretized via
```math
    \boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \nu \boldsymbol{\nabla} w \right )
    = \frac{1}{V} \left[
          \delta_x^{faa} ( \nu \overline{A_x}^{caa} \partial_x^{caa} w )
        + \delta_y^{afa} ( \nu \overline{A_y}^{aca} \partial_y^{aca} w )
        + \delta_z^{aaf} ( \nu \overline{A_z}^{aac} \partial_z^{aac} w )
    \right ] \, ,
```
where ``\nu`` is the kinematic viscosity.

An isotropic diffusion operator acting on a tracer ``c``, on the other hand, is discretized via
```math
   \boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \kappa \boldsymbol{\nabla} c \right )
    = \frac{1}{V} \left[ \vphantom{\overline{A_x}^{caa}}
        \delta_x^{caa} ( \kappa A_x \partial_x^{faa} c )
      + \delta_y^{aca} ( \kappa A_y \partial_y^{afa} c )
      + \delta_z^{aac} ( \kappa A_z \partial_z^{aaf} c )
    \right] \, .
```

## Vertical integrals
Vertical integrals are converted into sums along each column. For example, the hydrostatic pressure 
anomaly is
```math
    p_{HY}^\prime = \int_{-L_z}^0 b^\prime \, \mathrm{d} z \, ,
```
where ``b^\prime`` is the buoyancy perturbation. Converting it into a sum that we compute from 
the top downwards we get
```math
    \begin{equation}
    p_{HY}^\prime(k) =
        \begin{cases}
            - \overline{b_{N_z}^\prime}^{aaf} \Delta z^F_{N_z},               & \quad k = N_z \, , \\
            p_{HY}^\prime(k+1) - \overline{b_{k+1}^\prime}^{aaf} \Delta z^F_k, & \quad 1 \le k \le N_z - 1 \, ,
        \end{cases}
    \end{equation}
```
where we converted the sum into a recursive definition for ``p_{HY}^\prime(k)`` in terms of 
``p_{HY}^\prime(k+1)`` so that the integral may be computed with ``\mathcal{O}(N_z)`` operations 
by a single thread.

The vertical velocity ``w`` may be computed from ``u`` and ``v`` via the continuity equation
```math
    w = - \int_{-L_z}^0 (\partial_x u + \partial_y v) \, \mathrm{d} z \, ,
```
to satisfy the incompressibility condition ``\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{v} = 0``
to numerical precision. This also involves computing a vertical integral, in this case evaluated
from the bottom up
```math
    \begin{equation}
    w_k =
        \begin{cases}
            0, & \quad k = 1 \, , \\
            w_{k-1} - \left( \partial_x^{caa} u + \partial_y^{aca} v \right) \Delta z^C_k, & \quad 2 \le k \le N_z \, .
        \end{cases}
    \end{equation}
```
# Pressure decomposition

In the numerical implementation of the momentum equations, the kinematic pressure ``p`` 
is split into "hydrostatic anomaly" and "non-hydrostatic" parts via
```math
    \begin{equation}
    \label{eq:pressure}
    p(\boldsymbol{x}, t) = p_{\rm{hyd}}(\boldsymbol{x}, t) + p_{\rm{non}}(\boldsymbol{x}, t) \, .
    \end{equation}
```
The anomalous hydrostatic component of the kinematic pressure is defined by 
```math
    \begin{align}
    \label{eq:hydrostaticpressure}
    \partial_z p_{\rm{hyd}} \equiv -b \, ,
    \end{align}
```
such that the sum of the kinematic pressure and buoyancy perturbation becomes
```math
    \begin{align}
    -\boldsymbol{\nabla} p + b \boldsymbol{\hat z} = 
        - \boldsymbol{\nabla} p_{\rm{non}}
        - \boldsymbol{\nabla}_h p_{\rm{hyd}} \, ,
    \end{align}
```
where ``\boldsymbol{\nabla}_h \equiv \partial_x \boldsymbol{\hat x} + \partial_y \boldsymbol{\hat y}`` 
is the horizontal gradient. The hydrostatic pressure anomaly is so named because the "total" 
hydrostatic pressure contains additional components:
```math
\begin{align}
\partial_z p_{\text{total hydrostatic}} & = - g \left ( 1 + \frac{\rho_*}{\rho_0} + \frac{\rho'}{\rho_0} \right ) \, , \\
                                           & = \partial_z p_{\rm{hyd}} - g \left ( 1 + \frac{\rho_*}{\rho_0} \right ) \, .
\end{align}
```
Under this pressure decomposition the pressure gradient that appears in the momentum equations becomes
```math
   \boldsymbol{\nabla} p \mapsto \boldsymbol{\nabla} p_{\rm{non}} + \boldsymbol{\nabla}_h p_{\rm{hyd}}\, .
```
Mathematically, the non-hydrostatic pressure ``p_{\rm{non}}`` enforces the incompressibility constraint.
# Finite volume method on a staggered grid

The `Oceananigans.jl` staggered grid is defined by a rectilinear array of cuboids of horizontal dimensions 
``\Delta x_{i, j, k}, \Delta y_{i, j, k}`` and vertical dimension 
``\Delta z_{i, j, k}``, where ``(i, j, k)`` index the location of each cell in the staggered grid.
Note that the indices ``(i, j, k)`` increase with increasing coordinate ``(x, y, z)``.

![Schematic of staggered grid](assets/staggered_grid.png)
*A schematic of `Oceananigans.jl` finite volumes for a two-dimensional staggered grid in ``(x, z)``.
Tracers ``c`` and pressure ``p`` are defined at the center of the control volume. The ``u`` control volumes are 
centered on the left and right edges of the pressure control volume while the ``w`` control volumes are centered 
on the top and bottom edges of the pressure control volumes. The indexing convention places the ``i^{\rm{th}}`` 
``u``-node on cell ``x``-faces to the left of the ``i`` tracer point at cell centers.*

Dropping explicit indexing, the areas of cell faces are given by
```math
    A_x = \Delta y \Delta z, \quad A_y = \Delta x \Delta z, \quad A_z = \Delta x \Delta y \, ,
```
so that each cell encloses a volume ``V = \Delta x \Delta y \Delta z``.

A finite volume method discretizes a continuous quantity ``c`` by considering its average over a finite volume:
```math
    c_{i, j, k} \equiv \frac{1}{V_{i, j, k}} \int c(\boldsymbol{x}) \, \mathrm{d} V_{i, j, k} \, .
```
The finite volumes that discretize each of ``u``, ``v``, and ``w`` are located on a grid which is "staggered" 
with respect to the grid that defines tracer finite volumes. 
The nodes, or central points of the velocity finite volumes are co-located with the faces of the tracer 
finite volume.
In particular, the ``u``-nodes are located in the center of the "``x``-face" (east of the tracer point), 
``v``-nodes are located on ``y``-faces south of the tracer point, and ``w``-nodes are located on 
``z``-faces downwards from the tracer point.
# [Numerical implementation of boundary conditions](@id numerical_bcs)

We adopt a mixed approach for implementing boundary conditions that uses both halo regions and "direct"
imposition of boundary conditions, depending on the condition prescribed.

We illustrate how boundary conditions are implemented by considering the tracer equation \eqref{eq:tracer}.

See [Model setup: boundary conditions](@ref model_step_bcs) for how to create and use these
boundary conditions in Oceananigans.

## Gradient boundary conditions

Users impose gradient boundary conditions by prescribing the gradient ``\gamma`` of a field 
``c`` across an *external boundary* ``\partial \Omega_b``. The prescribed gradient ``\gamma`` 
may be a constant, discrete array of values, or an arbitrary function. The gradient boundary 
condition is enforced setting the value of halo points located outside the domain interior 
such that
```math
    \begin{equation}
    \label{eq:gradient-bc}
    \hat{\boldsymbol{n}} \boldsymbol{\cdot} \boldsymbol{\nabla} c |_{\partial \Omega_b} = \gamma \, .
    \end{equation}
```
where ``\hat{\boldsymbol{n}}`` is the vector normal to ``\partial \Omega_b``.

Across the bottom boundary in ``z``, for example, this requires that
```math
    \begin{equation}
    \label{eq:linear-extrapolation}
    c_{i, j, 0} = c_{i, j, 1} + \gamma_{i, j, 1} \tfrac{1}{2} \left ( \Delta z_{i, j, 1} + \Delta z_{i, j, 0} \right ) \, ,
    \end{equation}
```
where ``\Delta z_{i, j, 1} = \Delta z_{i, j, 0}`` are the heights of the finite volume at ``i, j`` and ``k=1`` and ``k=0``.
This prescription implies that the ``z``-derivative of ``c`` across the boundary at ``k=1`` is
```math
    \begin{equation}
    \partial_z c \, |_{i, j, 1} \equiv
        \frac{c_{i, j, 1} - c_{i, j, 0}}{\tfrac{1}{2} \left ( \Delta z_{i, j, 1} + \Delta z_{i, j, 0} \right )}
            = \gamma_{i, j, 1} \, ,
    \end{equation}
```
as prescribed by the user.

Gradient boundary conditions are represented by the [`Gradient`](@ref) type.

## Value boundary conditions

Users impose value boundary conditions by prescribing ``c^b``, the value of ``c`` on the external
boundary ``\partial \Omega_b``.
The value ``c^b`` may be a constant, array of discrete values, or an arbitrary function.
To enforce a value boundary condition, the gradient associated with the difference between
``c^b`` and ``c`` at boundary-adjacent nodes is diagnosed and used to set the value of the ``c`` halo point
located outside the boundary.

At the bottom boundary in ``z``, for example, this means that the gradient of ``c`` is determined by
```math
    \begin{equation}
    \gamma = \frac{c_{i, j, 1} - c^b_{i, j, 1}}{\tfrac{1}{2} \Delta z_{i, j, 1}} \, ,
    \end{equation}
```
which is then used to set the halo point ``c_{i, j, 0}`` via linear extrapolation.

Value boundary conditions are represented by the [`Value`](@ref) type.

## Flux boundary conditions

Users impose flux boundary conditions by prescribing the flux ``q_c \, |_b`` of ``c`` across
the external boundary ``\partial \Omega_b``. The flux ``q_c \, |_b`` may be a constant, array 
of discrete values, or arbitrary function. To explain how flux boundary conditions are imposed 
in `Oceananigans.jl`, we note that the average of the tracer conservation equation over a finite 
volume yields
```math
    \begin{equation}
    \label{eq:dc/dt}
    \partial_t c_{i, j, k} = - \frac{1}{V_{i, j, k}} \oint_{\partial \Omega_{i, j, k}} (\boldsymbol{v} c + \boldsymbol{q}_c) \, \mathrm{d} S
                             + \frac{1}{V_{i, j, k}} \int_{V_{i, j, k}} F_c \, \mathrm{d} V \, ,
    \end{equation}
```
where the surface integral over ``\partial \Omega_{i, j, k}`` averages the flux of ``c`` across 
the six faces of the finite volume. The right-hand-side of \eqref{eq:dc/dt} above is denoted as 
``G_c |_{i, j, k}``.


An external boundary of a finite volume is associated with a no-penetration condition such that
``\hat{\boldsymbol{n}} \boldsymbol{\cdot} \boldsymbol{v} \, |_{\partial \Omega_b} = 0``, where 
``\hat{\boldsymbol{n}}`` is the vector normal to ``\partial \Omega_b``. Furthermore, the closures 
currently available in `Oceananigans.jl` have the property that ``\boldsymbol{q}_c \propto \boldsymbol{\nabla} c``.
Thus setting ``\hat{\boldsymbol{n}} \boldsymbol{\cdot} \boldsymbol{\nabla} c \, |_{\partial \Omega_b} = 0`` 
on the external boundary implies that the total flux of ``c`` across the external boundary is
```math
    \begin{equation}
    \hat{\boldsymbol{n}} \boldsymbol{\cdot} \left ( \boldsymbol{v} c + \boldsymbol{q}_c \right ) |_{\partial \Omega_b} = 0 \, .
    \end{equation}
```
`Oceananigans.jl` exploits this fact to define algorithm that prescribe fluxes across external 
boundaries ``\partial \Omega_b``:

1. Impose a constant gradient ``\hat{\boldsymbol{n}} \boldsymbol{\cdot} \boldsymbol{\nabla} c 
   \, |_{\partial \Omega_b} = 0`` across external boundaries via using halo points (similar 
   to \eqref{eq:gradient-bc}), which ensures that the evaluation of ``G_c`` in boundary-adjacent
   cells does not include fluxes across the external boundary, and;
2. Add the prescribed flux to the boundary-adjacent volumes prior to calculating ``G_c``: 
   ``G_c \, |_b = G_c \, |_b - \frac{A_b}{V_b} q_c \, |_b \, \text{sign}(\hat{\boldsymbol{n}})``, 
   where ``G_c \, |_b`` denotes values of ``G_c`` in boundary-adjacent volumes, ``q_c \, |_b`` 
   is the flux prescribed along the boundary, ``V_b`` is the volume of the boundary-adjacent 
   cell, and ``A_b`` is the area of the external boundary of the boundary-adjacent cell.

   The factor ``\text{sign}(\hat{\boldsymbol{n}})`` is ``-``1 and ``+``1 on "left" and "right" 
   boundaries, and accounts for the fact that a positive flux on a left boundary where 
   ``\text{sign}(\hat{\boldsymbol{n}}) = -1`` implies an "inward" flux of ``c`` that increases 
   interior values of ``c``, whereas a positive flux on a right boundary where 
   ``\text{sign}(\hat{\boldsymbol{n}}) = 1`` implies an "outward" flux that decreases interior
   values of ``c``.

Flux boundary conditions are represented by the [`Flux`](@ref) type.
# [Turbulence closures](@id numerical_closures)

To truly simulate and resolve turbulence at high Reynolds number (so basically all interesting flows) would require
you resolve all motions down to the [Kolmogorov41](@cite) length scale ``\eta = (\nu^3 / \varepsilon)^{1/4}`` where
``\nu`` is the kinematic viscosity and ``\varepsilon`` the average rate of dissipation of turbulence kinetic energy per
unit mass.

As pointed out way back by [Corrsin61](@cite), to run a simulation on a horizontal domain about 10 times the size of an
"average eddy" with 100 vertical levels and where the grid spacing is given by ``\eta`` would require the computer to
store on the order of ``10^{14}`` variables.[^1] This is still impractical today, although may be within
reach in less than a decade. He ends by suggesting the use of an analog rather digital computer---a tank of water.

[^1]: And even then, ``\eta`` gives the *maximum* allowable grid spacing. There is significant flow structure
    smaller than ``\eta``.

To have any hope of simulating high Reynolds number flows we need some way of resolving the sub-grid scale motions.[^2]

[^2]: In reality there is no need to resolve all motions down to the Kolmogorov length scale to achieve
    acceptable accuracy. Perhaps good results can be achieved if 80\% of the kinetic energy is resolved
    (§13) [Pope00](@cite).


## Reynolds-averaged Navier–Stokes equations

Following [Reynolds1895](@cite) we can decompose flow variables such as velocity ``\boldsymbol{v}`` into the mean component
``\overline{\boldsymbol{v}}`` and the fluctuating component ``\boldsymbol{v}^\prime`` so that ``\boldsymbol{v} = \overline{\boldsymbol{v}} + \boldsymbol{v}^\prime``
[see §4 of [Pope00](@cite) for a modern discussion].

Expressing the Navier-Stokes equations in tensor notation
```math
\begin{align}
    \partial_i v_i &= 0  \, ,\\
    \partial_t v_i + v_j \partial_j v_i &= f_i - \alpha\partial_i p + \nu \partial_j \partial_j v_i \, ,
\end{align}
```
where ``\alpha = \rho^{-1}`` is the specific volume and ``f_i`` represents external forces. We can plug in the Reynolds
decomposition for ``\boldsymbol{v}`` and after some manipulation arrive at the following form for the *Reynolds-averaged
Navier-Stokes equations*
```math
\begin{align}
    \partial_i \overline{u}_i &= 0  \, ,\\
    \partial_t \overline{u}_i + \overline{u}_j \partial_j \overline{u}_i &= \overline{f}_i -
    \partial_j \left(-\alpha\overline{p}\delta_{ij} + 2\nu \overline{S}_{ij} - \overline{v_i^\prime v_j^\prime}\right) \, ,
\end{align}
```
where
```math
\overline{S}_{ij} = \frac{1}{2} ( \partial_j \overline{u}_i + \partial_i \overline{u}_j ) \, ,
```
is the mean rate of strain tensor.

Thanks to the non-linearity of the Navier-Stokes equations, even when averaged we are left with pesky fluctuation
terms which form the components of the *Reynolds stress tensor*
```math
\tau_{ij} = \rho \overline{v_i^\prime v_j^\prime} \, .
```
Attempting to close the equations leads to the *closure problem*: the time evolution of the Reynolds stresses
depends on  triple covariances ``\overline{v_i^\prime v_j^\prime v_k^\prime}`` and covariances with pressure, which depend
on quadruple covariances and so on [Chou45](@cite).

This is kind of hopeless so we will have to find some way to model the Reynolds stresses.

## Gradient-diffusion hypothesis and eddy viscosity models

The *gradient-diffusion hypothesis*, due to [Boussinesq1877](@cite), assumes that the transport of scalar fluxes
such as ``\overline{\boldsymbol{v}^\prime c^\prime}`` and ``\overline{v_i^\prime v_j^\prime}`` occurs down the mean scalar gradient
``\grad c`` as if they are being diffused (§4.4) [Pope00](@cite). This is in analogy with how momentum transfer by
molecular motion in a gas can be described by a molecular viscosity.

Taking this assumption we can express the Reynolds stresses and turbulent tracer fluxes in terms of the mean variables
and close the equations
```math
\overline{\boldsymbol{v}^\prime c^\prime} = -\kappa_e \boldsymbol{\nabla} \overline{c}
\quad \text{and} \quad
\overline{v_i^\prime v_j^\prime} = -2\nu_e \overline{S}_{ij} \, ,
```
where ``\nu_e = \nu_e(\boldsymbol{x}, t)`` is the turbulent or *eddy viscosity* and ``\kappa_e = \kappa_e(\boldsymbol{x}, t)``
is the *eddy diffusivity*.

The effective diffusivity ends up being the sum of the molecular and eddy diffusivities. So just by using an elevated
value for the viscosity and diffusivity, you are already using an eddy viscosity model.

The eddy viscosity model is simple and for that reason is very popular. It can work well even with a constant eddy
diffusivity. However, it does assume that the flux is aligned down gradient, which is not true even in simple turbulent
flows as the physics of turbulence is quite different from that of colliding molecules leading to the viscous stress law
(§4.4,10.1) [Pope00](@cite). So we might want something a little bit more sophisticated.
# The Boussinesq approximation

Oceananigans.jl employ often the Boussinesq approximation[^1]. In the Boussinesq approximation
the fluid density ``\rho`` is, in general, decomposed into three components:
```math
    \rho(\boldsymbol{x}, t) = \rho_0 + \rho_*(z) + \rho'(\boldsymbol{x}, t) \, ,
```
where ``\rho_0`` is a constant 'reference' density, ``\rho_*(z)`` is a background density
profile which, when non-zero, is typically associated with the hydrostatic compression
of seawater in the deep ocean, and ``\rho'(\boldsymbol{x}, t)`` is the dynamic component of density
corresponding to inhomogeneous distributions of a buoyant tracer such as temperature or salinity.

The fluid *buoyancy*, associated with the buoyant acceleration of fluid, is
defined in terms of ``\rho'`` as
```math
    b = - \frac{g \rho'}{\rho_0} \, ,
```
where ``g`` is gravitational acceleration.

The Boussinesq approximation is valid when ``\rho_* + \rho' \ll \rho_0``, which implies the
fluid is _approximately_ incompressible, and thus does not support acoustic waves. In this case, 
the mass conservation equation reduces to the continuity equation
```math
    \begin{equation}
    \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{v} = \partial_x u + \partial_y v + \partial_z w = 0 \, .
    \label{eq:continuity}
    \end{equation}
```

[^1]: Named after Boussinesq (1903) although used earlier by Oberbeck (1879), the Boussinesq
      approximation neglects density differences in the momentum equation except when associated
      with the gravitational term. It is an accurate approximation for many flows, and especially
      so for oceanic flows where density differences are very small. See Vallis (2017, section 2.4)
      for an oceanographic introduction to the Boussinesq equations and Vallis (2017, Section 2.A)
      for an asymptotic derivation. See Kundu (2015, Section 4.9) for an engineering
      introduction.

# Coriolis forces

The Coriolis model controls the manifestation of the term ``\boldsymbol{f} \times \boldsymbol{v}``
in the momentum equation.

## ``f``-plane approximation

Under an ``f``-plane approximation[^3] the reference frame in which
the momentum and tracer equations are solved rotates at a constant rate.

### The traditional ``f``-plane approximation

In the *traditional* ``f``-plane approximation, the coordinate system rotates around
a vertical axis such that
```math
    \boldsymbol{f} = f \boldsymbol{\hat z} \, ,
```
where ``f`` is constant and determined by the user.

## The arbitrary-axis constant-Coriolis approximation

In this approximation, the coordinate system rotates around an axis in the ``x,y,z``-plane, such
that
```math
    \boldsymbol{f} = f_x \boldsymbol{\hat x} + f_y \boldsymbol{\hat y} + f_z \boldsymbol{\hat z} \, ,
```
where ``f_x``, ``f_y`` and ``f_z`` are constants determined by the user.


[^3]: The ``f``-plane approximation is used to model the effects of Earth's rotation on anisotropic 
      fluid motion in a plane tangent to the Earth's surface. In this case, the projection of 
      the Earth's rotation vector at latitude ``\varphi`` and onto a coordinate system in which 
      ``x, y, z`` correspond to the directions east, north, and up is
      ``\boldsymbol{f} \approx \frac{4 \pi}{\text{day}} \left ( \cos \varphi \boldsymbol{\hat y} + \sin \varphi \boldsymbol{\hat z} \right ) \, ,``
      where the Earth's rotation rate is approximately ``2 \pi / \text{day}``. The *traditional* 
      ``f``-plane approximation neglects the ``y``-component of this projection, which is appropriate 
      for fluid motions with large horizontal-to-vertical aspect ratios.

## ``\beta``-plane approximation

### The traditional ``\beta``-plane approximation

Under the *traditional* ``\beta``-plane approximation, the rotation axis is vertical as for the
``f``-plane approximation, but ``f`` is expanded in a Taylor series around a central latitude 
such that
```math
    \boldsymbol{f} = \left ( f_0 + \beta y \right ) \boldsymbol{\hat z} \, ,
```
where ``f_0`` is the planetary vorticity at some central latitude, and ``\beta`` is the
planetary vorticity gradient.
The ``\beta``-plane model is not periodic in ``y`` and thus can be used only in domains that
are bounded in the ``y``-direction.

### The non-traditional ``\beta``-plane approximation

The *non-traditional* ``\beta``-plane approximation accounts for the latitudinal variation of both
the locally vertical and the locally horizontal components of the rotation vector
```math
    \boldsymbol{f} = \left[ 2\Omega\cos\varphi_0 \left( 1 -  \frac{z}{R} \right) + \gamma y \right] \boldsymbol{\hat y}
           + \left[ 2\Omega\sin\varphi_0 \left( 1 + 2\frac{z}{R} \right) + \beta  y \right] \boldsymbol{\hat z} \, ,
```
as can be found in the paper by [DellarJFM2011](@cite) where 
``\beta = 2 \Omega \cos \varphi_0 / R`` and ``\gamma = -4 \Omega \sin \varphi_0 / R``.
# Surface gravity waves and the Craik-Leibovich approximation

In Oceananiagns.jl, users model the effects of surface waves by specifying spatial and
temporal gradients of the Stokes drift velocity field.
At the moment, only uniform unidirectional Stokes drift fields are supported, in which case
```math
    \boldsymbol{u}^S = u^S(z, t) \hat{\boldsymbol{x}} + v^S(z, t) \hat{\boldsymbol{y}} \, .
```
Surface waves are modeled in Oceananigans.jl by the Craik-Leibovich approximation,
which governs interior motions under a surface gravity wave field that have been time- or
phase-averaged over the rapid oscillations of the surface waves.
The oscillatory vertical and horizontal motions associated with surface waves themselves,
therefore, are not present in the resolved velocity field ``\boldsymbol{v}``, and only the 
steady, averaged effect of surface waves that manifests over several or more wave oscillations 
are modeled.

In Oceananigans.jl with surface waves, the resolved velocity field ``\boldsymbol{v}`` is the 
Lagrangian-mean velocity field. The Lagrangian-mean velocity field at a particular location 
``(x, y, z)`` is average velocity of a fluid particle whose average position is ``(x, y, z)`` 
at time ``t``. The average position of a fluid particle ``\boldsymbol{\xi}(t) = (\xi, \eta, \zeta)`` 
is thus governed by
```math
    \partial_t \boldsymbol{\xi} + \boldsymbol{v}(\boldsymbol{\xi}, t) \boldsymbol{\cdot} \boldsymbol{\nabla} \boldsymbol{\xi} = \boldsymbol{v}(\boldsymbol{\xi}, t) \, ,
```
which is the same relationship that holds when surface waves are not present and ``\boldsymbol{v}`` 
ceases to be an averaged velocity field. The simplicity of the governing equations for Lagrangian-mean 
momentum is the main reason we use a Lagrangian-mean formulation in Oceananigans.jl, rather 
than an Eulerian-mean formulation: for example, the tracer conservation equation is unchanged 
by the inclusion of surface wave effects. Moreover, because the effect of surface waves manifests 
either as a bulk forcing of Lagrangian-mean momentum or as a modification to the effective background 
rotation rate of the interior fluid similar to any bulk forcing or Coriolis force, we do not 
explicitly include the effects of surface waves in turbulence closures that model the effects 
of subgrid turbulence. More specifically, the effect of steady surface waves does not effect 
the conservation of Lagrangian-mean turbulent kinetic energy.

The Lagrangian-mean velocity field ``\boldsymbol{v}`` contrasts with the Eulerian-mean velocity 
field ``\boldsymbol{v}^E``, which is the fluid velocity averaged at the fixed Eulerian position 
``(x, y, z)``. The surface wave Stokes drift field supplied by the user is, in fact, defined
by the difference between the Eulerian- and Lagrangian-mean velocity:
```math
    \boldsymbol{u}^S \equiv \boldsymbol{v} - \boldsymbol{v}^E \, .
```
The Stokes drift velocity field is typically prescribed for idealized scenarios, or determined
from a wave model for the evolution of surface waves under time-dependent atmospheric winds
in more realistic cases.
# Shallow water model

The `ShallowWaterModel` solves the shallow water dynamics for a fluid of constant density but 
with varying fluid depth ``h(x, y, t)``. The dynamics for the evolution of the two-dimensional 
flow ``\boldsymbol{u}(x, y, t) = u(x, y, t) \boldsymbol{\hat x} + v(x, y, t) \boldsymbol{\hat y}`` 
and the fluid's height ``h(x, y, t)`` is:
```math
  \begin{align}
    \partial_t \boldsymbol{u} + \boldsymbol{u} \boldsymbol{\cdot} \boldsymbol{\nabla} \boldsymbol{u} 
    + \boldsymbol{f} \times \boldsymbol{u} & = - g \boldsymbol{\nabla} h \, ,\\
    \partial_t h + \boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{u} h \right ) & = 0 \, .
  \end{align}
```

Using the transport along each direction ``\boldsymbol{u} h`` as our dynamical 
variables, we can express the shallow-water dynamics in conservative form:
```math
  \begin{align}
    \partial_t (\boldsymbol{u} h) + \boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{u} \boldsymbol{u} h \right ) + \boldsymbol{f} \times (\boldsymbol{u} h) & = - g \boldsymbol{\nabla} \left ( \frac1{2} h^2 \right ) \, ,\\
    \partial_t h + \boldsymbol{\nabla} \boldsymbol{\cdot} (\boldsymbol{u} h) & = 0 \, ,
  \end{align}
```
where ``\boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{u} \boldsymbol{u} h \right )`` 
denotes a vector whose components are ``[\boldsymbol{\nabla} \boldsymbol{\cdot} (\boldsymbol{u} \boldsymbol{u} h)]_i = \boldsymbol{\nabla} \boldsymbol{\cdot} (u_i \boldsymbol{u} h)``.

The `ShallowWaterModel` state variables are the transports, `uh` and `vh` and the fluid's 
height `h`. We can retrieve the flow velocities by dividing the corresponding transport by 
the fluid's height, e.g., `v = vh / h`.# Buoyancy model and equations of state

The buoyancy model determines the relationship between tracers and the buoyancy ``b`` in the momentum equation.

## Buoyancy tracer

The simplest buoyancy model uses buoyancy ``b`` itself as a tracer: ``b`` obeys the tracer
conservation equation and is used directly in the momentum equations in the momentum equation.

## Seawater buoyancy

For seawater buoyancy is, in general, modeled as a function of conservative temperature
``T``, absolute salinity ``S``, and depth below the ocean surface ``d`` via
```math
    \begin{equation}
    b = - \frac{g}{\rho_0} \rho' \left (T, S, d \right ) \, ,
    \label{eq:seawater-buoyancy}
    \end{equation}
```
where ``g`` is gravitational acceleration, ``\rho_0`` is the reference density.
The function ``\rho'(T, S, d)`` in the seawater buoyancy relationship that links conservative temperature,
salinity, and depth to the density perturbation is called the *equation of state*.
Both ``T`` and ``S`` obey the tracer conservation equation.

### Linear equation of state

Buoyancy is determined from a linear equation of state via
```math
    b = g \left ( \alpha T - \beta S \right ) \, ,
```
where ``g`` is gravitational acceleration, ``\alpha`` is the thermal expansion coefficient,
and ``\beta`` is the haline contraction coefficient.

### Nonlinear equation of state

Buoyancy is determined by the simplified equations of state introduced by [Roquet15TEOS](@cite).
# Boundary conditions

In Oceananigans.jl the user may impose \textit{no-penetration}, \textit{flux},
\textit{gradient} (Neumann), and \textit{value} (Dirichlet) boundary conditions in bounded,
non-periodic directions.
Note that the only boundary condition available for a velocity field normal to the bounded
direction is \textit{no-penetration}.

## Flux boundary conditions

A flux boundary condition prescribes flux of a quantity normal to the boundary.
  For a tracer ``c`` this corresponds to prescribing
```math
q_c \, |_b \equiv \boldsymbol{q}_c \boldsymbol{\cdot} \hat{\boldsymbol{n}} \, |_{\partial \Omega_b} \, ,
```
where ``\partial \Omega_b`` is an external boundary.

## Gradient (Neumann) boundary condition

A gradient boundary condition prescribes the gradient of a field normal to the boundary.
For a tracer ``c`` this prescribes
```math
\gamma \equiv \boldsymbol{\nabla} c \boldsymbol{\cdot} \hat{\boldsymbol{n}} \, |_{\partial \Omega_b} \, .
```

## Value (Dirichlet) boundary condition

A value boundary condition prescribes the value of a field on a boundary; for a tracer this
prescribes
```math
c_b \equiv c \, |_{\partial \Omega_b} \, .
```

## No penetration boundary condition

A no penetration boundary condition prescribes the velocity component normal to a boundary to be 0,
so that
```math
\boldsymbol{\hat{n}} \boldsymbol{\cdot} \boldsymbol{v} \, |_{\partial \Omega_b} = 0 \, .
```
# Turbulence closures

The turbulence closure selected by the user determines the form of stress divergence
``\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{\tau}`` and diffusive flux divergence
``\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{q}_c`` in the momentum and tracer conservation equations.

## Constant isotropic diffusivity

In a constant isotropic diffusivity model, the kinematic stress tensor is defined
```math
\tau_{ij} = - \nu \Sigma_{ij} \, ,
```
where ``\nu`` is a constant viscosity and
``\Sigma_{ij} \equiv \tfrac{1}{2} \left ( v_{i, j} + v_{j, i} \right )`` is the strain-rate
tensor. The divergence of ``\boldsymbol{\tau}`` is then
```math
\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{\tau} = -\nu \nabla^2 \boldsymbol{v} \, .
```
Similarly, the diffusive tracer flux is ``\boldsymbol{q}_c = - \kappa \boldsymbol{\nabla} c`` for tracer
diffusivity ``\kappa``, and the diffusive tracer flux divergence is
```math
\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{q}_c = - \kappa \nabla^2 c \, .
```
Each tracer may have a unique diffusivity ``\kappa``.

## Constant anisotropic diffusivity

In Oceananigans.jl, a constant anisotropic diffusivity implies a constant tensor
diffusivity ``\nu_{j k}`` and stress ``\boldsymbol{\tau}_{ij} = \nu_{j k} u_{i, k}`` with non-zero
components ``\nu_{11} = \nu_{22} = \nu_h`` and ``\nu_{33} = \nu_v``.
With this form the kinematic stress divergence becomes
```math
\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{\tau} = - \left [ \nu_h \left ( \partial_x^2 + \partial_y^2 \right )
                                    + \nu_v \partial_z^2 \right ] \boldsymbol{v} \, ,
```
and diffusive flux divergence
```math
\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{q}_c = - \left [ \kappa_{h} \left ( \partial_x^2 + \partial_y^2 \right )
                                    + \kappa_{v} \partial_z^2 \right ] c \, ,
```
in terms of the horizontal viscosities and diffusivities ``\nu_h`` and ``\kappa_{h}`` and the
vertical viscosity and diffusivities ``\nu_v`` and ``\kappa_{v}``.
Each tracer may have a unique diffusivity components ``\kappa_h`` and ``\kappa_v``.

## Constant anisotropic biharmonic diffusivity

In Oceananigans.jl, a constant anisotropic biharmonic diffusivity implies a constant tensor
diffusivity ``\nu_{j k}`` and stress ``\boldsymbol{\tau}_{ij} = \nu_{j k} \partial_k^3 u_i`` with non-zero
components ``\nu_{11} = \nu_{22} = \nu_h`` and ``\nu_{33} = \nu_v``.
With this form the kinematic stress divergence becomes
```math
\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{\tau} = - \left [ \nu_h \left ( \partial_x^2 + \partial_y^2 \right )^2
                                    + \nu_v \partial_z^4 \right ] \boldsymbol{v} \, ,
```
and diffusive flux divergence
```math
\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{q}_c = - \left [ \kappa_{h} \left ( \partial_x^2 + \partial_y^2 \right )^2
                                    + \kappa_{v} \partial_z^4 \right ] c \, ,
```
in terms of the horizontal biharmonic viscosities and diffusivities ``\nu_h`` and ``\kappa_{h}`` and the
vertical biharmonic viscosity and diffusivities ``\nu_v`` and ``\kappa_{v}``.
Each tracer may have a unique diffusivity components ``\kappa_h`` and ``\kappa_v``.

## Smagorinsky-Lilly turbulence closure

In the turbulence closure proposed by [Lilly62](@cite) and [Smagorinsky63](@cite),
the subgrid stress associated with unresolved turbulent motions is modeled diffusively via
```math
\tau_{ij} = \nu_e \Sigma_{ij} \, ,
```
where ``\Sigma_{ij} = \tfrac{1}{2} \left ( v_{i, j} + v_{j, i} \right )`` is the resolved
strain rate.
The eddy viscosity is given by
```math
    \begin{align}
    \nu_e = \left ( C \Delta_f \right )^2 \sqrt{ \Sigma^2 } \, \varsigma(N^2 / \Sigma^2) + \nu \, ,
    \label{eq:smagorinsky-viscosity}
    \end{align}
```
where ``\Delta_f`` is the "filter width" associated with the finite volume grid spacing,
``C`` is a user-specified model constant, ``\Sigma^2 \equiv \Sigma_{ij} \Sigma_{ij}``, and
``\nu`` is a constant isotropic background viscosity.
The factor ``\varsigma(N^2 / \Sigma^2)`` reduces ``\nu_e`` in regions of
strong stratification via
```math
    \varsigma(N^2 / \Sigma^2) = \sqrt{1 - \min \left ( 1, C_b N^2 / \Sigma^2 \right )} \, ,
```
where ``N^2 = \max \left (0, \partial_z b \right )`` is the squared buoyancy frequency for stable
stratification with ``\partial_z b > 0`` and ``C_b`` is a user-specified constant.  Lilly (1962)
proposed ``C_b = 1/Pr``, where ``Pr`` is a turbulent Prandtl number.
The filter width for the Smagorinsky-Lilly closure is
```math
\Delta_f(\boldsymbol{x}) = \left ( \Delta x \Delta y \Delta z \right)^{1/3} \, ,
```
where ``\Delta x``, ``\Delta y``, and ``\Delta z`` are the grid spacing in the
``\boldsymbol{\hat x}``, ``\boldsymbol{\hat y}``, and ``\boldsymbol{\hat z}`` directions at location ``\boldsymbol{x} = (x, y, z)``.

The effect of subgrid turbulence on tracer mixing is also modeled diffusively via
```math
\boldsymbol{q}_c = \kappa_e \boldsymbol{\nabla} c \, ,
```
where the eddy diffusivity ``\kappa_e`` is
```math
\kappa_e = \frac{\nu_e - \nu}{Pr} + \kappa \, ,
```
where ``\kappa`` is a constant isotropic background diffusivity.
Both ``Pr`` and ``\kappa`` may be set independently for each tracer.

## Anisotropic minimum dissipation (AMD) turbulence closure

Oceananigans.jl uses the anisotropic minimum dissipation (AMD) model proposed by
[Verstappen18](@cite) and described and tested by [Vreugdenhil18](@cite).
The AMD model uses an eddy diffusivity hypothesis similar the Smagorinsky-Lilly model.
In the AMD model, the eddy viscosity and diffusivity for each tracer are defined in terms
of eddy viscosity and diffusivity *predictors*
``\nu_e^\dagger`` and ``\kappa_e^\dagger``, such that
```math
    \nu_e = \max \left ( 0, \nu_e^\dagger \right ) + \nu
    \quad \text{and} \quad
    \kappa_e = \max \left ( 0, \kappa_e^\dagger \right ) + \kappa \, ,
```
to ensure that ``\nu_e \ge 0`` and ``\kappa_e \ge 0``, where ``\nu`` and ``\kappa`` are the
constant isotropic background viscosity and diffusivities for each tracer. The eddy viscosity 
predictor is
```math
    \begin{equation}
    \nu_e^\dagger = -C \Delta_f^2
    \frac
        {(\hat{\partial}_k \hat{v}_i) (\hat{\partial}_k \hat{v}_j) \hat{\Sigma}_{ij}
        + C_b \hat{\delta}_{i3} (\hat{\partial}_k \hat{v_i}) (\hat{\partial}_k b)}
        {(\hat{\partial}_l \hat{v}_m) (\hat{\partial}_l \hat{v}_m)} \, ,
    \label{eq:nu-dagger}
    \end{equation}
```
while the eddy diffusivity predictor for tracer ``c`` is
```math
    \begin{equation}
    \label{eq:kappa-dagger}
    \kappa_e^\dagger = -C \Delta_f^2
    \frac
        {(\hat{\partial}_k \hat{v}_i) (\hat{\partial}_k c) (\hat{\partial}_i c)}
        {(\hat{\partial}_l c) (\hat{\partial}_l c)} \, .
    \end{equation}
```
In the definitions of the eddy viscosity and eddy diffusivity predictor, ``C`` and ``C_b`` are
user-specified model constants, ``\Delta_f`` is a "filter width" associated with the finite volume
grid spacing, and the hat decorators on partial derivatives, velocities, and the Kronecker
delta ``\hat \delta_{i3}`` are defined such that
```math
    \hat \partial_i \equiv \Delta_i \partial_i, \qquad
    \hat{v}_i(x, t) \equiv \frac{v_i(x, t)}{\Delta_i}, \quad \text{and} \quad
    \hat{\delta}_{i3} \equiv \frac{\delta_{i3}}{\Delta_3} \, .
```
A velocity gradient, for example, is therefore
``\hat{\partial}_i \hat{v}_j(x, t) = \frac{\Delta_i}{\Delta_j} \partial_i v_j(x, t)``,
while the normalized strain tensor is
```math
    \hat{\Sigma}_{ij} =
        \frac{1}{2} \left[ \hat{\partial}_i \hat{v}_j(x, t) + \hat{\partial}_j \hat{v}_i(x, t) \right] \, .
```
The filter width ``\Delta_f`` in that appears in the viscosity and diffusivity predictors
is taken as the square root of the harmonic mean of the squares of the filter widths in
each direction:
```math
    \frac{1}{\Delta_f^2} = \frac{1}{3} \left(   \frac{1}{\Delta x^2}
                                              + \frac{1}{\Delta y^2}
                                              + \frac{1}{\Delta z^2} \right) \, .
```
The constant ``C_b`` permits the "buoyancy modification" term it multiplies to be omitted
from a calculation.
By default we use the model constants ``C=1/12`` and ``C_b=0``.

## Convective adjustment vertical diffusivity

This closure aims to model the enhanced mixing that occurs due to convection.
At every point and for every time instance, the closure diagnoses the gravitational stability of the fluid and applies the vertical diffusivities (i) `background_νz` to `u, v` and `background_κz` to all tracers if the fluid is gravitationally neutral or stable with `∂z(b) >= 0`, or (ii) `convective_νz` and `convective_κz` if `∂z(b) >= 0`.
This closure is a plausible model for convection if `convective_κz` ``\gg`` `background_κz` and `convective_νz` ``\gg`` `background_νz`.
# Hydrostatic model with a free surface

The `HydrostaticFreeSurfaceModel` solves the incompressible Navier-Stokes equations under the 
Boussinesq and hydrostatic approximations and with an arbitrary number of tracer conservation 
equations. Physics associated with individual terms in the momentum and tracer conservation
equations --- the background rotation rate of the equation's reference frame,
gravitational effects associated with buoyant tracers under the Boussinesq
approximation, generalized stresses and tracer fluxes associated with viscous and
diffusive physics, and arbitrary "forcing functions" --- are determined by the whims of the
user.

## Mass conservation and free surface evolution equation

The mass conservation equation is
```math
    0 = \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{u} + \partial_z w \, .  
```

The above is integrated from the bottom of the fluid up to ``z = 0`` to obtain ``w(x, y, z, t)``.

The free surface displacement ``\eta(x, y, t)`` satisfies the linearized kinematic boundary 
condition at the surface
```math
    \partial_t \eta = w(x, y, z=0, t) \, .
```

## The momentum conservation equation

The equations governing the conservation of momentum in a rotating fluid, including buoyancy
via the Boussinesq approximation are
```math
    \begin{align}
    \partial_t \boldsymbol{u} & = - \left ( \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} \right ) \boldsymbol{u}
                        - \boldsymbol{f} \times \boldsymbol{u} 
                        - \boldsymbol{\nabla}_h (p + g \eta)
                        - \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{\tau}
                        + \boldsymbol{F_u} \, , \label{eq:momentum}\\
                        0 & = b - \partial_z p \, , \label{eq:hydrostatic}
    \end{align}
```
where ``b`` the is buoyancy, ``\boldsymbol{\tau}`` is the hydrostatic kinematic stress tensor, 
``\boldsymbol{F_u}`` denotes an internal forcing of the velocity field ``\boldsymbol{u}``, 
``p`` is kinematic pressure, ``\eta`` is the free-surface displacement, and ``\boldsymbol{f}`` 
is the *Coriolis parameter*, or the background vorticity associated with the specified rate of 
rotation of the frame of reference.

Equation \eqref{eq:hydrostatic} above is the hydrostatic approximation and comes about as the 
dominant balance of terms in the Navier-Stokes vertical momentum equation under the Boussinesq 
approximation.

The terms that appear on the right-hand side of the momentum conservation equation are (in order):

* momentum advection: ``\left ( \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} \right ) 
  \boldsymbol{u}``,
* Coriolis: ``\boldsymbol{f} \times \boldsymbol{u}``,
* baroclinic kinematic pressure gradient: ``\boldsymbol{\nabla} p``,
* barotropic kinematic pressure gradient: ``\boldsymbol{\nabla} (g \eta)``,
* molecular or turbulence viscous stress: ``\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{\tau}``, and
* an arbitrary internal source of momentum: ``\boldsymbol{F_u}``.

## The tracer conservation equation

The conservation law for tracers in Oceananigans.jl is
```math
    \begin{align}
    \partial_t c = - \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} c
                   - \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{q}_c
                   + F_c \, ,
    \label{eq:tracer}
    \end{align}
```
where ``\boldsymbol{q}_c`` is the diffusive flux of ``c`` and ``F_c`` is an arbitrary source term.
Oceananigans.jl permits arbitrary tracers and thus an arbitrary number of tracer equations to 
be solved simultaneously with the momentum equations.

From left to right, the terms that appear on the right-hand side of the tracer conservation 
equation are

* tracer advection: ``\boldsymbol{u} \boldsymbol{\cdot} \boldsymbol{\nabla} c``,
* molecular or turbulent diffusion: ``\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{q}_c``, and
* an arbitrary internal source of tracer: ``F_c``.

The following subsections provide more details on the possible forms that each individual term 
in the momentum and tracer equations can take in Oceananigans.jl.
# Coordinate system and notation

Oceananigans.jl is formulated in a Cartesian coordinate system ``\boldsymbol{x} = (x, y, z)`` 
with unit vectors ``\boldsymbol{\hat x}``, ``\boldsymbol{\hat y}``, and ``\boldsymbol{\hat z}``, 
where ``\boldsymbol{\hat x}`` points east, ``\boldsymbol{\hat y}`` points north, and ``\boldsymbol{\hat z}`` 
points 'upward', opposite the direction of gravitational acceleration.

We denote time with ``t``, partial derivatives with respect to time ``t`` or a coordinate ``x`` 
with ``\partial_t`` or ``\partial_x``, and denote the gradient operator ``\boldsymbol{\nabla} \equiv 
\partial_x \boldsymbol{\hat x} + \partial_y \boldsymbol{\hat y} + \partial_z \boldsymbol{\hat z}``. 
Horizontal gradients are denoted with ``\boldsymbol{\nabla}_h \equiv \partial_x \boldsymbol{\hat x} + \partial_y \boldsymbol{\hat y}``.

We use ``u``, ``v``, and ``w`` to denote the east, north, and vertical velocity components,
such that ``\boldsymbol{v} = u \boldsymbol{\hat x} + v \boldsymbol{\hat y} + w \boldsymbol{\hat z}``.
We reserve ``\boldsymbol{v}`` for the three-dimensional velocity field and use ``\boldsymbol{u}``
to denote the horizontal components of flow, i.e., ``\boldsymbol{u} = u \boldsymbol{\hat x} + 
v \boldsymbol{\hat y}``.# Nonhydrostatic model

The `NonhydrostaticModel` solves the incompressible Navier-Stokes equations under the Boussinesq
approximation and an arbitrary number of tracer conservation equations.
Physics associated with individual terms in the momentum and tracer conservation
equations --- the background rotation rate of the equation's reference frame,
gravitational effects associated with buoyant tracers under the Boussinesq
approximation, generalized stresses and tracer fluxes associated with viscous and
diffusive physics, and arbitrary "forcing functions" --- are determined by the whims of the
user.

## The momentum conservation equation

The equations governing the conservation of momentum in a rotating fluid, including buoyancy
via the Boussinesq approximation and including the averaged effects of surface gravity waves
at the top of the domain via the Craik-Leibovich approximation are
```math
    \begin{align}
    \partial_t \boldsymbol{v} & = - \left ( \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} \right ) \boldsymbol{v}
                        - \left ( \boldsymbol{V} \boldsymbol{\cdot} \boldsymbol{\nabla} \right ) \boldsymbol{v}
                        - \left ( \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} \right ) \boldsymbol{V} \nonumber \\
                        & \qquad
                        - \left ( \boldsymbol{f} - \boldsymbol{\nabla} \times \boldsymbol{u}^S \right ) \times \boldsymbol{v} 
                        - \boldsymbol{\nabla} p
                        + b \boldsymbol{\hat g}
                        - \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{\tau}
                        + \partial_t \boldsymbol{u}^S
                        + \boldsymbol{F_v} \, ,
    \label{eq:momentum}
    \end{align}
```
where ``b \boldsymbol{\hat g}`` the is the buoyancy (a vector whose default direction is upward), 
``\boldsymbol{\tau}`` is the kinematic stress tensor, ``\boldsymbol{F_v}``
denotes an internal forcing of the velocity field ``\boldsymbol{v}``, ``p`` is the kinematic 
pressure, ``\boldsymbol{u}^S`` is the horizontal, two-dimensional 'Stokes drift' velocity field associated with surface gravity 
waves, and ``\boldsymbol{f}`` is the Coriolis parameter, or the background vorticity associated 
with the specified rate of rotation of the frame of reference.

The terms that appear on the right-hand side of the momentum conservation equation are (in order):

* momentum advection: ``\left ( \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} \right ) 
  \boldsymbol{v}``,
* advection of resolved momentum by the background velocity field ``\boldsymbol{V}``: 
  ``\left ( \boldsymbol{V} \boldsymbol{\cdot} \boldsymbol{\nabla} \right ) \boldsymbol{v}``,
* advection of background momentum by resolved velocity: ``\left ( \boldsymbol{v} \boldsymbol{\cdot} 
  \boldsymbol{\nabla} \right ) \boldsymbol{V}``,
* Coriolis: ``\boldsymbol{f} \times \boldsymbol{v}``,
* the effective background rotation rate due to surface waves: ``\left ( \boldsymbol{\nabla} \times 
  \boldsymbol{u}^S \right ) \times \boldsymbol{v}``,
* kinematic pressure gradient: ``\boldsymbol{\nabla} p``,
* buoyant acceleration: ``b``,
* vertical unit vector (pointing to the direction opposite to gravity): ``\boldsymbol{\hat g}``,
* molecular or turbulence viscous stress: ``\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{\tau}``,
* a source of momentum due to forcing or damping of surface waves: ``\partial_t \boldsymbol{u}^S``, and
* an arbitrary internal source of momentum: ``\boldsymbol{F_v}``.

## The tracer conservation equation

The conservation law for tracers in Oceananigans.jl is
```math
    \begin{align}
    \partial_t c = - \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} c
                   - \boldsymbol{V} \boldsymbol{\cdot} \boldsymbol{\nabla} c
                   - \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} C
                   - \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{q}_c
                   + F_c \, ,
    \label{eq:tracer}
    \end{align}
```
where ``\boldsymbol{q}_c`` is the diffusive flux of ``c`` and ``F_c`` is an arbitrary source term.
Oceananigans.jl permits arbitrary tracers and thus an arbitrary number of tracer equations to 
be solved simultaneously with the momentum equations.

From left to right, the terms that appear on the right-hand side of the tracer conservation equation are

* tracer advection: ``\boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} c``,
* tracer advection by the background velocity field, ``\boldsymbol{V}``: ``\boldsymbol{V} \boldsymbol{\cdot} \boldsymbol{\nabla} c``,
* advection of the background tracer field, ``C``, by the resolved velocity field: ``\boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} C``,
* molecular or turbulent diffusion: ``\boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{q}_c``, and
* an arbitrary internal source of tracer: ``F_c``.

The following subsections provide more details on the possible forms that each individual term 
in the momentum and tracer equations can take in Oceananigans.jl.
# Output writers

`AbstractOutputWriter`s save data to disk.
`Oceananigans` provides three ways to write output:

1. [`NetCDFOutputWriter`](@ref) for output of arrays and scalars that uses [NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl)
2. [`JLD2OutputWriter`](@ref) for arbitrary julia data structures that uses [JLD2.jl](https://github.com/JuliaIO/JLD2.jl)
3. [`Checkpointer`](@ref) that automatically saves as much model data as possible, using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl)

The `Checkpointer` is discussed on a separate documentation page.

## Basic usage

[`NetCDFOutputWriter`](@ref) and [`JLD2OutputWriter`](@ref) require four inputs:

1. The `model` from which output data is sourced (required to initialize the `OutputWriter`).
2. A key-value pairing of output "names" and "output" objects. `JLD2OutputWriter` accepts `NamedTuple`s and `Dict`s;
   `NetCDFOutputWriter` accepts `Dict`s with string-valued keys. Output objects are either `AbstractField`s or
   functions that return data when called via `func(model)`.
3. A `schedule` on which output is written. `TimeInterval`, `IterationInterval`, `WallTimeInterval` schedule
   periodic output according to the simulation time, simulation interval, or "wall time" (the physical time
   according to a clock on your wall). A fourth `schedule` called `AveragedTimeInterval` specifies
   periodic output that is time-averaged over a `window` prior to being written.
4. The filename and directory. Currently `NetCDFOutputWriter` accepts one `filepath` argument, while
   `JLD2OutputWriter` accepts a filename `prefix` and `dir`ectory.

Other important keyword arguments are

* `field_slicer::FieldSlicer` for outputting subregions, two- and one-dimensional slices of fields.
  By default a `FieldSlicer` is used to remove halo regions from fields so that only the physical
  portion of model data is saved to disk.

* `array_type` for specifying the type of the array that holds outputted field data. The default is
  `Array{Float32}`, or arrays of single-precision floating point numbers.

Once an `OutputWriter` is created, it can be used to write output by adding it the
ordered dictionary `simulation.output_writers`. prior to calling `run!(simulation)`.

More specific detail about the `NetCDFOutputWriter` and `JLD2OutputWriter` is given below.

!!! tip "Time step alignment and output writing"
    Oceananigans simulations will shorten the time step as needed to align model output with each
    output writer's schedule.

## NetCDF output writer

Model data can be saved to NetCDF files along with associated metadata. The NetCDF output writer is generally used by
passing it a dictionary of (label, field) pairs and any indices for slicing if you don't want to save the full 3D field.

### Examples

Saving the u velocity field and temperature fields, the full 3D fields and surface 2D slices
to separate NetCDF files:

```jldoctest netcdf1
using Oceananigans

grid = RectilinearGrid(size=(16, 16, 16), extent=(1, 1, 1));

model = NonhydrostaticModel(grid=grid, tracers=:c);

simulation = Simulation(model, Δt=12, stop_time=3600);

fields = Dict("u" => model.velocities.u, "c" => model.tracers.c);

simulation.output_writers[:field_writer] =
    NetCDFOutputWriter(model, fields, filepath="more_fields.nc", schedule=TimeInterval(60))

# output
NetCDFOutputWriter scheduled on TimeInterval(1 minute):
├── filepath: more_fields.nc
├── dimensions: zC(16), zF(17), xC(16), yF(16), xF(16), yC(16), time(0)
├── 2 outputs: ["c", "u"]
├── field slicer: FieldSlicer(:, :, :, with_halos=false)
└── array type: Array{Float32}
```

```jldoctest netcdf1
simulation.output_writers[:surface_slice_writer] =
    NetCDFOutputWriter(model, fields, filepath="another_surface_xy_slice.nc",
                       schedule=TimeInterval(60), field_slicer=FieldSlicer(k=grid.Nz))

# output
NetCDFOutputWriter scheduled on TimeInterval(1 minute):
├── filepath: another_surface_xy_slice.nc
├── dimensions: zC(1), zF(1), xC(16), yF(16), xF(16), yC(16), time(0)
├── 2 outputs: ["c", "u"]
├── field slicer: FieldSlicer(:, :, 16, with_halos=false)
└── array type: Array{Float32}
```

```jldoctest netcdf1
simulation.output_writers[:averaged_profile_writer] =
    NetCDFOutputWriter(model, fields,
                       filepath = "another_averaged_z_profile.nc",
                       schedule = AveragedTimeInterval(60, window=20),
                       field_slicer = FieldSlicer(i=1, j=1))

# output
NetCDFOutputWriter scheduled on TimeInterval(1 minute):
├── filepath: another_averaged_z_profile.nc
├── dimensions: zC(16), zF(17), xC(1), yF(1), xF(1), yC(1), time(0)
├── 2 outputs: ["c", "u"] averaged on AveragedTimeInterval(window=20 seconds, stride=1, interval=1 minute)
├── field slicer: FieldSlicer(1, 1, :, with_halos=false)
└── array type: Array{Float32}
```

`NetCDFOutputWriter` also accepts output functions that write scalars and arrays to disk,
provided that their `dimensions` are provided:

```jldoctest
using Oceananigans

grid = RectilinearGrid(size=(16, 16, 16), extent=(1, 2, 3));

model = NonhydrostaticModel(grid=grid);

simulation = Simulation(model, Δt=1.25, stop_iteration=3);

f(model) = model.clock.time^2; # scalar output
g(model) = model.clock.time .* exp.(znodes(Center, grid)); # vector/profile output
h(model) = model.clock.time .* (   sin.(xnodes(Center, grid, reshape=true)[:, :, 1])
                            .*     cos.(ynodes(Face, grid, reshape=true)[:, :, 1])); # xy slice output

outputs = Dict("scalar" => f, "profile" => g, "slice" => h);

dims = Dict("scalar" => (), "profile" => ("zC",), "slice" => ("xC", "yC"));

output_attributes = Dict(
    "scalar"  => Dict("longname" => "Some scalar", "units" => "bananas"),
    "profile" => Dict("longname" => "Some vertical profile", "units" => "watermelons"),
    "slice"   => Dict("longname" => "Some slice", "units" => "mushrooms")
);

global_attributes = Dict("location" => "Bay of Fundy", "onions" => 7);

simulation.output_writers[:things] =
    NetCDFOutputWriter(model, outputs,
                       schedule=IterationInterval(1), filepath="some_things.nc", dimensions=dims, verbose=true,
                       global_attributes=global_attributes, output_attributes=output_attributes)

# output
NetCDFOutputWriter scheduled on IterationInterval(1):
├── filepath: some_things.nc
├── dimensions: zC(16), zF(17), xC(16), yF(16), xF(16), yC(16), time(0)
├── 3 outputs: ["profile", "slice", "scalar"]
├── field slicer: FieldSlicer(:, :, :, with_halos=false)
└── array type: Array{Float32}
```

See [`NetCDFOutputWriter`](@ref) for more information.

## JLD2 output writer

JLD2 is a fast HDF5 compatible file format written in pure Julia.
JLD2 files can be opened in Julia with the [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) package
and in Python with the [h5py](https://www.h5py.org/) package.

The `JLD2OutputWriter` receives either a `Dict`ionary or `NamedTuple` containing
`name, output` pairs. The `name` can be a symbol or string. The `output` must either be
an `AbstractField` or a function called with `func(model)` that returns arbitrary output.
Whenever output needs to be written, the functions will be called and the output
of the function will be saved to the JLD2 file.

### Examples

Write out 3D fields for u, v, w, and a tracer c, along with a horizontal average:

```jldoctest jld2_output_writer
using Oceananigans
using Oceananigans.Utils: hour, minute

model = NonhydrostaticModel(grid=RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1)), tracers=(:c,))
simulation = Simulation(model, Δt=12, stop_time=1hour)

function init_save_some_metadata!(file, model)
    file["author"] = "Chim Riggles"
    file["parameters/coriolis_parameter"] = 1e-4
    file["parameters/density"] = 1027
    return nothing
end

c_avg = Field(Average(model.tracers.c, dims=(1, 2)))

# Note that model.velocities is NamedTuple
simulation.output_writers[:velocities] = JLD2OutputWriter(model, model.velocities,
                                                          prefix = "some_more_data",
                                                          schedule = TimeInterval(20minute),
                                                          init = init_save_some_metadata!)

# output
JLD2OutputWriter scheduled on TimeInterval(20 minutes):
├── filepath: ./some_more_data.jld2
├── 3 outputs: (:u, :v, :w)
├── field slicer: FieldSlicer(:, :, :, with_halos=false)
├── array type: Array{Float32}
├── including: [:grid, :coriolis, :buoyancy, :closure]
└── max filesize: Inf YiB
```

and a time- and horizontal-average of tracer `c` every 20 minutes of simulation time
to a file called `some_more_averaged_data.jld2`

```jldoctest jld2_output_writer
simulation.output_writers[:avg_c] = JLD2OutputWriter(model, (; c=c_avg),
                                                     prefix = "some_more_averaged_data",
                                                     schedule = AveragedTimeInterval(20minute, window=5minute))

# output
JLD2OutputWriter scheduled on TimeInterval(20 minutes):
├── filepath: ./some_more_averaged_data.jld2
├── 1 outputs: (:c,) averaged on AveragedTimeInterval(window=5 minutes, stride=1, interval=20 minutes)
├── field slicer: FieldSlicer(:, :, :, with_halos=false)
├── array type: Array{Float32}
├── including: [:grid, :coriolis, :buoyancy, :closure]
└── max filesize: Inf YiB
```


See [`JLD2OutputWriter`](@ref) for more information.

## Time-averaged output

Time-averaged output is specified by setting the `schedule` keyword argument for either `NetCDFOutputWriter` or
`JLD2OutputWriter` to [`AveragedTimeInterval`](@ref).

With `AveragedTimeInterval`, the time-average of ``a`` is taken as a left Riemann sum corresponding to

```math
\langle a \rangle = \frac{1}{T} \int_{t_i-T}^{t_i} a \, \mathrm{d} t \, ,
```

where ``\langle a \rangle`` is the time-average of ``a``, ``T`` is the time-`window` for averaging specified by
the `window` keyword argument to `AveragedTimeInterval`, and the ``t_i`` are discrete times separated by the
time `interval`. The ``t_i`` specify both the end of the averaging window and the time at which output is written.

### Example

Building an `AveragedTimeInterval` that averages over a 1 year window, every 4 years,

```jldoctest averaged_time_interval
using Oceananigans.OutputWriters: AveragedTimeInterval
using Oceananigans.Utils: year, years

schedule = AveragedTimeInterval(4years, window=1year)

# output
AveragedTimeInterval(window=1 year, stride=1, interval=4 years)
```

An `AveragedTimeInterval` schedule directs an output writer
to time-average its outputs before writing them to disk:

```jldoctest averaged_time_interval
using Oceananigans
using Oceananigans.OutputWriters: JLD2OutputWriter
using Oceananigans.Utils: minutes

model = NonhydrostaticModel(grid=RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1)))

simulation = Simulation(model, Δt=10minutes, stop_time=30years)

simulation.output_writers[:velocities] = JLD2OutputWriter(model, model.velocities,
                                                          prefix = "even_more_averaged_velocity_data",
                                                          schedule = AveragedTimeInterval(4years, window=1year, stride=2))

# output
JLD2OutputWriter scheduled on TimeInterval(4 years):
├── filepath: ./even_more_averaged_velocity_data.jld2
├── 3 outputs: (:u, :v, :w) averaged on AveragedTimeInterval(window=1 year, stride=2, interval=4 years)
├── field slicer: FieldSlicer(:, :, :, with_halos=false)
├── array type: Array{Float32}
├── including: [:grid, :coriolis, :buoyancy, :closure]
└── max filesize: Inf YiB
```
# Buoyancy models and equations of state

The buoyancy option selects how buoyancy is treated in `NonhydrostaticModel`s and
`HydrostaticFreeSurfaceModel`s (`ShallowWaterModel`s do not have that option given the physics of
the model). There are currently three alternatives:

1. No buoyancy (and no gravity).
2. Evolve buoyancy as a tracer.
3. _Seawater buoyancy_: evolve temperature ``T`` and salinity ``S`` as tracers with a value for the gravitational
   acceleration ``g`` and an equation of state of your choosing.

## No buoyancy

To turn off buoyancy (and gravity) you can simply pass `buoyancy = nothing` to the model
constructor. For example to create a `NonhydrostaticModel`:


```@meta
DocTestSetup = quote
    using Oceananigans
end
```


```jldoctest buoyancy
julia> grid = RectilinearGrid(size=(64, 64, 64), extent=(1, 1, 1));

julia> model = NonhydrostaticModel(grid=grid, buoyancy=nothing)
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: ()
├── closure: Nothing
├── buoyancy: Nothing
└── coriolis: Nothing
```

`buoyancy=nothing` is the default option for`NonhydrostaticModel`, so ommitting `buoyancy`
from the `NonhydrostaticModel` constructor yields an identical result:

```jldoctest buoyancy
julia> model = NonhydrostaticModel(grid=grid)
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: ()
├── closure: Nothing
├── buoyancy: Nothing
└── coriolis: Nothing
```

To create a `HydrostaticFreeSurfaceModel` without a buoyancy term we explicitly
specify `buoyancy=nothing` flag. The default tracers `T` and `S` for `HydrostaticFreeSurfaceModel`
may be eliminated when `buoyancy=nothing` by specifying `tracers=()`:

```jldoctest buoyancy
julia> model = HydrostaticFreeSurfaceModel(grid=grid, buoyancy=nothing, tracers=())
HydrostaticFreeSurfaceModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: ()
├── closure: Nothing
├── buoyancy: Nothing
├── free surface: ExplicitFreeSurface with gravitational acceleration 9.80665 m s⁻²
└── coriolis: Nothing
```

## Buoyancy as a tracer

Both `NonhydrostaticModel` and `HydrostaticFreeSurfaceModel` support evolving
a buoyancy tracer by including `:b` in `tracers` and specifying  `buoyancy = BuoyancyTracer()`:

```jldoctest buoyancy
julia> model = NonhydrostaticModel(grid=grid, buoyancy=BuoyancyTracer(), tracers=:b)
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: (:b,)
├── closure: Nothing
├── buoyancy: Buoyancy{BuoyancyTracer, Oceananigans.Grids.ZDirection}
└── coriolis: Nothing
```

We follow the same pattern to create a `HydrostaticFreeSurfaceModel` with buoyancy as a tracer:

```jldoctest buoyancy
julia> model = HydrostaticFreeSurfaceModel(grid=grid, buoyancy=BuoyancyTracer(), tracers=:b)
HydrostaticFreeSurfaceModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: (:b,)
├── closure: Nothing
├── buoyancy: Buoyancy{BuoyancyTracer, Oceananigans.Grids.ZDirection}
├── free surface: ExplicitFreeSurface with gravitational acceleration 9.80665 m s⁻²
└── coriolis: Nothing
```

## Seawater buoyancy

`NonhydrostaticModel` and `HydrostaticFreeSurfaceModel` support modeling the buoyancy of seawater
as a function of gravitational acceleration, conservative temperature ``T`` and absolute salinity ``S``.
The relationship between ``T``, ``S``, the geopotential height, and the density perturbation from
a reference value is called the `equation_of_state`.
Specifying `buoyancy = SeawaterBuoyancy()` (which uses a linear equation of state and
[Earth standard](https://en.wikipedia.org/wiki/Standard_gravity)
`gravitational_acceleration = 9.80665 \, \text{m}\,\text{s}^{-2}` by default)
requires the tracers `:T` and `:S`:

```jldoctest buoyancy
julia> model = NonhydrostaticModel(grid=grid, buoyancy=SeawaterBuoyancy(), tracers=(:T, :S))
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: (:T, :S)
├── closure: Nothing
├── buoyancy: Buoyancy{SeawaterBuoyancy{Float64, LinearEquationOfState{Float64}, Nothing, Nothing}, Oceananigans.Grids.ZDirection}
└── coriolis: Nothing
```

With `HydrostaticFreeSurfaceModel`,

```jldoctest buoyancy
julia> model = HydrostaticFreeSurfaceModel(grid=grid, buoyancy=SeawaterBuoyancy(), tracers=(:T, :S))
HydrostaticFreeSurfaceModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: (:T, :S)
├── closure: Nothing
├── buoyancy: Buoyancy{SeawaterBuoyancy{Float64, LinearEquationOfState{Float64}, Nothing, Nothing}, Oceananigans.Grids.ZDirection}
├── free surface: ExplicitFreeSurface with gravitational acceleration 9.80665 m s⁻²
└── coriolis: Nothing
```

is identical to the default,

```jldoctest buoyancy
julia> model = HydrostaticFreeSurfaceModel(grid=grid)
HydrostaticFreeSurfaceModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: (:T, :S)
├── closure: Nothing
├── buoyancy: Buoyancy{SeawaterBuoyancy{Float64, LinearEquationOfState{Float64}, Nothing, Nothing}, Oceananigans.Grids.ZDirection}
├── free surface: ExplicitFreeSurface with gravitational acceleration 9.80665 m s⁻²
└── coriolis: Nothing
```

To model flows near the surface of Europa where `gravitational_acceleration = 1.3 \, \text{m}\,\text{s}^{-2}`,
we might alternatively specify

```jldoctest buoyancy
julia> buoyancy = SeawaterBuoyancy(gravitational_acceleration=1.3)
SeawaterBuoyancy{Float64}: g = 1.3
└── equation of state: LinearEquationOfState{Float64}: α = 1.67e-04, β = 7.80e-04

julia> model = NonhydrostaticModel(grid=grid, buoyancy=buoyancy, tracers=(:T, :S))
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: (:T, :S)
├── closure: Nothing
├── buoyancy: Buoyancy{SeawaterBuoyancy{Float64, LinearEquationOfState{Float64}, Nothing, Nothing}, Oceananigans.Grids.ZDirection}
└── coriolis: Nothing
```

for example.

### Linear equation of state

To specify the thermal expansion and haline contraction coefficients
``\alpha = 2 \times 10^{-3} \; \text{K}^{-1}`` and ``\beta = 5 \times 10^{-4} \text{psu}^{-1}``,

```jldoctest
julia> buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=2e-3, β=5e-4))
SeawaterBuoyancy{Float64}: g = 9.80665
└── equation of state: LinearEquationOfState{Float64}: α = 2.00e-03, β = 5.00e-04
```

### Idealized nonlinear equations of state

Instead of a linear equation of state, five idealized (second-order) nonlinear equation of state as described by
[Roquet15Idealized](@cite) may be used. These equations of state are provided via the
[SeawaterPolynomials.jl](https://github.com/CliMA/SeawaterPolynomials.jl) package.

```jldoctest buoyancy
julia> using SeawaterPolynomials.SecondOrderSeawaterPolynomials

julia> eos = RoquetSeawaterPolynomial(:Freezing)
SecondOrderSeawaterPolynomial{Float64}(0.7718, -0.0491, 0.0, -2.5681e-5, 0.0, -0.005027, 0.0)

julia> buoyancy = SeawaterBuoyancy(equation_of_state=eos)
SeawaterBuoyancy{Float64}: g = 9.80665
└── equation of state: SeawaterPolynomials.SecondOrderSeawaterPolynomials.SecondOrderSeawaterPolynomial{Float64}(0.7718, -0.0491, 0.0, -2.5681e-5, 0.0, -0.005027, 0.0)
```

### TEOS-10 equation of state

A high-accuracy 55-term polynomial approximation to the TEOS-10 equation of state suitable for use in
Boussinesq models as described by [Roquet15TEOS](@cite) is implemented in the
[SeawaterPolynomials.jl](https://github.com/CliMA/SeawaterPolynomials.jl) package and may be used.

```jldoctest buoyancy
julia> using SeawaterPolynomials.TEOS10

julia> eos = TEOS10EquationOfState()
SeawaterPolynomials.BoussinesqEquationOfState{TEOS10SeawaterPolynomial{Float64}, Int64}(TEOS10SeawaterPolynomial{Float64}(), 1020)
```

## The direction of gravitational acceleration

To simulate gravitational accelerations that don't align with the vertical (`z`) coordinate,
we wrap the buoyancy model in
`Buoyancy()` function call, which takes the keyword arguments `model` and `vertical_unit_vector`,

```jldoctest buoyancy
julia> θ = 45; # degrees

julia> g̃ = (0, sind(θ), cosd(θ));

julia> model = NonhydrostaticModel(grid=grid, 
                                   buoyancy=Buoyancy(model=BuoyancyTracer(), vertical_unit_vector=g̃), 
                                   tracers=:b)
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: (:b,)
├── closure: Nothing
├── buoyancy: Buoyancy{BuoyancyTracer, Tuple{Int64, Float64, Float64}}
└── coriolis: Nothing
```

# Lagrangian particle tracking

Models can keep track of the location and properties of neutrally buoyant particles. Particles are
advected with the flow field using forward Euler time-stepping at every model iteration.

## Simple particles

If you just need to keep of particle locations ``(x, y, z)`` then you can construct some Lagrangian particles
using the regular `LagrangianParticles` constructor

```@meta
DocTestSetup = quote
    using Oceananigans
end
```

```jldoctest particles
grid = RectilinearGrid(size=(10, 10, 10), extent=(1, 1, 1));

n_particles = 10;

x₀ = zeros(n_particles);

y₀ = rand(n_particles);

z₀ = -0.5 * ones(n_particles);

lagrangian_particles = LagrangianParticles(x=x₀, y=y₀, z=z₀)

# output
10 Lagrangian particles with
├── 3 properties: (:x, :y, :z)
└── 0 tracked fields: ()
```

then pass it to a model constructor

```jldoctest particles
model = NonhydrostaticModel(grid=grid, particles=lagrangian_particles)

# output
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 10×10×10 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: ()
├── closure: Nothing
├── buoyancy: Nothing
├── coriolis: Nothing
└── particles: 10 Lagrangian particles with 3 properties: (:x, :y, :z)
```

!!! warn "Lagrangian particles on GPUs"
    Remember to use `CuArray` instead of regular `Array` when storing particle locations and properties on the GPU.

## Custom particles

If you want to keep track of custom properties, such as the species or DNA of a Lagrangian particle
representing a microbe in an agent-based model, then you can create your own custom particle type
and pass a `StructArray` to the `LagrangianParticles` constructor.

```jldoctest particles
using Oceananigans
using StructArrays

struct LagrangianMicrobe{T, S, D}
    x :: T
    y :: T
    z :: T
    species :: S
    dna :: D
end

n_particles = 3;

x₀ = zeros(n_particles);

y₀ = rand(n_particles);

z₀ = -0.5 * ones(n_particles);

species = [:rock, :paper, :scissors]

dna = ["TATACCCC", "CCTAGGAC", "CGATTTAA"]

particles = StructArray{LagrangianMicrobe}((x₀, y₀, z₀, species, dna));

lagrangian_particles = LagrangianParticles(particles)

# output
3 Lagrangian particles with
├── 5 properties: (:x, :y, :z, :species, :dna)
└── 0 tracked fields: ()
```

!!! warn "Custom properties on GPUs"
    Not all data types can be passed to GPU kernels. If you intend to advect particles on the GPU make sure
    particle properties consist of only simple data types. The symbols and strings in this example won't
    work on the GPU.

## Writing particle properties to disk

Particle properties can be written to disk using JLD2 or NetCDF.

When writing to JLD2 you can pass `model.particles` as part of the named tuple of outputs.

```julia
JLD2OutputWriter(model, (particles=model.particles,), prefix="particles", schedule=TimeInterval(15))
```

When writing to NetCDF you should write particles to a separate file as the NetCDF dimensions differ for
particle trajectories. You can just pass `model.particles` straight to `NetCDFOutputWriter`:

```julia
NetCDFOutputWriter(model, model.particles, filepath="particles.nc", schedule=TimeInterval(15))
```

!!! warn "Outputting custom particle properties to NetCDF"
    NetCDF does not support arbitrary data types. If you need to write custom particle properties to disk
    that are not supported by NetCDF then you should use JLD2 (which should support almost any Julia data type).
# Number type

Passing `float_type=Float64` or `float_type=Float32` to the `Model` constructor causes the model to store all numbers
with 64-bit or 32-bit floating point precision.

!!! note "Avoiding mixed-precision operations"
    When not using `Float64` be careful to not mix different precisions as it could introduce implicit type conversions
    which can negatively effect performance. You can pass the number type desires to many constructors to enforce
    the type you want: e.g. `RectilinearGrid(CPU(), Float32; size=(16, 16, 16), extent=(1, 1, 1))` and
    `IsotropicDiffusivity(Float16; κ=1//7, ν=2//7)`.

!!! warning "Effect of floating point precision on simulation accuracy"
    While we run many tests with both `Float32` and `Float64` it is not clear whether `Float32` is precise enough to
    provide similar accuracy in all use cases. If accuracy is a concern, stick to `Float64`.

    We will be actively investigating the possibility of using lower precision floating point numbers such as `Float32`
    and `Float16` for fluid dynamics as well as the use of alternative number types such as Posits and Sonums.
# Turbulent diffusivity closures and large eddy simulation models

A turbulent diffusivity closure representing the effects of viscous dissipation and diffusion can be passed via the
`closure` keyword.

See [turbulence closures](@ref numerical_closures) and [large eddy simulation](@ref numerical_les) for more details
on turbulent diffusivity closures.

## Constant isotropic diffusivity

To use constant isotropic values for the viscosity ``\nu`` and diffusivity ``\kappa`` you can use [`IsotropicDiffusivity`](@ref)

```jldoctest
julia> using Oceananigans.TurbulenceClosures

julia> closure = IsotropicDiffusivity(ν=1e-2, κ=1e-2)
IsotropicDiffusivity: ν=0.01, κ=0.01
```

## Constant anisotropic diffusivity

To specify constant values for the horizontal and vertical viscosities, ``\nu_h`` and ``\nu_z``, and horizontal and vertical
diffusivities, ``\kappa_h`` and ``\kappa_z``, you can use [`AnisotropicDiffusivity`](@ref)

```jldoctest
julia> using Oceananigans.TurbulenceClosures

julia> closure = AnisotropicDiffusivity(νh=1e-3, νz=5e-2, κh=2e-3, κz=1e-1)
AnisotropicDiffusivity: (νx=0.001, νy=0.001, νz=0.05), (κx=0.002, κy=0.002, κz=0.1)
```

## Smagorinsky-Lilly

To use the Smagorinsky-Lilly LES closure, no parameters are required

```jldoctest
julia> using Oceananigans.TurbulenceClosures

julia> closure = SmagorinskyLilly()
SmagorinskyLilly: C=0.16, Cb=1.0, Pr=1.0, ν=0.0, κ=0.0
```

although they may be specified. By default, the background viscosity and diffusivity are assumed to be the molecular
values for seawater. For more details see [`SmagorinskyLilly`](@ref).

## Anisotropic minimum dissipation

To use the constant anisotropic minimum dissipation (AMD) LES closure,

```jldoctest
julia> using Oceananigans.TurbulenceClosures

julia> closure = AnisotropicMinimumDissipation()
AnisotropicMinimumDissipation{Float64} turbulence closure with:
           Poincaré constant for momentum eddy viscosity Cν: 0.08333333333333333
    Poincaré constant for tracer(s) eddy diffusivit(ies) Cκ: 0.08333333333333333
                        Buoyancy modification multiplier Cb: nothing
                Background diffusivit(ies) for tracer(s), κ: 0.0
             Background kinematic viscosity for momentum, ν: 0.0
```

no parameters are required although they may be specified. By default, the background viscosity and diffusivity
are assumed to be the molecular values for seawater. For more details see [`AnisotropicMinimumDissipation`](@ref).

## Convective Adjustment Vertical Diffusivity--Viscosity

To use the a convective adjustment scheme that applies enhanced values for vertical diffusivity ``\kappa_z`` and/or
viscosity ``\nu_z``, anytime and anywhere the background stratification becomes unstable.

```jldoctest
julia> using Oceananigans

julia> closure = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 1.0, background_κz = 1e-3)
ConvectiveAdjustmentVerticalDiffusivity: (background_κz=0.001, convective_κz=1.0, background_νz=0.0, convective_νz=0.0)
```
# Clock

The clock holds the current simulation time, iteration number, and time step stage.
The time step stage is relevant only for the multi-stage time-stepper `RungeKutta3TimeStepper`.

By default, `Clock`s are initialized at iteration 0, and stage 1,

```@meta
DocTestSetup = quote
    using Oceananigans
end
```

```jldoctest
julia> clock = Clock(time=0.0)
Clock{Float64}: time = 0 seconds, iteration = 0, stage = 1
```

but can be modified to start the model clock at some other time.
For example, passing

```jldoctest
julia> clock = Clock(time=3600.0)
Clock{Float64}: time = 1 hour, iteration = 0, stage = 1
```

to the constructor for `NonhydrostaticModel` causes the simulation
time to start at ``t = 3600`` seconds.

The type of the keyword argument `time` should be a float or date type.
To use the date type `TimeDate` from the `TimesDates.jl` package,
for example, pass

```jldoctest
julia> using TimesDates

julia> clock = Clock(time=TimeDate(2020))
Clock{TimesDates.TimeDate}: time = 2020-01-01T00:00:00, iteration = 0, stage = 1
```

to `NonhydrostaticModel`.
`TimeDate` supports nanosecond resolution and is thus recommended over `Base.Dates.DateTime`,
which is also supported but has only millisecond resolution.
# Grids

We currently support only `RectilinearGrid`s with either constant or variable grid spacings.
The spacings can be different for each dimension.

A `RectilinearGrid` is constructed by specifying the `size` of the grid (a `Tuple` specifying
the number of grid points in each direction) and either the `extent` (a `Tuple` specifying the
physical extent of the grid in each direction), or by prescribing `x`, `y`, and `z`. Keyword
arguments `x`, `y`, and `z` could be either *(i)* 2-`Tuple`s that define the the _end points_ in
each direction, or *(ii)* arrays or functions of the corresponding indices `i`, `j`, or `k` that
specify the locations of cell faces in the `x`-, `y`-, or `z`-direction, respectively.

A regular rectilinear grid with ``N_x \times N_y \times N_z = 32 \times 64 \times 256`` grid points
and an `extent` of ``L_x = 128`` meters, ``L_y = 256`` meters, and ``L_z = 512`` meters is constructed
by

```@meta
DocTestSetup = quote
    using Oceananigans
end
```

```jldoctest
julia> grid = RectilinearGrid(size=(32, 64, 256), extent=(128, 256, 512))
32×64×256 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── Periodic x ∈ [0.0, 128.0)  regularly spaced with Δx=4.0
├── Periodic y ∈ [0.0, 256.0)  regularly spaced with Δy=4.0
└── Bounded  z ∈ [-512.0, 0.0] regularly spaced with Δz=2.0
```

!!! info "Default domain"
    When using the `extent` keyword, e.g., `extent = (Lx, Ly, Lz)`, then the ``x \in [0, L_x]``,
    ``y \in [0, L_y]``, and ``z \in [-L_z, 0]`` -- a sensible choice for oceanographic applications.

## Specifying the grid's topology

Another crucial keyword is a 3-`Tuple` that specifies the grid's `topology`.
In each direction the grid may be `Periodic`, `Bounded` or `Flat`.
By default, both the `RectilinearGrid` and the `RectilinearGrid` constructors 
assume the grid topology is horizontally-periodic
and bounded in the vertical, such that `topology = (Periodic, Periodic, Bounded)`.

A "channel" model that is periodic in the ``x``-direction and wall-bounded
in the ``y``- and ``z``-dimensions is build with,

```jldoctest
julia> grid = RectilinearGrid(topology=(Periodic, Bounded, Bounded), size=(64, 64, 32), extent=(1e4, 1e4, 1e3))
64×64×32 RectilinearGrid{Float64, Periodic, Bounded, Bounded} on CPU with 1×1×1 halo
├── Periodic x ∈ [0.0, 10000.0) regularly spaced with Δx=156.25
├── Bounded  y ∈ [0.0, 10000.0] regularly spaced with Δy=156.25
└── Bounded  z ∈ [-1000.0, 0.0] regularly spaced with Δz=31.25
```

The `Flat` topology is useful when running problems with fewer than 3 dimensions. As an example,
to use a two-dimensional horizontal, doubly periodic domain the topology is `(Periodic, Periodic, Flat)`.


## Specifying domain end points

To specify a domain with a different origin than the default, the `x`, `y`, and `z` keyword arguments must be used.
For example, a grid with ``x \in [-100, 100]`` meters, ``y \in [0, 12.5]`` meters, and ``z \in [-\pi, \pi]`` meters
is constructed via

```jldoctest
julia> grid = RectilinearGrid(size=(32, 16, 256), x=(-100, 100), y=(0, 12.5), z=(-π, π))
32×16×256 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── Periodic x ∈ [-100.0, 100.0)     regularly spaced with Δx=6.25
├── Periodic y ∈ [0.0, 12.5)         regularly spaced with Δy=0.78125
└── Bounded  z ∈ [-3.14159, 3.14159] regularly spaced with Δz=0.0245437
```

## Grids with non-regular spacing in some of the directions

For a "channel" model, as the one we constructed above, one would probably like to have finer resolution near
the channel walls. We construct a grid that has non-regular spacing in the bounded dimensions, here ``y`` and ``z``
by by prescribing functions for `y` and `z` keyword arguments. 

For example, we can use the Chebychev nodes, which are more closely stacked near boundaries, to prescribe the
``y``- and ``z``-faces.

```jldoctest
julia> Nx, Ny, Nz = 64, 64, 32;

julia> Lx, Ly, Lz = 1e4, 1e4, 1e3;

julia> chebychev_spaced_y_faces(j) = - Ly/2 * cos(π * (j - 1) / Ny);

julia> chebychev_spaced_z_faces(k) = - Lz/2 - Lz/2 * cos(π * (k - 1) / Nz);

julia> grid = RectilinearGrid(size = (Nx, Ny, Nz),
                              topology=(Periodic, Bounded, Bounded),
                              x = (0, Lx),
                              y = chebychev_spaced_y_faces,
                              z = chebychev_spaced_z_faces)
64×64×32 RectilinearGrid{Float64, Periodic, Bounded, Bounded} on CPU with 1×1×1 halo
├── Periodic x ∈ [0.0, 10000.0)    regularly spaced with Δx=156.25
├── Bounded  y ∈ [-5000.0, 5000.0] variably spaced with min(Δy)=6.02272, max(Δy)=245.338
└── Bounded  z ∈ [-1000.0, 0.0]    variably spaced with min(Δz)=2.40764, max(Δz)=49.0086
```

```@setup 1
using Oceananigans
using Plots
Plots.scalefontsizes(1.25)
Plots.default(lw=3)
Nx, Ny, Nz = 64, 64, 32
Lx, Ly, Lz = 1e4, 1e4, 1e3
chebychev_spaced_y_faces(j) = - Ly/2 * cos(π * (j - 1) / Ny);
chebychev_spaced_z_faces(k) = - Lz/2 - Lz/2 * cos(π * (k - 1) / Nz);
grid = RectilinearGrid(size = (Nx, Ny, Nz),
                              topology=(Periodic, Bounded, Bounded),
                              x = (0, Lx),
                              y = chebychev_spaced_y_faces,
                              z = chebychev_spaced_z_faces)
```

We can easily visualize the spacing of ``y`` and ``z`` directions.

```@example 1
using Plots

py = plot(grid.yᵃᶜᵃ[1:Ny],  grid.Δyᵃᶜᵃ[1:Ny],
           marker = :circle,
           ylabel = "y-spacing (m)",
           xlabel = "y (m)",
           legend = nothing,
           ylims = (0, 250))

pz = plot(grid.Δzᵃᵃᶜ[1:Nz], grid.zᵃᵃᶜ[1:Nz],
          marker = :square,
          ylabel = "z (m)",
          xlabel = "z-spacing (m)",
          legend = nothing,
          xlims = (0, 50))

plot(py, pz, layout=(2, 1), size=(800, 900))

savefig("plot_stretched_grid.svg"); nothing # hide
```

![](plot_stretched_grid.svg)
# Tracers

The tracers to be advected around can be specified via a list of symbols. By default the model doesn't evolve any
tracers.

```@meta
DocTestSetup = quote
    using Oceananigans
end
```

```jldoctest tracers
julia> grid = RectilinearGrid(size=(64, 64, 64), extent=(1, 1, 1));

julia> model = NonhydrostaticModel(grid=grid)
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: ()
├── closure: Nothing
├── buoyancy: Nothing
└── coriolis: Nothing
```

But tracers can be added with the `tracers` keyword.
For example, to add conservative temperature `T` and absolute salinity `S`:

```jldoctest tracers
julia> grid = RectilinearGrid(size=(64, 64, 64), extent=(1, 1, 1));

julia> model = NonhydrostaticModel(grid=grid, tracers=(:T, :S))
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: (:T, :S)
├── closure: Nothing
├── buoyancy: Nothing
└── coriolis: Nothing
```

whose fields can be accessed via `model.tracers.T` and `model.tracers.S`.

```jldoctest tracers
julia> model.tracers.T
64×64×64 Field{Center, Center, Center} on RectilinearGrid on CPU
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── boundary conditions: west=Periodic, east=Periodic, south=Periodic, north=Periodic, bottom=ZeroFlux, top=ZeroFlux, immersed=ZeroFlux
└── data: 66×66×66 OffsetArray(::Array{Float64, 3}, 0:65, 0:65, 0:65) with eltype Float64 with indices 0:65×0:65×0:65
    └── max=0.0, min=0.0, mean=0.0

julia> model.tracers.S
64×64×64 Field{Center, Center, Center} on RectilinearGrid on CPU
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── boundary conditions: west=Periodic, east=Periodic, south=Periodic, north=Periodic, bottom=ZeroFlux, top=ZeroFlux, immersed=ZeroFlux
└── data: 66×66×66 OffsetArray(::Array{Float64, 3}, 0:65, 0:65, 0:65) with eltype Float64 with indices 0:65×0:65×0:65
    └── max=0.0, min=0.0, mean=0.0
```

An arbitrary number of tracers may be simulated. For example, to simulate
``C_1``, ``CO₂``, and `nitrogen` as additional passive tracers,

```jldoctest tracers
julia> model = NonhydrostaticModel(grid=grid, tracers=(:T, :S, :C₁, :CO₂, :nitrogen))
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 64×64×64 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: (:T, :S, :C₁, :CO₂, :nitrogen)
├── closure: Nothing
├── buoyancy: Nothing
└── coriolis: Nothing
```

!!! info "Active versus passive tracers"
    An active tracer is a tracer whose distribution affects the evolution of momentum and other tracers.
    Typical ocean models evolve conservative temperature and absolute salinity as active tracers,
    which effect momentum through buoyancy forces.
    Passive tracers are "passive" in the sense that their distribution does not affect
    the evolution of other tracers or flow quantities.
# Coriolis

The Coriolis option determines whether the fluid experiences the effect of the Coriolis force, or rotation. Currently
three options are available: no rotation, ``f``-plane, and ``\beta``-plane.

!!! info "Coriolis vs. rotation"
    If you are wondering why this option is called "Coriolis" it is because rotational effects could include the
    Coriolis and centripetal forces, both of which arise in non-inertial reference frames. But here the model only
    considers the Coriolis force.

## No rotation

By default there is no rotation. This can be made explicit by passing `coriolis = nothing` to a model constructor.

## Traditional ``f``-plane

To set up an ``f``-plane with, for example, Coriolis parameter ``f = 10^{-4} \text{s}^{-1}``

```@meta
DocTestSetup = quote
    using Oceananigans
end
```

```jldoctest
julia> coriolis = FPlane(f=1e-4)
FPlane{Float64}: f = 1.00e-04
```

An ``f``-plane can also be specified at some latitude on a spherical planet with a planetary rotation rate. For example,
to specify an ``f``-plane at a latitude of ``\varphi = 45°\text{N}`` on Earth which has a rotation rate of
``\Omega = 7.292115 \times 10^{-5} \text{s}^{-1}``

```jldoctest
julia> coriolis = FPlane(rotation_rate=7.292115e-5, latitude=45)
FPlane{Float64}: f = 1.03e-04
```

in which case the value of ``f`` is given by ``2\Omega\sin\varphi``.

## Coriolis term for constant rotation in a Cartesian coordinate system

One can use `ConstantCartesianCoriolis` to set up a Coriolis acceleration term where the Coriolis parameter
is constant and the rotation axis is arbitrary. For example, with
``\boldsymbol{f} = (0, f_y, f_z) = (0, 2, 1) \times 10^{-4} \text{s}^{-1}``,

```jldoctest
julia> coriolis = ConstantCartesianCoriolis(fx=0, fy=2e-4, fz=1e-4)
ConstantCartesianCoriolis{Float64}: fx = 0.00e+00, fy = 2.00e-04, fz = 1.00e-04
```

Or alternatively, the same result can be achieved by specifying the magnitude of the Coriolis
frequency `f` and the `rotation_axis`. So another way to get a Coriolis acceleration with the same
values is:

```jldoctest
julia> rotation_axis = (0, 2e-4, 1e-4)./√(2e-4^2 + 1e-4^2) # rotation_axis has to be a unit vector
(0.0, 0.8944271909999159, 0.4472135954999579)

julia> coriolis = ConstantCartesianCoriolis(f=√(2e-4^2+1e-4^2), rotation_axis=rotation_axis)
ConstantCartesianCoriolis{Float64}: fx = 0.00e+00, fy = 2.00e-04, fz = 1.00e-04
```

An ``f``-plane with non-traditional Coriolis terms can also be specified at some latitude on a spherical planet
with a planetary rotation rate. For example, to specify an ``f``-plane at a latitude of ``\varphi = 45°\text{N}``
on Earth which has a rotation rate of ``\Omega = 7.292115 \times 10^{-5} \text{s}^{-1}``

```jldoctest
julia> coriolis = ConstantCartesianCoriolis(rotation_rate=7.292115e-5, latitude=45)
ConstantCartesianCoriolis{Float64}: fx = 0.00e+00, fy = 1.03e-04, fz = 1.03e-04
```

in which case ``f_z = 2\Omega\sin\varphi`` and ``f_y = 2\Omega\cos\varphi``.

## Traditional ``\beta``-plane

To set up a ``\beta``-plane the background rotation rate ``f_0`` and the ``\beta`` parameter must be specified. For example,
a ``\beta``-plane with ``f_0 = 10^{-4} \text{s}^{-1}`` and ``\beta = 1.5 \times 10^{-11} \text{s}^{-1}\text{m}^{-1}`` can be
set up with

```jldoctest
julia> coriolis = BetaPlane(f₀=1e-4, β=1.5e-11)
BetaPlane{Float64}: f₀ = 1.00e-04, β = 1.50e-11
```

Alternatively, a ``\beta``-plane can also be set up at some latitude on a spherical planet with a planetary rotation rate
and planetary radius. For example, to specify a ``\beta``-plane at a latitude of ``\varphi = 10^\circ{S}`` on Earth
which has a rotation rate of ``\Omega = 7.292115 \times 10^{-5} \text{s}^{-1}`` and a radius of ``R = 6,371 \text{km}``

```jldoctest
julia> coriolis = BetaPlane(rotation_rate=7.292115e-5, latitude=-10, radius=6371e3)
BetaPlane{Float64}: f₀ = -2.53e-05, β = 2.25e-11
```

in which case ``f_0 = 2\Omega\sin\varphi`` and ``\beta = 2\Omega\cos\varphi / R``.

## Non-traditional ``\beta``-plane

A non-traditional ``\beta``-plane requires either 5 parameters (by default Earth's radius and
rotation rate are used):

```jldoctest
julia> NonTraditionalBetaPlane(fz=1e-4, fy=2e-4, β=4e-11, γ=-8e-11)
NonTraditionalBetaPlane{Float64}: fz = 1.00e-04, fy = 2.00e-04, β = 4.00e-11, γ = -8.00e-11, R = 6.37e+06
```

or the rotation rate, radius, and latitude:

```jldoctest
julia> NonTraditionalBetaPlane(rotation_rate=5.31e-5, radius=252.1e3, latitude=10)
NonTraditionalBetaPlane{Float64}: fz = 1.84e-05, fy = 1.05e-04, β = 4.15e-10, γ = -1.46e-10, R = 2.52e+05
```
# [Boundary conditions](@id model_step_bcs)

A boundary condition is applied to each field, dimension, and endpoint. There are left and right boundary conditions
for each of the x, y, and z dimensions so each field has 6 boundary conditions. Each of these boundary conditions may
be specified individually. Each boundary condition can be specified via a constant value, an array, or a function.

The left and right boundary conditions associated with the x-dimension are called west and east, respectively. For the
y-dimension, left and right are called south and north. For the z-dimension, left and right are called bottom and top.

See [Numerical implementation of boundary conditions](@ref numerical_bcs) for more details.

## Boundary condition classifications

1. [`Periodic`](@ref)
2. [`Flux`](@ref)
3. [`Value`](@ref) (Dirchlet)
4. [`Gradient`](@ref) (Neumann)
5. [`Open`](@ref)

Boundary conditions are constructed using the classification as a prefix: `FluxBoundaryCondition`, `ValueBoundaryCondition`, and so on.

## Starter tips

Here's a short list of useful tips for defining and prescribing boundary conditions on a model:

1. Boundary conditions depend on the grid topology and can only be non-default or non-`Periodic` in `Bounded` directions.
   Tracer boundary conditions are no flux by default in `Bounded` directions.
   Momentum boundary conditions are free-slip for tangential components and impenetrable for wall-normal components in `Bounded` directions.
   
2. Another way to say point 1 is that you'll never need to set:
    * `Periodic` boundary conditions (default for `Periodic` directions);
    * Impenetrable / "no normal flow" boundary conditions (default for wall-normal momentum components in `Bounded` directions);
    * "No flux" or "free slip" boundary conditions (default for tracers and wall-tangential momentum components in `Bounded` directions).

3. `ValueBoundaryCondition` (aka "Dirichlet" boundary conditions) models boundary fluxes given a field's diffusive flux model, and assuming that a field has the prescribed value on the boundary.
   _Note_: You cannot use `ValueBoundaryCondition` on a wall-normal velocity component; you must use `Open` for that.
   Examples where you might use `ValueBoundaryCondition`:
   * Prescribe a surface to have a constant temperature, like 20 degrees.
     Heat will then flux in and out of the domain depending on the temperature difference between the surface and the interior, and the temperature diffusivity.
   * Prescribe a velocity tangent to a boundary as in a driven-cavity flow (for example), where the top boundary is moving.
     Momentum will flux into the domain do the difference between the top boundary velocity and the interior velocity, and the prescribed viscosity.

4. `FluxBoundaryCondition` _directly_ prescribes the flux of a quantity across a boundary rather than calculating it given a viscosity or diffusivity.
   For example, sunlight absorbed at the ocean surface imparts a temperature flux that heats near-surface fluid.
   If there is a known `diffusivity`, you can express `FluxBoundaryCondition(flux)` using `GradientBoundaryCondition(-flux / diffusivity)` (aka "Neumann" boundary condition).
   But when `diffusivity` is not known or is variable (as for large eddy simulation, for example), it's convenient and more straightforward to apply `FluxBoundaryCondition`.

## Default boundary conditions

By default, periodic boundary conditions are applied on all fields along periodic dimensions. Otherwise tracers
get no-flux boundary conditions and velocities get free-slip and no normal flow boundary conditions.

## Boundary condition structures

Oceananigans uses a hierarchical structure to express boundary conditions:

1. Each boundary has one [`BoundaryCondition`](@ref)
2. Each field has seven [`BoundaryCondition`](@ref) (`west`, `east`, `south`, `north`, `bottom`, `top` and
   and an additional experimental condition for `immersed` boundaries)
3. A set of `FieldBoundaryConditions`, up to one for each field, are grouped into a `NamedTuple` and passed
   to the model constructor.

## Specifying boundary conditions for a model

Boundary conditions are defined at model construction time by passing a `NamedTuple` of `FieldBoundaryConditions`
specifying non-default boundary conditions for fields such as velocities and tracers.

Fields for which boundary conditions are not specified are assigned a default boundary conditions.

A few illustrations are provided below. See the examples for
further illustrations of boundary condition specification.

## Creating individual boundary conditions with `BoundaryCondition`

```@meta
DocTestSetup = quote
   using Oceananigans

   using Random
   Random.seed!(1234)
end
```

Boundary conditions may be specified with constants, functions, or arrays.
In this section we illustrate usage of the different [`BoundaryCondition`](@ref) constructors.

### 1. Constant `Value` (Dirchlet) boundary condition

```jldoctest
julia> constant_T_bc = ValueBoundaryCondition(20.0)
BoundaryCondition: classification=Value, condition=20.0
```

A constant [`Value`](@ref) boundary condition can be used to specify constant tracer (such as temperature),
or a constant _tangential_ velocity component at a boundary. Note that boundary conditions on the
_normal_ velocity component must use the [`Open`](@ref) boundary condition type.

Finally, note that `ValueBoundaryCondition(condition)` is an alias for `BoundaryCondition(Value, condition)`.

### 2. Constant `Flux` boundary condition

```jldoctest
julia> ρ₀ = 1027;  # Reference density [kg/m³]

julia> τₓ = 0.08;  # Wind stress [N/m²]

julia> wind_stress_bc = FluxBoundaryCondition(-τₓ/ρ₀)
BoundaryCondition: classification=Flux, condition=-7.789678675754625e-5
```

A constant [`Flux`](@ref) boundary condition can be imposed on tracers and tangential velocity components
that can be used, for example, to specify cooling, heating, evaporation, or wind stress at the ocean surface.

!!! info "The flux convention in Oceananigans"
    `Oceananigans` uses the convention that positive fluxes produce transport in the
    _positive_ direction (east, north, and up for ``x``, ``y``, ``z``).
    This means, for example, that a _negative_ flux of momentum or velocity at a _top_
    boundary, such as in the above example, produces currents in the _positive_ direction,
    because it prescribes a downwards flux of momentum into the domain from the top.
    Likewise, a _positive_ temperature flux at the top boundary
    causes _cooling_, because it transports heat _upwards_, out of the domain.
    Conversely, a positive flux at a _bottom_ boundary acts to increase the interior
    values of a quantity.

### 3. Spatially- and temporally-varying flux

Boundary conditions may be specified by functions,

```jldoctest
julia> @inline surface_flux(x, y, t) = cos(2π * x) * cos(t);

julia> top_tracer_bc = FluxBoundaryCondition(surface_flux)
BoundaryCondition: classification=Flux, condition=surface_flux(x, y, t) in Main at none:1
```

!!! info "Boundary condition functions"
    By default, a function boundary condition is called with the signature
    ```julia
    f(ξ, η, t)
    ```
    where `t` is time and `ξ, η` are spatial coordinates that vary along the boundary:
    * `f(y, z, t)` on `x`-boundaries;
    * `f(x, z, t)` on `y`-boundaries;
    * `f(x, y, t)` on `z`-boundaries.
    Alternative function signatures are specified by keyword arguments to
    `BoundaryCondition`, as illustrated in subsequent examples.

### 4. Spatially- and temporally-varying flux with parameters

Boundary condition functions may be 'parameterized',

```jldoctest
julia> @inline wind_stress(x, y, t, p) = - p.τ * cos(p.k * x) * cos(p.ω * t); # function with parameters

julia> top_u_bc = FluxBoundaryCondition(wind_stress, parameters=(k=4π, ω=3.0, τ=1e-4))
BoundaryCondition: classification=Flux, condition=wind_stress(x, y, t, p) in Main at none:1
```

!!! info "Boundary condition functions with parameters"
    The keyword argument `parameters` above specifies that `wind_stress` is called
    with the signature `wind_stress(x, y, t, parameters)`. In principle, `parameters` is arbitrary.
    However, relatively simple objects such as floating point numbers or `NamedTuple`s must be used
    when running on the GPU.

### 5. 'Field-dependent' boundary conditions

Boundary conditions may also depend on model fields. For example, a linear drag boundary condition
is implemented with

```jldoctest
julia> @inline linear_drag(x, y, t, u) = - 0.2 * u
linear_drag (generic function with 1 method)

julia> u_bottom_bc = FluxBoundaryCondition(linear_drag, field_dependencies=:u)
BoundaryCondition: classification=Flux, condition=linear_drag(x, y, t, u) in Main at none:1
```

`field_dependencies` specifies the name of the dependent fields either with a `Symbol` or `Tuple` of `Symbol`s.

### 6. 'Field-dependent' boundary conditions with parameters

When boundary conditions depends on fields _and_ parameters, their functions take the form

```jldoctest
julia> @inline quadratic_drag(x, y, t, u, v, drag_coeff) = - drag_coeff * u * sqrt(u^2 + v^2)
quadratic_drag (generic function with 1 method)

julia> u_bottom_bc = FluxBoundaryCondition(quadratic_drag, field_dependencies=(:u, :v), parameters=1e-3)
BoundaryCondition: classification=Flux, condition=quadratic_drag(x, y, t, u, v, drag_coeff) in Main at none:1
```

Put differently, `ξ, η, t` come first in the function signature, followed by field dependencies,
followed by `parameters` is `!isnothing(parameters)`.

### 7. Discrete-form boundary condition with parameters

Discrete field data may also be accessed directly from boundary condition functions
using the `discrete_form`. For example:

```jldoctest
@inline filtered_drag(i, j, grid, clock, model_fields) =
   @inbounds - 0.05 * (model_fields.u[i-1, j, 1] + 2 * model_fields.u[i, j, 1] + model_fields.u[i-1, j, 1])

u_bottom_bc = FluxBoundaryCondition(filtered_drag, discrete_form=true)

# output
BoundaryCondition: classification=Flux, condition=filtered_drag(i, j, grid, clock, model_fields) in Main at none:1
```

!!! info "The 'discrete form' for boundary condition functions"
    The argument `discrete_form=true` indicates to [`BoundaryCondition`](@ref) that `filtered_drag`
    uses the 'discrete form'. Boundary condition functions that use the 'discrete form'
    are called with the signature
    ```julia
    f(i, j, grid, clock, model_fields)
    ```
    where `i, j` are grid indices that vary along the boundary, `grid` is `model.grid`,
    `clock` is the `model.clock`, and `model_fields` is a `NamedTuple`
    containing `u, v, w` and the fields in `model.tracers`.
    The signature is similar for ``x`` and ``y`` boundary conditions expect that `i, j` is replaced
    with `j, k` and `i, k` respectively.

### 8. Discrete-form boundary condition with parameters

```jldoctest
julia> Cd = 0.2;  # drag coefficient

julia> @inline linear_drag(i, j, grid, clock, model_fields, Cd) = @inbounds - Cd * model_fields.u[i, j, 1];

julia> u_bottom_bc = FluxBoundaryCondition(linear_drag, discrete_form=true, parameters=Cd)
BoundaryCondition: classification=Flux, condition=linear_drag(i, j, grid, clock, model_fields, Cd) in Main at none:1
```

!!! info "Inlining and avoiding bounds-checking in boundary condition functions"
    Boundary condition functions should be decorated with `@inline` when running on CPUs for performance reasons.
    On the GPU, all functions are force-inlined by default.
    In addition, the annotation `@inbounds` should be used when accessing the elements of an array
    in a boundary condition function (such as `model_fields.u[i, j, 1]` in the above example).
    Using `@inbounds` will avoid a relatively expensive check that the index `i, j, 1` is 'in bounds'.

### 9. A random, spatially-varying, constant-in-time temperature flux specified by an array

```jldoctest
julia> Nx = Ny = 16;  # Number of grid points.

julia> Q = randn(Nx, Ny); # temperature flux

julia> white_noise_T_bc = FluxBoundaryCondition(Q)
BoundaryCondition: classification=Flux, condition=16×16 Matrix{Float64}
```

When running on the GPU, `Q` must be converted to a `CuArray`.

## Building boundary conditions on a field

To create a set of [`FieldBoundaryConditions`](@ref) for a temperature field,
we write

```jldoctest
julia> T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(20),
                                       bottom = GradientBoundaryCondition(0.01))
Oceananigans.FieldBoundaryConditions, with boundary conditions
├── west: Oceananigans.BoundaryConditions.DefaultPrognosticFieldBoundaryCondition
├── east: Oceananigans.BoundaryConditions.DefaultPrognosticFieldBoundaryCondition
├── south: Oceananigans.BoundaryConditions.DefaultPrognosticFieldBoundaryCondition
├── north: Oceananigans.BoundaryConditions.DefaultPrognosticFieldBoundaryCondition
├── bottom: BoundaryCondition{Gradient, Float64}
├── top: BoundaryCondition{Value, Int64}
└── immersed: BoundaryCondition{Flux, Nothing}
```

If the grid is, e.g., horizontally-periodic, then each horizontal `DefaultPrognosticFieldBoundaryCondition`
is converted to `PeriodicBoundaryCondition` inside the model's constructor, before assigning the
boundary conditions to temperature `T`.

In general, boundary condition defaults are inferred from the field location and `topology(grid)`.

## Specifying model boundary conditions

To specify non-default boundary conditions, a named tuple of [`FieldBoundaryConditions`](@ref) objects is
passed to the keyword argument `boundary_conditions` in the [`NonhydrostaticModel`](@ref) constructor.
The keys of `boundary_conditions` indicate the field to which the boundary condition is applied.
Below, non-default boundary conditions are imposed on the ``u``-velocity and temperature.

```jldoctest
julia> topology = (Periodic, Periodic, Bounded);

julia> grid = RectilinearGrid(size=(16, 16, 16), extent=(1, 1, 1), topology=topology);

julia> u_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(+0.1),
                                       bottom = ValueBoundaryCondition(-0.1));

julia> c_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(20),
                                       bottom = GradientBoundaryCondition(0.01));

julia> model = NonhydrostaticModel(grid=grid, boundary_conditions=(u=u_bcs, c=c_bcs), tracers=:c)
NonhydrostaticModel{CPU, Float64}(time = 0 seconds, iteration = 0)
├── grid: 16×16×16 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── tracers: (:c,)
├── closure: Nothing
├── buoyancy: Nothing
└── coriolis: Nothing

julia> model.velocities.u
16×16×16 Field{Face, Center, Center} on RectilinearGrid on CPU
├── grid: 16×16×16 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── boundary conditions: west=Periodic, east=Periodic, south=Periodic, north=Periodic, bottom=Value, top=Value, immersed=ZeroFlux
└── data: 18×18×18 OffsetArray(::Array{Float64, 3}, 0:17, 0:17, 0:17) with eltype Float64 with indices 0:17×0:17×0:17
    └── max=0.0, min=0.0, mean=0.0

julia> model.tracers.c
16×16×16 Field{Center, Center, Center} on RectilinearGrid on CPU
├── grid: 16×16×16 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── boundary conditions: west=Periodic, east=Periodic, south=Periodic, north=Periodic, bottom=Gradient, top=Value, immersed=ZeroFlux
└── data: 18×18×18 OffsetArray(::Array{Float64, 3}, 0:17, 0:17, 0:17) with eltype Float64 with indices 0:17×0:17×0:17
    └── max=0.0, min=0.0, mean=0.0
```

Notice that the specified non-default boundary conditions have been applied at
top and bottom of both `model.velocities.u` and `model.tracers.c`.

# Setting initial conditions

Initial conditions are imposed after model construction. This can be easily done using the the `set!` function, which
allows the setting of initial conditions using constant values, arrays, or functions.

```@meta
DocTestSetup = quote
    using Oceananigans
end
```

```jldoctest
julia> grid = RectilinearGrid(size=(16, 16, 16), extent=(1, 1, 1));

julia> model = NonhydrostaticModel(grid=grid);

julia> set!(model, u=0.1, v=1.5)
```

```jldoctest
julia> grid = RectilinearGrid(size=(16, 16, 16), extent=(1, 1, 1));

julia> model = NonhydrostaticModel(grid=grid, buoyancy=SeawaterBuoyancy(), tracers=(:T, :S));

julia> ∂T∂z = 0.01;

julia> ϵ(σ) = σ * randn();

julia> T₀(x, y, z) = ∂T∂z * z + ϵ(1e-8);

julia> set!(model, T=T₀)
```

!!! tip "Divergence-free velocity fields"
    Note that as part of the time-stepping algorithm, the velocity field is made
    divergence-free at every time step. So if a model is not initialized with a
    divergence-free velocity field, it may change on the first time step. As a result
    tracers may not be conserved up to machine precision at the first time step.
# Architecture

Passing `architecture = CPU()` or `architecture = GPU()` to the `Model` constructor will determine whether the model
is time stepped on a CPU or GPU.

Ideally a set up or simulation script does not need to be modified to run on a GPU but we are still smoothing out
rough edges. Generally the CPU wants `Array` objects while the GPU wants `CuArray` objects.

!!! tip "Running on GPUs"
    If you are having issues with running Oceananigans on a GPU, please
    [open an issue](https://github.com/CLiMA/Oceananigans.jl/issues/new) and we'll do our best to help out.
# Background fields

`BackgroundField`s are velocity and tracer fields around which the resolved
velocity and tracer fields evolve. In `Oceananigans`, only the _advective_ terms
associated with the interaction between background and resolved fields are included.
For example, tracer advection is described by

```math
\boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{v} c \right ) \, ,
```

where ``\boldsymbol{v}`` is the resolved velocity field and ``c`` is the resolved
tracer field corresponding to `model.tracers.c`. 

When a background field ``C`` is provided, the tracer advection term becomes

```math
\boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{v} c \right ) 
    + \boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{v} C \right ) \, .
```

When both a background field velocity field ``\boldsymbol{U}`` and a background tracer field ``C``
are provided, then the tracer advection term becomes

```math
\boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{v} c \right ) 
    + \boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{v} C \right )
    + \boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{U} c \right ) \, .
```

Notice that the term ``\boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \boldsymbol{U} C \right )`` 
is neglected: only the terms describing the advection of resolved tracer by the background 
velocity field and the advection of background tracer by the resolved velocity field are included.
An analgous statement holds for the advection of background momentum by the resolved
velocity field.
Other possible terms associated with the Coriolis force, buoyancy, turbulence closures,
and surface waves acting on background fields are neglected.

## Specifying background fields

`BackgroundField`s are defined by functions of ``(x, y, z, t)`` and optional parameters. A 
simple example is

```jldoctest
using Oceananigans

U(x, y, z, t) = 0.2 * z

grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))

model = NonhydrostaticModel(grid = grid, background_fields = (u=U,))

model.background_fields.velocities.u

# output
FunctionField located at (Face, Center, Center)
├── func: U (generic function with 1 method)
├── grid: 1×1×1 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── clock: Clock(time=0 seconds, iteration=0)
└── parameters: nothing
```

`BackgroundField`s are specified by passing them to the kwarg `background_fields`
in the `NonhydrostaticModel` constructor. The kwarg `background_fields` expects
a `NamedTuple` of fields, which are internally sorted into `velocities` and `tracers`,
wrapped in `FunctionField`s, and assigned their appropriate locations.

`BackgroundField`s with parameters require using the `BackgroundField` wrapper:

```jldoctest moar_background
using Oceananigans

parameters = (α=3.14, N=1.0, f=0.1)

# Background fields are defined via function of x, y, z, t, and optional parameters
U(x, y, z, t, α) = α * z
B(x, y, z, t, p) = - p.α * p.f * y + p.N^2 * z 

U_field = BackgroundField(U, parameters=parameters.α)
B_field = BackgroundField(B, parameters=parameters)

# output
BackgroundField{typeof(B), NamedTuple{(:α, :N, :f), Tuple{Float64, Float64, Float64}}}
├── func: B (generic function with 1 method)
└── parameters: (α = 3.14, N = 1.0, f = 0.1)
```

When inserted into `NonhydrostaticModel`, we get

```jldoctest moar_background
grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))

model = NonhydrostaticModel(grid = grid, background_fields = (u=U_field, b=B_field),
                            tracers=:b, buoyancy=BuoyancyTracer())

model.background_fields.tracers.b

# output
FunctionField located at (Center, Center, Center)
├── func: B (generic function with 1 method)
├── grid: 1×1×1 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 1×1×1 halo
├── clock: Clock(time=0 seconds, iteration=0)
└── parameters: (α = 3.14, N = 1.0, f = 0.1)
```

# Forcing functions

"Forcings" are user-defined terms appended to right-hand side of
the momentum or tracer evolution equations. In `Oceananigans`, momentum
and tracer forcings are defined via julia functions. `Oceananigans` includes
an interface for implementing forcing functions that depend on spatial coordinates,
time, model velocity and tracer fields, and external parameters.

```@meta
DocTestSetup = quote
    using Oceananigans
end
```

Forcings are added to `Oceananigans` models by passing a `NamedTuple` of functions
or forcing objects to the `forcing` keyword argument in `NonhydrostaticModel`'s constructor.
By default, momentum and tracer forcing functions are assumed to be functions of
`x, y, z, t`. A basic example is

```jldoctest
u_forcing(x, y, z, t) = exp(z) * cos(x) * sin(t)

grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))
model = NonhydrostaticModel(grid=grid, forcing=(u=u_forcing,))

model.forcing.u

# output
ContinuousForcing{Nothing} at (Face, Center, Center)
├── func: u_forcing (generic function with 1 method)
├── parameters: nothing
└── field dependencies: ()
```

More general forcing functions are built via the `Forcing` constructor
described below. `Oceananigans` also provides a convenience type called `Relaxation`
that specifies "relaxation", or damping terms that restore a field to a
target distribution outside of a masked region of space. `Relaxation` can be
used to implement sponge layers near the boundaries of a domain.

## The `Forcing` constructor

The `Forcing` constructor provides an interface for specifying forcing functions that

1. Depend on external parameters; and
2. Depend on model fields at the `x, y, z` location that forcing is applied; and/or
3. Require access to discrete model data.

### Forcing functions with external parameters

Most forcings involve external, changeable parameters.
Here are two examples of `forcing_func`tions that depend on 
_(i)_ a single scalar parameter `s`, and _(ii)_ a `NamedTuple` of parameters, `p`:

```jldoctest parameterized_forcing
# Forcing that depends on a scalar parameter `s`
u_forcing_func(x, y, z, t, s) = s * z

u_forcing = Forcing(u_forcing_func, parameters=0.1)

# Forcing that depends on a `NamedTuple` of parameters `p`
T_forcing_func(x, y, z, t, p) = - p.μ * exp(z / p.λ) * cos(p.k * x) * sin(p.ω * t)

T_forcing = Forcing(T_forcing_func, parameters=(μ=1, λ=0.5, k=2π, ω=4π))

grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))
model = NonhydrostaticModel(grid=grid, forcing=(u=u_forcing, T=T_forcing), buoyancy=SeawaterBuoyancy(), tracers=(:T, :S))

model.forcing.T

# output
ContinuousForcing{NamedTuple{(:μ, :λ, :k, :ω), Tuple{Int64, Float64, Float64, Float64}}} at (Center, Center, Center)
├── func: T_forcing_func (generic function with 1 method)
├── parameters: (μ = 1, λ = 0.5, k = 6.283185307179586, ω = 12.566370614359172)
└── field dependencies: ()
```

```jldoctest parameterized_forcing
model.forcing.u

# output
ContinuousForcing{Float64} at (Face, Center, Center)
├── func: u_forcing_func (generic function with 1 method)
├── parameters: 0.1
└── field dependencies: ()
```

In this example, the objects passed to the `parameters` keyword in the construction of
`u_forcing` and `T_forcing` --- a floating point number for `u_forcing`, and a `NamedTuple`
of parameters for `T_forcing` --- are passed on to `u_forcing_func` and `T_forcing_func` when
they are called during time-stepping. The object passed to `parameters` is in principle arbitrary.
However, if using the GPU, then `typeof(parameters)` may be restricted by the requirements
of GPU-compiliability.

### Forcing functions that depend on model fields

Forcing functions may depend on model fields evaluated at the `x, y, z` where forcing is applied.
Here's a somewhat non-sensical example:

```jldoctest field_dependent_forcing
# Forcing that depends on the velocity fields `u`, `v`, and `w`
w_forcing_func(x, y, z, t, u, v, w) = - (u^2 + v^2 + w^2) / 2

w_forcing = Forcing(w_forcing_func, field_dependencies=(:u, :v, :w))

# Forcing that depends on salinity `S` and a scalar parameter
S_forcing_func(x, y, z, t, S, μ) = - μ * S

S_forcing = Forcing(S_forcing_func, parameters=0.01, field_dependencies=:S)

grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))
model = NonhydrostaticModel(grid=grid, forcing=(w=w_forcing, S=S_forcing), buoyancy=SeawaterBuoyancy(), tracers=(:T, :S))

model.forcing.w

# output
ContinuousForcing{Nothing} at (Center, Center, Face)
├── func: w_forcing_func (generic function with 1 method)
├── parameters: nothing
└── field dependencies: (:u, :v, :w)
```

```jldoctest field_dependent_forcing
model.forcing.S

# output
ContinuousForcing{Float64} at (Center, Center, Center)
├── func: S_forcing_func (generic function with 1 method)
├── parameters: 0.01
└── field dependencies: (:S,)
```

The `field_dependencies` arguments follow `x, y, z, t` in the forcing `func`tion in
the order they are specified in `Forcing`.
If both `field_dependencies` and `parameters` are specified, then the `field_dependencies`
arguments follow `x, y, z, t`, and `parameters` follow `field_dependencies`.

Model fields that arise in the arguments of continuous `Forcing` `func`tions are
automatically interpolated to the staggered grid location at which the forcing is applied.

### "Discrete form" forcing functions

"Discrete form" forcing functions are either called with the signature

```julia
func(i, j, k, grid, clock, model_fields)
```

or the parameterized form

```julia
func(i, j, k, grid, clock, model_fields, parameters)
```

Discrete form forcing functions can access the entirety of model field
data through the argument `model_fields`. The object `model_fields` is a `NamedTuple`
whose properties include the velocity fields `model_fields.u`, `model_fields.v`,
`model_fields.w` and all fields in `model.tracers`.

Using discrete forcing functions may require understanding the
staggered arrangement of velocity fields and tracers in `Oceananigans`.
Here's a slightly non-sensical example in which the vertical derivative of a buoyancy
tracer is used as a time-scale for damping the u-velocity field:

```jldoctest discrete_forcing
# A damping term that depends on a "local average":
local_average(i, j, k, grid, c) = @inbounds (c[i, j, k] + c[i-1, j, k] + c[i+1, j, k] +
                                                          c[i, j-1, k] + c[i, j+1, k] +
                                                          c[i, j, k-1] + c[i, j, k+1]) / 7

b_forcing_func(i, j, k, grid, clock, model_fields) = - local_average(i, j, k, grid, model_fields.b)

b_forcing = Forcing(b_forcing_func, discrete_form=true)

# A term that damps the local velocity field in the presence of stratification
using Oceananigans.Operators: ∂zᵃᵃᶠ, ℑxzᶠᵃᶜ

function u_forcing_func(i, j, k, grid, clock, model_fields, ε)
    # The vertical derivative of buoyancy, interpolated to the u-velocity location:
    N² = ℑxzᶠᵃᶜ(i, j, k, grid, ∂zᵃᵃᶠ, model_fields.b)

    # Set to zero in unstable stratification where N² < 0:
    N² = max(N², zero(typeof(N²)))

    return @inbounds - ε * sqrt(N²) * model_fields.u[i, j, k]
end

u_forcing = Forcing(u_forcing_func, discrete_form=true, parameters=1e-3)

grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))
model = NonhydrostaticModel(grid=grid, tracers=:b, buoyancy=BuoyancyTracer(), forcing=(u=u_forcing, b=b_forcing))

model.forcing.b

# output
DiscreteForcing{Nothing}
├── func: b_forcing_func (generic function with 1 method)
└── parameters: nothing
```

```jldoctest discrete_forcing
model.forcing.u

# output
DiscreteForcing{Float64}
├── func: u_forcing_func (generic function with 1 method)
└── parameters: 0.001
```

The annotation `@inbounds` is crucial for performance when accessing array indices
of the fields in `model_fields`.

## `Relaxation`

`Relaxation` defines a special forcing function that restores a field at a specified `rate` to
a `target` distribution, within a region uncovered by a `mask`ing function.
`Relaxation` is useful for implementing sponge layers, as shown in the second example.

The following code constructs a model in which all components
of the velocity field are damped to zero everywhere on a time-scale of 1000 seconds, or ~17 minutes:

```jldoctest
damping = Relaxation(rate = 1/1000)

grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1)) 
model = NonhydrostaticModel(grid=grid, forcing=(u=damping, v=damping, w=damping))

model.forcing.w

# output
ContinuousForcing{Nothing} at (Center, Center, Face)
├── func: Relaxation(rate=0.001, mask=1, target=0)
├── parameters: nothing
└── field dependencies: (:w,)
```

The constructor for `Relaxation` accepts the keyword arguments `mask`, and `target`,
which specify a `mask(x, y, z)` function that multiplies the forcing, and a `target(x, y, z)`
distribution for the quantity in question. By default, `mask` uncovered the whole domain
and `target` restores the field in question to 0

We illustrate usage of `mask` and `target` by implementing a sponge layer that relaxes
velocity fields to zero and restores temperature to a linear gradient in the bottom
1/10th of the domain:

```jldoctest sponge_layer
grid = RectilinearGrid(size=(1, 1, 1), x=(0, 1), y=(0, 1), z=(-1, 0))

        damping_rate = 1/100 # relax fields on a 100 second time-scale
temperature_gradient = 0.001 # ⁰C m⁻¹
 surface_temperature = 20    # ⁰C (at z=0)

target_temperature = LinearTarget{:z}(intercept=surface_temperature, gradient=temperature_gradient)
       bottom_mask = GaussianMask{:z}(center=-grid.Lz, width=grid.Lz/10)

uvw_sponge = Relaxation(rate=damping_rate, mask=bottom_mask)
  T_sponge = Relaxation(rate=damping_rate, mask=bottom_mask, target=target_temperature)

model = NonhydrostaticModel(grid=grid, forcing=(u=uvw_sponge, v=uvw_sponge, w=uvw_sponge, T=T_sponge), buoyancy=SeawaterBuoyancy(), tracers=(:T, :S))

model.forcing.u

# output
ContinuousForcing{Nothing} at (Face, Center, Center)
├── func: Relaxation(rate=0.01, mask=exp(-(z + 1.0)^2 / (2 * 0.1^2)), target=0)
├── parameters: nothing
└── field dependencies: (:u,)
```

```jldoctest sponge_layer
model.forcing.T

# output
ContinuousForcing{Nothing} at (Center, Center, Center)
├── func: Relaxation(rate=0.01, mask=exp(-(z + 1.0)^2 / (2 * 0.1^2)), target=20.0 + 0.001 * z)
├── parameters: nothing
└── field dependencies: (:T,)
```

# Diagnostics

Diagnostics are a set of general utilities that can be called on-demand during time-stepping to compute quantities of
interest you may want to save to disk, such as the horizontal average of the temperature, the maximum velocity, or to
produce a time series of salinity. They also include utilities for diagnosing model health, such as the CFL number or
to check for NaNs.

Diagnostics are stored as a list of diagnostics in `simulation.diagnostics`. Diagnostics can be specified at model creation
time or be specified at any later time and appended (or assigned with a key value pair) to `simulation.diagnostics`.

Most diagnostics can be run at specified frequencies (e.g. every 25 time steps) or specified intervals (e.g. every
15 minutes of simulation time). If you'd like to run a diagnostic on demand then do not specify any intervals
(and do not add it to `simulation.diagnostics`).
# Checkpointing

A checkpointer can be used to serialize the entire model state to a file from which the model
can be restored at any time. This is useful if you'd like to periodically checkpoint when running
long simulations in case of crashes or cluster time limits, but also if you'd like to restore
from a checkpoint and try out multiple scenarios.

For example, to periodically checkpoint the model state to disk every 1,000,000 seconds of simulation
time to files of the form `model_checkpoint_iteration12500.jld2` where `12500` is the iteration
number (automatically filled in)

```@meta
DocTestSetup = quote
    using Oceananigans
end
```

```@repl checkpointing
using Oceananigans, Oceananigans.Units

model = NonhydrostaticModel(grid=RectilinearGrid(size=(16, 16, 16), extent=(1, 1, 1)))

simulation = Simulation(model, Δt=1, stop_iteration=1)

simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=TimeInterval(5years), prefix="model_checkpoint")

run!(simulation)
```

The default options should provide checkpoint files that are easy to restore from in most cases.
For more advanced options and features, see [`Checkpointer`](@ref).

## Picking up a simulation from a checkpoint file

Picking up a simulation from a checkpoint requires the original script that was used to generate
the checkpoint data. Change the first instance of [`run!`](@ref) in the script to take `pickup=true`:

```@repl checkpointing
simulation.stop_iteration = 2

run!(simulation, pickup=true)
```

which finds the latest checkpoint file in the current working directory (in this trivial case,
this is the checkpoint associated with iteration 0), loads prognostic fields and their tendencies
from file, resets the model clock and iteration, and updates the model auxiliary state before
starting the time-stepping loop.

Use `pickup=iteration`, where `iteration` is an `Integer`, to pick up from a specific iteration.
Or, use `pickup=filepath`, where `filepath` is a string, to pickup from a specific file located
at `filepath`.
# Model setup

This section describes all the options and features that can be used to set up a model. For 
more detailed information consult the API documentation.

Each structure covered in this section can be constructed and passed to the models' constructors. 
For examples of model construction, see the examples. The validation experiments provide more 
advanced examples.

For reference, here are all the option or keyword arguments that can be passed to the models
currently implemented in Oceananigans.jl. See the different sections on the sidebar for more 
details and examples for each keyword argument.

### `NonhydrostaticModel`

```@docs
NonhydrostaticModel
```

### `HydrostaticFreeSurfaceModel`

```@docs
HydrostaticFreeSurfaceModel
```

### `ShallowWaterModel`

```@docs
ShallowWaterModel
```
# Index

```@index
```
# Convergence Tests

Convergence tests are implemented in `/validation/convergence_tests` and range
from zero-dimensional time-stepper tests to two-dimensional integration tests that
involve non-trivial pressure fields, advection, and diffusion.

For all tests except point exponential decay, we use the ``L_1`` norm,

```math
    L_1 \equiv \frac{\mathrm{mean} | \phi_\mathrm{sim} - \phi_\mathrm{exact} |}{\mathrm{mean} | \phi_\mathrm{exact} |}
```

and ``L_\infty`` norm,

```math
    L_\infty \equiv \frac{\max | \phi_\mathrm{sim} - \phi_\mathrm{exact} |}{\max | \phi_\mathrm{exact} |} \, ,
```

to compare simulated fields, ``\phi_\mathrm{sim}``, with exact, analytically-derived solutions
``\phi_\mathrm{exact}``.
The field ``\phi`` may be a tracer field or a velocity field.

## Point Exponential Decay

This test analyzes time-stepper convergence by simulating the zero-dimensional, or spatially-uniform equation

```math
    \partial_t c = - c \, ,
```

with the initial condition ``c = 1``, which has the analytical solution ``c = \mathrm{e}^{-t}``.

We find the expected first-order convergence with decreasing time-step ``\Delta t`` using our
first-order accurate, "modified second-order" Adams-Bashforth time-stepping method:

![Point exponential decay](https://raw.githubusercontent.com/CliMA/Oceananigans.jl/v0.35.0/docs/src/verification/convergence_plots/point_exponential_decay_time_stepper_convergence.png)

This result validates the correctness of the `Oceananigans` implementation of Adams-Bashforth time-stepping.

## One-dimensional advection and diffusion of a Gaussian

This and the following tests focus on convergence with grid spacing, ``\Delta x``.

In one dimension with constant diffusivity ``\kappa`` and in the presence of a
constant velocity ``U``, a Gaussian evolves according to

```math
    c = \frac{\mathrm{e}^{- (x - U t)^2 / 4 \kappa t}}{\sqrt{4 \pi \kappa t}} \, .
```

For this test we take the initial time as ``t=t_0``.
We simulate this problem with advection and diffusion, as well as with ``U=0`` and thus diffusion only, as well as with
``\kappa \approx 0`` and thus "advection only".
The solutions are

![Gaussian advection diffusion solutions](https://raw.githubusercontent.com/CliMA/Oceananigans.jl/v0.35.0/docs/src/verification/convergence_plots/gaussian_advection_diffusion_solutions.png)

which exhibit the expected second-order convergence with ``\Delta x^2 \propto 1 / N_x^2``:

![Gaussian advection diffusion convergence](https://raw.githubusercontent.com/CliMA/Oceananigans.jl/v0.35.0/docs/src/verification/convergence_plots/gaussian_advection_diffusion_error_convergence.png)

These results validate the correctness of time-stepping, constant diffusivity operators, and advection operators.

## One-dimensional advection and diffusion of a cosine

In one dimension with constant diffusivity ``\kappa`` and in the presence of a
constant velocity ``U``, a cosine evolves according to

```math
    c = \mathrm{e}^{-\kappa t} \cos (x - U t) \, .
```

The solutions are

![Cosine advection diffusion solutions](https://raw.githubusercontent.com/CliMA/Oceananigans.jl/v0.35.0/docs/src/verification/convergence_plots/cosine_advection_diffusion_solutions.png)

which exhibit the expected second-order convergence with ``\Delta x^2 \propto 1 / N_x^2``:

![Cosine advection diffusion convergence](https://raw.githubusercontent.com/CliMA/Oceananigans.jl/v0.35.0/docs/src/verification/convergence_plots/cosine_advection_diffusion_error_convergence.png)

These results validate the correctness of time-stepping, constant diffusivity operators, and advection operators.

## Two-dimensional diffusion

With zero velocity field and constant diffusivity ``\kappa``, the tracer field

```math
    c(x, y, t=0) = \cos(x) \cos(y) \, ,
```

decays according to

```math
    c(x, y, t) = \mathrm{e}^{-2 \kappa t} \cos(x) \cos(y) \, ,
```

with either periodic boundary conditions, or insulating boundary conditions in either ``x`` or ``y``.

The expected convergence with ``\Delta x^2 \propto 1 / N_x^2`` is observed:

![Two dimensional diffusion convergence](https://raw.githubusercontent.com/CliMA/Oceananigans.jl/v0.35.0/docs/src/verification/convergence_plots/two_dimensional_diffusion_convergence.png)

This validates the correctness of multi-dimensional diffusion operators.

## Decaying, advected Taylor-Green vortex

The velocity field

```math
\begin{align}
    u(x, y, t) & = U + \mathrm{e}^{-t} \cos(x - U t) \sin(y) \, , \\
    v(x, y, t) & =   - \mathrm{e}^{-t} \sin(x - U t) \cos(y) \, ,
\end{align}
```

is a solution to the Navier-Stokes equations with viscosity ``\nu = 1``.

The expected convergence with spatial resolution is observed:

![Decaying advected Taylor Green](https://raw.githubusercontent.com/CliMA/Oceananigans.jl/v0.35.0/docs/src/verification/convergence_plots/taylor_green_convergence.png)

This validates the correctness of the advection and diffusion of a velocity field.

## Forced two-dimensional flows

We introduce two convergence tests associated with forced flows in domains that are 
bounded in ``y``, and periodic in ``x`` with no tracers.

*Note: in this section, subscripts are used to denote derivatives to make reading 
and typing equations easier.*

In a two-dimensional flow in ``(x, y)``, the velocity field ``(u, v)`` can be expressed in terms
of a streamfunction ``\psi(x, y, t)`` such that

```math
    u \equiv - \psi_y \, , \quad \text{and} \quad v \equiv \psi_x \, ,
```
where subscript denote derivatives such that ``\psi_y \equiv \partial_y \psi``, for example.
With an isotropic Laplacian viscosity ``\nu = 1``, the momentum and continuity equations are
```math
    \begin{align}
    \boldsymbol{v}_t + \left ( \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} \right ) \boldsymbol{v} + \boldsymbol{\nabla} p & = \nabla^2 \boldsymbol{v} + \boldsymbol{F}_v \, , \\
    \boldsymbol{\nabla} \boldsymbol{\cdot} \boldsymbol{v} & = 0 \, ,
    \end{align}
```
while the equation for vorticity, ``\omega = v_x - u_y = \nabla^2 \psi``, is
```math
    \omega_t + \mathrm{J} \left ( \psi, \omega \right ) = \nabla^2 \omega + F_\omega \, .
```
Finally, taking the divergence of the momentum equation, we find a Poisson equation for pressure,
```math
    \nabla^2 p = - u_x^2 - v_y^2 - 2 u_y v_x + \partial_x F_v + \partial_y F_v \, .
```

To pose the problem, we first pick a streamfunction ``\psi``. This choice then yields the vorticity 
forcing ``F_{\omega}`` that satisfies the vorticity equation. We then determine ``F_u`` by solving 
``\partial_y F_v = - F_{\omega}``, and pick ``F_v`` so that we can solve the Poisson equation 
for pressure.

We restrict ourselves to a class of problems in which
```math
\psi(x, y, t) = - f(x, t) g(y) \, , \quad \text{with} \quad f \equiv \cos [x - \xi(t)] \, , \quad
\xi(t) \equiv 1 + \sin(t^2) \, .
```
Grinding through the algebra, this particular form implies that ``F_{\omega}`` is given by
```math
    F_{\omega} = -\xi^\prime f_x (g - g^{\prime\prime}) + f f_x (g g^{\prime\prime\prime} - g^\prime g^{\prime\prime}) + f (g - 2 g^{\prime\prime} + g^{\prime\prime\prime\prime}) \, ,
```
where primes denote derivatives of functions of a single argument. 
Setting ``\partial_y F_v = F_{\omega}``, we find that if ``F_v`` satisfies
```math
    \partial_y F_v = (g^\prime)^2 + g g^{\prime\prime} \, ,
```
then the pressure Poisson equation becomes
```math
    \nabla^2 p = \cos [2 (x - \xi)] [(g^\prime)^2 - g g^{\prime\prime}] + \partial_x F_v \, .
```
This completes the specification of the problem.

We set up the problem by imposing the time-dependent forcing functions ``F_u`` and ``F_v``
on ``u`` and ``v``, initializing the flow at ``t=0``, and integrating the problem forwards
in time using Oceananigans. We find the expected convergence of the numerical solution to the
analytical solution: the error between the numerical and analytical solutions
decreases with ``1/N_x^2 \sim \Delta x^2``, where ``N_x`` is the number of grid
points and ``\Delta x`` is the spatial resolution:

![Forced free slip convergence](https://raw.githubusercontent.com/CliMA/Oceananigans.jl/v0.35.0/docs/src/verification/convergence_plots/forced_free_slip_convergence.png)

The convergence tests are performed using both ``y`` and ``z`` as the bounded direction.

### Forced, free-slip flow

A forced flow satisfying free-slip conditions at ``y = 0`` and ``y = \pi`` has the streamfunction
```math
    \psi(x, y, t) = - \cos [x - \xi(t)] \sin (y) \, ,
```
and thus ``g(y) = \sin y``. The velocity field ``(u, v)`` is
```math
    u = \cos (x - \xi) \cos y \, , \quad \text{and} \quad v = \sin (x - \xi) \sin y \, ,
```
which satisfies the boundary conditions ``u_y |_{y=0} = u_y |_{y=\pi} = 0`` and
``v |_{y=0} = v |_{y=\pi} = 0``. The vorticity forcing is
```math
    F_{\omega} = - 2 \xi^\prime f_x \sin y + 4 f \sin y \, ,
```
which implies that
```math
    F_v = - 2 \xi^\prime f_x \cos y + 4 f \cos y \, ,
```
and ``F_v = \tfrac{1}{2} \sin 2 y``.

## Forced, fixed-slip flow

A forced flow satisfying "fixed-slip" boundary conditions at ``y=0`` and ``y=1`` has
the streamfunction
```math
    \psi(x, y, t) = - \cos [x - \xi(t)] (y^3 - y^2) \, ,
```
and thus ``g(y) = y^3 - y^2``. The velocity field ``(u, v)`` is
```math
    u = f (3y^2 - 2 y) \, , \quad \text{and} \quad v = - f_x (y^3 - y^2) \, ,
```
which satisfies the boundary conditions
```math
    u |_{y=0} = 0 \, , \quad u |_{y=1} = f \, , \quad \text{and} \quad v |_{y=0} = v |_{y=1} = 0 \, .
```
The vorticity forcing is
```math
    F_{\omega} = - \xi^\prime f_x (y^3 - y^2 - 6y + 2) - f f_x (12 y^3 - 12 y^2 + 4 y) + f (y^3 - y^2 - 12 y + 4) \, ,
```
which implies that
```math
    F_v = \xi^\prime f_x (\tfrac{1}{4} y^4 - \tfrac{1}{3} y^3 - 3 y^2 + 2y)
        + f f_x (3 y^4 - 4 y^3 + 2y^2 ) 
        - f (\tfrac{1}{4} y^4 - \tfrac{1}{3} y^3 - 6 y^2 + 4 y) \, ,
```
and
```math
    F_v = 3 y^5 - 5 y^4 + 2y^3 \, .
```

We set up the problem in the same manner as the forced, free-slip problem above. Note that we 
also must the no-slip boundary condition ``u |_{y=0} = 0`` and the time-dependent fixed-slip 
condition ``u |_{y=1} = f``. As for the free-slip problem, we find that the error between the 
numerical and analytical solutions decreases with ``1 / N_x^2 \sim \Delta x^2``, where ``N_x``
is the number of grid points and ``\Delta x`` is the spatial resolution:

![Forced fixed slip convergence](https://raw.githubusercontent.com/CliMA/Oceananigans.jl/v0.35.0/docs/src/verification/convergence_plots/forced_fixed_slip_convergence.png)

The convergence tests are performed using both ``y`` and ``z`` as the bounded direction.
# Fractional step method

In some models (e.g., `NonhydrostaticModel` or `HydrostaticFreeSurfaceModel`) solving the momentum 
coupled with the continuity equation can be cumbersome so instead we employ a fractional step 
method. To approximate the solution of the coupled system we first solve an approximation to 
the discretized momentum equation for an intermediate velocity field ``\boldsymbol{v}^\star`` 
without worrying about satisfying the incompressibility constraint. We then project ``\boldsymbol{v}^\star`` 
onto the space of divergence-free velocity fields to obtain a value for ``\boldsymbol{v}^{n+1}`` 
that satisfies continuity.

For example, for the `NonhydrostaticModel`, if we ignore the background velocity fields and the
surface waves, we thus discretize the momentum equation as
```math
  \frac{\boldsymbol{v}^\star - \boldsymbol{v}^n}{\Delta t}
    = - \left[ \boldsymbol{v} \boldsymbol{\cdot} \boldsymbol{\nabla} \boldsymbol{v} \right]^{n+\frac{1}{2}}
      - \boldsymbol{f} \times \boldsymbol{v}^{n+\frac{1}{2}}
      + \boldsymbol{\nabla} \boldsymbol{\cdot} \left ( \nu \boldsymbol{\nabla} \boldsymbol{v}^{n+\frac{1}{2}} \right )
      + \boldsymbol{F}_{\boldsymbol{v}}^{n+\frac{1}{2}} \, ,
```
where the superscript ``n + \frac{1}{2}`` indicates that these terms are evaluated at time step 
``n + \frac{1}{2}``, which we compute explicitly (see [Time Stepping section](../numerical_implementation/time_stepping)).

The projection is then performed
```math
   \boldsymbol{v}^{n+1} = \boldsymbol{v}^\star - \Delta t \, \boldsymbol{\nabla} p^{n+1} \, ,
```
to obtain a divergence-free velocity field ``\boldsymbol{v}^{n+1}``. Here the projection is performed by solving an elliptic
problem for the pressure ``p^{n+1}`` with the boundary condition
```math
  \boldsymbol{\hat{n}} \boldsymbol{\cdot} \boldsymbol{\nabla} p^{n+1} |_{\partial\Omega} = 0 \, .
```

[Orszag86](@cite) and [Brown01](@cite) raise an important issue regarding these fractional step 
methods, which is that "*while the velocity can be reliably computed to second-order accuracy 
in time and space, the pressure is typically only first-order accurate in the ``L_\infty``-norm.*" 
The numerical boundary conditions must be carefully accounted for to ensure the second-order 
accuracy promised by the fractional step methods.
# Staggered grid

Velocities ``u``, ``v``, and ``w`` are defined on the faces of the cells, which are coincident with three orthogonal
coordinate axes (the Cartesian axes in the case of Oceananigans). Pressure ``p`` and tracers ``c`` are stored at
the cell  centers as cell averages. See schematic below of the different control
volumes. Other quantities may be defined at other locations. For example, vorticity ``\boldsymbol{\omega} = \boldsymbol{\nabla} \times \boldsymbol{v}``
is defined at the cell edges.[^1]

[^1]: In 2D it would more correct to say the cell corners. In 3D, variables like vorticity lie at the same vertical
    levels as the cell-centered variables and so they really lie at the cell edges.

![Schematic of control volumes](../numerical_implementation/assets/staggered_grid.png)
*A schematic of `Oceananigans.jl` finite volumes for a two-dimensional staggered grid in ``(x, z)``.
Tracers ``c`` and pressure ``p`` are defined at the center of the control volume. The ``u`` control volumes are 
centered on the left and right edges of the pressure control volume while the ``w`` control volumes are centered 
on the top and bottom edges of the pressure control volumes. The indexing convention places the ``i^{\rm{th}}`` 
``u``-node on cell ``x``-faces to the left of the ``i`` tracer point at cell centers.*

This staggered arrangement of variables is more complicated than the collocated grid arrangement but is greatly
beneficial as it avoids the odd-even decoupling between the pressure and velocity if they are stored at the same
positions. §6.1 of [Patankar80](@cite) discusses this problem in the presence of a zigzag pressure field: on a 1D
collocated grid the velocity at the point ``i`` is influenced by the pressure at points ``i-1`` and ``i+1``, and a zigzag
pressure field will be felt as a uniform pressure, which is obviously wrong and would reduce the accuracy of the
solution. The pressure is effectively taken from a coarser grid than what is actually used. The basic problem is that
the momentum equations will use the pressure difference between two alternate points when it should be using two
adjacent points.

From the viewpoint of linear algebra, these spurious pressure modes correspond to solutions in the null space of the
pressure projection operator with eigenvalue zero and are thus indistinguishable from a uniform pressure field
[Sani81](@cite).

The staggered grid was first introduced by [Harlow65](@cite) with their *marker and cell* method. In meteorology
and oceanography, this particular staggered grid configuration is referred to as the Arakawa C-grid after [Arakawa77](@cite), who
investigated four different staggered grids and the unstaggered A-grid for use in an atmospheric model.

[Arakawa77](@cite) investigated the dispersion relation of inertia-gravity waves[^2] traveling in the ``x``-direction
```math
  \omega^2 = f^2 + gHk^2 \, ,
```
in the linearized rotating shallow-water equations for five grids. Here ``\omega`` is the angular frequency, ``H`` is the
height of the fluid and ``k`` is the wavenumber in the ``x``-direction. Looking at the effect of spatial discretization
error on the frequency of these waves they find that the B and C-grids reproduce the dispersion relation most closely
out of the five [Arakawa77](@cite) (Figure 5). In particular, the dispersion relation for the C-grid is given by
```math
  \omega^2 = f^2 \left[ \cos^2 \left( \frac{k\Delta}{2} \right)
             + 4 \left( \frac{\lambda}{\Delta} \right)^2 \sin^2 \left( \frac{k\Delta}{2} \right) \right] \, ,
```
where ``\lambda`` is the wavelength and ``\Delta`` is the grid spacing. Paraphrasing p. 184 of [Arakawa77](@cite): The
wavelength of the shortest resolvable wave is ``2\Delta`` with corresponding wavenumber ``k = \pi/\Delta`` so it is
sufficient to evaluate the dispersion relation over the range ``0 < k \Delta < \pi``. The frequency is monotonically
increasing for ``\lambda / \Delta > \frac{1}{2}`` and monotonically decreasing for ``\lambda / \Delta < \frac{1}{2}``. For the
fourth smallest wave ``\lambda / \Delta = \frac{1}{2}`` we get ``\omega^2 = f^2`` which matches the ``k = 0`` wave. Furthermore,
the group velocity is zero for all ``k``. On the other grids, waves with ``k \Delta = \pi`` can behave like pure inertial
oscillations or stationary waves, which is bad.

The B and C-grids are less oscillatory than the others and quite faithfully simulate geostrophic adjustment. However,
the C-grid is the only one that faithfully reproduces the two-dimensional dispersion relation ``\omega^2(k, \ell)``, all
the other grids have false maxima, and so [Arakawa77](@cite) conclude that the C-grid is best for simulating geostrophic
adjustment except for abnormal situations in which ``\lambda / \Delta`` is less than or close to 1. This seems to have held
true for most atmospheric and oceanographic simulations as the C-grid is popular and widely used.

[^2]: Apparently also called Poincaré waves, Sverdrup waves, and *rotational gravity waves* §13.9 of [Kundu15](@cite).
# Library

Documenting the public user interface.

## Oceananigans.jl

```@autodocs
Modules = [Oceananigans]
Private = false
Pages   = [
    "Oceananigans.jl"
]
```

## Abstract operations

```@autodocs
Modules = [Oceananigans.AbstractOperations]
Private = false
```

## Advection

```@autodocs
Modules = [Oceananigans.Advection]
Private = false
```

## Architectures

```@autodocs
Modules = [Oceananigans.Architectures]
Private = false
Pages   = ["Architectures.jl"]
```

## Boundary conditions

```@autodocs
Modules = [Oceananigans.BoundaryConditions]
Private = false
```

## BuoyancyModels

```@autodocs
Modules = [Oceananigans.BuoyancyModels]
Private = false
```

## Coriolis

```@autodocs
Modules = [Oceananigans.Coriolis]
Private = false
```

## Diagnostics

```@autodocs
Modules = [Oceananigans.Diagnostics]
Private = false
```

## Fields

```@autodocs
Modules = [Oceananigans.Fields]
Private = false
```

## Forcings

```@autodocs
Modules = [Oceananigans.Forcings]
Private = false
```

## Grids

```@autodocs
Modules = [Oceananigans.Grids]
Private = false
```

## Immersed boundaries

```@autodocs
Modules = [Oceananigans.ImmersedBoundaries]
Private = false
```

## Lagrangian particle tracking

```@autodocs
Modules = [Oceananigans.LagrangianParticleTracking]
Private = false
```

## Logger

```@autodocs
Modules = [Oceananigans.Logger]
Private = false
Pages   = ["Logger.jl"]
```

## Models

```@autodocs
Modules = [Oceananigans.Models]
Private = false
```

### Non-hydrostatic models

```@autodocs
Modules = [Oceananigans.Models.NonhydrostaticModels]
Private = false
```

### Hydrostatic free-surface models

```@autodocs
Modules = [Oceananigans.Models.HydrostaticFreeSurfaceModels]
Private = false
]
```

### Shallow-water models

```@autodocs
Modules = [Oceananigans.Models.ShallowWaterModels]
Private = false
Pages   = [
    "Models/ShallowWaterModels/shallow_water_model.jl"
]
```

## Output readers

```@autodocs
Modules = [Oceananigans.OutputReaders]
Private = false
]
```

## Output writers

```@autodocs
Modules = [Oceananigans.OutputWriters]
Private = false
]
```

## Simulations

```@autodocs
Modules = [Oceananigans.Simulations]
Private = false
```

## Solvers

```@autodocs
Modules = [Oceananigans.Solvers]
Private = false
```

## Stokes drift

```@autodocs
Modules = [Oceananigans.StokesDrift]
Private = false
```

## Time steppers

```@autodocs
Modules = [Oceananigans.TimeSteppers]
Private = false
```

## Turbulence closures

```@autodocs
Modules = [Oceananigans.TurbulenceClosures]
Private = false
```

## Utilities

```@autodocs
Modules = [Oceananigans.Utils]
Private = false
```

# [Performance benchmarks](@id performance_benchmarks)

The performance benchmarking scripts in the
[`benchmarks`](https://github.com/CliMA/Oceananigans.jl/tree/main/benchmark)
directory of the git repository can be run to benchmark Oceananigans.jl on your machine.
They use [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) to collect data and [PrettyTables.jl](https://github.com/ronisbr/PrettyTables.jl)  to nicely
format the benchmark results.



## Shallow Water Model

This benchmark tests the performance of the shallow water model run in a doubly periodic domain (`topology = (Periodic, Periodic, Flat)`)
on a CPU versus a GPU.  We find that with the `WENO5` advection scheme we get a maximum speedup of more than 400 times on a `16384^2` grid.

```
Oceananigans v0.58.1
Julia Version 1.6.0
Commit f9720dc2eb (2021-03-24 12:55 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) Silver 4216 CPU @ 2.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, cascadelake)
Environment:
  EBVERSIONJULIA = 1.6.0
  JULIA_DEPOT_PATH = :
  EBROOTJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0
  EBDEVELJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0/easybuild/avx2-Core-julia-1.6.0-easybuild-devel
  JULIA_LOAD_PATH = :
  GPU: Tesla V100-SXM2-32GB

                                              Shallow water model benchmarks
┌───────────────┬─────────────┬───────┬────────────┬────────────┬────────────┬────────────┬───────────┬────────┬─────────┐
│ Architectures │ Float_types │    Ns │        min │     median │       mean │        max │    memory │ allocs │ samples │
├───────────────┼─────────────┼───────┼────────────┼────────────┼────────────┼────────────┼───────────┼────────┼─────────┤
│           CPU │     Float64 │    32 │   2.677 ms │   2.876 ms │   3.047 ms │   4.806 ms │  1.36 MiB │   2253 │      10 │
│           CPU │     Float64 │    64 │   5.795 ms │   5.890 ms │   6.073 ms │   7.770 ms │  1.36 MiB │   2255 │      10 │
│           CPU │     Float64 │   128 │  16.979 ms │  17.350 ms │  17.578 ms │  19.993 ms │  1.36 MiB │   2255 │      10 │
│           CPU │     Float64 │   256 │  62.543 ms │  63.222 ms │  63.544 ms │  67.347 ms │  1.36 MiB │   2255 │      10 │
│           CPU │     Float64 │   512 │ 250.149 ms │ 251.023 ms │ 251.092 ms │ 252.389 ms │  1.36 MiB │   2315 │      10 │
│           CPU │     Float64 │  1024 │ 990.901 ms │ 993.115 ms │ 993.360 ms │ 996.091 ms │  1.36 MiB │   2315 │       6 │
│           CPU │     Float64 │  2048 │    4.002 s │    4.004 s │    4.004 s │    4.007 s │  1.36 MiB │   2315 │       2 │
│           CPU │     Float64 │  4096 │   16.371 s │   16.371 s │   16.371 s │   16.371 s │  1.36 MiB │   2315 │       1 │
│           CPU │     Float64 │  8192 │   64.657 s │   64.657 s │   64.657 s │   64.657 s │  1.36 MiB │   2315 │       1 │
│           CPU │     Float64 │ 16384 │  290.423 s │  290.423 s │  290.423 s │  290.423 s │  1.36 MiB │   2315 │       1 │
│           GPU │     Float64 │    32 │   3.468 ms │   3.656 ms │   3.745 ms │   4.695 ms │  1.82 MiB │   5687 │      10 │
│           GPU │     Float64 │    64 │   3.722 ms │   3.903 ms │   4.050 ms │   5.671 ms │  1.82 MiB │   5687 │      10 │
│           GPU │     Float64 │   128 │   3.519 ms │   3.808 ms │   4.042 ms │   6.372 ms │  1.82 MiB │   5687 │      10 │
│           GPU │     Float64 │   256 │   3.822 ms │   4.153 ms │   4.288 ms │   5.810 ms │  1.82 MiB │   5687 │      10 │
│           GPU │     Float64 │   512 │   4.637 ms │   4.932 ms │   4.961 ms │   5.728 ms │  1.82 MiB │   5765 │      10 │
│           GPU │     Float64 │  1024 │   3.240 ms │   3.424 ms │   3.527 ms │   4.553 ms │  1.82 MiB │   5799 │      10 │
│           GPU │     Float64 │  2048 │  10.783 ms │  10.800 ms │  11.498 ms │  17.824 ms │  1.98 MiB │  16305 │      10 │
│           GPU │     Float64 │  4096 │  41.880 ms │  41.911 ms │  42.485 ms │  47.627 ms │  2.67 MiB │  61033 │      10 │
│           GPU │     Float64 │  8192 │ 166.751 ms │ 166.800 ms │ 166.847 ms │ 167.129 ms │  5.21 MiB │ 227593 │      10 │
│           GPU │     Float64 │ 16384 │ 681.129 ms │ 681.249 ms │ 681.301 ms │ 681.583 ms │ 16.59 MiB │ 973627 │       8 │
└───────────────┴─────────────┴───────┴────────────┴────────────┴────────────┴────────────┴───────────┴────────┴─────────┘

        Shallow water model CPU to GPU speedup
┌─────────────┬───────┬──────────┬─────────┬─────────┐
│ Float_types │    Ns │  speedup │  memory │  allocs │
├─────────────┼───────┼──────────┼─────────┼─────────┤
│     Float64 │    32 │ 0.786715 │ 1.33777 │ 2.52419 │
│     Float64 │    64 │  1.50931 │ 1.33774 │ 2.52195 │
│     Float64 │   128 │  4.55587 │ 1.33774 │ 2.52195 │
│     Float64 │   256 │  15.2238 │ 1.33774 │ 2.52195 │
│     Float64 │   512 │  50.8995 │ 1.33771 │ 2.49028 │
│     Float64 │  1024 │  290.085 │ 1.33809 │ 2.50497 │
│     Float64 │  2048 │  370.777 │ 1.45575 │  7.0432 │
│     Float64 │  4096 │  390.617 │ 1.95667 │ 26.3641 │
│     Float64 │  8192 │  387.632 │ 3.82201 │ 98.3123 │
│     Float64 │ 16384 │   426.31 │  12.177 │ 420.573 │
└─────────────┴───────┴──────────┴─────────┴─────────┘
```
As shown in the graph below, speedups increase sharply starting at grid size `512^2` and then plateau off at around 400 times at grid size `4096^2` and beyond.

![shallow_water_speedup](https://user-images.githubusercontent.com/45054739/128793049-7bcbabaa-2d66-4209-a311-b02729fb93fa.png)

The time graph below shows that execution times on GPU are negligebly small up until grid size `1024^2` where it starts to scale similarly to times on CPU.

![shallow_water_times](https://user-images.githubusercontent.com/45054739/128793311-e4bbfd5a-aea8-4cdc-bee8-cb71128ff5fe.png)



## Nonhydrostatic Model

Similar to to shallow water model, the nonhydrostatic model benchmark tests for its performance on both a CPU and a GPU. It was also benchmarked with the `WENO5` advection scheme. The nonhydrostatic model is 3-dimensional unlike the 2-dimensional shallow water model. Total number of grid points is Ns cubed.

```
Oceananigans v0.58.8
Julia Version 1.6.1
Commit 6aaedecc44 (2021-04-23 05:59 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) Silver 4216 CPU @ 2.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, cascadelake)
Environment:
  EBVERSIONJULIA = 1.6.1
  JULIA_DEPOT_PATH = :
  EBROOTJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.1
  EBDEVELJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.1/easybuild/avx2-Core-julia-1.6.1-easybuild-devel
  JULIA_LOAD_PATH = :
  GPU: Tesla V100-SXM2-32GB

                                            Nonhydrostatic model benchmarks
┌───────────────┬─────────────┬─────┬────────────┬────────────┬────────────┬────────────┬──────────┬────────┬─────────┐
│ Architectures │ Float_types │  Ns │        min │     median │       mean │        max │   memory │ allocs │ samples │
├───────────────┼─────────────┼─────┼────────────┼────────────┼────────────┼────────────┼──────────┼────────┼─────────┤
│           CPU │     Float32 │  32 │  34.822 ms │  34.872 ms │  35.278 ms │  38.143 ms │ 1.38 MiB │   2302 │      10 │
│           CPU │     Float32 │  64 │ 265.408 ms │ 265.571 ms │ 265.768 ms │ 267.765 ms │ 1.38 MiB │   2302 │      10 │
│           CPU │     Float32 │ 128 │    2.135 s │    2.135 s │    2.136 s │    2.138 s │ 1.38 MiB │   2302 │       3 │
│           CPU │     Float32 │ 256 │   17.405 s │   17.405 s │   17.405 s │   17.405 s │ 1.38 MiB │   2302 │       1 │
│           CPU │     Float64 │  32 │  37.022 ms │  37.179 ms │  37.335 ms │  39.017 ms │ 1.77 MiB │   2302 │      10 │
│           CPU │     Float64 │  64 │ 287.944 ms │ 288.154 ms │ 288.469 ms │ 290.838 ms │ 1.77 MiB │   2302 │      10 │
│           CPU │     Float64 │ 128 │    2.326 s │    2.326 s │    2.326 s │    2.327 s │ 1.77 MiB │   2302 │       3 │
│           CPU │     Float64 │ 256 │   19.561 s │   19.561 s │   19.561 s │   19.561 s │ 1.77 MiB │   2302 │       1 │
│           GPU │     Float32 │  32 │   4.154 ms │   4.250 ms │   4.361 ms │   5.557 ms │ 2.13 MiB │   6033 │      10 │
│           GPU │     Float32 │  64 │   3.383 ms │   3.425 ms │   3.889 ms │   8.028 ms │ 2.13 MiB │   6077 │      10 │
│           GPU │     Float32 │ 128 │   5.564 ms │   5.580 ms │   6.095 ms │  10.725 ms │ 2.15 MiB │   7477 │      10 │
│           GPU │     Float32 │ 256 │  38.685 ms │  38.797 ms │  39.548 ms │  46.442 ms │ 2.46 MiB │  27721 │      10 │
│           GPU │     Float64 │  32 │   3.309 ms │   3.634 ms │   3.802 ms │   5.844 ms │ 2.68 MiB │   6033 │      10 │
│           GPU │     Float64 │  64 │   3.330 ms │   3.648 ms │   4.008 ms │   7.808 ms │ 2.68 MiB │   6071 │      10 │
│           GPU │     Float64 │ 128 │   7.209 ms │   7.323 ms │   8.313 ms │  17.259 ms │ 2.71 MiB │   8515 │      10 │
│           GPU │     Float64 │ 256 │  46.614 ms │  56.444 ms │  55.461 ms │  56.563 ms │ 3.17 MiB │  38253 │      10 │
└───────────────┴─────────────┴─────┴────────────┴────────────┴────────────┴────────────┴──────────┴────────┴─────────┘

      Nonhydrostatic model CPU to GPU speedup
┌─────────────┬─────┬─────────┬─────────┬─────────┐
│ Float_types │  Ns │ speedup │  memory │  allocs │
├─────────────┼─────┼─────────┼─────────┼─────────┤
│     Float32 │  32 │ 8.20434 │ 1.53786 │ 2.62076 │
│     Float32 │  64 │ 77.5308 │ 1.53835 │ 2.63988 │
│     Float32 │ 128 │ 382.591 │ 1.55378 │ 3.24805 │
│     Float32 │ 256 │ 448.619 │ 1.77688 │ 12.0421 │
│     Float64 │  32 │ 10.2308 │ 1.51613 │ 2.62076 │
│     Float64 │  64 │ 78.9952 │ 1.51646 │ 2.63727 │
│     Float64 │ 128 │ 317.663 │ 1.53759 │ 3.69896 │
│     Float64 │ 256 │ 346.554 │ 1.79466 │ 16.6173 │
└─────────────┴─────┴─────────┴─────────┴─────────┘
```

Like the shallow water model, it can be seen at grid size `64^3` that the GPU is beginning to be saturated as speedups rapidly increase. At grid sizes `128^3` and `256^3` we see the speedup stablise to around 400 times.

![incompressible_speedup](https://user-images.githubusercontent.com/45054739/129825248-adb8dfe5-e9ea-4321-bd11-fb415d81e2cb.png)

For both float types, the benchmarked GPU times of the nonhydrostatic model starts to scale like its CPU times when grid size reaches `128^3`.

![incompressible_times](https://user-images.githubusercontent.com/45054739/129825253-0d5739d9-f0a7-476e-8152-4ee462b71ad5.png)



## Distributed Shallow Water Model

By using `MPI.jl` the shallow water model can be run on multiple CPUs and multiple GPUs. For the benchmark results shown below, each rank is run on one CPU core and each uses a distinct GPU if applicable. 

### Weak Scaling Shallow Water Model
```
Oceananigans v0.58.2
Julia Version 1.6.0
Commit f9720dc2eb (2021-03-24 12:55 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2683 v4 @ 2.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
Environment:
  EBVERSIONJULIA = 1.6.0
  JULIA_DEPOT_PATH = :
  EBROOTJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0
  EBDEVELJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0/easybuild/avx2-Core-julia-1.6.0-easybuild-devel
  JULIA_LOAD_PATH = :

                                  Shallow water model weak scaling benchmark
┌───────────────┬──────────┬────────────┬────────────┬────────────┬────────────┬──────────┬────────┬─────────┐
│          size │    ranks │        min │     median │       mean │        max │   memory │ allocs │ samples │
├───────────────┼──────────┼────────────┼────────────┼────────────┼────────────┼──────────┼────────┼─────────┤
│   (4096, 256) │   (1, 1) │ 363.885 ms │ 364.185 ms │ 364.911 ms │ 370.414 ms │ 1.60 MiB │   2774 │      10 │
│   (4096, 512) │   (1, 2) │ 370.782 ms │ 375.032 ms │ 375.801 ms │ 394.781 ms │ 1.49 MiB │   3116 │      20 │
│  (4096, 1024) │   (1, 4) │ 369.648 ms │ 369.973 ms │ 371.613 ms │ 399.526 ms │ 1.49 MiB │   3116 │      40 │
│  (4096, 2048) │   (1, 8) │ 377.386 ms │ 379.982 ms │ 382.732 ms │ 432.787 ms │ 1.49 MiB │   3116 │      80 │
│  (4096, 4096) │  (1, 16) │ 388.336 ms │ 395.473 ms │ 400.079 ms │ 496.598 ms │ 1.49 MiB │   3116 │     160 │
│  (4096, 8192) │  (1, 32) │ 403.565 ms │ 447.136 ms │ 449.138 ms │ 545.945 ms │ 1.49 MiB │   3116 │     320 │
│ (4096, 16384) │  (1, 64) │ 397.965 ms │ 441.627 ms │ 453.465 ms │ 619.493 ms │ 1.49 MiB │   3125 │     640 │
│ (4096, 32768) │ (1, 128) │ 400.481 ms │ 447.789 ms │ 448.692 ms │ 590.028 ms │ 1.49 MiB │   3125 │    1280 │
└───────────────┴──────────┴────────────┴────────────┴────────────┴────────────┴──────────┴────────┴─────────┘

                Shallow water model weak scaling speedup
┌───────────────┬──────────┬──────────┬────────────┬──────────┬─────────┐
│          size │    ranks │ slowdown │ efficiency │   memory │  allocs │
├───────────────┼──────────┼──────────┼────────────┼──────────┼─────────┤
│   (4096, 256) │   (1, 1) │      1.0 │        1.0 │      1.0 │     1.0 │
│   (4096, 512) │   (1, 2) │  1.02978 │   0.971077 │ 0.930602 │ 1.12329 │
│  (4096, 1024) │   (1, 4) │  1.01589 │   0.984355 │ 0.930602 │ 1.12329 │
│  (4096, 2048) │   (1, 8) │  1.04338 │   0.958427 │ 0.930602 │ 1.12329 │
│  (4096, 4096) │  (1, 16) │  1.08591 │   0.920886 │ 0.930602 │ 1.12329 │
│  (4096, 8192) │  (1, 32) │  1.22777 │   0.814484 │ 0.930602 │ 1.12329 │
│ (4096, 16384) │  (1, 64) │  1.21264 │   0.824644 │ 0.930687 │ 1.12653 │
│ (4096, 32768) │ (1, 128) │  1.22957 │   0.813296 │ 0.930687 │ 1.12653 │
└───────────────┴──────────┴──────────┴────────────┴──────────┴─────────┘
```

As seen in the tables above and in the graph below, efficiency drops off to around 80% and remains as such from 16 to 128 ranks. GPUs are not used in this or the next benchmark setup. 

![ws_shallow_water_efficiency](https://user-images.githubusercontent.com/45054739/129826042-6ed4345b-b53a-49af-b375-6b7f11f53f31.png)

### Strong Scaling Shallow Water Model
```
Oceananigans v0.58.2
Julia Version 1.6.0
Commit f9720dc2eb (2021-03-24 12:55 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2683 v4 @ 2.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
Environment:
  EBVERSIONJULIA = 1.6.0
  JULIA_DEPOT_PATH = :
  EBROOTJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0
  EBDEVELJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0/easybuild/avx2-Core-julia-1.6.0-easybuild-devel
  JULIA_LOAD_PATH = :

                                Shallow water model strong scaling benchmark
┌──────────────┬──────────┬────────────┬────────────┬────────────┬────────────┬──────────┬────────┬─────────┐
│         size │    ranks │        min │     median │       mean │        max │   memory │ allocs │ samples │
├──────────────┼──────────┼────────────┼────────────┼────────────┼────────────┼──────────┼────────┼─────────┤
│ (4096, 4096) │   (1, 1) │    5.694 s │    5.694 s │    5.694 s │    5.694 s │ 1.60 MiB │   2804 │       1 │
│ (4096, 4096) │   (1, 2) │    2.865 s │    2.865 s │    2.866 s │    2.869 s │ 1.49 MiB │   3146 │       4 │
│ (4096, 4096) │   (1, 4) │    1.435 s │    1.437 s │    1.441 s │    1.475 s │ 1.49 MiB │   3146 │      16 │
│ (4096, 4096) │   (1, 8) │ 732.711 ms │ 736.394 ms │ 738.930 ms │ 776.773 ms │ 1.49 MiB │   3146 │      56 │
│ (4096, 4096) │  (1, 16) │ 389.211 ms │ 395.749 ms │ 396.813 ms │ 433.332 ms │ 1.49 MiB │   3116 │     160 │
│ (4096, 4096) │  (1, 32) │ 197.894 ms │ 219.211 ms │ 236.780 ms │ 367.188 ms │ 1.49 MiB │   3116 │     320 │
│ (4096, 4096) │  (1, 64) │ 101.520 ms │ 112.606 ms │ 116.809 ms │ 221.497 ms │ 1.49 MiB │   3125 │     640 │
│ (4096, 4096) │ (1, 128) │  51.452 ms │  60.256 ms │  70.959 ms │ 232.309 ms │ 1.49 MiB │   3125 │    1280 │
└──────────────┴──────────┴────────────┴────────────┴────────────┴────────────┴──────────┴────────┴─────────┘

              Shallow water model strong scaling speedup
┌──────────────┬──────────┬─────────┬────────────┬──────────┬─────────┐
│         size │    ranks │ speedup │ efficiency │   memory │  allocs │
├──────────────┼──────────┼─────────┼────────────┼──────────┼─────────┤
│ (4096, 4096) │   (1, 1) │     1.0 │        1.0 │      1.0 │     1.0 │
│ (4096, 4096) │   (1, 2) │ 1.98728 │   0.993641 │ 0.930621 │ 1.12197 │
│ (4096, 4096) │   (1, 4) │ 3.96338 │   0.990845 │ 0.930621 │ 1.12197 │
│ (4096, 4096) │   (1, 8) │ 7.73237 │   0.966547 │ 0.930621 │ 1.12197 │
│ (4096, 4096) │  (1, 16) │ 14.3881 │   0.899255 │ 0.930336 │ 1.11127 │
│ (4096, 4096) │  (1, 32) │ 25.9754 │   0.811731 │ 0.930336 │ 1.11127 │
│ (4096, 4096) │  (1, 64) │ 50.5666 │   0.790102 │ 0.930421 │ 1.11448 │
│ (4096, 4096) │ (1, 128) │ 94.4984 │   0.738269 │ 0.930421 │ 1.11448 │
└──────────────┴──────────┴─────────┴────────────┴──────────┴─────────┘
```

Slightly differing from the weak scaling results, efficiencies drop below 80% to around 74% at 128 ranks for the strong scaling distributed shallow water model benchmark. This is likely caused by the 128 CPU cores not being sufficiently saturated anymore by the constant `4096^2` grid size thus losing some efficiency overheads.

![ss_shallow_water_efficiency](https://user-images.githubusercontent.com/45054739/129826134-3c526b9f-efd1-436c-9dc1-bde376a035db.png)

### Multi-GPU Shallow Water Model

While still a work in progress, it is possible to use CUDA-aware MPI to run the shallow water model on multiple GPUs. Though efficiencies may not be as high as multi-CPU, the multi-GPU architecture is still worthwhile when keeping in mind the baseline speedups generated by using a single GPU. Note that though it is possible for multiple ranks to share the use of a single GPU, efficiencies would significantly decrease and memory may be insufficient. The results below show up to three ranks each using a separate GPU.

```
Julia Version 1.6.2
Commit 1b93d53fc4 (2021-07-14 15:36 UTC)
Platform Info:
  OS: Linux (powerpc64le-unknown-linux-gnu)
  CPU: unknown
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, pwr9)
Environment:
  JULIA_MPI_PATH = /home/software/spack/openmpi/3.1.4-nhjzelonyovxks5ydtrxehceqxsbf7ik
  JULIA_CUDA_USE_BINARYBUILDER = false
  JULIA_DEPOT_PATH = /nobackup/users/henryguo/projects/henry-test/Oceananigans.jl/benchmark/.julia
  GPU: Tesla V100-SXM2-32GB
  
                              Shallow water model weak scaling benchmark
┌───────────────┬──────────┬─────────┬─────────┬─────────┬─────────┬──────────┬────────┬─────────┐
│          size │    ranks │     min │  median │    mean │     max │   memory │ allocs │ samples │
├───────────────┼──────────┼─────────┼─────────┼─────────┼─────────┼──────────┼────────┼─────────┤
│   (4096, 256) │   (1, 1) │ 2.702 ms│ 2.728 ms│ 2.801 ms│ 3.446 ms│ 2.03 MiB │   5535 │      10 │
│   (4096, 512) │   (1, 2) │ 3.510 ms│ 3.612 ms│ 4.287 ms│16.546 ms│ 2.03 MiB │   5859 │      20 │
│   (4096, 768) │   (1, 3) │ 3.553 ms│ 3.653 ms│ 5.195 ms│39.152 ms│ 2.03 MiB │   5859 │      30 │
└───────────────┴──────────┴─────────┴─────────┴─────────┴─────────┴──────────┴────────┴─────────┘

                Shallow water model weak scaling speedup
┌───────────────┬──────────┬──────────┬────────────┬──────────┬─────────┐
│          size │    ranks │ slowdown │ efficiency │   memory │  allocs │
├───────────────┼──────────┼──────────┼────────────┼──────────┼─────────┤
│   (4096, 256) │   (1, 1) │      1.0 │        1.0 │      1.0 │     1.0 │
│   (4096, 512) │   (1, 2) │  1.32399 │   0.755293 │  1.00271 │ 1.05854 │
│   (4096, 768) │   (1, 3) │  1.33901 │   0.746818 │  1.00271 │ 1.05854 │
└───────────────┴──────────┴──────────┴────────────┴──────────┴─────────┘
```



## Distributed Nonhydrostatic Model

Similar to the distributed shallow water model benchmark results shown above, the distributed nonhydrostatic model was also benchmarked with the strong and weak scaling methods.

### Weak Scaling Nonhydrostatic Model

Weak scaling efficiencies can be improved for the nonhydrostatic model.

```
Oceananigans v0.60.1
Julia Version 1.6.1
Commit 6aaedecc44 (2021-04-23 05:59 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2683 v4 @ 2.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
Environment:
  JULIA_MPI_PATH = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/intel2020/openmpi/4.0.3
  EBVERSIONJULIA = 1.6.1
  JULIA_DEPOT_PATH = :
  JULIA_MPI_BINARY = system
  JULIA_MPI_LIBRARY = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/intel2020/openmpi/4.0.3/lib64/libmpi.so
  EBROOTJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.1
  JULIA_MPI_ABI = OpenMPI
  EBDEVELJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.1/easybuild/avx2-Core-julia-1.6.1-easybuild-devel
  JULIA_LOAD_PATH = :
  JULIA_MPIEXEC = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/intel2020/openmpi/4.0.3/bin/mpiexec

                                    Nonhydrostatic model weak scaling benchmark
┌──────────────────┬─────────────┬────────────┬────────────┬────────────┬────────────┬──────────┬────────┬─────────┐
│             size │       ranks │        min │     median │       mean │        max │   memory │ allocs │ samples │
├──────────────────┼─────────────┼────────────┼────────────┼────────────┼────────────┼──────────┼────────┼─────────┤
│   (128, 128, 16) │   (1, 1, 1) │  33.568 ms │  34.087 ms │  34.173 ms │  34.894 ms │ 2.05 MiB │   2762 │      10 │
│   (128, 128, 32) │   (1, 2, 1) │  36.650 ms │  37.161 ms │  37.393 ms │  42.411 ms │ 1.99 MiB │   3096 │      20 │
│   (128, 128, 64) │   (1, 4, 1) │  41.861 ms │  43.440 ms │  46.176 ms │  97.578 ms │ 1.99 MiB │   3136 │      40 │
│  (128, 128, 128) │   (1, 8, 1) │  59.995 ms │  64.110 ms │  68.021 ms │ 138.422 ms │ 1.99 MiB │   3216 │      80 │
│  (128, 128, 256) │  (1, 16, 1) │  62.633 ms │  71.266 ms │  74.775 ms │ 164.206 ms │ 2.01 MiB │   3376 │     160 │
│  (128, 128, 512) │  (1, 32, 1) │ 108.253 ms │ 135.611 ms │ 139.384 ms │ 225.336 ms │ 2.04 MiB │   3722 │     320 │
│ (128, 128, 1024) │  (1, 64, 1) │ 138.504 ms │ 181.043 ms │ 186.386 ms │ 335.170 ms │ 2.12 MiB │   4372 │     640 │
│ (128, 128, 2048) │ (1, 128, 1) │ 218.592 ms │ 285.293 ms │ 290.989 ms │ 434.878 ms │ 2.39 MiB │   5652 │    1280 │
└──────────────────┴─────────────┴────────────┴────────────┴────────────┴────────────┴──────────┴────────┴─────────┘

                   Nonhydrostatic model weak scaling speedup
┌──────────────────┬─────────────┬──────────┬────────────┬──────────┬─────────┐
│             size │       ranks │  speedup │ efficiency │   memory │  allocs │
├──────────────────┼─────────────┼──────────┼────────────┼──────────┼─────────┤
│   (128, 128, 16) │   (1, 1, 1) │      1.0 │        1.0 │      1.0 │     1.0 │
│   (128, 128, 32) │   (1, 2, 1) │ 0.917292 │   0.917292 │ 0.968543 │ 1.12093 │
│   (128, 128, 64) │   (1, 4, 1) │ 0.784698 │   0.784698 │ 0.969719 │ 1.13541 │
│  (128, 128, 128) │   (1, 8, 1) │ 0.531697 │   0.531697 │ 0.972279 │ 1.16437 │
│  (128, 128, 256) │  (1, 16, 1) │ 0.478315 │   0.478315 │ 0.978143 │  1.2223 │
│  (128, 128, 512) │  (1, 32, 1) │ 0.251361 │   0.251361 │ 0.992878 │ 1.34757 │
│ (128, 128, 1024) │  (1, 64, 1) │ 0.188283 │   0.188283 │  1.03539 │ 1.58291 │
│ (128, 128, 2048) │ (1, 128, 1) │ 0.119482 │   0.119482 │  1.16791 │ 2.04634 │
└──────────────────┴─────────────┴──────────┴────────────┴──────────┴─────────┘
```

![ws_nonhydrostatic_efficiency](https://user-images.githubusercontent.com/45054739/130146112-2dd7e24a-7a79-4000-a899-405362af0f2a.png)

### Strong Scaling Nonhydrostatic Model

Strong scaling efficiencies can also be improved for the nonhydrostatic model.

```
Oceananigans v0.60.1
Julia Version 1.6.1
Commit 6aaedecc44 (2021-04-23 05:59 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2683 v4 @ 2.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
Environment:
  JULIA_MPI_PATH = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/intel2020/openmpi/4.0.3
  EBVERSIONJULIA = 1.6.1
  JULIA_DEPOT_PATH = :
  JULIA_MPI_BINARY = system
  JULIA_MPI_LIBRARY = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/intel2020/openmpi/4.0.3/lib64/libmpi.so
  EBROOTJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.1
  JULIA_MPI_ABI = OpenMPI
  EBDEVELJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.1/easybuild/avx2-Core-julia-1.6.1-easybuild-devel
  JULIA_LOAD_PATH = :
  JULIA_MPIEXEC = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/intel2020/openmpi/4.0.3/bin/mpiexec

                                   Nonhydrostatic model strong scaling benchmark
┌─────────────────┬─────────────┬────────────┬────────────┬────────────┬────────────┬──────────┬────────┬─────────┐
│            size │       ranks │        min │     median │       mean │        max │   memory │ allocs │ samples │
├─────────────────┼─────────────┼────────────┼────────────┼────────────┼────────────┼──────────┼────────┼─────────┤
│ (256, 256, 256) │   (1, 1, 1) │    3.049 s │    3.053 s │    3.053 s │    3.057 s │ 2.05 MiB │   2762 │       2 │
│ (256, 256, 256) │   (1, 2, 1) │    1.609 s │    1.610 s │    1.611 s │    1.620 s │ 1.99 MiB │   3096 │       8 │
│ (256, 256, 256) │   (1, 4, 1) │ 814.290 ms │ 817.305 ms │ 818.685 ms │ 833.792 ms │ 1.99 MiB │   3136 │      28 │
│ (256, 256, 256) │   (1, 8, 1) │ 434.521 ms │ 439.352 ms │ 443.049 ms │ 508.913 ms │ 1.99 MiB │   3216 │      80 │
│ (256, 256, 256) │  (1, 16, 1) │ 251.632 ms │ 272.364 ms │ 277.555 ms │ 370.059 ms │ 2.01 MiB │   3376 │     160 │
│ (256, 256, 256) │  (1, 32, 1) │ 182.380 ms │ 233.322 ms │ 247.325 ms │ 441.971 ms │ 2.04 MiB │   3696 │     320 │
│ (256, 256, 256) │  (1, 64, 1) │ 119.546 ms │ 178.933 ms │ 204.036 ms │ 564.097 ms │ 2.12 MiB │   4346 │     640 │
│ (256, 256, 256) │ (1, 128, 1) │  73.802 ms │ 120.147 ms │ 136.395 ms │ 378.697 ms │ 2.39 MiB │   5626 │    1280 │
└─────────────────┴─────────────┴────────────┴────────────┴────────────┴────────────┴──────────┴────────┴─────────┘

                 Nonhydrostatic model strong scaling speedup
┌─────────────────┬─────────────┬─────────┬────────────┬──────────┬─────────┐
│            size │       ranks │ speedup │ efficiency │   memory │  allocs │
├─────────────────┼─────────────┼─────────┼────────────┼──────────┼─────────┤
│ (256, 256, 256) │   (1, 1, 1) │     1.0 │        1.0 │      1.0 │     1.0 │
│ (256, 256, 256) │   (1, 2, 1) │ 1.89655 │   0.948276 │ 0.968543 │ 1.12093 │
│ (256, 256, 256) │   (1, 4, 1) │ 3.73522 │   0.933804 │ 0.969719 │ 1.13541 │
│ (256, 256, 256) │   (1, 8, 1) │ 6.94845 │   0.868556 │ 0.972279 │ 1.16437 │
│ (256, 256, 256) │  (1, 16, 1) │ 11.2086 │   0.700536 │ 0.978143 │  1.2223 │
│ (256, 256, 256) │  (1, 32, 1) │ 13.0841 │   0.408879 │ 0.992685 │ 1.33816 │
│ (256, 256, 256) │  (1, 64, 1) │ 17.0612 │   0.266582 │  1.03519 │  1.5735 │
│ (256, 256, 256) │ (1, 128, 1) │  25.409 │   0.198508 │  1.16772 │ 2.03693 │
└─────────────────┴─────────────┴─────────┴────────────┴──────────┴─────────┘
```

![ss_nonhydrostatic_efficiency](https://user-images.githubusercontent.com/45054739/130146219-b354fa25-7d77-4206-8e7e-ec639b2250fa.png)



## Multithreading

Oceananigans can also achieve parallelism via multithreading. Though its efficiencies are less than that of the MPI distributed architectures, its simple setup still makes it a viable option for achieving speedups on simple systems.

### Weak Scaling Multithreaded Shallow Water Model

The initial drop and then rise in efficiencies going from 1 to 2 to 4 threads is likely caused by the 2 threads being automatically allocated onto only one physical CPU core. Though one physical CPU core may contain 2 logical cores each capable of running a separate thread, having 2 threads run on one core will still reduce efficiencies as many resources such as caches and buses must be shared by both threads. Note that there are as many CPU cores allocated as the maximum number of threads.

```
Oceananigans v0.58.9
Julia Version 1.6.0
Commit f9720dc2eb (2021-03-24 12:55 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2683 v4 @ 2.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
Environment:
  EBVERSIONJULIA = 1.6.0
  JULIA_DEPOT_PATH = :
  EBROOTJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0
  EBDEVELJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0/easybuild/avx2-Core-julia-1.6.0-easybuild-devel
  JULIA_LOAD_PATH = :

                  Shallow water model weak scaling with multithreading benchmark
┌───────────────┬─────────┬─────────┬─────────┬─────────┬─────────┬───────────┬─────────┬─────────┐
│          size │ threads │     min │  median │    mean │     max │    memory │  allocs │ samples │
├───────────────┼─────────┼─────────┼─────────┼─────────┼─────────┼───────────┼─────────┼─────────┤
│   (8192, 512) │       1 │ 1.458 s │ 1.458 s │ 1.458 s │ 1.458 s │  1.37 MiB │    2318 │       4 │
│  (8192, 1024) │       2 │ 2.925 s │ 2.989 s │ 2.989 s │ 3.052 s │ 18.06 MiB │ 1076944 │       2 │
│  (8192, 2048) │       4 │ 2.296 s │ 2.381 s │ 2.397 s │ 2.515 s │ 13.60 MiB │  760190 │       3 │
│  (8192, 4096) │       8 │ 2.347 s │ 2.369 s │ 2.377 s │ 2.415 s │ 16.36 MiB │  891860 │       3 │
│  (8192, 8192) │      16 │ 2.407 s │ 2.548 s │ 2.517 s │ 2.595 s │ 17.44 MiB │  863941 │       3 │
│ (8192, 16384) │      32 │ 3.023 s │ 3.069 s │ 3.069 s │ 3.115 s │ 23.03 MiB │ 1034063 │       2 │
└───────────────┴─────────┴─────────┴─────────┴─────────┴─────────┴───────────┴─────────┴─────────┘

        Shallow water model weak multithreading scaling speedup
┌───────────────┬─────────┬──────────┬────────────┬─────────┬─────────┐
│          size │ threads │ slowdown │ efficiency │  memory │  allocs │
├───────────────┼─────────┼──────────┼────────────┼─────────┼─────────┤
│   (8192, 512) │       1 │      1.0 │        1.0 │     1.0 │     1.0 │
│  (8192, 1024) │       2 │  2.04972 │   0.487872 │ 13.2156 │ 464.601 │
│  (8192, 2048) │       4 │  1.63302 │   0.612363 │ 9.95278 │ 327.951 │
│  (8192, 4096) │       8 │  1.62507 │   0.615359 │ 11.9706 │ 384.754 │
│  (8192, 8192) │      16 │  1.74747 │   0.572257 │  12.755 │  372.71 │
│ (8192, 16384) │      32 │  2.10486 │    0.47509 │  16.846 │ 446.101 │
└───────────────┴─────────┴──────────┴────────────┴─────────┴─────────┘
```

### Strong Scaling Multithreaded Nonhydrostatic Model

The notable and continuous decrease in efficiencies for the strong scaling nonhydrostatic model is likely caused by the `256^3` grid not sufficiently saturating 32 threads running on 32 CPUs. At the time this benchmark was produced, multithreading for both nonhydrostatic and shallow water models is still an active area of improvement. Please use the appropriate scripts found in [`benchmarks`](https://github.com/CliMA/Oceananigans.jl/tree/main/benchmark) to obtain more recent and hopefully ameliorated benchmark results.

```
Oceananigans v0.58.9
Julia Version 1.6.1
Commit 6aaedecc44 (2021-04-23 05:59 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2683 v4 @ 2.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
Environment:
  EBVERSIONJULIA = 1.6.1
  JULIA_DEPOT_PATH = :
  EBROOTJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.1
  EBDEVELJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.1/easybuild/avx2-Core-julia-1.6.1-easybuild-devel
  JULIA_LOAD_PATH = :

                                     Multithreading benchmarks
┌──────┬─────────┬────────────┬────────────┬────────────┬────────────┬──────────┬────────┬─────────┐
│ size │ threads │        min │     median │       mean │        max │   memory │ allocs │ samples │
├──────┼─────────┼────────────┼────────────┼────────────┼────────────┼──────────┼────────┼─────────┤
│  256 │       1 │    2.496 s │    2.637 s │    2.637 s │    2.777 s │ 1.70 MiB │   2251 │       2 │
│  256 │       2 │    2.385 s │    2.618 s │    2.618 s │    2.851 s │ 7.03 MiB │ 342397 │       2 │
│  256 │       4 │    1.320 s │    1.321 s │    1.333 s │    1.371 s │ 3.69 MiB │ 113120 │       4 │
│  256 │       8 │ 850.438 ms │ 855.292 ms │ 855.952 ms │ 861.966 ms │ 3.31 MiB │  65709 │       6 │
│  256 │      16 │ 642.225 ms │ 645.458 ms │ 648.685 ms │ 674.259 ms │ 3.60 MiB │  40992 │       8 │
│  256 │      32 │ 680.938 ms │ 694.376 ms │ 701.272 ms │ 746.599 ms │ 4.88 MiB │  36729 │       8 │
└──────┴─────────┴────────────┴────────────┴────────────┴────────────┴──────────┴────────┴─────────┘

     Nonhydrostatic Strong Scaling Multithreading speedup
┌──────┬─────────┬──────────┬────────────┬─────────┬─────────┐
│ size │ threads │ slowdown │ efficiency │  memory │  allocs │
├──────┼─────────┼──────────┼────────────┼─────────┼─────────┤
│  256 │       1 │      1.0 │        1.0 │     1.0 │     1.0 │
│  256 │       2 │ 0.992966 │   0.503542 │ 4.14014 │ 152.109 │
│  256 │       4 │ 0.501089 │   0.498913 │ 2.17724 │ 50.2532 │
│  256 │       8 │ 0.324366 │   0.385367 │ 1.94899 │  29.191 │
│  256 │      16 │ 0.244788 │   0.255323 │ 2.12262 │ 18.2106 │
│  256 │      32 │ 0.263339 │   0.118668 │ 2.87624 │ 16.3167 │
└──────┴─────────┴──────────┴────────────┴─────────┴─────────┘
```



## Tracers

This benchmark tests the performance impacts of running with various amounts of active
and passive tracers and compares the difference in speedup going from CPU to GPU. Number of tracers are listed in the tracers column as (active, passive). 

```
Oceananigans v0.58.1
Julia Version 1.6.0
Commit f9720dc2eb (2021-03-24 12:55 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) Silver 4216 CPU @ 2.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, cascadelake)
Environment:
  EBVERSIONJULIA = 1.6.0
  JULIA_DEPOT_PATH = :
  EBROOTJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0
  EBDEVELJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0/easybuild/avx2-Core-julia-1.6.0-easybuild-devel
  JULIA_LOAD_PATH = :
  GPU: Tesla V100-SXM2-32GB

                                       Arbitrary tracers benchmarks
┌───────────────┬─────────┬───────────┬───────────┬───────────┬───────────┬────────────┬────────┬─────────┐
│ Architectures │ tracers │       min │    median │      mean │       max │     memory │ allocs │ samples │
├───────────────┼─────────┼───────────┼───────────┼───────────┼───────────┼────────────┼────────┼─────────┤
│           CPU │  (0, 0) │   1.439 s │   1.440 s │   1.440 s │   1.441 s │ 908.03 KiB │   1656 │       4 │
│           CPU │  (0, 1) │   1.539 s │   1.574 s │   1.575 s │   1.613 s │   1.24 MiB │   1942 │       4 │
│           CPU │  (0, 2) │   1.668 s │   1.669 s │   1.670 s │   1.671 s │   1.76 MiB │   2291 │       3 │
│           CPU │  (1, 0) │   1.527 s │   1.532 s │   1.532 s │   1.536 s │   1.24 MiB │   1942 │       4 │
│           CPU │  (2, 0) │   1.690 s │   1.697 s │   1.695 s │   1.698 s │   1.77 MiB │   2301 │       3 │
│           CPU │  (2, 3) │   2.234 s │   2.239 s │   2.241 s │   2.251 s │   3.59 MiB │   3928 │       3 │
│           CPU │  (2, 5) │   2.755 s │   2.838 s │   2.838 s │   2.921 s │   5.18 MiB │   4908 │       2 │
│           CPU │ (2, 10) │   3.588 s │   3.748 s │   3.748 s │   3.908 s │  10.39 MiB │   7682 │       2 │
│           GPU │  (0, 0) │  9.702 ms │ 12.755 ms │ 12.458 ms │ 12.894 ms │   1.59 MiB │  12321 │      10 │
│           GPU │  (0, 1) │ 13.863 ms │ 13.956 ms │ 14.184 ms │ 16.297 ms │   2.20 MiB │  14294 │      10 │
│           GPU │  (0, 2) │ 15.166 ms │ 15.230 ms │ 15.700 ms │ 19.893 ms │   2.93 MiB │  15967 │      10 │
│           GPU │  (1, 0) │ 13.740 ms │ 13.838 ms │ 14.740 ms │ 22.940 ms │   2.20 MiB │  14278 │      10 │
│           GPU │  (2, 0) │ 15.103 ms │ 15.199 ms │ 16.265 ms │ 25.906 ms │   2.93 MiB │  15913 │      10 │
│           GPU │  (2, 3) │ 13.981 ms │ 18.856 ms │ 18.520 ms │ 20.519 ms │   5.56 MiB │  17974 │      10 │
│           GPU │  (2, 5) │ 15.824 ms │ 21.211 ms │ 21.064 ms │ 24.897 ms │   7.86 MiB │  23938 │      10 │
│           GPU │ (2, 10) │ 22.085 ms │ 27.236 ms │ 28.231 ms │ 38.295 ms │  15.02 MiB │  31086 │      10 │
└───────────────┴─────────┴───────────┴───────────┴───────────┴───────────┴────────────┴────────┴─────────┘

  Arbitrary tracers CPU to GPU speedup
┌─────────┬─────────┬─────────┬─────────┐
│ tracers │ speedup │  memory │  allocs │
├─────────┼─────────┼─────────┼─────────┤
│  (0, 0) │ 112.881 │ 1.78792 │ 7.44022 │
│  (0, 1) │ 112.761 │ 1.77743 │ 7.36045 │
│  (0, 2) │ 109.618 │  1.6627 │ 6.96945 │
│  (1, 0) │ 110.717 │ 1.77723 │ 7.35221 │
│  (2, 0) │ 111.678 │ 1.66267 │ 6.91569 │
│  (2, 3) │ 118.737 │ 1.55043 │ 4.57587 │
│  (2, 5) │ 133.803 │  1.5155 │ 4.87734 │
│ (2, 10) │ 137.615 │ 1.44535 │  4.0466 │
└─────────┴─────────┴─────────┴─────────┘

       Arbitrary tracers relative performance (CPU)
┌───────────────┬─────────┬──────────┬─────────┬─────────┐
│ Architectures │ tracers │ slowdown │  memory │  allocs │
├───────────────┼─────────┼──────────┼─────────┼─────────┤
│           CPU │  (0, 0) │      1.0 │     1.0 │     1.0 │
│           CPU │  (0, 1) │  1.09293 │ 1.39873 │ 1.17271 │
│           CPU │  (0, 2) │  1.15948 │ 1.99019 │ 1.38345 │
│           CPU │  (1, 0) │  1.06409 │ 1.39873 │ 1.17271 │
│           CPU │  (2, 0) │  1.17887 │ 1.99054 │ 1.38949 │
│           CPU │  (2, 3) │  1.55493 │ 4.04677 │ 2.37198 │
│           CPU │  (2, 5) │  1.97115 │ 5.84537 │ 2.96377 │
│           CPU │ (2, 10) │   2.6031 │ 11.7179 │ 4.63889 │
└───────────────┴─────────┴──────────┴─────────┴─────────┘

       Arbitrary tracers relative performance (GPU)
┌───────────────┬─────────┬──────────┬─────────┬─────────┐
│ Architectures │ tracers │ slowdown │  memory │  allocs │
├───────────────┼─────────┼──────────┼─────────┼─────────┤
│           GPU │  (0, 0) │      1.0 │     1.0 │     1.0 │
│           GPU │  (0, 1) │   1.0941 │ 1.39053 │ 1.16013 │
│           GPU │  (0, 2) │  1.19399 │ 1.85081 │ 1.29592 │
│           GPU │  (1, 0) │  1.08489 │ 1.39037 │ 1.15883 │
│           GPU │  (2, 0) │  1.19157 │ 1.85109 │ 1.29153 │
│           GPU │  (2, 3) │  1.47824 │ 3.50924 │ 1.45881 │
│           GPU │  (2, 5) │  1.66293 │ 4.95474 │ 1.94286 │
│           GPU │ (2, 10) │  2.13524 │ 9.47276 │ 2.52301 │
└───────────────┴─────────┴──────────┴─────────┴─────────┘
```



## Turbulence closures

This benchmark tests the performance impacts of various turbulent diffusivity closures
and large eddy simulation (LES) models as well as how much speedup they experience going from CPU to GPU.

```
Oceananigans v0.58.1
Julia Version 1.6.0
Commit f9720dc2eb (2021-03-24 12:55 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) Silver 4216 CPU @ 2.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, cascadelake)
Environment:
  EBVERSIONJULIA = 1.6.0
  JULIA_DEPOT_PATH = :
  EBROOTJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0
  EBDEVELJULIA = /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/julia/1.6.0/easybuild/avx2-Core-julia-1.6.0-easybuild-devel
  JULIA_LOAD_PATH = :
  GPU: Tesla V100-SXM2-32GB

                                                  Turbulence closure benchmarks
┌───────────────┬──────────────────────────────────┬───────────┬───────────┬───────────┬───────────┬──────────┬────────┬─────────┐
│ Architectures │                         Closures │       min │    median │      mean │       max │   memory │ allocs │ samples │
├───────────────┼──────────────────────────────────┼───────────┼───────────┼───────────┼───────────┼──────────┼────────┼─────────┤
│           CPU │ AnisotropicBiharmonicDiffusivity │   3.634 s │   3.637 s │   3.637 s │   3.639 s │ 1.77 MiB │   2316 │       2 │
│           CPU │           AnisotropicDiffusivity │   2.045 s │   2.052 s │   2.059 s │   2.079 s │ 1.77 MiB │   2316 │       3 │
│           CPU │    AnisotropicMinimumDissipation │   3.240 s │   3.240 s │   3.240 s │   3.241 s │ 2.09 MiB │   2763 │       2 │
│           CPU │             IsotropicDiffusivity │   2.342 s │   2.344 s │   2.344 s │   2.345 s │ 1.77 MiB │   2316 │       3 │
│           CPU │                 SmagorinskyLilly │   3.501 s │   3.504 s │   3.504 s │   3.507 s │ 2.03 MiB │   2486 │       2 │
│           CPU │              TwoDimensionalLeith │   4.813 s │   4.820 s │   4.820 s │   4.828 s │ 1.88 MiB │   2481 │       2 │
│           GPU │ AnisotropicBiharmonicDiffusivity │ 24.699 ms │ 24.837 ms │ 26.946 ms │ 46.029 ms │ 3.16 MiB │  29911 │      10 │
│           GPU │           AnisotropicDiffusivity │ 16.115 ms │ 16.184 ms │ 16.454 ms │ 18.978 ms │ 2.97 MiB │  17169 │      10 │
│           GPU │    AnisotropicMinimumDissipation │ 15.858 ms │ 25.856 ms │ 24.874 ms │ 26.014 ms │ 3.57 MiB │  24574 │      10 │
│           GPU │             IsotropicDiffusivity │ 14.442 ms │ 17.415 ms │ 17.134 ms │ 17.513 ms │ 2.99 MiB │  19135 │      10 │
│           GPU │                 SmagorinskyLilly │ 16.315 ms │ 23.969 ms │ 23.213 ms │ 24.059 ms │ 3.86 MiB │  24514 │      10 │
│           GPU │              TwoDimensionalLeith │ 34.470 ms │ 34.628 ms │ 35.535 ms │ 43.798 ms │ 3.56 MiB │  45291 │      10 │
└───────────────┴──────────────────────────────────┴───────────┴───────────┴───────────┴───────────┴──────────┴────────┴─────────┘

              Turbulence closure CPU to GPU speedup
┌──────────────────────────────────┬─────────┬─────────┬─────────┐
│                         Closures │ speedup │  memory │  allocs │
├──────────────────────────────────┼─────────┼─────────┼─────────┤
│ AnisotropicBiharmonicDiffusivity │ 146.428 │ 1.78781 │ 12.9149 │
│           AnisotropicDiffusivity │ 126.804 │ 1.67787 │ 7.41321 │
│    AnisotropicMinimumDissipation │ 125.324 │ 1.70856 │ 8.89396 │
│             IsotropicDiffusivity │ 134.607 │ 1.69269 │ 8.26209 │
│                 SmagorinskyLilly │ 146.187 │ 1.89602 │ 9.86082 │
│              TwoDimensionalLeith │ 139.196 │ 1.89218 │ 18.2551 │
└──────────────────────────────────┴─────────┴─────────┴─────────┘

```



## Older Benchmarks
The following benchmark results are generated from an older version of Oceananigans and with deprecated benchmarking scripts. These legacy benchmark results can still be resonably used as a reference for gauging performance changes across versions of Oceananigans.

### Static ocean

This is a benchmark of a simple "static ocean" configuration. The time stepping and Poisson
solver still takes the same amount of time whether the ocean is static or active, so it is
indicative of actual performance. It tests the performance of a bare-bones
horizontally-periodic model with `topology = (Periodic, Periodic, Bounded)`.

```
Oceananigans v0.34.0 (DEVELOPMENT BRANCH)
Julia Version 1.4.2
Commit 44fa15b150* (2020-05-23 18:35 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) Silver 4214 CPU @ 2.20GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, skylake)
  GPU: TITAN V

 ──────────────────────────────────────────────────────────────────────────────────────
        Static ocean benchmarks                Time                   Allocations      
                                       ──────────────────────   ───────────────────────
           Tot / % measured:                 291s / 29.6%           27.7GiB / 0.50%    

 Section                       ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────────────────
  16× 16× 16  [CPU, Float32]       10   15.6ms  0.02%  1.56ms   2.61MiB  1.84%   267KiB
  16× 16× 16  [CPU, Float64]       10   16.9ms  0.02%  1.69ms   2.61MiB  1.84%   267KiB
  16× 16× 16  [GPU, Float32]       10   53.4ms  0.06%  5.34ms   11.5MiB  8.14%  1.15MiB
  16× 16× 16  [GPU, Float64]       10   69.7ms  0.08%  6.97ms   11.5MiB  8.14%  1.15MiB
  32× 32× 32  [CPU, Float32]       10   54.6ms  0.06%  5.46ms   2.61MiB  1.84%   267KiB
  32× 32× 32  [CPU, Float64]       10   57.1ms  0.07%  5.71ms   2.61MiB  1.84%   267KiB
  32× 32× 32  [GPU, Float32]       10   57.5ms  0.07%  5.75ms   11.6MiB  8.15%  1.16MiB
  32× 32× 32  [GPU, Float64]       10   75.0ms  0.09%  7.50ms   11.6MiB  8.16%  1.16MiB
  64× 64× 64  [CPU, Float32]       10    424ms  0.49%  42.4ms   2.61MiB  1.84%   267KiB
  64× 64× 64  [CPU, Float64]       10    425ms  0.49%  42.5ms   2.61MiB  1.84%   267KiB
  64× 64× 64  [GPU, Float32]       10   61.7ms  0.07%  6.17ms   11.6MiB  8.16%  1.16MiB
  64× 64× 64  [GPU, Float64]       10   82.4ms  0.10%  8.24ms   11.6MiB  8.17%  1.16MiB
 128×128×128  [CPU, Float32]       10    3.67s  4.26%   367ms   2.61MiB  1.84%   267KiB
 128×128×128  [CPU, Float64]       10    3.64s  4.23%   364ms   2.61MiB  1.84%   267KiB
 128×128×128  [GPU, Float32]       10   74.8ms  0.09%  7.48ms   11.6MiB  8.16%  1.16MiB
 128×128×128  [GPU, Float64]       10   94.0ms  0.11%  9.40ms   11.6MiB  8.17%  1.16MiB
 256×256×256  [CPU, Float32]       10    38.5s  44.8%   3.85s   2.61MiB  1.84%   267KiB
 256×256×256  [CPU, Float64]       10    37.9s  44.1%   3.79s   2.61MiB  1.84%   267KiB
 256×256×256  [GPU, Float32]       10    350ms  0.41%  35.0ms   11.6MiB  8.18%  1.16MiB
 256×256×256  [GPU, Float64]       10    352ms  0.41%  35.2ms   11.6MiB  8.17%  1.16MiB
 ──────────────────────────────────────────────────────────────────────────────────────

CPU Float64 -> Float32 speedup:
 16× 16× 16 : 1.084
 32× 32× 32 : 1.046
 64× 64× 64 : 1.000
128×128×128 : 0.993
256×256×256 : 0.986

GPU Float64 -> Float32 speedup:
 16× 16× 16 : 1.304
 32× 32× 32 : 1.303
 64× 64× 64 : 1.335
128×128×128 : 1.257
256×256×256 : 1.004

CPU -> GPU speedup:
 16× 16× 16  [Float32]: 0.291
 16× 16× 16  [Float64]: 0.242
 32× 32× 32  [Float32]: 0.949
 32× 32× 32  [Float64]: 0.762
 64× 64× 64  [Float32]: 6.876
 64× 64× 64  [Float64]: 5.152
128×128×128  [Float32]: 49.036
128×128×128  [Float64]: 38.730
256×256×256  [Float32]: 109.868
256×256×256  [Float64]: 107.863
```

### Channel

This benchmark tests the channel model (`topology = (Periodic, Bounded, Bounded)`)
configuration which can be slower due to the use of a more complicated algorithm
(involving 2D cosine transforms) for the pressure solver in the listed version
of Oceananigans.

```
Oceananigans v0.34.0 (DEVELOPMENT BRANCH)
Julia Version 1.4.2
Commit 44fa15b150* (2020-05-23 18:35 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) Silver 4214 CPU @ 2.20GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, skylake)
  GPU: TITAN V

 ──────────────────────────────────────────────────────────────────────────────────────
           Channel benchmarks                  Time                   Allocations      
                                       ──────────────────────   ───────────────────────
           Tot / % measured:                 453s / 19.5%           26.3GiB / 0.48%    

 Section                       ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────────────────
  32× 32× 32  [CPU, Float32]       10   58.5ms  0.07%  5.85ms   2.84MiB  2.22%   291KiB
  32× 32× 32  [CPU, Float64]       10   60.8ms  0.07%  6.08ms   2.85MiB  2.22%   291KiB
  32× 32× 32  [GPU, Float32]       10   68.7ms  0.08%  6.87ms   12.6MiB  9.85%  1.26MiB
  32× 32× 32  [GPU, Float64]       10   88.2ms  0.10%  8.82ms   12.6MiB  9.85%  1.26MiB
  64× 64× 64  [CPU, Float32]       10    459ms  0.52%  45.9ms   2.84MiB  2.22%   291KiB
  64× 64× 64  [CPU, Float64]       10    442ms  0.50%  44.2ms   2.85MiB  2.22%   291KiB
  64× 64× 64  [GPU, Float32]       10   91.0ms  0.10%  9.10ms   12.8MiB  10.0%  1.28MiB
  64× 64× 64  [GPU, Float64]       10    108ms  0.12%  10.8ms   12.8MiB  10.0%  1.28MiB
 128×128×128  [CPU, Float32]       10    3.87s  4.38%   387ms   2.84MiB  2.22%   291KiB
 128×128×128  [CPU, Float64]       10    3.92s  4.44%   392ms   2.85MiB  2.22%   291KiB
 128×128×128  [GPU, Float32]       10    145ms  0.16%  14.5ms   13.2MiB  10.3%  1.32MiB
 128×128×128  [GPU, Float64]       10    163ms  0.18%  16.3ms   13.2MiB  10.3%  1.32MiB
 256×256×256  [CPU, Float32]       10    38.6s  43.6%   3.86s   2.85MiB  2.22%   292KiB
 256×256×256  [CPU, Float64]       10    38.7s  43.8%   3.87s   2.85MiB  2.22%   292KiB
 256×256×256  [GPU, Float32]       10    805ms  0.91%  80.5ms   14.0MiB  10.9%  1.40MiB
 256×256×256  [GPU, Float64]       10    805ms  0.91%  80.5ms   14.0MiB  10.9%  1.40MiB
 ──────────────────────────────────────────────────────────────────────────────────────

CPU Float64 -> Float32 speedup:
 32× 32× 32 : 1.040
 64× 64× 64 : 0.963
128×128×128 : 1.015
256×256×256 : 1.004

GPU Float64 -> Float32 speedup:
 32× 32× 32 : 1.283
 64× 64× 64 : 1.188
128×128×128 : 1.120
256×256×256 : 0.999

CPU -> GPU speedup:
 32× 32× 32  [Float32]: 0.851
 32× 32× 32  [Float64]: 0.689
 64× 64× 64  [Float32]: 5.044
 64× 64× 64  [Float64]: 4.088
128×128×128  [Float32]: 26.602
128×128×128  [Float64]: 24.097
256×256×256  [Float32]: 47.891
256×256×256  [Float64]: 48.116
```
