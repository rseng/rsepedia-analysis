<div align="center">
  <h1>Kinetic.jl</h1>
  <img
    src="https://i.postimg.cc/ncXfgjXd/dancing-circles.gif"
    alt="Kinetic Logo" width="300">
  </img>

  [![version](https://juliahub.com/docs/Kinetic/version.svg)](https://juliahub.com/ui/Packages/Kinetic/wrVmu)
  [![](https://img.shields.io/badge/docs-latest-cornflowerblue)](https://xiaotianbai.com/Kinetic.jl/dev/)
  [![](https://img.shields.io/badge/docs-stable-blue)](https://xiaotianbai.com/Kinetic.jl/stable/)
  [![status](https://joss.theoj.org/papers/65d56efef938caf92c2cc942d2c25ea4/status.svg?style=flat-square)](https://joss.theoj.org/papers/65d56efef938caf92c2cc942d2c25ea4)
  [![Visits Badge](https://badges.pufler.dev/visits/vavrines/Kinetic.jl)](https://badges.pufler.dev)
</div>

<!--
![](https://img.shields.io/github/v/tag/vavrines/Kinetic.jl?include_prereleases&label=latest%20version&logo=github&sort=semver)
![](https://img.shields.io/badge/License-MIT-yellow.svg)
![](https://zenodo.org/badge/243490351.svg?style=flat-square)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![GitHub commits since tagged version](https://img.shields.io/github/commits-since/vavrines/Kinetic.jl/v0.7.0.svg?style=social&logo=github)](https://github.com/vavrines/Kinetic.jl)
-->

<!--<div align="center"> <img
  src="https://i.postimg.cc/ncXfgjXd/dancing-circles.gif"
  alt="Kinetic Logo" width="300"></img>
</div>-->
<!--
# Kinetic.jl
<img src="https://i.postimg.cc/ncXfgjXd/dancing-circles.gif" width="300"/>
-->

**Kinetic** is a lightweight [Julia](https://julialang.org) toolbox for the study of computational fluid dynamics.
The main module is split into portable components:

- [KitBase.jl](https://github.com/vavrines/KitBase.jl): basic physics and numerical schemes
- [KitML.jl](https://github.com/vavrines/KitML.jl): neural dynamics and machine learning methods
- [KitFort.jl](https://github.com/vavrines/KitFort.jl): alternative high-performance Fortran backend

As an optional module, the alternative Fortran backend can be manually imported into the current ecosystem when the ultimate executing efficiency is pursued.
A Python wrapper [kineticpy](https://github.com/vavrines/kineticpy) has been built as well to call the structs and methods here through [pyjulia](https://github.com/JuliaPy/pyjulia).

| [Kinetic](https://github.com/vavrines/Kinetic.jl) | [KitBase](https://github.com/vavrines/KitBase.jl) | [KitML](https://github.com/vavrines/KitML.jl) | [KitFort](https://github.com/vavrines/KitFort.jl) |
| ---------- | --------- | ---------------- | ------ |
| ![CI](https://img.shields.io/github/workflow/status/vavrines/Kinetic.jl/CI?style=flat-square) | ![CI](https://img.shields.io/github/workflow/status/vavrines/KitBase.jl/CI?style=flat-square) | ![CI](https://img.shields.io/github/workflow/status/vavrines/KitML.jl/CI?style=flat-square) | ![CI](https://img.shields.io/github/workflow/status/vavrines/KitFort.jl/CI?style=flat-square) |
| [![codecov](https://img.shields.io/codecov/c/github/vavrines/Kinetic.jl?style=flat-square)](https://codecov.io/gh/vavrines/Kinetic.jl) | [![codecov](https://img.shields.io/codecov/c/github/vavrines/KitBase.jl?style=flat-square)](https://codecov.io/gh/vavrines/KitBase.jl) | [![codecov](https://img.shields.io/codecov/c/github/vavrines/KitML.jl?style=flat-square)](https://codecov.io/gh/vavrines/KitML.jl) | [![codecov](https://img.shields.io/codecov/c/github/vavrines/KitFort.jl?style=flat-square)](https://codecov.io/gh/vavrines/KitFort.jl) |

## Installation

Kinetic.jl is a registered package in the official [Julia package registry](https://github.com/JuliaRegistries/General).
We recommend installing it with the Julia package manager. 
From the Julia REPL, you can get in the package manager (by pressing `]`) and add the package

```julia
julia> ]
(v1.7) pkg> add Kinetic
```
This will automatically install a currently stable release and all its dependencies.
Similarly, the previously installed versions can be updated to the latest tagged release by

```julia
(v1.7) pkg> update Kinetic
```

## Physics

Kinetic.jl focuses on theoretical and numerical studies of many-particle systems of gases, photons, plasmas, neutrons, etc.
It employs the finite volume method (FVM) to conduct 1-3 dimensional numerical simulations on CPUs and GPUs.
Any advection-diffusion-type equation can be solved within the framework.
Special attentions have been paid on Hilbert's sixth problem, i.e. to build the numerical passage between kinetic theory of gases, e.g. the Boltzmann equation, and continuum mechanics, e.g. the Euler and Navier-Stokes equations.
A partial list of current supported models and equations include:
- Boltzmann equation
- radiative transfer equation
- Fokker-Planck-Landau equation
- direct simulation Monte Carlo
- advection-diffusion equation
- Burgers equation
- Euler equations
- Navier-Stokes equations
- Magnetohydrodynamical equations
- Maxwell's equations

## Documentation

For the detailed information on the implementation and usage of the package,
[check the documentation](https://xiaotianbai.com/Kinetic.jl/dev/).

## Citing

If you benefit from Kinetic.jl in your research, teaching, or other activities, we would be happy if you could mention or cite it:

```
@article{Xiao2021,
  doi = {10.21105/joss.03060},
  url = {https://doi.org/10.21105/joss.03060},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {62},
  pages = {3060},
  author = {Tianbai Xiao},
  title = {Kinetic.jl: A portable finite volume toolbox for scientific and neural computing},
  journal = {Journal of Open Source Software}
}
```

## Contributing

If you have further questions regarding Kinetic.jl or have got an idea on improving it, please feel free to get in touch. Open an issue or pull request if you'd like to work on a new feature or even if you're new to open-source and want to find a cool little project or issue to work on that fits your interests. We're more than happy to help along the way.
# Known Issues

- FFT doesn't work well with offset arrays---
title: 'Kinetic.jl: A portable finite volume toolbox for scientific and neural computing'
tags:
  - kinetic theory
  - computational fluid dynamics
  - scientific machine learning
  - julia
authors:
  - name: Tianbai Xiao
    orcid: 0000-0001-9127-9497
    affiliation: 1
affiliations:
 - name: Karlsruhe Institute of Technology, 76131 Karlsruhe, Germany
   index: 1
date: 06 January 2021
bibliography: paper.bib
---

# Summary

Kinetic.jl is a lightweight finite volume toolbox written in the Julia programming language for the study of computational physics and scientific machine learning.
It is an open-source project hosted on GitHub and distributed under the MIT license.
The main module consists of KitBase.jl for basic physics and KitML.jl for neural dynamics.
The function library provides a rich set of numerical fluxes and source terms for differential and integral equations.
Any advection-diffusion type mechanical or neural equation can be set up and solved within the framework.
Machine learning methods can be seamlessly integrated to build data-driven closure models and accelerate the calculation of nonlinear terms.
The package is designed to balance programming flexibility for scientific research, algorithmic efficiency for applications, the simplicity for educational usage.

# Statement of need

A physical system can perform a wonderfully diverse set of acts on different characteristic scales.
It is challenging to propose a universal theory that can be applied to describing multi-scale physical evolutions quantitatively.
For example, particle transport can be depicted statistically by fluid mechanics at a macroscopic level [@batchelor2000], but needs to be followed in more detail by the kinetic theory of gases at the molecular mean free path scale [@chapman1990].
With rapidly advancing computing power, the finite volume method (FVM) provides a prevalent method to conduct direct numerical simulations based on first physical principles.

Most existing FVM libraries, e.g., OpenFOAM [@jasak2007], are dedicated to solving the Euler and the Navier-Stokes equations.
Very limited work has been done for phase-field models [@zhu2017; @krause2021].
Since classical fluid dynamics basically requires an one-shot simulation process from initial to final solution fields, these libraries are mostly written in compiled languages (C/C++ and Fortran).
Such approaches enjoy good execution efficiency but sacrifice the flexibility of secondary development.
This makes it cumbersome to integrate existing numerical solvers with scientific machine learning (SciML) packages, as interactive programming is becoming a mainstream practice in data science.
This also causes considerable difficulties to general or educational users who are not familiar with the package in configuring environments and compiling binaries.

One compromise can be made by using a combination of static and dynamic languages [@clawpack2020], where the high-level front-ends and the low-level computational back-ends are split.
This methodology benefits general users, while researchers still need to work on the back-end if a new feature is required. 
The so-called two-language problem introduces additional tradeoffs in both development and execution.
For example, a two-tiered system brings unavoidable challenges for type domain transition and memory management.
Special attention needs to be paid on optimizing the high-level codes, e.g., the vectorization of massive computation part, which can be unnatural in a physical simulation and might generate additional temporary objects. 
In addition, interfacing between layers may add significant overhead and makes whole-program optimization much more difficult [@bezanson2012].
Unlike these packages, Kinetic.jl is built upon the Julia programming language [@bezanson2017], which is dynamically typed and designed for high performance computing for a broad range of devices. 
Based on type inference and multiple dispatch, it is a promising choice to solve the two-language problem.

Kinetic.jl focuses on the theoretical and numerical studies of many-particle systems of gases, photons, plasmas, neutrons, etc. [@xiao2017; @xiao2020a]
A hierarchy of abstractions is implemented in the library.
At the highest level, it is feasible to model and simulate a fluid dynamic problem within ten lines of code. 
At the lowest level, we designed methods for general numbers and arrays so that it is possible to cooperate with existing packages in Julia ecosystem.
For example, Flux.jl [@Flux2018] can be used to create and train scientific machine learning models.
Innovations of the package are:

- 100% Julia stack that encounters no two-language problem

- Comprehensive support for kinetic theory and phase-space equations

- Lightweight design to ensure the flexibility for secondary development

- Close coupling with scientific machine learning

# KitBase.jl

The main module of Kinetic.jl is split into two pieces to reduce the just-in-time (JIT) compilation time for domain specific applications.
The basic physical laws and finite volume method are implemented in KitBase.jl.
It provides a variety of solvers for the Boltzmann equation, Maxwell's equations, advection-diffusion equation, Burgers' equation, Euler and Navier-Stokes equations, etc.
Different parallel computing techniques are provided, e.g., multi-threading, distributed computing, and CUDA programming.

In the following, we present an illustrative example of solving a lid-driven cavity problem with the Boltzmann equation. 
Two initialization methods, i.e., configuration text and Julia script, are available for setting up the solver.
With the configuration file `config.toml` set as below,
```toml
# setup
matter = gas # material
case = cavity # case
space = 2d2f2v # phase
flux = kfvs # flux function
collision = bgk # intermolecular collision
nSpecies = 1 # number of species
interpOrder = 2 # interpolation order of accuracy
limiter = vanleer # limiter function
boundary = maxwell # boundary condition
cfl = 0.8 # CFL number
maxTime = 5.0 # maximal simulation time

# physical space
x0 = 0.0 # starting point in x
x1 = 1.0 # ending point in x
nx = 45 # number of cells in x
y0 = 0.0 # starting point in y
y1 = 1.0 # ending point in y
ny = 45 # number of cells in y
pMeshType = uniform # mesh type
nxg = 0 # number of ghost cell in x
nyg = 0 # number of ghost cell in y

# velocity space
umin = -5.0 # starting point in u
umax = 5.0 # ending point in u
nu = 28 # number of cells in u
vmin = -5.0 # starting point in v
vmax = 5.0 # ending point in v
nv = 28 # number of cells in v
vMeshType = rectangle # mesh type
nug = 0 # number of ghost cell in u
nvg = 0 # number of ghost cell in v

# gas property
knudsen = 0.075 # Knudsen number
mach = 0.0 # Mach number
prandtl = 1.0 # Prandtl number
inK = 1.0 # molecular inner degree of freedom
omega = 0.72 # viscosity index of hard-sphere gas
alphaRef = 1.0 # viscosity index of hard-sphere gas in reference state
omegaRef = 0.5 # viscosity index of hard-sphere gas ub reference state

# boundary condition
uLid = 0.15 # U-velocity of moving wall
vLid = 0.0 # V-velocity of moving wall
tLid = 1.0 # temperature of wall
```

we can execute the following codes
```julia
using Kinetic
set, ctr, xface, yface, t = initialize("config.toml")
t = solve!(set, ctr, xface, yface, t)
plot_contour(set, ctr)
```

In the above codes, the computational setup is stored in `set`. 
The solutions over control volumes are represented in an array `ctr`, while `xface` and `yface` record the interface fluxes along x and y directions.
In this example, the structured mesh is generated automatically by Kinetic.jl, while a non-structured mesh file can also be imported and used for computation.
The result is visualized with built-in function `plot_contour`, which presents the distributions of gas density, velocity, and temperature inside the cavity.

![Fig. 1](cavity.png)
Fig. 1: macroscopic variables in the lid-driven cavity (top left: density, top right: U-velocity, bottom left: V-velocity, bottom right: temperature).

# KitML.jl

Machine learning has increasing momentum in scientific computing.
Given the nonlinear structure of differential and integral equations, it is promising to incorporate the universal function approximators from machine learning surrogate models into the governing equations and achieve a better balance between efficiency and accuracy.
In KitML.jl, we implement strategies to construct hybrid mechanical-neural differential operators and form structure-preserving data-driven closure models.
The detailed background can be found in @xiao2020b.

# Extension

Numerical simulations of nonlinear models and differential equations are essentially connected with supercomputers and high-performance computing (HPC). 
Considering that some existing hardware architecture, e.g., Sunway TaihuLight with Chinese-designed SW26010 processors, only provides optimization for specific languages, we have developed an accompanying package KitFort.jl.
This is not a default component of Kinetic.jl but can be manually imported.
In addition, a wrapper, kineticpy, has been built to locate structures and methods from the Python ecosystem. 

# Acknowledgements

The current work is funded by the Alexander von Humboldt Foundation (Ref3.5-CHN-1210132-HFST-P).

# References
# Reconstruction

```@docs
reconstruct!
```

The reconstruction solver interpolates piecewise solutions with the desirable order of accuracy.
The reconstruction stencils can be based on 2 or 3 cells
```@docs
reconstruct2
reconstruct2!
reconstruct3
reconstruct3!
```

The available schemes are
```@docs
vanleer
minmod
superbee
vanalbaba
weno5
```
# Guide for Contributors

Thank you for considering contributing to Kinetic! This short guide will
give you ideas on how you can contribute and help you make a contribution.
Please feel free to ask us questions and chat with us at any time if you're
unsure about anything.

## What can I do?

* Tackle an existing issue.
* Try to run Kinetic and play around with it to simulate your favorite
  fluid and kinetic physics. If you run into any problems or find it difficult
  to use or understand, please open an issue!
* Write up an example or tutorial on how to do something useful with
  Kinetic, like how to set up a new physical configuration.
* Improve documentation or comments if you found something hard to use.
* Implement a new feature if you need it to use Kinetic.

If you're interested in working on something, let us know by commenting on
existing issues or by opening a new issue if. This is to make sure no one else
is working on the same issue and so we can help and guide you in case there
is anything you need to know beforehand.

## Philosophy

* Each pull request should consist of a logical collection of changes. You can
  include multiple bug fixes in a single pull request, but they should be related.
  For unrelated changes, please submit multiple pull requests.
* Do not commit changes to files that are irrelevant to your feature or bugfix
  (eg: .gitignore).
* Be willing to accept criticism and work on improving your code; we don't want
  to break other users' code, so care must be taken not to introduce bugs. We
  discuss pull requests and keep working on them until we believe we've done a
  good job.
* Be aware that the pull request review process is not immediate, and is
  generally proportional to the size of the pull request.

## Reporting a bug

The easiest way to get involved is to report issues you encounter when using
Kinetic or by requesting something you think is missing.

* Head over to the issues in [KitBase](https://github.com/vavrines/KitBase.jl/issues) or [KitML](https://github.com/vavrines/KitML.jl/issues) page.
* Search to see if your issue already exists or has even been solved previously.
* If you indeed have a new issue or request, click the "New Issue" button.
* Please be as specific as possible. Include the version of the code you were using, as
  well as what operating system you are running. The output of Julia's `versioninfo()`
  and `] status` is helpful to include. If possible, include complete, minimal example
  code that reproduces the problem.

## Setting up your development environment

* Install [Julia](https://julialang.org/) on your system.
* Install git on your system if it is not already there (install XCode command line tools on
  a Mac or git bash on Windows).
* Login to your GitHub account and make a fork of the
  [KitBase](https://github.com/vavrines/KitBase.jl) or [KitML](https://github.com/vavrines/KitML.jl) by
  clicking the "Fork" button.
* Clone your fork of the Kinetic repository (in terminal on Mac/Linux or git shell/
  GUI on Windows) in the location you'd like to keep it.
  ```
  git clone https://github.com/your-user-name/KitBase.jl.git or https://github.com/your-user-name/KitML.jl.git
  ```
* Navigate to that folder in the terminal or in Anaconda Prompt if you're on Windows.
* Connect your repository to the upstream (main project).
  ```git remote add KitBase https://github.com/vavrines/KitBase.jl.git``` or 
  ```git remote add KitML https://github.com/vavrines/KitML.jl.git```
* Create the development environment by opening Julia via `julia --project` then
  typing in `] instantiate`. This will install all the dependencies in the Project.toml
  file.
* You can test to make sure Kinetic works by typing in `] test` which will run all
  the tests (this can take a while).

Your development environment is now ready!

## Pull requests

Changes and contributions should be made via GitHub pull requests against the ``master`` branch.

When you're done making changes, commit the changes you made. Chris Beams has
written a [guide](https://chris.beams.io/posts/git-commit/) on how to write
good commit messages.

When you think your changes are ready to be merged into the main repository,
push to your fork and submit a pull request in https://github.com/vavrines/KitBase.jl/compare/ and https://github.com/vavrines/KitML.jl/compare/.

**Working on your first Pull Request?** You can learn how from the video series
[How to Contribute to an Open Source Project on GitHub](https://egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github), Aaron Meurer's [tutorial on the git workflow](https://www.asmeurer.com/git-workflow/), or the guide [“How to Contribute to Open Source"](https://opensource.guide/how-to-contribute/).

## Documentation

Now that you've made your awesome contribution, it's time to tell the world how to use it.
Writing documentation strings is really important to make sure others use your functionality
properly. Didn't write new functions? That's fine, but be sure that the documentation for
the code you touched is still in great shape. It is not uncommon to find some strange wording
or clarification that you can take care of while you are here.

## Credits

This contributor's guide is based on the [MetPy contributor's guide](https://github.com/Unidata/MetPy/blob/master/CONTRIBUTING.md) and [Oceananigans](https://github.com/CliMA/Oceananigans.jl).# Illustrative examples

Thanks to the brilliant expressiveness and low-overhead abstraction in Julia, we provide different levels of solution algorithm for modeling and simulating advection-diffusion dynamics.
The high-level solver is able to solve complex physics in a few lines, while the low-level APIs keep all the detailed implementations and benefit the secondary development.
The low-level methods are easy to be called from Python and C.
In the following, we present some quick tutorials to illustrate the usage of Kinetic.
For more examples, please refer the example directories in [Kinetic.jl](https://github.com/vavrines/Kinetic.jl/tree/master/example), [KitBase.jl](https://github.com/vavrines/KitBase.jl/tree/main/example) and [KitML.jl](https://github.com/vavrines/KitML.jl/tree/main/example).
# KitFort and high performance computing

Numerical simulations of nonlinear models and differential equations are essentially connected with supercomputers and high-performance computing (HPC).
The performance of a supercomputer or a software program is commonly measured in floating-point operations per second (FLOPS).
Through the milestone astronomy research of [Celeste](https://juliacomputing.com/case-studies/celeste/), Julia has entered the PetaFLOPS club (together with C/C++ and Fortran) since 2017.
Julia is experiencing a dramatic Rise in HPC and elsewhere, and that is why we use Julia to organize the Kinetic.
However, compared with the mature C/C++ ecosystem, the equivalent execution efficiency isn't going to happen in all time and situations.
Some existing hardware architecture, e.g. [Sunway TaihuLight](https://en.wikipedia.org/wiki/Sunway_TaihuLight), the previou fastest supercomputer in [TOP500](https://www.top500.org/) list, is built upon 40,960 Chinese-designed SW26010 manycore 64-bit RISC processors, which is not specifically optimized for Julia.
Therefore, we've develop an accompanying package [KitFort.jl](https://github.com/vavrines/KitFort.jl).
The Fortran codes have been linked to the Julia syntax with the built-in `ccall` function.
It's not a default submodule of Kinetic since we believe the Julia codes are sufficient for general users and developers and encounter no two-language problem.
However, it can be manually imported when the executing efficiency becomes the first priority by executing
```julia
julia> ]
(v1.5) pkg> add KitFort
```
After that, using/import the package.
```julia
julia> using KitFort
```

It can be updated to the latest tagged release from the package manager by executing
```julia
(v1.5) pkg> update KitFort
```
# Math

```@docs
linspace
heaviside
fortsign
mat_split
central_diff
central_diff!
upwind_diff
upwind_diff!
unstruct_diff
KitBase.lgwt
KitBase.extract_last
```# Kinetic.jl

Kinetic is a portable [Julia](https://julialang.org/) toolbox for the study of computational fluid dynamics and scientific machine learning.
The default module consists of [KitBase.jl](https://github.com/vavrines/KitBase.jl) with basic physics and [KitML.jl](https://github.com/vavrines/KitML.jl) with neural dynamics. 
The high-performance Fortran library [KitFort.jl](https://github.com/vavrines/KitFort.jl) can be manually imported when the ultimate efficiency is pursued.
Besides, a wrapper [kineticpy](https://github.com/vavrines/kineticpy) is built to locate the data hierarchies and methods in Python.

## Scope of application

Kinetic is interested in the evolution of many-particle systems, e.g. gases, photons, plasmas, neutrons, electrons, etc.
Based on the finite volume method (FVM), it provides an efficient tool where 1-3 dimensional theoretical modeling and numerical simulation can be conducted.
Any advection-diffusion type equation can be hooked within the framework.
Special attention has been paid to the [kinetic theory](https://en.wikipedia.org/wiki/Kinetic_theory_of_gases) and the Boltzmann-type equations,
which depicts the time-space evolution of particles via ensemble averaging at the mesoscopic level.
A partial list of current supported models and equations is as follows.
- linear Boltzmann equation
- nonlinear Boltzmann equation
- multi-component Boltzmann equation
- Fokker-Planck-Landau equation
- direct simulation Monte Carlo
- advection-diffusion equation
- Burgers' equation
- Euler equations
- Navier-Stokes equations
- Extended hydrodynamical equations from asymptotic expansion
- Magnetohydrodynamical equations
- Maxwell's equations

## Design philosophy

The code hierarchy is designed as intuitive and neat as possible.
It's dedicated to providing a friendly interface for educational usage in kinetic theory and rich functionality for scientific research.
Benefiting from the brilliant expressiveness and low-overhead abstraction provided by the [Julia programming language](https://julialang.org/), 
we provide different levels of APIs to allow the users to focus on physics and to cooperate with the existing packages in the Julia ecosystem.

## What is new?

Finite volume method is a proven approach for simulating conservation laws.
Compared with the existing open-source softwares, e.g. [OpenFOAM](https://openfoam.org/), [SU2](https://su2code.github.io/) and [Clawpack](https://www.clawpack.org/), 
Kinetic holds the novelty through the following points:
- 100% Julia stack that encounters no two-language problem
- Comprehensive support for kinetic theory and phase-space equations
- Lightweight design to ensure the flexibility for secondary development
- Closely coupling with scientific machine learning

## How to get help?

If you are interested in using Kinetic.jl or are trying to figure out how to use it, please feel free to get in touch and raise questions.
Do open an issue or pull request if you have questions, suggestions or solutions.# Multiple threading

The multi-threading computation is built upon Julia's `@threads` macro.
```julia
Base.Threads.@threads for ... end
```
It provides an OpenMP type parallelization.
The iteration space is splitted among multiple tasks and those tasks are parallelized on threads according to a scheduling policy.
A barrier is placed at the end of the loop which waits for all tasks to finish execution.

In Kinetic, `@threads` is set in front of the loops for reconstruction, evolution and update.
For example, the evaluation of fluxes is conducted as follows.
```julia
@inbounds Threads.@threads for i = idx0:idx1
    flux_gks!(
        face[i].fw,
        ctr[i-1].w .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sw,
        ctr[i].w .- 0.5 .* ctr[i].dx .* ctr[i].sw,
        KS.gas.γ,
        KS.gas.K,
        KS.gas.μᵣ,
        KS.gas.ω,
        dt,
        0.5 * ctr[i-1].dx,
        0.5 * ctr[i].dx,
        ctr[i-1].sw,
        ctr[i].sw,
    )
end
```
It automatically makes use of multiple threading if Julia is initialized with
```bash
julia -t n
```

Besides of `@threads`, finer dispatch can be made with `@spawn` macro.
```julia
Base.Threads.@spawn
```
It creates and runs a task on any available thread.
To wait for the task to finish, call `wait` on the result of this macro, or call fetch to wait and then obtain its return value.
Values can be interpolated into `@spawn` via `$`, which copies the value directly into the constructed underlying closure.
It allows user to insert the value of a variable, isolating the aysnchronous code from changes to the variable's value in the current task.
This can be conduct with the low-level reconstruction, flux and step functions.# Installation Instructions

Kinetic is a registered Julia package in the official entry.
We recommend installing it with the built-in Julia package manager.
It automatically installs a currently stable and tagged release. 
From the Julia REPL, you can add the package.
```julia
julia> ]
(v1.6) pkg> add Kinetic
```

This will automatically install Kinetic and all its dependencies, and it's not needed to build the package manually.
You can also build the dependencies if some of them were removed by mistake.
```julia
julia> ]
(v1.6) pkg> build Kinetic
```
After that, we can `using` or `import` the package.
`using` will load the module and make its exported names available for direct use.
```julia
julia> using Kinetic
julia> linspace(0, 1, 5)
5-element Vector{Float64}:
 0.0
 0.25
 0.5
 0.75
 1.0
```
Correspondingly, `import` only loads the module while the names needs to be accessed with dot syntax.
```julia
julia> import Kinetic
julia> Kinetic.linspace(0, 1, 5)
5-element Vector{Float64}:
 0.0
 0.25
 0.5
 0.75
 1.0
```

Kinetic can be updated to the latest tagged release from the package manager.
```julia
(v1.5) pkg> update Kinetic
```

!!! tip "Use Julia 1.3 or newer"
    Kinetic matches perfectly with Julia 1.3 and newer versions.
    Installing it with an older version of Julia will locate incomplete functionality.# Index of Types and Methods

```@index
```# Configuration of solver

Kinetic is organized with the data structures and methods of both generality and convenience. 
While the low-level methods can be applied to multi-dimensional arrays directly, we provide a set of domain-specific structs that handles multiple dispatch in an elegant way.

For a solver pending for execution, its configurations can be handled in a `SolverSet <: AbstractSolverSet` struct.
```@docs
SolverSet
```
It contains six fields:
- set: general setup of a simulation
- pSpace: physical space settings
- vSpace: particle velocity space settings
- gas: properties of the simulated substance
- ib: initial and boundary conditions
- outputFolder: file directory for the output results

This struct plays an key role in the solution algorithm.# Benchmark

Here we provide a benchmark to identity the performance variation between Julia and Fortran implementations.
For brevity, we direct make use of the dynamic library [kitmod.so](https://github.com/vavrines/KitFort.jl/tree/main/src/fortran) by `ccall` function in Julia, and
compare the efficiency of computing numerical fluxes by [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl).

```julia
using Kinetic, BenchmarkTools

begin
    u = collect(-5.0:0.1:5.0)
    nu = length(u)
    weights = ones(nu) .* 0.5

    fw = zeros(3)
    fh = zeros(nu)
    fb = zeros(nu)

    inK = 2
    γ = 5.0 / 3.0
    primL = [1., 0., 1.]
    wL = prim_conserve(primL, γ)
    hL = maxwellian(u, primL) |> Array;
    bL = hL .* 2 ./ (2.)
    shL = zeros(nu)
    sbL = zeros(nu)
    lenL = 0.1

    primR = [0.5, 0., 1.]
    wR = prim_conserve(primR, γ)
    hR = maxwellian(u, primR) |> Array;
    bR = hR .* 2 ./ (2.)
    shR = zeros(nu)
    sbR = zeros(nu)
    lenR = 0.1

    muref = 0.001
    omega = 0.72
    prandtl = 1.0
    dt = 1e-4
end

#--- kfvs ---#
@btime ccall(
    (:__kinetic_MOD_flux_kfvs_2f1v, "kitmod.so"),
    Nothing,
    (
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64},
        Ref{Float64}, 
        Ref{Int}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64},
        Ref{Float64},
    ),
    fw,
    fh,
    fb,
    hL,
    bL,
    hR,
    bR,
    u,
    weights,
    nu,
    dt,
    shL,
    sbL,
    shR,
    sbR,
)

@btime flux_kfvs!(fw, fh, fb, hL, bL, hR, bR, u, weights, dt, shL, sbL, shR, sbR)

#--- ugks ---#
@btime ccall(
    (:__kinetic_MOD_flux_ugks_2f1v, "kitmod.so"),
    Nothing,
    (
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64},
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Int}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64},
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64},
    ),
    fw,
    fh,
    fb,
    wL,
    hL,
    bL,
    wR,
    hR,
    bR,
    u,
    weights,
    nu,
    inK,
    γ,
    muref,
    omega,
    prandtl,
    dt,
    lenL,
    lenR,
    shL,
    sbL,
    shR,
    sbR,
)

@btime flux_ugks!(fw, fh, fb, wL, hL, bL, wR, hR, bR, u, weights, inK, γ, muref, omega, prandtl, dt, lenL, lenR, shL, sbL, shR, sbR)

```

The results on a intel NUC8i7BEH with i7-8559U with 101 velocity points is as follows

Kinetic.jl
- KFVS flux ~ 6.747 μs (13 allocations: 11.38 KiB)
- UGKS flux ~ 13.344 μs (123 allocations: 20.94 KiB)

KitFort.jl
- KFVS flux ~ 5.421 μs (37 allocations: 800 bytes)
- UGKS flux ~ 11.413 μs (55 allocations: 1.09 KiB)

As presented, there is an improvement on efficiency by around 15%.# Advection diffusion

The first example is the scalar advection-diffusion equation.
It's a one dimensional problem in spatial domain ``x``.
Let's first configure the solver setup.
```julia
using Kinetic, Plots

set = Setup(
    matter = "scalar", # material
    case = "advection", # test case
    space = "1d0f0v", # phase space
    flux = "gks", # flux
    collision = "", # collision: for scalar conservation laws there are none
    interpOrder = 1, # interpolation order
    boundary = "period", # boundary condition
    cfl = 0.5, # cfl
    maxTime = 1.0, # simulation time
)
```

Then we generate the computational mesh.
Since we solve the macroscopic transport equation, the phase space is set to be nothing.
```julia
ps = PSpace1D(0.0, 1.0, 100, 1)
vs = nothing
```

And we define the physical property of material.
For the advection-diffusion equation, the two fields are the advection speed and viscosity respectively.
```julia
property = Scalar(1.0, 1e-6)
```

A sine wave is used as the initial condition.
```julia
ib = IB(x -> sin(2π * x), property)
```

For brevity, the above setups can be integrated into a single structure.
We also allocate the structures for cell-centered solutions and interface fluxes.
```julia
ks = SolverSet(set, ps, vs, property, ib)
ctr, face = init_fvm(ks)
```

The solution algorithm can be processed together with visualization.
```julia
t = 0.0
dt = KitBase.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int

anim = @animate for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, 0.0)

    plot(ks, ctr, xlabel="x", label="u", ylims=[-1, 1])
end

gif(anim, "advection.gif", fps = 45)
```

![](./assets/advection.gif)
# Physical space

A struct `set <: AbstractPhysicalSpace` defines the geometric setup of a simulation.
For the structured topology, structs for 1 and 2 dimensional physical space are built.
```@docs
PSpace1D
PSpace2D
```
It contains:
- x0 (y0): location of starting point
- x1 (y1): location of ending point
- nx (ny): number of cells in one direction
- x (y): locations of middle points of all cells
- dx (dy): intervals of all cell points

Besides, a unstrctured mesh struct is built, which supports 1-3 dimensional geometries.
```@docs
UnstructPSpace
```
It can be created by the built-in mesh reader.
```@docs
read_mesh
```
# Finite volume data

In the finite volume method, the data is stored separately throughout the cells.
Therefore, we provide `AbstractControlVolume` and `AbstractInterface` structs for processing in-cell and edge information,
which are used as arrays of structs (AoS) in numerical simulations.
Considering one-dimensional physical space ``x``, we provide the following control volume structs.
The structs differs from the number of particle distribution functions.
```@docs
ControlVolume1D
ControlVolume1D1F
ControlVolume1D2F
ControlVolume1D3F
ControlVolume1D4F
```

Within each cell, different numbers of particle distribution function can be tracked.
The interface data is stored correspondingly.
```@docs
Interface1D
Interface1D1F
Interface1D2F
Interface1D3F
Interface1D4F
```

The 2D control volume structs are implemented as well.
```@docs
ControlVolume2D
ControlVolume2D1F
ControlVolume2D2F
ControlVolume2D3F
```

The numerical fluxes are evaluated through `AbstractInterface` structs.
```@docs
Interface2D
Interface2D1F
Interface2D2F
```
# Lid-driven cavity

We then show the lid-driven cavity.
It's a four dimensional problem, with two in physical domain ``(x,y)`` and another in particle velocity domain ``(u,v)``.
Similarly, we prepare the configuration file as
```
# setup
matter = gas
case = cavity
space = 2d2f2v
flux = kfvs
collision = bgk
nSpecies = 1
interpOrder = 2
limiter = vanleer
boundary = maxwell
cfl = 0.8
maxTime = 5.0

# phase space
x0 = 0.0
x1 = 1.0
nx = 45
y0 = 0.0
y1 = 1.0
ny = 45
pMeshType = uniform
nxg = 0
nyg = 0

# velocity space
umin = -5.0
umax = 5.0
nu = 28
vmin = -5.0
vmax = 5.0
nv = 28
vMeshType = rectangle
nug = 0
nvg = 0

# gas
knudsen = 0.075
mach = 0.0
prandtl = 1.0
inK = 1.0
omega = 0.72
alphaRef = 1.0
omegaRef = 0.5

# boundary
uLid = 0.15
vLid = 0.0
tLid = 1.0
```

We then execute the following codes to conduct a simulation
```julia
using Kinetic
ks, ctr, a1face, a2face, t = initialize("config.txt")
t = solve!(ks, ctr, a1face, a2face, t)
```

The high-level solver `solve!` is equivalent as the following low-level procedures
```julia
using ProgressMeter
res = zeros(4)
dt = timestep(ks, ctr, t)
nt = floor(ks.set.maxTime / dt) |> Int
@showprogress for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, a1face, a2face, dt; mode = Symbol(ks.set.flux), bc = Symbol(ks.set.boundary))
    update!(ks, ctr, a1face, a2face, dt, res; coll = Symbol(ks.set.collision), bc = Symbol(ks.set.boundary))
end
```

It can be further expanded into the lower-level backend.
```julia
# lower-level backend 
@showprogress for iter = 1:nt
    # horizontal flux
    @inbounds Threads.@threads for j = 1:ks.pSpace.ny
        for i = 2:ks.pSpace.nx
            KitBase.flux_kfvs!(
                a1face[i, j].fw,
                a1face[i, j].fh,
                a1face[i, j].fb,
                ctr[i-1, j].h,
                ctr[i-1, j].b,
                ctr[i, j].h,
                ctr[i, j].b,
                ks.vSpace.u,
                ks.vSpace.v,
                ks.vSpace.weights,
                dt,
                a1face[i, j].len,
            )
        end
    end
    
    # vertical flux
    vn = ks.vSpace.v
    vt = -ks.vSpace.u
    @inbounds Threads.@threads for j = 2:ks.pSpace.ny
        for i = 1:ks.pSpace.nx
            KitBase.flux_kfvs!(
                a2face[i, j].fw,
                a2face[i, j].fh,
                a2face[i, j].fb,
                ctr[i, j-1].h,
                ctr[i, j-1].b,
                ctr[i, j].h,
                ctr[i, j].b,
                vn,
                vt,
                ks.vSpace.weights,
                dt,
                a2face[i, j].len,
            )
            a2face[i, j].fw .= KitBase.global_frame(a2face[i, j].fw, 0., 1.)
        end
    end
    
    # boundary flux
    @inbounds Threads.@threads for j = 1:ks.pSpace.ny
        KitBase.flux_boundary_maxwell!(
            a1face[1, j].fw,
            a1face[1, j].fh,
            a1face[1, j].fb,
            ks.ib.bcL,
            ctr[1, j].h,
            ctr[1, j].b,
            ks.vSpace.u,
            ks.vSpace.v,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[1, j].dy,
            1.,
        )

        KitBase.flux_boundary_maxwell!(
            a1face[ks.pSpace.nx+1, j].fw,
            a1face[ks.pSpace.nx+1, j].fh,
            a1face[ks.pSpace.nx+1, j].fb,
            ks.ib.bcR,
            ctr[ks.pSpace.nx, j].h,
            ctr[ks.pSpace.nx, j].b,
            ks.vSpace.u,
            ks.vSpace.v,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[ks.pSpace.nx, j].dy,
            -1.,
        )
    end
    
    @inbounds Threads.@threads for i = 1:ks.pSpace.nx
        KitBase.flux_boundary_maxwell!(
            a2face[i, 1].fw,
            a2face[i, 1].fh,
            a2face[i, 1].fb,
            ks.ib.bcD,
            ctr[i, 1].h,
            ctr[i, 1].b,
            vn,
            vt,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[i, 1].dx,
            1,
        )
        a2face[i, 1].fw .= KitBase.global_frame(a2face[i, 1].fw, 0., 1.)
        
        KitBase.flux_boundary_maxwell!(
            a2face[i, ks.pSpace.ny+1].fw,
            a2face[i, ks.pSpace.ny+1].fh,
            a2face[i, ks.pSpace.ny+1].fb,
            [1., 0.0, -0.15, 1.0],
            ctr[i, ks.pSpace.ny].h,
            ctr[i, ks.pSpace.ny].b,
            vn,
            vt,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[i, ks.pSpace.ny].dy,
            -1,
        )
        a2face[i, ks.pSpace.ny+1].fw .= KitBase.global_frame(
            a2face[i, ks.pSpace.ny+1].fw,
            0.,
            1.,
        )
    end

    # update
    @inbounds for j = 1:ks.pSpace.ny
        for i = 1:ks.pSpace.nx
            KitBase.step!(
                ctr[i, j].w,
                ctr[i, j].prim,
                ctr[i, j].h,
                ctr[i, j].b,
                a1face[i, j].fw,
                a1face[i, j].fh,
                a1face[i, j].fb,
                a1face[i+1, j].fw,
                a1face[i+1, j].fh,
                a1face[i+1, j].fb,
                a2face[i, j].fw,
                a2face[i, j].fh,
                a2face[i, j].fb,
                a2face[i, j+1].fw,
                a2face[i, j+1].fh,
                a2face[i, j+1].fb,
                ks.vSpace.u,
                ks.vSpace.v,
                ks.vSpace.weights,
                ks.gas.K,
                ks.gas.γ,
                ks.gas.μᵣ,
                ks.gas.ω,
                ks.gas.Pr,
                ctr[i, j].dx * ctr[i, j].dy,
                dt,
                zeros(4),
                zeros(4),
                :bgk,
            )
        end
    end
end
```

The result can be visualized with built-in function `plot_contour`, which presents the contours of gas density, U-velocity, V-velocity and temperature inside the cavity.
```julia
KitBase.plot_contour(ks, ctr)
```

![](./assets/cavity.png)

It is equivalent as the following low-level backend.
```julia
begin
    using Plots
    sol = zeros(4, ks.pSpace.nx, ks.pSpace.ny)
    for i in axes(sol, 2)
        for j in axes(sol, 3)
            sol[1:3, i, j] .= ctr[i, j].prim[1:3]
            sol[4, i, j] = 1.0 / ctr[i, j].prim[4]
        end
    end
    contourf(ks.pSpace.x[1:ks.pSpace.nx, 1], ks.pSpace.y[1, 1:ks.pSpace.ny], sol[3, :, :]')
end
```
# Phase Space

```@docs
newton_cotes
legendre_quadrature
octa_quadrature
quadrature_weights
```# Shock tube problem

We then use the Boltzmann equation to solve the shock tube problem in gas dynamics.
It's a two dimensional problem, with one in physical domain ``x`` and another in particle velocity domain ``u``.
First let us prepare the configuration file as
```
# case
matter = gas
case = sod
space = 1d2f1v
nSpecies = 1
flux = kfvs
collision = bgk
interpOrder = 2
limiter = vanleer
boundary = fix
cfl = 0.5
maxTime = 0.2

# physical space
x0 = 0
x1 = 1
nx = 200
pMeshType = uniform
nxg = 1

# velocity space
vMeshType = rectangle
umin = -5
umax = 5
nu = 28
nug = 0

# gas
knudsen = 0.0001
mach = 0.0
prandtl = 1
inK = 2
omega = 0.81
alphaRef = 1.0
omegaRef = 0.5
```

The configuration file can be understood as follows:
- The simulation case is the standard Sod shock tube
- A phase space in 1D physical and 1D velocity space is created with two particle distribution functions inside
- The numerical flux function is the kinetic flux vector splitting method and the collision term uses the BGK relaxation
- The reconstruction step employs van Leer limiter to create 2nd-order interpolation
- The two boundaries are fixed with Dirichlet boundary condition
- The timestep is determined with a CFL number of 0.5
- The maximum simulation time is 0.2
- The physical space spans in [0, 1] with 200 uniform cells
- The velocity space spans in [-5, 5] with 28 uniform cells
- The reference Knudsen number is set as 1e-4
- The reference Mach number is absent
- The reference Prandtl number is 1
- The gas molecule contains two internal degrees of freedom
- The viscosity is evaluated with the following formulas
```math
\mu = \mu_{ref} \left(\frac{T}{T_{ref}}\right)^{\omega}
```
```math
\mu_{ref}=\frac{5(\alpha+1)(\alpha+2) \sqrt{\pi}}{4 \alpha(5-2 \omega)(7-2 \omega)} Kn_{ref}
```

The configuration file directly generate variables during runtime via meta-programming in Julia,
and it can be stored in any text format (txt, toml, cfg, etc.). 
For example, if `config.txt` is created, 
we then execute the following codes to conduct a simulation
```julia
using Kinetic
set, ctr, face, t = initialize("config.txt")
t = solve!(set, ctr, face, t)
```

The computational setup is stored in `set` and the control volume solutions are stored in `ctr` and `face`. 
The high-level solver `solve!` is equivalent as the following low-level procedures
```julia
dt = timestep(ks, ctr, t)
nt = Int(floor(ks.set.maxTime / dt))
res = zeros(3)
for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, res)
end
```

The result can be visualized with built-in function `plot_line`, which presents the profiles of gas density, velocity and temperature inside the tube.
```julia
plot_line(set, ctr)
```
![](./assets/sod.png)
# Update

The update solver calculate the variables at n+1 step based on numerical fluxes and in-cell collisions.
```@docs
update!
```

The current solver supports different collision models, for example:
- `:bgk`: BGK relaxation model
- `:shakhov`: Shakhov relaxation model
- `:boltzmann`: Boltzmann: original Boltzmann collision integral

The boundary conditions vary.
- `:fix`: fixed Dirichlet boundary
- `:period`: periodic boundary
- `:extra`: extrapolation
- `:maxwell`: Maxwell's diffusive boundary

The current solver adopts implicit-explicit (IMEX) uniformly.
Further Multi-step time integrators can be used in conjunction with method of lines in [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
# I / O 

```@docs
read_dict
```# Configuration

```@docs
ib_rh
ib_sod
ib_briowu
ib_cavity
```# General framework

Kinetic employs the finite volume method (FVM) for modeling and simulation. 
The general solution algorithm can be conclude as follows, where both explicit and implicit methods are implemented.

![](./assets/solver_process.png)

The high-level solver function is 
```@docs
solve!
```

The detailed solution procedures can be concluded as follows
- pre-process
- timestep calculation
- reconstruction
- evolution
- update
- post-process
# Initial and boundary conditions

A struct `ib <: AbstractCondition` defines the initial and boundary conditions of a simulation.
It contains the values of conservative and primitive variables, and particle distribution functions at left and right (up and down) domain for both initial and boundary conditions.
It is set this way to easily deal with discontinuous initial conditions.

```@docs
IB
IB1F
IB2F
IB3F
IB4F
```# KitML and scientific machine learning

Machine learning is building its momentum in scientific computing.
Given the nonlinear structure of differential and integral equations, it is promising to incorporate the universal function approximator from machine learning models into the governing equations and achieve the balance between efficiency and accuracy.
KitML is designed as a scientific machine learning toolbox, which devotes to fusing mechanical and neural models.
For example, the Boltzmann collision operator can be divided into a combination of relaxation model and neural network, i.e. the so-called universal Boltzmann equation.
```math
\frac{df}{dt} = \int_{\mathcal{R}^{3}} \int_{\mathcal{S}^{2}} \mathcal{B}(\cos \beta, g)\left[f\left(\mathbf{u}^{\prime}\right) f\left(\mathbf{u}_{*}^{\prime}\right)-f(\mathbf{u}) f\left(\mathbf{u}_{*}\right)\right] d \mathbf{\Omega} d \mathbf{u}_{*} \simeq \nu(\mathcal{M}-f)+\mathrm{NN}_{\theta}(\mathcal{M}-f)
```
The UBE has the following benefits. 
First, it automatically ensures the asymptotic limits. 
Let's consider the Chapman-Enskog method for solving Boltzmann equation, where the distribution function is approximated with expansion series.
```math
f \simeq f^{(0)}+f^{(1)}+f^{(2)}+\cdots, \quad f^{(0)}=\mathcal{M}
```
Take the zeroth order truncation, and consider an illustrative multi-layer perceptron.
```math
\mathrm{NN}_{\theta}(x)=\operatorname{layer}_{n}\left(\ldots \text { layer }_{2}\left({\sigma}\left(\text { layer }_{1}(x)\right)\right)\right), \quad \operatorname{layer}(x)=w x
```
Given the zero input from ``M − f``, the contribution from collision term is absent, and the moment equation naturally leads to the Euler equations.
```math
\frac{\partial}{\partial t}\left(\begin{array}{c}
\rho \\
\rho \mathbf{U} \\
\rho E
\end{array}\right)+\nabla_{\mathbf{x}} \cdot\left(\begin{array}{c}
\rho \mathbf{U} \\
\rho \mathbf{U} \otimes \mathbf{U} \\
\mathbf{U}(\rho E+p)
\end{array}\right)=\int\left(\begin{array}{c}
1 \\
\mathbf{u} \\
\frac{1}{2} \mathbf{u}^{2}
\end{array}\right)\left(\mathcal{M}_{t}+\mathbf{u} \cdot \nabla_{\mathbf{x}} \mathcal{M}\right) d \mathbf{u}=0
```

KitML provides two functions to construct universal Boltzmann equation, and it works seamlessly with any modern ODE solver in [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
```@docs
ube_dfdt
ube_dfdt!
```

Besides, we provide an input convex neural network (ICNN) developed by Amos et al.

The neural network parameters are constrained such that the output of the network is a convex function of the inputs. 
The structure of the ICNN is shown as follows, and it allows for efficient inference via optimization over some inputs to the network given others, and can be applied to settings including structured prediction, data imputation, reinforcement learning, and others. 
It is important for entropy-based modelling, since the minimization principle works exclusively with convex function.

![](./assets/icnn.png)

```@docs
ICNNLayer
ICNNChain
```
Besides, we also provide scientific machine learning training interfaces and I/O methods.
They are consistent with both [Flux.jl](https://github.com/FluxML/Flux.jl) and [DiffEqFlux.jl](https://github.com/SciML/DiffEqFlux.jl) ecosystem.

```@docs
sci_train
sci_train!
load_data
save_model
```# Particle velocity space

A struct `vSpace <: AbstractSetup` defines the particle velocity setup of a simulation.
Structs for 1-3 dimensional particle velocity space are built.
```@docs
VSpace1D
VSpace2D
VSpace3D
```
It contains
- u0 (v0, w0): location of starting point
- u1 (v1, w1): location of ending point
- nu (nv, nw): number of cells in one direction
- u (v, w): locations of middle points of all cells
- du (dv, dw): intervals of all cell points
- weights: quadrature weights for numerical integral

Note that the one-dimensional velocity space can be used to handle 1-3 dimensional unstructured topology as well.
In addition, velocity space structs for multi-component substance are implemented.
```@docs
MVSpace1D
MVSpace2D
MVSpace3D
```

For the simulation cases where no phase-space evolution is involved, `vSpace` can be set as `nothing` directly.
```julia
vSpace = nothing
```# Parallel computing

Julia supports different categories of parallel programming natively:
- Asynchronous "tasks", or coroutines
- Multi-threading
- Distributed computing

Kinetic integrates the the latter two mechanism along with the CUDA-based GPU computing.
An initialization function is built in Kinetic.
```julia
function __init__()
    np = nworkers()
    nt = Threads.nthreads()
    if nt > 1 || np > 1
        @info "Kinetic will run with $np processors and $nt threads"
    else
        @info "Kinetic will run serially"
    end

    if has_cuda()
        @info "Kinetic will run with CUDA"
        for (i, dev) in enumerate(CUDA.devices())
            @info "$i: $(CUDA.name(dev))"
        end
        @info "Scalar operation is disabled in CUDA"
        CUDA.allowscalar(false)
    end
end
```
As the package is imported, it will report the computational resources (processors, threads and CUDA devices) that are going to be utilized.# Universal Boltzmann equation

In the following, we present a universal differential equation strategy to construct the neural network enhanced Boltzmann equation.
The complicated fivefold integral operator is replaced by a combination of relaxation and neural models.
It promises a completely differential structure and thus the neural ODE type training and computing becomes possible.
The approach reduces the computational cost up to three orders of magnitude and preserves the perfect accuracy.
The detailed theory and implementation can be found in [Tianbai Xiao and Martin Frank, Using neural networks to accelerate the solution of the Boltzmann equation](https://arxiv.org/pdf/2010.13649.pdf).

First we load all the packages needed and set up the configurations.
```julia
using OrdinaryDiffEq, Flux, DiffEqFlux, Plots
using KitBase, KitML

begin
    case = "homogeneous"
    maxTime = 3
    tlen = 16
    u0 = -5
    u1 = 5
    nu = 80
    nug = 0
    v0 = -5
    v1 = 5
    nv = 28
    nvg = 0
    w0 = -5
    w1 = 5
    nw = 28
    nwg = 0
    vMeshType = "rectangle"
    nm = 5
    knudsen = 1
    inK = 0
    alpha = 1.0
    omega = 0.5
    nh = 8
end
```

The dataset is produced by the fast spectral method, which solves the nonlinear Boltzmann integral with fast Fourier transformation.
```julia
begin
    tspan = (0.0, maxTime)
    tran = linspace(tspan[1], tspan[2], tlen)
    γ = heat_capacity_ratio(inK, 3)
    vSpace = VSpace3D(u0, u1, nu, v0, v1, nv, w0, w1, nw, vMeshType)

    f0 =
        Float32.(
            0.5 * (1 / π)^1.5 .*
            (exp.(-(vSpace.u .- 0.99) .^ 2) .+ exp.(-(vSpace.u .+ 0.99) .^ 2)) .*
            exp.(-vSpace.v .^ 2) .* exp.(-vSpace.w .^ 2),
        ) |> Array
    prim0 =
        conserve_prim(moments_conserve(f0, vSpace.u, vSpace.v, vSpace.w, vSpace.weights), γ)
    M0 = Float32.(maxwellian(vSpace.u, vSpace.v, vSpace.w, prim0)) |> Array

    mu_ref = ref_vhs_vis(knudsen, alpha, omega)
    kn_bzm = hs_boltz_kn(mu_ref, 1.0)
    τ0 = mu_ref * 2.0 * prim0[end]^(0.5) / prim0[1]

    phi, psi, phipsi = kernel_mode(
        nm,
        vSpace.u1,
        vSpace.v1,
        vSpace.w1,
        vSpace.du[1, 1, 1],
        vSpace.dv[1, 1, 1],
        vSpace.dw[1, 1, 1],
        vSpace.nu,
        vSpace.nv,
        vSpace.nw,
        alpha,
    )

    # Boltzmann
    prob = ODEProblem(boltzmann_ode!, f0, tspan, [kn_bzm, nm, phi, psi, phipsi])
    data_boltz = solve(prob, Tsit5(), saveat = tran) |> Array

    # BGK
    prob1 = ODEProblem(bgk_ode!, f0, tspan, [M0, τ0])
    data_bgk = solve(prob1, Tsit5(), saveat = tran) |> Array


    data_boltz_1D = zeros(Float32, axes(data_boltz, 1), axes(data_boltz, 4))
    data_bgk_1D = zeros(Float32, axes(data_bgk, 1), axes(data_bgk, 4))
    for j in axes(data_boltz_1D, 2)
        data_boltz_1D[:, j] .=
            reduce_distribution(data_boltz[:, :, :, j], vSpace.weights[1, :, :])
        data_bgk_1D[:, j] .=
            reduce_distribution(data_bgk[:, :, :, j], vSpace.weights[1, :, :])
    end
    f0_1D = reduce_distribution(f0, vSpace.weights[1, :, :])
    M0_1D = reduce_distribution(M0, vSpace.weights[1, :, :])

    X = Array{Float32}(undef, vSpace.nu, 1)
    for i in axes(X, 2)
        X[:, i] .= f0_1D
    end
    Y = Array{Float32}(undef, vSpace.nu, 1, tlen)
    for i in axes(Y, 2)
        Y[:, i, :] .= data_boltz_1D
    end
    M = Array{Float32}(undef, nu, size(X, 2))
    for i in axes(M, 2)
        M[:, i] .= M0_1D
    end
    τ = Array{Float32}(undef, 1, size(X, 2))
    for i in axes(τ, 2)
        τ[1, i] = τ0
    end
end
```

Then we define the neural network and construct the unified model with mechanical and neural parts.
The training is conducted by DiffEqFlux.jl with ADAM optimizer.
```julia
begin
    model_univ = DiffEqFlux.FastChain(
        DiffEqFlux.FastDense(nu, nu * nh, tanh),
        DiffEqFlux.FastDense(nu * nh, nu),
    )
    p_model = DiffEqFlux.initial_params(model_univ)

    function dfdt(f, p, t)
        df = (M .- f) ./ τ .+ model_univ(M .- f, p)
    end
    prob_ube = ODEProblem(dfdt, X, tspan, p_model)

    function loss(p)
        sol_ube = solve(prob_ube, Midpoint(), u0 = X, p = p, saveat = tran)
        loss = sum(abs2, Array(sol_ube) .- Y)

        return loss
    end

    his = []
    cb = function (p, l)
        display(l)
        push!(his, l)
        return false
    end
end

res = DiffEqFlux.sciml_train(loss, p_model, ADAM(), cb = cb, maxiters = 200)
res = DiffEqFlux.sciml_train(loss, res.minimizer, ADAM(), cb = cb, maxiters = 200)
```

Once we have trained a hybrid Boltzmann collision term, we could solve it as a normal differential equation with any desirable solvers.
Consider the Midpoint rule as an example, the solution algorithm and visualization can be organized.
```julia
ube = ODEProblem(KitML.ube_dfdt, f0_1D, tspan, [M0_1D, τ0, (model_univ, res.minimizer)]);
sol = solve(
    ube,
    Midpoint(),
    u0 = f0_1D,
    p = [M0_1D, τ0, (model_univ, res.minimizer)],
    saveat = tran,
);

plot(
    vSpace.u[:, vSpace.nv÷2, vSpace.nw÷2],
    data_boltz_1D[:, 1],
    lw = 2,
    label = "Initial",
    color = :gray32,
    xlabel = "u",
    ylabel = "particle distribution",
)
plot!(
    vSpace.u[:, vSpace.nv÷2, vSpace.nw÷2],
    data_boltz_1D[:, 2],
    lw = 2,
    label = "Boltzmann",
    color = 1,
)
plot!(
    vSpace.u[:, vSpace.nv÷2, vSpace.nw÷2],
    data_bgk_1D[:, 2],
    lw = 2,
    line = :dash,
    label = "BGK",
    color = 2,
)
plot!(
    vSpace.u[:, vSpace.nv÷2, vSpace.nw÷2],
    M0_1D,
    lw = 2,
    label = "Maxwellian",
    color = 10,
)
scatter!(vSpace.u[:, vSpace.nv÷2, vSpace.nw÷2], sol.u[2], lw = 2, label = "UBE", color = 3)
```

![](./assets/ube.png)# GPU computing

The thriving development of GPUs provides an alternative choice for scientific computing.
Kinetic enables computation on the graphical architecture on the basis of [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl).
It provides the main programming interface for working with NVIDIA CUDA GPUs. 
It features a user-friendly array abstraction, a compiler for writing CUDA kernels in Julia, and wrappers for various CUDA libraries.
In the following, we present an illustrative test of kinetic flux vector splitting method to evaluate upwind flux of the Boltzmann equation.
The test is conducted on a Tesla K80 GPU on [nextjournal.com](nextjournal.com).
We first load all the modules, and do a CPU-based computation.

```julia
import Pkg
Pkg.add("Revise")
Pkg.add("KitBase")
Pkg.add("CUDA")
Pkg.add("BenchmarkTools")

using Revise, CUDA, BenchmarkTools, KitBase

dt = 1e-3

primL = [1.0, 0.0, 0.5]
primR = [0.125, 0.0, 0.625]

u = collect(-5.0:0.01:5.0)
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)
```
The benchmark result on a Intel NUC8i7BEH desktop is around `5.244 μs (3 allocations: 24.00 KiB)`.
Then let's turn to GPU.
```julia
u = collect(-5.0:0.01:5.0) |> CuArray
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)
```
The benchmark result is around `32.965 μs (187 allocations: 10.73 KiB)`.
As can be seen, due to the relative small input size, the GPU threads aren't fully occupied, and therefore CPU is more efficient in this case.

Then let's increase the input vector size, i.e. to consider more discrete particle velocity points for distribution functions.
```julia
u = collect(-5.0:0.001:5.0)
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)

u = collect(-5.0:0.001:5.0) |> CuArray
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)
```
The results become around `50.011 μs (6 allocations: 234.80 KiB)` for CPU and `33.640 μs (187 allocations: 10.73 KiB)` for GPU.

We could further increase the computation size.
```julia
u = collect(-5.0:0.0001:5.0)
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)

u = collect(-5.0:0.0001:5.0) |> CuArray
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)
```
The results become around `507.960 μs (6 allocations: 2.29 MiB)` for CPU and `32.021 μs (187 allocations: 10.73 KiB)` for GPU.
Under this size of computation, the GPU brings about 16x efficiency increment.# Evolution

```@docs
evolve!
```

The evolution solver calculate the interface numerical fluxes based on two neighbor cells.
Different flux functions can be used with the option `model`.
- macroscopic: Godunov, Lax, Roe, HLL, wave-propagation
- mesoscopic: upwind, central-upwind, gas-kinetic scheme

The available flux solvers are listed as follows.
```@docs
flux_lax!
flux_hll!
flux_roe!
flux_gks
flux_gks!
flux_kfvs!
flux_kcu!
flux_ugks!
flux_boundary_maxwell!
flux_boundary_specular!
flux_em!
flux_emx!
flux_emy!
```
# Reference

- Bezanson, J., Edelman, A., Karpinski, S., & Shah, V. B. (2017). Julia: A fresh approach to numerical computing. SIAM review, 59(1), 65-98.
- Chapman, S., Cowling, T. G., & Burnett, D. (1990). The mathematical theory of non-uniform gases: an account of the kinetic theory of viscosity, thermal conduction and diffusion in gases. Cambridge University Press.
- Landau, L.D., & Lifshitz E. M. (1959). Fluid mechanics, London: Pergamon Press.
- Blazek, J. (2015). Computational fluid dynamics: principles and applications. Butterworth-Heinemann.
- Xu, K., & Huang, J. C. (2010). A unified gas-kinetic scheme for continuum and rarefied flows. Journal of Computational Physics, 229(20), 7747-7764.
- Bird, G. A. (1994). Molecular gas dynamics and the direct simulation of gas flows. Molecular gas dynamics and the direct simulation of gas flows.
- Xiao, T., Cai, Q., & Xu, K. (2017). A well-balanced unified gas-kinetic scheme for multiscale flow transport under gravitational field. Journal of Computational Physics, 332, 475-491.
- Xiao, T., Xu, K., & Cai, Q. (2019). A unified gas-kinetic scheme for multiscale and multicomponent flow transport. Applied Mathematics and Mechanics, 40(3), 355-372.
- Xiao, T., Liu, C., Xu, K., & Cai, Q. (2020). A velocity-space adaptive unified gas kinetic scheme for continuum and rarefied flows. Journal of Computational Physics, 415, 109535.
- Xiao, T., & Frank, M. (2020). Using neural networks to accelerate the solution of the Boltzmann equation. arXiv:2010.13649.
- Amos, B., Xu, L., & Kolter, J. Z. (2017, July). Input convex neural networks. In International Conference on Machine Learning (pp. 146-155). PMLR.# Preprocess

```@docs
initialize
```

The pre-process solver initializes the simulation that returns solver set, control volumes, interfaces, and current time.
It could be a new simulation or restart of an interrupted one.
- new run: .txt / .cfg / .toml / etc.
- restart: .jld2
# Parameter Settings

A struct `set <: AbstractSetup` defines the general setup of a simulation.
```@docs
Setup
```
It contains
- matter: fluid substance
- case: simulation case name
- space: ``n_1 d n_2 f n_3 v``, which denotes the physical dimensionality, numbers of particle distribution functions and velocity dimensionality
- flux: numerical flux function name
- collision: collision operator of kinetic equation
- nSpecies: number of species
- interpOrder: order of accuracy for reconstruction
- limiter: limiter function name
- boundary: boundary condition
- cfl: Courant-Friedrichs-Lewy number for determining time step
- maxTime: maximum simulation time# Postprocess

The post-process solver handles the simulation data and visualization.

```@docs
plot_line
plot_contour
write_jld
```
# Distributed computing

The distributed computation is built upon Julia's `@distributed` macro in the `Distributed` module.
```julia
using Distributed
@distributed [reducer] for var = range
    body
end
```
It provides a MPI-type parallelization with a leaner code size.
The specified range is partitioned and locally executed across all workers. 
In case an optional reducer function is specified, `@distributed` performs local reductions on each worker with a final reduction on the calling process.
Without a reducer function, @distributed will execute asynchronously, i.e. it spawns independent tasks on all available workers and returns immediately without waiting for completion. 
To make it wait for completion, prefix the call with `@sync` like :
```julia
@sync @distributed for var = range
    body
end
```

In the following, we present an example to conduct distributed computing with the help Julia's `SharedArrays` module, which creates arrays shared by all the processors.
More massive computing can be made with [DistributedArrays](https://github.com/JuliaParallel/DistributedArrays.jl).
First, we consider a distributed computing.

```julia
using Distributed, SharedArrays
addprocs(3)
@everywhere using KitBase

begin
    vars = Dict{Symbol,Any}()
    vars[:matter] = "gas"
    vars[:case] = "sod"
    vars[:space] = "1d0f0v"
    vars[:flux] = "kfvs"
    vars[:collision] = "bgk"
    vars[:nSpecies] = 1
    vars[:interpOrder] = 1
    vars[:limiter] = "vanleer"
    vars[:boundary] = "fix"
    vars[:cfl] = 0.5
    vars[:maxTime] = 0.2
    vars[:x0] = 0.0
    vars[:x1] = 1.0
    vars[:nx] = 2000
    vars[:pMeshType] = "uniform"
    vars[:nxg] = 0
    vars[:knudsen] = 0.001
    vars[:mach] = 0.0
    vars[:prandtl] = 1.0
    vars[:inK] = 0.0
    vars[:omega] = 0.81
    vars[:alphaRef] = 1.0
    vars[:omegaRef] = 0.5
end

set = KitBase.set_setup(vars)
pSpace = KitBase.set_geometry(vars)
vSpace = KitBase.set_velocity(vars)
gas = KitBase.set_property(vars)
ib = KitBase.set_ib(vars, set, vSpace, gas)
folder = @__DIR__
ks = KitBase.SolverSet(set, pSpace, vSpace, gas, ib, folder)

dt = ks.pSpace.dx[1] / (5.0 + KitBase.sound_speed(ks.ib.primL, ks.gas.γ))
nt = floor(ks.set.maxTime / dt) |> Int

wp = SharedArray{Float64}((ks.pSpace.nx, 3), init=A->(A=zeros(ks.pSpace.nx, 3)))
for i in 1:ks.pSpace.nx
    if i <= ks.pSpace.nx ÷ 2
        wp[i,:] .= ks.ib.wL
    else
        wp[i,:] .= ks.ib.wR
    end
end     
fwp = SharedArray{Float64}((ks.pSpace.nx+1, 3), init=A->(A=zeros(ks.pSpace.nx+1, 3)))

@time for iter = 1:nt÷3
    @sync @distributed for i in 2:ks.pSpace.nx
        flux = @view fwp[i,:]
        KitBase.flux_gks!(
            flux,
            wp[i-1,:],
            wp[i,:],
            ks.gas.γ,
            ks.gas.K,
            ks.gas.μᵣ,
            ks.gas.ω,
            dt,
            0.5 * ks.pSpace.dx[i-1],
            0.5 * ks.pSpace.dx[i],
        )
    end
    
    @sync @distributed for i in 2:ks.pSpace.nx-1
        for j in 1:3
            wp[i,j] += (fwp[i,j] - fwp[i+1,j]) / ks.pSpace.dx[i]
        end
    end
end
```

The benchmark result on a Intel NUC8i7BEH desktop is around `13.620491 seconds (2.26 M allocations: 101.219 MiB, 0.22% gc time)`.
Then, we compare the efficiency with a serial execution.

```julia
w = zeros(ks.pSpace.nx, 3)
for i in 1:ks.pSpace.nx
    if i <= ks.pSpace.nx ÷ 2
        w[i,:] .= ks.ib.wL
    else
        w[i,:] .= ks.ib.wR
    end
end     
fw = zeros(ks.pSpace.nx+1, 3)

@time for iter = 1:nt÷3
    for i in 2:ks.pSpace.nx
        flux = @view fw[i,:]
        KitBase.flux_gks!(
            flux,
            w[i-1,:],
            w[i,:],
            ks.gas.γ,
            ks.gas.K,
            ks.gas.μᵣ,
            ks.gas.ω,
            dt,
            0.5 * ks.pSpace.dx[i-1],
            0.5 * ks.pSpace.dx[i],
        )
    end
    
    for i in 2:ks.pSpace.nx-1
        for j in 1:3
            w[i,j] += (fw[i,j] - fw[i+1,j]) / ks.pSpace.dx[i]
        end
    end
end
```

The result on the same desktop is around `20.830331 seconds (323.96 M allocations: 24.472 GiB, 16.89% gc time)`.
With more grid cells being used, the performance deviation is expected to be more significant.# Calling from Python 

For maximum convenience, a wrapper [kineticpy](https://github.com/vavrines/kineticpy) has been built to locate all the methods from Python.

## How to use?

Let's start by cloning the repository and changing into the directory.
```bash
git clone https://github.com/vavrines/kineticpy.git
cd kineticpy
```

Next, we start `python`.
The Julia main module can be installed and initialized by
```python
>>> import kineticpy
>>> kineticpy.install()
```

The basic structs and methods are stored in the base module, and can be imported via
```python
>>> from kineticpy import base
```

## Example

We provide some quick tutorial here for kineticpy.

```python
>>> from kineticpy import base
>>> import numpy as np

>>> u = np.linspace(-5, 5, 28) # velocity space
>>> prim_var = np.array([1.0, 0.0, 1.0]) # primitive flow variables
>>> M = base.maxwellian(u, prim_var) # compute Maxwellian distribution
>>> M.view()

array([7.83543327e-12, 2.77323769e-10, 7.46041809e-09, 1.52542631e-07,
       2.37067103e-06, 2.80029217e-05, 2.51412806e-04, 1.71562923e-03,
       8.89839075e-03, 3.50793472e-02, 1.05109877e-01, 2.39379825e-01,
       4.14365469e-01, 5.45169515e-01, 5.45169515e-01, 4.14365469e-01,
       2.39379825e-01, 1.05109877e-01, 3.50793472e-02, 8.89839075e-03,
       1.71562923e-03, 2.51412806e-04, 2.80029217e-05, 2.37067103e-06,
       1.52542631e-07, 7.46041809e-09, 2.77323769e-10, 7.83543327e-12])
```# Timestep

```@docs
timestep
```

The timestep solver returns the time interval used for the upcoming solution loop based on the current variables.
# Physical Space

```@docs
global_frame
local_frame
uniform_mesh
meshgrid
mesh_connectivity_2D
mesh_center_2D
mesh_area_2D
```# Stepper

```@docs
KitBase.step!
```
# Basic Physics

## Microscopic formulation

The physical world shows a diverse set of behaviors on different characteristic scales.
Consider the molecular motion of gases as an example.
Down to the finest scale of a many-particle system, the Newton's second law depicts particle motions via
```math
\mathbf{F} = m \mathbf{a}.
```

As a first order system it reads
```math
\frac{d \mathbf x}{dt} = \mathbf v, \ \frac{d \mathbf v}{dt} = \frac{\mathbf F}{m},
```
where ``\mathbf F`` is external force and ``m`` is particle mass.

An intuitive numerical algorithm is to get the numerous particles on board and track the trajectories of them.
A typical example is the [Molecular Dynamics](https://en.wikipedia.org/wiki/Molecular_dynamics).
This is not going to be efficient since there are more than `2e25` molecules per cubic meter in normal atmosphere, and things get even more complicated when you count on the N-body interactions all the time.
Some methods have been proposed to simplify the computation.
As an example, the [Direct simulation Monte Carlo](https://en.wikipedia.org/wiki/Direct_simulation_Monte_Carlo) employs certain molecular models and conduct the intermolecular collisions in a stochastic manner.
It significantly reduces the computational cost, while the trade-off is the artificial fluctuations.
Many realizations must be simulated successively to average the solutions and reduce the errors.

## Mesoscopic formulation

An alternative strategy is made from ensemble averaging, where the coarse-grained modeling is used to provide a bottom-up view.
At the mean free path and collision time scale of molecules, particles travel freely during most of time with mild intermolecular collisions.
Such dynamics can be described with an operator splitting approach, i.e. the kinetic transport equation
```math
\frac{\partial f}{\partial t}+ \mathbf v \cdot \nabla_\mathbf x f + \mathbf a \cdot \nabla_\mathbf v f = Q(f),
```
where ``f`` denotes the probability of finding a particle at certain location in phase space.
The left-hand side of the equation above model the transport phenomena due to the inhomogeneous distribution of particles and external force field,
while the right-hand side depicts intermolecular collision.
Different collision models can be inserted into such equation.
If the particles only collide with a background material one obtains linear Boltzmann collision operator
```math
Q(f)=\int_{\mathbb R^3} \mathcal B(\mathbf v_*, \mathbf v) \left[ f(\mathbf v_*)-f(\mathbf v)\right] d\mathbf v_*,
```
where the collision kernel ``\mathcal B`` models the strength of collisions at different velocities. 
If the interactions among particles are considered, the collision operator becomes nonlinear. 
For example, the two-body collision results in nonlinear Boltzmann equation
```math
Q(f)=\int_{\mathbb R^3} \int_{\mathcal S^2} \mathcal B(\cos \beta, |\mathbf{v}-\mathbf{v_*}|) \left[ f(\mathbf v')f(\mathbf v_*')-f(\mathbf v)f(\mathbf v_*)\right] d\mathbf \Omega d\mathbf v_*.
```
To solve the Boltzmann equation, a discretized phase space needs to be introduced and the solution algorithm is called [discrete ordinates method](https://en.wikipedia.org/wiki/Discrete_ordinates_method) or discrete velocity method.
Due to the complicated fivefold integral in the nonlinear Boltzmann collision operator, sometimes it is replaced by the simplified models in the discrete velocity method, e.g. the relaxation model
```math
Q(f) = \nu (\mathcal M - f).
```
From the [H-theorem](https://en.wikipedia.org/wiki/H-theorem), we learn that an isolated system evolves in the direction with entropy increment.
The maximal entropy status corresponds to the well-known Maxwellian distribution
```math
\mathcal M = n\left(\frac{m}{2\pi k T}\right)^{D/2}\exp \left( -\frac{m}{2kT} (\mathbf v - \mathbf V)^2 \right),
```
where ``k`` is the Boltzmann constant, ``\mathbf V`` and ``T`` are the bulk velocity and temperature.
The Boltzmann dynamics can be projected onto lower dimensionality.
For example, with one-dimensional velocity space formulation, the high-dimensional particle distribution can be integrated with respect to the rest coordinates as
```math
h_0 = \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} f d v dw, \ h_1 = \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} (v^2+w^2) f dv dw,
```
where ``h_0`` and ``h_1`` are called reduced distribution functions and form a so-called `1d2f1v` system.

## Macroscopic formulation

Meanwhile, with the enlargement of modeling scale to a macroscopic hydrodynamic level, the accumulating effect of particle collisions results in an equalization of local temperature and velocity,
where the moderate non-equilibrium effects can be well described by viscous transport, heat conduction and mass diffusion,
i.e., the so called transport phenomena. 
Large-scale dynamics presents the property of waves, and the macroscopic transport equations can be constructed to describe the bulk behaviors of fluids.
Typical examples are the Euler and Navier-Stokes equations
```math
\frac{\partial \mathbf W}{\partial t} + \nabla_\mathbf x \cdot \mathbf F = \mathbf S
```
From microscopic particle transport to macroscopic fluid motion, there is a continuous variation of flow dynamics. 
We pay special attentions to Hilbert's sixth problem, i.e. building the numerical passage between the kinetic theory of gases and continuum mechanics. 
# Burgers

Now we could turn to the Burgers equation.
It's a typical hyperbolic conservation law, where discontinuous solution can emerge in a self-evolving system.
Let's consider the same initial configuration as advection-diffusion example.
```julia
using Kinetic, Plots

set = Setup(
    matter = "scalar", # material
    case = "burgers", # test case
    space = "1d0f0v", # phase space
    flux = "gks", # flux
    collision = "", # collision: for scalar conservation laws there are none
    interpOrder = 1, # interpolation order
    boundary = "period", # boundary condition
    cfl = 0.5, # cfl
    maxTime = 1.0, # simulation time
)
ps = PSpace1D(0.0, 1.0, 100, 1)
vs = nothing
property = Scalar(0, 1e-6)
ib = IB(x -> sin(2π * x), property)

ks = SolverSet(set, ps, vs, property, ib)
ctr, face = init_fvm(ks)
```

The solution algorithm can be processed together with visualization.
```julia
t = 0.0
dt = KitBase.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int

anim = @animate for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, 0.0)

    plot(ks, ctr, xlabel="x", label="u", ylims=[-1, 1])
end

gif(anim, "burgers.gif", fps = 45)
```

![](./assets/burgers.gif)
# Particle properties

A struct `gas <: AbstractProperty` defines the properties of particle model.
It currently supports the following models:
- scalar
- gas-type molecule
- plasma

```@docs
Scalar
Gas
Mixture
Plasma1D
Plasma2D
```

The fields denote, for example:
- Kn: reference Knudsen number
- Ma: reference Mach number
- Pr: reference Prandtl number
- K: internal degree of freedom of molecule
- γ: adiabatic index
- ω: viscosity index
- αᵣ: reference ``\alpha`` in viscosity evaluation
- ωᵣ: reference ``\omega`` in viscosity evaluation
- μᵣ: reference viscosity
- m: mass of each particle
- np: number of particles

The viscosity is evaluated the following hard-sphere model.
```math
\mu = \mu_{ref} \left(\frac{T}{T_{ref}}\right)^{\omega}
```
```math
\mu_{ref}=\frac{5(\alpha+1)(\alpha+2) \sqrt{\pi}}{4 \alpha(5-2 \omega)(7-2 \omega)} Kn_{ref}
```
# Theory

```@docs
prim_conserve
conserve_prim
mixture_prim_conserve
mixture_conserve_prim
em_coefficients
advection_flux
burgers_flux
euler_flux
euler_jacobi
gauss_moments
mixture_gauss_moments
moments_conserve
mixture_moments_conserve
pdf_slope
mixture_pdf_slope
moments_conserve_slope
mixture_moments_conserve_slope
discrete_moments
stress
heat_flux
maxwellian
mixture_maxwellian
shakhov
reduce_distribution
full_distribution
ref_vhs_vis
vhs_collision_time
aap_hs_collision_time
aap_hs_prim
aap_hs_diffeq!
shift_pdf!
hs_boltz_kn
kernel_mode
boltzmann_fft
boltzmann_fft!
heat_capacity_ratio
sound_speed
```