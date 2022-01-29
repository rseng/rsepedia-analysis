# GeophysicalFlows.jl

<!-- description -->
<p>
  <strong>💨🌏🌊 Geophysical fluid dynamics pseudospectral solvers with Julia and <a href="http://github.com/FourierFlows/FourierFlows.jl">FourierFlows.jl</a>. https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable</strong>
</p>

<!-- Badges -->
<p align="left">
    <a href="https://buildkite.com/julialang/geophysicalflows-dot-jl">
        <img alt="Buildkite CPU+GPU build status" src="https://img.shields.io/buildkite/4d921fc17b95341ea5477fb62df0e6d9364b61b154e050a123/main?logo=buildkite&label=Buildkite%20CPU%2BGPU">
    </a>
    <a href="https://ci.appveyor.com/project/navidcy/geophysicalflows-jl">
        <img alt="Build Status for Window" src="https://img.shields.io/appveyor/ci/navidcy/geophysicalflows-jl/main?label=Window&logo=appveyor&logoColor=white">
    </a>
    <a href="https://FourierFlows.github.io/GeophysicalFlowsDocumentation/stable">
        <img alt="stable docs" src="https://img.shields.io/badge/documentation-stable%20release-blue">
    </a>
    <a href="https://FourierFlows.github.io/GeophysicalFlowsDocumentation/dev">
        <img alt="latest docs" src="https://img.shields.io/badge/documentation-in%20development-orange">
    </a>
    <a href="https://codecov.io/gh/FourierFlows/GeophysicalFlows.jl">
        <img src="https://codecov.io/gh/FourierFlows/GeophysicalFlows.jl/branch/main/graph/badge.svg" />
    </a>
    <a href="https://doi.org/10.5281/zenodo.1463809">
        <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1463809.svg" alt="DOI">
    </a>
    <a href="https://github.com/SciML/ColPrac">
      <img alt="ColPrac: Contributor's Guide on Collaborative Practices for Community Packages" src="https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet">
    </a>
    <a href="https://doi.org/10.21105/joss.03053">
      <img src="https://joss.theoj.org/papers/10.21105/joss.03053/status.svg" alt="DOI badge" >
    </a>
</p>

This package leverages the [FourierFlows.jl] framework to provide modules for solving problems in
Geophysical Fluid Dynamics on periodic domains using Fourier-based pseudospectral methods.


## Installation

To install, use Julia's  built-in package manager (accessed by pressing `]` in the Julia REPL command prompt) to add the package and also to instantiate/build all the required dependencies

```julia
julia>]
(v1.6) pkg> add GeophysicalFlows
(v1.6) pkg> instantiate
```

The most recent version of GeophysicalFlows.jl requires Julia v1.6  (the current long-term-release) or later. _We strongly urge you to use this version._

The latest version that is compatible with Julia v1.5 is GeophysicalFlows.jl v0.13.1.


## Examples

See `examples/` for example scripts. These examples are best viewed by browsing them within 
the package's [documentation]. 

Some animations created with GeophysicalFlows.jl are [online @ youtube].


## Modules

* `TwoDNavierStokes`: the two-dimensional vorticity equation.
* `SingleLayerQG`: the barotropic or equivalent-barotropic quasi-geostrophic equation, which 
  generalizes `TwoDNavierStokes` to cases with topography, Coriolis parameters of the form 
  `f = f₀ + βy`, and finite Rossby radius of deformation.
* `MultiLayerQG`: a multi-layer quasi-geostrophic model over topography and with the ability 
  to impose a zonal flow `U_n(y)` in each layer.
* `SurfaceQG`: a surface quasi-geostrophic model.
* `BarotropicQGQL`: the quasi-linear barotropic quasi-geostrophic equation.


## Scalability

For now, GeophysicalFlows.jl is restricted to run on either a single CPU or single GPU. These
restrictions come from FourierFlows.jl. Multi-threading can enhance performance for the Fourier
transforms. By default, FourierFlows.jl will use the maximum number of threads available on 
your machine. You can set the number of threads used by FourierFlows.jl by setting the 
environment variable, e.g.,

```
$ export JULIA_NUM_THREADS=4
```

For more information on multi-threading users are directed to the [Julia Documentation](https://docs.julialang.org/en/v1/manual/multi-threading/).

If your machine has more than one GPU available, then functionality within CUDA.jl package 
enables the user to choose the GPU device that FourierFlows.jl should use. The user is referred
to the [CUDA.jl Documentation](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#Device-Management);
in particular, [`CUDA.devices`](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#CUDA.devices) 
and [`CUDA.CuDevice`](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#CUDA.CuDevice). 
The user is also referred to the [GPU section](https://fourierflows.github.io/FourierFlowsDocumentation/stable/gpu/) in the FourierFlows.jl documentation.


## Getting help

Interested in GeophysicalFlows.jl or trying to figure out how to use it? Please feel free to 
ask us questions and get in touch! Check out the 
[examples](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/main/examples) and 
[open an issue](https://github.com/FourierFlows/GeophysicalFlows.jl/issues/new) or 
[start a discussion](https://github.com/FourierFlows/GeophysicalFlows.jl/discussions/new) 
if you have any questions, comments, suggestions, etc.


## Citing

If you use GeophysicalFlows.jl in research, teaching, or other activities, we would be grateful 
if you could mention GeophysicalFlows.jl and cite our paper in JOSS:

> Constantinou et al., (2021). GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs. _Journal of Open Source Software_, **6(60)**, 3053, doi:[10.21105/joss.03053](https://doi.org/10.21105/joss.03053).

The bibtex entry for the paper is:

```bibtex
@article{GeophysicalFlowsJOSS,
  doi = {10.21105/joss.03053},
  url = {https://doi.org/10.21105/joss.03053},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {60},
  pages = {3053},
  author = {Navid C. Constantinou and Gregory LeClaire Wagner and Lia Siegelman and Brodie C. Pearson and André Palóczy},
  title = {GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs \& GPUs},
  journal = {Journal of Open Source Software}
}
```


## Contributing

If you're interested in contributing to the development of GeophysicalFlows.jl we are excited 
to get your help, no matter how big or small a contribution you make! It's always great to have 
new people look at the code with fresh eyes: you will see errors that other developers have missed.

Let us know by [open an issue](https://github.com/FourierFlows/GeophysicalFlows.jl/issues/new) 
or [start a discussion](https://github.com/FourierFlows/GeophysicalFlows.jl/discussions/new) 
if you'd like to work on a new feature or implement a new module, if you're new to open-source 
and want to find a cool little project or issue to work on that fits your interests! We're more 
than happy to help along the way.

For more information, check out our [contributor's guide](https://github.com/FourierFlows/GeophysicalFlows.jl/blob/main/CONTRIBUTING.md).


[FourierFlows.jl]: https://github.com/FourierFlows/FourierFlows.jl
[documentation]: https://fourierflows.github.io/GeophysicalFlowsDocumentation/dev/
[online @ youtube]: https://www.youtube.com/channel/UCO_0ugkNUwCsFUMtepwYTqw
# Contributors Guide

Thank you for considering contributing to GeophysicalFlows.jl! 

Please feel free to ask us questions and chat with us at any time if you're
unsure about anything.

Best way to get in touch is to either just open a GitHub [issue](https://github.com/FourierFlows/GeophysicalFlows.jl/issues) 
(don't be shy!) or starting a [discussion](https://github.com/FourierFlows/GeophysicalFlows.jl/discussions).

We follow the [ColPrac guide](https://github.com/SciML/ColPrac) for collaborative
practices. New contributors should make sure to read that guide.

## What can I do?

* Tackle an existing [issue](https://github.com/FourierFlows/GeophysicalFlows.jl/issues).

* Try to run your favorite GeophysicalFlows.jl module and play around with it to simulate 
  your particular favorite setup. If you run into any problems or find it difficult
  to use, modify, or understand, please [open an issue](https://github.com/FourierFlows/GeophysicalFlows.jl/issues)!

* Write up an example or tutorial on how to do something useful with one of the current modules
  in GeophysicalFlows.jl, like how to set up a new physical configuration.

* Improve documentation, docstrings, or comments if you found something is hard to use.

* Implement a new feature (e.g., a new diagnostic into a module).

* Implement a new module from scratch to solve your favorite partial differential equation with
  periodic boundary conditions.

If you're interested in working on something, let us know by commenting on existing issues or 
by opening a new issue if. This is to make sure no one else is working on the same issue and 
so we can help and guide you in case there is anything you need to know beforehand.
---
title: 'GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs'
tags:
  - geophysical fluid dynamics
  - computational fluid dynamics
  - Fourier methods
  - pseudospectral
  - Julia
  - gpu
authors:
  - name: Navid C. Constantinou
    orcid: 0000-0002-8149-4094
    affiliation: "1, 2"
  - name: Gregory LeClaire Wagner
    orcid: 0000-0001-5317-2445
    affiliation: 3
  - name: Lia Siegelman
    orcid: 0000-0003-3330-082X
    affiliation: 4
  - name: Brodie C. Pearson
    orcid: 0000-0002-0202-0481
    affiliation: 5
  - name: André Palóczy
    orcid: 0000-0001-8231-8298
    affiliation: 6
affiliations:
 - name: Australian National University
   index: 1
 - name: ARC Centre of Excellence for Climate Extremes
   index: 2
 - name: Massachussetts Institute of Technology
   index: 3
 - name: University of California San Diego
   index: 4
 - name: Oregon State University
   index: 5
 - name: University of Oslo
   index: 6
date: 7 April 2021
bibliography: paper.bib
---


# Summary

`GeophysicalFlows.jl` is a Julia [@Bezanson2017] package that contains partial differential 
equations solvers for a collection of geophysical fluid systems in periodic domains. All 
modules use Fourier-based pseudospectral numerical methods and leverage the framework provided 
by the `FourierFlows.jl` [@FourierFlows] Julia package for time-stepping, custom diagnostics, 
and saving output.


# Statement of need

Conceptual models in simple domains often provide stepping stones for better understanding geophysical and astrophysical systems, particularly the atmospheres and oceans of Earth and other planets. These conceptual models are used in research but also are of great value for helping students in class to grasp new concepts and phenomena. Oftentimes people end up coding their own versions of solvers for the same partial differential equations for research or classwork. `GeophysicalFlows.jl` package is designed to be easily utilized and adaptable for a wide variety of both research and pedagogical purposes.

On top of the above-mentioned needs, the recent explosion of machine-learning applications in atmospheric and oceanic sciences advocates for the need that solvers for partial differential equations can be run on GPUs. 

`GeophysicalFlows.jl` provides a collection of modules for solving sets of partial differential equations often used as conceptual models. These modules are continuously tested (unit tests and tests for the physics involved) and are well-documented. `GeophysicalFlows.jl` utilizes Julia's functionality and abstraction to enable all modules to run on CPUs or GPUs, and to provide a high level of customizability within modules. The abstractions allow simulations to be tailored for specific research questions, via the choice of parameters, domain properties, and schemes for damping, forcing, time-stepping etc. Simulations can easily be carried out on different computing architectures. Selection of the architecture on which equations are solved is done by providing the argument `CPU()` or `GPU()` during the construction of a particular problem.
 
Documented examples for each geophysical system (module) appear in the package's documentation, 
providing a starting point for new users and for the development of new or customized modules. 
Current modules include two-dimensional flow and a variety of quasi-geostrophic (QG) dynamical 
systems, which provide analogues to the large-scale dynamics of atmospheres and oceans. The QG 
systems currently in `GeophysicalFlows.jl` extend two-dimensional dynamics to include the leading
order effects of a third dimension through planetary rotation, topography, surface boundary 
conditions, stratification and quasi-two-dimensional layering. A community-based collection 
of diagnostics throughout the modules are used to compute quantities like energy, enstrophy, 
dissipation, etc.

![Potential vorticity snapshots from a nonlinearly equilibrated simulation of the Eady instability 
over a meridional ridge. Simulation used `MultiLayerQG` module of `GeophysicalFlows.jl`. The Eady 
problem was approximated here using 5 fluid layers stacked up in the vertical. Each layer was
simulated with 512² grid-points. Plots were made with the `Plots.jl` Julia package, which 
utilizes the `cmocean` colormaps collection [@Thyng2016]. Scripts to reproduce the simulation 
reside in the repository `github.com/FourierFlows/MultilayerQG-example`.
\label{fig1}](PV_eady_nlayers5.png)


# State of the field

`GeophysicalFlows.jl` is a unique Julia package that shares some features and similarities with 
other packages. In particular:

- `pyqg` [@pyqg] (Python)

  Beyond their base language, the major differences between `GeophysicalFlows.jl` and `pyqg` 
  is that `GeophysicalFlows.jl` can be run on GPUs or CPUs and leverages a separate package (`FourierFlows.jl`; which is continuously developed) to solve differential equations and compute diagnostics, while `pyqg` can only be run on CPUs and uses a self-contained kernel. 
  
- Dedalus [@Burns2020] (Python)
  
  Dedalus is a Python package with an intuitive script-based interface that uses spectral methods 
  to solve general partial differential equations, such as the ones within `GeophysicalFlows.jl`.
  Dedalus allows for more general boundary conditions in one of the dimensions. It only runs on 
  CPUs (not on GPUs) but can be MPI-parallelized.
  
- `Oceananigans.jl` [@Oceananigans] (Julia)
  
  `Oceananigans.jl` is a fluid solver focussed on the Navier-Stokes equations under the Boussinesq
  approximation. `Oceananigans.jl` also runs on GPUs, and it allows for more variety of boundary
  conditions but it does not have spectral accuracy as it uses finite-volume discretization methods.
  
- `MAOOAM` [@MAOOAM] (Fortran, Python, and Lua) and its expanded Python implementation `qgs` [@qgs]

  `MAOOAM` and `qgs` simulate two atmospheric layers with QG dynamics, above either land or 
  an oceanic fluid layer with reduced-gravity QG dynamics. The dynamics of individual layers 
  have overlap with the `MultiLayerQG` and `SingleLayerQG` modules, however the layer configuration 
  of `MOAAM` and `qgs` is specifically designed to study the dynamics of Earth's mid-latitude 
  atmosphere. Neither `MAOOAM` nor `qgs` can run on GPUs.
  
- Isolated codes/scripts 

  Several codes/scripts exist in personal websites and in open-source public repositories with
  similar functionality as some `GeophysicalFlows.jl` modules (e.g., `TwoDNavierStokes` or 
  `SingleLayerQG`). Usually, though, these codes come without any or poor documentation and 
  typically they are not continuously tested.

`GeophysicalFlows.jl` can be used to investigate a variety of scientific research questions 
thanks to its various modules and high customizability, and its ease-of-use makes it an ideal 
teaching tool for fluids courses [@GeophysicalFlows-Examples; @CLExWinterSchool2020]. 
`GeophysicalFlows.jl` has been used in developing Lagrangian vortices identification algorithms 
[@Karrasch2020] and to test new theories for diagnosing turbulent energy transfers in geophysical 
flows [@Pearson2021]. Currently, `GeophysicalFlows.jl` is being used, e.g., (i) to compare 
different observational sampling techniques in these flows, (ii) to study the bifurcation properties 
of Kolmogorov flows [@KolmogorovFlow], (iii) to study the genesis and persistence of the polygons 
of vortices present at Jovian high latitudes (Siegelman, Young, and Ingersoll; in prep), and 
(iv) to study how mesoscale macroturbulence affects mixing of tracers [@QG_tracer_advection].


# Acknowledgements

We acknowledge discussions with Keaton Burns, Valentin Churavy, Theodore Drivas, Cesar Rocha, 
and William Young. B. C. P. was supported by the National Science Foundation under Grant 
No. 2023721. We would also like to take a moment to remember our friend and colleague 
Sean R. Haney (February 1987 - January 2021) who left us a bit too early.


# References
# GeophysicalFlows.jl/examples


These are some basic examples of the various modules included in GeophysicalFlows.jl. The best way to go through the examples by browsing them within the package's <a href="https://fourierflows.github.io/GeophysicalFlowsDocumentation/dev/">documentation <img src="https://img.shields.io/badge/docs-dev-blue.svg"></a>.


## Run the examples

You can run the examples in an executable environment via [Binder](https://mybinder.org) by clicking on the [![badge](https://img.shields.io/badge/binder-badge-579ACA.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFkAAABZCAMAAABi1XidAAAB8lBMVEX///9XmsrmZYH1olJXmsr1olJXmsrmZYH1olJXmsr1olJXmsrmZYH1olL1olJXmsr1olJXmsrmZYH1olL1olJXmsrmZYH1olJXmsr1olL1olJXmsrmZYH1olL1olJXmsrmZYH1olL1olL0nFf1olJXmsrmZYH1olJXmsq8dZb1olJXmsrmZYH1olJXmspXmspXmsr1olL1olJXmsrmZYH1olJXmsr1olL1olJXmsrmZYH1olL1olLeaIVXmsrmZYH1olL1olL1olJXmsrmZYH1olLna31Xmsr1olJXmsr1olJXmsrmZYH1olLqoVr1olJXmsr1olJXmsrmZYH1olL1olKkfaPobXvviGabgadXmsqThKuofKHmZ4Dobnr1olJXmsr1olJXmspXmsr1olJXmsrfZ4TuhWn1olL1olJXmsqBi7X1olJXmspZmslbmMhbmsdemsVfl8ZgmsNim8Jpk8F0m7R4m7F5nLB6jbh7jbiDirOEibOGnKaMhq+PnaCVg6qWg6qegKaff6WhnpKofKGtnomxeZy3noG6dZi+n3vCcpPDcpPGn3bLb4/Mb47UbIrVa4rYoGjdaIbeaIXhoWHmZYHobXvpcHjqdHXreHLroVrsfG/uhGnuh2bwj2Hxk17yl1vzmljzm1j0nlX1olL3AJXWAAAAbXRSTlMAEBAQHx8gICAuLjAwMDw9PUBAQEpQUFBXV1hgYGBkcHBwcXl8gICAgoiIkJCQlJicnJ2goKCmqK+wsLC4usDAwMjP0NDQ1NbW3Nzg4ODi5+3v8PDw8/T09PX29vb39/f5+fr7+/z8/Pz9/v7+zczCxgAABC5JREFUeAHN1ul3k0UUBvCb1CTVpmpaitAGSLSpSuKCLWpbTKNJFGlcSMAFF63iUmRccNG6gLbuxkXU66JAUef/9LSpmXnyLr3T5AO/rzl5zj137p136BISy44fKJXuGN/d19PUfYeO67Znqtf2KH33Id1psXoFdW30sPZ1sMvs2D060AHqws4FHeJojLZqnw53cmfvg+XR8mC0OEjuxrXEkX5ydeVJLVIlV0e10PXk5k7dYeHu7Cj1j+49uKg7uLU61tGLw1lq27ugQYlclHC4bgv7VQ+TAyj5Zc/UjsPvs1sd5cWryWObtvWT2EPa4rtnWW3JkpjggEpbOsPr7F7EyNewtpBIslA7p43HCsnwooXTEc3UmPmCNn5lrqTJxy6nRmcavGZVt/3Da2pD5NHvsOHJCrdc1G2r3DITpU7yic7w/7Rxnjc0kt5GC4djiv2Sz3Fb2iEZg41/ddsFDoyuYrIkmFehz0HR2thPgQqMyQYb2OtB0WxsZ3BeG3+wpRb1vzl2UYBog8FfGhttFKjtAclnZYrRo9ryG9uG/FZQU4AEg8ZE9LjGMzTmqKXPLnlWVnIlQQTvxJf8ip7VgjZjyVPrjw1te5otM7RmP7xm+sK2Gv9I8Gi++BRbEkR9EBw8zRUcKxwp73xkaLiqQb+kGduJTNHG72zcW9LoJgqQxpP3/Tj//c3yB0tqzaml05/+orHLksVO+95kX7/7qgJvnjlrfr2Ggsyx0eoy9uPzN5SPd86aXggOsEKW2Prz7du3VID3/tzs/sSRs2w7ovVHKtjrX2pd7ZMlTxAYfBAL9jiDwfLkq55Tm7ifhMlTGPyCAs7RFRhn47JnlcB9RM5T97ASuZXIcVNuUDIndpDbdsfrqsOppeXl5Y+XVKdjFCTh+zGaVuj0d9zy05PPK3QzBamxdwtTCrzyg/2Rvf2EstUjordGwa/kx9mSJLr8mLLtCW8HHGJc2R5hS219IiF6PnTusOqcMl57gm0Z8kanKMAQg0qSyuZfn7zItsbGyO9QlnxY0eCuD1XL2ys/MsrQhltE7Ug0uFOzufJFE2PxBo/YAx8XPPdDwWN0MrDRYIZF0mSMKCNHgaIVFoBbNoLJ7tEQDKxGF0kcLQimojCZopv0OkNOyWCCg9XMVAi7ARJzQdM2QUh0gmBozjc3Skg6dSBRqDGYSUOu66Zg+I2fNZs/M3/f/Grl/XnyF1Gw3VKCez0PN5IUfFLqvgUN4C0qNqYs5YhPL+aVZYDE4IpUk57oSFnJm4FyCqqOE0jhY2SMyLFoo56zyo6becOS5UVDdj7Vih0zp+tcMhwRpBeLyqtIjlJKAIZSbI8SGSF3k0pA3mR5tHuwPFoa7N7reoq2bqCsAk1HqCu5uvI1n6JuRXI+S1Mco54YmYTwcn6Aeic+kssXi8XpXC4V3t7/ADuTNKaQJdScAAAAAElFTkSuQmCC)](https://mybinder.org)
 at the top of each example page in the documentation.

Alternatively, you can run these scripts directly using the `.toml` files in this directory.
To do that, first open julia and activate this directory's project, e.g.,
```
$ julia --project="path/to/examples/directory"
```
Then instantiate the project in this directory, i.e.,
```julia
] instantiate
```
to install dependencies. Then run any of the examples by
```julia
include("path/to/examples/directory/example_script.jl")
```
# GeophysicalFlows.jl Documentation

## Overview

`GeophysicalFlows.jl` is a collection of modules which leverage the 
[FourierFlows.jl](https://github.com/FourierFlows/FourierFlows.jl) framework to provide
solvers for problems in Geophysical Fluid Dynamics, on periodic domains using Fourier-based pseudospectral methods.


## Examples

Examples aim to demonstrate the main functionalities of each module. Have a look at our Examples collection!


!!! note "Fourier transforms normalization"
    
    Fourier-based pseudospectral methods rely on Fourier expansions. Throughout the 
    documentation we denote symbols with hat, e.g., ``\hat{u}``, to be the Fourier transform 
    of ``u`` like, e.g.,
    
    ```math
    u(x) = \sum_{k_x} \hat{u}(k_x) \, e^{i k_x x} .
    ```
    
    The convention used in the modules is that the Fourier transform of a variable, e.g., `u` 
    is denoted with `uh` (where the trailing `h` is there to imply "hat"). Note, however, 
    that `uh` is obtained via a FFT of `u` and due to different normalization factors that the 
    FFT algorithm uses, `uh` _is not_ exactly the same as ``\hat{u}`` above. Instead,
    
    ```math
    \hat{u}(k_x) = \frac{𝚞𝚑}{n_x e^{- i k_x x_0}} ,
    ```
    
    where ``n_x`` is the total number of grid points in ``x`` and ``x_0`` is the left-most 
    point of our ``x``-grid.
    
    Read more in the FourierFlows.jl Documentation; see 
    [Grids](https://fourierflows.github.io/FourierFlowsDocumentation/stable/grids/) section.


!!! info "Unicode"
    Oftentimes unicode symbols are used in modules for certain variables or parameters. For 
    example, `ψ` is commonly used to denote the  streamfunction of the flow, or `∂` is used 
    to denote partial differentiation. Unicode symbols can be entered in the Julia REPL by 
    typing, e.g., `\psi` or `\partial` followed by the `tab` key.
    
    Read more about Unicode symbols in the 
    [Julia Documentation](https://docs.julialang.org/en/v1/manual/unicode-input/).


## Developers

The development of GeophysicalFlows.jl started by [Navid C. Constantinou](http://www.navidconstantinou.com) and [Gregory L. Wagner](https://glwagner.github.io) during the 21st AOFD Meeting 2017. During the 
course of time various people have contributed to GeophysicalFlows.jl, including 
[Lia Siegelman](https://scholar.google.com/citations?user=BQJtj6sAAAAJ), [Brodie Pearson](https://brodiepearson.github.io), and [André Palóczy](https://scholar.google.com/citations?user=o4tYEH8AAAAJ) (see the [example in FourierFlows.jl](https://fourierflows.github.io/FourierFlowsDocumentation/stable/literated/OneDShallowWaterGeostrophicAdjustment/)).


## Citing

If you use GeophysicalFlows.jl in research, teaching, or other activities, we would be grateful 
if you could mention GeophysicalFlows.jl and cite our paper in JOSS:

Constantinou et al., (2021). GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs. _Journal of Open Source Software_, **6(60)**, 3053, doi:[10.21105/joss.03053](https://doi.org/10.21105/joss.03053).

The bibtex entry for the paper is:

```bibtex
@article{GeophysicalFlowsJOSS,
  doi = {10.21105/joss.03053},
  url = {https://doi.org/10.21105/joss.03053},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {60},
  pages = {3053},
  author = {Navid C. Constantinou and Gregory LeClaire Wagner and Lia Siegelman and Brodie C. Pearson and André Palóczy},
  title = {GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs \& GPUs},
  journal = {Journal of Open Source Software}
}
```
# GPU

GPU-functionality is enabled via `FourierFlows.jl`. For more information on how `FourierFlows.jl`
handled with GPUs we urge you to the corresponding [`FourierFlows.jl` documentation section ](https://fourierflows.github.io/FourierFlowsDocumentation/stable/gpu/).

All `GeophysicalFlows.jl` modules can be run on GPU by providing `GPU()` as the device (`dev`) 
argument in the problem constructors. For example,

```julia
julia> GeophysicalFlows.TwoDNavierStokes.Problem(GPU())
Problem
  ├─────────── grid: grid (on GPU)
  ├───── parameters: params
  ├────── variables: vars
  ├─── state vector: sol
  ├─────── equation: eqn
  ├────────── clock: clock
  └──── timestepper: RK4TimeStepper
```

## Selecting GPU device

`FourierFlows.jl` can only utilize a single GPU. If your machine has more than one GPU available, 
then using functionality within `CUDA.jl` package enables you can choose the GPU device that 
`FourierFlows.jl` should use. The user is referred to the [`CUDA.jl` Documentation](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#Device-Management); in particular, [`CUDA.devices`](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#CUDA.devices) and [`CUDA.CuDevice`](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#CUDA.CuDevice).
# Stochastic Forcing

Forcing terms are implemented in various modules. Forcing can be either deterministic or 
stochastic (random). For deterministic forcing the implementation is straightforward; for 
stochastic forcing there are two main train of thoughts: Itô calculus and Stratonovich calculus.

Both stochastic calculi give the same results. But once we decide to use one of the two calculi 
we have to remain consistent and use that calculus throughout. There can be a lot of confusion 
and oftentimes confusion stems from mixing the two different stochastic calculi in a single 
computation instead of using one of the two consistently all along.

!!! note "Itô or Stratonovich in GeophysicalFlows.jl?"
    All modules included in GeophysicalFlows.jl use **Stratonovich calculus**.
		
The choice of Stratonovich calculus for GeophysicalFlows.jl was made since this calculus "works 
the same" with both stochastic and deterministic forcing, i.e. with Stratonovich calculus we 
have the same chain rules for differentiation for stochastic functions as the chain rules we 
learn in normal-deterministic calculus). Therefore, with the Stratonovich calculus the code 
does not really "care" whether the user implement deterministic or stochastic forcing.

If you are interested in learning more regarding the two stochastic calculi and  how they are 
numerically implemented then read on; otherwise you can skip this section of the documentation.

## Stochastic Differential Equations (SDEs)

A differential equation:

```math
	\frac{\mathrm{d} x}{\mathrm{d} t} = f(x) , \quad x(t_0) = 0,
```

can also be equivalently written in an integral form:

```math
	x(t) = \int_{t_0}^t f(x(s)) \, \mathrm{d} s.
```

In a similar manner, a stochastic differential equation (SDE),

```math
	\mathrm{d} x = f(x) \, \mathrm{d} t + g(x) \, \mathrm{d} W_t , \quad x(t_0) = 0 ,
```

with ``W_t`` a Brownian motion or Wiener process, can be written in an integral form as:

```math
	x(t) = \int_{t_0}^{t} f(x(s)) \, \mathrm{d} s + \int_{t_0}^{t} g(x(s)) \, \mathrm{d} W_s .
```

!!! tip "Wiener process"
    A Wiener process is a random variable ``W_t`` that depends continuously on ``t \ge 0`` and satisfies the following properties:
    1. Independence. For ``0 \le s \le t`` the increment ``W_t - W_s`` is independent of any prior values, i.e., independent of all ``W_\tau``, ``\tau \le s``.
    2. Stationarity. The statistical distribution of the increment ``W_{t+s} − W_s`` does not depend on  ``s`` (and so is identical in distribution to ``W_t``).
    3. Gaussianity. ``W_t`` is a Gaussian process with mean ``\langle W_t \rangle = 0`` and covariance ``\langle W_t W_s \rangle = \min(t, s)``.

!!! tip Notation, e.g., ``x_t``
    It's common to use notation ``x_t`` to denote explicit ``t``-dependence of variable ``x``. Not to be confused with the other common usage of subscripts for denoting partial differentiation.

The last integral in the integral representation of a SDE expression above is a stochastic integral 
(it involves a stochastic differential, `` \mathrm{d} W_t``). There is not a single straight-forward 
way for computing the value of a stochastic integral. The various ways we can approximate the 
value of a stochastic integral as a Riemannian sum each lead to a different answer. The two most 
popular ways for computing such stochastic integrals are:

```math
\begin{aligned}
{\color{Green} \text{Itô}} &: {\color{Green}\int_{t_0}^{t} g(x(s)) \, \mathrm{d} W_s \approx \sum_{j} g \left ( x(t_j) \right )(W_{j+1} - W_j)} , \\
{\color{Magenta} \text{Stratonovich}} &: {\color{Magenta} \int_{t_0}^{t} g(x(s)) \, \mathrm{d} W_s \approx \sum_{j} g \left (x \left (\tfrac{1}{2}(t_j + t_{j+1}) \right ) \right)(W_{j+1} - W_j)} .
\end{aligned}
```

The difference in the two calculi above lies in the point at which we choose to evaluate ``g(x)``:
we take the start of the time-interval for ``{\color{Green} \text{Itô, } t_j}``, while we use
the mid-point for ``{\color{Magenta}\text{Stratonovich, } \tfrac{1}{2}(t_j+t_{j+1})}``. In the case the 
stochastic noise is additive, i.e., its prefactor ``g`` does not depend on the state ``x_t``,
then the two interpretations coincide. When the noise does depend on the state of the system, 
i.e., ``g=g(x(t))``, then the two interpretations above give thoroughly different results. This
happens because the white noise process is not continuous and, therefore, the two interpretations
of the stochastic integrals above do not converge to the same result.

To overcome this apparent inconsistency, the two choices above come together with different 
chain rules, i.e., chain rules that are not necessarily the same as those in plain old calculus.
Let us see how different choices for computing the stochastic integrals bring about the need 
for different chain rules.

An SDE can be written also in differential form. Because we cannot formally form the derivative
``\mathrm{d} W / \mathrm{d} t``, since ``W`` is nowhere differentiable, we write an SDE in 
differential form as:

```math
\begin{aligned}
{\color{Green}\text{Itô}} &: {\color{Green}\mathrm{d} x_t = f(x_t) \mathrm{d} t + g(x_t) \mathrm{d} W_t} , \\
{\color{Magenta}\text{Stratonovich}} &: {\color{Magenta}\mathrm{d} x_t = f(x_t) \mathrm{d} t + g(x_t) \circ \mathrm{d} W_t} .
\end{aligned}
```

The circle in the term ``{\color{Magenta}g(x_t) \circ \mathrm{d} W_t}`` is used to differentiate 
between Itô and Stratonovich calculus.

Let's now assume we perform a variable change ``y = G(x)``. It turns out that according to 
which interpretation of the stochastic integral one chooses to use, then the following chain 
rule must be used:

```math
\begin{aligned}
{\color{Green}\text{Itô}} &: {\color{Green}\mathrm{d} y_t = \frac{\mathrm{d} G}{\mathrm{d} x} \mathrm{d} x_t + \frac{1}{2} g(x_t)^2 \frac{\mathrm{d}^2 G}{\mathrm{d} x^2} \mathrm{d} t = \left[ \frac{\mathrm{d} G}{\mathrm{d} x} f(x_t) + \frac{1}{2} g(x_t)^2 \frac{\mathrm{d}^2 G}{\mathrm{d} x^2} \right] \mathrm{d} t + \frac{\mathrm{d} G}{\mathrm{d} x} g(x_t) \mathrm{d} W_t} , \\
{\color{Magenta}\text{Stratonovich}} &: {\color{Magenta}\mathrm{d} y_t  = \frac{\mathrm{d} G}{\mathrm{d} x} \mathrm{d} x_t = \frac{\mathrm{d} G}{\mathrm{d} x} f(x_t) \mathrm{d} t + \frac{\mathrm{d} G}{\mathrm{d} x} g(x_t) \circ \mathrm{d} W_t} .
\end{aligned}
```

The above are the so-called *stochastic chain rules*. All derivatives of ``G`` are evaluated 
at ``x_t``. For Stratonovich calculus, the chain rule resembles the usual chain rule one learns
in calculus; for Itô there exists an additional term, often referred to as the "drift-term": 
``{\color{Green}\tfrac{1}{2} g^2 \, \mathrm{d}^2G / \mathrm{d} x^2}``.

It's easy to see that the extra drift-term in Itô's interpretation of the stochastic integral, 
is *exactly* equal to the ensemble mean  over forcing realizations of the Stratonovich 
stochastic integral. That's because the Itô stochastic integral has, by construction, zero 
ensemble mean since at every instant the noise is multiplied with ``g`` which is  evaluated 
at a time instance before the action of the noise; ``g`` and ``\mathrm{d} W`` are uncorrelated 
and thus:

```math
{\color{Green} \left \langle g(x_t) \mathrm{d} W_t \right \rangle = 0} \quad \text{while} \quad {\color{Magenta} \left \langle g(x_t) \circ \mathrm{d} W_t \right \rangle \ne 0} .
```

The above is demonstrated by evaluating the simple stochastic integral:

```math
\begin{aligned}
{\color{Green} \text{Itô}} &: {\color{Green} \left \langle \int_{t_0}^{t} W_s \, \mathrm{d} W_s \right \rangle \approx \sum_{j} \left \langle W_j (W_{j+1} - W_j) \right \rangle} \\
& \hspace{7.3em} {\color{Green} = \sum_j \left \langle W_j W_{j+1} \right \rangle - \left \langle W_j W_j \right \rangle \sim \sum_{j} t_j - t_j = 0} , \\
{\color{Magenta}\text{Stratonovich}} &: {\color{Magenta}\left \langle \int_{t_0}^{t} W_s \circ \mathrm{d} W_s \right \rangle \approx \sum_j \left \langle \frac1{2}(W_j + W_{j+1}) (W_{j+1} - W_j)\right \rangle} \\
& \hspace{7.3em} {\color{Magenta} = \frac1{2} \sum_j \left \langle W_{j+1} W_{j+1} \right \rangle - \left \langle W_j W_j \right \rangle  \sim \frac1{2} \sum_j t_{j+1} - t_j = \frac{t}{2}} .
\end{aligned}
```

SDEs rarely can be solved in closed form; most often numerical solution of SDEs is brought to 
the rescue. Itô calculus has the advantage that is very easily implemented numerically. On 
the other hand, the chain rule in Stratonovich calculus coincides with that in normal calculus. 
This stems from the fact that in the Stratonovich interpretation the white noise process is as 
a series of colored noise processes with the de-correlation time tending to zero. This made 
Stratonovich calculus more popular in the physics community. A nice discussion on the differences 
and similarities between the two calculi is given by [van Kampen](https://doi.org/10.1007/BF01007642).

## A simple Stochastic Differential Equation: the Ornstein--Uhlenbeck process

One of the simplest SDEs is the Ornstein--Uhlenbeck process, a variation of which is:

```math
x(t) = - \int_{t_0}^{t} \mu x(s) \, \mathrm{d} s + \int_{t_0}^{t} \sqrt{\sigma} \, \mathrm{d} W_s . \tag{1}
```

Note that in differential form (1) is written as:

```math
\mathrm{d} x_t = - \mu x_t \, \mathrm{d} t + \sqrt{\sigma} \, \mathrm{d} W_t . \tag{2}
```

Luckily, for (2) we don't need to distinguish between Itô and Stratonovich, since ``g`` is 
independent of ``x(t)``. But note that oftentimes this is not the case; that ``g`` is independent 
of ``x(t)`` is only a fortuitous coincident for this particular SDE.

How do we time-step SDE (2) numerically? Let us assume a discretization of time into time-steps
of duration ``\tau``, i.e., ``t_j = (j-1) \tau``, ``j=1, 2, \dots``. (What follows is easily 
generalized to non-uniform time discretizations.) With that in mind, we denote ``x_j \equiv x(t_j)``. 
Then the Euler--Mayorama time-stepping scheme for (2) is

```math
	x_{j+1} = x_j + (-\mu x_j) \tau + \sqrt{\sigma} (W_{j+1} - W_j) .
```

Now let us ask the following question: How can we compute the work done by the noise?
In other words, if we are interested in the evolution of the "energy", defined as
``E \equiv \tfrac{1}{2} x^2``, then how does the noise term attribute in the growth of ``E``? 
To answer that we first have to find the SDE that energy ``E`` obeys. But, in doing so, it 
is important to adopt a single interpretation for computing stochastic integrals as now a 
transformation of variables is needed. That is, depending on whether we choose to interpret 
the stochastic integrals according to Itô or to Stratonovich calculus, ``E`` evolves according
to:

```math
\hspace{3.35em} {\color{Green} \text{Itô}} : {\color{Green} \mathrm{d} E_t = \left ( -2 \mu E_t + \tfrac{1}{2} \sigma \right ) \mathrm{d} t + x_t \sqrt{\sigma} \mathrm{d} W_t} , \tag{3}
```
```math
\hspace{-3.35em} {\color{Magenta} \text{Stratonovich}} : {\color{Magenta} \mathrm{d} E_t = -2 \mu E_t \mathrm{d} t + x_t \circ \sqrt{\sigma} \mathrm{d} W_t} . \tag{4}
```

The term ``-2 \mu E_t`` in both cases is the dissipation of energy by the ``\mu`` term; the 
rest of the terms involve the noise. How do we compute the work ``P`` done by the noise? 
Well, it follows that:

```math
\begin{aligned}
{\color{Green} \text{Itô}} &: {\color{Green} P_t = \tfrac{1}{2} \sigma \mathrm{d} t + \sqrt{\sigma} x_t \mathrm{d} W_t \approx \tfrac{1}{2} \sigma \, \mathrm{d}t + \sqrt{\sigma} x_j (W_{j+1} - W_j)} , \\
{\color{Magenta} \text{Stratonovich}} &: {\color{Magenta} P_t =  x_t \circ \sqrt{\sigma} \mathrm{d} W_t \approx \sqrt{\sigma} x \left ( \tfrac{1}{2} (t_j + t_{j+1}) \right ) (W_{j+1} - W_j)} .
\end{aligned}
```

Now let's assume for a moment that we didn't know the rules for transforming Stratonovich to 
Itô and we were wondering what is the extra drift term we have to include in the Itô formulations, 
i.e., the ``\tfrac{1}{2} \sigma`` term. We can compute the Itô's drift-term using the fact that 
it is exactly equal to ``\langle x_t \circ \sqrt{\sigma} \mathrm{d} W_t \rangle``; and for the 
latter we can use the "usual" calculus. That is, we rewrite (1) as:

```math
\dot{x} = -\mu x + \xi , \tag{5}
```

where ``\xi(t)`` is understood to be the "continuous" version of the white-noise process (which 
is formally only understood in terms of distributions). The forcing ``\xi`` has the properties:

```math
\left \langle \xi(t) \right \rangle = 0 \quad \text{and} \quad \left \langle \xi(t) \xi(t') \right \rangle = \sigma \delta(t - t') .
```

Thus we need to compute ``\langle P_t \rangle = \langle x(t) \xi(t) \rangle``. But (5) formally
has the solution:

```math
x(t) = e^{-\mu t} x(0) + \int_0^t e^{-\mu (t - s)} \xi(s) \, \mathrm{d} s .
```

and using this solution we get

```math
\langle P_t \rangle = \langle x(t) \xi(t) \rangle =  e^{-\mu t} \underbrace{\langle x(0) \xi(t) \rangle}_{=0} + \int_0^t e^{-\mu (t - s)} \langle \xi(t) \xi(s) \rangle \, \mathrm{d} s = \sigma \int_0^t e^{- \mu (t - s)} \delta(t - s) \, \mathrm{d} s = \frac{\sigma}{2} .
```

Above we used that ``\int_0^t \delta(t - s) \mathrm{d} s = \tfrac{1}{2}``, which is consistent 
with the Stratonovich symmetric interpretation of stochastic integrals.

### Numerical implementation

How do we time-step the equation for ``E``? In the case of Itô's interpretation, (3), we use 
the Euler--Maruyama time-stepping scheme:

```math
	E_{j+1} = E_j + \left ( -2 \mu E_j + \frac{\sigma}{2} \right ) \tau + \sqrt{\sigma} x_j (W_{j+1} - W_j).
```
However, we cannot use Euler--Maruyama for time-stepping the corresponding Stratonovich 
version of (4), since the Euler--Maruyama scheme involves "Itô"-thinking. To time-step (4) we 
have to approximate ``g`` in the middle of the time-step. There are many ways to do that, one 
of which is the, so called, Euler--Heun method:

```math
\begin{aligned}
\widetilde{E}_{j+1} &= E_j + (- 2\mu E_j) \tau + \sqrt{\sigma} x_j (W_{j+1} - W_j), \\
E_{j+1} &= E_j + \left( -2 \mu \frac{E_j + \widetilde{E}_{j + 1}}{2} \right)\tau + \sqrt{\sigma}\frac{x_j + x_{j+1}}{2} (W_{j+1} - W_j) .
\end{aligned}
```

Let's apply not Euler--Maruyama and Euler--Heun schemes to time-step (3) and (4) respectively
and compare the results with those obtained from time-stepping (2) and computing ``E`` a 
posteriori. 

Figure below compares the energy evolution as predicted by:
- direct computation from the ``x_t`` time-series: ``\tfrac{1}{2} x_t^2``,
- time-integration of (3) using Euler--Maruyama, and
- time-integration of (4) using Euler--Heun.


```@setup 1
using Plots
Plots.default(lw=2)
```

```@example 1
using Plots
using Statistics: mean
using Random: randn, seed!
seed!(1234) # for reproducing the same plots

                μ = 0.2
                σ = 0.2    # noise strength
               dt = 0.01   # timestep
           nsteps = 2001   # total timesteps
   n_realizations = 1000   # how many forcing realizations
some_realizations = 20     # used for plotting to illustrate convergence

t = 0:dt:(nsteps-1)*dt 	# time

ΔW = randn(nsteps, n_realizations) * sqrt(dt) # noise

# Numerical calculation
x = zeros(size(ΔW))
E_ito = zeros(size(ΔW))
E_str = zeros(size(ΔW))
E_numerical = zeros(size(ΔW))

for j = 2:nsteps # time step the equations
	
  # time-step dx = - μ x dt + √σ dW
  @. x[j, :] = x[j-1, :] - μ * x[j-1, :] * dt + sqrt(σ) * ΔW[j-1, :]

  # time-step dE = (- 2μ E + ½σ) dt + √σ x dW
  @. E_ito[j, :] = E_ito[j-1, :] + (-2μ * E_ito[j-1, :]
	                   + σ/2) * dt + sqrt(σ) * x[j-1, :] * ΔW[j-1, :]

  # time-step dE = - 2μ E dt + √σ x ∘ dW
  xbar = @. x[j-1, :] - μ * x[j-1, :] * dt + sqrt(σ) * ΔW[j-1, :]
  Ebar = @. E_str[j-1, :] + (-2μ * E_str[j-1, :]) * dt + sqrt(σ) * x[j-1, :] * ΔW[j-1, :]
  @. E_str[j, :] = E_str[j-1, :] + (-2μ * (E_str[j-1, :]
		+ Ebar) / 2) * dt + sqrt(σ) * (x[j-1, :] + xbar) / 2 * ΔW[j-1, :]
end

# direct computation of E from x
@. E_numerical = 0.5 * x^2

# compare the three E(t) solutions
plot(μ * t, [E_numerical[:, 1] E_ito[:, 1] E_str[:, 1]],
          linewidth = [3 2 1],
              label = ["½ xₜ²" "Eₜ (Ito)" "Eₜ (Stratonovich)"],
          linestyle = [:solid :dash :dashdot],
             xlabel = "μ t",
             ylabel = "E",
             legend = :topleft,
              title = "comparison of E(t) for single realization")

savefig("assets/energy_comparison.svg"); nothing # hide
```

![energy_comparison](assets/energy_comparison.svg)

Now we can further compute the "energy" budgets, i.e., the work done by the noise versus the
energy loss by the ``μ`` term, using Itô and Stratonovich formalisms. Figures below show 
the ensemble mean energy budgets (using 1000 ensemble members) as computed using Itô and
Stratonovich calculus. For the energy budget to close we have to be consistent: if we time-step 
the  energy equation based on Stratonovich calculus then we must compute the work also according 
to Stratonovich and vice versa.

```@example 1
# theoretical results for ⟨E⟩ and d⟨E⟩/dt
   E_theory = @. σ/4μ * (1 - exp(-2μ * t))
dEdt_theory = @. σ/2  * exp(-2μ * t)

# compute d⟨E⟩/dt numerically
dEdt_ito = mean(E_ito[2:nsteps, :] .- E_ito[1:nsteps-1, :], dims=2) / dt

# compute the work and dissipation
work_ito = mean(sqrt(σ) * ΔW[1:nsteps-1, :] / dt .* x[1:nsteps-1, :] .+ σ/2, dims=2)
diss_ito = 2*μ * (mean(E_ito[1:nsteps-1, :], dims=2))

# Ensemble mean energy budgets from the Itô integration

plot_E = plot(μ * t, [E_theory mean(E_ito[:, 1:some_realizations], dims=2) mean(E_ito, dims=2)],
                linewidth = [3 2],
	                  label = ["theoretical ⟨E⟩" "⟨E⟩ from $some_realizations ensemble members" "⟨E⟩ from $n_realizations ensemble members"],
	                 xlabel = "μ t",
	                 ylabel = "E",
	                 legend = :bottomright,
	                  title = "Ito: 𝖽Eₜ = (-2μ Eₜ + ½σ) 𝖽t + xₜ √σ 𝖽Wₜ")

plot_Ebudget = plot(μ * t[1:nsteps-1], [dEdt_ito work_ito.-diss_ito dEdt_theory[1:nsteps-1]],
                linestyle = [:dash :dashdot :solid],
                linewidth = [2 1 3],
                    label = ["numerical 𝖽⟨E⟩/𝖽t" "⟨work - dissipation⟩" "theoretical 𝖽⟨E⟩/𝖽t"],
                   legend = :topright,
                   xlabel = "μ t")

plot(plot_E, plot_Ebudget, layout=grid(2, 1, heights=[0.65 ,0.35]), size=(600, 525))

savefig("assets/energy_budgets_Ito.svg"); nothing # hide
```

![energy_budgets_Ito](assets/energy_budgets_Ito.svg)


```@example 1
# compute d⟨E⟩/dt numerically
dEdt_str = mean(E_str[2:nsteps, :] .- E_str[1:nsteps-1, :], dims=2) / dt

# compute the work and dissipation
work_str = mean(sqrt(σ) * ΔW[1:nsteps-1, :] / dt .* (x[1:nsteps-1, :] .+ x[2:nsteps, :])/2, dims=2)
diss_str = 2*μ * (mean(E_str[1:nsteps-1, :], dims=2))

plot_E = plot(μ * t, [E_theory mean(E_str[:, 1:some_realizations], dims=2) mean(E_str, dims=2)],
                linewidth = [3 2],
                    label = ["theoretical ⟨E⟩" "⟨E⟩ from $some_realizations ensemble members" "⟨E⟩ from $n_realizations ensemble members"],
                   xlabel = "μ t",
                   ylabel = "E",
                   legend = :bottomright,
                    title = "Stratonovich: 𝖽Eₜ = -2μ Eₜ 𝖽t + xₜ ∘ √σ 𝖽Wₜ")

plot_Ebudget = plot(μ * t[1:nsteps-1], [dEdt_str[1:nsteps-1] work_str[1:nsteps-1].-diss_str[1:nsteps-1] dEdt_theory[1:nsteps-1]],
                linestyle = [:dash :dashdot :solid],
                linewidth = [2 1 3],
                    label = ["numerical 𝖽⟨E⟩/𝖽t" "⟨work - dissipation⟩" "theoretical 𝖽⟨E⟩/𝖽t"],
                   legend = :bottomleft,
                   xlabel = "μ t")

plot(plot_E, plot_Ebudget, layout=grid(2, 1, heights=[0.65 ,0.35]), size=(600, 525))

savefig("assets/energy_budgets_Stratonovich.svg"); nothing # hide
```

![energy_budgets_Stratonovich](assets/energy_budgets_Stratonovich.svg)


## A simple Stochastic Partial Differential Equation (SPDE)

We would like now to transfer all the knowledge we got from the previous sections to PDEs. 
In particular we'll start by focussing on the simple SPDE:

```math
\partial_t \nabla^2 \psi(\bm{x}, t) =  - \mu \nabla^2 \psi(\bm{x}, t) + \xi(\bm{x}, t) , \tag{6}
```

with periodic boundary conditions in both ``x`` and ``y``. SPDE (6) is also equivalently written as:

```math
\mathrm{d} \nabla^2 \psi_{t}(\bm{x}) = - \mu \nabla^2 \psi_{t} (\bm{x}) \mathrm{d} t + \sqrt{\sigma} \mathrm{d} W_{t} (\bm{x}) .
```

The form (6) is the continuous version, similar to (2). In this SPDE, since the forcing is 
purely additive, i.e., it does not depend on the state of the system, both Itô and Stratonovich 
interpretations coincide.

The forcing ``\xi`` obeys:

```math
\langle \xi(\bm{x}, t) \rangle = 0 \quad \text{and} \quad \langle \xi(\bm{x}, t) \xi(\bm{x}', t') \rangle = Q(\bm{x} - \bm{x}') \delta(t - t') ,
```

that is, the forcing is white in time but spatially correlated; its spatial correlation is 
prescribed by the function ``Q`` which is, necessarily, homogeneous in all its arguments
(see discussion by [Constantinou (2015)](http://arxiv.org/abs/1503.07644); Appendix A).

Equation (6) above describes the vorticity evolution of a two-dimensional fluid ``\nabla^2 \psi`` 
that is stochastically forced while dissipated by linear drag ``\mu``. The energy of the 
fluid is:

```math
E = \tfrac{1}{2} \overline{|\bm{\nabla} \psi|^2}^{x, y} = -\tfrac{1}{2} \overline{\psi \nabla^2 \psi}^{x, y} ,
```

where the overbar denotes average over ``x`` and ``y`` and an integration-by-parts was carried
through in the last equality. To obtain the energy equation we multiply (6) with ``-\psi`` and 
average over the whole domain. Thus, the work done by the forcing is given by:

```math
P = - \, \overline{\psi \, \xi}^{x, y} ,
```

but the above is a stochastic integral and it is meaningless without a rule for computing the stochastic integral.

Numerically, the work done by the forcing at the ``j``-th timestep can be obtained 
Stratonovich-wise via:

```math
\begin{aligned}
P_j = - \, \overline{\frac{\psi(\bm{x}, t_j) + \psi(\bm{x}, t_{j+1})}{2}  \xi(\bm{x}, t_{j+1}) }^{x,y} ,
\end{aligned}
```

or Itô-wise as

```math
\begin{aligned}
P_j = -\, \overline{ \psi(\bm{x}, t_j) \xi(\bm{x}, t_{j+1}) }^{x,y} + \text{drift} .
\end{aligned}
```

But how much is the Itô drift term in this case? As in the previous section, the drift is 
*precisely* the ensemble mean of the Stratonovich work, i.e.:

```math
\textrm{Ito drift}= - \overline{\langle \underbrace{\psi(\bm{x}, t) \circ  \xi(\bm{x}, t)}_{\textrm{Stratonovich}} \rangle}^{x, y} .
```

But again, the above can be computed using the "formal" solution of (6):

```math
\psi(\bm{x}, t) = e^{-\mu t} \psi(\bm{x}, 0) + \int_0^t e^{- \mu (t - s)} \nabla^{-2} \xi(\bm{x}, s) \, \mathrm{d} s ,
```

which implies

```math
\begin{aligned}
\text{drift} & = -\overline{e^{- \mu t} \underbrace{\left \langle \psi(\bm{x}, 0) \xi(\bm{x}, t) \right \rangle}_{=0}}^{x, y} - \int_0^t e^{- \mu (t - s)} \overline{\nabla^{-2} \left \langle \xi(\bm{x}, s) \xi(\bm{x}, t) \right\rangle}^{x, y} \, \mathrm{d} s \\
& = - \int_0^t e^{-\mu(t - s)} \overline{\underbrace{\left [ \nabla^{-2} Q (\bm{x}) \right ] \big|_{\bm{x}=0}}_{\text{independent of }x, y} \, \delta(t - s)}^{x,y} \, \mathrm{d} s \\
& = - \frac1{2} \nabla^{-2} Q(\bm{x}) \big|_{\bm{x}=0} \\
& = - \frac1{2} \left [ \nabla^{-2} \int \frac{\mathrm{d}^2 \bm{k}}{(2\pi)^2} \widehat{Q}(\bm{k}) \, e^{i \bm{k} \bm{\cdot} \bm{x}} \right ]_{\bm{x}=0} \\
& = \int \frac{\mathrm{d}^2 \bm{k}}{(2\pi)^2} \frac{\widehat{Q}(\bm{k})}{2 |\bm{k}|^2} .
\end{aligned}
```

Thus, the drift, or in this case the mean energy input rate by the stochastic forcing, is 
precisely determined from the spatial correlation of the forcing, ``Q``. Let us denote the 
drift as:

```math
\varepsilon \equiv \int \frac{\mathrm{d}^2 \bm{k}}{(2\pi)^2} \frac{\widehat{Q}(\bm{k})}{2 |\bm{k}|^2} . \tag{7}
```

Using the above, the work for a single forcing realization at the ``j``-th timestep is numerically 
computed as:

```math
{\color{Green} \text{Itô}} : {\color{Green} P_j = -\overline{\psi(\bm{x}, t_j) \xi(\bm{x}, t_{j+1})}^{x, y} + \varepsilon} , \tag{8}
```
```math
{\color{Magenta} \text{Stratonovich}} : {\color{Magenta} P_j = -\overline{\frac{\psi(\bm{x}, t_j) + \psi(\bm{x}, t_{j+1})}{2} \xi(\bm{x}, t_{j+1})}^{x, y}} . \tag{9}
```

Remember, previously the work done by the stochastic forcing was:
```math
\mathrm{d} P_t = {\color{Green} \frac{\sigma}{2}\mathrm{d} t + \sqrt{\sigma} x_t \mathrm{d} W_t} = {\color{Magenta} \sqrt{\sigma} x_t \circ \mathrm{d} W_t} ,
```
and by sampling over various forcing realizations:
```math
\langle \mathrm{d} P_t \rangle = \frac{\sigma}{2} \mathrm{d} t = \langle \sqrt{\sigma} x_t \circ \mathrm{d} W_t \rangle .
```

All modules in GeophysicalFlows.jl use Stratonovich calculus. For example, the domain-averaged 
energy injected per unit time by the forcing in the `TwoDNavierStokes` module is computed 
using (9) via the [`energy_work`](@ref GeophysicalFlows.TwoDNavierStokes.energy_work) function.

## A bit more elaborate SPDE

It turns out everything carries through if in our SPDE above for the 2D vorticity equation we 
also include the nonlinear advection terms:

```math
\partial_t \nabla^2 \psi(\bm{x}, t) + \mathsf{J}(\psi, \nabla^2 \psi) = - \mu \nabla^2 \psi(\bm{x}, t) + \xi(\bm{x}, t) . \tag{10}
```

The nonlinearity does not alter the Itô drift; thus the ensemble mean energy input by the 
stochastic forcing, remains the same. We can easily verify this from the "formal" solution 
of (10):

```math
\psi(\bm{x}, t) = e^{- \mu t} \psi(\bm{x}, 0) + \int_0^t e^{- \mu (t - s)} \nabla^{-2} \xi(\bm{x}, s) \, \mathrm{d} s - \int_0^t \nabla^{-2} \mathsf{J} \left ( \psi(\bm{x}, s), \nabla^2 \psi(\bm{x}, s) \right ) \mathrm{d} s .
```

When multiplied with ``\xi(\bm{x}, t)`` the last term vanishes since its only non-zero 
contribution comes from the point ``s = t``, which is of measure zero (in the integrated sense). 

A demonstration of how the energy budgets can be computed when we have stochastic forcing is 
illustrated in an [example of the TwoDNavierStokes](../literated/twodnavierstokes_stochasticforcing_budgets/) 
module.
# Contributors Guide

This is a short guide for potential GeophysicalFlows.jl contributors.

Please feel free to ask us questions and chat, either by raising an [issue](https://github.com/FourierFlows/GeophysicalFlows.jl/issues) or starting a [discussion](https://github.com/FourierFlows/GeophysicalFlows.jl/discussions).

We follow the [ColPrac guide](https://github.com/SciML/ColPrac) for collaborative practices. 
New contributors should make sure to read that guide.

## What can I do?

* Tackle an existing [issue](https://github.com/FourierFlows/GeophysicalFlows.jl/issues).

* Try to run your favorite GeophysicalFlows.jl module and play around with it to simulate 
  your favorite setup. If you run into any problems or find it difficult to use, modify, or 
  understand, please [open an issue](https://github.com/FourierFlows/GeophysicalFlows.jl/issues)!

* Write up an example or tutorial on how to do something useful with one of the current modules
  in GeophysicalFlows.jl, like how to set up a new physical configuration.

* Improve documentation, docstrings, or comments if you found something is hard to use.

* Implement a new feature (e.g., a new diagnostic into a module).

* Implement a new module from scratch to solve your favorite partial differential equation with
  periodic boundary conditions.

If you're interested in working on something, let us know by commenting on an existing issue 
or by opening a new issue. This is to make sure no one else is working on the same issue and 
so we can help and guide you in case there is anything you need to know beforehand.

## Ground Rules

* Each pull request should consist of a logical collection of changes. You can
  include multiple bug fixes in a single pull request, but they should be related.
  For unrelated changes, please submit multiple pull requests.
* Do not commit changes to files that are irrelevant to your feature or bugfix
  (e.g., `.gitignore`).
* Be willing to accept criticism and work on improving your code; we don't want
  to break other users' code, so care must be taken not to introduce bugs. We
  discuss pull requests and keep working on them until we believe we've done a
  good job.
* Be aware that the pull request review process is not immediate, and is
  generally proportional to the size of the pull request.

## Reporting a bug

The easiest way to get involved is to report issues you encounter when using GeophysicalFlows.jl 
or by requesting something you think is missing.

* Head over to the [issues](https://github.com/FourierFlows/GeophysicalFlows.jl/issues) page.
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
  [GeophysicalFlows.jl repository](https://github.com/FourierFlows/GeophysicalFlows.jl) by
  clicking the "Fork" button.
* Clone your fork of the GeophysicalFlows.jl repository (in terminal on Mac/Linux or git shell/
  GUI on Windows) in the location you'd like to keep it.
  ```
  git clone https://github.com/your-user-name/GeophysicalFlows.jl.git
  ```
* Navigate to that folder in the terminal or in Anaconda Prompt if you're on Windows.
* Connect your repository to the upstream (main project).
  ```
  git remote add geophysicalflows https://github.com/FourierFlows/GeophysicalFlows.jl.git
  ```
* Create the development environment by opening Julia via `julia --project` then
  typing in `] instantiate`. This will install all the dependencies in the `Project.toml`
  file.
* You can test to make sure GeophysicalFlows.jl works by typing in `] test` which will run all
  the tests (this can take a while). In an ideal world you should run the tests on a machine
  with a GPU capability but if that's not a possibility that is available to you then don't 
  worry -- simply comment in a PR that you didn't test on GPU.

Your development environment is now ready!

## Pull Requests

Changes and contributions should be made via GitHub pull requests against the `main` branch.

When you're done making changes, commit the changes you made. Chris Beams has written 
a [guide](https://chris.beams.io/posts/git-commit/) on how to write good commit messages.

When you think your changes are ready to be merged into the main repository,
push to your fork and [submit a pull request](https://github.com/FourierFlows/GeophysicalFlows.jl/compare/).

**Working on your first Pull Request?** You can learn how from this _free_ video series
[How to Contribute to an Open Source Project on GitHub](https://egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github), Aaron Meurer's [tutorial on the git workflow](https://www.asmeurer.com/git-workflow/), 
or the guide [“How to Contribute to Open Source"](https://opensource.guide/how-to-contribute/).

## Documentation

All PRs that introduce new features or new modules should be accompanied with appropriate 
docstrings and documentation. Writing documentation strings is really important to make sure 
others use your functionality properly. Didn't write new functions? That's fine, but be sure 
that the documentation for the code you touched is still in great shape. It is not uncommon 
to find some strange wording or clarification that you can take care of while you are here.

We encourage using [unicode](https://docs.julialang.org/en/v1/manual/unicode-input/) characters 
when writing docstrings, e.g., use `α` instead of `\alpha`. This makes the rendering of the 
docstrings in the Documentation and in the Julia REPL's `help?>` mode as similar as possible.

You can preview how the Documentation will look like after merging by building the documentation 
locally. To do that, from the main directory of your local repository call

```
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs/ docs/make.jl
```
 
and then open `docs/build/index.html` in your favorite browser.

## Credits

This contributor's guide is heavily based on the [MetPy contributor's guide](https://github.com/Unidata/MetPy/blob/master/CONTRIBUTING.md) 
and on its "cover" made by [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/contributing/).
# Installation instructions

You can install the latest version of GeophysicalFlows.jl via the built-in package manager 
(accessed by pressing `]` in the Julia REPL command prompt) to add the package and also to 
instantiate/build all the required dependencies

```julia
julia>]
(v1.6) pkg> add GeophysicalFlows
(v1.6) pkg> instantiate
```

We recommend installing GeophysicalFlows.jl with the built-in Julia package manager, because 
this installs a stable, tagged release. Later on, you can update GeophysicalFlows.jl to the 
latest tagged release again via the package manager by typing

```julia
(v1.6) pkg> update GeophysicalFlows
```

Note that some releases might induce breaking changes to certain modules. If after anything 
happens or your code stops working, please open an [issue](https://github.com/FourierFlows/GeophysicalFlows.jl/issues) 
or start a [discussion](https://github.com/FourierFlows/GeophysicalFlows.jl/discussions). We're 
more than happy to help with getting your simulations up and running.

!!! warn "Use Julia 1.6 or newer"
    The latest GeophysicalFlows.jl requires at least Julia v1.6 to run.
    Installing GeophysicalFlows with an older version of Julia will install an older version 
    of GeophysicalFlows.jl (the latest version compatible with your version of Julia).
    
    Last version compatible with Julia v1.5: GeophysicalFlows.jl v0.13.1

    Last version compatible with Julia v1.0.5 (the current long-term-release): GeophysicalFlows.jl v0.5.1


# Aliasing


In pseudospectral methods, computing nonlinear terms results in aliasing errors. (Read more about
aliasing errors in the [FourierFlows.jl Documentation](https://fourierflows.github.io/FourierFlowsDocumentation/stable/aliasing/).) To avoid aliasing errors, we need to apply some dealiasing to our fields 
in Fourier space before transforming to physical space to compute nonlinear terms.

!!! info "De-aliasing scheme"
    FourierFlows.jl curently implements dealiasing by zeroing out the highest-`aliased_fraction` 
    wavenumber components on a `grid`. By default in FourierFlows.jl, `aliased_fraction=1/3`.
    Users can construct a `grid` with different `aliased_fraction` via
    
    ```julia
    julia> grid = OneDGrid(64, 2π, aliased_fraction=1/2)
    
    julia> OneDimensionalGrid
             ├─────────── Device: CPU
             ├──────── FloatType: Float64
             ├────────── size Lx: 6.283185307179586
             ├──── resolution nx: 64
             ├── grid spacing dx: 0.09817477042468103
             ├─────────── domain: x ∈ [-3.141592653589793, 3.0434178831651124]
             └─ aliased fraction: 0.5
    ```
    or provide the keyword argument `aliased_fraction` to the `Problem()` constructor of each
    module, e.g.,
    
    ```julia
    julia> prob = GeophysicalFlows.TwoDNavierStokes.Problem(; aliased_fraction=1/2)
    Problem
      ├─────────── grid: grid (on CPU)
      ├───── parameters: params
      ├────── variables: vars
      ├─── state vector: sol
      ├─────── equation: eqn
      ├────────── clock: clock
      └──── timestepper: RK4TimeStepper
      
    julia> prob.grid.aliased_fraction
    0.5
    ```

Currently, all nonlinearities in all modules included in GeophysicalFlows.jl modules are quadratic 
nonlinearities. Therefore, the default `aliased_fraction` of 1/3 is appropriate.

All modules apply de-aliasing by calling, e.g., `dealias!(prob.sol, prob.grid)` both before
computing any nonlinear terms and also during updating all variable, i.e., within `updatevars!`.

To disable de-aliasing you need to create a problem with a grid that has been constructed with 
the keyword `aliased_fraction=0`.
# Functions

## `GeophysicalFlows`

### Exported functions

```@docs
GeophysicalFlows.lambdipole
GeophysicalFlows.peakedisotropicspectrum
```


## `TwoDNavierStokes`

### Exported functions

```@docs
GeophysicalFlows.TwoDNavierStokes.Problem
GeophysicalFlows.TwoDNavierStokes.energy_dissipation_hyperviscosity
GeophysicalFlows.TwoDNavierStokes.energy_dissipation_hypoviscosity
GeophysicalFlows.TwoDNavierStokes.energy_work
GeophysicalFlows.TwoDNavierStokes.enstrophy_dissipation_hyperviscosity
GeophysicalFlows.TwoDNavierStokes.enstrophy_dissipation_hypoviscosity
GeophysicalFlows.TwoDNavierStokes.enstrophy_work
```

### Private functions

```@docs
GeophysicalFlows.TwoDNavierStokes.calcN_advection!
GeophysicalFlows.TwoDNavierStokes.addforcing!
GeophysicalFlows.TwoDNavierStokes.energy_dissipation
GeophysicalFlows.TwoDNavierStokes.enstrophy_dissipation
```


## `SingleLayerQG`

### Exported functions

```@docs
GeophysicalFlows.SingleLayerQG.Problem
GeophysicalFlows.SingleLayerQG.streamfunctionfrompv!
GeophysicalFlows.SingleLayerQG.energy_dissipation
GeophysicalFlows.SingleLayerQG.energy_work
GeophysicalFlows.SingleLayerQG.energy_drag
GeophysicalFlows.SingleLayerQG.enstrophy
GeophysicalFlows.SingleLayerQG.enstrophy_dissipation
GeophysicalFlows.SingleLayerQG.enstrophy_work
GeophysicalFlows.SingleLayerQG.enstrophy_drag
```

### Private functions

```@docs
GeophysicalFlows.SingleLayerQG.calcN_advection!
GeophysicalFlows.SingleLayerQG.addforcing!
```


## `MultiLayerQG`

### Exported functions

```@docs
GeophysicalFlows.MultiLayerQG.Problem
GeophysicalFlows.MultiLayerQG.fwdtransform!
GeophysicalFlows.MultiLayerQG.invtransform!
GeophysicalFlows.MultiLayerQG.streamfunctionfrompv!
GeophysicalFlows.MultiLayerQG.pvfromstreamfunction!
```

### Private functions

```@docs
GeophysicalFlows.MultiLayerQG.LinearEquation
GeophysicalFlows.MultiLayerQG.calcNlinear!
GeophysicalFlows.MultiLayerQG.calcN_advection!
GeophysicalFlows.MultiLayerQG.calcN_linearadvection!
GeophysicalFlows.MultiLayerQG.addforcing!
```


## `SurfaceQG`

### Exported functions

```@docs
GeophysicalFlows.SurfaceQG.Problem
GeophysicalFlows.SurfaceQG.buoyancy_dissipation
GeophysicalFlows.SurfaceQG.buoyancy_work
```

### Private functions

```@docs
GeophysicalFlows.SurfaceQG.calcN_advection!
GeophysicalFlows.SurfaceQG.addforcing!
```


## `BarotropicQGQL`

### Exported functions

```@docs
GeophysicalFlows.BarotropicQGQL.Problem
GeophysicalFlows.BarotropicQGQL.dissipation
GeophysicalFlows.BarotropicQGQL.work
GeophysicalFlows.BarotropicQGQL.drag
```

### Private functions

```@docs
GeophysicalFlows.BarotropicQGQL.calcN_advection!
GeophysicalFlows.BarotropicQGQL.addforcing!
```
# Private types

## `TwoDNavierStokes`

```@docs
GeophysicalFlows.TwoDNavierStokes.Params
GeophysicalFlows.TwoDNavierStokes.Vars
GeophysicalFlows.TwoDNavierStokes.DecayingVars
GeophysicalFlows.TwoDNavierStokes.ForcedVars
GeophysicalFlows.TwoDNavierStokes.StochasticForcedVars
```


## `SingleLayerQG`

```@docs
GeophysicalFlows.SingleLayerQG.Params
GeophysicalFlows.SingleLayerQG.BarotropicQGParams
GeophysicalFlows.SingleLayerQG.EquivalentBarotropicQGParams
GeophysicalFlows.SingleLayerQG.Vars
GeophysicalFlows.SingleLayerQG.DecayingVars
GeophysicalFlows.SingleLayerQG.ForcedVars
GeophysicalFlows.SingleLayerQG.StochasticForcedVars
```


## `MultiLayerQG`

```@docs
GeophysicalFlows.MultiLayerQG.Params
GeophysicalFlows.MultiLayerQG.SingleLayerParams
GeophysicalFlows.MultiLayerQG.Vars
GeophysicalFlows.MultiLayerQG.DecayingVars
GeophysicalFlows.MultiLayerQG.ForcedVars
GeophysicalFlows.MultiLayerQG.StochasticForcedVars
```


## `SurfaceQG`

```@docs
GeophysicalFlows.SurfaceQG.Params
GeophysicalFlows.SurfaceQG.Vars
GeophysicalFlows.SurfaceQG.DecayingVars
GeophysicalFlows.SurfaceQG.ForcedVars
GeophysicalFlows.SurfaceQG.StochasticForcedVars
```


## `BarotropicQGQL`

```@docs
GeophysicalFlows.BarotropicQGQL.Params
GeophysicalFlows.BarotropicQGQL.Vars
GeophysicalFlows.BarotropicQGQL.DecayingVars
GeophysicalFlows.BarotropicQGQL.ForcedVars
GeophysicalFlows.BarotropicQGQL.StochasticForcedVars
```
### a placeholder directory for output generated by Docs# BarotropicQGQL

### Basic Equations

This module solves the *quasi-linear* quasi-geostrophic barotropic vorticity equation on a beta 
plane of variable fluid depth ``H - h(x, y)``. Quasi-linear refers to the dynamics that *neglects* 
the eddy--eddy interactions in the eddy evolution equation after an eddy--mean flow decomposition, 
e.g., 

```math
\phi(x, y, t) = \overline{\phi}(y, t) + \phi'(x, y, t) ,
```

where overline above denotes a zonal mean, ``\overline{\phi}(y, t) = \int \phi(x, y, t) \, 𝖽x / L_x``, and prime denotes deviations from the zonal mean. This approximation is used in many process-model studies of zonation, e.g., 

- Farrell, B. F. and Ioannou, P. J. (2003). [Structural stability of turbulent jets.](http://doi.org/10.1175/1520-0469(2003)060<2101:SSOTJ>2.0.CO;2) *J. Atmos. Sci.*, **60**, 2101-2118.
- Srinivasan, K. and Young, W. R. (2012). [Zonostrophic instability.](http://doi.org/10.1175/JAS-D-11-0200.1) *J. Atmos. Sci.*, **69 (5)**, 1633-1656.
- Constantinou, N. C., Farrell, B. F., and Ioannou, P. J. (2014). [Emergence and equilibration of jets in beta-plane turbulence: applications of Stochastic Structural Stability Theory.](http://doi.org/10.1175/JAS-D-13-076.1) *J. Atmos. Sci.*, **71 (5)**, 1818-1842.

As in the [SingleLayerQG module](singlelayerqg.md), the flow is obtained through a 
streamfunction ``\psi`` as ``(u, v) = (-\partial_y \psi, \partial_x \psi)``. All flow fields 
can be obtained from the quasi-geostrophic potential vorticity (QGPV). Here, the QGPV is

```math
\underbrace{f_0 + \beta y}_{\text{planetary PV}} + \underbrace{\partial_x v
	- \partial_y u}_{\text{relative vorticity}} + \underbrace{\frac{f_0 h}{H}}_{\text{topographic PV}} .
```

The dynamical variable is the component of the vorticity of the flow normal to the plane of 
motion, ``\zeta \equiv \partial_x v - \partial_y u = \nabla^2 \psi``. Also, we denote the 
topographic PV with ``\eta \equiv f_0 h / H``. After we apply the eddy-mean flow decomposition 
above, the QGPV dynamics are:

```math
\begin{aligned}
	\partial_t \overline{\zeta} + \mathsf{J}(\overline{\psi}, \overline{\zeta} + \overline{\eta}) + \overline{\mathsf{J}(\psi', \zeta' + \eta')} & = \underbrace{- \left[\mu + \nu(-1)^{n_\nu} \nabla^{2n_\nu}
	\right] \overline{\zeta} }_{\textrm{dissipation}} , \\
	\partial_t \zeta'  + \mathsf{J}(\psi', \overline{\zeta} + \overline{\eta}) + \mathsf{J}(\overline{\psi}, \zeta' + \eta') + & \underbrace{\mathsf{J}(\psi', \zeta' + \eta') - \overline{\mathsf{J}(\psi', \zeta' + \eta')}}_{\textrm{EENL}} + \beta \partial_x \psi' = \\
	& = \underbrace{-\left[\mu + \nu(-1)^{n_\nu} \nabla^{2n_\nu} \right] \zeta'}_{\textrm{dissipation}} + F .
\end{aligned}
```

where ``\mathsf{J}(a, b) = (\partial_x a)(\partial_y b) - (\partial_y a)(\partial_x b)``. On 
the right hand side, ``F(x, y, t)`` is forcing (which is assumed to have zero zonal mean, 
``\overline{F} = 0``), ``\mu`` is linear drag, and ``\nu`` is hyperviscosity. Plain old 
viscosity corresponds to ``n_{\nu} = 1``.

*Quasi-linear* dynamics **neglects the term eddy-eddy nonlinearity (EENL) term** above.


### Implementation

The equation is time-stepped forward in Fourier space:

```math
\partial_t \widehat{\zeta} = - \widehat{\mathsf{J}(\psi, \zeta + \eta)}^{\textrm{QL}} + \beta \frac{i k_x}{|𝐤|^2} \widehat{\zeta} - \left ( \mu + \nu |𝐤|^{2n_\nu} \right ) \widehat{\zeta} + \widehat{F} .
```

The state variable `sol` is the Fourier transform of vorticity, [`ζh`](@ref GeophysicalFlows.BarotropicQGQL.Vars).

The Jacobian is computed in the conservative form: ``\mathsf{J}(f, g) = \partial_y 
[ (\partial_x f) g] - \partial_x [ (\partial_y f) g]``. The superscript QL on the Jacobian term 
above denotes that triad interactions that correspond to the EENL term are removed.

The linear operator is constructed in `Equation`

```@docs
GeophysicalFlows.BarotropicQGQL.Equation
```

and the nonlinear terms are computed via

```@docs
GeophysicalFlows.BarotropicQGQL.calcN!
```

which in turn calls [`calcN_advection!`](@ref GeophysicalFlows.BarotropicQGQL.calcN_advection!) 
and [`addforcing!`](@ref GeophysicalFlows.BarotropicQGQL.addforcing!).


### Parameters and Variables

All required parameters are included inside [`Params`](@ref GeophysicalFlows.BarotropicQGQL.Params)
and all module variables are included inside [`Vars`](@ref GeophysicalFlows.BarotropicQGQL.Vars).

For decaying case (no forcing, ``F = 0``), `vars` can be constructed with [`DecayingVars`](@ref GeophysicalFlows.BarotropicQGQL.DecayingVars). 
For the forced case (``F \ne 0``) the `vars` struct is with [`ForcedVars`](@ref GeophysicalFlows.BarotropicQGQL.ForcedVars) or [`StochasticForcedVars`](@ref GeophysicalFlows.BarotropicQGQL.StochasticForcedVars).


### Helper functions

```@docs
GeophysicalFlows.BarotropicQGQL.updatevars!
GeophysicalFlows.BarotropicQGQL.set_zeta!
```


### Diagnostics

The kinetic energy of the fluid is obtained via:

```@docs
GeophysicalFlows.BarotropicQGQL.energy
```

while the enstrophy via:

```@docs
GeophysicalFlows.BarotropicQGQL.enstrophy
```

Other diagnostic include: [`dissipation`](@ref GeophysicalFlows.BarotropicQGQL.dissipation), 
[`drag`](@ref GeophysicalFlows.BarotropicQGQL.drag), and [`work`](@ref GeophysicalFlows.BarotropicQGQL.work).


## Examples

- [`examples/barotropicqgql_betaforced.jl`](../../literated/barotropicqgql_betaforced/): A script that simulates forced-dissipative quasi-linear quasi-geostrophic flow on a beta plane demonstrating zonation. The forcing is temporally delta-correlated and its spatial structure is isotropic with power in a narrow annulus of total radius ``k_f`` in wavenumber space. This example demonstrates that the anisotropic inverse energy cascade is not required for zonation.
# SurfaceQG

### Basic Equations

This module solves the non-dimensional surface quasi-geostrophic (SQG) equation for surface 
buoyancy ``b_s = b(x, y, z=0)``, as described in Capet et al., 2008. The buoyancy and the fluid 
velocity at the surface are related through a streamfunction ``\psi`` via:

```math
(u_s, v_s, b_s) = (-\partial_y \psi, \partial_x \psi, -\partial_z \psi) .
```

The SQG model evolves the surface buoyancy,

```math
\partial_t b_s + \mathsf{J}(\psi, b_s) = \underbrace{-\nu(-1)^{n_\nu} \nabla^{2n_\nu} b_s}_{\textrm{buoyancy diffusion}} + \underbrace{F}_{\textrm{forcing}} .
```

Above, ``\mathsf{J}(\psi, b) = (\partial_x \psi)(\partial_y b) - (\partial_y \psi)(\partial_x b)`` 
is the two-dimensional Jacobian. The evolution of buoyancy is only solved for the surface 
layer, but ``b_s`` is a function of the vertical gradient of ``\psi``. In the SQG system, the 
potential vorticity in the interior of the flow is identically zero. That is, relative vorticity 
is equal and opposite to the vertical stretching of the buoyancy layers,

```math
\underbrace{\left(\partial_x^2 + \partial_y^2 \right) \psi}_{\textrm{relative vorticity}} + \underbrace{\partial_z^2 \psi}_{\textrm{stretching term}} = 0 ,
```

with the boundary conditions ``b_s = - \partial_z \psi|_{z=0}`` and ``\psi \rightarrow 0`` as ``z \rightarrow -\infty``. (We take here the oceanographic convention: ``z \le 0``.)

These equations describe a system where the streamfunction (and hence the dynamics) at all depths is prescribed entirely by the surface buoyancy. By taking the Fourier transform in the horizontal (``x`` and ``y``), the streamfunction-buoyancy relation is:

```math
\widehat{\psi}(k_x, k_y, z, t) = - \frac{\widehat{b_s}}{|𝐤|} \, e^{|𝐤|z} , 
```

where ``|𝐤| = \sqrt{k_x^2 + k_y^2}`` is the total horizontal wavenumber.

### Implementation

The buoyancy equation is time-stepped forward in Fourier space:

```math
\partial_t \widehat{b_s} = - \widehat{\mathsf{J}(\psi, b_s)} - \nu |𝐤|^{2 n_\nu} \widehat{b_s} + \widehat{F} .
```

The surface buoyancy is [`b`](@ref GeophysicalFlows.SurfaceQG.Vars). The state variable 
`sol` is the Fourier transform of the surface buoyancy, [`bh`](@ref GeophysicalFlows.SurfaceQG.Vars).

The Jacobian is computed in the conservative form: ``\mathsf{J}(f, g) =
\partial_y [ (\partial_x f) g] -\partial_x[ (\partial_y f) g]``.

The linear operator is constructed in `Equation`

```@docs
GeophysicalFlows.SurfaceQG.Equation
```

while the nonlinear terms via 

```@docs
GeophysicalFlows.SurfaceQG.calcN!
```

which in turn calls [`calcN_advection!`](@ref GeophysicalFlows.SurfaceQG.calcN_advection!) 
and [`addforcing!`](@ref GeophysicalFlows.SurfaceQG.addforcing!).


### Parameters and Variables

All required parameters are included inside [`Params`](@ref GeophysicalFlows.SurfaceQG.Params)
and all module variables are included inside [`Vars`](@ref GeophysicalFlows.SurfaceQG.Vars).

For decaying case (no forcing, ``F = 0``), `vars` can be constructed with [`DecayingVars`](@ref GeophysicalFlows.SurfaceQG.DecayingVars). 
For the forced case (``F \ne 0``) the `vars` struct is with [`ForcedVars`](@ref GeophysicalFlows.SurfaceQG.ForcedVars) or [`StochasticForcedVars`](@ref GeophysicalFlows.SurfaceQG.StochasticForcedVars).


### Helper functions

```@docs
GeophysicalFlows.SurfaceQG.updatevars!
GeophysicalFlows.SurfaceQG.set_b!
```


### Diagnostics

```@docs
GeophysicalFlows.SurfaceQG.kinetic_energy
GeophysicalFlows.SurfaceQG.buoyancy_variance
```

Other diagnostic include: [`buoyancy_dissipation`](@ref GeophysicalFlows.SurfaceQG.buoyancy_dissipation) and
[`buoyancy_work`](@ref GeophysicalFlows.SurfaceQG.buoyancy_work).


## Examples

- [`examples/surfaceqg_decaying.jl`](../../literated/surfaceqg_decaying/): A script that simulates decaying surface quasi-geostrophic flow with a prescribed initial buoyancy field, producing an animation of the evolution of the surface buoyancy.

  > Capet, X. et al., (2008). Surface kinetic energy transfer in surface quasi-geostrophic flows. *J. Fluid Mech.*, **604**, 165-174.
# SingleLayerQG

### Basic Equations

This module solves the barotropic or equivalent barotropic quasi-geostrophic vorticity equation 
on a beta plane of variable fluid depth ``H - h(x, y)``. The flow is obtained through a 
streamfunction ``\psi`` as ``(u, v) = (-\partial_y \psi, \partial_x \psi)``. All flow fields 
can be obtained from the quasi-geostrophic potential vorticity (QGPV). Here the QGPV is

```math
	\underbrace{f_0 + \beta y}_{\text{planetary PV}} + \underbrace{\partial_x v
	- \partial_y u}_{\text{relative vorticity}}
	\underbrace{ - \frac{1}{\ell^2} \psi}_{\text{vortex stretching}} + 
	\underbrace{\frac{f_0 h}{H}}_{\text{topographic PV}} ,
```

where ``\ell`` is the Rossby radius of deformation. Purely barotropic dynamics corresponds to 
infinite Rossby radius of deformation (``\ell = \infty``), while a flow with a finite Rossby 
radius follows is said to obey equivalent-barotropic dynamics. We denote the sum of the relative
vorticity and the vortex stretching contributions to the QGPV with ``q \equiv \nabla^2 \psi - \psi / \ell^2``.
Also, we denote the topographic PV with ``\eta \equiv f_0 h / H``.

The dynamical variable is ``q``.  Thus, the equation solved by the module is:

```math
\partial_t q + \mathsf{J}(\psi, q + \eta) + \beta \partial_x \psi = 
\underbrace{-\left[\mu + \nu(-1)^{n_\nu} \nabla^{2n_\nu} \right] q}_{\textrm{dissipation}} + F .
```

where ``\mathsf{J}(a, b) = (\partial_x a)(\partial_y b)-(\partial_y a)(\partial_x b)`` is the 
two-dimensional Jacobian. On the right hand side, ``F(x, y, t)`` is forcing, ``\mu`` is 
linear drag, and ``\nu`` is hyperviscosity of order ``n_\nu``. Plain old viscosity corresponds 
to ``n_\nu = 1``.


### Implementation

The equation is time-stepped forward in Fourier space:

```math
\partial_t \widehat{q} = - \widehat{\mathsf{J}(\psi, q + \eta)} + \beta \frac{i k_x}{|𝐤|^2 + 1/\ell^2} \widehat{q} - \left(\mu + \nu |𝐤|^{2n_\nu} \right) \widehat{q} + \widehat{F} .
```

The state variable `sol` is the Fourier transform of the sum of relative vorticity and vortex 
stretching (when the latter is applicable), [`qh`](@ref GeophysicalFlows.SingleLayerQG.Vars).

The Jacobian is computed in the conservative form: ``\mathsf{J}(f, g) =
\partial_y [ (\partial_x f) g] - \partial_x[ (\partial_y f) g]``.

The linear operator is constructed in `Equation`

```@docs
GeophysicalFlows.SingleLayerQG.Equation
```

The nonlinear terms are computed via

```@docs
GeophysicalFlows.SingleLayerQG.calcN!
```

which in turn calls [`calcN_advection!`](@ref GeophysicalFlows.SingleLayerQG.calcN_advection!) 
and [`addforcing!`](@ref GeophysicalFlows.SingleLayerQG.addforcing!).


### Parameters and Variables

All required parameters are included inside [`Params`](@ref GeophysicalFlows.SingleLayerQG.Params)
and all module variables are included inside [`Vars`](@ref GeophysicalFlows.SingleLayerQG.Vars).

For decaying case (no forcing, ``F=0``), `vars` can be constructed with [`DecayingVars`](@ref GeophysicalFlows.SingleLayerQG.DecayingVars). 
For the forced case (``F \ne 0``) the `vars` struct is with [`ForcedVars`](@ref GeophysicalFlows.SingleLayerQG.ForcedVars) or [`StochasticForcedVars`](@ref GeophysicalFlows.SingleLayerQG.StochasticForcedVars).


### Helper functions

Some helper functions included in the module are:

```@docs
GeophysicalFlows.SingleLayerQG.updatevars!
GeophysicalFlows.SingleLayerQG.set_q!
```


### Diagnostics

The kinetic energy of the fluid is computed via:

```@docs
GeophysicalFlows.SingleLayerQG.kinetic_energy
```

while the potential energy, for an equivalent barotropic fluid, is computed via:

```@docs
GeophysicalFlows.SingleLayerQG.potential_energy
```

The total energy is:

```@docs
GeophysicalFlows.SingleLayerQG.energy
```

Other diagnostic include: [`energy_dissipation`](@ref GeophysicalFlows.SingleLayerQG.energy_dissipation), 
[`energy_drag`](@ref GeophysicalFlows.SingleLayerQG.energy_drag), [`energy_work`](@ref GeophysicalFlows.SingleLayerQG.energy_work), 
[`enstrophy_dissipation`](@ref GeophysicalFlows.SingleLayerQG.enstrophy_dissipation), and
[`enstrophy_drag`](@ref GeophysicalFlows.SingleLayerQG.enstrophy_drag), [`enstrophy_work`](@ref GeophysicalFlows.SingleLayerQG.enstrophy_work).


## Examples

- [`examples/singlelayerqg_betadecay.jl`](../../literated/singlelayerqg_betadecay/): A script that simulates decaying quasi-geostrophic flow on a beta plane demonstrating zonation.

- [`examples/singlelayerqg_betaforced.jl`](../../literated/singlelayerqg_betaforced/): A script that simulates forced-dissipative quasi-geostrophic flow on a beta plane demonstrating zonation. The forcing is temporally delta-correlated with isotropic spatial structure with power in a narrow annulus in wavenumber space with total wavenumber ``k_f``.

- [`examples/singlelayerqg_decay_topography.jl`](../../literated/singlelayerqg_decay_topography/): A script that simulates two dimensional turbulence (barotropic quasi-geostrophic flow with ``\beta=0``) above topography.

- [`examples/singlelayerqg_decaying_barotropic_equivalentbarotropic.jl`](../../literated singlelayerqg_decaying_barotropic_equivalentbarotropic/): A script that simulates two dimensional turbulence (``\beta=0``) with both infinite and finite Rossby radius of deformation and compares the evolution of the two.# TwoDNavierStokes


### Basic Equations

This module solves two-dimensional incompressible Navier-Stokes equations using the 
vorticity-streamfunction formulation. The flow ``\bm{u} = (u, v)`` is obtained through a 
streamfunction ``\psi`` as ``(u, v) = (-\partial_y \psi, \partial_x \psi)``. The only non-zero 
component of vorticity is that normal to the plane of motion, 
``\partial_x v - \partial_y u = \nabla^2 \psi``. The module solves the two-dimensional 
vorticity equation:

```math
\partial_t \zeta + \mathsf{J}(\psi, \zeta) = \underbrace{-\left [ \mu (-\nabla^2)^{n_\mu}
+ \nu (-\nabla^2)^{n_\nu} \right ] \zeta}_{\textrm{dissipation}} + F ,
```

where ``\mathsf{J}(\psi, \zeta) = (\partial_x \psi)(\partial_y \zeta) - (\partial_y \psi)(\partial_x \zeta)`` 
is the two-dimensional Jacobian and ``F(x, y, t)`` is forcing. The Jacobian term is the advection
of relative vorticity, ``\mathsf{J}(ψ, ζ) = \bm{u \cdot \nabla} \zeta``. Both ``ν`` and ``μ`` 
terms are viscosities; typically the former is chosen to act at small scales (``n_ν ≥ 1``), 
while the latter at large scales (``n_ν ≤ 0``). Plain old viscocity corresponds to ``n_ν=1`` 
while ``n_μ=0`` corresponds to linear drag. Values of ``n_ν ≥ 2`` or ``n_μ ≤ -1`` are referred 
to as hyper- or hypo-viscosities, respectively.


### Implementation

The equation is time-stepped forward in Fourier space:

```math
\partial_t \widehat{\zeta} = - \widehat{\mathsf{J}(\psi, \zeta)} - \left ( \mu |𝐤|^{2n_\mu}
+ \nu |𝐤|^{2n_\nu} \right ) \widehat{\zeta} + \widehat{F} .
```

The state variable `sol` is the Fourier transform of vorticity, [`ζh`](@ref GeophysicalFlows.TwoDNavierStokes.Vars).

The Jacobian is computed in the conservative form: ``\mathsf{J}(a, b) = 
\partial_y [(\partial_x a) b] - \partial_x[(\partial_y a) b]``.

The linear operator is constructed in `Equation`

```@docs
GeophysicalFlows.TwoDNavierStokes.Equation
```

The nonlinear terms are computed via

```@docs
GeophysicalFlows.TwoDNavierStokes.calcN!
```

which in turn calls [`calcN_advection!`](@ref GeophysicalFlows.TwoDNavierStokes.calcN_advection!) 
and [`addforcing!`](@ref GeophysicalFlows.TwoDNavierStokes.addforcing!).


### Parameters and Variables

All required parameters are included inside [`Params`](@ref GeophysicalFlows.TwoDNavierStokes.Params)
and all module variables are included inside [`Vars`](@ref GeophysicalFlows.TwoDNavierStokes.Vars).

For decaying case (no forcing, ``F=0``), `vars` can be constructed with [`Vars`](@ref GeophysicalFlows.TwoDNavierStokes.Vars). 
For the forced case (``F \ne 0``) the `vars` struct is with [`ForcedVars`](@ref GeophysicalFlows.TwoDNavierStokes.ForcedVars) or [`StochasticForcedVars`](@ref GeophysicalFlows.TwoDNavierStokes.StochasticForcedVars).


### Helper functions

Some helper functions included in the module are:

```@docs
GeophysicalFlows.TwoDNavierStokes.updatevars!
GeophysicalFlows.TwoDNavierStokes.set_ζ!
```


### Diagnostics

Some useful diagnostics are:

```@docs
GeophysicalFlows.TwoDNavierStokes.energy
GeophysicalFlows.TwoDNavierStokes.enstrophy
```

Other diagnostic include: [`energy_dissipation`](@ref GeophysicalFlows.TwoDNavierStokes.energy_dissipation), 
[`energy_work`](@ref GeophysicalFlows.TwoDNavierStokes.energy_work), 
[`enstrophy_dissipation`](@ref GeophysicalFlows.TwoDNavierStokes.enstrophy_dissipation), and
[`enstrophy_work`](@ref GeophysicalFlows.TwoDNavierStokes.enstrophy_work).


## Examples

- [`examples/twodnavierstokes_decaying.jl`](../../literated/twodnavierstokes_decaying/): A script that simulates decaying two-dimensional turbulence reproducing the results by

  > McWilliams, J. C. (1984). The emergence of isolated coherent vortices in turbulent flow. *J. Fluid Mech.*, **146**, 21-43.

- [`examples/twodnavierstokes_stochasticforcing.jl`](../../literated/twodnavierstokes_stochasticforcing/): A script that simulates forced-dissipative two-dimensional turbulence with isotropic temporally delta-correlated stochastic forcing.

- [`examples/twodnavierstokes_stochasticforcing_budgets.jl`](../../literated/twodnavierstokes_stochasticforcing_budgets/): A script that simulates forced-dissipative two-dimensional turbulence demonstrating how we can compute the energy and enstrophy budgets.
# MultiLayerQG

### Basic Equations

This module solves the layered quasi-geostrophic equations on a beta plane of variable fluid 
depth ``H - h(x, y)``. The flow in each layer is obtained through a streamfunction ``\psi_j`` as 
``(u_j, v_j) = (-\partial_y \psi_j, \partial_x \psi_j)``, ``j = 1, \dots, n``, where ``n`` 
is the number of fluid layers.

The QGPV in each layer is

```math
\mathrm{QGPV}_j = q_j + \underbrace{f_0 + \beta y}_{\textrm{planetary PV}} + \delta_{j, n} \underbrace{\frac{f_0 h}{H_n}}_{\textrm{topographic PV}}, \quad j = 1, \dots, n .
```

where ``q_j`` incorporates the relative vorticity in each layer ``\nabla^2 \psi_j`` and the 
vortex stretching terms:

```math
q_1 = \nabla^2 \psi_1 + F_{3/2, 1} (\psi_2 - \psi_1) ,\\
q_j = \nabla^2 \psi_j + F_{j-1/2, j} (\psi_{j-1} - \psi_j) + F_{j+1/2, j} (\psi_{j+1} - \psi_j) , \quad j = 2, \dots, n-1 ,\\
q_n = \nabla^2 \psi_n + F_{n-1/2, n} (\psi_{n-1} - \psi_n) .
```

with

```math
F_{j+1/2, k} = \frac{f_0^2}{g'_{j+1/2} H_k} \quad \text{and} \quad
g'_{j+1/2} = g \frac{\rho_{j+1} - \rho_j}{\rho_{j+1}} .
```

In view of the relationships above, when we convert to Fourier space ``q``'s and ``\psi``'s are 
related via the matrix equation

```math
\begin{pmatrix} \widehat{q}_{𝐤, 1}\\\vdots\\\widehat{q}_{𝐤, n} \end{pmatrix} =
\underbrace{\left(-|𝐤|^2 \mathbb{1} + \mathbb{F} \right)}_{\equiv \mathbb{S}_{𝐤}}
\begin{pmatrix} \widehat{\psi}_{𝐤, 1}\\\vdots\\\widehat{\psi}_{𝐤, n} \end{pmatrix}
```

where

```math
\mathbb{F} \equiv \begin{pmatrix}
 -F_{3/2, 1} &              F_{3/2, 1}  &   0   &  \cdots    & 0\\
  F_{3/2, 2} & -(F_{3/2, 2}+F_{5/2, 2}) & F_{5/2, 2} &       & \vdots\\
 0           &                  \ddots  &   \ddots   & \ddots & \\
 \vdots      &                          &            &        &  0 \\
 0           &       \cdots             &   0   & F_{n-1/2, n} & -F_{n-1/2, n}
\end{pmatrix} .
```

Including an imposed zonal flow ``U_j(y)`` in each layer, the equations of motion are:

```math
\partial_t q_j + \mathsf{J}(\psi_j, q_j ) + (U_j - \partial_y\psi_j) \partial_x Q_j +  U_j \partial_x q_j  + (\partial_y Q_j)(\partial_x \psi_j) = -\delta_{j, n} \mu \nabla^2 \psi_n - \nu (-1)^{n_\nu} \nabla^{2 n_\nu} q_j ,
```

with

```math
\partial_y Q_j \equiv \beta - \partial_y^2 U_j - (1-\delta_{j,1}) F_{j-1/2, j} (U_{j-1} - U_j) - (1 - \delta_{j,n}) F_{j+1/2, j} (U_{j+1} - U_j) + \delta_{j,n} \partial_y \eta , \\
\partial_x Q_j \equiv \delta_{j, n} \partial_x \eta .
```


### Implementation

Matrices ``\mathbb{S}_{𝐤}`` as well as ``\mathbb{S}^{-1}_{𝐤}`` are included in `params` as 
`params.S` and `params.S⁻¹` respectively. Additionally, the background PV gradients 
``\partial_x Q`` and ``\partial_y Q`` are also included in the `params` as `params.Qx` and 
`params.Qy`.

One can get the ``\widehat{\psi}_j`` from ``\widehat{q}_j`` via 
`streamfunctionfrompv!(psih, qh, params, grid)`, while the inverse, i.e. obtain ``\widehat{q}_j`` from ``\widehat{\psi}_j``, is done via  `pvfromstreamfunction!(qh, psih, params, grid)`.

The equations of motion are time-stepped forward in Fourier space:

```math
\partial_t \widehat{q}_j = - \widehat{\mathsf{J}(\psi_j, q_j)}  - \widehat{U_j \partial_x Q_j} - \widehat{U_j \partial_x q_j}
+ \widehat{(\partial_y \psi_j) \partial_x Q_j}  - \widehat{(\partial_x \psi_j)(\partial_y Q_j)} + \delta_{j, n} \mu |𝐤|^{2} \widehat{\psi}_n - \nu |𝐤|^{2n_\nu} \widehat{q}_j .
```

In doing so the Jacobian is computed in the conservative form: ``\mathsf{J}(f,g) =
\partial_y [ (\partial_x f) g] - \partial_x[ (\partial_y f) g]``.

The state variable `sol` consists of the Fourier transforms of ``q_j`` at each layer, i.e., 
[`qh`](@ref GeophysicalFlows.MultiLayerQG.Vars).

The linear operator is constructed in `Equation`

```@docs
GeophysicalFlows.MultiLayerQG.Equation
GeophysicalFlows.MultiLayerQG.hyperviscosity
```

The nonlinear terms are computed via

```@docs
GeophysicalFlows.MultiLayerQG.calcN!
```

which in turn calls [`calcN_advection!`](@ref GeophysicalFlows.MultiLayerQG.calcN_advection!) 
and [`addforcing!`](@ref GeophysicalFlows.MultiLayerQG.addforcing!).

!!! tip "Linearized MultiLayerQG dynamics "
    The `MultiLayerQG` module includes also a linearized version of the dynamics about a base
    flow ``U_j(y)``, ``j = 1, \dots, n``; see [`LinearEquation`](@ref GeophysicalFlows.MultiLayerQG.LinearEquation), 
    [`calcNlinear!`](@ref GeophysicalFlows.MultiLayerQG.calcNlinear!), and 
    [`calcN_linearadvection!`](@ref GeophysicalFlows.MultiLayerQG.calcN_linearadvection!).


### Parameters and Variables

All required parameters are included inside [`Params`](@ref GeophysicalFlows.MultiLayerQG.Params)
and all module variables are included inside [`Vars`](@ref GeophysicalFlows.MultiLayerQG.Vars).

For decaying case (no forcing, ``F=0``), `vars` can be constructed with [`DecayingVars`](@ref GeophysicalFlows.MultiLayerQG.DecayingVars). 
For the forced case (``F \ne 0``) the `vars` struct is with [`ForcedVars`](@ref GeophysicalFlows.MultiLayerQG.ForcedVars) or [`StochasticForcedVars`](@ref GeophysicalFlows.MultiLayerQG.StochasticForcedVars).


### Helper functions

```@docs
GeophysicalFlows.MultiLayerQG.set_q!
GeophysicalFlows.MultiLayerQG.set_ψ!
GeophysicalFlows.MultiLayerQG.updatevars!
```

### Diagnostics

The eddy kinetic energy in each layer and the eddy potential energy that corresponds to each 
fluid interface is computed via `energies()`:

```@docs
GeophysicalFlows.MultiLayerQG.energies
```

The lateral eddy fluxes in each layer and the vertical fluxes across fluid interfaces are
computed via `fluxes()`:

```@docs
GeophysicalFlows.MultiLayerQG.fluxes
```


## Examples

 - [`examples/multilayerqg_2layer.jl`](../../literated/multilayerqg_2layer/): A script that simulates the growth and equilibration of baroclinic eddy turbulence in the Phillips 2-layer model.
