<a href="https://github.com/JuliaOcean/AIBECS.jl">
  <img src="https://user-images.githubusercontent.com/4486578/60554111-8fc27400-9d79-11e9-9ca7-6d78ee89ea70.png" alt="logo" title="The AIBECS logo: It represents three global marine biogeochemical cycles, where each element affects the others" align="center" width="50%"/>
</a>

# AIBECS.jl

*The ideal tool for exploring global marine biogeochemical cycles.*

<p>
  <a href="https://JuliaOcean.github.io/AIBECS.jl/stable/">
    <img src="https://img.shields.io/github/workflow/status/JuliaOcean/AIBECS.jl/Documentation?style=for-the-badge&label=Documentation&logo=Read%20the%20Docs&logoColor=white">
  </a>
  <a href="https://doi.org/10.21105/joss.03814">
    <img src="https://img.shields.io/static/v1?label=JOSS&message=10.21105/joss.03814&color=9cf&style=flat-square" alt="DOI badge">
  </a>
  <a href="https://www.bpasquier.com/talk/osm_sandiego_2020/OSM_SanDiego_2020.pdf">
    <img src=https://img.shields.io/static/v1?label=Poster&message=OSM2020&color=9cf&style=flat-square>
  </a>
</p>

<p>
  <a href="https://doi.org/10.5281/zenodo.2864051">
    <img src="http://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.2864051-blue.svg?&style=flat-square">
  </a>
  <a href="https://github.com/JuliaOcean/AIBECS.jl/blob/master/LICENSE">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-blue.svg?&style=flat-square">
  </a>
</p>

<p>
  <a href="https://github.com/JuliaOcean/AIBECS.jl/actions">
    <img src="https://img.shields.io/github/workflow/status/JuliaOcean/AIBECS.jl/Mac%20OS%20X?label=OSX&logo=Apple&logoColor=white&style=flat-square">
  </a>
  <a href="https://github.com/JuliaOcean/AIBECS.jl/actions">
    <img src="https://img.shields.io/github/workflow/status/JuliaOcean/AIBECS.jl/Linux?label=Linux&logo=Linux&logoColor=white&style=flat-square">
  </a>
  <a href="https://github.com/JuliaOcean/AIBECS.jl/actions">
    <img src="https://img.shields.io/github/workflow/status/JuliaOcean/AIBECS.jl/Windows?label=Windows&logo=Windows&logoColor=white&style=flat-square">
  </a>
  <a href="https://codecov.io/gh/JuliaOcean/AIBECS.jl">
    <img src="https://img.shields.io/codecov/c/github/JuliaOcean/AIBECS.jl/master?label=Codecov&logo=codecov&logoColor=white&style=flat-square">
  </a>
</p>




**AIBECS** (for **A**lgebraic **I**mplicit **B**iogeochemical **E**lemental **C**ycling **S**ystem, pronounced like the cool [ibex](https://en.wikipedia.org/wiki/Ibex)) is a Julia package that provides ocean biogeochemistry modellers with an easy-to-use interface for creating and running models of the ocean system.

AIBECS is a system because it allows you to choose some biogeochemical tracers, define their interactions, select an ocean circulation and *Voilà!* — your model is ready to run.

## Getting started

If you are new to AIBECS, head over to the [documentation](https://JuliaOcean.github.io/AIBECS.jl/stable/) and look for the tutorials.
(You can also click on the big "Documentation" badge above.)

## Concept

This package was developed to exploit linear-algebra tools and algorithms in Julia to efficiently simulate marine tracers.
AIBECS represents global biogeochemical cycles with a discretized system of nonlinear ordinary differential equations that takes the generic form

∂***x***/∂*t* + **T*****x*** = ***G***(***x***)

where ***x*** represents the model state variables, i.e., the marine tracer(s) concentration.
For a single tracer, ***x*** can be interpreted as the 3D field of its concentration.
In AIBECS, ***x*** is represented as a column vector (that's why it's **bold** and *italic*).

The operator **T** is a spatial differential operator that represents the transport of tracers.
For example, for a single tracer transported by ocean circulation,

**T** = ∇⋅(***u*** - **K**∇)

represents the effects of advection and eddy-diffusion.
(***u*** is the 3D vector of the marine currents and **K** is a 3×3 eddy-diffusivity matrix.)
Thus, **T** "acts" on ***x*** such that **T*****x*** is the flux divergence of that tracer.
In AIBECS, **T** is represented by a matrix (that's why it's **bold** and upstraight).

Lastly, the right-hand-side, ***G***(***x***), represents the local sources minus sinks of each tracer, whcih must be provided as functions of the tracer(s) ***x***.

To simulate tracers using the AIBECS, you just need to define the transport operators ***T*** and the net sources and sinks ***G***.
That's pretty much the whole concept!

## References

If you use this package, please cite it.

If you use data provided by this package (like the ocean circulation from the OCIM), please cite them as well.

For convenience, all the references are available in [BibTeX](https://en.wikipedia.org/wiki/BibTeX) format in the [CITATION.bib](./CITATION.bib) file.

Also, if you want to do research using the AIBECS, and you think I could help, do not hesitate to contact me directly (contacts on my [website](www.bpasquier.com)), I would be happy to contribute and collaborate!

<img src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png" alt="NSF" title="NSF_logo" align="right" height="50"/>

The authors acknowledge funding from the Department of Energy grant DE-SC0016539 and from the National Science Foundation grant 1658380.
---
title: 'AIBECS.jl: A tool for exploring global marine biogeochemical cycles.'
tags:
  - Julia
  - Oceanography
  - Biogeochemical Cycles
  - Earth Sciences
authors:
  - name: Benoît Pasquier^[corresponding author]
    orcid: 0000-0002-3838-5976
    affiliation: "1,3" # (Multiple affiliations must be quoted)
  - name: François W. Primeau
    orcid: 0000-0001-7452-9415
    affiliation: "2"
  - name: Seth G. John
    orcid: 0000-0002-8257-626X
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
  - name: Department of Earth Sciences, University of Southern California
    index: 1
  - name: Department of Earth System Science, University of California, Irvine
    index: 2
  - name: (now at) School of Mathematics and Statistics, University of New South Wales, Sydney
    index: 3
date: 7 Dec 2021
bibliography: paper.bib

---

# Summary




The ocean transports, mixes, and transforms chemical constituents on a multitude of time and length scales.
Observations and models are both essential for making sense of this complex system.
However, the days of just publishing data without any quantitative modelling are over, increasing pressure on sea-going oceanographers that are expected to be proficient in the use of biogeochemical models.
And while the simplest of these models consist of only a few boxes [e.g., the two-box model of @Archer_etal_GBC_2000], with solutions that are easily obtained by hand or on a simple desktop computer, the most advanced models have high-resolution three-dimensional meshes that require high-performance computing (HPC) clusters and a considerable amount of computational-science expertise [e.g., the MITgcm, @Campin_etal_2021].
A high barrier to entry for modelling in oceanography hinders advances in the field.
[AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl) [@Pasquier:2021] is a [JuliaOcean](https://github.com/JuliaOcean/)-affiliated package written in [Julia](https://julialang.org/) [@Bezanson_etal:2017], which aims to lower this barrier by providing an easy-to-use, open-source, and modular framework for simulating global marine tracers in steady-state and biogeochemical parameter fitting.






A conceptual model of the cycle of any marine tracer essentially requires two components.
(i) A model of how the tracer is transported (be it the ocean currents and eddies, gravitational settling, or a combination of those), and (ii) a model of the local sources and sinks at any location.
AIBECS.jl is built on this concept and allows users to build numerical models of marine tracers by selecting a circulation and/or vertical transport of the tracer with particles and local sources and sinks.
Tools for generating the steady-state equations, solving them, and then diagnosing and plotting the simulated tracers are also provided, either directly by AIBECS.jl, by its dependencies, or by satellite packages in the AIBECS and Julia ecosystem.



@Box_1979 has remarked that, by definition, all models are wrong because they make simplifying assumptions.
These assumptions can give rise to model parameters that can capture fundamental characteristics of the system.
Once the model of the cycle of a marine tracer is implemented, it is thus natural to try to estimate biogeochemical parameters by minimizing model–observation mismatches, which may improves the skill of the model.
Parameter optimization is generally more efficiently performed when first and second derivatives are available, and these derivatives can also be used to estimate parameter sensitivity [@Thacker_JGRO_1989].
AIBECS.jl was designed with parameter fitting in mind and provides, along with satellite packages, tools for generating first- and second-order derivatives automatically and efficiently.



The AIBECS.jl package builds on the vision of the [AWESOME OCIM](https://github.com/MTEL-USC/AWESOME-OCIM) [A Working Environment for Simulating Ocean Movement and Elemental cycling within the Ocean Circulation Inverse Model, or AO, @John_etal_ChemGeo_2020], which was designed to make three-dimensional element-cycling models more broadly accessible.
Written in [Julia](https://julialang.org/) — chosen for of its combined expressive power and efficiency and offering a truly open-source solution — the AIBECS.jl framework improves on the AO on a number of fronts to target education and research.
Simple AIBECS.jl simulations can be produced in minutes in interactive notebooks while resource-hungry projects, e.g., for optimizing models with multiple tracers interacting nonlinearly within a high-resolution mesh, can be easily version-controlled, hosted, and run on HPC clusters.




Through a simple user interface, AIBECS.jl provides access to a variety of steady-state ocean circulation models.
These currently include
the Ocean Circulation Inverse Model (OCIM) v0.1 [@Primeau_etal_JGRO_2013],
v1.0 [@DeVries_Primeau_JPO_2011; @DeVries_GBC_2014],
and v2.0 [@DeVries_Holzer_JGRO_2019],
and the MITgcm-built Ocean Comprehensible Atlas (OCCA) ocean-state estimate model [@Forget:2010], of which the downloads are handled by the [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl) package [@White_etal:2019].
AIBECS.jl also offers classic two-box and three-box models [@Sarmiento_Gruber:2006; @Archer_etal_GBC_2000].
The [OceanGrids.jl](https://github.com/briochemc/OceanGrids.jl) package, on which AIBECS depends, provides the underlying grid configuration types as well as regridding and interpolating routines.
Swapping the underlying circulation model and grid requires a single-line-of-code change, facilitating intercomparison projects.
As new circulation models that are represented in matrix form are made publicly available, they will be added to the collection.
These could include past- and future-ocean circulation models, for paleoceanographic or climate-change studies, for example.





AIBECS.jl also provides extra functionality to facilitate the generation of numerical models.
Tooling to simulate gravitational settling of tracers with non-buoyant particles is provided by the `transportoperator` function.
In addition, AIBECS provides access to a number of predefined fields that can be used to generate source and sink processes.
Fine-resolution (1-arc-minute) topography from the ETOPO1 dataset [@Amante_Eakins_2009] can be used for a refined interception of particulate fluxes by subgrid topographic features not captured by coarser circulation models.
For aeolian deposition, AIBECS.jl includes aerosol-type- and region-of-origin-partitioned dust deposition fields [@Chien_etal_2016; @Kok_etal_2021b].
Datasets for global river discharge [@Dai_2017; @Dai_Trenberth_2002] and surface groundwater discharge [@Luijendijk_etal_2019; @Luijendijk_etal_2020] are included.
For hydrothermal-sourced tracers, the helium fluxes from the Earth's mantle computed with the OCIM v1.0 and v2.0 are available when loading the corresponding circulation models [@DeVries_GBC_2014; @DeVries_Holzer_JGRO_2019].
AIBECS.jl also provides access to the data included with the AWESOME OCIM framework [@John_etal_ChemGeo_2020], namely data from the Global Ocean Data Analysis Project [GLODAP, @Lauvset_etal_2016; @Olsen_etal_2016], P-cycling modelled fields from @Weber_etal_Science_2018, nepheloid layers [@Gardner_etal_EPSL_2018; @Gardner_etal_ProgOcn_2018; @Taburet_etal_OcnSci_2019], as well as other data already present within AIBECS or satellite packages.
Also useful to global biogeochemistry modelling are data from the World Ocean Atlas [@WOA_2018_nut] that can be downloaded, assisted by external package [WorldOceanAtlasTools.jl](https://github.com/briochemc/WorldOceanAtlasTools.jl) [@WorldOceanAtlasTools.jl-2019].
Similarly, GEOTRACES data [@Schlitzer_etal_ChemGeo_2018] can be handled by the [GEOTRACES.jl](https://github.com/briochemc/GEOTRACES.jl) package (although GEOTRACES requires manual download of the data).
More advanced usage such as optimization is facilitated by the [F1Method.jl](https://github.com/briochemc/F1Method.jl) package [@F1Method], which provides efficient gradient and Hessian computations of objective functions defined through AIBECS.jl, which can then be directly fed to optimization routines from, e.g., the [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) package [@Optim.jl-2018].
Finally, plotting recipes for the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package [@Plots.jl], are available.






Internally, AIBECS.jl uses a quasi-Newton solver [@Kelley_2003_1] translated from MATLAB to Julia and tailored to the context of marine tracers to solve/simulate tracers.
AIBECS uses forward-mode auto-differentiation from the [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) package [@RevelsLubinPapamarkou2016] for the nonlinear parts of the system of ordinary differential equations to generate the Jacobian required for the solver.
Metadata such as units and prior distributions can be attached to model parameters, which are handled with the help of the [UnPack.jl](https://github.com/mauro3/Unpack.jl), [FieldMetadata.jl](https://github.com/rafaqz/FieldMetadata.jl), [Flatten.jl](https://github.com/rafaqz/Flatten.jl), [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl), [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) [@Besancon_etal_2021; @Distributions.jl-2019], [Unitful.jl](https://github.com/PainterQubits/Unitful.jl), and [Bijectors.jl](https://github.com/TuringLang/Bijectors.jl) dependencies.






The [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl) package is registered with Julia's [default package registry (GENERAL)](https://github.com/JuliaRegistries/General), such that installation takes a single line of code from within Julia.
The package [documentation](https://juliaocean.github.io/AIBECS.jl/stable/), which is built through continuous integration (CI), includes [tutorials](https://juliaocean.github.io/AIBECS.jl/stable/#.-Tutorials) and [how-to guides](https://juliaocean.github.io/AIBECS.jl/stable/#.-How-to-guides) that are both available online for consultation and as runnable Jupyter notebooks.
Continuous integration through GitHub actions also includes a fairly complete suite of tests.







# Statement of need




Pioneered by @Schlitzer_JPO_1993 to study the ocean circulation using ventilation tracers such as radiocarbon [see also @Schlitzer_GeoMonoAGU_2000], the low computational costs of steady-state circulation models allow for efficient optimization and inference/estimation of biogeochemical parameters.
Today, despite reduced resolution and steady-state assumption, data-constrained matrix-transport models such as the OCIM are at the forefront of oceanographic research.
This is evidenced for example by the growing number of high-profile studies that use such models for parameter estimation published in recent years [e.g., @DeVries_etal_Nature_2017; @DeVries_etal_Biogeosciences_2013; @DeVries_etal_GRL_2012; @Weber_Deutsch_Nature_2010; @DeVries_GBC_2014; @DeVries_Weber_GBC_2017; @Teng_etal_NatGeosci_2014; @DeVries_Deutsch_NatGeosci_2014; @Weber_etal_PNAS_2016; @Roshan_DeVries_NatCom_2017; @Wang_etal_Nature_2019].
However, although the steady-state assumption and matrix representation simplify the simulation of tracers compared to traditional Ocean General Circulation Models (OGCMs), most studies that employ a steady-state matrix representation of marine cycling remain difficult to reproduce without significant computer-science and modelling expertise because they are built on private implementations.
Comparisons between different circulation models are moreover complicated by the lack of standardization across models.



Hence, there is a need to facilitate the use of steady-state ocean-circulation models by providing:
(i) an integrated framework for handling a number of different ocean-circulation models with tools for swapping circulations (including interpolating from one model grid to another),
(ii) a user-friendly interface for translating mathematical models of biogeochemical cycles into the corresponding code (e.g., for sources, sinks, and vertical transport of tracers), and
(iii) solvers for efficient simulations, optimization, diagnosis, and statistical analysis.




AIBECS.jl provides a free, open-source, unified framework for biogeochemical-tracer-modelling studies that use steady-state circulation models.
Among other advantages over existing solutions (i.e., the AO), AIBECS.jl offers better computational efficiency, enhanced versatility, composability with other Julia packages, and ease of reproducibility (granted by version control and Julia's package manager) and improved syntax, which are pillars of modern scientific dissemination.
Thus, AIBECS users may include sea-going oceanographers and educators who will benefit from its simplicity, as well as more experienced modellers who can leverage its computational advantages.
AIBECS.jl has been used for teaching and is currently used for research focused on marine trace metals.



The publication of circulation models as transport matrices from existing general circulation models [@Khatiwala_etal_OM_2005; @Khatiwala_GBC_2007; @Bardin_etal_OM_2014; @Bardin_etal_OM_2016; @Kvale_etal_GMD_2017; @Zanna_etal_PNAS_2019], and hopefully, future publication of transport matrix estimates from standard ocean general circulation models [e.g., @Chamberlain_etal_OM_2019] will increase the collection of circulations available from AIBECS.jl.
Of particular interest to the broader community to facilitate the simulation of past and future marine biogeochemical states would be transport matrix models extracted from the Climate Model Intercomparison Project (CMIP), which includes past and future simulations of the ocean circulation.





Further devepment could include exposing advanced Grren-function-based diagnostic tools [e.g., @Holzer_etal_GBC_2021; @Pasquier_Holzer_Biogeosciences_2018], coupling tracers on different grids, Newton–Krylov solvers for cyclo-stationary states [e.g., CYCLOCIM: @Huang_etal_OM_2021], or time-dependent solvers for transient biogeochemical simualtions provided, e.g., by the SciML ecosystem [@Rackauckas_Nie_JORS_2017].
Bridging packages could be implemented for improved composability with statistical packages [e.g., [Turing.jl](https://github.com/TuringLang/Turing.jl): @ge2018t], optimization tools, and plotting software [e.g., [Makie.jl](https://github.com/JuliaPlots/Makie.jl): @Makie.jl; @Danisch_Krumbiegel_2021].



# Acknowledgements

FWP and BP acknowledge funding from the Department of Energy (grant DE-SC0016539) and the National Science Foundation (grant 1658380).
BP and SGJ acknowledge funding provided by the Simons Foundation (Award #426570SP to SGJ) and the National Science Foundation (grant 1736896).
The authors thank Dr. Zhen Wu and Prof. Dan Kelley for their insightful reviews that have helped improve the software and this manuscript.

# ReferencesThis ReadMe is **not** part of the AIBECS documentation.
It is merely for myself to remember how to edit the documentation.

To edit the documentation, it is better to avoid running the full suite of scripts and examples.
In the `generate.jl` script, comment/uncomment the Literate functions that generate the markdown and` notebooks files to not execute them.

Then, from the root of AIBECS, run
```bash
julia --project=docs --color=yes -L docs/live.jl
```
which should start a local website at `http://localhost:8000/` that will host the documentation and refresh on any save of the examples or other documentation source files.

And don't forget to turn the execution on for a last run before you push!
```@raw html
<img src="https://user-images.githubusercontent.com/4486578/60554111-8fc27400-9d79-11e9-9ca7-6d78ee89ea70.png" alt="logo" title="AIBECS_logo" align="middle" width="50%"/>
```

# [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl)

The **Algebraic Implicit Biogeochemistry Elemental Cycling System**

!!! note
    If you are using the AIBECS for the first time, you must add it to your Julia environment, by typing
    ```
    ]add AIBECS
    ```
    at the REPL.

This documentation is organized in 4 parts:

#### 1. Tutorials

If you want to try AIBECS for the first time, this is where you should start.

- The [ideal age tutorial](@ref ideal-age) is a good place to start.
    Generate a simple linear model of an idealized tracer in a few lines of code.
- The [radiocarbon tutorial](@ref radiocarbon) is a little bit more involved, with some nonlinearities and more advanced use of AIBECS features and syntax.
- The [coupled PO₄–POP model tutorial](@ref P-model) couples 2 interacting tracers, dissolved phosphate and particulate organic phosphorus (POP).
- The [dust model tutorial](@ref dust-model) simulates some fictitious metals being injected at the ocean–atmosphere interface and being reversibly scavenged.
- The [river discharge tutorial](@ref river-discharge) similarly simulates another fictitious radioactive tracer that is injected in the ocean by rivers and that decays away as time passes.
- The [groundwater discharge tutorial](@ref groundwater-discharge) is almost identical to the river-discharge tutorial, except it uses groundwater discharge data.

#### 2. How-to guides

Here you will find goal-oriented walk-through's.

- [Parameters guide](@ref parameters)
- [Plotting basic things](@ref plots)
- [Plotting cruise/transects data](@ref cruiseplots)
- [Estimate fluxes](@ref fluxes)
- How to simulate, i.e., solve or timestep your model (coming soon!)
- How to optimize parameters (coming soon!)
- How to simulate sinking particles (coming soon!)

#### 3. Explanation/discussion

Here you will find more general discussions and explanations surrounding the AIBECS.

- [The concept of the AIBECS](@ref concept)
- [Tracer transport operators](@ref tracer-transport-operators)

#### 4. Reference

This section contains the docstrings of (almost all) the functions available in AIBECS.

----

!!! note
    The AIBECS is being developed primarily by Benoît Pasquier (postdoc with Seth John at the University of Southern California).
    If you use the AIBECS in your research, [please cite it](https://doi.org/10.5281/zenodo.2864051)!
    Similarly, if you access data from within AIBECS (like the OCIM or OCCA ocean circulations) please cite them too.

!!! warning
    This package is in active development, so nothing is set in stone, and things may be broken sometimes.
    And if you have any suggestions or feature requests, do not hesitate to [start an issue](https://github.com/JuliaOcean/AIBECS.jl/issues), or better even, [submit a pull request](https://github.com/JuliaOcean/AIBECS.jl/pulls)!
# [Concept](@id concept)

The AIBECS (pronounced like the cool [ibex](https://en.wikipedia.org/wiki/Ibex) if you have a french accent) is a new software written in [Julia](https://julialang.org) to easily create some marine biogeochemistry models in just a few commands.


AIBECS is not just a single model.
It's a **system** that allows you to create a global steady-state biogeochemistry model with just a few simple commands.
Basically, you just need to tell AIBECS which ocean circulation to use (from simple toy models of just a few boxes to more complicated global models of the circulation), what elements you want to model/track, and how they are converted into other (the net local sources and sinks of your model).
Once the model is set up, chose some parameter values and you can run simulations.

AIBECS relies on many tools from linear algebra to run simulations and perform optimizations really fast.
AIBECS-generated models are described by a state function, denoted $\boldsymbol{F}$, which defines how the concentrations of elements in the ocean evolve with time.
In mathematical terms, this translates to a system of nonlinear differential equations with the generic form 

$$\frac{\partial \boldsymbol{x}}{\partial t} = \boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}),$$
where $\boldsymbol{x}$ is the state of the model (i.e., the concentrations of the tracers), and $\boldsymbol{p}$ are model parameters.
With AIBECS, you can efficiently find the equilibrium of the system (AKA the steady-state).
That is when the time-derivative is $0$, so that

$$\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0,$$

and $\boldsymbol{x}$ does not change with time.
Instead of simulating the evolution of $\boldsymbol{x}$ with time and waiting for the system to reach equilibrium — like most biogeochemistry models do — AIBECS uses linear algebra techniques, like Newton's method in multiple dimensions, or Krylov spaces, to implicitly solve for the steady-state solution, hence the "algebraic" and "implicit" names.
This makes AIBECS much faster than the competition!


# [Tracer transport](@id tracer-transport-operators)

## Transport operator

To model marine biogeochemical tracers on a global scale we need to be able to account for their movement within the 3D ocean.
We do this with tracer **transport operators**, generically denoted by $\mathcal{T}$.
These operators *act* on a tracer field to give its divergence.
In other words, take a tracer with concentration $x(\boldsymbol{r})$ at location $\boldsymbol{r}$ that is transported by some mechanism represented by the operator $\mathcal{T}$.
Its divergence at $\boldsymbol{r}$ is then $(\mathcal{T} x)(\boldsymbol{r})$.

In the case of the ocean circulation, i.e., the currents and eddies that transport marine tracers floating around in sea water, the transport operator can be represented by

$$\mathcal{T} = \nabla \cdot \left[ \boldsymbol{u}(\boldsymbol{r}) - \mathbf{K}(\boldsymbol{r}) \nabla \right]$$

The $\boldsymbol{u}(\boldsymbol{r})$ term represents the mean marine current velocity at location $\boldsymbol{r}$.
Thus, $\boldsymbol{u}(\boldsymbol{r})$ is a 3D vector aligned with the ocean currents and whose amplitude, in m/s, gives the velocity of the moving sea water.
The $\mathbf{K}(\boldsymbol{r})$ term is a 3×3 matrix that represents eddy diffusivity.
That is, it reproduces the effective mixing effect of unresolved eddies, the turbulent vortices that are too small relative to the model grid to be explicitly captured.
In earlier models of the ocean circulation, like OCIM's ancestor (denoted OCIM0 in AIBECS), $\mathbf{K}$ was diagonal.
Since OCIM1, $\mathbf{K}$ contains non-diagonal to orient mixing preferentially along isopycnals.

In the case of sinking particles, one can assume that they are only transported downwards with some terminal settling velocity $\boldsymbol{w}(\boldsymbol{r}).$
The corresponding transport operator is then simply

$$\mathcal{T} = \nabla \cdot \boldsymbol{w}(\boldsymbol{r}).$$

## Discretization

In order to represent a marine tracer on a computer one needs to **discretize** the 3D ocean into a discrete grid, i.e., a grid with a finite number of boxes.
One can go a long way towards understanding what a tracer transport operator is by playing with a simple model with only a few boxes, which is the goal of this piece of documentation.

The simple box model we consider is embedded in a 2×2×2 "shoebox".
It has 5 *wet* boxes and 3 *dry* boxes, as illustrated below:

```@raw html
<img src="https://user-images.githubusercontent.com/4486578/58314610-3b130b80-7e53-11e9-9fe8-9527cdcca2d0.png" width =800>
```

An example of discretized advection, $\boldsymbol{u}$, is shown on the image above, and consists of
- a meridional overturning circulation flowing in a cycle through boxes 1 → 2 → 6 → 5 → 1 (shown in the meridional section 1 panel)
- a zonal current in a reentrant cycling through boxes 1 → 3 → 1 (shown in the layer 1 panel)

!!! note
    This circulation is available as the `Primeau_2x2x2` model.
    You can load it in AIBECS via
    ```julia
    grd, T = Primeau_2x2x2.load()
    ```
    Note that this circulation also contains vertical mixing representing deep convection between boxes 2 ↔ 6 (not shown on the image)

## Vectorization

In AIBECS, tracers are represented by column vectors.
That is, the 3D tracer field, $x(\boldsymbol{r})$, is **vectorized**, in the sense that the concentrations in each box are rearranged into a column vector.

```@raw html
<img src="https://user-images.githubusercontent.com/4486578/61757212-fe3ba480-ae02-11e9-8d17-d72866eaafb5.gif" width =800>
```

The continuous transport operator $\mathcal{T}$ can then be represented by a matrix, denoted by $\mathbf{T}$, and sometimes called the *transport matrix*.
It turns out that in most cases, this matrix is sparse.

!!! note
    A sparse matrix behaves the same way as a regular matrix.
    The only difference is that in a sparse matrix the majority of the entries are zeros.
    These zeros are not stored explicitly to save computer memory making it possible to deal with fairly high resolution ocean models.

Mathematically, the discretization and vectorization convert an expression with partial derivatives into a matrix vector product.
In summary, for the ocean circulation, we do the following conversion

$$(\mathcal{T} x)(\boldsymbol{r}) = \nabla \cdot \left[ \boldsymbol{u}(\boldsymbol{r}) - \mathbf{K}(\boldsymbol{r}) \nabla \right] x(\boldsymbol{r}) \longrightarrow \mathbf{T} \, \boldsymbol{x}$$

## Ocean circulations in AIBECS

In AIBECS, there are currently a few available circulations that you can directly load with the AIBECS:

- `Primeau_2x2x2`, the circulation described in this documentation page
- `TwoBoxModel`, the 2-box model of the *Sarmiento and Gruber* (2006) book
- `Archer_etal_2000`, the 3-box model of *Archer et al*. (2000) (also used in the *Sarmiento and Gruber* (2006) book)
- `OCIM0`, the unnamed precursor of the OCIM1 (and which should actually be named OCIM0.1)
- `OCIM1`
- `OCIM2`

To load any of these, you just need to do

```julia
grd, T = Circulation.load()
```

where `Circulation` is one of the circulations listed above.

If you are adventurous, you can create your own circulations.
AIBECS provides some tools for this.
(This is how the `Primeau_2x2x2` and `Archer_etal_2000` circulations were created.)

## Sinking particles in AIBECS

There are no loadable transport operators for sinking particles, because there are too many ways to represent sinking particles.
However, the AIBECS provides tools to create your own transport operators by providing a scalar or a vector of the downward settling velocities, `w`.
For example, you can create the transport operator for particles sinking at 100 m/d everywhere simply via

```julia
T = transportoperator(grd, w=100u"m/d")
``` 

# AIBECS functions

## Circulations

```@docs
OCIM2.load
OCIM1.load
OCIM0.load
OCCA.load
Primeau_2x2x2.load
Archer_etal_2000.load
TwoBoxModel.load
```

## Plotting

```@docs
plothorizontalslice
surfacemap
plot∫dz
plotverticalmean
minimap
plotmeridionalslice
plotzonalmean
plot∫dx
plotzonalslice
plotmeridionalmean
plot∫dy
plot∫dxdy
plothorizontalmean
plotdepthprofile
plottransect
ratioatstation
plotparameter
plotparameters
```

## Parameters

```@docs
AbstractParameters
unpack
length(<:AbstractParameters)
size(<:AbstractParameters)
values(<:AbstractParameters)
symbols
flattenable_values
flattenable_symbols
latex
table
vec
```

## For simulations

```@docs
state_function_and_Jacobian
split_state_function_and_Jacobian
SteadyStateProblem
solve
```

## For particles

```@docs
transportoperator
PFDO
DIVO
FATO
```

## For optimization

```@docs
mismatch
∇mismatch
```
