# Contributing to landlab

Thank you for contributing to Landlab! We appreciate
your help as this is largely as volunteer effort! :heart: :heart: :heart:

# How to contribute

## Reporting Bugs

Before creating a bug report, please do at least a cursory check that the
bug has not already been reported by searching the Issues portion of the
GitHub repository. If it has, add a comment to the existing issue instead of
opening a new one.

### Submitting a Bug Report

Bugs are tracked as
[GitHub issues](https://guides.github.com/features/issues/). After you've
determined you've found a new bug, please open a
[new issue](https://github.com/landlab/landlab/issues).

Explain the problem and include additional details to help maintainers
reproduce the problem. Here are some items that will make it easier
to track down the source of the problem.

*  **Use a clear and descriptive title** for the issue that identifies the
   problem.
*  **Describe the exact steps that reproduce the problem**.
*  **Provide a [minimal example](https://stackoverflow.com/help/minimal-reproducible-example)
   that demonstrates the steps** as, for example, a bash script
   along with input files. This example should reproduce your
   problem with as few lines of code as possible and easily
   reproducible my another person. Such an example almost certainly will not
   include an input file or any dependencies beyond those required by the
   `landlab_dev` conda environment.
*  **Describe the behavior you are seeing after these steps**.
*  **Describe the behavior you expect to see after these steps**.

Additionally, the answers to the following questions about your run
environment will be helpful.

*  **Which version of landlab are you using?** This could be a specific
   git sha or a release number. The best way to find this information is to
   import landlab and evaluate `landlab.__version__`
*  **What is he name and version of you OS?**
*  **What compiler are you using?**
*  **How did you build landlab (if using the development version)?**


## Submitting Changes

:tada: Whoa! This is great! We love it when folks contibute code! :tada:

Changes to landlab should be submitted as
[pull requests](http://help.github.com/pull-requests/)).

*  Create a GitHub issue that describes what you propose to do.
*  Create a topic branch that contains your changes.
*  Open a new [GitHub pull request](https://github.com/landlab/landlab/pull/new/master).
*  Ensure the pull request description clearly describes the problem
   and solution. Include the relevant issue number.

### Git Commit Messages

* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally
* For fun, consider starting the commit message with an applicable emoji:
    * :art: `:art:` when improving the format/structure of the code
    * :racehorse: `:racehorse:` when improving performance
    * :non-potable_water: `:non-potable_water:` when plugging memory leaks
    * :memo: `:memo:` when writing docs
    * :penguin: `:penguin:` when fixing something on Linux
    * :apple: `:apple:` when fixing something on macOS
    * :checkered_flag: `:checkered_flag:` when fixing something on Windows
    * :bug: `:bug:` when fixing a bug
    * :fire: `:fire:` when removing code or files
    * :green_heart: `:green_heart:` when fixing the CI build
    * :white_check_mark: `:white_check_mark:` when adding tests
    * :shirt: `:shirt:` when removing linter warnings

### Pull Request Messages

  * Rename the pull request and provide a comment that synthesizes what
    the pull request changes or adds. This helps us synthesize what
    changes have occured between Landlab releases.

## Adding new components

If you would like to create a new component, we a few conventions that we would
like you to follow.

Please visit [this part](https://landlab.readthedocs.io/en/master/development/index.html)
of the main Landlab documentation page to read about developer installation,
guidelines to contributing code, and our software development practices.

**Landlab 2 is Python >=3.6 only.**

Thanks! :heart: :heart: :heart:

The Landlab team
# Agent-based modeling with Landlab and Mesa

The tutorials in this collection illustrate how to combine Landlab components
and capabilities with the agent-based modeling (ABM) package Mesa to create
integrated simulations that combine ABMs and grid-based continuum-type models.
The examples are deliberately simple, and are designed to show how integrated
models can be built rather than to demonstrate any particular application.


---
title: 'NetworkSedimentTransporter: A Landlab component for bed material transport through river networks'
tags:
  - Python
  - Landlab
authors:
  - name: Allison M. Pfeiffer
    orcid: 0000-0002-3974-132X
    affiliation: 1

  - name: Katherine R. Barnhart
    orcid: 0000-0001-5682-455X
    affiliation: 2, 3, 4

  - name: Jonathan A. Czuba
    orcid: 0000-0002-9485-2604
    affiliation: 5

  - name: Eric W. H. Hutton
    orcid: 0000-0002-5864-6459
    affiliation: 6

affiliations:
  - name: Western Washington University, Geology Department
    index: 1
  - name: University of Colorado at Boulder, Department of Geological Sciences
    index: 2
  - name: University of Colorado at Boulder, Cooperative Institute for Research in Environmental Sciences
    index: 3
  - name: "Present affiliation: U.S. Geological Survey, Landslide Hazards Program, 1711 Illinois St., Golden, CO 80401"
    index: 4
  - name: Virginia Tech, Department of Biological Systems Engineering and The Global Change Center
    index: 5
  - name: University of Colorado at Boulder, Community Surface Dynamics Modeling System Integration Facility
    index: 6
date: 27 Apr 2020
bibliography: papers.bib
---

# Summary

Coarse sediment (sand, gravel, and cobbles) moves downstream through river networks. The transport rate of any particular sediment grain on the river bed surface is a function of both the hydraulics of that reach of river and the size distribution of the other grains in the reach. As sediment moves through a river system, grains may be deposited or eroded, burying and exposing other grains, and in the process changing the elevation and slope of each segment of river. This process of river channel evolution through the process of sediment transport is referred to as morphodynamics [@ParkerEbook]. Computational morphodynamic models allow for the prediction of sediment pulse transport, such as that which occurs after dam removal [@Cuietal2006a; @Cuietal2006b; @Cui2007a] or landsliding events [@An2017; @Benda&Dunne1997], as well as the prediction of changes in river channel bed surface grain size [@Fergusonetal2015].

Most computational morphodynamic models take an Eulerian approach, which tracks changes in bed elevation through time as a function of the spatial gradient in sediment flux [e.g., @ParkerEbook]. These models directly compute bed elevation change and sediment flux throughout the domain. One of the major drawbacks with Eulerian morphodynamic models is the difficulty in being able to 'tag' individual sediment particles to answer questions about how an individual sediment particle/input may move, when it might arrive, and what affect it will have on river morphology when it arrives downstream. To overcome this drawback and to more easily extend morphodynamic models to entire river networks, recent work has focused on developing river-network based Lagrangian sediment transport models, which track the locations of individual sediment units on a river network.  

A more comprehensive overview of river-network based sediment transport models is described by @Czubaetal2017. Of most relevance to the work described herein, @Czuba2018 introduced a network-based, Lagrangian bed material morphodynamic model that tracks the motion of individual units (referred to as “parcels”) of sediment through a river network. The model presented by @Czuba2018 has been applied to post-wildfire debris-flow sediment movement through a river network in Utah [@Murphyetal2019]. Czuba's approach improves on the existing morphodynamic models by: (1) accounting for the full river network, rather than a single longitudinal profile, (2) allowing the user to ‘tag’ particular sediment inputs and track their fate through time. Despite its advances, this existing network sediment transport model, however, has two notable drawbacks: 1) it written in a proprietary scripting language (MATLAB), and 2) it is not explicitly designed to be interoperable with other Earth-surface models, such as streamflow or landslide models.

Here, we present software that overcomes these two drawbacks, translating and expanding upon the network sediment transport model of @Czuba2018 in Landlab, a modular, Python-based package for the modeling of Earth-surface dynamics. Landlab is an open-source Python package for modeling Earth-surface processes [@Hobley2017Creative; @Barnhart2020Short]. It was designed as a modular framework, hosting a variety of process components such as flow routing, hillslope diffusion, and stream power erosion that function on a common set of landscape model grids. The ``NetworkSedimentTransporter`` is the newest of these components. We first describe computational infrastructure built in order to create the ``NetworkSedimentTransporter`` and then describe the new component itself.

The creation of the ``NetworkSedimentTransporter`` required the addition of two new data structures in the Landlab framework. First, the ``NetworkModelGrid``, which represents the model domain as connected nodes and links. Second, the ``DataRecord``, which stores a generic set of items in time and on the model grid. It is used here to store all attributes associated with the sediment parcels that move into, through, and out of the model domain.

In the ``NetworkSedimentTransporter``, sediment is represented as "parcels"-a quantity of sediment grains with common attributes such as grain diameter, lithology, and density. Each parcel is transported, buried, and eroded as a coherent unit. The river network is represented as a series of links and nodes on a ``NetworkModelGrid``. Each time the ``NetworkSedimentTransporter`` is run forward in time, the set of parcels that are in active transport is identified based on the flow conditions and bed surface grain size in each link, transport distances are calculated for all active parcels based on the @WilcockCrowe2003 equations, and parcels move through links on the network by updating their locations based on their transport distances [@Czuba2018]. As a result of parcel redistribution, the elevation of nodes and slope of the links evolves [@Czubaetal2017; @Czuba2018].

Our implementation is not a direct translation of the model implemented in MATLAB and described in @Czuba2018. Here we add three new elements to the model: sediment density that varies across parcels, downstream bed-material abrasion, and enhanced capabilities for specifying the active layer thickness.

The use of the ``DataRecord`` attributes to store density and the abrasion-rate coefficient permits different values for each sediment parcel. The density influences which parcels are mobile and how far they move each timestep. Variable (rather than constant) density permits better representing study sites with lithologic variation. Similarly, different rock types may abrade at different rates. The abrasion-rate is calculated as the loss of particle mass (or volume, because density is constant within each parcel) during transport downstream as:

$W_x = W_0 \exp \left(\alpha x \right)$

Where $x$ is the downstream transport distance, $\alpha$ is the abrasion rate (for mass loss), and $W_x$ and $W_0$ are the resulting and original sediment parcel masses, respectively. The model tracks parcel volumes (not masses) so the actual implementation replaces $W_x$ and $W_0$ with volumes (e.g., $W_0=V_0\rho_s$, where $V_0$ is the original sediment parcel volume and $\rho_s$ is the rock density of the sediment in the parcel); however, the form of the equation for mass or volume is equivalent for a parcel with a constant sediment density (i.e., the $\rho_s$ on both sides of the equation cancel out). Furthermore, once a volume reduction of each parcel is computed, the model also updates the associated reduction in parcel sediment grain size as:

$D_x = D_0 \left(\frac{V_x}{V_0}\right)^{1/3}$

Where $D_x$ and $D_0$ are the resulting and original sediment parcel diameters, respectively.

Our final modification to @Czuba2018 is enhancing the methods used for calculating variable active layer thickness. Many sediment transport models [e.g., @Cui2007b; @Czuba2018] represent the mobile portion of the grains on the riverbed at any given time as an "active layer" of constant thickness. All grains in this layer are transported, whereas all grains below this layer are immobile. Within ``NetworkSedimentTransporter`` the user has the option to specify active layer thickness as a constant value or a multiple of the mean grain size in each link. Alternatively, we incorporated the formulation of @Wongetal2007 to calculate an active layer thickness for each link in the network at each timestep as a function of Shields stress and median grain diameter.

The ``NetworkSedimentTransporter`` component of Landlab is capable of routing mixed grain size sediment through river networks to answer questions about how sediment pulses move through river networks and when, where, and how they affect downstream reaches. The accessibility of this code within the Landlab framework will make it easier for future users to modify and contribute to its continual evolution.

Source code for ``NetworkSedimentTransporter`` is available as part of the [Landlab Python package](https://github.com/landlab/landlab) and can be found in
the [``NetworkSedimentTransporter`` component](https://github.com/landlab/landlab/tree/release/landlab/components/network_sediment_transporter). The first release version of Landlab that includes the ``NetworkSedimentTransporter`` component is tagged as v2.1.0.

The Landlab project maintains a separate repository containing tutorials that introduce core concepts and the use of individual components. In addition to the
source code, a set of Jupyter Notebooks introducing the use of NetworkSedimentTransporter
are now part of the Landlab tutorials repository:
- [Part  1: Introduction with a synthetic network](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/network_sediment_transporter/network_sediment_transporter.ipynb)
- [Part  2: Using a shapefile-based river network](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/network_sediment_transporter/network_sediment_transporter_shapefile_network.ipynb)
- [Part  3: Plotting options](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/network_sediment_transporter/network_plotting_examples.ipynb)

# Acknowledgements

Barnhart supported by an NSF EAR Postdoctoral Fellowship (NSF Award Number 1725774). Czuba was partially supported by NSF-EAR (1848672), Virginia Agricultural Experiment Station, and USDA Hatch program (1017457). Pfeiffer was supported by the NCED II Synthesis Postdoctoral program and NSF-PREEVENTS (NSF Award Number 1663859 to PI Istanbulluoglu). Landlab is supported by the National Science Foundation (NSF Award Numbers 1147454, 1148305, 1450409, 1450338, and 1450412) and by the Community Surface Dynamics Modeling System (NSF Award Numbers 1226297 and 1831623). The authors thank Associate Editor Kristen Thyng, along with Zoltán Sylvester and Evan Goldstein for their thorough review of this contribution in the midst of a pandemic.

# References
---
title: 'Lithology: A Landlab submodule for spatially variable rock properties'
tags:
  - Python
  - Landlab
authors:
  - name: Katherine R. Barnhart
    orcid: 0000-0001-5682-455X
    affiliation: 1, 2
  - name: Eric Hutton
    orcid: 0000-0002-5864-6459
    affiliation: 3, 4
  - name: Nicole M. Gasparini
    orcid: 0000-0002-0803-3697
    affiliation: 5
  - name: Gregory E. Tucker
    orcid: 0000-0003-0364-5800
    affiliation: 1, 2, 3
affiliations:
  - name: University of Colorado at Boulder, Department of Geological Sciences
    index: 1
  - name: University of Colorado at Boulder, Cooperative Institute for Research in Environmental Sciences
    index: 2
  - name: University of Colorado at Boulder, Community Surface Dynamics Modeling System Integration Facility
    index: 3
  - name: University of Colorado at Boulder, Institute for Arctic and Alpine Research
    index: 4
  - name: Tulane University, Department of Earth and Environmental Sciences
    index: 5
date: 16 August 2018
bibliography: papers.bib
---

# Summary

The surface of the Earth reflects the competing advection of rock by tectonic
processes and the erosion of rock by wind, water, and ice. Rock properties
influence erosion rates by changing the processes responsible for erosion and
the rate at which rock is weathered, detached, and turned into mobile sediment.
Variations in the rock properties over space and with depth reflect the legacy
of sedimentary deposition and tectonic deformation. Long-term landscape
evolution modeling experiments that include the impact of spatially and
temporally variable rock characteristics can be used to identify the impact of
rock strength patterns on other geologic observables such as topography, erosion
rates, and detrital mineral records [e.g., @Forte2016Complexites;
@Perne2017Steady]. Identifying these relationships allows for better
interpretations of the geologic record.

Landlab is an Open Source Python package that provides a framework for the
development of 2D numerical models, typically in Earth surface dynamics
[@Hobley2017Creative]. Landlab was designed as a  modular framework in which
different process components can be mixed and matched to construct a model based
on a user's needs. Prior work on spatially variable lithology in landscape
evolution has been done using a modified version of the channel-hillslope
integrated landscape development [CHILD, @tucker2001channel] model [e.g.,
@Forte2016Complexites] and the FastScape V5 model [@braun2013very;
@Perne2017Steady]. To provide these capabilities within the Landlab framework,
there is a need for a Landlab submodule that can treat spatial variations in
rock materials.  

This contribution describes ``Lithology``, a Landlab submodule designed to
support the representation of 3D variations in rock material properties within
the Landlab framework. It includes two classes: ``Lithology`` is a generic
representation of spatially varying rock material, and ``LithoLayers`` is a
derived class that treats parallel layers of material variations. In both
classes, each rock type may have multiple attributes. Rock layers may be removed
through erosion, or added to through deposition. Two options for the underlying
datastructure are supported: "event layers", in which the data structure stores
each time-step as an event even if there is no material in the layer, or
"material layers", in which entries in the datastructure represent contiguous
material of the same property, but not necessarily the same age. This second
option is more memory efficient but does not record the transient dynamics of
erosion and deposition.

Source code for ``Lithology`` and ``Litholayers`` is available as part of the
[Landlab python package](https://github.com/landlab/landlab) and can be found in
the [``Lithology``
submodule](https://github.com/landlab/landlab/tree/release/landlab/components/lithology).
The ``Lithology`` submodule is documented using Docstrings, and the
documentation can be found on the Landlab ReadTheDocs site. One page exists for
the [Lithology
component](https://landlab.readthedocs.io/en/release/landlab.components.lithology.html)
and a second for the [LithoLayers
component](https://landlab.readthedocs.io/en/release/landlab.components.litholayers.html).
Unit and docstring tests provide 100% coverage of this submodule. [Pull Request #
674](https://github.com/landlab/landlab/pull/674) brought the ``Lithology``
submodule into the core Landlab source code. The first release version of
Landlab that includes the ``Lithology`` submodule is tagged as v1.5.4. The
concept DOI for Landlab is archived in Zenodo with the linked DOI: [@ZenodoLithologySourceCode]
and the archive for this manuscript points to the Zenodo archive of v1.5.4.

The Landlab project maintains a separate repository containing tutorials that
introduce core concepts and the use of individual submodules. In addition to the
source code, a [Jupyter Notebook introducing the use of Lithology and
Litholayers](https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/lithology/lithology_and_litholayers.ipynb)
is now part of the Landlab tutorials repository. This tutorial was brought into
the repository with [Pull Request #
19](https://github.com/landlab/tutorials/pull/19). The first release version of
the Landlab tutorials that includes this notebooks is tagged as v1.5.4 and is
archived in Zenodo with the linked DOI: [@ZenodoLithologyNotebook].

# Acknowledgements

The authors thank Adam Forte, Matt Rossi, and Brian Yanites for helpful
discussions during the development of this code and the accompanying Jupyter
notebooks.

# References
---
title: 'GroundwaterDupuitPercolator: A Landlab component for groundwater flow'
tags:
  - Python
  - Landlab
authors:
  - name: David G. Litwin
    orcid: 0000-0002-8097-4029
    affiliation: 1
  - name: Gregory E. Tucker
    orcid: 0000-0003-0364-5800
    affiliation: 2, 3, 4
  - name: Katherine R. Barnhart
    orcid: 0000-0001-5682-455X
    affiliation: 2, 3
  - name: Ciaran J. Harman
    orcid: 0000-0002-3185-002X
    affiliation: 1, 5
affiliations:
  - name: Johns Hopkins University, Department of Environmental Health and Engineering
  - index: 1
  - name: University of Colorado at Boulder, Department of Geological Sciences
    index: 2
  - name: University of Colorado at Boulder, Cooperative Institute for Research in Environmental Sciences
    index: 3
  - name: University of Colorado at Boulder, Community Surface Dynamics Modeling System Integration Facility
    index: 4
  - name: Johns Hopkins University, Department of Earth and Planetary Science
  - index: 5  
date: 18 November 2019
bibliography: papers.bib
---

# Summary
A large portion of the water that enters a catchment as precipitation percolates through soil and rock before exiting to water bodies or returning to the atmosphere as evapotranspiration. In many places, the discharge of water stored in the subsurface is a primary source of streamflow, and thus controls the ways in which catchments respond to stochastic variations in precipitation and climate [@beck_global_2013; @jasechko_pronounced_2014]. Previous studies have shown the importance of groundwater for a diverse range of processes, from transpiration to solute export [@maxwell_connections_2016; @van_verseveld_role_2009], and across diverse timescales, from rainfall-runoff response to landscape evolution [@huang_modelling_2006; @sklash_role_1979]. Of particular relevance to landscape evolution, groundwater can be an important control on the occurrence of overland flow, as the interaction of the water table with the ground surface controls the spatial extent of saturation and groundwater return flow [@dunne_partial_1970].
Variably-saturated groundwater flow is often assumed to be governed by the Richards equation, which describes how water content and/or total energy potential evolve in an idealized porous medium due to fluxes of water driven by gradients in total potential. Numerical solutions to the Richards equation are computationally expensive [e.g. @kirkland_algorithms_1992], often limiting their applications. For computational efficiency, we use the widely applied Dupuit-Forcheimer approximation, which simplifies the Richards equation when aquifers are laterally extensive in comparison to their thickness, and the capillary fringe above the water table is relatively thin [e.g. @childs_drainage_1971; @troch_hillslope-storage_2003]. In this case, the water table is modeled as a free surface with groundwater flow driven by gradients in water table elevation. This formulation is known as the Boussinesq model of an unconfined aquifer. When the model assumptions are valid, this greatly reduces the model complexity while still producing water table elevations and discharges comparable to Richards equation solutions [@hilberts_hillslope-storage_2004].

Groundwater models of varying complexity are available for different purposes. Fully coupled groundwater and surface water models such as PARFLOW [@kollet_integrated_2006], CATHY [@camporese_surface-subsurface_2010], and PIHM [@qu_semidiscrete_2007] solve the three-dimensional Richards equation for variably saturated flow, and couple this with precipitation and runoff components. MODFLOW [@langevin_modflow_2019] also solves the three dimensional Richards equation, and may be coupled to precipitation and runoff models, as in GSFLOW [@regan_gsflow_2018]. More parsimonious models (with fewer necessary parameters and calculations) are also available, such as those that solve the hillslope storage Boussinesq equation for one-dimensional groundwater flow in hillslopes of non-constant width [@marcais_dynamic_2017, @broda_low-dimensional_2012]. Two-dimensional implementations of the Boussinesq model are common, but do not appear to be widely available as open source packages. The simplicity and computational efficiency of this method is advantageous for capturing the first-order effects of groundwater flow on other Earth surface processes. More sophisticated groundwater models may be necessary depending on the hydrological features that the user intends to capture.

Implementations of the Boussinesq model often encounter numerical instabilities where the water table intersects the surface and groundwater return flow occurs through a seepage face. This is due to the presence of a discontinuity in the energy gradient from inside the hillslope (where it is determined by the water table) to the seepage face (where it is determined by the topography). @marcais_dynamic_2017 introduced a regularization that smooths the transition between subsurface flow and surface flow. Although introduced for numerical stability and not based on physical principles, this smoothing may reproduce the effect of subgrid heterogeneity where saturation within a grid cell is unlikely to be homogenous or well-reproduced by a binary condition (saturated vs unsaturated).

The ``GroundwaterDupuitPercolator`` solves the governing groundwater flow equations with an explicit, forward in time finite volume method, using the @marcais_dynamic_2017 regularization at seepage faces. While the explicit method limits the maximum timestep that can be used without jeopardizing stability, the model includes an adaptive timestep solver that subdivides the user-provided timestep in order to satisfy a Courant–Friedrichs–Lewy stability criterion.

The ``GroundwaterDupuitPercolator`` can be implemented on both regular (e.g. rectangular and hexagonal) and irregular grids determined by the user. Recharge, hydraulic conductivity, and porosity may be specified as single values uniform over the model domain, or as vectors on the nodes (recharge, porosity) or links (hydraulic conductivity) of the grid. Link hydraulic conductivity can also be specified from a two-dimensional hydraulic conductivity tensor using an included function. For mass balance calculations, the model includes methods to determine the total groundwater storage on the grid domain, the total recharge flux in, and total groundwater and surface water fluxes leaving through the boundaries.

The ``GroundwaterDupuitPercolator`` is implemented in Landlab, a Python-based open source Earth surface modeling toolkit [@hobley_creative_2017]. Landlab has a modular framework, which allows for easy coupling of different process components to meet the needs of the modeler. For example, the surface water flux from the ``GroundwaterDupuitPercolator`` can be passed to the ``FlowAccumulator`` module to route overland flow and calculate discharge at nodes. A summary of links to the documentation and example Jupyter notebooks is provided by the submodule [README](https://github.com/landlab/landlab/tree/master/landlab/components/groundwater). A diverse array of components are available, yielding many possibilities for model coupling that have not yet been explored. Given the importance of groundwater for many Earth surface processes, this component is an important contribution to the Landlab environment.  


# Acknowledgements

D. Litwin was supported in part by a Horton Research Award from the American Geophysical Union. K. Barnhart was supported by NSF EAR Postdoctoral Fellowship (NSF 1725774).

# References
---
title: 'SpeciesEvolver: A Landlab component to evolve life in simulated landscapes'
tags:
  - Python
  - Jupyter Notebook
  - Landlab
  - Landscape evolution
  - Phylogeography
  - Biodiversity
authors:
  - name: Nathan J. Lyons
    orcid: 0000-0001-6965-3374
    affiliation: 1
  - name: James S. Albert
    orcid: 0000-0001-5477-1749
    affiliation: 2
  - name: Nicole M. Gasparini
    orcid: 0000-0002-0803-3697
    affiliation: 1
affiliations:
  - name: Tulane University, Department of Earth and Environmental Sciences
    index: 1
  - name: University of Louisiana Lafayette, Department of Biology
    index: 2
date: 15 January 2020
bibliography: paper.bib
---

# Introduction and Statement of Need

The surface of the Earth and its biota evolve together. Climate and tectonics ultimately drive the physical and chemical surface processes that evolve landscape structure, including the connectivity of landscape portions that facilitate or impede movement of organismal populations [@Stanley:1979; @Antonelli:2018]. Impeded organismal movement reduces gene flow among populations and genetic diversity within populations, increasing the probability of species extinction [@Bohonak:1999]. Long-term geographic separation of populations (i.e., allopatry) is a mechanism of speciation as populations genetically diverge due to reproductive isolation [@Coyne:1992]. Speciation within the same area (i.e., sympatry) can emerge as subpopulations specialize in different resources, and along environmental gradients (e.g., in surface air temperature) where a continuum of reproductive isolation develops, among other proposed mechanisms [@Dieckmann:1999; @Doebeli:2003]. These macroevolutionary processes—dispersal, extinction, and speciation—determine regional biodiversity [@Taylor:1993].

Biodiversity is influenced by changes in landscape connectivity, nevertheless linking the evolution of a landscape with its biota is complicated by limited landform and genetic preservation, disparate timescales of surface and macroevolutionary processes, and landscape heterogeneity (e.g., species-dense assemblages, lithological variability), among numerous other factors. Discovery of surface dynamics has been aided by numerical landscape evolution models [@Tucker:2010]. These models act as digital laboratories where researchers can apply current theory to simulate complex surface dynamics while quantitatively exploring new ideas. Yet to be developed and widely shared is a modeling tool that integrates macroevolutionary processes with numerical representations of surface processes. Such a tool can help research communities overcome the complications of linking the evolution of life and landscapes.

We built ``SpeciesEvolver`` to simulate biotic evolution at geologic and macroevolutionary timescales. This software is adapted from SEAMLESS (Spatially Explicit Area Model of Landscape Evolution by SimulationS) that models organismal diversification in one-dimensional space [@Albert:2016]. ``SpeciesEvolver`` operates in two-dimensional landscapes, is built for extension, and is a component of the ``Landlab`` modeling toolkit.

``Landlab`` is an open source Python package that provides tools to build numerical models of surface dynamics [@Hobley:2017]. A landscape is represented by a model grid with configurable spatial dimensions. The surface processes components of ``Landlab`` drive change in landscape attributes (e.g., topographic elevation). The use of multiple components in a model effectively couples the processes because the processes work with the same grid and its landscape attributes. Building ``SpeciesEvolver`` into ``Landlab`` allows its users to build landscape-life evolution models to examine biological evolution alongside surface processes-driven landscape evolution.

# Software Extensibility

``SpeciesEvolver`` is adaptable to various modeling approaches and taxa, defined as a group of organisms (e.g., a population, species, or broader taxonomic level). The simulated taxa are implemented as classes with methods for macroevolutionary processes. The base class, ``Taxon`` provides abstract methods and properties that can be expanded or overridden. Software users can readily subclass ``Taxon`` designed for their model, including specialty properties (e.g., body size), behaviors (e.g., dispersal as a function of least cost paths), and taxon composition (e.g., composed of individual organisms).

The built-in implementation of taxon, ``ZoneTaxon`` evaluates macroevolutionary processes using a concept of landscape connectivity analyzed with ``Zone`` objects that manage the location of taxa in the model grid. The software user creates a function that returns the total extent of zones in the grid. Individual zones are automatically identified as the adjacent grid nodes within the total zone extent. In the example in Fig. 1a, zones are created where a grid variable—a variable selected by the user—is greater than 1. Two zones are created in this example because the nodes that meet this condition cluster in two discrete areas.

Landscape connectivity of ``ZoneTaxon`` is determined by the spatiotemporal relationships of zones. At a model time step, the relationship is determined for each zone. This relationship describes the spatial intersection of zones as ‘x-to-y’ where ‘x’ and ‘y’ are the descriptive counts of zones belonging to the intersection in the prior and current time steps, respectively. Zone connectivity relationships are illustrated in Fig. 1b. The relationships and their impact on ``ZoneTaxon`` objects are

* one-to-one: This relationship occurs where a zone in the prior time step either (a) precisely overlaps one zone in the current time step, (b) partially overlaps one zone in the current step, or (c) intersects one zone in the current step with a different shape and/or size. This relationship triggers only dispersal of the ``ZoneTaxon`` to the zone location updated in the current time step.
* one-to-many: A zone in the prior time step is overlapped by multiple zones in the current time step. The taxa in the prior step zone disperse across the multiple zones in the current step. Speciation is triggered following a duration set by the ``allopatric_wait_time`` parameter set for instances of ``ZoneTaxon``, as described in the documentation of this class.
* many-to-one: Multiple zones in the prior time step are overlapped by a zone in the current time step. Taxa extant across the multiple prior step zones exist within the same zone of the current step. By default, taxa in a zone do not affect each other, although predatory and resource-limited dynamics, for example, can be implemented in subclasses of ``ZoneTaxon``.
* many-to-many: Multiple zones in the prior step are overlapped by multiple zones in the current step. Zone counts in both the prior and current steps must be greater than one for the connectivity to be assigned this relationship. Zone counts in these steps do not need to be the same for the connectivity to be assigned this relationship. Taxa extant in prior step zones are relocated to current step zones. Speciation occurs following the ``allopatric_wait_time`` parameter set for the taxon.
* one-to-none: A zone in the prior step overlaps no zones in the current step. Taxa in the zone of the prior step become extinct as of the current time step.

@Lyons:2019 used the built-in taxon type, ``ZoneTaxon`` to investigate how changes in stream network connectivity impacted the diversity of simulated riverine species in this first application of ``SpeciesEvolver``. The species were populated to stream grid nodes and diversification emerged where stream connectivity changed. The flexibility of ``SpeciesEvolver`` with the growing library of surface processes in ``Landlab`` provides ample opportunities to discover links between landscapes and its biota. Links to ``SpeciesEvolver`` documentation and Jupyter Notebook tutorials are provided in the component [README](https://github.com/landlab/landlab/tree/master/landlab/components/species_evolution).

# Figures

![Schematics of zone creation and connectivity. Zone creation and connectivity types are explained in the text.](fig_zones.png)

# Acknowledgements

Reviewers Fiona Clubb and Evan Goldstein helped to improve this paper and code documentation.

# References
[![status](http://joss.theoj.org/papers/74487c5a6820fb2fe2898960ad6d2ea0/status.svg)](http://joss.theoj.org/papers/74487c5a6820fb2fe2898960ad6d2ea0)

Welcome to the README for the Lithology component submodule.

This submodule was published as a [Journal of Open Source
Software](http://joss.theoj.org) publication. Click the badge above to be
directed to the paper.

There is a [jupyter notebook in the Landlab Tutorials
repository](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/lithology/lithology_and_litholayers.ipynb)
that describes the use of this submodule.

If you have any questions, comments, issues, or bugs related to this submodule,
please [open an Issue](https://github.com/landlab/landlab/issues) so we can
respond.
## Welcome to the README for the GroundwaterDupuitPercolator component submodule.

[![status](https://joss.theoj.org/papers/6936ca6851c622de48b2c5f6cf45a7bd/status.svg)](https://joss.theoj.org/papers/6936ca6851c622de48b2c5f6cf45a7bd)

The GroundwaterDupuitPercolator is a component in Landlab for simulating shallow
subsurface flow. A [paper describing it](https://joss.theoj.org/papers/6936ca6851c622de48b2c5f6cf45a7bd)
was published in February 2020 in the Journal of Open Source Software. Here we
summarize installation, documentation, tutorials, tests, and getting help with
this component.

As this component lives within the larger Landlab package ecosystem, most of the
information below provides links into the [main Landlab documentation](https://landlab.readthedocs.io/).

### Installation
To use this component, you will need to install Landlab. Two options for
installation are available:
[a pre-packaged binary](https://landlab.readthedocs.io/en/master/install/index.html)
distributed through PyPI or conda-forge and a
[source code installation](https://landlab.readthedocs.io/en/master/development/install/index.html#developer-install).

The dependencies of the Landlab package are described [here](https://landlab.readthedocs.io/en/master/development/practices/dependencies.html).  

### Documentation
The documentation specific to this component is housed within the Landlab
documentation. There are two pages in the documentation that are most relevant
to this component:
- [The component API](https://landlab.readthedocs.io/en/master/reference/components/groundwater.html).
- [A page](https://landlab.readthedocs.io/en/master/reference/components/dupuit_theory.html#dupuit-theory)
describing the theory and numerical implementation of this component.

If you are new to Landlab and components, we recommend that you also look at the
[User Guide](https://landlab.readthedocs.io/en/master/user_guide/index.html),
in particular, the page on the [model grid](https://landlab.readthedocs.io/en/master/user_guide/grid.html), and [components](https://landlab.readthedocs.io/en/master/user_guide/components.html).

### Tutorials
There is a [Jupyter notebook in the Landlab Tutorials repository](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/groundwater/groundwater_flow.ipynb)
that describes the use of the `GroundwaterDupuitPercolator`.
The link takes you to a binder instance of this notebook. Its filepath within
the repository is `notebooks/tutorials/groundwater/groundwater_flow.ipynb`

A directory of all Landlab notebooks can be found (as a binder instance) [here](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/welcome.ipynb)

### Tests of this Component
Along with the rest of the Landlab package, this component uses
[`pytest`](https://docs.pytest.org/en/latest/)
to  discover and run its tests. General information about running the Landlab
tests can be found [here](https://landlab.readthedocs.io/en/master/development/install/test.html#testing).

If you want to run the tests locally, you will need to use a
[source code installation](https://landlab.readthedocs.io/en/master/development/install/index.html#developer-install).

### Getting Help
If you have any questions, comments, issues, or bugs related to this submodule,
please [open an Issue](https://github.com/landlab/landlab/issues/new) so we can
respond.
# NetworkSedimentTransporter: move sediment parcels in a river network.

The NetworkSedimentTransporter is a lagrangian model for sediment transport on a river network. This component is the subject of a forthcoming Journal of Open Source Software submission.

## Documentation and installation

Landlab documentation is hosted on this [ReadTheDocs page](https://landlab.readthedocs.io/en/release),
including instructions to install Landlab. NetworkSedimentTransporter is installed with
Landlab.

NetworkSedimentTransporter documentation is located [here](https://landlab.readthedocs.io/en/master/reference/components/network_sediment_transporter.html).

## NetworkSedimentTransporter tutorial

Three tutorials exist on NetworkSedimentTransporter. They are Jupyter notebooks accessible in the Landlab notebooks. The following are links to Binder instances of the notebooks. If instead you want to run the notebooks locally, you can clone the landlab repository and find them in the directory `landlab/notebooks/tutorials/network_sediment_transporter`.

- [Part  1: Introduction with a synthetic network](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/network_sediment_transporter/network_sediment_transporter.ipynb)
- [Part  2: Using a shapefile-based river network](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/network_sediment_transporter/network_sediment_transporter_shapefile_network.ipynb)
- [Part  3: Plotting options](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/network_sediment_transporter/network_plotting_examples.ipynb)

The index of all Landlab tutorials on Binder can bee found [here](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/welcome.ipynb) using Binder.

## Get or give help

[Open an Issue here](https://github.com/landlab/landlab/issues) where we can
respond to your questions, comments, issues, ideas, or any identified bugs
related to Landlab including NetworkSedimentTransporter.
# FlowDirAccPf: efficient filling and flow routing



``FlowDirAccPf`` is a ``Landlab`` component that provides an alternative and efficent approach to fill or breach DEMs, calculate flow directions and update flow accumulations. The component is restricted to structured grids and contains a wrapper for the RichDEM python package [@barnes2016parallel,@barnes2017parallel]. [``RichDEM``](https://richdem.readthedocs.io/en/latest/intro.html) is a set of hydrologic analysis tools using parallel processing to process large DEMs and calculate hydrologic properties.

#TODO

FlowDirAccPf is introduced [in this paper]()


## Documentation and installation

Landlab documentation is hosted on this [ReadTheDocs page](https://landlab.readthedocs.io/en/release),
including instructions to install Landlab. ``FlowDirAccPf`` is installed with
Landlab.

#TODO

``FlowDirAccPf`` documentation is located [here](https://landlab.readthedocs.io/en/release/reference/components/FlowDirAccPf.html).

## FlowDirAccPf tutorial

A ``FlowDirAccPf`` tutorial exists in the form of a Jupyter Notebooks accessible
through the following links:

#TODO adjust links when in main branch. [This](https://github.com/BCampforts/landlab/blob/bc/priority_flood/notebooks/tutorials/flow_direction_and_accumulation/PriorityFlood_realDEMs.ipynb) is a direct link to the notebook. 

- [Launch the tutorial](https://mybinder.org/v2/gh/BCampforts/landlab/blob/bc/priority_flood/notebooks/tutorials/PriorityFlood/PriorityFlood_realDEMs.ipynb)
as interactive notebook in your browser, with no need to install software,
launched using Binder.
- [A static version of the same tutorial](https://nbviewer.jupyter.org/github/BCampforts/landlab/blob/bc/priority_flood/notebooks/tutorials/PriorityFlood/PriorityFlood_realDEMs.ipynb)
- All Landlab tutorials can be launched from [this directory](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/welcome.ipynb) using Binder.

## Get or give help

[Open an Issue here](https://github.com/landlab/landlab/issues) where we can
respond to your questions, comments, issues, ideas, or any identified bugs
related to Landlab including ``FlowDirAccPf``.
# SpeciesEvolver: evolve life in simulated landscapes

[![status](https://joss.theoj.org/papers/446f3d17d642682b234ffed2b53198f6/status.svg)](https://joss.theoj.org/papers/446f3d17d642682b234ffed2b53198f6)

Life evolves alongside landscapes by biotic and abiotic processes under complex
dynamics at Earth's surface. Researchers who wish to explore these dynamics can
use this component as a tool for them to build landscape-life evolution models.
Landlab components, including SpeciesEvolver are designed to work with a shared
model grid. Researchers can build novel models using plug-and-play surface
process components to evolve the grid's landscape alongside the life tracked by
SpeciesEvolver. The simulated life evolves following customizable processes.

SpeciesEvolver is introduced [in this paper](https://doi.org/10.21105/joss.02066)
published February 2020 by the Journal of Open Source Software.

## Documentation and installation

Landlab documentation is hosted on this [ReadTheDocs page](https://landlab.readthedocs.io/en/release),
including instructions to install Landlab. SpeciesEvolver is installed with
Landlab.

SpeciesEvolver documentation is located [here](https://landlab.readthedocs.io/en/release/reference/components/species_evolution.html).

## SpeciesEvolver tutorial

A SpeciesEvolver tutorial exists in the form of a Jupyter Notebook accessible
through the following links:
- [Launch the tutorial](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/species_evolution/Introduction_to_SpeciesEvolver.ipynb)
as interactive notebook in your browser, with no need to install software,
launched using Binder.
- [A static version of the same tutorial](https://nbviewer.jupyter.org/github/landlab/landlab/blob/master/notebooks/tutorials/species_evolution/Introduction_to_SpeciesEvolver.ipynb)
- All Landlab tutorials can be launched from [this directory](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/welcome.ipynb) using Binder.

## Get or give help

[Open an Issue here](https://github.com/landlab/landlab/issues) where we can
respond to your questions, comments, issues, ideas, or any identified bugs
related to Landlab including SpeciesEvolver.
