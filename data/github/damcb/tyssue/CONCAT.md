---
title: 'Tyssue: an epithelium simulation library'
tags:
  - Python
  - Developmental biology
  - Epithelium
authors:
  - name: Sophie Theis
    orcid: 0000-0003-4206-8153
    affiliation: "1, 2"
  - name: Magali Suzanne
    orcid: 0000-0003-0801-3706
    affiliation: "2"
  - name: Guillaume Gay
    orcid: 0000-0002-3957-2474
    affiliation: "1, 3"
affiliations:
 - name: Morphogénie Logiciels, 32110 St Martin d’Armagnac, France.
   index: 1
 - name: LBCMCP, Centre de Biologie Intégrative (CBI), Université de Toulouse, CNRS, UPS, Toulouse 31062, France.
   index: 2
 - name: Turing Center for Living Systems, Marseille, France
   index: 3
date: 2 december 2020
bibliography: paper.bib

---
# Summary

The `tyssue` Python library seeks to provide a unified interface to implement bio-mechanical models of living tissues. Its main focus is on vertex based epithelium models. `tyssue` allows to model the mechanical behavior of 2D, apical 3D or full 3D epithelia based on the numerical resolution of the equations of motion for the mesh vertices. Biological processes are modeled through changes in the topological and dynamical properties of the mesh. `tyssue` is a modular library. Starting with the same tissue geometry, the choice of constraints, energy potential terms and parameters increases the possibility to answer different biological questions and easily explore mechanical hypotheses.


# Statement of Need


Tissue remodeling is a complex process integrating a large number of inputs such as gene expression pattern, cell adherent properties or cell mechanics. It can be difficult to manipulate specific aspects genetically. It can even be hard to simply capture, when the process takes only few minutes. Furthermore, morphogenesis is inherently a mechanical process. To execute complex morphogenetic movements, epithelia are driven by in-plane forces, like constriction of apical cell surface [@Heer:2017], and/or out-of-plane forces, such as the apico-basal cable in apoptotic cell [@Monier:2015; @Gracia:2019] or lateral tension [@Sherrard:2010; @Sui:2018]. Modeling these processes help us to understand how tissues acquire their shape, in complement of the experimental systems, and beyond their limitations. Several vertex models have been developed in the past few years to describe the physics of epithelia (for a review, see @Alt:2017), and common features can be identified. Several kinds of models have already been published. The apical vertex model has been used several times to study topology changes during morphogenetic movement in _Drosophila, Hydra and Xenopus_ [@Staple:2010; @Farhadifar:2007; @Aegerter:2012]. Associated with protein dynamics, it has been used to study the effect of protein position on tissue organisation in zebrafish retina [@Salbreux:2012]. 3D vertex model have been used to study epithelium deformation due to normal development or to cancer development [@Okuda:2015; @Eritano:2020]. Most of the time, models are developed for a specific biological question and are difficult to adapt to other systems, for several reasons. However, there is some exception like Chaste [@Cooper:2020], which propose an open source C++ library to model cell populations or how specific events arise at the system level. With the `tyssue` library, we propose models which are adaptable and scalable with the field of research and the biological question. Topology and mechanics are implemented independently to improve the versatility of models.

The `tyssue` library defines epithelium as meshes. A vertex model defines a tissue as an assembly of vertices and edges, which can form polygonal face (in 2D) or polyhedron (in 3D). For now, we assume that cell junctions are straight lines. In `tyssue`, each edge is split, so that every face is limited by oriented "half-edges" (see figure 1 A), in a structure identical to the [Linear Cell Complex](https://doc.cgal.org/latest/Linear_cell_complex/index.html) in the CGAL library. The `tyssue` library  allows to produce different kinds of tissue, from 2D to 3D tissue (see figure 1 B). The library implements concepts and mechanisms common to all vertex models, for both topological and mechanical aspects. 


![figure 1: overview of epithelium representations in tyssue](doc/illus/Figure1.png  "Figure 1")
**Figure 1: Description of tyssue geometry.** _A-Composition of a tissue from a vertex to a cell. B-Three kind of geometry that can be used in tyssue. C-Example of cell dynamics usable in 2D and 3D._

# Features of tyssue

### Topology

Common cellular processes are implemented in our library such as cell elimination, division or rearrangements. We implemented those processes based on previous works.

Cell division is modeled as the splitting of a cell by a straight line (or plane in 3D) [@Brodland:2002], the angle and position of the division plane can be decided (see Figure 1 C, right panel).

Changes in cell neighbors - also called rearrangements - happen when the size of the boundary between two neighboring cells passes below a certain threshold length in 2D (type 1 transition), or area in 3D (I-H or H-I transition)  [@Okuda:2015]. In that case, the linked vertices fuse and are separated again, which can lead to a change in the local topology  (see Figure 1 C, left panel).

Cell elimination happens when a cell area (volume) reaches a low threshold. When this happens, cell starts to lose contact with neighboring cells through series of rearrangements. Once the cell is reduced to a triangle (in 2D) or a tetrahedron (in 3D) the remaining vertices are merged to create a new vertex.

Although it was customary to assume the neighbor exchange to be a single-step process, we follow the work by Finegan et al. which describes cell exchange as a multistep, stochastic process [@Finegan:2019]. As a consequence, in `tyssue`, vertices are not limited to 3 (in 2D) or 4 (in 3D) linked edges, but can form "rosettes" - see [type1 and rosette](https://github.com/DamCB/tyssue-demo/blob/master/05-Rearangments.ipynb) examples.


### Mechanics

In `tyssue`, the dynamical behavior of epithelium is described by solving the equation of motions following Newton's principle. At the scales of the studied processes, the inertia is negligible compared to other forces such as friction, adhesion or contraction of the actin cytoskeleton.

Honda et al. assume that cell movements respond to mechanical forces in an overdamped manner and the vertices are driven by the sum of interfacial tension on cell boundaries and the resistance force against the deformation of cells [@Honda:1978; @Honda:1983]. The `EulerSolver` class in `tyssue` allows to simulate such an over-damped movement.

Interactions in the epithelium are described as potentials depending on the mesh geometry, as described in Farhadifar et al., who showed that a 2D epithelium geometry and topology can be faithfully reproduced by finding the quasi-static equilibrium of energy depending on cell areas and junction lengths [@Farhadifar:2007]. The `QSSolver` class allows to solve this gradient descent problem.

More recently, Bi et al. focused his work on tissue rigidity which allows or not cell displacement in an epithelium, based on the relation between area and perimeter of a cell [@Bi:2015]. In `tyssue`, it is easy to define custom terms of the potential, through an object oriented model "factory" design, and use them to solve either the over-damped or gradient descent problem.

This way, it is easy to test various combinations of energy terms and find those that best fit the observed _in vivo_ dynamics.

![figure 2: Organisation of tyssue](doc/illus/Figure2.png  "Figure 2")
**Figure 2: Organisation of different part of tyssue**


Documentation of the `tyssue` Python library can be found [here](http://tyssue.io/). Notebook introduction on how to use `tyssue` library can be found [here](https://github.com/DamCB/tyssue-demo).

The `tyssue` library has already been used in several studies with different context of epithelia morphogenesis, such as leg folding and mesoderm invagination in *Drosophila melanogaster* [@Monier:2015; @Gracia:2019; @Martin:2021]. Github repository from those publications can be found [here](https://github.com/glyg/leg-joint), [here](https://github.com/suzannelab/invagination) and [here](https://github.com/suzannelab/polarity) respectively.


# Acknowledgements

The development of `tyssue` was supported by grants from the European Research Council (ERC) under the European Union Horizon 2020 research and innovation program (grant number EPAF: 648001), and from the Association Nationale de la recherche et de la Technologie (ANRT). `tyssue` has benefited from the contributions of [Hadrien Mary](https://github.com/hadim), [George Courcoubetis](https://github.com/gcourcou), [Bertrand Caré](https://github.com/bcare) and [Félix Quinton](https://github.com/felixquinton).

We wish to thank all past and present members of Magali Suzanne's team for their continuous support and for providing valuable insight on epithelium biology, Cyprien Gay for the discussions on the physics of epithelia, and the scientific python community for the core tools [@Virtanen:2020; @Harris:2020; @Hunter:2007; @mckinney:2010] we use in this project.


# Correspondence
Please contact guillaume@damcb.com

# Code

tyssue is written in Python 3. Code and detailed installation instructions can be found [here](https://github.com/DamCB/tyssue/blob/master/INSTALL.md). Continuous integration is performed with [Travis](https://travis-ci.org/github/DamCB/tyssue). The associated code coverage can be found at [CodeCov](https://codecov.io/gh/DamCB/tyssue).


# References
# What's new in 0.8.0

##

## Core

- Set an `"id"` column by default
- Add force inference algorithm in sheet class, based on method describe in Brodland et al. 2014.
- New `get_simple_index` to "quickly" get a non-oriented, full edge representation of the epithelium
- New `euler_characteristic` function to asses wether polyhedrons are closed


### Geometry
- Add two new geometry classes.

### History

- Fix drop index name in retrieve method.
- Make HistoryHdf5 record every columns by default


## Behaviors

- Removes the `sheet_events.relax` deprecated function

## Dynamics

- More coherent use of `is_active`

## Topology

- close_cell now raises an error if the cell has more than one hole
- changed the algorythm for 3D cell division to have a more robust
  way to attribute faces to the daughter cells (see `return_verts` argument in `get_division_edges`)
- do not check for condition_4 on division


## Generation

- Adds a Lloyd relaxation step for spherical sheet generation (needs to be generalized)

# What's new in 0.7.0

## Coordinate system and vertex ordering

- Added `"rx", "ry", "rz"` columns, reprensenting the source vertex positions relative to face center
- The `face_projected_pos` method has been rewritten with a much faster implementation.
- There is now a `"phi"` angle representing the angle in the
face plane.
- There is a new `order` optional parameter to `reset_index` and `sanitize` that re-indexes `edge_df` such that the vertices for each face are contiguous and ordered clockwize in the face plane.
- The `face_polygons` method benefits from this rewrite, so does the 2D `sheet_view`.


## Quasi-static solver

- Added periodic boundary conditions support

## Generation

- adds spherical sheet and monolayers from CGAL

## Install

- CGAL-less version installable with pip

## Geometry & Dynamics

- New geometry classes
- BorderElasticity, RadialTension, and  BarrierElasticity effectors
- Effectors now take nested dictionnaries as specs, instead of dictionnary of sets:
```py
    specs = {
        "face": {
            "is_alive": 1,
            "perimeter": 1.0,
            "perimeter_elasticity": 0.1,
            "prefered_perimeter": 3.81,
        }
    }
```
Instead of the previous:

```py
    specs = {
        "face": {
            "is_alive",
            "perimeter",
            "perimeter_elasticity",
            "prefered_perimeter",
        }
    }
```

This allows to give default values to the paramters directly in the effector class definition, and to have a `specs` attribute for the model directly.

A warning is raised if the effector definition uses the old style.


- Cell division on periodic boundary conditions



# What's new in 0.6.0

## Topology

Some bug fixes for 3D rearangements

## Collisions

Puts an exact kernel in c_collisions.cpp


## We switched to CodeCoV for coverage reports (purely from hype/style motivations)

## SolverIVP disapeared

In the absence of a clear way to deal with rearangement, we had to let this go for now, it may come back later...

## Behaviors

- We add two basics function in actions for sheet tissue : `increase` and `decrease`. In the near future, we will removed deprecated function that `increase` and `decrease` will replace (such as `growth`, `shrink`, `contract` and `relax`).

## History and HistoryHdf5

- new `HistoryHdf5` class that records each time point in a hdf file instead of in memory.

- new `browse_history` function that creates a widget to slide through the different time points with an ipyvolume 3D view

## Draw

- the `color` entries in edge and face specs can now be functions that take a sheet object as sole argument:
```py
specs = {
    "edge":{
        'color':'lightslategrey',
        'visible':True
    },
    "face":{
        'color': lambda sheet : sheet.face_df["apoptosis"],
        'colormap':'Reds',
        'visible':True
    }
}
```

This way, the color is updated at each function call, without having to define a new function.


## Utils

- new `get_next` function returns the indexes of the next half-edge for every edge (e.g the edge whose `srce` is the `trgt` of the current edge)


# What's new in 0.5

## Major rewrite of the rearangements

We now allow for rosettes to form, and split type1 transition in two steps: merging of edges bellow the critical length and spliting more than rank 3 vertices (or more than rank 4 in 3D). The splitting frequency is governed by two settings `p_4` and `p5p`.This follows Finegan et al 2019. See  doc/notebooks/Rosettes.ipynb for a demo.

A look a the diffs in `sheet_topology` and `bulk_topology` should convince the reader that this should result in a major increase in stability.

Automated reconnection is now treated as an event (treated by an `EventManager` instance), see `tyssue.behavior.base_events.reconnect`.

In EulerSolver, this avoids raising `TopologyChangeError` at least in explicit Euler. Topology changes in IVPSolver are not supported for now.


## Viscous solver

- added a `bounds` attribute to EulerSolver. This simply clips the displacement to avoid runaway conditions.


## Core and topology

- A new `update_rank` method allows to compute the rank of a vertex (as the number of _flat_ edges connected to it). This required to move the `connectivity` module to utils to avoid circular imports.

- We explicitly allow two sided faces to be created by `collapse_edge` or `remove_face`, they are directly eliminated.


# What's new in 0.4

##  Time dependant solvers

- Added two time-dependant solvers `solvers.viscous.EulerSolver` and `solvers.viscous.IVPSolver` -- see [the demo](doc/notebooks/SimpleTimeDependent.ipynb)

- Added a connectivity matrices module

- removed layer type 1 transition (maybe a bit harshly)

- rewrote effectors specification with default values

- added a `merge_border` function to put single edges at a 2D sheet border, and a `trim_borders` option to the `sheet.remove()` and `sheet.sanitize()` methods

## Pipes

- collision detection should now be an optional dependency
- `doc/notebooks` is now synched with the [tyssue-demo](https://github.com/damcb.tyssue-demo) repo
- new `make nbtest` uses nbval to run the notebooks under pytest
(provinding some kind of integration tests)



# What's new in 0.3.1

- Collision detection also works for the outer layers of bulk tissues, i.e. collisions of the apical or basal surfaces are avoided for a monolayer.

- Added `get_neighbors` and `get_neighborhood` method to `Epithelium` to allow local patch queries in bulk epithelia (was initially only possible for 2D sheets).


# What's new in 0.3

## Solvers

The `solvers.quasistatic.QSSolver` class provides a refactored solver that includes automatic Type 1, Type 3 and collision detection solving after each function evaluation. Use it with:

```
solver = QSSolver(with_t1=True, with_t3=True)
solver.find_energy_min(sheet, **minimize_kwargs)
```

The function signature is a bit different from the previous `sheet_vertex_solver.Solver` as key-word arguments are directly passed to scipy `minimize` function. You thus need to replace:

```python
solver_kw = {'minimize': {'method': 'L-BFGS-B',
                          'options': {'ftol': 1e-8,
                                      'gtol': 1e-8}}}
solver.find_energy_min(sheet, **solver_kw)
```

by:

```python
solver_kw = {'method': 'L-BFGS-B',
             'options': {'ftol': 1e-8,
                         'gtol': 1e-8}}}
solver.find_energy_min(sheet, **solver_kw)
```

to use the new solver.
Note that `sheet_vertex_solver.Solver` is still available.

## Behavior

###  Event management refactoring

We refactored event management with a keyword arguments only design to make passing complex parameter dictionnaries  easier.


Actions and events where added for monolayer objects.

There is now an option in the manager `append` methods kwargs to add unique event or not.

## Licence

We switched to GPL to be able to use CGAL without worrying. If this is
a problem to you, we can offer accomodations.

## Vizualisation

The use of the top level `draw.sheet_view` function is encouraged. It is now possible to specify visibility at the single face level with a `"visible"` column in the face DataFrame.


## Core

* Added a `History` class to handle time series of sheet movements

## Geometry

* Lumen volume calculation on a new geometry class (#110)
* Create a new segment vertex category : lateral in Monolayer
* adds `finally` statement to scale_unscale utils
* Change 'sagittal' key word by 'lateral' key word


## Dynamics

### New quasitatic solver class.

### New effectors

* Add LumenVolumeElasticity effector
* added SurfaceTension effector

## Bug fixes

* reset catched ValueError to Exception waiting for pandas to publish 0.24
* Better opposite management and validation for Sheet, closes #72
* Correction of color face (#85)
* fixes reset_specs warning formatting bug
* Correction of segment category for new faces create in IH transition

## Misc

The codebase now uses [black](https://github.com/ambv/black) to format all the code base.

## Pruning

* removed old isotropic model
* removes multisheet (#105)
# The `tyssue` roadmap

## Version 0.2
# tyssue : An epithelium simulation library

![A nice banner](doc/illus/banner.png)

<hr/>

[![Build Status](https://travis-ci.org/DamCB/tyssue.svg?branch=master)](https://travis-ci.org/DamCB/tyssue)

[![Coverage Status](https://coveralls.io/repos/DamCB/tyssue/badge.svg)](https://coveralls.io/r/DamCB/tyssue)

[![Doc Status](https://readthedocs.org/projects/tyssue/badge/?version=latest)](http://tyssue.readthedocs.io/en/latest/
)

[![DOI](https://zenodo.org/badge/32533164.svg)](https://zenodo.org/badge/latestdoi/32533164) [![Join the chat at https://gitter.im/DamCB/tyssue](https://badges.gitter.im/DamCB/tyssue.svg)](https://gitter.im/DamCB/tyssue?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)



| Name | Downloads | Version | Platforms |
| --- | --- | --- | --- |
| [![Conda Recipe](https://img.shields.io/badge/recipe-tyssue-green.svg)](https://anaconda.org/conda-forge/tyssue) | [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/tyssue.svg)](https://anaconda.org/conda-forge/tyssue) | [![Conda Version](https://img.shields.io/conda/vn/conda-forge/tyssue.svg)](https://anaconda.org/conda-forge/tyssue) | [![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/tyssue.svg)](https://anaconda.org/conda-forge/tyssue) |

# tyssue is now published in the Journal of Open Source Software!

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02973/status.svg)](https://doi.org/10.21105/joss.02973)



The `tyssue` library seeks to provide a unified interface to implement
bio-mechanical models of living tissues.
It's main focus is on **vertex based epithelium models**.

## Overview

### What kind of Models does it implement?

The first model implemented is the one described in
Monier et al. [monier2015apico]. It is an example of a vertex model,
where the interactions are only evaluated on the apical surface sheet
of the epithelium. The second class of models is still at an
stage. They implement a description of the tissue's rheology, within a
dissipation function formalism.

![The two models considered](doc/illus/two_models.png)

### General Structure of the modeling API

#### Design principles

> [APIs not apps](https://opensource.com/education/15/9/apis-not-apps)

Each biological question, be it in morphogenesis or cancer studies is
unique, and requires tweeking of the models developed by the
physicists. Most of the modelling softwares follow an architecture
based on a core C++ engine with a combinaison of markup or scripting
capacities to run specific simulation.

In `tyssue`, we rather try to expose an API that simplifies the
building of tissue models and running simulations, while keeping the
possibilities as open as possible.

> Separate structure, geometry and models

We seek to have a design as modular as possible, to allow the same
epithlium mesh to be fed to different physical models.

> Accessible, easy to use data structures

The core of the tyssue library rests on two structures: a set of
`pandas DataFrame` holding the tissue geometry and associated data,
and nested dictionnaries holding the model parameters, variables and
default values.

![Tyssue data structure](doc/illus/tyssue_data_management.png)

The API thus defines an `Epithelium` class. An instance of this class
is a container for the datasets and the specifications, and implements
methods to manipulate indexing of the dataframes to ease calculations.

The mesh structure is heavily inspired by
[CGAL Linear Cell Complexes](http://doc.cgal.org/latest/Linear_cell_complex/index.html),
most importantly, in the case of a 2D vertex sheet for example, each
junction edge between the cells is "splitted" between two oriented **half
edges**.


#### Creating an Epithelium

```python
## Core object
from tyssue import Sheet
## Simple 2D geometry
from tyssue import PlanarGeometry
## Visualisation (matplotlib based)
from tyssue.draw import sheet_view

sheet = Sheet.planar_sheet_2d('basic2D', nx=6, ny=7,
                              distx=1, disty=1)
PlanarGeometry.update_all(sheet)
sheet.sanitize()
fig, ax = sheet_view(sheet)
```

### Features

* Easy data manipulation.
* Multiple geometries (Sheets in 2D and 3D, monolayers, bulk).
* Easy to extend.
* 2D (matplotlib) and 3D (ipyvolume) customisable visualisation.
* Easy quasistatic model definition.
* Self collision detection. **new in 0.3**


### Documentation

* The documentation is now browsable on [tyssue.io](http://tyssue.io)
* The old documentation is still browsable online [here](http://tyssue.readthedocs.io/en/latest/)
* Introduction notebooks are available [here](https://github.com/DamCB/tyssue-demo).

### Mailing list:

tyssue@framaliste.org - https://framalistes.org/sympa/info/tyssue

Subscribe ➙ https://framalistes.org/sympa/subscribe/tyssue
Unsubscribe ➙ https://framalistes.org/sympa/sigrequest/tyssue


### Authors

* Bertrand Caré - @bcare
* Cyprien Gay - @cypriengay
* Guillaume Gay (maintainer) - @glyg
* Hadrien Mary - @hadim
* François Molino
* Magali Suzanne
* Sophie Theis - @sophietheis

## Dependencies

As all the dependencies are already completely supported in
python 3.x, **we won't be maintaining a python 2.x version**, because
it's time to move on...

### Core

- CGAL > 4.7
- Python >= 3.6
- numpy
- scipy
- matplotlib
- pandas
- pytables
- jupyter
- notebook
- quantities
- ipywidgets
- pythreejs
- ipyvolume
- vispy

### Tests

- pytest
- coverage
- pytest-cov

## Install

You can install the library with the conda package manager


```bash
conda install -c conda-forge tyssue
```


### Through PyPi

You can also install tyssue from PyPi, this is a CGAL-less version (pure python), lacking some features:

`python -m pip install --user --upgrade tyssue`

### From source

See [INSTALL.md](INSTALL.md) for a step by step install, including the necessary python environment.


## Licence

Since version 0.3, this project is distributed under the terms of the [General Public Licence](https://www.gnu.org/licenses/gpl.html).


Versions 2.4 and earlier were distributed under the [Mozilla Public Licence](https://www.mozilla.org/en-US/MPL/2.0/).

If GPL licencing is too restrictive for your intended usage, please contact the maintainer.

## Bibliography

* There is a [Bibtex file here](doc/bibliography/tyssue.bib) with collected relevant publications.

The tyssue library stemed from a refactoring of the `leg-joint` code used in [monier2015apico].


[monier2015apico]: Monier, B. et al. Apico-basal forces exerted by
  apoptotic cells drive epithelium folding. Nature 518, 245–248 (2015).

[Tamulonis2013]: Tamulonis, C. Cell-based models. (Universiteit ven Amsterdam, 2013). doi:10.1177/1745691612459060.

[Tlili2013]: Tlili,S. et al. Mechanical formalism for tissue dynamics. 6, 23 (2013).

## Research notice
Please note that this repository is participating in a study into sustainability
 of open source projects. Data will be gathered about this repository for
 approximately the next 12 months, starting from June 2021.

Data collected will include number of contributors, number of PRs, time taken to
 close/merge these PRs, and issues closed.

For more information, please visit
[our informational page](https://sustainable-open-science-and-software.github.io/) or download our [participant information sheet](https://sustainable-open-science-and-software.github.io/assets/PIS_sustainable_software.pdf).
# DamCB Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

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
reported by contacting the project team at contact@damcb.com . All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
Since version 0.3, tyssue depends on CGAL for collision detection, and thus a c++ compiler toolchain. It is thus advised to use conda for a simple installation procedure.

## Installing tyssue with conda

If you have a conda environment ready:
```
conda install -c conda-forge tyssue
```

This will install tyssue and all its dependencies, with the pre-compiled binary parts.

## Install tyssue using pip

This install a cross-platform, pure python version of tyssue.
Some advanced features are not available, namely:

- Collision detection
- Periodic boundary sheet generation

```sh
python -m pip install --user --upgrade tyssue
```

## Installing from source

Those are the instructions to install the package from source on a
debian-like linux distribution. If you allready have a basic
scientific python stack, use it, don't install anaconda.

### Install a C++ compiler

With an Debian like system, this is achieved by:

```bash
sudo apt install build-essential cmake g++
```

### Download and install `tyssue` from source

If you want to do that, I assume you allready know how to manage
dependencies on your platform. The simplest way to manage dependencies is to use [`conda`](https://docs.conda.io/en/latest/miniconda.html) to manage the dependencies (you can use [`mamba`](https://github.com/mamba-org/mamba) as a faster alternative to conda).

Start by cloning tyssue recursively to also grab pybind11:

```bash
git clone --recursive https://github.com/damcb/tyssue.git
cd tyssue
```

Then create a virtual environement:

```bash
conda env create -f environment.yml
```

Then install python:
```
python setup.py install
```



If all went well, you have successfully installed tyssue.

### Install testing utilities

```sh
pip install pytest pytest-cov nbval
```

A `Makefile` provides some utility function. Try :

```sh
make tests  # Run tests with nose
make coverage  # Run tests with coverage
make flake8  # Check PEP8 on the code
make nbtest # Tests all  the demo notebooks - requires nbval
```


### Building the documentation

The documentation uses
[nbsphinx](http://nbsphinx.readthedocs.io/en/0.2.9/index.html) to
convert the jupyter notebooks in doc/notebooks to html with sphinx.


```sh
pip install sphinx nbsphinx sphinx-autobuild
cd tyssue/doc
make html
```
# Welcome to Tyssue's documentation!


The `tyssue` library seeks to provide a unified interface to implement
bio-mechanical models of living tissues.
It's main focus is on **vertex based epithelium models**.
