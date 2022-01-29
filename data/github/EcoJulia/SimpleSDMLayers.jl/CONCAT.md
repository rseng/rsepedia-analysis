## Simple Layers for Species Distributions Modelling

**What is this package**? `SimpleSDMLayers` offers a series of types, methods,
and additional helper functions to build species distribution models. It does
*not* implement any species distribution models, although there are a few
examples of how this can be done in the documentation.

**Who is developping this package**? This package is primarily maintained by the
[Quantitative & Computational Ecology][qce] group at Université de Montréal, and
is part of the broader EcoJulia organisation.

[qce]: https://poisotlab.io/

**Can I sponsor this project**? Sure! There is a link in the sidebar on the right.
Any money raised this way will go towards the snacks and coffee fund for students,
or any charitable cause we like to support.

**How can I cite this package**? This repository itself can be cited through its
Zenodo archive ([`4902317`][zendoi]; this will generate a DOI for every
release), and there is a manuscript in *Journal of Open Science Software*
describing the package as well ([`10.21105/joss.02872`][jossdoi]).

[zendoi]: https://zenodo.org/record/4902317
[jossdoi]: https://joss.theoj.org/papers/10.21105/joss.02872

**Is there a manual to help with the package**? Yes. You can read the
documentation for the current [stable release][stable], which includes help on
the functions, as well as a series of tutorials and vignettes ranging from
simple analyses to full-fledged mini-studies.

[stable]: https://ecojulia.github.io/SimpleSDMLayers.jl/stable/

**Don't you have some swanky badges to display**? We do. They are listed at the
very end of this README.

**Can I contribute to this project**? Absolutely. The most immediate way to
contribute is to *use* the package, see what breaks, or where the documentation
is incomplete, and [open an issue]. If you have a more general question, you can
also [start a discussion]. Please read the [Code of Conduct][CoC] and the
[contributing guidelines][contr].

[CoC]: https://github.com/EcoJulia/SimpleSDMLayers.jl/blob/master/CODE_OF_CONDUCT.md
[contr]: https://github.com/EcoJulia/SimpleSDMLayers.jl/blob/master/CONTRIBUTING.md
[open an issue]: https://github.com/EcoJulia/SimpleSDMLayers.jl/issues
[start a discussion]: https://github.com/EcoJulia/SimpleSDMLayers.jl/discussions

**How do I install the package**? The latest tagged released can be installed
just like any Julia package: `]add SimpleSDMLayers`. To get the most of it, we
strongly suggest to also add `StatsPlots` and `GBIF`.

**Why are there no code examples in this README**? In short, because keeping the
code in the README up to date with what the package actually does is tedious;
the documentation is built around many case studies, with richer text, and with
a more narrative style. This is where you will find the code examples and the
figures you are looking for!

---


![GitHub](https://img.shields.io/github/license/EcoJulia/SimpleSDMLayers.jl)

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) ![GitHub contributors](https://img.shields.io/github/contributors/EcoJulia/SimpleSDMLayers.jl) ![GitHub commit activity](https://img.shields.io/github/commit-activity/m/EcoJulia/SimpleSDMLayers.jl)

![GitHub last commit](https://img.shields.io/github/last-commit/EcoJulia/SimpleSDMLayers.jl) ![GitHub issues](https://img.shields.io/github/issues-raw/EcoJulia/SimpleSDMLayers.jl) ![GitHub pull requests](https://img.shields.io/github/issues-pr-raw/EcoJulia/SimpleSDMLayers.jl)

![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/EcoJulia/SimpleSDMLayers.jl?sort=semver) ![GitHub Release Date](https://img.shields.io/github/release-date/EcoJulia/SimpleSDMLayers.jl)

![GitHub Workflow Status](https://img.shields.io/github/workflow/status/EcoJulia/SimpleSDMLayers.jl/CI?label=CI%20workflow) ![Codecov](https://img.shields.io/codecov/c/github/EcoJulia/SimpleSDMLayers.jl) ![GitHub Workflow Status](https://img.shields.io/github/workflow/status/EcoJulia/SimpleSDMLayers.jl/Documentation?label=Documentation%20workflow)
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at timothee.poisot@umontreal.ca. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Welcome!

We are looking forward to changes and improvement to this package. Before you
contribute, please read the [Code of Conduct][CoC].

[CoC]: https://github.com/EcoJulia/SimpleSDMLayers.jl/blob/master/CODE_OF_CONDUCT.md

## Repository structure

- `src` contains all the functions, types, and methods
- `tests` contains the unit tests and some integration tests

## Don't know where to start?

- report a *bug* or suggest an *improvement* -- open an [issue] on *GitHub*
- write a *vignette* -- used `SimpleSDMLayers.jl` to do something? You can add it to `docs/src/examples`!
- improve the *documentation* -- all functions have a `docstring` where they are declared, and improving them is a great way to get started

[issue]: https://github.com/EcoJulia/SimpleSDMLayers.jl/issues

## Setting up your environment

Have a look at the current [Julia documentation][pkgdoc].

[pkgdoc]: https://docs.julialang.org/en/stable/manual/packages/#Making-changes-to-an-existing-package-1

## Workflow

This section describes the general steps to make sure that your contribution is
integrated rapidly. The general workflow is as follows:

1. Before you do anything, open an issue -- this is where we discuss the potential changes to the package
1. Fork the repository (see *Branches, etc.* below)
2. Create an *explicitly named branch* from `next` (if present) or `main`
3. Create a pull request *as soon as you make the first commit*
4. Be as explicit as possible on your goals
5. Do not squash / rebase commits while you work -- we will do so when merging

### Pull requests

Creating a pull request *before* you push any code will signal that you are
interested in contributing to the project. Once this is done, push often, and be
explicit about what the commits do (see commits, below). This gives the
opportunity for feedback during your work, and allow for tweaks in what you are
doing.

A *good* pull request (in addition to satisfying to all of the criteria below)
is:

1. Single purpose - it should do one thing, and one thing only
2. Short - it should ideally involve less than 250 lines of code
3. Limited in scope - it should ideally not span more than a handful of files
4. Well tested and well documented
5. Written in a style similar to the rest of the codebase

This will ensure that your contribution is rapidly reviewed and evaluated.

### Branches, etc.

The *tagged* versions of anything on `main` are stable releases. The `main`
branch itself is the latest version, but it *must* always work (after the first
tagged release). For more intensive development, use the `next` branch, or
feature-specific branches. All significant branches are under continuous
integration *and* code coverage analysis.

### Versioning

We use [semantic versioning][sv] (`major`.`minor`.`patch`). Anything that adds
no new feature should increase the `patch` number, new non-API-breaking changes
should increase `minor`, and major changes should increase `major`. Any increase
of a higher level resets the number after it (*e.g*, `0.3.1` becomes `1.0.0` for
the first major release). It is highly recommended that you do not start working
on API-breaking changes without having received a go from us first.

[sv]: http://semver.org/
**What the pull request does**   

Explain in a few words what the pull request does.

**Type of change**   

Please indicate the relevant option(s)

- [ ] :bug: Bug fix (non-breaking change which fixes an issue)
- [ ] :sparkle: New feature (non-breaking change which adds functionality)
- [ ] :boom: Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] :book: This change requires a documentation update

**Checklist**

- [ ] The changes are documented
  - [ ] The docstrings of the different functions describe the arguments, purpose, and behavior 
  - [ ] There are examples in the documentation website
- [ ] The changes are tested
- [ ] The changes **do not** modify the behavior of the previously existing functions
  - If they **do**, please explain why and how in the introduction paragraph
- [ ] For **new contributors** - my name and information are added to `.zenodo.json`
- [ ] The `Project.toml` field `version` has been updated
  - Change the *last* number for a `v0.0.x` series release, a bug fix, or a performance improvement
  - Change the *middle* number for additional features that *do not* break the current test suite (unless you fix a bug in the test suite)
  - Change the *first* number for changes that break the current test suite
# Simple SDM Layers in *Julia*

The `SimpleSDMLayers` provides an interface to facilitate the manipulation of
raster data for species distributions modeling in *Julia*.

The two core types of the package are `SimpleSDMPredictor` and
`SimpleSDMResponse`. The only difference between the two is that predictors are
immutable, but responses are. All types belong to the abstract `SimpleSDMLayer`,
and are organised in the same way: a `grid` field storing a matrix of data (of
any type!), and the `left`, `right`, `bottom`, and `top` coordinates (as
floating point values). Of course these details are largely irrelevant, since we
have overloaded a large number of methods from `Base`, to make indexing,
converting, and modifying data as easy as possible.

The aim of the package is to deliver (i) a series of methods to manipulate
raster data, and (ii) facilitated access to common datasets used to train
species distribution models. Despite what the name may suggest, this package
does *not* implement SDMs, but is instead intended as a library usable for this
purpose. Nevertheless, the documentation contains a few example of building
models, and integrating this package with the GBIF API.

The documentations is split in three sections. The **manual** is a fairly
exhaustive documentation of the functions, methods, types, and interfaces, and
it is a good idea to skim it first, and come back for a focused reading when you
have a specific use-case in mind. The **general examples** section is a
collection of mostly disconnected workflows, intended to show how
`SimpleSDMLayers` interacts with other packages. It should give you a better
understanding of you can actually use the package.

Finally, the **SDM case studies** are a more linear series of vignettes,
covering occurrence data, variable selection, bulding a presence-only model,
generating pseudo-absences, and using a machine learning approach to do range
forecasting under climate change. This last section can be used as a template to
develop new analyses, and will use almost all the features in the package. All
of the SDM vignettes use the same species throughout - *Hypomyces lactifluorum*
is a fungus of moderate commerical importance in North America, whose
distribution is probably going to be affected by climate change.# Reading and writing files

```@docs
SimpleSDMLayers.ascii
geotiff
```# Indexing interface

In versions prior to `0.8`, the indexing syntax (`x[y]`) had been a little bit
abused to mean two different things: getting data out of a raster, and getting a
raster out of a raster. Starting from `0.8`, the indexing interface (*i.e.*
anything relying on `getindex` and `setindex!`) is used to act on values. The
resizing of rasters is now handled by the `clip` function.

## Getting values out of a raster

```@docs
getindex
```

## Writing values in a raster

```@docs
setindex!
```# Clipping and stitching rasters

```@docs
clip
```

```@docs 
Base.vcat
Base.hcat
```

```@docs
mosaic
```# Other operations

```@docs
coarsen
slidingwindow
mask
```
# Operations on values

```@docs
rescale!
rescale
replace
replace!
```

```@docs
broadcast
```# Datasets

The package offers access to bioclimatic and other datasets - they are
downloaded, saved to the disk, and then read locally. Please note that some of
them require a lot of memory, so make sure your machine can handle them.

By default, the layers are stored in the `assets` subfolder of the current
project. This being said, the prefered solution is to define a `SDMLAYERS_PATH`
environment variable pointing to a specific path, where the layers will live.
This will ensure that they are re-used between projects.

All layers are returned as `SimpleSDMPredictor`, and therefore constructed by
calling the `SimpleSDMPredictor` function on a `LayerProvider` and a
`LayerDataset`, possibly with a future climate model and scenario. In all cases,
the method accepts either a single layer, or an array of layers.

| Data provider | Dataset                | Layers | Future models    | Future scenarios                                                   |
| ------------- | ---------------------- | ------ | ---------------- | ------------------------------------------------------------------ |
| `EarthEnv`    | `Landcover`            | 12     |                  |                                                                    |
| `EarthEnv`    | `HabitatHeterogeneity` | 14     |                  |                                                                    |
| `WorldClim`   | `BioClim`              | 19     | `CMIP6`          | `SharedSocioeconomicPathway`                                       |
| `WorldClim`   | `Elevation`            | 1      | `Elevation`      |                                                                    |
| `CHELSA`      | `BioClim`              | 12     | `CMIP5`, `CMIP6` | `RepresentativeConcentrationPathway`, `SharedSocioeconomicPathway` |

## Providers and datasets

The `layernames` method (inputs are a provider and a dataset) will return a
tuple with the name of the layers.

### Data providers

```@docs
SimpleSDMLayers.LayerProvider
```

```@docs
WorldClim
CHELSA
EarthEnv
```

### Datasets

```@docs
SimpleSDMLayers.LayerDataset
```

```@docs
BioClim
LandCover
HabitatHeterogeneity
Elevation
```

## Future data

### CMIP5

```@docs
CMIP5
RepresentativeConcentrationPathway
```

### CMIP6

```@docs
CIMP6
SharedSocioeconomicPathway
```
# Types

Layers are represented by a grid, storing the content of cells in a `Matrix`,
and a bounding box indicated by the floating point coordinates of its limits.

## Implemented types

```@docs
SimpleSDMLayer
SimpleSDMResponse
SimpleSDMPredictor
```

## Getting the coordinates

```@docs
latitudes
longitudes
boundingbox
```

## Type conversion

```@docs
convert
```# Methods overloaded

To facilitate writing julian code, we have overloaded a number of methods from
`Base`. These methods should remove the need to interact with the `grid` field
directly, and also allow to set and get values using the geographic coordinates
(as opposed to the grid positions).

## From `Base`

```@docs
copy
collect
eltype
size
stride
eachindex
similar
Base.sum
Base.maximum
Base.minimum
Base.extrema
Base.max
Base.min
+
-
*
/
==
isequal
```

## From `Broadcast`

## From `Statistics`

```@docs
Statistics.mean
Statistics.median
Statistics.std
Statistics.quantile
```
