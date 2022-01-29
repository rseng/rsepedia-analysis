# v0.16.9

Replace LightGraps with Graphs.

# v0.16.4

Add [MLJ.jl](https://github.com/alan-turing-institute/MLJ.jl) support.

# v0.16.3

New function: `midlife`.

# v0.16.0

External changes:

* `progress` keyword argument renamed to `verbose`, `field_type` keyword argument renamed to
  `field`.
* New interface: `ripserer(::Type{AbstractFiltration}, args...; kwargs...)`.
* Added `CircularCoordinates`.

Interface changes:

* All filtration constructors now have to take `verbose` as a keyword argument.
* Replaced vectors of `ChainElement`s with `Chain`s.
* Added `AbstractCell`.
* `Cube` is now an `AbstractCell`, `AbstractSimplex` is reserved for actual simplices.
* Simplices are no longer `Array`s.
* `simplex`, `unsafe_simplex`, and `unsafe_cofacet` no longer take a `sign` argument.

# v0.15.4

* Use `PersistenceDiagrams` v0.8.

# v0.15.3

* Fix type instability in zeroth interval generation.

# v0.15.2

* Update compat with Distances.jl.

# v0.15.1

* Representative cocycles are computed for infinite intervals.

# v0.15.0

* `SparseRips(...)` is deprecated. Use `Rips(...; sparse=true)`.
* `AbstractFiltration`s now need to define `births` instead of `birth`.
* Results are now sorted by persistence instead of birth time.
* Homology is now computed with the `alg=:homology` keyword argument.
* Involuted homology can be computed with the `alg=:involuted` keyword argument.
* The `reps` keyword argument can be set to a collection of integers, finding
  representatives only for specified dimensions.
* Implicit or explicit reduction can be set with the `implicit` keyword argument.
* Improved progress printing.
* New function: `find_apparent_pairs`.
<div align="center">
  <img src="https://raw.githubusercontent.com/mtsch/Ripserer.jl/master/docs/src/assets/logo-title.svg" alt="Ripserer.jl" width="480">

_Flexible and efficient persistent homology computation._

[![Coverage Status](https://coveralls.io/repos/github/mtsch/Ripserer.jl/badge.svg?branch=master)](https://coveralls.io/github/mtsch/Ripserer.jl?branch=master)
[![Build Status](https://github.com/mtsch/Ripserer.jl/workflows/Test/badge.svg)](https://github.com/mtsch/Ripserer.jl/actions?query=workflow%3ATest)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://mtsch.github.io/Ripserer.jl/dev)
[![status](https://joss.theoj.org/papers/0c8b6abead759ba068ee178fedc998a9/status.svg)](https://joss.theoj.org/papers/0c8b6abead759ba068ee178fedc998a9)

</div>

## Introduction

Ripserer is a pure Julia implementation of the [Ripser](https://github.com/Ripser/ripser)
algorithm for computing persistent homology. Its aims are to be easy to use, generic, and
fast.

See the [documentation](https://mtsch.github.io/Ripserer.jl/dev) for more information and
usage examples.

If you're looking for persistence diagram-related functionality such as Wasserstein or
bottleneck distances, persistence images, or persistence curves, please see
[PersistenceDiagrams.jl](https://github.com/mtsch/PersistenceDiagrams.jl).

## Quick start

This package is registered. To install it, run the following.

```julia
julia> using Pkg
julia> Pkg.add("Ripserer")
```

Now, generate some data.

```julia
julia> data = [(rand(), rand(), rand()) for _ in 1:200]
```

The main exported function in this package is
[`ripserer`](https://mtsch.github.io/Ripserer.jl/dev/api/ripserer/#Ripserer.ripserer). By
default, it computes Vietoris-Rips persistent homology on point cloud data and distance
matrices.

```julia
julia> ripserer(data)
# 2-element Vector{PersistenceDiagramsBase.PersistenceDiagram}:
#  200-element 0-dimensional PersistenceDiagram
#  84-element 1-dimensional PersistenceDiagram
```

[Several other filtration
types](https://mtsch.github.io/Ripserer.jl/dev/api/ripserer/#Filtrations) are supported. We
tell `ripserer` to use them by passing them as the first argument.

```julia
julia> ripserer(EdgeCollapsedRips, data)
# 2-element Vector{PersistenceDiagramsBase.PersistenceDiagram}:
#  200-element 0-dimensional PersistenceDiagram
#  84-element 1-dimensional PersistenceDiagram
```

Sometimes you may want to initialize a filtration in advance.

```julia
julia> rips = EdgeCollapsedRips(data, threshold=1)
# EdgeCollapsedRips{Int64, Float64}(nv=200)
```
```julia
julia> ripserer(rips, dim_max=2)
# 3-element Vector{PersistenceDiagramsBase.PersistenceDiagram}:
#  200-element 0-dimensional PersistenceDiagram
#  84-element 1-dimensional PersistenceDiagram
#  16-element 2-dimensional PersistenceDiagram
```

Ripserer supports plotting with
[Plots.jl](https://github.com/JuliaPlots/Plots.jl). Experimental
[Makie.jl](https://github.com/JuliaPlots/Makie.jl) support is also available
[here](https://github.com/mtsch/MakieRipserer.jl).

Plotting persistence diagrams and barcodes is straightforward:

```julia
using Plots
result = ripserer(data, dim_max=2)
plot(plot(result), barcode(result)
```
![](docs/src/assets/readme-plot-1.svg)

```julia
barcode(result)
```
![](docs/src/assets/readme-plot-2.svg)
# Contributor Guide

All contribution are welcome! Please read these guidelines before starting to work on this
project. Following these guidelines will reduce friction and improve the speed at which your
code gets merged.

## Bug Reports

If you notice code that crashes, is incorrect, or is too slow, please file a bug report. The
report should be raised as a GitHub issue with a minimal working example that reproduces the
condition. The example should include any data needed. If the problem is incorrectness, then
please post the correct result along with an incorrect result.

Please include version numbers of all relevant libraries and Julia itself.

## Documentation Fixes and Improvements

If you find typos, weird wording or missing information in the docs, don't hesitate to open
a pull request or issue. When editing docs, I recommend clicking "Edit on GitHub" and using
GitHub's online editor.

Ideas for examples or tutorials are more than welcome as well!

## Code Style

This project attempts to follow the [BlueStyle](https://github.com/invenia/BlueStyle). In
addition to that, I would like to highlight the following:

* Use descriptive variable names. Exceptions can be made for variables with a short scope.
* Avoid long lines. Line length should be at most 92 characters.
* Avoid introducing new dependencies if you can.
* Don't over-optimize non-critical code, but avoid doing unnecessary work if you can.
* Open an issue before opening starting work on a complex PR so we can discuss the changes.
---
title: 'Ripserer.jl: flexible and efficient persistent homology computation in Julia'
tags:
  - Julia
  - Persistent Homology
  - Topological Data Analysis
authors:
  - name: Matija Čufar
    affiliation: 1
affiliations:
 - name: Independent Researcher
   index: 1
date: 27 August 2020
bibliography: paper.bib
---

# Introduction

Persistent homology [@edelsbrunner2008persistent] is a relatively recent computational
technique that extracts topological information from various kinds of datasets. This
topological information gives us a good overview of the global shape of the data as well as
giving us a description of its local geometry. Since its introduction, it has been used in a
diverse range of applications, including biology [@bernoff2016biological], material science
[@lee2017quantifying], signal processing [@tralie2016high], and computer vision
[@asaad2017topological]. A problem persistent homology faces is the very large size of
combinatorial structures it has to work with. Recent algorithmic advances employ various
computational shortcuts to overcome this problem.

Among the most successful implementations of persistent homology is Ripser
[@bauer2019ripser]. With its speed and low memory usage, it makes persistent homology
practical for larger datasets, even in higher dimensions. The introduction of Ripser has
spawned a whole cottage industry of extensions and wrappers. Some examples include Ripser++
[@zhang2020gpu], Lock-free Ripser [@morozov2020towards], Ripser.py [@tralie2018ripser],
Cubical Ripser [@kaji2020cubical], and Flagser [@lutgehetmann2020computing].

In the Julia [@bezanson2017julia] space, there are few persistent homology packages
available. The ones we were able to find include Eirene.jl[^1] [@henselman2016matroid],
ComputationalHomology.jl[^2], Sparips.jl[^3] [@brehm2018sparips], and
PersistentCohomology.jl[^4], of which only Eirene.jl is available through the Julia package
manager.

# Statement of Need

A significant hurdle in developing new approaches to persistent homology stems from the fact
that developing an efficient implementation of its matrix reduction algorithm is nontrivial.

To solve this problem, we introduce Ripserer.jl, a pure Julia implementation of persistent
homology based on the algorithm that powers Ripser. It provides users with an intuitive
user interface and is readily useful as a topological data analysis framework. The other
main feature Ripserer.jl provides is the ability to hook into its algorithm through an
API. This allows researchers to experiment with different approaches to persistent homology
without having to reimplement the algorithm from scratch or forking an existing repository.

![Example visualizations. The plot on the left shows the three main representative cocycles
in the data. The right plot shows the persistence diagram.](figure.png)

# Summary

Along with its companion package, PersistenceDiagrams.jl[^5], Ripserer.jl provides a
featureful environment for computing persistent homology and integrating it with the rest of
Julia's data science stack. At the time of writing, it offers the following features.

* Fast Vietoris-Rips, alpha complex, and cubical persistent homology computation.
* Representative cocycle and critical simplex computation.
* Support for coefficients in any, possibly user-defined, field.
* Convenient persistence diagram and representative cocycle visualization via Plots.jl[^6]
  recipes.
* Bottleneck and Wasserstein matching and distance computation.
* Various persistence diagram vectorization functions, implemented with persistence
  images [@adams2017persistence] and persistence curves [@chung2019persistence].
* Easy extensibility through a documented API.

Our benchmarks[^7] show that Ripserer's performance is very close to that of Ripser. It
tends to be slightly slower for dense inputs and slightly faster for very sparse inputs. In
the cubical case, we compared it to Cubical Ripser. There the performance was worse, taking
up to 3 times as long to compute some results. This is expected as Cubical Ripser is much
more specialized for its use case and even splits its code into different repositories for
different dimensions.

We have not compared performance with newer, parallel implementations such as Ripser++ or
lock-free Ripser. Judging from the benchmarks they provide, we expect them to perform much
better. Their downside, however, is that they require powerful hardware, such as GPUs or
large numbers of processors.

[^1]: https://github.com/Eetion/Eirene.jl
[^2]: https://github.com/wildart/ComputationalHomology.jl
[^3]: https://github.com/bbrehm/Sparips.jl
[^4]: https://github.com/piever/PersistentCohomology.jl
[^5]: https://github.com/mtsch/PersistenceDiagrams.jl
[^6]: https://github.com/JuliaPlots/Plots.jl
[^7]: https://mtsch.github.io/Ripserer.jl/dev/benchmarks/

# Acknowledgments

We would like to thank Žiga Virk for comments, suggestions and ideas, and Ulrich Bauer for
making the source code of Ripser freely available.


# References
![](assets/logo-title.svg)

_Flexible and efficient persistent homology computation._

Author: Matija Čufar ([@mtsch](https://github.com/mtsch/))

## Introduction

Ripserer is a pure Julia library for computing persistent homology based on the
[Ripser](https://github.com/Ripser/ripser) algorithm. Roughly speaking, persistent homology
detects the global topological and local geometric structure of data in a noise-resistant,
stable way. If you are unfamiliar with persistent homology, I recommend reading this
[excellent
introduction](https://towardsdatascience.com/persistent-homology-with-examples-1974d4b9c3d0).

Please see the [Usage Guide](@ref) for a quick introduction, and the [API](@ref) page for
detailed descriptions of Ripserer's functionality.

While this package is fully functional, it is still in development and should not be
considered stable. I try to disrupt the public interface as little as possible, but breaking
changes might still occur from time to time.

## Installation

This package is registered. To install it, simply run the following and everything should
just work.

```julia
julia> import Pkg
julia> Pkg.add("Ripserer")
```

All versions of Julia from 1.0 onward are supported, but I recommend using the latest
version of Julia for optimal performance.

## Features

Ripserer and its companion package
[PersistenceDiagrams.jl](https://github.com/mtsch/PersistenceDiagrams.jl) currently support

* Fast Vietoris-Rips and cubical, and alpha complex persistent homology computation.
* Representative cocycle, cycle, and critical simplex computation.
* Convenient persistence diagram and representative cocycle visualization via
  [Plots.jl](https://github.com/JuliaPlots/Plots.jl). Experimental
  [Makie.jl](https://github.com/JuliaPlots/Makie.jl) is also available
  [here](https://github.com/mtsch/MakieRipserer.jl).
* Bottleneck and Wasserstein matching and distance computation.
* Various persistence diagram vectorization functions, implemented with persistence images
  and persistence curves.
* Easy extensibility through a documented API.
* Integration with [MLJ.jl](https://github.com/alan-turing-institute/MLJ.jl).
* Experimental shortest representative cycle computation.
* Experimental sparse circular coordinate computation.

To access some of the features, you need to use the PersistenceDiagrams.jl package.

## Performance

Much like Ripser, Ripserer uses several computational tricks to achieve its speed. Among
others, these include an implicit simplicial complex representation and the clearing
optimization. For a more detailed overview of these optimizations, check out [Ulrich Bauer's
article on Ripser](https://arxiv.org/abs/1908.02518).

In general, the performance of Ripserer is very close to
[Ripser](https://github.com/Ripser/ripser), usually within around 30%. Ripserer's strength
performance-wise is very sparse inputs, where it can sometimes outperform Ripser. It also
computes some things Ripser skips, like the critical simplices.

Ripserer's Cubical homology is up to 3× slower than that of [Cubical
Ripser](https://github.com/CubicalRipser/), which uses a more specialized
algorithm. Ripserer is still a good choice for small 3d images and large 2d images. Unlike
Cubical Ripser, it also supports computations on images of dimensions higher than 4.

See the [Benchmarks](@ref) section for more detailed benchmarks.

## Extending

Ripserer is designed to be easily extended with new simplex or filtration types. See the
[Abstract Types and Interfaces](@ref) API section for more information.

If you have written an extension or are having trouble implementing one, please feel free to
open a pull request or an issue. You may also contact me directly.

## Contributing

All contributions are welcome, even small things like typo fixes and ideas! See the
[contribution guidelines](https://github.com/mtsch/Ripserer.jl/blob/master/CONTRIBUTING.md)
for more information.

If you used this software in a cool project, or if you have any comments, questions, or
suggestions, feel free to contact me at
[matijacufar@gmail.com](mailto:matijacufar@gmail.com).

## Citing

If you used Ripserer in your work, consider citing the [JOSS
paper](https://joss.theoj.org/papers/10.21105/joss.02614).

A bibtex entry is provided in
[CITATION.bib](https://github.com/mtsch/Ripserer.jl/blob/master/CITATION.bib).
# Acknowledgments

I would like to thank:

* [@ubauer](https://github.com/ubauer) for creating the original
  [Ripser](https://github.com/Ripser/ripser) on which this project is based.
* [@ctralie](https://github.com/ctralie) and [@sauln](https://github.com/sauln) for creating
  [ripser.py](https://github.com/scikit-tda/ripser.py/) which has been a source of
  inspiration.
* Žiga Virk, for giving ideas and helping with the theoretical side of things.

# References

Bauer, U. (2019). Ripser: efficient computation of Vietoris-Rips persistence barcodes. arXiv
preprint [arXiv:1908.02518](https://arxiv.org/abs/1908.02518).

Kaji, S., Sudo, T., & Ahara, K. (2020). Cubical Ripser: Software for computing persistent
homology of image and volume data. arXiv preprint
[arXiv:2005.12692](https://arxiv.org/pdf/2005.12692).

Wagner, H., Chen, C., & Vuçini, E. (2012). Efficient computation of persistent homology for
cubical data. In Topological methods in data analysis and visualization II
(pp. 91-106). Springer, Berlin, Heidelberg.

Chen, C., & Kerber, M. (2011, March). Persistent homology computation with a twist. In
Proceedings 27th European Workshop on Computational Geometry (Vol. 11, pp. 197-200).

De Silva, V., Morozov, D., & Vejdemo-Johansson, M. (2011). Persistent cohomology and
circular coordinates. Discrete & Computational Geometry, 45(4), 737-759.

Zomorodian, A., & Carlsson, G. (2005). Computing persistent homology. Discrete &
Computational Geometry, 33(2), 249-274.

Edelsbrunner, H. (1993, July). The union of balls and its dual shape. In Proceedings of the
ninth annual symposium on Computational geometry (pp. 218-231).

Čufar, M. & Virk, Ž. (2021). Fast computation of persistent homology representatives with
involuted persistent homology. arXiv preprint
[arxiv:2105.03629](https://arxiv.org/abs/2105.03629)
# API

## Ripserer

```@docs
ripserer
```

## Filtrations

```@docs
Rips
```

```@docs
Cubical
```

```@docs
Custom
```

```@docs
Alpha
```

```@docs
EdgeCollapsedRips
```

## Persistence Diagrams

Persistence diagrams live in a separate package,
[PersistenceDiagrams.jl](https://github.com/mtsch/PersistenceDiagrams.jl). The package is
documented in detail [here](https://mtsch.github.io/PersistenceDiagrams.jl/dev/).

If you are looking for
Wasserstein or bottleneck distances, persistence images, betti curves, landscapes, and
similar, you will need to run `using PersistenceDiagrams`.

For convenience, the following basic functionality is reexported by Ripserer:

```@docs
PersistenceDiagrams.PersistenceDiagram
```

```@docs
PersistenceDiagrams.PersistenceInterval
```

```@docs
birth(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
death(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
persistence(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
midlife
```

```@docs
representative(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
birth_simplex(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
death_simplex(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
barcode
```

## Simplices and Representatives

```@docs
Ripserer.Simplex
```

```@docs
Ripserer.Cube
```

```@docs
Ripserer.dim(::Ripserer.AbstractCell)
```

```@docs
Ripserer.birth(::Ripserer.AbstractCell)
```

```@docs
Ripserer.index(::Ripserer.AbstractCell)
```

```@docs
Ripserer.vertices(::Ripserer.AbstractCell)
```

```@docs
Ripserer.Chain
```

```@docs
Mod
```

## MLJ.jl Interface

```@docs
Ripserer.RipsPersistentHomology
```

```@docs
Ripserer.AlphaPersistentHomology
```

```@docs
Ripserer.CubicalPersistentHomology
```

## Experimental Features

```@docs
Ripserer.reconstruct_cycle
```

```@docs
Ripserer.Partition
```

```@docs
Ripserer.CircularCoordinates
```

## Abstract Types and Interfaces

```@docs
Ripserer.AbstractFiltration
```

```@docs
Ripserer.nv(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.births(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.vertices(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.edges(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.simplex_type
```

```@docs
Ripserer.simplex
```

```@docs
Ripserer.unsafe_simplex
```

```@docs
Ripserer.unsafe_cofacet
```

```@docs
Ripserer.threshold(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.columns_to_reduce
```

```@docs
Ripserer.emergent_pairs
```

```@docs
Ripserer.postprocess_diagram
```

```@docs
Ripserer.distance_matrix
```

```@docs
Ripserer.AbstractRipsFiltration
```

```@docs
Ripserer.adjacency_matrix(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.AbstractCustomFiltration
```

```@docs
Ripserer.simplex_dicts
```

```@docs
Ripserer.AbstractCell
```

```@docs
Ripserer.AbstractSimplex
```

```@docs
Base.sign(::Ripserer.AbstractCell)
```

```@docs
Base.:-(::Ripserer.AbstractCell)
```

```@docs
Ripserer.coboundary
```

```@docs
Ripserer.boundary
```
# Benchmarks

The following tables show benchmarks that compare Ripserer's performance with
[Ripser](https://github.com/Ripser/ripser), [Cubical
Ripser](https://github.com/CubicalRipser/), and
[Eirene.jl](https://github.com/Eetion/Eirene.jl). The benchmarking code and more info about
the datasets are available [here](https://github.com/mtsch/RipsererBenchmarks.jl).

All benchmarks were performed on a laptop with an Intel(R) Core(TM) i5-4200U CPU @ 1.60GHz
with 8GB of RAM.

We used [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl/) to perform the
timing benchmarks and [Valgrind's Massif
tool](https://www.valgrind.org/docs/manual/ms-manual.html) to measure peak heap sizes
(i.e. total memory footprint).

The benchmarks were performed with Ripserer v0.15, `master` versions of Ripser (commit
hash `286d369`) and Cubical Ripser (commit hashes `6edb9c5` for 2D and `a063dac` for 3D),
and Eirene v1.3.5.

The timings show the minimum time taken among five runs of the benchmark.

The heap sizes for Ripserer include the Julia runtime.

## Comparison with Ripser

In this experiment, we performed benchmarks with the datasets presented in the [Ripser
article](https://arxiv.org/abs/1908.02518). We only used the datasets that we were able to
run with less than 8GB memory. All datasets were parsed as `Float32` as that is what Ripser
supports. The time it takes to parse a file is included for both Ripser and Ripserer.

### Dense results

|dataset       |size|dim|threshold|Ripserer|Ripser   |ratio|Ripserer heap|Ripser heap|
|:-------------|:---|:--|:--------|:-------|:--------|:----|:------------|:----------|
|`o3_1024`     |1024|3  |1.8      |4.576 s |3.057 s  |1.497|374.1 MiB    |151.0 MiB  |
|`o3_4096`     |4096|3  |1.4      |151.527 s|76.177 s|1.989|4.7 GiB      |4.1 GiB    |
|`dragon2000`  |2000|1  |         |3.133 s |2.833 s  |1.106|316.7 MiB    |296.8 MiB  |
|`fract-r`     |512 |2  |         |22.807 s|19.482 s |1.171|2.2 GiB      |2.0 GiB    |
|`random16`    |50  |2  |         |8 ms    |10 ms    |0.803|111.1 MiB    |1.1 MiB    |
|`sphere_3_192`|192 |2  |         |1.549 s |1.491 s  |1.039|287.0 MiB    |209.5 MiB  |

### Sparse results

These benchmarks were performed with the `sparse=true` keyword argument.

|dataset       |size|dim|threshold|Ripserer|Ripser   |ratio|Ripserer heap|Ripser heap|
|:-------------|:---|:--|:--------|:-------|:--------|:----|:------------|:----------|
|`o3_1024`     |1024|3  |1.8      |3.036 s |3.057 s  |0.993|418.2 MiB    |151.0 MiB  |
|`o3_4096`     |4096|3  |1.4      |76.052 s|76.177 s |0.998|4.9 GiB      |4.1 GiB    |
|`dragon2000`  |2000|1  |         |3.588 s |2.833 s  |1.267|350.4 MiB    |296.8 MiB  |
|`fract-r`     |512 |2  |         |25.399 s|19.482 s |1.304|2.2 GiB      |2.0 GiB    |
|`random16`    |50  |2  |         |9 ms    |10 ms    |0.932|111.1 MiB    |1.1 MiB    |
|`sphere_3_192`|192 |2  |         |1.734 s |1.491 s  |1.163|288.5 MiB    |209.5 MiB  |

### Alpha-Rips

These benchmarks were performed on sparse matrices that correspond to the 1-skeleta of
Delaunay triangulations. The purpose of these is to show performance with very sparse
inputs.

|dataset              |size |dim|Ripserer  |Ripser    |ratio|Ripserer heap|Ripser heap|
|:--------------------|:----|:--|:---------|:---------|:----|:------------|:----------|
|`alpha_3_sphere_3000`|3000 |3  |636 ms    |789 ms    |0.807|138.4 MiB    |33.2 MiB   |
|`alpha_torus_10_000` |10000|2  |872 ms    |1.179 s   |0.741|130.0 MiB    |27.7 MiB   |
|`alpha_5_sphere_1000`|1000 |5  |49.431 s  |46.707 s  |1.058|387.2 MiB    |202.0 MiB  |
|`alpha_dragon_2000`  |2000 |2  |56 ms     |76 ms     |0.744|2.4 GiB      |1.5 GiB    |
|`alpha_4_sphere_2000`|2000 |4  |5.844 s   |6.203 s   |0.942|110.9 MiB    |33.2 MiB   |

## Comparison with Cubical Ripser

In these benchmarks, we used some of the datasets presented in the [Cubical
Ripser](https://arxiv.org/abs/2005.12692) article. We limited the 2D image size to 1999×999
as the current `master` (commit hash `6edb9c5`) version of 2D Cubical Ripser throws an
assertion error for anything larger. We were also unable to perform 3D 256×256×256 image
benchmarks due to Ripserer running out of memory. The `eltype` of all datasets is `Float64`,
because that is what Cubical Ripser supports. When running Ripserer in the real world, it's
a good idea to use the image's native data types. This will _slightly_ reduce the memory
footprint and increase performance.

|dataset       |size   |dim|Ripserer  |Cubical Ripser|ratio|Ripserer heap|Cubical Ripser heap|
|:-------------|:------|:--|:---------|:-------------|:----|:------------|:------------------|
|`lena512`     |262144 |1  |787 ms    |299 ms        |2.631|145.0 MiB    |49.3 MiB           |
|`lena1999x999`|1997001|1  |2.87 s    |2.009 s       |1.429|514.4 MiB    |186.7 MiB          |
|`bonsai64`    |262144 |2  |2.875 s   |2.996 s       |0.96 |280.6 MiB    |1.3 GiB            |
|`bonsai128`   |2097152|2  |31.151 s  |14.733 s      |2.114|1.5 GiB      |1.9 GiB            |
|`head128`     |2097152|2  |24.102 s  |12.434 s      |1.938|1.5 GiB      |1.9 GiB            |

## Comparison with Eirene

In these benchmarks, we compare Ripserer to
[Eirene.jl](https://github.com/Eetion/Eirene.jl). Ripserer benchmarks were run with
`alg=:involuted`, so this measures the time it takes to compute representative cycles.

|dataset     |size|dim|threshold|Ripserer  |Eirene  |ratio|
|:-----------|:---|:--|:--------|:---------|:-------|:----|
|`gcycle`    |100 |3  |         |6.231 s   |24.158 s|0.258|
|`hiv`       |1088|1  |         |1.824 s   |7.774 s |0.235|
|`dragon1000`|1000|1  |         |575 ms    |8.441 s |0.068|
|`celegans`  |297 |2  |         |4.217 s   |4.588 s |0.919|
|`o3_1024`   |1024|3  |1.8      |5.735 s   |8.314 s |0.69 |
|`random16`  |50  |7  |         |8.577 s   |7.688 s |1.116|
# Related Julia Packages

This section attempts to provide an overview of Julia packages implementing persistent
homology. If you're trying to do something Ripserer is missing, one of these might have what
you're looking for.

This list is incomplete and only includes things I'm aware of. The descriptions are based
on my (limited) experience with these packages.

* [ComputationalHomology.jl](https://github.com/wildart/ComputationalHomology.jl) uses a
  different, slower algorithm, but can do some things Ripserer can't such as Čech
  persistent homology and homology of CW complexes. The package is a part of
  [TDA.jl](https://github.com/wildart/TDA.jl), which also offers other topological data
  analysis tools like Mapper.

* [Eirene.jl](https://github.com/Eetion/Eirene.jl) uses a [different algorithm based on
  matroids](https://arxiv.org/abs/1606.00199). A benefit of this algorithm is that it can
  recover persistent homology generators and a bunch of other data.

* [Ripser.jl](https://github.com/mtsch/Ripser.jl) my deprecated wrapper of the original
  C++ program. Very bare-bones and outdated. Runs a bit slower than Ripserer.

* [Sparips.jl](https://github.com/mtsch/Ripser.jl) this is a preprocessor that comes with
  a different wrapper of Ripser. The preprocessor allows you to compute persistent homology
  of very large datasets. Integrating it with Ripserer is on my TODO list.

* [PersistentCohomology.jl](https://github.com/piever/PersistentCohomology.jl) does
  essentially the same thing, but with a different algorithm, that is orders of magnitude
  slower.

If you are a developer of a persistent homology package not on this list, or any of the
information here is incorrect, let me know or open a PR.
