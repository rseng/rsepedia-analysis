[![DOI](https://joss.theoj.org/papers/10.21105/joss.02538/status.svg)](https://doi.org/10.21105/joss.02538)
[![linux-test Status](https://github.com/morphometry/papaya2/workflows/linux-test/badge.svg)](https://github.com/morphometry/papaya2/actions)
[![mac-test Status](https://github.com/morphometry/papaya2/workflows/mac-test/badge.svg)](https://github.com/morphometry/papaya2/actions)

# Overview

papaya2 is a small header-only C++ library for computing irreducible 2D Minkowski tensors of image and polygonal data.
The library needs a C++ 11 compliant compiler, and no external dependencies.
More detailed information can be found at <https://morphometry.org/software/papaya2/>.

If you're using this work in published research, please cite

[Schaller et al., (2020). papaya2: 2D Irreducible Minkowski Tensor computation. Journal of Open Source Software, 5(54), 2538](https://doi.org/10.21105/joss.02538)

# Installation

papaya2 needs no installation, it is an header-only library.


# Demos

papaya2 includes several demos which can be found in the demo folder,
see the [README](https://github.com/morphometry/papaya2/blob/master/demos/README.md) file.

An interactive demo can be found at <https://morphometry.org/morphometer>.
It is based on the JavaScript binding of papaya2.


# Tests

papaya2 inclues a test suite based on catch2. It can be run with

    cd test
    make


# Other language bindings

papaya2 includes bindings to Python, Matlab and JavaScript.
See READMEs in correspoding folders for more details or <https://morphometry.org/software/papaya2/>.


# Contributing

If you have some contribution to the software, please write an email to <info@morphometry.org>.
Bugs can be filed on [github](https://github.com/morphometry/papaya2/issues).
---
title: 'papaya2: 2D Irreducible Minkowski Tensor computation'
tags:
  - Minkowski Tensors
  - Morphometry
  - Image analysis
  - Voronoi diagram
  - C++
  - Python
  - Matlab
  - JavaScript
authors:
  - name: Fabian M. Schaller
    orcid: 0000-0003-2609-9988
    affiliation: "1, 2"
  - name: Jenny Wagner
    orcid: 0000-0002-4999-3838
    affiliation: "3"
  - name: Sebastian C. Kapfer
    orcid: 0000-0002-7591-2739
    affiliation: "1"
affiliations:
 - name: Theoretische Physik 1, FAU Erlangen-Nürnberg, Germany
   index: 1
 - name: Institut für Stochastik, Karlsruhe Institute for Technology, Germany
   index: 2
 - name: Zentrum für Astronomie, Universität Heidelberg, Germany
   index: 3
date: 24 May 2020
bibliography: paper.bib

---

# Summary

A common challenge in scientific and technical domains is the quantitative
description of geometries and shapes, e.g. in the analysis of microscope
imagery or astronomical observation data.  Frequently, it is desirable to
go beyond scalar shape metrics such as porosity and surface to volume ratios
because the samples are anisotropic or because direction-dependent quantities
such as conductances or elasticity are of interest.  Popular analysis software
such as [ImageJ](https://imagej.nih.gov/ij/) and [SExtractor](https://imagej.nih.gov/ij/)
provide only limited tooling for higher-order anisotropy characterization;
usually only the tensor of inertia (rank 2) is available.

Minkowski Tensors are a systematic family of versatile and robust higher-order
shape descriptors, originating in integral geometry, see @bib:AdvMatReview for an introduction and detailed references.  They
allow for shape characterization of arbitrary order and promise a path to
systematic structure-function relationships for direction-dependent properties.
Minkowski Tensors have previously been applied to data as diverse as ice grain
microstructure [@bib:SchroederMicro2010],
granular packing geometries [@bib:BeadPacksAnisotropic2010; @bib:Schaller2015],
astronomical data [@bib:Kerscher2001; @bib:Joby2019; @bib:Klatt2020],
neuronal data [@bib:Beisbart2006],
foams [@bib:Saadatfar2012; @bib:Evans2017]
and random sets, tessellations and point patterns [@bib:AnisoFluids2010; @bib:Springer2017].
An accessible introduction to Minkowski Tensors can be found on
[www.morphometry.org](https://morphometry.org/theory/anisotropy-analysis-by-imt/).

Here, we present `papaya2`, a C++ library which facilitates computation of
Irreducible Minkowski Tensors for two-dimensional geometries and shapes, including planar
objects bounded by polygonal contours, collections of points (point patterns)
and greyscale pixel data.
This library is accompanied by example programs and
bindings for Python, Matlab, and the JavaScript language.

`Papaya2` is a rewrite of [`papaya`](https://github.com/skapfer/papaya) with a
library interface, support for Irreducible Minkowski Tensors and interpolated marching squares, and
extensions to Matlab, JavaScript and Python provided.  While the tensor of inertia is computed
by many tools, we are not aware of other open-source software which provides
higher-rank shape characterization in 2D.

For the analysis of the examples in this paper, we employ our interactive online resource
[Morphometer](https://morphometry.org/morphometer/) which uses `papaya2` for its computations.

# C++ library papaya2

The C++ 11 library `papaya2` contains the core algorithms to compute Irreducible
Minkowski Tensors of two-dimensional geometries.  It processes both polygonal
and 2D image input data.

`papaya2` is a header-only template library designed to operate on user data structures.
We bundle several example programs which can be adapted to user requirements,
or employed directly for simple analyses (see section *Demos*).

The main components of the library are defined in the header file `<papaya2.hpp>`.
Analysis results are returned in a `MinkowskiAccumulator` object, which offers
accessors to retrieve common morphometric data, including the following:

- `area()`  The 2D volume (area) enclosed by the geometry

- `perimeter()`  The perimeter (boundary length) of the geometry

- `msm(s)`  The $s$-th Minkowski structure metric $q_s$,
see [Morphometry page](https://morphometry.org/theory/anisotropy-analysis-by-imt/) and @bib:Mickel2013 for details

- `imt(s)`  The $s$-th Irreducible Minkowski Tensor $\Psi_s$,
see previous item for details

The library provides convenient wrapper functions which encapsulate common analysis tasks.
In general, these functions are C++ function templates which operate on user data structures.
User-supplied data structures need to include some required methods and operators as documented in the headers.
The most important entrypoints are

- `papaya2::imt_polygon`:
compute the Irreducible Minkowski Tensors of closed simple polygons, specified as a sequence
of vertices in counterclockwise order.

- `papaya2::imt_interpolated_marching_squares`:
computes the Irreducible Minkowski Tensors of an excursion set of a single channel of a raster
graphics image (bitmap).  An extended version of the Marching Squares algorithm is
used which computes interpolated contours from 2x2 neighborhoods, see @bib:Mantz2008 for details.
The input data is passed to `papaya2` by reference via a suitable adapter class to avoid copies.
There are several examples of adapter classes provided, as well as a copying container (`BasicPhoto`).

- `papaya2::minkowski_map_interpolated_marching_squares`:
implements the Minkowski map algorithm [@bib:SchroederMicro2010] for a space-resolved anisotropy analysis.

The supplementary header `<papaya2/voronoi.hpp>` implements the Minkowski Tensor analysis of point
patterns via the Voronoi tessellation approach [@bib:AnisoFluids2010].  The demo
`ppanalysis` exemplifies how to use this header file.

# Application Examples

Here we show some examples analyzed in the [Morphometer web application](https://morphometry.org/morphometer/),
which uses `papaya2.js`, the JavaScript version of the `papaya2`.
Morphometer provides rapid analysis of small amounts of data (up to 1000 points, or 500x500 pixels).
For routine analysis we recommend using the `ppanalysis` and `imganalysis` demos or Python/Matlab bindings.

![Minkowski Tensor analysis of a polygon in Morphometer.\label{fig:morpho-ui}](morphometer-single-polygon.png)

Minkowski Tensors can be applied to different types of data:

- Single polygons: $s$-fold symmetric polygons are characterized by high values of $q_s$.
\autoref{fig:morpho-ui} shows a polygon with approximates an equilateral triangle.
Therefore, we find high values of $q_3$, $q_6$, $q_9$, etc.
The distinguished directions of each $\Psi_s$ are depicted on the right of the $q_s$ bar diagram.

![Minkowski Tensor analysis of a point pattern induced by a granular crystal cluster.\label{fig:morpho-pp-mode}](morphometer-granular-cryst-cluster.png)

- Point patterns can be, for instance, realizations of abstract point processes or data of physical particle systems.
For the Minkowski Tensor analysis, a Voronoi tessellation of the points is constructed and
Minkowski Tensors of the individual Voronoi cells are computed.
\autoref{fig:morpho-pp-mode} (left) shows a hexagonal crystal cluster surrounded by an amorphous background.
The Minkowski structure metric $q_6$ (indicated by the color of the Voronoi cells) is very well suited to detect hexagonal crystalline structures.
The presence of ideal hexagonal cells is demonstrated by the peak at $q_6 = 1$ in the histogram on the right-hand side.

![Minkowski Tensor analysis of a greyscale image: a Gaussian random field.\label{fig:morpho-image-mode}](morphometer-image-analysis.png)
 
- Greyscale images can also be analyzed in terms of Minkowski Tensors.
\autoref{fig:morpho-image-mode} (left) shows a detail of an anisotropic Gaussian random field,
which is converted into a binary image by thresholding (center) and analyzed using Minkowski Tensors (right).
The significant $q_2$ value in the Minkowski analysis (right) shows that the random field has a preferred direction,
which is also reflected by the distinguished direction marker (red color).

# Demos and language bindings

In the directory `demos`, we provide a number of example programs which use the library
for data analysis.  These are meant to be modified and adapted to user needs as required.
For simple analyses, they can be used directly, see the 
[README file](https://github.com/morphometry/papaya2/blob/master/demos/README.md) in the `demos` folder.

We also provide bindings of the library for Python, Matlab and JavaScript.

# Acknowledgements

We acknowledge funding by Deutsche Forschungsgemeinschaft as part of the [Forschergruppe GPSRS](http://gpsrs.de).

We are grateful to Daniel Hug, Günther Last, Klaus Mecke and Gerd Schröder-Turk for guidance and support,
to Michael Klatt for many discussions and for his contributions to the scientific concept and content of the
Morphometer, to Dennis Müller and Thomas Schindler for example data, and to Simon Weis for technical advice.

We use [picopng](https://lodev.org/lodepng/) for loading PNG images,
[emscripten](https://emscripten.org/) for compiling to JavaScript,
[CGAL](https://cgal.org/) for Voronoi diagrams,
and
[Catch2](https://github.com/catchorg/Catch2) for unit tests.

# References
# Demos

There are four demos explaining how to use the `papaya2` library in this directory.
They are meant to be modified and adapted to your needs.
It should be pretty easy to get them to compile on Linux/BSD/Mac systems.

More detailed information can be found at <https://morphometry.org/software/papaya2/>.

## `imganalysis`

`imganalysis` demonstrates how to analyze pixel data.  Type `make imganalysis` to build it;
it needs no external dependencies.

`imganalysis` reads images in PNG format.

    ./imganalysis in example_inputs/GRF_matern_C5.png out outdata.txt

See `./imganalysis help` and <https://morphometry.org/software/papaya2/> for further information.

## `ppanalysis`

`ppanalysis` analyzes point patterns.  For building it the CGAL library is required,
see the end of this document.  Type `make ppanalysis` to compile.

An example file can be run by executing

    ./ppanalysis in example_inputs/granular-cryst-cluster.txt out out.txt boxL 500    

See `./ppanalysis help` and <https://morphometry.org/software/papaya2/> for further information.

## `banana`

`banana` analyzes astrophysics data in FITS format.  The `CCfits` library is required to build,
see the end of this document.  After successful compilation with `make banana` it can be used like

    ./banana in example_inputs/SIE_detail2.fits out SIE_detail2_out.txt mint 0.0003 maxt 0.003 numt 100

See `./banana help` and <https://morphometry.org/software/papaya2/> for further information.

## `sersic`

`sersic` is another example analyzing pixel data.  It samples 
[Sérsic density profiles](https://en.wikipedia.org/wiki/Sersic_profile)
and computes Minkowski Tensors of their excursion sets.

    make sersic
    ./sersic scan_angle threshold 1.5 aspect 0.3 resolution 300 interpolated_marching_squares >sersic0.3out.txt

The `sersic` example does not require external dependencies.

## Installing dependencies

On Debian/Ubuntu-based systems:

    sudo apt-get install libcgal-dev libccfits-dev

On MacOS with [Homebrew](https://docs.brew.sh/)

    brew install boost cgal ccfits

Once you have installed any missing dependencies, type `make clean` to have the Makefile detect them.

If any required libraries are not found in standard paths, those paths must be added to your compiler's configuration by extending `features.mk` with the following lines:

        CXXFLAGS += -I /extra/directory/to/include -I /even/more/directories/to/include
        LDFLAGS += -L /extra/directory/which/has/the/lib

Run `make clean` to have the Makefile detect the libraries in the new paths.

### Manual download of CGAL

If your system does not come with the CGAL library, you will have to download it for compiling `ppanalysis`.
With a recent version of your C++ compiler, there is no reason to install CGAL at all as it supports a header-only mode.
The following commands will download and compile `ppanalysis` with CGAL 5.1:

        curl -L https://dl.bintray.com/boostorg/release/1.74.0/source/boost_1_74_0.tar.bz2 | tar xj
        curl -L https://github.com/CGAL/cgal/releases/download/v5.1/CGAL-5.1.tar.xz | tar xJ
        echo CXXFLAGS += -I CGAL-5.1/include -I boost_1_74_0 >features.mk
        make clean
        make ppanalysis
# Usage

Install [emscripten](https://github.com/emscripten-core/emscripten).
(papaya2.js is tested with emscripten Version 2.0.1)

To compile type `make`.

# Usage

# Installation

1. Install the header files for Python 3, numpy and CGAL (see the end of this document).
On Debian-based systems: `sudo apt-get install libcgal-dev python3-dev python3-numpy`.

2. To compile type `make`.

   If you get the error message

        Could not find the numpy C headers. Is the numpy Python module installed properly?

   you do not have the [`numpy`](https://numpy.org/) Python package installed.

3. You can put the resulting `pypaya2.so` in the same directory as your analysis
Python scripts, or on your PYTHONPATH.

# Usage

Computing the Minkowski Tensors of a polygon:

    import pypaya2
    vertices = [[2., 3.], [1., 1.], [-3., 4.], [-2., 0.], [2., 0.]
    minkval = pypaya2.imt_for_polygon(vertices)
    print(minkval['psi2'])

Performing Voronoi-Minkowski analysis of a point pattern:

    import pypaya2
    seeds = [[2., 3.], [1., 1.], [-3., 4.], [-2., 0.], [2., 0.]
    minkval = pypaya2.imt_for_pointpattern(seeds)
    print(minkval['psi2'])

Computing the Minkowski tensors of an image, assumed to be a 2D array:

    minkval = pypaya2.imt_for_image(image)
    print(minkval['psi2'])

    minkval = pypaya2.imt_for_image(image, threshold=[1, 3, 10, 30])
    print(minkval['psi2'])


## Running the tests

Run `make test`.

## Installing dependencies

On Debian/Ubuntu-based systems:

    sudo apt-get install libcgal-dev python3-dev python3-numpy

On MacOS with [Homebrew](https://docs.brew.sh/)

    brew install boost cgal numpy

**Once you have installed any missing dependencies, type `make clean` to have the Makefile detect them.**

## Compiling for Python 2

If you still have Python 2, edit the file `features.mk` and add the line `PYTHON_VERSION = 2`.
You will need `python2-dev` and `python-numpy` packages.
# Matlab examples

These extensions were tested with Matlab R2018b.  They require the data API introduced in this version.

See `imganalysis.m` for an example script that builds and uses the image analysis Matlab extension.

See `ppanalysis.m` for an example script that builds and uses the point pattern analysis Matlab extension.
