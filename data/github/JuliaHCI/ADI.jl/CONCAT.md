# ADI.jl

[![Build Status](https://github.com/juliahci/ADI.jl/workflows/CI/badge.svg?branch=main)](https://github.com/juliahci/ADI.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/A/ADI.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliahci/ADI.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliahci/ADI.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliahci.github.io/ADI.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliahci.github.io/ADI.jl/dev)
[![JOSS](https://joss.theoj.org/papers/32605be405e024fcbd15cd81dfdf9985/status.svg)](https://joss.theoj.org/papers/32605be405e024fcbd15cd81dfdf9985)
[![DOI](https://zenodo.org/badge/250468435.svg)](https://zenodo.org/badge/latestdoi/250468435)

A package for angular differential imaging (ADI) post-processing algorithms.

## Installation

ADI.jl is a registered package and can be installed using the Julia package manager. From the Julia REPL, enter Pkg mode (by pressing `]`)

```julia
julia>]

(@v1.5) pkg> add ADI
```

To exit Pkg mode, just backspace. Once the package is installed it can be imported with

```julia
julia> using ADI
```

To exit Pkg mode, just backspace. Once the package is installed it can be imported with
For more information, see the [Pkg documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/).

## Citations

If you use ADI.jl or derivatives in your work, please consider citing both the JOSS paper and the code record. The JOSS paper citation can be found in [`CITATION.bib`](CITATION.bib). The code will have a unique reference for each released version, so visit the [Zenodo record](https://doi.org/10.5281/zenodo.3977789) to grab the BibTeX for whichever version you used.

## Usage

The following is an extremely brief PCA reduction of an ADI cube. Please see the [documentation](https://juliahci.github.io/ADI.jl/dev/) for further usage, tutorials, and api reference.

```julia
julia> using ADI

julia> cube, angles = # load data

julia> alg = PCA(ncomps=10)

julia> flat_resid = alg(cube, angles) # ADI

julia> flat_resid_rdi = alg(cube, angles; ref=cube_ref) # flexible RDI
```

get the S/N and significance

```julia
julia> fwhm = # PSF fwhm in pixels

julia> snmap = detectionmap(snr, flat_residual, fwhm)

julia> sigmap = detectionmap(significance, flat_residual, fwhm)
```

get the contrast curve

```julia
julia> psf = # load psf or choose from HCIToolbox.Kernels

julia> cc = contrast_curve(alg, cube, angles, psf; fwhm=fwhm)
```

which can be easily loaded into a `DataFrame` or anything adopting the Tables.jl interface.

```julia
julia> using DataFrames

julia> df = DataFrame(cc)

julia> first(df, 5)
```

## Contributing and Support

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

In general contributions should follow [ColPrac](https://github.com/SciML/ColPrac). If you are interested in extending/improving ADI.jl, head to the [discussions](https://github.com/JuliaHCI/ADI.jl/discussions) to reach out. For support with using ADI.jl, please open an [issue](https://github.com/JuliaHCI/ADI.jl/issues/new/) describing the problem and steps to reproduce it.

## License

This package is licensed under the MIT Expat license. See [LICENSE](LICENSE) for more information.

---

**Author's Note**: This package is still under active development and is subject to change. Anything from minor behind-the-scenes details to large-scale design can change as I incorporate more methods into ADI.jl. I don't plan on spending much time with deprecation warnings throughout this process, since that limits my ability to experiment with implementation ideas and design goals. This package follows [semantic versioning](https://semver.org/), so an upgrade from `0.6` to `0.7` may be breaking and I recommend anybody using this package to browse the release notes for changes. Once ADI.jl is somewhat stable, I'll release a version `1.0`, at which point I'll worry about deprecations and other long-term usability considerations.
---
title: 'ADI.jl: A Julia Package for High-Contrast Imaging'
tags:
  - Julia
  - astronomy
  - high-contrast imaging
  - direct imaging
  - image processing
authors:
  - name: Miles Lucas
    orcid: 0000-0001-6341-310X
    affiliation: 1
  - name: Michael Bottom
    orcid: 0000-0003-1341-5531
    affiliation: 1
affiliations:
 - name: Institute for Astronomy, University of Hawai'i
   index: 1
date: 11/02/2020
bibliography: paper.bib
---

# Summary

High-contrast imaging (HCI) is a powerful technique for discovering and characterizing exoplanets. Being able to probe the architecture, formation, and atmospheres of planets *directly* is necessary for advancing companion formation and evolution theory [@bowler_imaging_2016]. The process required to image an exoplanet is daunting, however, due to the brightness and proximity of their host stars. One part of making such a difficult detection is the image processing of HCI data to remove systematic signals and attenuate noise. The size of HCI data and complexity of processing algorithms requires an efficient numerical framework that simultaneously offers modularity for rapid experimentation and discovery.

Angular differential imaging (ADI) is an observational technique for HCI that utilizes the Earth's rotation throughout a night of observing [@liu:2004; @marois_angular_2006]. Normally telescopes have optics to counter this rotation, but disabling these optics and taking images throughout the night will give a sequence of frames where the sky appears to rotate. The telescope optics produce a systematic noise that will not rotate with the sky, although it can slowly vary over time. The sequence of frames are pre-processed to calibrate the images and fix defects like bad pixels. The pre-processed images are then co-aligned and concatenated together into a data cube. ADI algorithms exploit the difference in rotation to approximate and subtract the systematics without substantially overfitting the signals from potential companions.

ADI algorithms differ in how they estimate the noise model. For example instead of estimating the systematics from the *target* star, a *reference* star can be used (reference differential imaging, RDI). Instead of using the entire frame, the cube can be processed in annuli corresponding to different circumstellar regions. These geometric techniques are independent of the underlying description of the algorithms. Similarly, GPU programming or out-of-core processing are computational techniques which appear like implementation details in comparison to the algorithms or how they are applied. Creating a *modular* and *generic* framework for HCI enables scientists to explore different algorithms and techniques flexibly, allowing more thorough and deeper investigations into the capabilities of HCI for finding exoplanets.

# Statement of need

`ADI.jl` is a Julia framework for post-processing high-contrast imaging (HCI) data. By organizing algorithms separately from their application, `ADI.jl` offers a modular API that benefits both observers and algorithm designers. Observers can rapidly experiment with multiple post-processing algorithms and techniques to optimally reduce their data. Julia's dynamic just-in-time (JIT) LLVM compiler [@Julia-2017; @2018arXiv180803370B] means this experimentation comes at a low runtime cost to the observer, enabling broader experimentation or higher throughput, for example, in the case of survey pipelines.

Algorithm designers will find that Julia is highly composable, so extending or adding a new algorithm only requires writing the code that is *unique to that algorithm*. Julia's language interoperability also means the algorithm can be implemented in Python or C, for example. In other words, to be able to fully use the post-processing capabilities of `ADI.jl` a new algorithm only needs to implement one or two methods. Furthermore, computational techniques like GPU programming are available *generically* through packages like `CUDA.jl` [@besard2018juliagpu].

Currently `ADI.jl` supports full-frame ADI and RDI processing, with experimental support for spectral differential imaging (SDI). The algorithms that are currently implemented are median subtraction [@marois_angular_2006], principal component analysis [PCA/KLIP; @soummer_detection_2012], non-negative matrix factorization [NMF; @ren_non-negative_2018], and fixed-point greedy disk subtraction [GreeDS; @pairet_reference-less_2019; @pairet_mayonnaise_2020]. In addition, common metrics such as S/N maps and contrast curves are available for posterior analysis. Forward modeling is being built in a separate Julia package, `Firefly.jl`, as part of active research.

# Comparisons with existing software

High-contrast imaging as a field predominantly utilizes Python for data reduction. We break down some of the necessary computations into *pre-processing*, which includes raw calibration of data, the centering and stacking of the data cube, bad-pixel removal, etc., *post-processing*, which includes the PSF approximation and subtraction, *detection metrics* which includes methods for analyzing post-processed data to make detections or find limits of detections, and finally *forward modeling* which includes various statistical models for companions and disks that can be used with post-processing algorithms in a maximum likelihood framework. `ADI.jl` primarily focuses on post-processing and detection metrics.

Some notable libraries for HCI tasks include the Vortex Imaging Pipeline (`VIP`) [@gomez_gonzalez_vip_2016], `pyKLIP` [@2015ascl.soft06001W], and `PynPoint` [@pynpoint:2019]. A table of the feature sets of these packages alongside `ADI.jl` is presented in the [online documentation](https://juliahci.github.io/ADI.jl/dev/gettingstarted/#Feature-comparison). In particular, `VIP` has served as a useful source of information regarding HCI image-processing as well as detailed implementations of common ADI algorithms. This has been indispensable in the development of `ADI.jl`, although this package is not a direct translation. In our small benchmark suite, `ADI.jl` is roughly the same speed or slightly quicker than `VIP` for the algorithms we tested, except for NMF, while detection maps were ~2 orders of magnitude quicker.

In general `VIP` offers the most diversity in algorithms and their applications, but not all algorithms are as feature-complete as the PCA implementation. `VIP` also contains many useful utilities for pre-processing and a pipeline framework. `pyKLIP` primarily uses the PCA (KLIP) algorithm, but offers many forward modeling implementations. `PynPoint` has a highly modular pre-processing module that is focused on pipelines.

# Acknowledgments

We acknowledge the `VIP` team for creating a package that has been an indispensable resource in the creation of `ADI.jl`. We also acknowledge the ongoing research using `ADI.jl` for forward modeling of exoplanets with the package `Firefly.jl`.

# References
```@meta
CurrentModule = ADI
```

# ADI.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliahci/ADI.jl)
[![Build Status](https://github.com/juliahci/ADI.jl/workflows/CI/badge.svg?branch=main)](https://github.com/juliahci/ADI.jl/actions)
[![Coverage](https://codecov.io/gh/juliahci/ADI.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliahci/ADI.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![JOSS](https://joss.theoj.org/papers/32605be405e024fcbd15cd81dfdf9985/status.svg)](https://joss.theoj.org/papers/32605be405e024fcbd15cd81dfdf9985)
[![DOI](https://zenodo.org/badge/250468435.svg)](https://zenodo.org/badge/latestdoi/250468435)

A package for angular differential imaging (ADI) along with its variants, such as reference differential imaging (RDI) and spectral differential imaging (SDI).

## Installation and Setup

ADI.jl is a registered package and can be installed using the Julia package manager. From the Julia REPL, enter Pkg mode (by pressing `]`)

```julia
julia>]

(@v1.5) pkg> add ADI
```

To exit Pkg mode, just backspace. Once the package is installed it can be imported with

```julia
julia> using ADI
```

For more information, see the [Pkg documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/).

## Citations

If you use ADI.jl or derivatives in your work, please consider citing both the JOSS paper and the code record. The JOSS paper citation can be found in [`CITATION.bib`](https://github.com/juliahci/ADI.jl/blob/master/CITATION.bib). The code will have a unique reference for each released version, so visit the [Zenodo record](https://doi.org/10.5281/zenodo.3977789) to grab the BibTeX for whichever version you used.

## Contributing and Support

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

In general contributions should follow [ColPrac](https://github.com/SciML/ColPrac). If you are interested in extending/improving ADI.jl, head to the [discussions](https://github.com/JuliaHCI/ADI.jl/discussions) to reach out. For support with using ADI.jl, please open an [issue](https://github.com/JuliaHCI/ADI.jl/issues/new/) describing the problem and steps to reproduce it.

## License

This work is distributed under the MIT "expat" license. See [`LICENSE`](https://github.com/juliahci/ADI.jl/blob/main/LICENSE) for more information.
# Index

```@index
```
# Introduction to High-Contrast Imaging

This will serve as a brief primer to high-contrast imaging (HCI) to illustrate the concepts and themes prevalent within this package. This is not meant to be an exhaustive lecture note on the topic, but more of a gentle introduction.

If you are comfortable with HCI topics, consider skipping to the [Getting Started](@ref gettingstarted) section to get introduced to ADI.jl, browse the Examples to see sample workflows, or browse through the API to start exploring the capabilities of ADI.jl.

## What is HCI

HCI is an advanced imaging technique comprising modern instrumentation, clever observational techniques, and post-processing algorithms. The goal of HCI is to probe the circumstellar regions of a star in the search for companions, debris disks, and more. The use of large aperture telescopes, adaptive optics (AO), and coronagraphs are the basis of HCI. An example of an image taken with these instruments is shown below.

![](assets/speckles.png)

What is notable in this image is that there is a lot of structured noise in the image that is overwhelming any potential companion signal. The center of the image is particularly noisy, which is precisely where we are most interested in searching for exoplanets. This noise is the effect of quasi-static speckles in the focal-plane. These speckles occur from non-common path aberrations in the AO system and are a fundamental part of the data received by these instruments. Improving the quality of instrumentation is an active topic of HCI research, but it is beyond the scope of this introduction.

## [Angular Differential Imaging (ADI)](@id adi)

Because this noise is fundamental to the data, post-processing must be performed in order to see any circumstellar signal. This post-processing requires us to fit the PSF of the speckles and then remove it. An example of the above frame with the speckles removed is shown below.

![](assets/S-1.png)

Unfortunately, there is no companion evident; but the speckles have been removed, so what is left? The exoplanet is still sitting below the statistical noise in this frame, but the noise can be averaged out by combining many frames together. Since we are concerned with subtracting the speckles, we need to be careful and consider **how do we fit and subtract the speckles without removing potential companion signal**?

This is where angular differential imaging (ADI) comes in. ADI is an observational technique pioneered in the early 2000s as an extension of *roll deconvolution* for ground-based telescopes. The core of this process is that the quasi-static speckles are a function of the optical system, not the astrophysical source. Throughout a night of observing we can leverage the rotation of the Earth to make the field-of-view (FOV) appear to rotate (on an Alt-Az mounted telescope with the field rotator disabled). Even though the sky appears to rotate, because the speckles are due to the telescope optics they will not appear to rotate. The animation below shows a cube of data with a bright fake companion that illustrates the sky rotation typical of ADI.

![](assets/fake_cube.gif)

By taking this sequence of images (commonly referred to as a cube) we can more easily model and fit the speckle signal separate from any companions. If you median combine the cube as-is, the non-stationary companion signal will attenuate leaving just the speckles. If we derotate the sequence according to the parallactic angles for each frame we align the sky to a common heading. Now we can collapse the derotated sequence and the planet will constructively interfere while the now-rotating speckles will attenuate. The figure below shows these two competing reductions.

![](assets/adi_example.png)


## Post-Processing Algorithms

Using data cubes (as described in the [ADI section](@ref adi)), we are tasked with fitting the speckles without capturing the rotating companion signal. Quite a few algorithms have been proposed and a thorough discussion of them is beyond the scope of this introduction. For now, let's assume the algorithms are a black-box that produce speckle approximation cubes.

If we have this cube, all we need to post-process the data is

1. Retrieve a speckle estimate cube
2. Subtract the speckle estimate from the target cube and form a residual cube
3. Derotate the residual cube according to the parallactic angles of the target
4. Collapse the derotated residual cube

Steps 2-4 are shown in the following figure

![](assets/adi_process.png)

After all this processing, finally the substellar companion HR8799e is evident. Hopefully this shows the difficulty of HCI and builds up part of the process that occurs outside of the reduction you'll be doing with ADI.jl.

## References

Here is a selection of further reading for information about high-contrast imaging, ADI, and similar techniques

* [Traub, Oppenheimer 2010, "Direct Imaging of Exoplanets"](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwj4rKn8_a3tAhVpvFkKHcepDoEQFjAKegQIBRAC&url=https%3A%2F%2Fwww.amnh.org%2Fcontent%2Fdownload%2F53052%2F796511%2Ffile%2FDirectImagingChapter.pdf&usg=AOvVaw0JT9cGTkuFGknAsfvyMxkY)
* [Bowler 2016, "Imaging Extrasolar Giant Planets"](https://ui.adsabs.harvard.edu/abs/2016PASP..128j2001B/abstract)
* [Pueyo 2018, "Direct Imaging as a Detection Technique for Exoplanets"](https://link.springer.com/referenceworkentry/10.1007%2F978-3-319-55333-7_10)
* [Marois et al. 2006, "Angular Differential Imaging: A Powerful High-Contrast Imaging Technique"](https://ui.adsabs.harvard.edu/abs/2006ApJ...641..556M/abstract)
# [Getting Started](@id gettingstarted)

Here is a quick-start guide for people familiar with ADI and experience using tools like [VIP](https://github.com/vortex-exoplanet/VIP) or [pyKLIP](https://pyklip.readthedocs.io/en/latest/). For installation and setup information, see the [Installation and Setup](@ref) section.

## Expected Data Formats

### ADI Cube

For standard ADI data, we store the values in a 3-dimensional array, where the first dimension is temporal, and the remaining dimensions are pixel coordinates. This is how most ADI data are stored on disk (typically in FITS files) and allow specifying operations like a tensor. This cube should already be registered with the star in the center of the frames (note the center is only well-defined for odd-sized frames, even though even-sized frames will work fine).

### Parallactic Angles

The parallactic angles should be stored as *degrees* in a vector. The parallactic angle `X[i]` will result in rotating frame `i` `X[i]` degrees counter-clockwise.

### SDI Cube/Tensor

For standard SDI data, we store the values in a 4-dimensional array, where the first dimension is spectral, the second is temporal, and the remaining dimensions are pixel coordinates. This is how *some* SDI data are stored on disk (typically in FITS files) and allow specifying operations like a tensor. For SDI data that is stored with the temporal axis first, the dimensions should be permuted before processing (see `permutedims`). This cube should also be registered with the star in the center of the frame.

In addition to the SDI tensor and parallactic angles, the list of wavelengths are required (for scaling speckles) and a spectral template can be used. To create a scale list from the wavelengths or a template, use [`scale_list`](@ref).

## Algorithms

The following algorithms are implemented:
* [Classic Subtraction](@ref classic)
* [PCA](@ref pca)
* [NMF](@ref nmf)
* [GreeDS](@ref greeds)

## Processing Patterns

### Full Frame ADI Reduction

Given an algorithm `alg`, we can fully process ADI data by calling `alg` like a function, or using the [`process`](@ref) method

```julia
julia> using ADI

julia> alg = PCA(ncomps=5)

julia> resid = alg(cube, angles)

julia> resid === process(alg, cube, angles)
true
```

### Full Frame RDI Reduction

The only difference here is the inclusion of a reference cube.

```julia
julia> alg = PCA(ncomps=5)

julia> resid = alg(cube, angles; ref=cube_ref)
```

### Reduction Process

The process for producing the flat, residual frame follows this general workflow

1. Create a cube of the speckle approximation, `S`
2. Subtract `S` from the data cube to create the residual cube `R`
3. Derotate `R` frame-by-frame according to the parallactic angle
4. Collapse the derotated `R`

In ADI.jl this process looks like this:

```julia
cube, angles = # load data
S = reconstruct(PCA(10), cube)
R = cube .- S
R_derotate = derotate(R, angles)
resid = collapse(R_derotate)

# or, more succinctly
R = subtract(PCA(10), cube)
resid = collapse(R, angles)
```

Notice how the only part of this specific to the algorithm is [`reconstruct`](@ref)? This lets us have the super-compact functional form from above without having to copy the common code for each algorithm.

### Altering the Geometry

[HCIToolbox.jl](https://github.com/JuliaHCI/HCIToolbox.jl) has utilities for geometrically filtering the input data, such as only taking an annulus of the input cube or iterating over many annuli. This is exactly the purpose of [`AnnulusView`](@ref) and [`MultiAnnulusView`](@ref), which use indexing tricks to retrieve the pixels *only* within the spatial region of interest without having to copy the input data.

If you wrap a cube in one of these views, ADI.jl will handle it automatically (if the algorithm supports it). Since these views filter the pixels, the runtime performance will generally be faster than the full-frame equivalents.

```julia
ann = AnnulusView(cube; inner=15, outer=25)
res = PCA(10)(ann, angles)
```

```julia
# annuli of width 5 starting at 5 pixels and ending at the edge of the cube
anns = MultiAnnulusView(cube, 5; inner=5)
res = PCA(10)(anns, angles)

# use different algorithms for each annulus
N_ann = length(anns.indices)
algs = [PCA(10), PCA(9), PCA(8), ...]
res = process(algs, anns, angles)
```

## Comparison to VIP

ADI.jl took a lot of ideas from VIP and expanded them using the power of Julia. To begin with, Julia typically has smaller and more self-contained packages, so most of the basic image-processing that is used here is actually written in the [HCIToolbox.jl](https://github.com/JuliaHCI/HCIToolbox.jl) package. In the future, I have plans to incorporate forward-modeling distributions in [Firefly.jl](https://github.com/JuliaHCI/Firefly.jl), which currently is an active research project.

Some technical distinctions to VIP

* Julia is 1-indexed. This means all the positions for apertures, bounds, images, etc. start at 1. This is distinct from 0-based indexing in python, but is equivalent to the indexing in DS9 and IRAF.
* Julia's `std` uses the sample statistic (`n-1` degrees of freedom) while numpy's `std` uses the population statistic (`n` degrees of freedom). This may cause very slight differences in measurements that rely on this.
* Aperture mapping - many of the [`Metrics`](@ref) are derived by measuring statistics in an annulus of apertures. In VIP, this ring is not equally distributed- the angle between apertures is based on the exact number of apertures rather than the integral number of apertures that are actually measured. In ADI.jl the angle between apertures is evenly distributed. The same number of pixels are discarded in both packages, but in VIP they all end up in the same region of the image (see [this figure](assets/aperture_masks.png)).
* Collapsing - by default VIP collapses a cube by derotating it then finding the median frame. In ADI.jl, the default collapse method is a weighted sum using the inverse of the temporal variance for weighting. This is documented in `HCIToolbox.collapse` and can be overridden by passing the keyword argument `method=median` or whichever statistical function you want to use.
* Image interpolation - by default VIP uses a `lanczos4` interpolator from opencv, by default ADI.jl uses a bilinear b-spline interpolator through Interpolations.jl
* Annular and framewise processing - some of the VIP algorithms allow you to go annulus-by-annulus and optionally filter the frames using parallactic angle thresholds. ADI.jl does not bake these options in using keyword arguments; instead, the geometric filtering is achieved through [`AnnulusView`](@ref) and [`MultiAnnulusView`](@ref). Parallactic angle thresholds are implemented in the [`Framewise`](@ref) algorithm wrapper. I've separated these techniques because they are fundamentally independent and because it greatly increases the composability of the algorithms.

The biggest difference, though, is Julia's multiple-dispatch system and how that allows ADI.jl to *do more with less code*. For example, the [`GreeDS`](@ref) algorithm was designed explicitly for [`PCA`](@ref), but the formalism of it is more generic than that. Rather than hard-coding in PCA, the GreeDS algorithm was written generically, and Julia's multiple-dispatch  allows the use of, say, [`NMF`](@ref) instead of PCA. By making the code *generic* and *modular*, ADI.jl enables rapid experimentation with different post-processing algorithms and techniques as well as minimizing the code required to implement a new algorithm and be able to fully use the ADI.jl API.

## Feature comparison

Some notable libraries for HCI tasks include [VIP](https://github.com/vortex-exoplanet/VIP), [pyKLIP](https://pyklip.readthedocs.io/en/latest/), and [PynPoint](https://github.com/PynPoint/PynPoint). A table of the feature sets of these packages alongside ADI.jl is presented below. In general VIP offers the most diversity in algorithms and their applications, but not all algorithms are as feature-complete as the PCA implementation. VIP also contains many useful utilities for pre-processing and a pipeline framework. pyKLIP primarily uses the PCA (KLIP) algorithm, but offers many forward modeling implementations. PynPoint has a highly modular pre-processing module that is focused on pipelines.

| - | Pre. | Algs. | Techs. | D.M. | F.M. |
|:---:|:---:|:---:|:---:|:---:|:---:|
| ADI.jl | ✗ | median, LOCI, PCA, NMF, fixed-point GreeDS | Full-frame ADI/RDI, SDI (experimental), annular ADI* | detection maps, STIM, SLIMask, contrast curve | ✗ |
| VIP | ✓ | median, LOCI, PCA, NMF, LLSG, ANDROMEDA, pairwise frame differencing | Full-frame ADI/RDI, SDI, annular ADI/RDI* | detection maps, blob detection, STIM, ROC, contrast curve | NegFC |
| pyKLIP | ✗ | PCA, NMF, weighted PCA | Full-frame ADI/RDI, SDI, annular ADI/RDI | detection maps, blob detection, contrast curve, cross-correlation | KLIP-FM, Planet Evidence, matched filter (FMMF), spectrum fitting, DiskFM |
| PynPoint | ✓ | median, PCA | Full-frame ADI/RDI, SDI | detection maps, contrast curve | ✗ |

**Column labels:** Pre-processing, Algorithms, Techniques, Detection Metrics, Forward Modeling.

*Techniques marked with * indicate partial support, meaning that not all algorithms are supported.*
# SDI

!!! warning
    SDI should be considered _very_ experimental.

```@docs
ADI.SDIAlgorithm
```

## API/Reference

```@docs
SingleSDI
DoubleSDI
SliceSDI
scale_list
```
# Metrics

```@index
Modules = [ADI.Metrics]
```

```@docs
ADI.Metrics
```

## Detection maps

```@docs
detectionmap
snr
significance
noise
stimmap
stim_threshold
Metrics.stim
```

## Ensemble methods

The following methods utilize the results of multiple ADI reductions, in some form.

```@docs
slimmap
```

## Throughput

```@docs
throughput
```

## Contrast curve

```@docs
contrast_curve
Metrics.subsample_contrast
Metrics.estimate_starphot
```
# Framewise

For an example of [`Framewise`](@ref) in use, along with spatial filtering, see [this example](@ref ex2).

```@docs
Framewise
```

## API/Reference

```@docs
AnnulusView
MultiAnnulusView
eachannulus
inverse
inverse!
```
# Benchmarks

The large scale image-processing required for ADI algorithms can lead to concerns about runtime efficiency. To this end, ADI.jl (and the associated JuliaHCI packages) are developed with performance in mind. These packages do not aim to be as fast as possible; rather they focus on being as fast as *is convenient* (for the users and the devs).

The [Vortex Imaging Pipeline](https://github.com/vortex-exoplanet/vip) (VIP) is the inspiration for ADI.jl. It is one of the major Python HCI packages and it offers many more features than ADI.jl. Some of the common uses for both packages include full-frame ADI processing, S/N maps, and contrast curves.

### System/Setup Information

The benchmarks here can be found in the [`bench/`](https://github.com/JuliaHCI/ADI.jl/blob/main/bench/) folder organized into Julia files. The benchmarks utilize BenchmarkTools.jl, PyCall.jl with `virtualenv`, and CSV.jl for accuracy, reproducibility, and organization.

```
Julia Version 1.6.0-beta1
Commit b84990e1ac* (2021-01-08 12:42 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin19.6.0)
  CPU: Intel(R) Core(TM) i5-8259U CPU @ 2.30GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.0 (ORCJIT, skylake)
Environment:
  JULIA_NUM_THREADS = 4
```

For the python code, there is a `requirements.txt` file in `bench/`. To reproduce this environment, (optionally) activate a virtual environment then install from the requirements file.

```
(venv) $ pip install -r requirements.txt
```

For reproducibility, there is a `Manifest.toml` file in `bench/`. To reproduce this environment, first activate it, then instantiate it

```
$ julia --project=bench -e 'using Pkg; Pkg.instantiate()'
```

!!! warning "PyCall.jl and virtual environments"
    The interface between Julia and python is handled by [PyCall.jl](https://github.com/juliapy/PyCall.jl). When using a virtual environment, PyCall may not use the correct python library. Before running the benchmarks, please read [this reference](https://github.com/juliapy/PyCall.jl#python-virtual-environments).

!!! tip "Multi-threading"
    Some of the image-processing methods in ADI.jl and HCIToolbox.jl are multi-threaded, and will lead to a noticable difference in some benchmarks. To take advantage of this, set the environment variable `JULIA_NUM_THREADS` before starting your runtime. [Multi-Threading documentation](https://docs.julialang.org/en/v1/manual/multi-threading/).

```@setup bench
using CSV
using DataFrames
using StatsPlots
benchdir(args...) = joinpath("..", ".." ,"bench", args...);
```

## ADI Reduction

These benchmarks show the duration to fully reduce ADI data for various algorithms. The data used are $\beta$ Pictoris and HR8799 from [HCIDatasets.jl](https://github.com/JuliaHCI/HCIDatasets.jl).

```@example bench
adi_data = CSV.File(benchdir("adi_benchmarks.csv")) |> DataFrame |> sort!
cube_labels = @. ifelse(adi_data.N == 622261, "Beta Pictoris", "HR8799")
insertcols!(adi_data, 4, :cube => cube_labels)
adi_groups = groupby(adi_data, :framework)
```

```@example bench
cube_groups = groupby(adi_data, :cube)
plot(
    @df(cube_groups[1], groupedbar(:alg, :time, group=:framework, yscale=:log10)),
    @df(cube_groups[2], groupedbar(:alg, :time, group=:framework)),
    size=(700, 350),
    leg=:topleft,
    ylabel="time (s)",
    title=["Beta Pictoris" "HR8799"]
)
```

*Please note the log-scale for the left figure.*

## Detection Maps

This benchmark measures the duration to produce a signal-to-noise ratio (S/N) map. Rather than test exact cubes, these benchmarks test randomly generated frames of various sizes. The FWHM is fixed at 5.

```@example bench
snrmap_data = CSV.File(benchdir("snrmap_benchmarks.csv")) |> DataFrame |> sort!
snrmap_groups = groupby(snrmap_data, :framework)
```

```@example bench
@df snrmap_data scatter(
    :N,
    :time,
    group=:framework,
    ms=6,
    xlabel="number of pixels",
    ylabel="time (s)"
)
```

## Contrast Curves

Finally, this benchmark measures the duration to generate a contrast curve for analyzing the algorithmic throughput of an ADI algorithm. For both benchmarks 3 azimuthal branches are used for throughput injections and a FWHM of 8. A Gaussian PSF function is evaluated in a `(21, 21)` grid for the injections. The data used are $\beta$ Pictoris and HR8799 from [HCIDatasets.jl](https://github.com/JuliaHCI/HCIDatasets.jl).

```@example bench
contrast_data = CSV.File(benchdir("contrast_benchmarks.csv")) |> DataFrame |> sort!
cube_labels = @. ifelse(contrast_data.N == 622261, "Beta Pictoris", "HR8799")
insertcols!(contrast_data, 4, :cube => cube_labels)
contrast_groups = groupby(contrast_data, :framework)
```

```@example bench
@df contrast_data groupedbar(
    :cube,
    :time,
    group=:framework,
    leg=:topleft,
    ylabel="time (s)",
    yscale=:log10,
)
```

*Please note the log-scale.*
# [GreeDS](@id greeds)

## API/Reference

```@docs
GreeDS
```
# API/Reference

```@index
Pages = ["algorithms/api.md",
         "algorithms/classic.md",
         "algorithms/greeds.md",
         "algorithms/loci.md",
         "algorithms/nmf.md",
         "algorithms/pca.md"]
```

```@docs
ADI.ADIAlgorithm
reconstruct
subtract
process
ADI.fit
ADI.design
ADI.ADIDesign
ADI.LinearDesign
```
# [PCA](@id pca)


## API/Reference

```@docs
PCA
TPCA
```# [NMF](@id nmf)

## API/Reference

```@docs
NMF
```# [LOCI](@id loci)

## API/Reference

```@docs
LOCI
```
# [Classic](@id classic)

## API/Reference

```@docs
Classic
ADI.ClassicDesign
```
