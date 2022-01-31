---
title: 'deeplenstronomy: A dataset simulation package for strong gravitational lensing'
tags:
  - Python
  - astronomy
  - strong lensing
  - simulation
authors:
  - name: Robert Morgan^[Corresponding author]
    orcid: 0000-0002-7016-5471
    affiliation: "1, 2" 
  - name: Brian Nord
    orcid: 0000-0001-6706-8972
    affiliation: "3, 4"
  - name: Simon Birrer
    orcid: 0000-0003-3195-5507
    affiliation: 5
  - name: Joshua Yao-Yu Lin
    orcid: 0000-0003-0680-4838
    affiliation: 6
  - name: Jason Poh
    orcid: 0000-0002-5040-093X
    affiliation: 4
affiliations:
  - name: University of Wisconsin-Madison
    index: 1
  - name: Legacy Survey of Space and Time Data Science Fellowship Program
    index: 2
  - name: Fermi National Accelerator Laboratory
    index: 3
  - name: University of Chicago
    index: 4
  - name: Stanford University
    index: 5
  - name: University of Illinois Urbana-Champaign
    index: 6
date: 28 October 2020
bibliography: paper.bib
---

# Background

Astronomical observations and statistical modeling permit the high-fidelity analysis of strong gravitational lensing (SL) systems, which display an astronomical phenomenon in which light from a distant object is deflected by the gravitational field of another object along its path to the observer.
These systems are of great scientific interest because they provide information about multiple astrophysical and cosmological phenomena, including the nature of dark matter, the expansion rate of the Universe, and characteristics of galaxy populations. 
They also serve as standing tests of the theory of General Relativity and modified theories of gravity. 

Traditional searches for SL systems have involved time- and effort-intensive visual or manual inspection of images by humans to identify characteristic features --- like arcs, particular color combinations, and  object orientations. 
However, a comprehensive search using the traditional approach is prohibitively expensive for large numbers of images, like those in cosmological surveys --- e.g., the Sloan Digital Sky Survey [@sdss], the Dark Energy Survey [@des], and the Legacy Survey of Space and Time (LSST) [@lsst]. 
To automate the SL detection process, techniques based on machine learning (ML) are beginning to overtake traditional approaches for scanning  astronomical images. 
In particular, deep learning techniques have been the focus, but they require large sets of labeled images to train these models. 
Because of the relatively low number of observed SL systems, simulated datasets of images are often needed. 
Thus, the composition and production of these simulated datasets have become integral parts of the SL detection process.

One of the premier tools for simulating and analyzing SL systems, `lenstronomy` [@lenstronomy], works by the user specifying the properties of the physical systems, as well as how they are observed (e.g., telescope and camera) through a `python`-based application programming interface (API) to generate a single image. 
Generating populations of SL systems that are fit for neural network training requires additional infrastructure. 

# Statement of need 

Due to the inherent dependence of the performance of ML approaches on their training data, the deep learning approach to SL detection is in tension with scientific reproducibility without a clear prescription for the simulation of the training data. 
There is a critical need for a tool that simulates full datasets in an efficient and reproducible manner, while enabling the use of all the features of the `lenstronomy` simulation API. 
Additionally, this tool should  simplify user interaction with `lenstronomy` and organize the simulations and associated metadata into convenient data structures for deep learning problems.
Multiple packages have been developed to generate realistic training data by wrapping around `lenstronomy`: `baobab` [@baobab] generates training sets for lens modeling and hierarchical inference and the LSST Dark Energy Science Collaboration's `SLSprinkler` [@lsstdescsl] adds strongly lensed variable objects into catalogs and images. 
Nonetheless, the need for a simple, general tool capable of efficiently simulating any astronomical system in a reproducible manner while giving the user complete freedom to set the properties of objects remains. 


# Summary

`deeplenstronomy` generates SL datasets by organizing and expediting user interaction with `lenstronomy`. 
The user creates a single yaml-style configuration file that describes the aspects of the dataset: number of images, properties of the telescope and camera, cosmological parameters, observing conditions, properties of the physical objects, and geometry of the SL systems. 
`deeplenstronomy` parses the configuration file and generates the dataset, producing both the images and the parameters that led to the production of each image as outputs. 
The configuration files can easily be shared, enabling users to easily reproduce each other's training datasets.

The premier objective of `deeplenstronomy` is to help astronomers make their training datasets as realistic as possible. 
To that end, `deeplenstronomy` contains built-in features for the following functionalities: use any stellar light profile or mass profile in `lenstronomy`; simulate a variety of astronomical systems such as single galaxies, foreground stars, galaxy clusters, supernovae, and kilonovae, as well as any combination of those systems; fully control the placements of objects in the simulations; use observing conditions of real astronomical surveys; draw any parameter from any probability distribution; introduce any correlation; and incorporate real images into the simulation.
Furthermore, `deeplenstronomy` facilitates realistic time-domain studies by providing access to public spectral energy distributions of observed supernovae and kilonovae and incorporating the transient objects into time series of simulated images.
Finally, `deeplenstronomy` provides data visualization functions to enable users to inspect their simulation outputs.
These features and the path from configuration file to full data set are shown in \autoref{fig:flowchart}.

`deeplenstronomy` makes use of multiple open-source software packages: `lenstronomy` is used for all gravitational lensing calculations and image simulation; `numpy` [@numpy] `Array`s are used internally to store image data and perform vectorized calculations; `pandas` [@pandas] `DataFrame`s are utilized for storing simulation metadata and file reading and writing; `scipy` [@scipy] is used for integration and interpolation; `matplotlib` [@matplotlib] functions are used for image visualization; `astropy` [@astropy] is used for cosmological calculations and color image production; `h5py` [@h5py] is utilized for saving images; and `PyYAML` [@pyyaml] is used to manage the configuration file. 
While not used directly, some `python-benedict` [@benedict] functionalities helped to create `deeplenstronomy`'s data structures and internal search algorithms. 

`deeplenstronomy` is packaged and disseminated via [PyPI](https://pypi.org/project/deeplenstronomy/). 
Documentation and example notebooks are available on the [`deeplenstronomy` website](https://deepskies.github.io/deeplenstronomy/). 
Any bugs or feature requests can be opened as issues in the [GitHub
repository](https://github.com/deepskies/deeplenstronomy/issues) [@deeplenstronomy].

![The `deeplenstronomy` process. Dataset properties, camera and telescope properties, observing conditions, object properties (e.g., `lenstronomy` light and mass profiles, point sources, and temporal behavior), the geometry of the SL systems, and optional supplemental input files (e.g., probability distributions, covariance matrices, and image backgrounds) are specified in the main configuration file. `deeplenstronomy` then intreprets the configuration file, calls `lenstronomy` simulation functionalities, and organizes the resulting images and metadata.\label{fig:flowchart}](flowchart.png)

# Acknowledgements

R. Morgan thanks the LSSTC Data Science Fellowship Program, which is funded by LSSTC, NSF Cybertraining Grant #1829740, the Brinson Foundation, and the Moore Foundation; his participation in the program has benefited this work. 
R. Morgan also thanks the Universities Research Association Fermilab Visiting Scholar Program for funding his work on this project.

We acknowledge the Deep Skies Lab as a community of multi-domain experts and collaborators who’ve facilitated an environment of open discussion, idea-generation, and collaboration. 
This community was important for the development of this project.
We acknowledge contributions from Joao Caldeira during the early stages of this project.

Work supported by the Fermi National Accelerator Laboratory, managed and operated by Fermi Research Alliance, LLC under Contract No. DE-AC02-07CH11359 with the U.S. Department of Energy. 
The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a non-exclusive, paid-up, irrevocable, world-wide license to publish or reproduce the published form of this manuscript, or allow others to do so, for U.S. Government purposes.


# References
# Welcome to `deeplenstronomy`!

[![status](https://joss.theoj.org/papers/e978dd566d1f290055a02d76288e95e1/status.svg)](https://joss.theoj.org/papers/e978dd566d1f290055a02d76288e95e1)
[![status](https://img.shields.io/badge/arXiv-2102.02830-red)](http://arxiv.org/abs/2102.02830)
[![status](https://img.shields.io/badge/PyPi-0.0.2.3-blue)](https://pypi.org/project/deeplenstronomy/)
[![status](https://img.shields.io/badge/License-MIT-lightgrey)](https://github.com/deepskies/deeplenstronomy/blob/master/LICENSE)

`deeplenstronomy` is a tool for simulating large datasets for applying deep learning to strong gravitational lensing. 
It works by wrapping the functionalities of [`lenstronomy`](https://github.com/sibirrer/lenstronomy) in a convenient yaml-style interface, allowing users to embrace the astronomer part of their brain rather than their programmer part when generating training datasets.

## Installation

**With conda (Recommended)**

- Step 0: Set up an environment. This can be done straightforwardly with a `conda` installation:

```
conda create -n deeplens python=3.7 jupyter scipy pandas numpy matplotlib astropy h5py PyYAML mpmath future
conda activate deeplens
```

- Step 1: `pip install lenstronomy`
- Step 2: `pip install deeplenstronomy`

**With pip**

- Step 1: `pip install deeplenstronomy`

## [Getting Started and Example Notebooks](https://deepskies.github.io/deeplenstronomy/Notebooks/)

Start by reading the [Getting Started Guide](https://deepskies.github.io/deeplenstronomy/Notebooks/GettingStarted.html) to familiarize yourself with the `deeplenstronomy` style.

After that, check out the example notebooks below:

### Notebooks for `deeplenstronomy` Utilities
- [Creating `deeplenstronomy` Configuration Files](https://deepskies.github.io/deeplenstronomy/Notebooks/ConfigFiles.html)
- [Generating Datasets](https://deepskies.github.io/deeplenstronomy/Notebooks/DeepLenstronomyDemo.html)
- [Visualizing `deeplenstronomy` Images](https://deepskies.github.io/deeplenstronomy/Notebooks/Visualization.html)
- [Utilizing Astronomical Surveys](https://deepskies.github.io/deeplenstronomy/Notebooks/Surveys.html)
- [Defining Your Own Probability Distributions](https://deepskies.github.io/deeplenstronomy/Notebooks/UserDistributions.html)
- [Using Your Own Images as Backgrounds](https://deepskies.github.io/deeplenstronomy/Notebooks/BackgroundsDemo.html)
- [Simulating Time-Series Datasets](https://deepskies.github.io/deeplenstronomy/Notebooks/TimeSeriesDemo.html)

### Notebooks for Applying `deeplenstronomy` to Machine Learning Analyses
- [Using `deeplenstronomy` for Active Learning](https://deepskies.github.io/deeplenstronomy/Notebooks/ActiveUpdateDemo.html)
- [Using `deeplenstronomy` for Classification and Regression](https://deepskies.github.io/deeplenstronomy/Notebooks/Metrics.html)

### Notebooks for Suggested Science Cases
- [A Walkthrough of Using `deeplenstronomy` for Science](https://deepskies.github.io/deeplenstronomy/Notebooks/FullExample.html)


## API Documentation

`deeplenstronomy` is designed so that users only need to work with their personal configuration files and the dataset generatation and image visualization functions.
However, if you would like to view the full API documentation, you can visit the [docs](https://deepskies.github.io/deeplenstronomy/docs/) page.

## Citation

If you use `deeplenstronomy` in your work, please include the following citations:
```
@article{deeplenstronomy,
  doi = {10.21105/joss.02854},
  url = {https://doi.org/10.21105/joss.02854},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {58},
  pages = {2854},
  author = {Robert Morgan and Brian Nord and Simon Birrer and Joshua Yao-Yu Lin and Jason Poh},
  title = {deeplenstronomy: A dataset simulation package for strong gravitational lensing},
  journal = {Journal of Open Source Software}
}

@article{lenstronomy,
    title     =   "lenstronomy: Multi-purpose gravitational lens modelling software package",
    journal   =   "Physics of the Dark Universe",
    volume    =   "22",
    pages     =   "189 - 201",
    year      =   "2018",
    issn      =   "2212-6864",
    doi       =   "10.1016/j.dark.2018.11.002",
    url       =   "http://www.sciencedirect.com/science/article/pii/S2212686418301869",
    author    =   "Simon Birrer and Adam Amara",
    keywords  =   "Gravitational lensing, Software, Image simulations"
}
```

## Contact

If you have any questions or run into any errors with the beta release of `deeplenstronomy`, please don't hesitate to reach out:

Rob Morgan 
<br>
robert [dot] morgan [at] wisc.edu

You can also message me on the DES, DELVE, LSSTC, deepskies, or lenstronomers Slack workspaces





<!---
.. image:: https://badge.fury.io/py/deeplenstronomy.png
    :target: http://badge.fury.io/py/deeplenstronomy

.. image:: https://travis-ci.org/bnord/deeplenstronomy.png?branch=master
    :target: https://travis-ci.org/bnord/deeplenstronomy
--->



# DeepLenstronomy Configuration Files

An overview of how to structure the main configuration file for `deeplenstronomy`.

## Introduction

The purpose of `deeplenstronomy` is to wrap the strong-lensing image simulation utilities of `lenstronomy` into a convenient and reproducible framework that enables various science analyses.
This is achieved in `deeplenstronomy` by allowing the user to generate an entire large-scale dataset from a single, human-readable configuration file.

This tutorial will help you to create the main configuration file for problem you are working on.

### Contents
- Sections in the Config File
- Drawing from Probability Distributions
- Advanced Features
- Troubleshooting

## Sections in the Config File

The main configuration file for `deeplenstronomy` is a yaml-style document.
The top-level keywords in this document will tell `deeplenstronomy` everything it needs to know to build you a dataset. These keywords are **mandatory** for `deeplenstronomy` to function properly.

_Note_: You must use spaces (not tabs) for indentation in `deeplenstronomy` configuration files.

### The DATASET Section

This section will contain information regarding the dataset as a whole, rather than the individual images or systems you want to generate. 
Here you must specify the `NAME` and `SIZE` of the dataset. 
If you plan on saving the dataset, you must also specify the `OUTDIR` where you want the simulated images and associated metadata to be located.
This is done with a block like this

```
DATASET:
    NAME: MyDeeplenstronomyDataset
    PARAMETERS:
        SIZE: 100
        OUTDIR: MySimulationResults
	SEED: 6     # optional
```
This block will set the name of the dataset to `MyDeeplenstronomyDataset`, specify that the dataset will contain 100 images, set a random seed for the simulation, and specify that the simulated images and associated metadata will be saved in a directory named `MySimulationResults`. 
If this directory does not already exist on your system, `deeplenstronomy` will create it automatically.

### The COSMOLOGY Section

This section will specify the cosmology with which `lenstronomy` will perform all calculations. In the config file, the bare-minimum cosmology section would look like this

```
COSMOLOGY:
    PARAMETERS:
        H0: 70
        Om0: 0.3
```
Internally, all that happens within `deeplenstronomy` is the instantiation of an `astropy.cosmology.FlatLambdaCDM` object. Thus, you can specify the value of `H0`, `Om0`, `Tcmb0`, `Neff`, `m_nu`, and `Ob0`. You can also view the [astropy documentation](https://docs.astropy.org/en/stable/api/astropy.cosmology.FlatLambdaCDM.html#flatlambdacdm) for more information on these parameters.

### The IMAGE Section

This section will specify all the information needed to make simulated images look realistic. The following parameters are required for `deeplenstronomy` to generate images:

- `exposure_time`: the number of seconds used to capture the image
- `numPix`: the width/height of the images in pixels
- `pixel_scale`: the number of arcseconds in each pixel
- `psf_type`: describes the model for the point-spread-function of the instrument
- `read_noise`: standard deviation of noise generated by read-out (in units of electrons)
- `ccd_gain`: the number of photo-electrons required for each analog-to-digital unit

A typical `IMAGE` block will look like this

```
IMAGE:
    PARAMETERS:
        exposure_time: 90
        numPix: 100
        pixel_scale: 0.263
        psf_type: 'GAUSSIAN'
        read_noise: 7
        ccd_gain: 6.083
```

### The SURVEY Section

This section specifies the properties of the survey that collected the dataset you are seeking to simulate. This could be the characteristics of DES, LSST, SDSS, ZTF, etc. The parameters `deeplenstronomy` needs in this section are:

- `BANDS`: a comma-separated list of the optical filters used
- `seeing`: the FWHM in arcseconds of the measured PSF from the observations
- `magnitude_zero_point`: the magnitude at which one photo-electron is measured each second
- `sky_brightness`: measured light contamination from the atmosphere in magnitude per square arcsecond
- `num_exposures`: the number of individual exposures in each coadded image

Thus, an example `SURVEY` block could look like this

```
SURVEY:
    PARAMETERS:
        BANDS: g,r,i,z,Y
        seeing: 0.9
        magnitude_zero_point: 30.0
        sky_brightness: 23.5
        num_exposures: 10
```
Often the parameters in the `SURVEY` section will have a different value for each band. This relationship is technically a correlation between the band and the parameter, and therefore is covered in the **Advanced Features** section of this document.

### Half-Way-There Summary

Just to recap what we've talked about, let's put everything together. So far, we've covered the top-level aspects of the simulation: the properties of the dataset, the cosmological model for calculations, the image properties, and the survey conditions. Our in-progress config file currently looks like this:

```
DATASET:
    NAME: MyDeeplenstronomyDataset
    PARAMETERS:
        SIZE: 100
        OUTDIR: MySimulationResults
	SEED: 6

COSMOLOGY:
    PARAMETERS:
        H0: 70
        Om0: 0.3

IMAGE:
    PARAMETERS:
        exposure_time: 90
        numPix: 100
        pixel_scale: 0.263
        psf_type: 'GAUSSIAN'
        read_noise: 7
        ccd_gain: 6.083

SURVEY:
    PARAMETERS:
        BANDS: g,r,i,z,Y
        seeing: 0.9
        magnitude_zero_point: 30.0
        sky_brightness: 23.5
        num_exposures: 10
```

Next, we'll look at how to specify the properties of the physical objects in the images as well as the geometry of the strong-lensing systems you're interested in.

### The SPECIES Section

This section specifies the properties of each object or noise source that could potentially be in one of your images. 
There are currently three species of object supported in `deeplenstronomy`: `GALAXY`, `POINTSOURCE`, and `NOISE`. The most information-rich species is `GALAXY`, so let's start there.

**GALAXY.**
This species has up to 5 attributes for which you can provide information. 
The first of these is the `NAME` keyword, which `deeplenstronomy` will use internally to group the properties of the object. Thus, the `NAME` attribute is mandatory and must be unique. 
The next three attributes (light profiles, mass profiles, and shear profiles) are the core of the `lenstronomy ` simulation.

Light profiles are specified with the `LIGHT_PROFILE_X` keyword and describe how light is distributed around the object.
The `X` represents an index to differentiate multiple light profiles in the case where you would like to utilize a superposition of profiles and must be indexed incrementally starting at 1.
The available profiles are any that are available in `lenstronomy`, and you can look at the `lenstronomy` documentation to see all the options and associated parameters.

As an example, let's create a `GALAXY` species that we'll name "LENS" with a "SERSIC_ELLIPSE" light profile.

```
SPECIES:
    GALAXY_1:
        NAME: LENS
        LIGHT_PROFILE_1:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 19.5
                center_x: 0.0
                center_y: 0.0
                R_sersic: 10
                n_sersic: 4
                e1: 0.2
                e2: -0.1
```

The parameters you need to specify will depend on the `lenstronomy` light profile you select, so check the `lenstronomy` documentation for what parameters you need.
`deeplenstronomy` will warn you if you supply an incorrect parameter. `deeplenstronomy` light profiles use `magnitude` as the defualt way of specifying the brightness; in `lenstronomy` you also have the `amp` parameter available, but that it not utilized in `deeplenstronomy`.
If we wanted to add a second light profile to this `GALAXY`, perhaps to add more light in the center, we can do so by adding a second section like this:

```
SPECIES:
    GALAXY_1:
        NAME: LENS
        LIGHT_PROFILE_1:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 19.5
                center_x: 0.0
                center_y: 0.0
                R_sersic: 10
                n_sersic: 4
                e1: 0.2
                e2: -0.1
        LIGHT_PROFILE_2:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 18.0
                center_x: 0.0
                center_y: 0.0
                R_sersic: 3
                n_sersic: 8
                e1: 0.05
                e2: -0.05
```

The next bit of information `lenstronomy` needs to calculate how light is lensed is the mass distibution.
Since we have conveniently named this `GALAXY` "LENS", it's not too farfetched to image that it will act as the lensing galaxy when this is all said and done, so we'll give it what `deeplenstronomy` calls a mass profile.

Mass profiles will also draw directly from the `lenstronomy` allowed models, so check the `lenstronomy` documentation for all your options. In `deeplenstronomy`, you will adopt the same indexing scheme as the light profiles. Let's add a mass profile to "LENS".

```
SPECIES:
    GALAXY_1:
        NAME: LENS
        LIGHT_PROFILE_1:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 19.5
                center_x: 0.0
                center_y: 0.0
                R_sersic: 10
                n_sersic: 4
                e1: 0.2
                e2: -0.1
        LIGHT_PROFILE_2:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 18.0
                center_x: 0.0
                center_y: 0.0
                R_sersic: 3
                n_sersic: 8
                e1: 0.05
                e2: -0.05
        MASS_PROFILE_1:
            NAME: SIE 
            PARAMETERS:
                theta_E: 2.0
                e1: 0.1
                e2: -0.1
                center_x: 0.0
                center_y: 0.0
```

Here we've used `lenstronomy`'s "SIE" model and specified the necessary parameters for it. 
Let's finish up describing "LENS" by giving it a shear profile.
A shear profile describes how _all_ the light in the system is lensed by an external source, perhaps a dark matter halo.
You can add this type of profile to "LENS" as follows:

```
SPECIES:
    GALAXY_1:
        NAME: LENS
        LIGHT_PROFILE_1:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 19.5
                center_x: 0.0
                center_y: 0.0
                R_sersic: 10
                n_sersic: 4
                e1: 0.2
                e2: -0.1
        LIGHT_PROFILE_2:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 18.0
                center_x: 0.0
                center_y: 0.0
                R_sersic: 3
                n_sersic: 8
                e1: 0.05
                e2: -0.05
        MASS_PROFILE_1:
            NAME: SIE 
            PARAMETERS:
                theta_E: 2.0
                e1: 0.1
                e2: -0.1
                center_x: 0.0
                center_y: 0.0
        SHEAR_PROFILE_1:
            NAME: SHEAR
            PARAMETERS:
                gamma1: 0.08
                gamma2: 0.01
```

Here we've used `lenstronomy`'s "SHEAR" model and specified the necessary parameters for it to work. 

To specify `GALAXY` objects with different properties, just create a new entry in the SPECIES section. Since we will want to do some lensing eventually, lets create another galaxy that will act the light that is actually being lensed.
Fittingly, we might choose to name this `GALAXY` object "SOURCE".
The SPECIES section would then look like this

```
SPECIES:
    GALAXY_1:
        NAME: LENS
        LIGHT_PROFILE_1:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 19.5
                center_x: 0.0
                center_y: 0.0
                R_sersic: 10
                n_sersic: 4
                e1: 0.2
                e2: -0.1
        LIGHT_PROFILE_2:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 18.0
                center_x: 0.0
                center_y: 0.0
                R_sersic: 3
                n_sersic: 8
                e1: 0.05
                e2: -0.05
        MASS_PROFILE_1:
            NAME: SIE 
            PARAMETERS:
                theta_E: 2.0
                e1: 0.1
                e2: -0.1
                center_x: 0.0
                center_y: 0.0
        SHEAR_PROFILE_1:
            NAME: SHEAR
            PARAMETERS:
                gamma1: 0.08
                gamma2: 0.01
    GALAXY_2:
        NAME: SOURCE
        LIGHT_PROFILE_1:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 21.5
                center_x: 0.0
                center_y: 0.0
                R_sersic: 6
                n_sersic: 5
                e1: 0.2
                e2: -0.1 
        SHEAR_PROFILE_1:
            NAME: SHEAR
            PARAMETERS:
                gamma1: 0.08
                gamma2: 0.01               
```

Since we'll be using this galaxy as the so-called source in our lensing system, we don't need to worry about specifying a mass profile.
We have also given it the same shear profile as the lens, since they will be part of the same system.

While the `GALAXY` objects will perhaps require the most thought when constructing the config file, the `POINTSOURCE` and `NOISE` objects will be much more straightforward.

**POINTSOURCE.**
This object will allow you to add a point-like light source to the images.
You may wish to utilize this feature to include foreground stars, to make a galaxy look like an AGN, or to include a supernova or other transient in the image.
I'll show all three of those possibilities here.

Let's start with an AGN.
In similar fashion to the way we specified the different `GALAXY` objects, we'll describe a point source object as follows:

```
<Still inside SPECIES section>
    POINTSOURCE_1:
        NAME: AGN
        HOST: SOURCE
        PARAMETERS:
            magnitude: 16
```
Each `POINTSOURCE` object is characterized by a `NAME`, a `HOST`, and `magnitude`, `sep`, `sep_unit`, and `angle` parameters. The `NAME` attribute will be used to track all properties of this object, so make sure all names are unique. 
The `HOST` attribute will tell `deeplenstronomy` where your object will be in the image. 

The additional parameters determine how the `POINTSOURCE` will appear.
The `sep` parameter will specify how far to distance the `POINTSOURCE` from its host.
If this parameter is not specified, the `POINTSOURCE` will be placed exactly on the center of the `HOST`.
If you specify `sep`, you also need to sepcify the `sep_unit` (either `kpc` or `arcsec`) to indicate the physical units you wish to use.
The `angle` parameter can also be specified if you wish to place the `POINTSOURCE` in a certain orientation with respect to its host. The units of `angle` are radians and the value must be between 0 and 2$\pi$.
In this case, we have selected to place this point source on the center of the `GALAXY` object that we named "SOURCE".
Lastly, the magnitude is the apparent magnitude of the point source.

If we wanted, we could also add a supernova-like point source by specifying a nonzero `sep` parameter. The config entry would then look like this:

```
<Still inside SPECIES section>
    POINTSOURCE_2:
        NAME: SUPERNOVA
        HOST: SOURCE
        PARAMETERS:
            magnitude: 21.0
            sep: 2.0
	    sep_unit: arcsec
```
Specifying this second `POINTSOURCE` for the `GALAXY` we named "SOURCE" **will not** by itself put both `POINTSOURCES` in "SOURCE" at the same time.
That aspect is not determined until the `GEOMETRY` section of the config file.
Here, you are just preparing `deeplenstronomy` to put your `POINTSOURCE`s in "SOURCE" later on.

The final type way you may choose to use the `POINTSOURCE` object is to include foreground stars. This is done by setting the `HOST` keyword to 'Foreground' as follows:

```
<Still inside SPECIES section>
    POINTSOURCE_3:
        NAME: STAR
        HOST: Foreground
        PARAMETERS:
            magnitude: 14.0
```
This setting will place the point source randomly in your image, and it will not contribute to or be affected by any of the lensing calculations.

**NOISE.**
The final type of object that `deeplenstronomy` allows you to include in your images is the `NOISE` object.
Currently, only `POISSON_NOISE` is built in to `deeplenstronomy`, but other types of noise are on the way.
If you would like a specific type of noise, just ping me and I'll make it happen.

You can define a `NOISE` species like this:

```
<Still inside SPECIES section>
    NOISE_1:
        NAME: POISSON_NOISE
        PARAMETERS:
            mean: 2.0
```
In this implementation, each pixel in the image will have a noise value added to it, and that noise value will be drawn from a Poisson distribution with mean 2.0.
The values assigned to different pixels are assumed to be independent.

### Summary of the SPECIES Section

That certainly was a lot of material.
The `SPECIES` section carries all the properties of all the objects in your entire dataset, so it is important to create it with care.
To summarize this component of the config file, let's look at the `SEPCIES` section in our example as a whole.

```
SPECIES:
    GALAXY_1:
        NAME: LENS
        LIGHT_PROFILE_1:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 19.5
                center_x: 0.0
                center_y: 0.0
                R_sersic: 10
                n_sersic: 4
                e1: 0.2
                e2: -0.1
        LIGHT_PROFILE_2:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 18.0
                center_x: 0.0
                center_y: 0.0
                R_sersic: 3
                n_sersic: 8
                e1: 0.05
                e2: -0.05
        MASS_PROFILE_1:
            NAME: SIE 
            PARAMETERS:
                theta_E: 2.0
                e1: 0.1
                e2: -0.1
                center_x: 0.0
                center_y: 0.0
        SHEAR_PROFILE_1:
            NAME: SHEAR
            PARAMETERS:
                gamma1: 0.08
                gamma2: 0.01
    GALAXY_2:
        NAME: SOURCE
        LIGHT_PROFILE_1:
            NAME: SERSIC_ELLIPSE
            PARAMETERS:
                magnitude: 21.5
                center_x: 0.0
                center_y: 0.0
                R_sersic: 6
                n_sersic: 5
                e1: 0.2
                e2: -0.1 
        SHEAR_PROFILE_1:
            NAME: SHEAR
            PARAMETERS:
                gamma1: 0.08
                gamma2: 0.01               
    POINTSOURCE_1:
        NAME: AGN
        HOST: SOURCE
        PARAMETERS:
            magnitude: 16
    POINTSOURCE_2:
        NAME: SUPERNOVA
        HOST: SOURCE
        PARAMETERS:
            magnitude: 21.0
            sep: 2.0
            sep_unit: arcsec
    POINTSOURCE_3:
        NAME: STAR
        HOST: Foreground
        PARAMETERS:
            magnitude: 14.0
    NOISE_1:
        NAME: POISSON_NOISE
        PARAMETERS:
            mean: 2.0
```
In making the dataset, this is the collection of objects you will have at your disposal. 
We've created two different `GALAXY`s, three different `POINTSOURCE`s, and a single `NOISE` source. 
Now we'll see how to specify the geometry of the systems that these objects will constitute.

### The GEOMETRY Section
This is the final required section in the config file, and will contain information about the relative orientations of the objects in the `SPECIES` section.
The geometry section is organized by the user specifying different `CONFIGURATION`s, which will follow the same indexing rules as the objects in the `SPECIES` section.

Each configuration will have a `NAME`, a `FRACTION`, and then informataion about the physical systems. 
The `NAME` will be used to track the association of properties within `deeplenstronomy`. 
The `FRACTION` is how much of the dataset you would like to be in the specific configuration.
Internally, `deeplenstronomy` just multiplies the `FRACTION` by the `SIZE` in the `DATASET` section to obtain a number of images to simulate, but it is necessary that the `FRACTION` values across your `CONFIGURATION`s sum to one.

Next, you will defind the `PLANE`s in your physical systems. Let's look at an example to work through this.

```
GEOMETRY:
    CONFIGURATION_1:
        NAME: GALAXY_AGN
        FRACTION: 0.25
        PLANE_1:
            OBJECT_1: LENS
            PARAMETERS:
                REDSHIFT: 0.2                    
        PLANE_2:
            OBJECT_1: SOURCE
            OBJECT_2: AGN
            PARAMETERS:
                REDSHIFT: 0.7                  
        NOISE_SOURCE_1: POISSON_NOISE
```                   

After the `NAME` and `FRACTION` attributes, we can see a `PLANE_1` and a `PLANE_2` have been defined. 
`PLANE`s require you to list the objects you want in the plane and to define the redshift as a parameter.
The defining characteristic of a `PLANE` in `deeplenstronomy` is that all objects place in that plane are assigned the same redshift.

Here, we've chosen to put a single object in `PLANE_1`, which is the plane closest to Earth. 
We've specified which `SPECIES` to serve as `OBJECT_1` by using the name "LENS", which recall was the name of `GALAXY_1` in the `SPECIES` section.
We've also set this plane to have a redshift of 0.2.

In `PLANE_2`, we've place two objects: "SOURCE" which refers to `GALAXY_2` and "AGN" which refers to `POINTSOURCE_1`. 
We've set this plane to have redshift 0.7 so that it will be behind the objects in `PLANE_1`.

Lastly, we specify that we want to simulate extra noise on top of the image by specifying a single noise source: `NOISE_SOURCE_1`.
We tell `deeplenstronomy` that we want "POISSON\_NOISE" to be used which will reference `NOISE_1` in the `SPECIES` section since it was named "POISSON\_NOISE".

Now. 
That seems like a lot, but let's say you also want to simulate a set of images with no noise to evaluate the performance of some neural network on high quality images. 
All you need to change to the config file is to add another configuration:

```
<Still inside GEOMETRY section>
    CONFIGURATION_2:
        NAME: GALAXY_AGN_NOISELESS
        FRACTION: 0.25
        PLANE_1:
            OBJECT_1: LENS
            PARAMETERS:
                REDSHIFT: 0.2                    
        PLANE_2:
            OBJECT_1: SOURCE
            OBJECT_2: AGN
            PARAMETERS:
                REDSHIFT: 0.7                  
```

What if we want to simulate lensed supernovae instead of lensed AGN?
Just add another configuration:

```
<Still inside GEOMETRY section>
    CONFIGURATION_3:
        NAME: LENSED_SNE
        FRACTION: 0.25
        PLANE_1:
            OBJECT_1: LENS
            PARAMETERS:
                REDSHIFT: 0.2                    
        PLANE_2:
            OBJECT_1: SOURCE
            OBJECT_2: SUPERNOVA
            PARAMETERS:
                REDSHIFT: 0.7      
        NOISE_SOURCE_1: POISSON_NOISE        
```

What if you're looking to add some spice to your life, and you want to include in this dataset a double source plane system, with a supernova and AGN in the source galaxy, two foreground stars, and Poisson noise?
Just add another configuration:

```
<Still inside GEOMETRY section>
    CONFIGURATION_4:
        NAME: SPICY_LIFE
        FRACTION: 0.25
        PLANE_1:
            OBJECT_1: LENS
            OBJECT_2: STAR
            OBJECT_3: STAR
            PARAMETERS:
                REDSHIFT: 0.2                    
        PLANE_2:
            OBJECT_1: LENS
            PARAMETERS:
                REDSHIFT: 0.7
        PLANE_3:
            OBJECT_1: SOURCE
            OBJECT_2: SUPERNOVA
            OBJECT_3: AGN
            PARAMETERS:
                REDSHIFT: 1.3      
        NOISE_SOURCE_1: POISSON_NOISE        
```
Whatever you want to do, the answer will likely be "Just add another configuration."

One of the main ideas behind `deeplenstronomy` is that once you have defined all the properties of the objects, survey, and images, you are free to simulate any physical system you want.

### Summary 
That's the basics of how to organize the main config file. You will need `DATASET`, `COSMOLOGY`, `IMAGE`, `SURVEY`, `SPECIES`, and `GEOMETRY` sections.
These sections make it possible to separate different aspects of the simulation, hopefully in a way that seems logical to astronomers.

Next we'll look at how to draw parameter values from distributions instead of hard-coding fixed values into the main config file.

## Drawing from Probability Distributions

It doesn't do much good to have to specify all the light profile, mass profile, redshift, survey, image (and many more) parameters directly in the config file. 
You can, but often you will want to draw those parameters from a probability distribution. 
This action is communicated to `deeplenstronomy` by using the `DISTRIBUTION` keyword in place of the actual parameter value.

As an example, maybe you want to have a range of magnitudes for an object.
Furthermore, let's say you want this range to be a uniform distribution from 16.0 to 23.0.
In the main config file, replace

```
magnitude: 20.0
```

with

```
magnitude:
    DISTRIBUTION:
        NAME: uniform
        PARAMETERS:
            minimum: 16.0
            maximum: 23.0
```

and `deeplenstronomy` will do the rest.
This approach will work for any parameter.

### Available Distributions

The name "uniform" targets a specific function `uniform()` in the file `deeplenstronomy/distributions.py`.
This file contains several distributions that can also be utilized by specifying their name.
The main functions are:
- `uniform` with parameters `minimum` and `maximum
- `uniform_int` with parameters `minimum` and `maximum` for ensuring that the drawn parameter is an integer
- `normal` with parameters `mean` and `std`
- `lognormal` with parameters `mean` and `sigma`
- `delta_function` with parameters `value`
- `symmetric_uniform_annulus` with parameters `r1` and `r2` for drawing values on the interval [-r2,-r1]U[r1,r2]

You are welcome to edit the source code to add any distribution you need, or open a pull request if you think it will be useful for others.

## Advanced Features

Above is all the information you need to begin constructing configuration files and building `deeplenstronomy` datasets. 
To sculpt your simulated datasets with a specific science case in mind, you may find some of these advanced features useful.

### Forcing Correlations

Here we have looked at how to draw parameter values from one-dimensional distributions.
In practice, you may find it useful to enforce physically-motivated correlations, say between an object's redshift and magnitude, or have color-dependent properties of an object.
In either case,  this can be accomplished using "User-Specified Distributions".
In the configuration file, you can add an entry (at the same level as `IMAGE`, `SURVEY`, `GEOMETRY`, etc.) like this:
```
DISTRIBUTIONS:
    USERDIST_1:
        FILENAME: distribution_file.txt
        MODE: interpolate
```
The properties of the file `distribution_file.txt` are described in the "UserDistributions" Notebook.

### Including Personal Images

You may find it useful to include real astronomical images in your simulations, either to have real-data like noise or to use actual object pictures in your dataset.
In the configuration file, you can add an entry (at the same level as `IMAGE`, `SURVEY`, `GEOMETRY`, etc.) like this:
```
BACKGROUNDS: 
    PATH: <path/to/image_directory_name> #(no trailing '/')
    CONFIGURATIONS: <configuration list, e.g. ['CONFIGURATION_1'] or ['CONFIGURATION_1', 'CONFIGURATION_3']>
```
The `PATH` and `CONFIGURATIONS` keys are described in detail in the "UserImages" Notebook.

### Supplemental Input Files

As you may have noticed by now, this main config file can get pretty long.
To make staying organized when working in this framework a little easier, try out the `INPUT` keyword.

As an example, let's clean up the `SURVEY` section.
This section currently looks like this:

```
SURVEY:
    PARAMETERS:
        BANDS: g,r,i,z,Y
        seeing: 0.9
        magnitude_zero_point: 30.0
        sky_brightness: 27.5
        num_exposures: 10
```

Copy the parameters subsection into a new file, which can be named whatever you want, but we'll name the file `survey.yaml`.
The file will look like this:

```
PARAMETERS:
    BANDS: g,r,i,z,Y
    seeing: 0.9
    magnitude_zero_point: 30.0
    sky_brightness: 27.5
    num_exposures: 10
```

Now in the main config file, change the `SURVEY` section to this:

```
SURVEY:
    INPUT: survey.yaml
```

Thus, the `INPUT` keyword will substitue the contents of whatever file it is used on in place of itself.
Hopefully, you will find this tool useful in modularizing and organizing your `deeplenstronomy` workflow.


## Troubleshooting

This software is new and certainly breakable. 
In an ideal world, the only time you will encounter an error will be typos in the configuration file. 
`deeplenstronomy` automatically parses your configuration file for areas that are entered incorrectly and will alert you before simulating your dataset. 
However, the checks it performs are by no means exhaustive. 
If you do encounter an error, please open an issue so that it can be addressed in the list of checks `deeplenstronomy` performs at runtime.


# Notebooks

Start by reading the [Getting Started Guide](https://deepskies.github.io/deeplenstronomy/Notebooks/GettingStarted.html) to familiarize yourself with the `deeplenstronomy` style.

After that, check out the example notebooks below:

### Notebooks for `deeplenstronomy` Utilities
- [Creating `deeplenstronomy` Configuration Files](https://deepskies.github.io/deeplenstronomy/Notebooks/ConfigFiles.html)
- [Generating Datasets](https://deepskies.github.io/deeplenstronomy/Notebooks/DeepLenstronomyDemo.html)
- [Visualizing `deeplenstronomy` Images](https://deepskies.github.io/deeplenstronomy/Notebooks/Visualization.html)
- [Utilizing Astronomical Surveys](https://deepskies.github.io/deeplenstronomy/Notebooks/Surveys.html)
- [Defining Your Own Probability Distributions](https://deepskies.github.io/deeplenstronomy/Notebooks/UserDistributions.html)
- [Using Your Own Images as Backgrounds](https://deepskies.github.io/deeplenstronomy/Notebooks/BackgroundsDemo.html)
- [Simulating Time-Series Datasets](https://deepskies.github.io/deeplenstronomy/Notebooks/TimeSeriesDemo.html)

### Notebooks for Applying `deeplenstronomy` to Machine Learning Analyses
- [Using `deeplenstronomy` for Active Learning](https://deepskies.github.io/deeplenstronomy/Notebooks/ActiveUpdateDemo.html)
- [Using `deeplenstronomy` for Classification and Regression](https://deepskies.github.io/deeplenstronomy/Notebooks/Metrics.html)

### Notebooks for Suggested Science Cases
- [A Walkthrough of Using `deeplenstronomy` for Science](https://deepskies.github.io/deeplenstronomy/Notebooks/FullExample.html)

# Simulator

Scrits to simulate DES-like image cutouts of lensed quasars. 
Quasar color information is based on SNANA simulated AGN light curves.

### Writen by Rob Morgan and Josh Yao-Yu Lin

## Data

All light curves are located in the [google bucket](https://console.cloud.google.com/storage/browser/deepskies-strong-lenses/quasar_sims).
Download and untar the file `lcs_plus_gal_param.tar.gz` into the current working directory.

## Main Interface: run_simulator.py

`run_simulator.py` takes 2 command line arguments:

1. The number of griz images to be simulated. Max is ~115,000 at present.
2. The object class to be simulated. Choose from `['lenses', 'foregrounds', 'galaxies']`.

### Usage Example

To simulate 10,000 double/quad lensed quasars, the command would be:

`python run_simulator.py 10000 lenses`

## Inner Workings: simulator.py

`run_simulator.py` calls the image simulation functions in `simulator.py` and stores all images in an array.
The functions in `simulator.py` were adapted from `SImulation_full_pipeline_double_quad_non_lens.ipynb`.
Changes include passing light curve, psf, and geometric information to from the input data files directly to the image generation.

## Notes:

This code likely needs to be proofed and made more efficient.





============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given. 

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/deepskies/deeplenstronomy/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.
* Your main configuration file

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

deeplenstronomy could always use more documentation, whether as part of the 
official deeplenstronomy docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/deepskies/deeplenstronomy/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `deeplenstronomy` for
local development.

1. Fork_ the `deeplenstronomy` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/deeplenstronomy.git

3. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

Now you can make your changes locally.

4. When you're done making changes, make sure to also add a test.

See the test directory for examples

5. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

6. Submit a pull request through the GitHub website.

.. _Fork: https://github.com/deepskies/deeplenstronomy/fork

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add an
   example demonstrating your feature.
3. The pull request should work for Python >= 3.7, and for PyPy.

Tips
----

To run a subset of tests::

	 $ cd test && pytest
=======
Credits
=======

Development Lead
----------------

* Rob Morgan <robert.morgan@wisc.edu>

Contributors
------------

* Brian Nord <nord@fnal.gov>
* Simon Birrer <sibirrer@gmail.com>
* Jason Poh <jasonpoh@uchicago.edu>
* Joshua Yao-Yu Lin <joshualin24@gmail.com>

.. :changelog:

History
-------

0.0.2.3 (2021-12-20)
+++++++++++++++++++++
* Bug fixes for ITERATE option in image BACKGROUNDS

0.0.2.2 (2021-12-13)
+++++++++++++++++++++
* Bug fixes for ITERATE option in image BACKGROUNDS
* Bug fix for time delay calculations (very minor)
* New PARAMS option to specify USERDIST columns

0.0.2.1 (2021-11-15)
+++++++++++++++++++++
* ITERATE option for image BACKGROUNDS
* Drop any systems with unphysical redshifts from simulations

0.0.2.0 (2021-05-07)
+++++++++++++++++++++
* Verified stability of all new time series features

0.0.1.9 (became version 0.0.2.0)
+++++++++++++++++++++
* Fix bug in number of lensed point sources

* Remove a couple CC SNe SEDs that led to erroneous magnitudes

0.0.1.8 (2021-04-26)
+++++++++++++++++++++
* Fix bug in extrapolating nites outside of SEDs for time series

* Require each galaxy to have at least one mass profile

* Require each configuration to have at least two planes

0.0.1.7 (2021-04-01)
+++++++++++++++++++++
* Write simulation input dicts to disk to limit memory usage

* Improve accuracy of calculated magnitudes from SEDs (SNe, KN)

* Track cosmographic time delays in the metadata

* Make image backgrounds compatible with time series

0.0.1.6 (2021-03-16)
+++++++++++++++++++++
* Fix bug in calculation of K-Correction

* Add DES deep field distributiions
  
0.0.1.5 (2021-03-10)
+++++++++++++++++++++
* Fix bug in the number of times a USERDIST gets sampled

* Fix bug in lsst survey mode

* Fix bug in the redshifting calculations for supernovae to prevent NaNs

0.0.1.4 (2021-03-03)
+++++++++++++++++++++
* Fix bug in checking configuration file geometry section

* Speed improvements for timeseries functionalities

* Corner plot functionality for metadata visualization

0.0.1.3 (2021-02-02)
+++++++++++++++++++++

* Introducing the static model for timeseries

* Introducing the peakshift parameter for timeseries

* More accurate treatment of noise for timeseries

0.0.1.2 (2021-01-29)
+++++++++++++++++++++

* Fix bug in saving both sigma_v and theta_E 

* Full API documentation

0.0.1.0 (2020-11-09)
+++++++++++++++++++++

* First official public release

0.0.0.14 (2020-11-09)
+++++++++++++++++++++

* Bug fixes in distributions

* Unit tests

0.0.0.11 (2020-10-27)
+++++++++++++++++++++

* Bug fixes in image backgrounds

* Random seeds

* Search for dataset parameter names

0.0.0.10 (2020-09-30)
+++++++++++++++++++++

* Beta Release

0.0.0.9 (2020-09-30)
++++++++++++++++++++

* Image Backgrounds

* User Distributions

0.0.0.6 (2020-08-17)
++++++++++++++++++++

* Implement time-series functinalities

0.0.0.1 (2020-01-24)
++++++++++++++++++++

* Rebrand to yaml-style configuration file

0.1.0 (2019-01-03)
++++++++++++++++++

* First release on PyPI.
.. include:: ../CONTRIBUTING.rst
.. include:: ../HISTORY.rst
========
Usage
========

To use deeplenstronomy in a project::

	import deeplenstronomy
.. include:: ../AUTHORS.rst
============
Installation
============

At the command line either via easy_install or pip::

    $ easy_install deeplenstronomy
    $ pip install deeplenstronomy

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv deeplenstronomy
    $ pip install deeplenstronomy
.. complexity documentation master file, created by
   sphinx-quickstart on Tue Jul  9 22:26:36 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: ../README.rst

Contents:
=========

.. toctree::
   :maxdepth: 2

   installation
   usage
   contributing
   authors
   history

Feedback
========

If you have any suggestions or questions about **deeplenstronomy** feel free to email me
at nord@fnal.gov.

If you encounter any errors or problems with **deeplenstronomy**, please let me know!
Open an Issue at the GitHub http://github.com/bnord/deeplenstronomy main repository.
