---
title: 'tacmagic: Positron emission tomography analysis in R'
authors:
- affiliation: "1, 2"
  name: Eric E. Brown
  orcid: 0000-0002-1575-2606
date: "25 January 2019"
bibliography: paper.bib
tags:
- positron emission tomography
- biomedical imaging
- neuroimaging
- neuroscience
- neuroinformatics
affiliations:
- index: 1
  name: Department of Psychiatry and Institute of Medical Science, University of Toronto, Toronto, Canada
- index: 2
  name: Centre for Addiction and Mental Health, Toronto, Canada
---

# Background

Positron emission tomography (PET) is a research and clinical imaging modality that uses radioactive tracers that bind to target molecules of interest. A PET scanner identifies the tracer location by virtue of the tracer's radioactive decay, providing information about the distribution of target in the body.  Analysis pipelines are used to calculate radiotracer activity over time within a spatial region of interest (ROI). The resulting time-activity curves (TAC) are analyzed to answer important clinical and research questions using kinetic models [@Dierckx:2014].

# The ``tacmagic`` R package

By supporting multiple source data formats, ``tacmagic`` provides an open R [@R] platform for the analysis of PET TAC data that has been produced by existing image analysis pipelines. The data loading functions provide a common format for subsequent analysis in R. We have also implemented basic non-invasive models commonly used in PET research [@Lopresti:2005; @Logan:1996], which have been tested against existing tools [@tpcclib]. The goal is to facilitate open, explicit and reproducible research.
 
The major features of ``tacmagic`` are documented in a walkthrough vignette that is included with the package. The features include:
 
1. loading TAC and volume data to analyze in R,
2. merging regional TAC data into larger ROIs weighted by volume,
3. basic TAC plotting,
4. calculation of standardized uptake value ratio (SUVR) [@Lopresti:2005;@Dierckx:2014],
5. calculation and plotting of the non-invasive reference region Logan DVR model [@Logan:1996;@tpcclib] and
6. calculation of cut-off values for dichotomizing data [@Aizenstein:2008].
 
The package is published with an open source licence, enabling future collaboration and expansion of the package's functions, which may include future support for additional data formats, kinetic models, plotting and cut-off calculation.

# Acknowledgments

Many thanks are due to the kind mentorship of Ariel Graff-Guerrero, Philip Gerretsen and Bruce Pollock, as well as to Fernando Caravaggio, Jun Chung, and Tiffany Chow for their guidance in PET analysis techniques prior to the development of this package. 

Development of parts of this package involved work supported by the Canadian Institute for Health Research Canada Graduate Scholarship, the Ontario Graduate Scholarship, and the Clinician Scientist Program of the University of Toronto's Department of Psychiatry.

# References
# Contributing

A primary purpose of tacmagic is to foster open and reproducible in PET neuroimaging research. Therefore, the package is released under a GPL licence and contributions are welcome.

At this relatively early stage in development, it would be extremely helpful to hear feedback about what is working well compared to existing workflows, and what needs work.

Please suggest changes or report bugs using the Issues page of this repository. 

Expanded funcitonality could be considered in the following areas:

1. Additional file format support for loading TAC and related data.
2. Implementation of additional kinetic models
3. Implementation or improvement of plotting functions
4. Implementation or improvmeent of cutoff calculation algorithms

### File formats

To propose additional file format support, a sample file is required for testing. Preferably, the file could be generated from the same PET data used for other formats. If a file format specification document is available, it should be referenced.

### Kinetic models

At a minimum, proposals to expand the supported models should include the R code, a refence describing the model, and the expected results using an exisiting tool for testing purposes. If possible, the existing example tac data would be used if possible.
# tacmagic: PET Analysis in R


[![DOI](https://zenodo.org/badge/131427691.svg)](https://zenodo.org/badge/latestdoi/131427691) [![Build Status](https://travis-ci.org/ropensci/tacmagic.svg?branch=master)](https://travis-ci.org/ropensci/tacmagic) [![Coverage status](https://codecov.io/gh/ropensci/tacmagic/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/tacmagic?branch=master) [![](https://badges.ropensci.org/280_status.svg)](https://github.com/ropensci/software-review/issues/280)
 [![JOSS](http://joss.theoj.org/papers/10.21105/joss.01281/status.svg)](https://doi.org/10.21105/joss.01281)


To foster openness, replicability, and efficiency, `tacmagic` facilitates loading and analysis of positron emission tomography data in R.

As a `tacmagic` is a new package, we strongly recommend checking all work against existing analyses to confirm the results are as expected, and welcome any feedback.

## Installation

The stable version of the package can be installed from CRAN, and the more recent development version can be installed with the devtools package. 

Use the following R commands to download the version you would like: for the CRAN release,  `install.packages("tacmagic")`, for the github release version that may not yet be available on CRAN, `devtools::install_github("ropensci/tacmagic")`, and for the very latest in-development version that is more likely to have bugs or errors and thus is not suitable for production use, use `devtools::install_github("ropensci/tacmagic", ref="devel")`.

## Features

The features of `tacmagic` are demonstrated in the package's walkthrough vignette, which is highly recommended for first-time uses.

### Data loading and weighted-averages

Time-activity curve (TAC) and/or region of interest (ROI) volume data can be loaded from various file formats including [PMOD](https://www.pmod.com/web/) .tac and .voistat files, a .mat file from the [magia](http://aivo.utu.fi/magia/) pipeline, and [Turku PET Centre's](http://turkupetcentre.fi) .DFT format. 

There is support for converting the radioactivity units in TAC data.

This package is not affiliated with any of the above pipelines.

### Time-activity curve plotting

Basic plotting of one or more TAC from one or more participants is available.

### Binding potential models

Non-invasive models are implemented including the standardized uptake volume (SUV), SUV ratio (SUVR), and the non-invasive Logan reference method.

### Batch and group-wise analysis

Loading and analysis functions can be run as a batch or by individual participant.

## Licence

    Copyright (C) 2018 Eric E. Brown

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

We also note specifically that this package is not intended for clinical use, and may contain bugs or errors, so any results should be verified. As above, we provide no warranty and assume no liability.

## Citation

Please cite this software package if you use it in your analyses. 

`Brown, E. E. (2019). tacmagic: PET Analysis in R. Journal of Open Source Software, 4(34), 1281. doi:10.21105/joss.01281.`

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# tacmagic News

## tacmagic 0.3.0 (2019-06-06)

### Major changes

- Radioactivity unit conversion can now be done for tac objects and numeric objects
- Conversion of a tac object to standardized uptake value (SUV) by controlling for weight and injected dose
- Vignette updated with new features

## tacmagic 0.2.1 (2019-03-06)

### Minor changes

This version has no changes to features but addresses minor, mostly stylistic, issues for the initial CRAN release.

## tacmagic 0.2.0 (2019-02-26)

### Major changes

This is the first release as a complete R package, and the version that has undergone open review with rOpenSci. The main features of tacmagic 0.2.0 are implemented and tested: loading PET time activity curve (tac) data from multiple formats, merging ROIs weighted for volume, calculating binding potential models including SUVR and DVR, basic plotting, and calculation of cut-off values. As tacmagic is still a new package and testing to date is based on few example files, it is essential to verify that all functions are behaving as expected when applied to your own data.
