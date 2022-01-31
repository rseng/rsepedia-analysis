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
---
title: "Analysis with tacmagic"
author: "Eric Brown"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{Analysis with tacmagic}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(tacmagic)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Background

[Positron emission tomography](https://en.wikipedia.org/wiki/Positron_emission_tomography) (PET) is a research and clinical imaging modality that uses radioactive tracers that bind to target molecules of interest. A PET scanner identifies the tracer location by virtue of the tracer's radioactive decay, providing information to determine the location of the target in the body. As the spatial resolution of PET is relatively poor, analysis is frequently combined with higher resolution imaging such as magnetic resonance imaging (MRI), which can be spatially co-registered to the PET image. Subsequently, radiotracer activity (over time) can be identified by spatial region of interest (ROI).

An image analysis pipeline is required to extract regional time activity curves (TACs) from a dynamic PET image. There are various pipelines available including widely-used commercial solutions (e.g. [PMOD](https://www.pmod.com/web/)) and newer open-source options (e.g. magia^[Tomi Karjalainen, Severi Santavirta, Tatu Kantonen, Jouni Tuisku, Lauri Tuominen, Jussi Hirvonen, Jarmo Hietala, Juha Rinne, Lauri Nummenmaa. Magia: Robust automated modeling and image processing toolbox for PET neuroinformatics. bioRxiv 604835; <https://doi.org/10.1101/604835>]). Pipelines generally implement the following steps:

* Dynamic PET pre-processing (e.g. motion correction, decay-correction)
* PET image co-registration with structural MRI
* MRI segmentation and normalization to atlas

Various pipelines save TAC, volume and related data in various formats. This package enables the loading and analysis of TAC and ROI volume data from image analysis pipelines for further analysis in R. 

## Vignette data

The sample data for this vignette uses an anonymized scans of a participant with Alzheimer's dementia, data from http://www.gaain.org which was generously made available for unrestricted use.^[Klunk, William E., Robert A. Koeppe, Julie C. Price, Tammie L. Benzinger, Michael D. Devous, William J. Jagust, Keith A. Johnson, et al. “The Centiloid Project: Standardizing Quantitative Amyloid Plaque Estimation by PET.” Alzheimer’s & Dementia 11, no. 1 (January 2015): 1-15.e4. https://doi.org/10.1016/j.jalz.2014.07.003.] The radiotracer used is Pittsburgh Compound B (PIB) which binds to beta-amyloid, a protein found in high concentration in the brains of individuals with Alzheimer's dementia.

There are two approaches to using the **tacmagic** package to analyze PET time-activity curve data: either by loading participant data individually and using the various functions to analyze it, or via the batch functions to list and analyze data from multiple participants. Here, we illustrate the main features of tacmagic, by walking through the analysis of a single participant. We provide an explanation of the batch mode at the end.

## Time-activity curve operations

### Data loading

Time-activity curve (TAC) data is loaded via `load_tac()`, which is a wrapper for format-specific functions. To specify which file format the TAC data is stored as, use the `format` parameter. Supported formats can be viewed in `help(load_tac)`.

The minimal amount of information required is the TAC data for one or more ROI, including the start and stop times of each frame, the time units and the activity units. This information may be in 1 or more files depending on the format and software that created it. 

For example, PMOD's `.tac` files contain all of the information, but the TAC .voistat files do not contain start and stop times, but this information could be specified using a `.acqtimes` file. Support is also available for DFT format, which contains both TAC and volume data.

We processed the PIB PET and T1 MRI data with the **PMOD PNEURO** software suite to produce a `.tac` file with TACs for all ROIs in the Hammer's atlas. The `.tac` file can be loaded with `load_tac()`:

```{r}
# Filename is a character string of the file's path on your computer.
filename <- system.file("extdata", "AD06.tac", package="tacmagic")
# Note: This file can also serve as a template if the TAC data is in some other 
# format that is not yet supported.

AD06_tac <- load_tac(filename, format="PMOD")
```

A TAC object is a data frame with extra attributes including time and activity units. A summary can be printed with the generic `print()` function.

```{r}
summary(AD06_tac) 

AD06_tac[1:5,1:5] # the first 5 frames of the first 3 ROIs
```

PMOD's suite also produces `.voistat` and `.acqtimes` formats, than can be used to produce the same data if you do not have `.tac` files:

```{r}
filename_acq <- system.file("extdata", "AD06.acqtimes", package="tacmagic")
filename_voistat <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")

tac2 <- load_tac(filename_voistat, format="voistat", acqtimes=filename_acq)

all.equal(AD06_tac, tac2)
```

We also used Turku's **magia** pipeline to process the same data. It can be loaded similarly, though with units explicitly entered because the information is not encoded in the .mat file:

```{r}
f_magia <- system.file("extdata", "AD06_tac_magia.mat", package="tacmagic")

AD06_tac_magia <- load_tac(f_magia, format="magia", 
                           time_unit="seconds", activity_unit="kBq/cc")

AD06_tac_magia[1:5,1:5]
```

#### Manually-created TAC objects

For other data sources, **tacmagic** TAC objects can be created from data.frame objects with `as.tac()`. The time and activity units must be specified as arguments if not already set as attributes in the data.frame. The columns of the data.frame are the regional TACs, with the column names the names of the ROIs. 

```{r}
manual <- data.frame(start=c(0:4), end=c(2:6), ROI1=c(10.1:14.2), ROI2=c(11:15))
manual_tac <- as.tac(manual, time_unit="minutes", activity_unit="kBq/cc")

summary(manual_tac)
```

### Radioactivity unit conversion

Most modern PET tools use kBq/cc as the standard activity units required by the software (e.g. TPCCLIB, PMOD). Most often data will be in this format and not require conversion. However, conversion of TAC objects to these units, or to other radioactivity units for that matter, is possible with `change_units()`. This is a generic function that works on `tac` objects as well as `numeric` objects. For `numeric` objects, both "to" and "from" units need to be specified. The function works regardless of whether the units are per volume (i.e. kBq is treated the same was as kBq/cc or kBq/mL).

```{r}
change_units(5, to_unit = "kBq", from_unit = "nCi")
change_units(0.5, to_unit = "nCi/cc", from_unit = "kBq/cc")
```

With `tac` objects, as the activity units are stored in the object, they should not be provided to `change_units()`:

```{r}
AD06_nCi <- change_units(AD06_tac, to_unit = "nCi/cc")
summary(AD06_nCi)
```

### ROI merging

Often it is desirable to combine TAC ROIs into larger ROIs. For example, if the PET analysis pipeline created TACs for each atlas ROI, your analysis may call for merging these atomic ROIs into larger regions, such as merging all of the atlas ROIs that make up the frontal lobe into a single frontal lobe ROI.

If this is done, the means should be weighted for the relative volumes of the atomic ROIs. If volume information is available, `tac_roi()` provides this functionality.

In PMOD's software, volume information is available in `.voistat` files. Units do not matter because it is the relative volume information that is needed.

In addition to TAC and volume information, we must specify which atomic ROIs make up the merged ROI. This is done by providing a named list, where the names are the merged ROIs and the list items are themselves lists of the atomic ROIs that make up each merged ROI. For the Hammer's atlas, and as an example, typical data is provided in `roi_ham_stand()`, `roi_ham_full()`, or `roi_ham_pib()`. 

```{r}
AD06_volume <- load_vol(filename_voistat, format="voistat")

roi_ham_pib()[1:2] # The first 2 definitions of merged ROIs, as an example.

AD06 <- tac_roi(tac=AD06_tac,           # The TAC file we loaded above.
                volumes=AD06_volume,    # Volume information loaded.
                ROI_def=roi_ham_pib(),  # ROI definitions for the Hammers atlas
                merge=F,                # T to also return atomic ROIs
                PVC=F                   # to use _C ROIs (PMOD convention)            
                )

AD06[1:5,1:5]
```

### Plotting

Basic TAC plotting can be done by calling `plot`, which accepts two TAC objects, e.g. from 2 participants or group means. The ROIs to plot are specified as a vector of ROI names as they appear in the TAC object. As the TAC object contains time unit information, the plot can convert to desired units, which can be specified with the `time` argument.

```{r, fig.show='hold', fig.height=4.5, fig.width=6.5, fig.align='center'}
plot(AD06,                                                    # TAC data
     ROIs=c("frontal", "temporal", "parietal", "cerebellum"), # ROIs to plot
     time="minutes",                   # Convert x axis from seconds to minutes
     title="PIB time activity curves for AD06"        # A title for the plot
     )
```


## Model calculation

### Standardized uptake value (SUV)

As the activity in the TAC is impacted by the dose of the radiotracer administered and the participant's body weight, a value adjusted for these factors is sometimes used, the [SUV](http://www.turkupetcentre.net/petanalysis/model_suv.html):

$$SUV = \frac{Ct}{\frac{Dose}{Weight}}$$

Where activity is measured in _kBq/mL (kBq/cc)_, the dose is in _MBq_, and the weight is in _kg_, the radioactivity units cancel and the SUV units are g/mL.

With the `tac_suv()` function, a `tac` object can be converted to SUV values with units _g/mL_, as demonstrated below (note the weight and dose are fabricated here for the demonstration).

```{r}
AD06_suv_tac <- tac_suv(AD06, dose = 8.5, dose_unit = "mCi", weight_kg = 70)
```

More often, an everage value over a certain time period, or a maximum value may be desired. This can be calculated with the `suv()` function.

```{r}
AD06_suv_calc <- suv(AD06, SUV_def = c(3000, 3300, 3600), dose = 8.5, dose_unit = "mCi", weight_kg = 70)
AD06_suv_calc["frontal",]
AD06_suv_max <- suv(AD06, SUV_def = "max", dose = 8.5, dose_unit = "mCi", weight_kg = 70)
AD06_suv_max["frontal",]
```

### SUV ratio (SUVR)

The standardized uptake value ratio ($SUVR$) is a simple quantification of PET activity that is commonly used from many tracers including PIB. It is the ratio of the tracer activity over a specified time period ($Ct$) in a target ROI to a reference region. Using a ratio allows factors that are normally required to calculate an $SUV$ to cancel out, namely tracer dose and patient body weight, and therefore $SUVR$ can be calculated from TAC data alone, i.e. without the need to specify tracer dose or body weight: 

$$SUVR = \frac{SUV_{TARGET}}{SUV_{REF}} = \frac{Ct_{TARGET}}{Ct_{REF}}$$

In the literature, SUVR is variably described and calculated using the mean of activity for the frames of the specified time period, or the area under the curve. For PIB, the mean/summed activity has been used, and the time windows have varied from starting at 40-50 minutes and ending at 60-90 minutes.^[Lopresti, B. J., W. E. Klunk, C. A. Mathis, J. A. Hoge, S. K. Ziolko, X. Lu, C. C. Meltzer, et al. “Simplified Quantification of Pittsburgh Compound B Amyloid Imaging PET Studies: A Comparative Analysis.” J Nucl Med 46 (2005): 1959–72.]

The `suvr()` function calculates SUVR for all regions in a TAC file based on the provided time information (as a vector of frame start times) and the specified reference region (a string). If the frames used are of different durations, the weighted mean is used.

```{r}
AD06_SUVR <- suvr(AD06,                       # TAC data
                  SUVR_def=c(3000,3300,3600), # = 50-70 minute window
                  ref="cerebellum"            # reference region in TAC data
                  )

AD06_SUVR

```
An alternative method, using the area under the curve with the mid-frame times as the x-axis is available with `suvr_auc()` and should provide very similar results.

```{r}
AD06_altSUVR <- suvr_auc(AD06, SUVR_def=c(3000,3300,3600), ref="cerebellum")

all.equal(AD06_SUVR, AD06_altSUVR) # Should be similar but not exact

```

### DVR

The Distribution Volume Ratio (DVR) is a method of quantifying tracer uptake that is used as an alternative to the SUVR in PIB studies, for example. Like SUVR, it can be calculated from TAC data without the need for arterial blood sampling, by making use of a reference region. In this case, it is called the _non-invasive_ Logan plot method. It is calculated with a graphical analysis technique described by Logan et al.^[Logan, J., Fowler, J. S., Volkow, N. D., Wang, G.-J., Ding, Y.-S., & Alexoff, D. L. (1996). Distribution Volume Ratios without Blood Sampling from Graphical Analysis of PET Data. Journal of Cerebral Blood Flow & Metabolism, 16(5), 834-840. https://doi.org/10.1097/00004647-199609000-00008] 

In addition to the TAC data, depending on the tracer, a value for k2' may need to be specified. For PIB, this has limited effect on the result, but can be specified, and a value of 0.2 has been recommended.^[http://www.turkupetcentre.net/petanalysis/analysis_11c-pib.html]

The non-invasive Logan plot works by finding the slope of the line of the following equation after time $t*$ where linearity has been reached:

$$\frac{\int_0^{T}C_{roi}(t)dt}{C_{roi}(t)} = DVR[\frac{\int_0^{T}C_{cer}(t)dt + C_{cer}(t) / k2`}{C_{roi}(T)}] + int  $$

#### Find t*

The time, $t*$ (`t_star`), after which the relationship is linear can be found by testing the point after which the error is below a certain threshold (default is 10%). If `t_star=0`, then tacmagic can find the suitable value.

```{r}
AD06_DVR_fr <- DVR_ref_Logan(AD06, 
                             target="frontal", # target ROI
                             ref="cerebellum", # reference region
                             k2prime=0.2,      # suitable k2' for tracer
                             t_star=0,        # 0 to find, or can specify frame
                             )

AD06_DVR_fr$DVR

``` 

To visually confirm that the model behaved as expected with linearity, there is a plotting function:

```{r, fig.show='hold', fig.height=4.5, fig.width=6.5, fig.align='center'}
plot(AD06_DVR_fr)

``` 

The right plot shows the Logan model, with the vertical line representing the identified $t*$, and the linear model fitted to the points after that time. In this case, the line after $t*$ can be seen to fit well. The slope of that line is the DVR.

Similarly, DVR can be calculated for all ROIs, either by setting `t_star` manually or to 0 as before. If 0, a different value will be identified for each ROI.

```{r}
AD06_DVR <- DVR_all_ref_Logan(AD06, ref="cerebellum", k2prime=0.2, t_star=23)

AD06_DVR

```

For this data, the DVR calculation has been shown to produce equivalent results as an existing tool.^[https://gitlab.utu.fi/vesoik/tpcclib] 

A wrapper function `dvr()` is available to conveniently calculate DVR for a target ROI or all ROIs, and currently defaults to using the Logan reference method:

```{r}
ADO6_frontal_DVR <- dvr(AD06, target="frontal", ref="cerebellum", k2prime=0.2, 
                        t_star=23)

```

## Batch analysis

In most cases, a project will involve the analysis of multiple participants. The above workflow can be used to test and visualize an analysis, but a batch workflow will likely be preferred to analyze multiple participants.

All analyses can be run using 2 steps: a batch data loading step and a batch analysis step. 

### Batch loading

Data loading is done by `batch_load()`. See `help(batch_load)` for the required arguments. 

The first argument is a vector of participant IDs that corresponds to file names, e.g.:

`participants <- c("participant01", "participant02")` if the files are located e.g. `/mypath/participant01.tac` and `/mypath/participant01_TAC.voistat`. In this case, the function call might look like: 

`my_data <- batch_load(participants, dir="/mypath/", tac_format="PMOD", roi_m=T, vol_file_suffix="_TAC.voistat", vol_format="voistat", ROI_def=roi_ham_stand(), merge=F)`

The above would load the appropriate TAC and voistat files, perform the ROI merging specified by `ROI_def`, because `roi_m = TRUE`, and would return a list where each element represents a participants, e.g. the first participant would be `my_data$participant1`.

To calculate SUV in batch, the participants dose and weight must be specified when loading with `batch_load()`, as it is then added to the respective `tac` objects.

### Batch analysis

Once the TAC data is loaded, all analyses can be run using `batch_tm()`. The output from `batch_load()` is the first argument for `batch_tm()`. The models implemented in tacmagic can be specified using the `models` argument, e.g. `models = c("SUVR", "Logan")` to calculate both SUVR and Logan DVR. The relevant model parameters will also need to be specified, so see `help(batch_tm)` for all possible arguments.

### Batch example

For the purpose of the vignette, the list of participants will be a list of the full TAC filenames (hence tac_file_suffix=""). In real-world data, the participants parameter can be a list of participant IDs that correspond to the actual filenames, i.e. the filename is made up of dir + participant + tac_file_suffix.

We will also choose not to use the roi_m option in batch_load(), which could be used to combine ROIs as outlined above.

```{r}

participants <- c(system.file("extdata", "AD06.tac", package="tacmagic"),
                   system.file("extdata", "AD07.tac", package="tacmagic"),
                   system.file("extdata", "AD08.tac", package="tacmagic"))

tacs <- batch_load(participants, dir="", tac_file_suffix="")

# Since the PMOD TAC files used here have 2 copies of ROIs, with and without 
# PVC, we can use split_pvc to keep the PVC-corrected verions. If we had used 
# roi_m here to combine ROIs, we could have specified to use the PVC versions 
# in batch_load() with PVC = TRUE.
tacs <- lapply(tacs, split_pvc, PVC=TRUE)
 
batch <- batch_tm(tacs, models=c("SUVR", "Logan"), ref="Cerebellum_r_C",
                  SUVR_def=c(3000,3300,3600), k2prime=0.2, t_star=23)

```

## Cut-off calculations

In the analysis of PIB/amyloid PET data, often researchers want to dichotomize patients into PIB+ vs. PIB-, i.e. to identify those with significant AD-related amyloid pathology (PIB+).

There are a number of approaches to this depending on the available data. We have implemented a method described by Aizenstein et al.^[Aizenstein HJ, Nebes RD, Saxton JA, et al. 2008. Frequent amyloid deposition without significant cognitive impairment among the elderly. Arch Neurol 65: 1509-1517.] which uses a group of participants with normal cognition to establish a cutoff value above which participants are unlikely to have minimal amyloid pathology.

The method identifies a group of participants out of the normal cognition group with higher-PIB outliers removed. An outlier is a participant with any ROI with a DVR higher than the upper inner fence, from a set of ROIs known to be associated with amyloid deposition. Such participants are removed from the group, and this process is done iteratively until no more outliers are removed. Then, cutoff values are determined from this new group for each ROI, set again as the upper inner fence. Then these cutoff values are applied to all participants, and a participant is deemed PIB+ if they have at least 1 ROI above its cutoff.

To demonstrate, a fake dataset of DVR values for 50 fake participants was generated and is available as `fake_DVR`. This would be equivalent to using `batch_tm()` on a group of participants with the `"Logan"` model specified.

```{r}
fake_DVR[1:5,]

```

To calculate the cutoff values using this iterative method, `cutoff_aiz()` takes 2 arguments: the DVR data, and the names of the variables of the ROI DVRs to use (and there must be at least 2 for this method).

```{r}
cutoffs <- cutoff_aiz(fake_DVR, c("ROI1_DVR", "ROI2_DVR", "ROI3_DVR", "ROI4_DVR"))

cutoffs

```

The final step is to apply the cutoffs to the full set of participants. We will use the same sample data:

```{r}
positivity_table <- pos_anyroi(fake_DVR, cutoffs)

positivity_table

```

The algorithm identified 11 PIB+ participants. In the generation of the sample data, the DVRs from the first 10 participants were drawn from a normal distribution with mean 1.9, sd 0.6 and for the latter 40 participants, from mean 1.3, sd 0.3; thus this pattern is in line with what we would expect: all 10 of the first participants are PIB+, and just 1 of the latter 40 was (by chance).
---
title: "Analysis with tacmagic"
author: "Eric Brown"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{Analysis with tacmagic}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(tacmagic)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Background

[Positron emission tomography](https://en.wikipedia.org/wiki/Positron_emission_tomography) (PET) is a research and clinical imaging modality that uses radioactive tracers that bind to target molecules of interest. A PET scanner identifies the tracer location by virtue of the tracer's radioactive decay, providing information to determine the location of the target in the body. As the spatial resolution of PET is relatively poor, analysis is frequently combined with higher resolution imaging such as magnetic resonance imaging (MRI), which can be spatially co-registered to the PET image. Subsequently, radiotracer activity (over time) can be identified by spatial region of interest (ROI).

An image analysis pipeline is required to extract regional time activity curves (TACs) from a dynamic PET image. There are various pipelines available including widely-used commercial solutions (e.g. [PMOD](https://www.pmod.com/web/)) and newer open-source options (e.g. magia^[Tomi Karjalainen, Severi Santavirta, Tatu Kantonen, Jouni Tuisku, Lauri Tuominen, Jussi Hirvonen, Jarmo Hietala, Juha Rinne, Lauri Nummenmaa. Magia: Robust automated modeling and image processing toolbox for PET neuroinformatics. bioRxiv 604835; <https://doi.org/10.1101/604835>]). Pipelines generally implement the following steps:

* Dynamic PET pre-processing (e.g. motion correction, decay-correction)
* PET image co-registration with structural MRI
* MRI segmentation and normalization to atlas

Various pipelines save TAC, volume and related data in various formats. This package enables the loading and analysis of TAC and ROI volume data from image analysis pipelines for further analysis in R. 

## Vignette data

The sample data for this vignette uses an anonymized scans of a participant with Alzheimer's dementia, data from http://www.gaain.org which was generously made available for unrestricted use.^[Klunk, William E., Robert A. Koeppe, Julie C. Price, Tammie L. Benzinger, Michael D. Devous, William J. Jagust, Keith A. Johnson, et al. “The Centiloid Project: Standardizing Quantitative Amyloid Plaque Estimation by PET.” Alzheimer’s & Dementia 11, no. 1 (January 2015): 1-15.e4. https://doi.org/10.1016/j.jalz.2014.07.003.] The radiotracer used is Pittsburgh Compound B (PIB) which binds to beta-amyloid, a protein found in high concentration in the brains of individuals with Alzheimer's dementia.

There are two approaches to using the **tacmagic** package to analyze PET time-activity curve data: either by loading participant data individually and using the various functions to analyze it, or via the batch functions to list and analyze data from multiple participants. Here, we illustrate the main features of tacmagic, by walking through the analysis of a single participant. We provide an explanation of the batch mode at the end.

## Time-activity curve operations

### Data loading

Time-activity curve (TAC) data is loaded via `load_tac()`, which is a wrapper for format-specific functions. To specify which file format the TAC data is stored as, use the `format` parameter. Supported formats can be viewed in `help(load_tac)`.

The minimal amount of information required is the TAC data for one or more ROI, including the start and stop times of each frame, the time units and the activity units. This information may be in 1 or more files depending on the format and software that created it. 

For example, PMOD's `.tac` files contain all of the information, but the TAC .voistat files do not contain start and stop times, but this information could be specified using a `.acqtimes` file. Support is also available for DFT format, which contains both TAC and volume data.

We processed the PIB PET and T1 MRI data with the **PMOD PNEURO** software suite to produce a `.tac` file with TACs for all ROIs in the Hammer's atlas. The `.tac` file can be loaded with `load_tac()`:

```{r}
# Filename is a character string of the file's path on your computer.
filename <- system.file("extdata", "AD06.tac", package="tacmagic")
# Note: This file can also serve as a template if the TAC data is in some other 
# format that is not yet supported.

AD06_tac <- load_tac(filename, format="PMOD")
```

A TAC object is a data frame with extra attributes including time and activity units. A summary can be printed with the generic `print()` function.

```{r}
summary(AD06_tac) 

AD06_tac[1:5,1:5] # the first 5 frames of the first 3 ROIs
```

PMOD's suite also produces `.voistat` and `.acqtimes` formats, than can be used to produce the same data if you do not have `.tac` files:

```{r}
filename_acq <- system.file("extdata", "AD06.acqtimes", package="tacmagic")
filename_voistat <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")

tac2 <- load_tac(filename_voistat, format="voistat", acqtimes=filename_acq)

all.equal(AD06_tac, tac2)
```

We also used Turku's **magia** pipeline to process the same data. It can be loaded similarly, though with units explicitly entered because the information is not encoded in the .mat file:

```{r}
f_magia <- system.file("extdata", "AD06_tac_magia.mat", package="tacmagic")

AD06_tac_magia <- load_tac(f_magia, format="magia", 
                           time_unit="seconds", activity_unit="kBq/cc")

AD06_tac_magia[1:5,1:5]
```

#### Manually-created TAC objects

For other data sources, **tacmagic** TAC objects can be created from data.frame objects with `as.tac()`. The time and activity units must be specified as arguments if not already set as attributes in the data.frame. The columns of the data.frame are the regional TACs, with the column names the names of the ROIs. 

```{r}
manual <- data.frame(start=c(0:4), end=c(2:6), ROI1=c(10.1:14.2), ROI2=c(11:15))
manual_tac <- as.tac(manual, time_unit="minutes", activity_unit="kBq/cc")

summary(manual_tac)
```

### Radioactivity unit conversion

Most modern PET tools use kBq/cc as the standard activity units required by the software (e.g. TPCCLIB, PMOD). Most often data will be in this format and not require conversion. However, conversion of TAC objects to these units, or to other radioactivity units for that matter, is possible with `change_units()`. This is a generic function that works on `tac` objects as well as `numeric` objects. For `numeric` objects, both "to" and "from" units need to be specified. The function works regardless of whether the units are per volume (i.e. kBq is treated the same was as kBq/cc or kBq/mL).

```{r}
change_units(5, to_unit = "kBq", from_unit = "nCi")
change_units(0.5, to_unit = "nCi/cc", from_unit = "kBq/cc")
```

With `tac` objects, as the activity units are stored in the object, they should not be provided to `change_units()`:

```{r}
AD06_nCi <- change_units(AD06_tac, to_unit = "nCi/cc")
summary(AD06_nCi)
```

### ROI merging

Often it is desirable to combine TAC ROIs into larger ROIs. For example, if the PET analysis pipeline created TACs for each atlas ROI, your analysis may call for merging these atomic ROIs into larger regions, such as merging all of the atlas ROIs that make up the frontal lobe into a single frontal lobe ROI.

If this is done, the means should be weighted for the relative volumes of the atomic ROIs. If volume information is available, `tac_roi()` provides this functionality.

In PMOD's software, volume information is available in `.voistat` files. Units do not matter because it is the relative volume information that is needed.

In addition to TAC and volume information, we must specify which atomic ROIs make up the merged ROI. This is done by providing a named list, where the names are the merged ROIs and the list items are themselves lists of the atomic ROIs that make up each merged ROI. For the Hammer's atlas, and as an example, typical data is provided in `roi_ham_stand()`, `roi_ham_full()`, or `roi_ham_pib()`. 

```{r}
AD06_volume <- load_vol(filename_voistat, format="voistat")

roi_ham_pib()[1:2] # The first 2 definitions of merged ROIs, as an example.

AD06 <- tac_roi(tac=AD06_tac,           # The TAC file we loaded above.
                volumes=AD06_volume,    # Volume information loaded.
                ROI_def=roi_ham_pib(),  # ROI definitions for the Hammers atlas
                merge=F,                # T to also return atomic ROIs
                PVC=F                   # to use _C ROIs (PMOD convention)            
                )

AD06[1:5,1:5]
```

### Plotting

Basic TAC plotting can be done by calling `plot`, which accepts two TAC objects, e.g. from 2 participants or group means. The ROIs to plot are specified as a vector of ROI names as they appear in the TAC object. As the TAC object contains time unit information, the plot can convert to desired units, which can be specified with the `time` argument.

```{r, fig.show='hold', fig.height=4.5, fig.width=6.5, fig.align='center'}
plot(AD06,                                                    # TAC data
     ROIs=c("frontal", "temporal", "parietal", "cerebellum"), # ROIs to plot
     time="minutes",                   # Convert x axis from seconds to minutes
     title="PIB time activity curves for AD06"        # A title for the plot
     )
```


## Model calculation

### Standardized uptake value (SUV)

As the activity in the TAC is impacted by the dose of the radiotracer administered and the participant's body weight, a value adjusted for these factors is sometimes used, the [SUV](http://www.turkupetcentre.net/petanalysis/model_suv.html):

$$SUV = \frac{Ct}{\frac{Dose}{Weight}}$$

Where activity is measured in _kBq/mL (kBq/cc)_, the dose is in _MBq_, and the weight is in _kg_, the radioactivity units cancel and the SUV units are g/mL.

With the `tac_suv()` function, a `tac` object can be converted to SUV values with units _g/mL_, as demonstrated below (note the weight and dose are fabricated here for the demonstration).

```{r}
AD06_suv_tac <- tac_suv(AD06, dose = 8.5, dose_unit = "mCi", weight_kg = 70)
```

More often, an everage value over a certain time period, or a maximum value may be desired. This can be calculated with the `suv()` function.

```{r}
AD06_suv_calc <- suv(AD06, SUV_def = c(3000, 3300, 3600), dose = 8.5, dose_unit = "mCi", weight_kg = 70)
AD06_suv_calc["frontal",]
AD06_suv_max <- suv(AD06, SUV_def = "max", dose = 8.5, dose_unit = "mCi", weight_kg = 70)
AD06_suv_max["frontal",]
```

### SUV ratio (SUVR)

The standardized uptake value ratio ($SUVR$) is a simple quantification of PET activity that is commonly used from many tracers including PIB. It is the ratio of the tracer activity over a specified time period ($Ct$) in a target ROI to a reference region. Using a ratio allows factors that are normally required to calculate an $SUV$ to cancel out, namely tracer dose and patient body weight, and therefore $SUVR$ can be calculated from TAC data alone, i.e. without the need to specify tracer dose or body weight: 

$$SUVR = \frac{SUV_{TARGET}}{SUV_{REF}} = \frac{Ct_{TARGET}}{Ct_{REF}}$$

In the literature, SUVR is variably described and calculated using the mean of activity for the frames of the specified time period, or the area under the curve. For PIB, the mean/summed activity has been used, and the time windows have varied from starting at 40-50 minutes and ending at 60-90 minutes.^[Lopresti, B. J., W. E. Klunk, C. A. Mathis, J. A. Hoge, S. K. Ziolko, X. Lu, C. C. Meltzer, et al. “Simplified Quantification of Pittsburgh Compound B Amyloid Imaging PET Studies: A Comparative Analysis.” J Nucl Med 46 (2005): 1959–72.]

The `suvr()` function calculates SUVR for all regions in a TAC file based on the provided time information (as a vector of frame start times) and the specified reference region (a string). If the frames used are of different durations, the weighted mean is used.

```{r}
AD06_SUVR <- suvr(AD06,                       # TAC data
                  SUVR_def=c(3000,3300,3600), # = 50-70 minute window
                  ref="cerebellum"            # reference region in TAC data
                  )

AD06_SUVR

```
An alternative method, using the area under the curve with the mid-frame times as the x-axis is available with `suvr_auc()` and should provide very similar results.

```{r}
AD06_altSUVR <- suvr_auc(AD06, SUVR_def=c(3000,3300,3600), ref="cerebellum")

all.equal(AD06_SUVR, AD06_altSUVR) # Should be similar but not exact

```

### DVR

The Distribution Volume Ratio (DVR) is a method of quantifying tracer uptake that is used as an alternative to the SUVR in PIB studies, for example. Like SUVR, it can be calculated from TAC data without the need for arterial blood sampling, by making use of a reference region. In this case, it is called the _non-invasive_ Logan plot method. It is calculated with a graphical analysis technique described by Logan et al.^[Logan, J., Fowler, J. S., Volkow, N. D., Wang, G.-J., Ding, Y.-S., & Alexoff, D. L. (1996). Distribution Volume Ratios without Blood Sampling from Graphical Analysis of PET Data. Journal of Cerebral Blood Flow & Metabolism, 16(5), 834-840. https://doi.org/10.1097/00004647-199609000-00008] 

In addition to the TAC data, depending on the tracer, a value for k2' may need to be specified. For PIB, this has limited effect on the result, but can be specified, and a value of 0.2 has been recommended.^[http://www.turkupetcentre.net/petanalysis/analysis_11c-pib.html]

The non-invasive Logan plot works by finding the slope of the line of the following equation after time $t*$ where linearity has been reached:

$$\frac{\int_0^{T}C_{roi}(t)dt}{C_{roi}(t)} = DVR[\frac{\int_0^{T}C_{cer}(t)dt + C_{cer}(t) / k2`}{C_{roi}(T)}] + int  $$

#### Find t*

The time, $t*$ (`t_star`), after which the relationship is linear can be found by testing the point after which the error is below a certain threshold (default is 10%). If `t_star=0`, then tacmagic can find the suitable value.

```{r}
AD06_DVR_fr <- DVR_ref_Logan(AD06, 
                             target="frontal", # target ROI
                             ref="cerebellum", # reference region
                             k2prime=0.2,      # suitable k2' for tracer
                             t_star=0,        # 0 to find, or can specify frame
                             )

AD06_DVR_fr$DVR

``` 

To visually confirm that the model behaved as expected with linearity, there is a plotting function:

```{r, fig.show='hold', fig.height=4.5, fig.width=6.5, fig.align='center'}
plot(AD06_DVR_fr)

``` 

The right plot shows the Logan model, with the vertical line representing the identified $t*$, and the linear model fitted to the points after that time. In this case, the line after $t*$ can be seen to fit well. The slope of that line is the DVR.

Similarly, DVR can be calculated for all ROIs, either by setting `t_star` manually or to 0 as before. If 0, a different value will be identified for each ROI.

```{r}
AD06_DVR <- DVR_all_ref_Logan(AD06, ref="cerebellum", k2prime=0.2, t_star=23)

AD06_DVR

```

For this data, the DVR calculation has been shown to produce equivalent results as an existing tool.^[https://gitlab.utu.fi/vesoik/tpcclib] 

A wrapper function `dvr()` is available to conveniently calculate DVR for a target ROI or all ROIs, and currently defaults to using the Logan reference method:

```{r}
ADO6_frontal_DVR <- dvr(AD06, target="frontal", ref="cerebellum", k2prime=0.2, 
                        t_star=23)

```

## Batch analysis

In most cases, a project will involve the analysis of multiple participants. The above workflow can be used to test and visualize an analysis, but a batch workflow will likely be preferred to analyze multiple participants.

All analyses can be run using 2 steps: a batch data loading step and a batch analysis step. 

### Batch loading

Data loading is done by `batch_load()`. See `help(batch_load)` for the required arguments. 

The first argument is a vector of participant IDs that corresponds to file names, e.g.:

`participants <- c("participant01", "participant02")` if the files are located e.g. `/mypath/participant01.tac` and `/mypath/participant01_TAC.voistat`. In this case, the function call might look like: 

`my_data <- batch_load(participants, dir="/mypath/", tac_format="PMOD", roi_m=T, vol_file_suffix="_TAC.voistat", vol_format="voistat", ROI_def=roi_ham_stand(), merge=F)`

The above would load the appropriate TAC and voistat files, perform the ROI merging specified by `ROI_def`, because `roi_m = TRUE`, and would return a list where each element represents a participants, e.g. the first participant would be `my_data$participant1`.

To calculate SUV in batch, the participants dose and weight must be specified when loading with `batch_load()`, as it is then added to the respective `tac` objects.

### Batch analysis

Once the TAC data is loaded, all analyses can be run using `batch_tm()`. The output from `batch_load()` is the first argument for `batch_tm()`. The models implemented in tacmagic can be specified using the `models` argument, e.g. `models = c("SUVR", "Logan")` to calculate both SUVR and Logan DVR. The relevant model parameters will also need to be specified, so see `help(batch_tm)` for all possible arguments.

### Batch example

For the purpose of the vignette, the list of participants will be a list of the full TAC filenames (hence tac_file_suffix=""). In real-world data, the participants parameter can be a list of participant IDs that correspond to the actual filenames, i.e. the filename is made up of dir + participant + tac_file_suffix.

We will also choose not to use the roi_m option in batch_load(), which could be used to combine ROIs as outlined above.

```{r}

participants <- c(system.file("extdata", "AD06.tac", package="tacmagic"),
                   system.file("extdata", "AD07.tac", package="tacmagic"),
                   system.file("extdata", "AD08.tac", package="tacmagic"))

tacs <- batch_load(participants, dir="", tac_file_suffix="")

# Since the PMOD TAC files used here have 2 copies of ROIs, with and without 
# PVC, we can use split_pvc to keep the PVC-corrected verions. If we had used 
# roi_m here to combine ROIs, we could have specified to use the PVC versions 
# in batch_load() with PVC = TRUE.
tacs <- lapply(tacs, split_pvc, PVC=TRUE)
 
batch <- batch_tm(tacs, models=c("SUVR", "Logan"), ref="Cerebellum_r_C",
                  SUVR_def=c(3000,3300,3600), k2prime=0.2, t_star=23)

```

## Cut-off calculations

In the analysis of PIB/amyloid PET data, often researchers want to dichotomize patients into PIB+ vs. PIB-, i.e. to identify those with significant AD-related amyloid pathology (PIB+).

There are a number of approaches to this depending on the available data. We have implemented a method described by Aizenstein et al.^[Aizenstein HJ, Nebes RD, Saxton JA, et al. 2008. Frequent amyloid deposition without significant cognitive impairment among the elderly. Arch Neurol 65: 1509-1517.] which uses a group of participants with normal cognition to establish a cutoff value above which participants are unlikely to have minimal amyloid pathology.

The method identifies a group of participants out of the normal cognition group with higher-PIB outliers removed. An outlier is a participant with any ROI with a DVR higher than the upper inner fence, from a set of ROIs known to be associated with amyloid deposition. Such participants are removed from the group, and this process is done iteratively until no more outliers are removed. Then, cutoff values are determined from this new group for each ROI, set again as the upper inner fence. Then these cutoff values are applied to all participants, and a participant is deemed PIB+ if they have at least 1 ROI above its cutoff.

To demonstrate, a fake dataset of DVR values for 50 fake participants was generated and is available as `fake_DVR`. This would be equivalent to using `batch_tm()` on a group of participants with the `"Logan"` model specified.

```{r}
fake_DVR[1:5,]

```

To calculate the cutoff values using this iterative method, `cutoff_aiz()` takes 2 arguments: the DVR data, and the names of the variables of the ROI DVRs to use (and there must be at least 2 for this method).

```{r}
cutoffs <- cutoff_aiz(fake_DVR, c("ROI1_DVR", "ROI2_DVR", "ROI3_DVR", "ROI4_DVR"))

cutoffs

```

The final step is to apply the cutoffs to the full set of participants. We will use the same sample data:

```{r}
positivity_table <- pos_anyroi(fake_DVR, cutoffs)

positivity_table

```

The algorithm identified 11 PIB+ participants. In the generation of the sample data, the DVRs from the first 10 participants were drawn from a normal distribution with mean 1.9, sd 0.6 and for the latter 40 participants, from mean 1.3, sd 0.3; thus this pattern is in line with what we would expect: all 10 of the first participants are PIB+, and just 1 of the latter 40 was (by chance).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROI_definitions.R
\name{roi_ham_stand}
\alias{roi_ham_stand}
\title{Return a list of merged ROIs made up of the atomic ROIs in the Hammer's 
atlas.}
\usage{
roi_ham_stand()
}
\value{
A list of lists, where each list is an ROI (e.g.) frontal lobe that 
specifies the atomic ROIs from the atlas that make it up.
}
\description{
Return a list of merged ROIs made up of the atomic ROIs in the Hammer's 
atlas.
}
\examples{
roi_ham_stand()
}
\references{
Hammers, Alexander, Richard Allom, Matthias J. Koepp, Samantha L. Free, 
Ralph Myers, Louis Lemieux, Tejal N. Mitchell, David J. Brooks, and John S. 
Duncan. 2003. Three-dimensional Maximum Probability Atlas of the Human 
Brain, with Particular Reference to the Temporal Lobe. Human Brain Mapping 
19 (4): 224-247. doi:10.1002/hbm.10123
}
\seealso{
Other ROI definitions: \code{\link{roi_ham_full}},
  \code{\link{roi_ham_pib}}
}
\concept{ROI definitions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tac_methods.R
\name{as.tac}
\alias{as.tac}
\title{Creates a tac object from a data.frame}
\usage{
as.tac(x, time_unit = NULL, activity_unit = NULL)
}
\arguments{
\item{x}{data.frame with start, end time and tac data}

\item{time_unit}{NULL if in data.frame or set to "seconds" or "minutes"}

\item{activity_unit}{NULL if in data.frame or set to "kBq/cc", "Bq/cc", 
"nCi/cc"}
}
\value{
tac object
}
\description{
tac objects can be created from data.frame objects with `as.tac()`. The time 
and activity units must be specified as arguments if not already set as 
attributes in the data.frame. The columns of the data frame are the regional
time activity curves, with the column names the names of the ROIs.
}
\details{
If the time_unit and activity_unit attributes are already in the data.frame,
they do not need to be set again, but otherwise they will need to be
specified in the input parameters.
}
\examples{
manual <- data.frame(start=c(0:4), end=c(2:6), 
                     ROI1=c(10.1:14.2), ROI2=c(11:15))
manual_tac <- as.tac(manual, time_unit="minutes", activity_unit="kBq/cc")
}
\seealso{
Other Loading functions: \code{\link{load_tac}},
  \code{\link{load_voistat}}, \code{\link{load_vol}}
}
\concept{Loading functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batches.R
\name{batch_load}
\alias{batch_load}
\title{Load (+/- merge) ROIs for batch of participants}
\usage{
batch_load(participants, dir = "", tac_file_suffix = ".tac",
  tac_format = "PMOD", roi_m = FALSE, PVC = NULL,
  vol_file_suffix = NULL, vol_format = NULL, merge = NULL,
  ROI_def = NULL, tracer_dose = NULL, dose_unit = NULL,
  weight_kg = NULL)
}
\arguments{
\item{participants}{A vector of participant IDs}

\item{dir}{A directory and/or file name prefix for the tac/volume files}

\item{tac_file_suffix}{How participant IDs corresponds to the TAC files}

\item{tac_format}{Format of tac files provided: See load_tac()}

\item{roi_m}{TRUE if you want to merge atomic ROIs into larger ROIs (and if 
not, the following parameters are not used)}

\item{PVC}{For PVC, true where the data is stored as _C in same tac file}

\item{vol_file_suffix}{How participant IDs correspond to volume files}

\item{vol_format}{The file format that includes volumes: See load_vol()}

\item{merge}{Passes value to tac_roi(); T to also incl. original atomic ROIs}

\item{ROI_def}{Object that defines combined ROIs, see ROI_definitions.R}

\item{tracer_dose}{optionally, a vector of tracer doses (in the same order as
participants), for SUV}

\item{dose_unit}{if tracer_dose is specified, note the unit (e.g "MBq")}

\item{weight_kg}{optionally, a vector of participant weights in kg, for SUV}
}
\value{
A list of data.frames, each is a participant's TACs
}
\description{
For a vector of participant IDs and correspondingly named tac files,
this loads the tac files. If roi_m = T, then can also merge ROIs into 
larger ROIs based on the optional parameters that follow.
}
\details{
See load_tac() for specifics.
}
\examples{
# For the working example, the participants are full filenames.
participants <- c(system.file("extdata", "AD06.tac", package="tacmagic"),
                  system.file("extdata", "AD07.tac", package="tacmagic"),
                  system.file("extdata", "AD08.tac", package="tacmagic"))

tacs <- batch_load(participants, tac_file_suffix="")
}
\seealso{
Other Batch functions: \code{\link{batch_tm}},
  \code{\link{batch_voistat}}
}
\concept{Batch functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batches.R
\name{batch_voistat}
\alias{batch_voistat}
\title{Obtain values from voistat files (using load_voistat() for a batch.}
\usage{
batch_voistat(participants, ROI_def, dir = "", filesuffix = ".voistat",
  varname = "VALUE")
}
\arguments{
\item{participants}{A vector of participant IDs}

\item{ROI_def}{Object that defines combined ROIs, see ROI_definitions.R}

\item{dir}{Directory and/or filename prefix of the files}

\item{filesuffix}{Optional filename characters between ID and ".voistat"}

\item{varname}{The name of the variable being extracted, e.g. "SRTM"}
}
\value{
A table of values for the specified ROIs for all participants
}
\description{
For a vector of participant IDs and correspondingly named .voistat files,
this extracts the value from the files for the specified ROIs.
participants can also be a vector of filenames, in which case set dir="" and
filesuffix="", as in the example.
}
\details{
See load_voistat() for specifics.
}
\examples{
participants <- c(system.file("extdata", "AD06_BPnd_BPnd_Logan.voistat", 
                              package="tacmagic"),
                   system.file("extdata", "AD07_BPnd_BPnd_Logan.voistat", 
                               package="tacmagic"),
                   system.file("extdata", "AD08_BPnd_BPnd_Logan.voistat", 
                               package="tacmagic"))

batchtest <- batch_voistat(participants=participants, ROI_def=roi_ham_pib(), 
                           dir="", filesuffix="", varname="Logan") 

}
\seealso{
Other Batch functions: \code{\link{batch_load}},
  \code{\link{batch_tm}}
}
\concept{Batch functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tac.R
\name{split_pvc}
\alias{split_pvc}
\title{Subset PMOD tacs with or without PVC}
\usage{
split_pvc(tac, PVC = TRUE)
}
\arguments{
\item{tac}{The time-activity curve data from loading function (PMOD)}

\item{PVC}{If TRUE, includes columns with "_C", if FALSE, ones without "_C"}
}
\value{
Time-activity curve object
}
\description{
When partial volume correction (PVC) is used in PMOD, the saved tac files
have ROIs with and without PVC. When loaded with load_tac()) it may be 
desirable to keep only either the PVC or non-PVC tacs. This returns a tac 
object that is a subset of the input tac object with only the PVC or non-PVC
tacs. This relies on PMOD's convention of labelling tac columns with "_C".
}
\examples{
# f_raw_tac and f_raw_vol are the filenames of PMOD-generated files
f_raw_tac <- system.file("extdata", "AD06.tac", package="tacmagic") 

tac <- load_tac(f_raw_tac)
tac_pvc <- split_pvc(tac, TRUE)
tac_nc <- split_pvc(tac, FALSE)
}
\seealso{
Other tac functions: \code{\link{plot.tac}},
  \code{\link{save_tac}}, \code{\link{tac_roi}}
}
\concept{tac functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cutoff.R
\name{pos_anyroi}
\alias{pos_anyroi}
\title{Dichotomize participants based on ROI cutoff values}
\usage{
pos_anyroi(modelstats, cutoff)
}
\arguments{
\item{modelstats}{SUVR or DVR data for group of participants from batch_tm()}

\item{cutoff}{cutoffs for ROIs as from cutoff_aiz()}
}
\value{
data.frame of participants and positive/negative status
}
\description{
Aizenstein et al. (2008) proposed a standardized method of calculating PIB+ 
cutoff values to classify participants as PIB+ or PIB-. They used the DVR 
from 7 ROIs associated with amyloid deposition. This function takes the 
ROI-based cutoff values, e.g. from cutoff_aiz(), and returns a table 
specifying which participants are positive, i.e. which have at least one ROI
greater than the cutoff.
}
\references{
Aizenstein HJ, Nebes RD, Saxton JA, et al. 2008. Frequent amyloid 
deposition without significant cognitive impairment among the elderly. 
Arch Neurol 65: 1509-1517.
}
\seealso{
Other Cutoff functions: \code{\link{cutoff_aiz}}
}
\concept{Cutoff functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tacmagic.R
\docType{package}
\name{tacmagic}
\alias{tacmagic}
\alias{tacmagic-package}
\title{tacmagic: PET Analysis in R}
\description{
The main features of tacmagic are to load PET time activity curve (tac) data 
from multiple formats, merge ROIs weighted for volume, calculate binding 
potential models including SUVR and DVR, basic plotting, and calculation of 
cut-off values. Please see the walkthrough vignette for a detailed overview.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{change_units}
\alias{change_units}
\title{Convert radioactivity units}
\usage{
change_units(x, to_unit, from_unit)
}
\arguments{
\item{x}{time-activity curve or numeric object}

\item{to_unit}{the desired unit (e.g. "kBq")}

\item{from_unit}{not used for tac object (it is in the tac object), but for 
numeric objects, must be specified (e.g. "nCi")}
}
\value{
the converted object, same type as x
}
\description{
Change the radioactivity units of a tac or numeric object to the specified
desired units (e.g. Bq, kBq, MBq, nCi, uCi, mCi, Ci). For convenience, if the
unit is per volume ("x/cc" or "x/mL"), the "/cc" part is ignored for the 
conversion.
}
\examples{
f <- system.file("extdata", "AD06.tac", package="tacmagic")
AD06_tac <- load_tac(f, format="PMOD")
AD06_tac_nCicc <- change_units(AD06_tac, to_unit = "nCi/cc")

change_units(5, to_unit = "kBq", from_unit = "nCi")
change_units(0.185, to_unit = "nCi", from_unit = "kBq")
}
\concept{unit functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loading.R
\name{load_voistat}
\alias{load_voistat}
\title{Reads PMOD .voistat files and optionally merges volume-weighted ROIs}
\usage{
load_voistat(filename, ROI_def = NULL, model = "VALUE")
}
\arguments{
\item{filename}{(e.g. participant_logan.voistat)}

\item{ROI_def}{Optional ROI definitions to combine ROIs (e.g. roi_ham_pib())}

\item{model}{A string to name the variable being extracted, e.g. "Logan_DVR"}
}
\value{
data.frame with loaded model data in specified combined weighted ROIs
}
\description{
PMOD can produce .voistat files with the average model values by ROI for 
its voxelwise binding potential (BPnd) models, such as Logan, SRTM, etc.
This function reads the .voistat file and returns a data.frame with the
ROI as rows and the model value as the column. Optionally, the ROIs can be
combined into larger ROIs if ROI_def is specified, just as with TAC loading.
}
\examples{
f <- system.file("extdata", "AD06_BPnd_BPnd_Logan.voistat", 
                 package="tacmagic")
vs <- load_voistat(f, ROI_def=roi_ham_pib(), model="Logan")
}
\seealso{
Other Loading functions: \code{\link{as.tac}},
  \code{\link{load_tac}}, \code{\link{load_vol}}
}
\concept{Loading functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.ref_Logan}
\alias{plot.ref_Logan}
\title{Non-invasive reference Logan plot}
\usage{
\method{plot}{ref_Logan}(x, ...)
}
\arguments{
\item{x}{Reference Logan model data object from DVR_ref_Logan()}

\item{...}{Additional parameters than can be passed to plotting function}
}
\value{
No return
f <- system.file("extdata", "AD06.tac", package="tacmagic")
fv <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")
AD06_tac <- load_tac(f, format="PMOD")
AD06_volume <- load_vol(fv, format="voistat")
AD06 <- tac_roi(tac=AD06_tac, volumes=AD06_volume, ROI_def=roi_ham_pib(),  
                merge=FALSE, PVC=FALSE)  
AD06_DVR_fr <- DVR_ref_Logan(AD06, target="frontal", ref="cerebellum",
                             k2prime=0.2, t_star=0) 
plot(AD06_DVR_fr)
}
\description{
This plots the non-invasive Logan plot.
}
\seealso{
Other Logan plot functions: \code{\link{DVR_all_ref_Logan}},
  \code{\link{DVR_ref_Logan}}, \code{\link{dvr}}
}
\concept{Logan plot functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saving.R
\name{save_tac}
\alias{save_tac}
\title{Save a tac object as a .tac file}
\usage{
save_tac(tac, outfile)
}
\arguments{
\item{tac}{The time-activity curve data, e.g. from load_tac() or tac_roi()}

\item{outfile}{The output filename}
}
\value{
Does not return an object, only saves a file
}
\description{
Saves a tac object, created by load_tac(), tac_roi() or manually, and 
saves it as a PMOD-formatted tac file. Using the .tac extension in the 
file name is recommended.
}
\seealso{
Other tac functions: \code{\link{plot.tac}},
  \code{\link{split_pvc}}, \code{\link{tac_roi}}
}
\concept{tac functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logan.R
\name{DVR_all_ref_Logan}
\alias{DVR_all_ref_Logan}
\title{Non-invasive reference Logan method for all ROIs in tac data}
\usage{
DVR_all_ref_Logan(tac_data, ref, k2prime, t_star, error = 0.1,
  method = "trapz", ...)
}
\arguments{
\item{tac_data}{The time-activity curve data from tac_roi()}

\item{ref}{Required -- The reference region, e.g. "cerebellum"}

\item{k2prime}{Required -- A fixed value for k2' must be specified (e.g. 0.2)}

\item{t_star}{Required -- If 0, t* will be calculated using find_t_star()}

\item{error}{For find_t_star()}

\item{method}{Method of integration, "trapz" or "integrate"}

\item{...}{When called from tm_batch, unused parameters may be supplied}
}
\value{
Data frame with calculated DVRs for all ROIs
}
\description{
This calculates the DVR using the non-invasive reference Logan method for
all TACs in a supplied tac file. It uses DVR_ref_Logan.
}
\examples{
f <- system.file("extdata", "AD06.tac", package="tacmagic")
fv <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")
AD06_tac <- load_tac(f, format="PMOD")
AD06_volume <- load_vol(fv, format="voistat")
AD06 <- tac_roi(tac=AD06_tac, volumes=AD06_volume, ROI_def=roi_ham_pib(),  
                merge=FALSE, PVC=FALSE)  

AD06_DVR <- DVR_all_ref_Logan(AD06, ref="cerebellum", k2prime=0.2, t_star=23)

}
\references{
Logan, J., Fowler, J. S., Volkow, N. D., Wang, G.-J., 
Ding, Y.-S., & Alexoff, D. L. (1996). Distribution Volume Ratios without 
Blood Sampling from Graphical Analysis of PET Data. Journal of Cerebral 
Blood Flow & Metabolism, 16(5), 834-840. 
https://doi.org/10.1097/00004647-199609000-00008
}
\seealso{
Other Logan plot functions: \code{\link{DVR_ref_Logan}},
  \code{\link{dvr}}, \code{\link{plot.ref_Logan}}
}
\concept{Logan plot functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logan.R
\name{DVR_ref_Logan}
\alias{DVR_ref_Logan}
\title{Non-invasive reference Logan method}
\usage{
DVR_ref_Logan(tac_data, target, ref, k2prime, t_star, error = 0.1,
  method = "trapz")
}
\arguments{
\item{tac_data}{The time-activity curve data from tac_roi()}

\item{target}{The name of the target ROI, e.g. "frontal"}

\item{ref}{The reference region, e.g. "cerebellum"}

\item{k2prime}{A fixed value for k2' must be specified (e.g. 0.2)}

\item{t_star}{If 0, t* will be calculated using find_t_star()}

\item{error}{For find_t_star()}

\item{method}{Method of integration, "trapz" or "integrate"}
}
\value{
Data frame with calculate DVRs for all ROIs
}
\description{
This calculates the coefficient from the non-invasive Logan method, which
is equal to DVR. Works for a single tac (target).
}
\examples{
f <- system.file("extdata", "AD06.tac", package="tacmagic")
fv <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")
AD06_tac <- load_tac(f, format="PMOD")
AD06_volume <- load_vol(fv, format="voistat")
AD06 <- tac_roi(tac=AD06_tac, volumes=AD06_volume, ROI_def=roi_ham_pib(),  
                merge=FALSE, PVC=FALSE)                             
               
AD06_DVR_fr <- DVR_ref_Logan(AD06, target="frontal", ref="cerebellum",
                             k2prime=0.2, t_star=0) 
                            
}
\references{
Logan, J., Fowler, J. S., Volkow, N. D., Wang, G.-J., 
Ding, Y.-S., & Alexoff, D. L. (1996). Distribution Volume Ratios without 
Blood Sampling from Graphical Analysis of PET Data. Journal of Cerebral 
Blood Flow & Metabolism, 16(5), 834-840. 
https://doi.org/10.1097/00004647-199609000-00008
}
\seealso{
Other Logan plot functions: \code{\link{DVR_all_ref_Logan}},
  \code{\link{dvr}}, \code{\link{plot.ref_Logan}}
}
\concept{Logan plot functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROI_definitions.R
\name{roi_ham_full}
\alias{roi_ham_full}
\title{Return a list of larger ROIs made up of the ROIs in the Hammer's atlas.}
\usage{
roi_ham_full()
}
\value{
A list of lists, where each list is an ROI (e.g.) frontal lobe that 
specifies the atomic ROIs from the atlas that make it up.
}
\description{
This includes the cortical regions of roi_ham_stand() but also other regions.
It can be modified to suit the user's needs.
}
\examples{
roi_ham_full()
}
\references{
Hammers, Alexander, Richard Allom, Matthias J. Koepp, Samantha L. Free, 
Ralph Myers, Louis Lemieux, Tejal N. Mitchell, David J. Brooks, and John S. 
Duncan. 2003. Three-dimensional Maximum Probability Atlas of the Human 
Brain, with Particular Reference to the Temporal Lobe. Human Brain Mapping 
19 (4): 224-247. doi:10.1002/hbm.10123
}
\seealso{
Other ROI definitions: \code{\link{roi_ham_pib}},
  \code{\link{roi_ham_stand}}
}
\concept{ROI definitions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SUV.R
\name{tac_suv}
\alias{tac_suv}
\title{Calculate SUV from TAC}
\usage{
tac_suv(tac, dose = NULL, dose_unit = NULL, weight_kg = NULL)
}
\arguments{
\item{tac}{time-activity curve object (decay-corrected)}

\item{dose}{the injected tracer dose}

\item{dose_unit}{unit of tracer dose (e.g. "MBq", "kBq", "mCi"...)}

\item{weight_kg}{the participant's weight in kg}
}
\value{
tac object with SUV values
}
\description{
Calculate the standardized uptake value (SUV) time-activity curve from a tac
object, the participant's weight, and the tracer dose. The weight must be in
kg, and the tracer dose must be specified. The dose is converted to MBq, the
tac is converted to kBq/cc, and the final SUV units are thus in g/cc. Aside
from the tac object, the remaining parameters should be left NULL if the
required data is in the tac object attributes (as can be done with
batch_load().
}
\examples{
f <- system.file("extdata", "AD06.tac", package="tacmagic")
fv <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")
AD06_tac <- load_tac(f, format="PMOD")
AD06_volume <- load_vol(fv, format="voistat")
AD06 <- tac_roi(tac=AD06_tac, volumes=AD06_volume, ROI_def=roi_ham_pib(),
                merge=FALSE, PVC=FALSE)
# dose and weight are fabricated for the example
AD06_suv <- tac_suv(AD06, dose = 9.0, dose_unit = "mCi", weight_kg = 70)
}
\seealso{
Other SUV functions: \code{\link{suv}}
}
\concept{SUV functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batches.R
\name{batch_tm}
\alias{batch_tm}
\title{Calculate one or more models for a batch of participants}
\usage{
batch_tm(all_tacs, models, custom_model = NULL, ...)
}
\arguments{
\item{all_tacs}{A list by participant, of tac data (load_batch())}

\item{models}{A vector of names of the models to calculate}

\item{custom_model}{A function that can be run like other models (advanced)}

\item{...}{The arguments that get passed to the specified models/custom model,
many are required; please check with model desired.}
}
\value{
A table of SUVR values for the specified ROIs for all participants
}
\description{
For a list of tac data (from load_batch) this calculates specified models
and saves in a tidy data.frame. Current model options are "SUVR", "Logan".
}
\details{
For further details about how the models are calculated, see the individual
functions that they rely on. "SUVR" uses suvr(), "Logan" uses
DVR_all_ref_Logan().
}
\examples{
participants <- c(system.file("extdata", "AD06.tac", package="tacmagic"),
                  system.file("extdata", "AD07.tac", package="tacmagic"),
                  system.file("extdata", "AD08.tac", package="tacmagic"))

tacs <- batch_load(participants, tac_file_suffix="")

# Keeps only the ROIs without partial-volume correction (PMOD convention)
tacs <- lapply(tacs, split_pvc, FALSE)

batch <- batch_tm(tacs, models=c("SUVR", "Logan"), ref="Cerebellum_r",
                  SUVR_def=c(3000,3300,3600), k2prime=0.2, t_star=23)

}
\seealso{
Other Batch functions: \code{\link{batch_load}},
  \code{\link{batch_voistat}}
}
\concept{Batch functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SUVR.R
\name{suvr}
\alias{suvr}
\title{Calculate weighted SUVRs for specified regions of interest}
\usage{
suvr(tac, SUVR_def, ref, ...)
}
\arguments{
\item{tac}{The time-activity curve data from tac_roi()}

\item{SUVR_def}{a vector of start times for window to be used in SUVR}

\item{ref}{a string, e.g. "cerebellum", to specify reference region}

\item{...}{When called from tm_batch, unused parameters may be supplied}
}
\value{
A data.frame of SUVR values for the specified ROIs
}
\description{
Calculate the standardized uptake value ratio (SUVR) for all ROIs in the
provided tac data, using the specified reference region.
}
\examples{
f <- system.file("extdata", "AD06.tac", package="tacmagic")
fv <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")
AD06_tac <- load_tac(f, format="PMOD")
AD06_volume <- load_vol(fv, format="voistat")
AD06 <- tac_roi(tac=AD06_tac, volumes=AD06_volume, ROI_def=roi_ham_pib(),  
                merge=FALSE, PVC=FALSE)

AD06_SUVR <- suvr(AD06, SUVR_def=c(3000,3300,3600), ref="cerebellum")

}
\seealso{
Other SUVR functions: \code{\link{suvr_auc}}
}
\concept{SUVR functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loading.R
\name{load_tac}
\alias{load_tac}
\title{Loads TAC from file for use by other functions (default is PMOD .tac format)}
\usage{
load_tac(filename, format = "PMOD", acqtimes = NULL,
  time_unit = NULL, activity_unit = NULL)
}
\arguments{
\item{filename}{(e.g. "participant01.tac")}

\item{format}{A character string, with options listed above (e.g. "PMOD")}

\item{acqtimes}{Filename for a .acqtimes file (as in PMOD), required for 
format="voistat"}

\item{time_unit}{NULL if in file (e.g. PMOD .tac), or set to "seconds" or 
"minutes" if not in file or to override file}

\item{activity_unit}{NULL if in file (e.g. PMOD .tac), or set to "kBq/cc", 
"Bq/cc", "nCi/cc"}
}
\value{
tac object
}
\description{
This is the main function for loading an individual participant's TAC data.
The minimal required information within the supplied files is the start and 
stop times and a time unit (either seconds or minutes), as well as the 
activity values for 1 or more ROIs, and units for activity. The currently 
supported formats (with the corresponding format argument), include:
\itemize{
  \item "PMOD": PMOD .tac files
  \item "voistat": PMOD TAC .voistat files used in combination with PMOD 
         .acqtimes file for start/stop times.
  \item "magia": magia pipeline .mat tac file
  \item "DFT": Turku PET Centre's DFT format
}
}
\examples{
f_raw_tac <- system.file("extdata", "AD06.tac", package="tacmagic") 
tac <- load_tac(f_raw_tac)

}
\seealso{
Other Loading functions: \code{\link{as.tac}},
  \code{\link{load_voistat}}, \code{\link{load_vol}}
}
\concept{Loading functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tac.R
\name{tac_roi}
\alias{tac_roi}
\title{Calculate weighted time-activity curves for specified regions of interest}
\usage{
tac_roi(tac, volumes, ROI_def, merge, PVC)
}
\arguments{
\item{tac}{The time-activity curve data from loading function}

\item{volumes}{The ROI volume data from loading function}

\item{ROI_def}{The definition of ROIs by combining smaller ROIs from TAC file}

\item{merge}{If TRUE, includes the original ROIs in the output data}

\item{PVC}{If TRUE, appends "_C" to ROI name header (as in PMOD TAC files)}
}
\value{
Time-activity curves for the specified ROIs
}
\description{
Calculate weighted time-activity curves for specified regions of interest
}
\examples{
# f_raw_tac and f_raw_vol are the filenames of PMOD-generated files
f_raw_tac <- system.file("extdata", "AD06.tac", package="tacmagic") 
f_raw_vol <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")

tac <- load_tac(f_raw_tac)
vol <- load_vol(f_raw_vol)
AD06_tac_nc <- tac_roi(tac, vol, roi_ham_full(), merge=FALSE, PVC=FALSE)
}
\seealso{
Other tac functions: \code{\link{plot.tac}},
  \code{\link{save_tac}}, \code{\link{split_pvc}}
}
\concept{tac functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fake_DVR}
\alias{fake_DVR}
\title{Fake DVR data for vignette and package testing}
\format{A data frame with 50 rows and 4 variables representing ROIs}
\usage{
fake_DVR
}
\description{
A fake dataset of 50 simulated participants in the format that the function
tm_batch() would be expected to produce with the "Logan" model specified.
The data itself was generated as follows: \cr \cr
#higher <- matrix(rnorm(40, 1.9, 0.6), ncol=4, nrow=10) \cr
#lower <- matrix(rnorm(160, 1.3, 0.3), ncol=4, nrow=40) \cr
#fake_data <- as.data.frame(rbind(higher, lower)) \cr
#row.names(fake_data) <- paste0("p", 1:50) \cr 
#colnames(fake_data) <- c("ROI1_DVR", "ROI2_DVR", "ROI3_DVR", "ROI4_DVR")\cr
#save(fake_data, "fake_DVR.Rda") \cr
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SUVR.R
\name{suvr_auc}
\alias{suvr_auc}
\title{Calculate SUVRs for regions of interest with AUC from mid-frame times}
\usage{
suvr_auc(tac, SUVR_def, ref, ...)
}
\arguments{
\item{tac}{The time-activity curve data from tac_roi()}

\item{SUVR_def}{a vector of start times for window to be used in SUVR}

\item{ref}{is a string, e.g. "cerebellum", to specify reference region}

\item{...}{When called from tm_batch, unused parameters may be supplied}
}
\value{
A data.frame of SUVR values for the specified ROIs
#' f <- system.file("extdata", "AD06.tac", package="tacmagic")
fv <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")
AD06_tac <- load_tac(f, format="PMOD")
AD06_volume <- load_vol(fv, format="voistat")
AD06 <- tac_roi(tac=AD06_tac, volumes=AD06_volume, ROI_def=roi_ham_pib(),  
                merge=FALSE, PVC=FALSE)

AD06_SUVR <- suvr_auc(AD06, SUVR_def=c(3000,3300,3600), ref="cerebellum")
}
\description{
Calculate the standardized uptake value ratio (SUVR) for all ROIs in the
provided tac data, using the specified reference region. This is an 
alternate to suvr() which should provide very similar values.
}
\seealso{
Other SUVR functions: \code{\link{suvr}}
}
\concept{SUVR functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cutoff.R
\name{cutoff_aiz}
\alias{cutoff_aiz}
\title{Cutoff value calculation using method described in Aizenstein et al. 2008}
\usage{
cutoff_aiz(modelstats, ROIs)
}
\arguments{
\item{modelstats}{SUVR or DVR data for group of participants from batch_tm()}

\item{ROIs}{list of variables (ROIs) to use for cutoff detection}
}
\value{
Cutoff values for each ROI based on the above method
}
\description{
See the reference below and the tacmagic walkthrough vignette. Aizenstein et
al. (2008) proposed a standardized method of calculating Pittsburgh Compound
B (PIB) cutoff values to classify participants as PIB+ or PIB-. They used the
distribution volume ratio (DVR) from several ROIs associated with amyloid 
deposition. The steps are summarized below. cutoff_aiz() implements 1-3,
returning cutoff values for each ROI. It can be used to dichotomize
participants, with pos_anyroi().
}
\details{
1. Remove outliers from a group of cognitively normal individuals. An outlier
is defined as having any ROI with DVR > upper inner fence of that ROI (= 3rd
quartile + (1.5 * IQR).
2. Iterate step 1 as needed until there are no more outlying participants.
3. From this subset of the group with outliers removed, the cutoff value for 
each ROI is set as the upper inner fence. 
4. For all participants, if there is any ROI above the cutoff for that 
region, then the participant is deemed to be PIB+.
}
\examples{
cutoff_aiz(fake_DVR, c("ROI1_DVR", "ROI2_DVR", "ROI3_DVR", "ROI4_DVR"))
}
\references{
Aizenstein HJ, Nebes RD, Saxton JA, et al. 2008. Frequent amyloid 
deposition without significant cognitive impairment among the elderly. 
Arch Neurol 65: 1509-1517.
}
\seealso{
Other Cutoff functions: \code{\link{pos_anyroi}}
}
\concept{Cutoff functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loading.R
\name{load_vol}
\alias{load_vol}
\title{Loads ROI volumes from file for use by other functions}
\usage{
load_vol(filename, format = "voistat")
}
\arguments{
\item{filename}{(e.g. participant.voistat)}

\item{format}{(default is the TAC .voistat format from PMOD, also accepts 
"DFT and "BPndPaste")}
}
\value{
data.frame with loaded TAC data
}
\description{
Loads ROI volumes from file for use by other functions
}
\examples{
f_raw_vol <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")

vol <- load_vol(f_raw_vol)
}
\seealso{
Other Loading functions: \code{\link{as.tac}},
  \code{\link{load_tac}}, \code{\link{load_voistat}}
}
\concept{Loading functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.tac}
\alias{plot.tac}
\title{Plots time activity curves from 1 or 2 participants or groups.}
\usage{
\method{plot}{tac}(x, tac2 = NULL, ROIs, ymax = 25, time = "minutes",
  title = "", colors = rainbow, ...)
}
\arguments{
\item{x}{A tac object containing time-activity curves to plot, e.g. from 
tac_roi() or load_tac()}

\item{tac2}{An optional, second TAC, to plot for comparison}

\item{ROIs}{A vector of ROIs to plot, names matching the TAC headers}

\item{ymax}{The maximum value on the y-axis}

\item{time}{"seconds" or "minutes" depending on desired x-axis, converts tac}

\item{title}{A title for the plot}

\item{colors}{If null, rainbow palette is used, otherwise another palette can
be specified (heat.colors, terrain.colors, topo.colors, cm.colors}

\item{...}{Additional arguments}
}
\value{
Creates a plot
}
\description{
Plots time activity curves from 1 or 2 participants or groups.
}
\examples{
# f_raw_tac and f_raw_vol are the filenames of PMOD-generated files
f_raw_tac <- system.file("extdata", "AD06.tac", package="tacmagic") 
f_raw_vol <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")

tac <- load_tac(f_raw_tac)
vol <- load_vol(f_raw_vol)
AD06_tac_nc <- tac_roi(tac, vol, roi_ham_full(), merge=FALSE, PVC=FALSE)
plot(AD06_tac_nc, ROIs=c("frontal", "cerebellum"), title="Example Plot")
}
\seealso{
Other tac functions: \code{\link{save_tac}},
  \code{\link{split_pvc}}, \code{\link{tac_roi}}
}
\concept{tac functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROI_definitions.R
\name{roi_ham_pib}
\alias{roi_ham_pib}
\title{Return a list of merged ROIs made up of atomic ROIs in the Hammer's atlas.}
\usage{
roi_ham_pib()
}
\value{
A list of lists, where each list is an ROI (e.g.) frontal lobe that 
specifies the atomic ROIs from the atlas that make it up.
}
\description{
This includes the ROIs from roi_ham_full and also the PIB cortical composite
ROI as defined in the PMOD documentation and as widely used in PIB studies.
See PMOD Neuro Tool (PNEURO) (Version 4.0) documentation.
}
\examples{
roi_ham_pib()
}
\references{
Hammers, Alexander, Richard Allom, Matthias J. Koepp, Samantha L. Free, 
Ralph Myers, Louis Lemieux, Tejal N. Mitchell, David J. Brooks, and John S. 
Duncan. 2003. Three-dimensional Maximum Probability Atlas of the Human 
Brain, with Particular Reference to the Temporal Lobe. Human Brain Mapping 
19 (4): 224-247. doi:10.1002/hbm.10123
}
\seealso{
Other ROI definitions: \code{\link{roi_ham_full}},
  \code{\link{roi_ham_stand}}
}
\concept{ROI definitions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dvr.R
\name{dvr}
\alias{dvr}
\title{Distribution volume ratio (DVR) for one or more ROIs}
\usage{
dvr(tac, model = "logan", target = NULL, ref, k2prime, t_star,
  error = 0.1, method = "trapz")
}
\arguments{
\item{tac}{The time-activity curve data from load_tac() or tac_roi()}

\item{model}{Only model currently available is "logan"}

\item{target}{Optional - otherwise will calculate DVR for all regions}

\item{ref}{Required -- The reference region, e.g. "cerebellum"}

\item{k2prime}{Required -- A fixed value for k2' must be specified (e.g. 0.2)}

\item{t_star}{Required -- If 0, t* will be calculated using find_t_star()}

\item{error}{For find_t_star()}

\item{method}{Method of integration, "trapz" or "integrate"}
}
\value{
Data frame with calculated DVRs
}
\description{
This calculates the DVR using the non-invasive reference Logan method for
all TACs in a supplied tac file. It uses DVR_ref_Logan if a target ROI is 
specified, otherwise will calculate DVR for all ROIs with DVR_ref_all_Logan()
}
\details{
For other model parameters, directly call DVR_ref_Logan().
}
\examples{
f <- system.file("extdata", "AD06.tac", package="tacmagic")
fv <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")
AD06_tac <- load_tac(f, format="PMOD")
AD06_volume <- load_vol(fv, format="voistat")
AD06 <- tac_roi(tac=AD06_tac, volumes=AD06_volume, ROI_def=roi_ham_pib(),  
                merge=FALSE, PVC=FALSE)  

AD06_DVRs <- dvr(AD06, ref="cerebellum", k2prime=0.2, t_star=23)

AD06_DVR <- dvr(AD06, target="frontal", ref="cerebellum", 
             k2prime=0.2, t_star=23)
}
\references{
Logan, J., Fowler, J. S., Volkow, N. D., Wang, G.-J., 
Ding, Y.-S., & Alexoff, D. L. (1996). Distribution Volume Ratios without 
Blood Sampling from Graphical Analysis of PET Data. Journal of Cerebral 
Blood Flow & Metabolism, 16(5), 834-840. 
https://doi.org/10.1097/00004647-199609000-00008
}
\seealso{
Other Logan plot functions: \code{\link{DVR_all_ref_Logan}},
  \code{\link{DVR_ref_Logan}}, \code{\link{plot.ref_Logan}}
}
\concept{Logan plot functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SUV.R
\name{suv}
\alias{suv}
\title{Calculate average SUV over time window, or maximum SUV}
\usage{
suv(tac, SUV_def, dose = NULL, dose_unit = NULL, weight_kg = NULL,
  ...)
}
\arguments{
\item{tac}{time-activity curve object (decay-corrected)}

\item{SUV_def}{vector of start times for window for SUV weighted average, or
alternatively, "max" for the maximum ROI SUV value}

\item{dose}{the injected tracer dose}

\item{dose_unit}{unit of tracer dose (e.g. "MBq", "kBq", "mCi"...)}

\item{weight_kg}{the participant's weight in kg}

\item{...}{When called from tm_batch, unused parameters may be supplied}
}
\value{
table of SUV values
}
\description{
Calculate the standardized uptake value (SUV) from a tac object, the
participant's weight, and the tracer dose. These values may be in the tac
object or manually supplied. The weight must be in kg, and the tracer units
must be specified. The dose is converted to MBq, the tac is converted to
kBq/cc, and the final SUV units are thus in g/cc. Aside from the tac object,
the remaining parameters should be left NULL if the required data is in the
tac object attributes (as can be done with batch_load()).
}
\examples{
f <- system.file("extdata", "AD06.tac", package="tacmagic")
fv <- system.file("extdata", "AD06_TAC.voistat", package="tacmagic")
AD06_tac <- load_tac(f, format="PMOD")
AD06_volume <- load_vol(fv, format="voistat")
AD06 <- tac_roi(tac=AD06_tac, volumes=AD06_volume, ROI_def=roi_ham_pib(),
                merge=FALSE, PVC=FALSE)
# dose and weight are fabricated for the example
AD06_suvmax <- suv(AD06, "max", dose = 9.0, dose_unit = "mCi",
                     weight_kg = 70)
AD06_suv <- suv(AD06, c(3000, 3300, 3600), dose = 9.0, dose_unit = "mCi",
                  weight_kg = 70)
}
\seealso{
Other SUV functions: \code{\link{tac_suv}}
}
\concept{SUV functions}
