# RHRT
[RHRT](/RHRT/README.md) is an R package to assess Heart Rate Turbulence from RR interval data.

This repository includes the source files of RHRT and its publication.
It is sorted into:

* the files of the [RHRT package](/RHRT)
* scripts used to create dummy RR interval data in [testdata-scripts](/testdata-scripts)
* all files needed for the [publication](/Paper) in JOSS

All contents of the repository are licensed under the [GNU GPL v2](http://www.gnu.org/licenses/old-licenses/gpl-2.0.html).

Please feel free to submit pull requests and use the issue tracker for issues or problems with the package.
For any other support or questions please contact [Valeria Blesius](mailto:rhrt@blesius.eu).# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement via mail to
[Valeria Blesius](mailto:rhrt@blesius.eu).
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.1, available at
[https://www.contributor-covenant.org/version/2/1/code_of_conduct.html][v2.1].

Community Impact Guidelines were inspired by
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.1]: https://www.contributor-covenant.org/version/2/1/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations

# Changelog

All notable changes of RHRT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.1] 2021-09-17
### Added
- added parameters to the plot methods to hide and adjust legend
### Changed
- updated vignette and broke it down into smaller parts
### Fixed
- corrected the example pipelines in the vignette

## [1.0.1] 2021-06-28

### Changed
- the reliability check returns now 0 as p-values (formerly NA) if parameter values are all identical (and there is more than one HRT in the HRTList)
### Fixed
- corrected testdata which was prone for round-off errors---
title: 'RHRT: An R package to assess Heart Rate Turbulence'
tags:
  - Heart Rate Turbulence
  - R package
  - data analysis
  - Cardiology
  - arrhythmia
  - ventricular premature beats
authors:
  - name: Valeria Blesius^[corresponding author]
    orcid: 0000-0002-2391-242X
    affiliation: "1"
  - name: Andreas Dominik^[co-author]
    orcid: 0000-0002-9368-0812
    affiliation: "1"
affiliations:
 - name: THM University of Applied Sciences, Giessen, Germany
   index: 1
date: 17 June 2021
bibliography: paper.bib

---

[](250 - 1000 words)
[](Mention if applicable a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it)

# Summary

[](describing the high-level functionality and purpose of the software for a diverse, non-specialist audience)

Heart Rate Turbulence (HRT) is a naturally occurring phenomenon of the heart that was first described by @schmidt_heart-rate_1999.
After a premature ventricular contraction (VPC) – an atypical heart beat that originates in the ventricles and occurs with a lower latency than a regular beat – the heart rate fluctuates.
This fluctuation consists of a fast increase and later decrease in heart rate.
This variability of interval lengths is an indicator for the condition of the autonomic nervous system:
While no reaction to the premature ventricular beat and therefore no variation suggests an underlying pathology, a distinctive turbulence is considered healthy.
Therefore, HRT can be used in medicine to determine the risk status of a person especially with certain diseases or conditions, e.g. after a myocardial infarction. [@bauer_heart_2008]

To assess HRT in R we wrote the package ``RHRT`` that can be found on [CRAN](https://cran.r-project.org/package=RHRT) and [GitHub](https://github.com/VBlesius/RHRT).
It finds occurrences of HRT in heart beat interval data, calculates the most used parameters (Turbulence Onset (TO), Turbulence Slope (TS), Turbulence Timing (TT) and normalised Turbulence Slope (nTS)) and plots the results.
The package works best and fastest when given annotation data, but can also find HRT based on commonly accepted filtering rules that were first published in @grimm_heart_2003.
Most filtering parameters and calculation methods can be freely adjusted to enable research on the methodology itself.
In addition to parameter calculation, ``RHRT`` can classify the data into common risk categories (HRT0-2 and HRTA-C) and estimate the reliability of the results based on the number and parameter values of the HRTs.

# Statement of need

[](clearly illustrates the research purpose of the software)

Since it reflects the status of the autonomic nervous system, HRT is a feasible method to estimate the health risk of a person [@lombardi_origin_2011].
Together with other autonomic markers as heart rate variability, HRT analysis is already used for risk stratification in the clinical practice.
Several tools for HRT analysis have been published until now, but are not available anymore (an HRT program available on request on the discontinued [www.h-r-t.org](www.h-r-t.org) and a software tool published in @kudrynski_computer_2011) or are restricted to specific platforms like HRVAnalysis to Windows [@pichot_hrvanalysis_2016].
Furthermore, they focus on risk stratification and therefore implement the standard HRT assessment workflows and parameters.

However, analysis of the methodology of HRT assessment is still needed [@blesius_hrt_2020], because optimal filtering and calculation parameters have not been systematically assessed yet.
Variations in methodology can lead to less comparable or even conflicting data and reduce the validity of HRT.
To tackle this, a tool like ``RHRT`` is required with batch processing and the possibility to alter the used methodology.
For example, the package enabled us to study the optimal number of intervals needed for TS calculation, which is a continous issue that leads to reduced comparability of research data [@blesius_comparability_2021].

# Minimal example

To install ``RHRT`` use

`install.packages("RHRT")`

for the version on CRAN.
To install the continuosly developed version on GitHub you can use the devtools package:

`devtools::install_github("VBlesius/RHRT/RHRT")`



The most straightforward use of `RHRT` is to scan your interval data for valid HRTs and analyse the results:

``` r
library(RHRT)
## scan your interval data and save the results as an HRTList
### the input should be a numeric vector consisting of RR interval data
### testdataLong is dummy data included in the package
hrtl <- vectorToHRT(testdataLong)

## get the HRT class of your data
getResults(hrtl, type = "class")
#> [1] "HRT0"

## have a look at the data and the parameters
plot(hrtl)
```
![Plot of an HRTList Object: the plot resembles the standard visualisation of HRT – a tachogram – in which the indices of the intervals are plotted against their lengths. The tachogram and HRT parameters TO and TS are drawn in black, red and blue, respectively, while the tachograms of all underlying HRTs are drawn in grey in the background. The plot is zoomed in by default to show the HRT parameter values more precisely, so the intervals before and after the VPC are outside the plot range. RR interval: interval between two heartbeats measured between the R-peaks in the ECG, couplRR: coupling interval (interval between last sinus induced contraction and VPC), compRR: compensatory interval (interval between VPC and following sinus induced contraction), TO: Turbulence Onset, TS: Turbulence Slope. \label{fig:plot}](../RHRT/man/figures/README-example-1.png)

More examples and a detailed description of the objects and functions can be found in the [vignette](https://github.com/VBlesius/RHRT/blob/main/RHRT/vignettes/rhrt-vignette.md) of the package.

# Acknowledgements
We are sincerely grateful to [Christopher Schölzel](https://orcid.org/0000-0001-8627-0594) for testing the package and providing excellent input through the course of its development.

# References
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RHRT

RHRT provides tools for a Heart Rate Turbulence analysis of RR interval
data. It can either find the ventricular premature complexes (VPCs) via
a set of filter rules or can use annotation data to only check the beats
with the correct annotation. VPC snippets are filtered for validity and
HRT parameters calculated.

In addition to standard calculation methods the package allows to modify
the filter and calculation parameters. It is therefore not only helpful
to identify HRT classes of measurements for risk assessment but also for
assessment of the methodology itself.

For more information please check the vignette of the package. For more
information about HRT have a look into the [original publication by
Schmidt et al.](https://doi.org/10.1016/S0140-6736(98)08428-1) or our
[review](https://doi.org/10.1088/1361-6579/ab98b3) with focus on the
methodology.

## Installation

You can install the released version of RHRT from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RHRT")
```

And the development version from [GitHub](https://github.com/) with:

``` r
install.packages("devtools") # if not already done
devtools::install_github("VBlesius/RHRT/RHRT")
```

## Example

The general workflow of RHRT is to scan your interval data for HRT and
check the results via HRT class and plot:

``` r
library(RHRT)
## scan your interval data and save the results as an HRTList
hrtl <- vectorToHRT(testdataLong)
## get the HRT class of your data
getResults(hrtl, type = "class")
#> [1] "HRT0"
## have a look at the data and the parameters
plot(hrtl)
```

<img src="man/figures/README-example-1.png" width="100%" />

## Data

Data to test the package can be found on
[Physionet](https://physionet.org/). Via the [WFDB
Toolkit](https://physionet.org/content/wfdb/10.6.2/) ECG data can be
downloaded and/or converted, for example:

``` bash
ann2rr -r chf2db/chf201 -a ecg -i s3 -w > ~/some/path/chf201.csv
```

Then load the data and use RHRT to find VPCSs:

``` r
chf201 <- read.table("~/some/path/chf201.csv")
hrtl <- RHRT::vectorToHRT(chf201[[1]]*1000, ann = chf201[[2]])
```

More example workflows can be found in the vignettes: [Quick start
guide](vignettes/synopsis.md), [Objects &
Functions](vignettes/objects_functions.md), [Scientific
Background](vignettes/background.md) and [Example
Pipelines](vignettes/examples.md).

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->
---
title: "RHRT: Synopsis for the hasty (Quick-start guide)"
author: "Valeria Blesius"
date: "2021-08-12"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{RHRT: Quick-start guide}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The RHRT package helps you assess **Heart Rate Turbulence** (HRT) in RR intervals and calculate **turbulence onset** (TO), **slope** (TS) and **timing** (TT). It can plot the tachograms and checks the results for reliability. The **ventricular premature beats** (VPCs) with **coupling** (CPI) and **compensatory interval** (CMI) can either be given with annotations or found on the basis of the filter rules as first described  by [Grimm et al. 2003](https://doi.org/10.1046/j.1542-474X.2003.08206.x). The type of average and order of calculation for all parameters can be set.

This vignette sums up the most common functions and parameters needed when using RHRT.

--------

## Loading package and data


```r
library("RHRT")
# testdataLong is a numeric vector of RR intervals in msec
data("testdataLong", package = "RHRT")
ints <- testdataLong
# testdataLong_Ann is a character vector of annotations corresponding to testdataLong
data("testdataLong_Ann", package = "RHRT")
ann <- testdataLong_Ann
```

## Checking interval data for HRTs

The **core function** of RHRT is `vectorToHRT` that finds valid VPCs in RR intervals and returns an `HRTList` object (see *HRTList object* in [Objects & Functions](objects_functions.md) for more information):


```r
hrtl <- vectorToHRT(ints) 
```

Every RR interval sequence that matches the needed interval lengths is considered to be a coupling and compensatory interval of a VPC, which can lead to wrong matches. If your data is annotated, you can provide the **annotation data** with the parameters `annotations` and `PVCAnn`.


```r
hrtl <- vectorToHRT(ints, annotations = ann, PVCAnn = "V")
```

Other parameters are:

* `numPreRRs` & `numPostRRs` are used to modify the **filter rules** to find HRTs (number of intervals before and after the VPC that have to match the filter criteria).
* `minHRT` is the **minimal number of HRTs** needed to calculate HRT / create a HRTList
* `normHallstrom` defines whether TS should be **normalised** with the method of Hallstrom et al. (see the chapter *Normalisation of Turbulence Slope* in the [scientific background](background.md) for more information). 

## Getting HRT parameters or class

```r
getResults(hrtl) # get the HRT class of the data
```

```
## [1] "HRT0"
```

Per default `getResults` checks whether all needed HRT parameters can be calculated reliably. This is done via a t-test per parameter value (for more information see chapter *Reliability Check* in the [scientific background](background.md) vignette). If any of the parameter values is **not reliable** `getResults` returns NR (not reliable). 


```r
getResults(hrtl, safe = FALSE) # get the HRT class without safety check
```

```
## [1] "HRT0"
```

In addition to the classification system HRT0-2 RHRT implements **HRTA-C** that is based on the three parameters TO, TS and TT. 


```r
getResults(hrtl, safe = FALSE, TT = TRUE) # include TT
```

```
## [1] "HRTA"
```

With the parameter `type` you can choose between getting only the HRT **class**, all **parameter values** or the parameter values with the corresponding **p-values** (types "class", "parameter" or "full", respectively).


```r
getResults(hrtl, type = "parameter", TT = TRUE) # get the averaged HRT parameters
```

```
##        TO        TS        TT 
## -10.57551  37.61417   3.00000
```

Other parameters are:

* `nTS`: the **normalised TS** is returned or used for classification instead of TS.
* `num`: forces the function to return **numerics** when using `type = parameter`. Depending on the results and your setting of `type` the `getResults` returns characters or numerics.
* `pmax`: changes the needed **significance level** for the parameters to be reliable.

## Plotting


```r
plot(hrtl, TT = TRUE) # plots the averaged VPCS and all underlying VPCSs in background
```

![](synopsis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

Per default the VPCS is **zoomed** in. If you want to also see the CPI and CMI use `cropped = FALSE`.


```r
plot(hrtl, cropped = FALSE) # shows also coupling and compensatory interval
```

![](synopsis_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

--------

Further information can be found in the other vignettes about the [objects and functions](objects_functions.md), the [scientific background](background.md) or with [example pipelines](examples.md).

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->
---
title: "RHRT: Objects and Functions"
author: "Valeria Blesius"
date: "2021-09-17"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{RHRT: Objects and Functions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This vignette gives a detailed overview of all objects and functions provided by RHRT.

--------

## HRT Object

### Slots

An HRT object saves the data of one VPCS and its HRT results. It consists of the following slots:

- **Intervals**
  - `preRRs`: preceding regular intervals, 5 per default
  - `couplRR`: CPI
  - `compRR`: CMI
  - `postRRs`: following regular intervals, 15 per default

- **HRT Parameters**
  - `TO`
  - `TS`
  - `TT`
  - `nTS`: normalised TS, for more Information see chapter *Normalisation of Turbulence Slope* in [Scientific Background](background.md)

- **Line coefficients**
  - `intercept`: The intercept of the TS regression line that is needed for plotting
  - `nintercept`: Analogously, the intercept of `nTS`

### Functions

- **`getRRs`** &nbsp;&nbsp; This function returns all intervals saved in the HRT object. This if helpful for i.e. calculating an averaged VPCS out of a set of HRT objects.

- **`plot`** &nbsp;&nbsp; Per default `plot` displays a zoomed in plot of the VPCS with highlighted TO and TS and a legend given the rounded HRT parameter values. In addition to the common parameters of graphics.plot the method accepts the following parameters:
  - `cropped`: switches between showing a zoomed in version and the full VPCS including CPI and CMI, the default is `TRUE`
  - `add`: adds the plot to the current one
  - `TT`: highlights TT and includes it to the legend
  - `paramsLegend`: switches between showing and hiding the parameter values in the legend, the default is `TRUE`
  - `colTO`, `colTS` and `colTT` determine the colour in which the different parameters are highlighted in the plot (red, blue and purple per default).

## HRTList Object

An HRTList object is created when `vectorToHRT` is called. All HRT parameters are calculated automatically in this process and can be called with the `getHRTParams`-methods.

### Slots

The HRTList object sums up all HRTs found in the given dataset. The slots are:

- `name`: The name of the object if given to `vectorToHRT`
- `IL`: the average length of all intervals in the given cleaned vector (see *vectorToHRT: Cleaning Input* for more information) 
- `pos`: the indices of the CPIs as found in the given vector
- `HRTs`: list of all HRT objects found
- `avHRT`: an avHRT object averaged from all HRTs
- `RMSSD`: the HRV parameter RMSSD calculated from all intervals in the given cleaned vector (see *vectorToHRT: Cleaning Input* for more information): RMSSD is the square root of the mean of the squared respective differences of the successive RR intervals, in R calculated as `sqrt(mean(diff(intervals)^2))`

### Functions

- **`calcAvHRT`** &nbsp;&nbsp; The function calculates the parameters of the averaged HRT. This is called automatically when using `vectorToHRT` with default parameters. If the avHRT should be calculated differently the following options are available:
  - `av`: The function with which averaging should be done: `mean` or `median`.
  - `orTO` and `orTS`: sets the order in which TO and TS should be calculated. With `avAfter` the parameter is assessed separately for every HRT and averaged afterwards, with `avBefore` the intervals of all VPCSs are averaged first and the parameter is assessed afterwards. The default is `avAfter` for `orTO` and `avBefore` for `orTS`.
  - `IL`: the average interval length that is needed to calculate nTS. The default value is automatically calculated from the whole cleaned vector (see *VectorToHRT: Cleaning Input* for more information) when calling `vectorToHRT` and saved in the `IL` slot of the HRTList object. 
  - `normIL`: the interval length to which the other parameters should be normalised. The default is 800 ms.
  - `normHallstrom`: Should nTS be normalised with the method by Hallstrom et al. or just based on the interval length? The default is `TRUE`.
  - `coTO`, `coTS` and `coTT`: The cut-off that should be used to calculate the reliability check for the different parameters. The default is `coTO` 0, `coTS` 2.5 and `coTT` 10.
- **`getResults`** &nbsp;&nbsp; This function returns either the HRT class or the parameter values. You can determine the output with
  - `type`: "class" returns the HRT class (system HRT0-2 or HRTA-C depending on `TT`), "parameter" returns the HRT parameters and "full" returns the HRT parameters and the p-values of the reliability check.
  - `TT`: Should TT be included in the return? The default is `TRUE`.
  - `nTS`: Switches between giving TS (default) and nTS.
  - `safe`: Per default `safe` is `TRUE` so only results that are reliable are returned. For not reliable results the function returns "NR" or, if `num` is `TRUE`, `NA`.
  - `pmax`: The cut-off of the p-value to determine reliability. Per default this is 0.05.
  - `num`: Forces the function to return numerics. Keep in mind that this is inapplicable when using `type` "class", in this case the function gives a warning and returns `NA`. With `type` "full" `num` is ignored, because in that case the result is already numeric.
  - `coTO`, `coTS` and `coTT`: Analogously to `calcAvHRT` the cut-off that should be used to determine whether the parameter values are normal. The default is `coTO` 0, `coTS` 2.5 and `coTT` 10. Be sure to give the same cut-offs to `calcAvHRT` and `getResults` if you don't use the result, otherwise the p-values won't match the results.
- **`getHRTParams`** &nbsp;&nbsp; Returns the values of the given slot of all HRT objects in the `HRTList.` This can be used to quickly list all separate HRT parameters of an `HRTList`. Although the function name focuses on the HRT parameters, it can return any other slot of the HRT objects.
- **`getPositions`** &nbsp;&nbsp; Returns the positions of the couplRRs which is identical to `HRTList@pos`.
- **`plot`** &nbsp;&nbsp; Analogously to the plot function of the HRT object HRTList objects can be plotted. This function plots the avHRT and adds the VPCSs of all HRTs as grey lines in the background.

## avHRT Object

An avHRT object is stored in an HRTList and inherits from the HRT object. It is averaged from the other HRTs in the HRTList automatically and can be recalculated with calcAvHRT.

### Slots

In addition to the HRT slots avHRT stores data about its calculation and the parameter validity:

- `av`, `orTO` and `orTS`: for more information see *HRTList: Functions (calcAvHRT)*
- `pTO`, `pTS`, `pTT` and `pnTS`: p-values from the reliability check, for more information see *Reliability Check* in [Scientific Background](background.md)
- `nRMSSD`: the RMSSD normalised to the a given heart rate, per default to 75 bpm

## vectorToHRT

This is the core function of the package. It finds VPCs, checks the respective VPCS for validity and saves the results in an HRTList object. Its parameters are

- `input`: RR interval data that should be searched for HRT. Data formatted as timestamps should be converted before using `vectorToHRT`.
- `annotations`: If no annotations are given `vectorToHRT` searches for matching patterns of interval lengths in the given vector regardless of any other information in respect to the type of beats. Therefore, the function could also save HRTs based on atrial premature complexes or other arrhythmia if the surrounding intervals match the filter rules (for more information see *Methods & Background: Filter Rules*). If annotations are given the function only checks the intervals marked to stem from ventricular beats. This leads to more accurate results and speeds up the runtime considerably. The annotations should match the beats *at the end* of the intervals. 
- `PVCAnn`: A character or string with which the VPCs are marked in the annotation vector. The default is "V".
- `normIL`: The interval length to which the other parameters should be normalised. The default is 800 ms.
- `normHallstrom`: Should nTS be normalised with the method by Hallstrom et al. or just based on the interval length? The default is `TRUE`.
- `numPreRRs` and `numPostRRs`: The number of regular intervals before and after the CPI and CMI, respectively, on which the filter rules are applied and from which TS and nTS are being calculated. The default is 5 and 15, respectively.
- `inputName`: You can give a name to the `HRTList` to easier organise your data. If no name is given, the slot is set with `NA`. 
- `minHRT`: This sets the minimal amount of HRTs that have to be found. Per default an `HRTList` is only created if the vector contains 5 or more HRTs.
- `cleaning`: To calculate `IL` and `RMSSD` the data is cleaned per default. (for more information see *VectorToHRT: Cleaning Input*).

### Cleaning Input
The `IL` and `RMSSD` can be highly biased through outliers. Since ECG data can include artefacts, especially at the end and the beginning, RHRT cleans the data before calculating these parameters. Intervals are removed if they

- are greater than 2000 or less than 300 ms or
- differ more than 20 % of their own value from the next interval.

--------

Further information can be found in the other vignettes: [synopsis](synopsis.md), [scientific background](background.md) and [example pipelines](examples.md).

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->
---
title: "RHRT: Scientific Background"
author: "Valeria Blesius"
date: "2021-08-12"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{RHRT: Scientific Background}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->

This vignette gives descriptions and references regarding the scientific background of the package. RHRT includes the possibility for non-standard methodology in order to further analyse HRT and its optimal application. The respective sources and explanations of these variations can be found here. 

--------

## Filter Rules

To ensure snippets free of any bias and containing effective VPCs, the VPCSs are filtered based on their interval lengths. The first publication to mention filter rules was [Grimm et al.][FR1]. With little variations these are the criteria that are used in the package as possible VPCSs are only saved as HRT objects if they match the following criteria:

1) Filter rules for CPI and CMI:

    * CPI must have a maximal length of 80 % 
    * CMI must have a minimal length of 120 %

    Both intervals are compared to the **reference interval (RFI)**. This interval is calculated as the mean of the preceding intervals before the coupling interval.
<br><br>

2) Filter rules for regular intervals:

    * The length has to be between 300 ms and 2000 ms
    * They must not differ more than 20 % from RFI
    * or more than 200 ms from the preceding interval

    How many preceding and following intervals of CPI and CMI are checked is based on `numPreRRs` and `numPostRRs` of `vectorToHRT`. The default is 5 and 15, respectively.  If any of the intervals do not fit the rules, the complete set is neglected.

[FR1]: https://doi.org/10.1046/j.1542-474X.2003.08206.x "Grimm et al., Heart rate turbulence following ventricular premature beats in healthy controls, 2003, Ann. Noninvas. Electro. 8 127–31"

## Normalisation of Turbulence Slope

HRT is influenced by the heart rate. While there is no clear conclusion for TO, TS values clearly positively correlate with the RR interval length (reviewed in [Blesius et al. 2020][NTS1]). Therefore, RHRT calculates `nTS` that is normalised to a fixed interval length (800 ms per default) in addition to the common TS.

Beside the heart rate, TS is biased by the number of HRTs used to calculate it (reviewed in [Blesius et al. 2020][NTS1]). While physiological reasons were suggested for this phenomenon ([Cygankiewicz et al. 2004][NTS2] and [Chen 2009][NTS3]), [Hallstrom et al. 2004][NTS4] reasoned it to be a mathematically induced relation based on the number of VPCSs as well as the number of postRRs to determine TS. They proposed a method to normalise TS in which, firstly, TS is normalised to a HR of 75 bpm (which is 800 ms interval length). Here, it makes no mathematical difference whether TS is normalised or the intervals themselves before assessing TS. Secondly, the following formula is used:

&nbsp;&nbsp;&nbsp;&nbsp; nTS = TS - ( 0.02475 * (numPostRRs-2)^0.9449 * (RMSSD / √#VPCSs) )
    
RHRT uses this normalisation per default. This can be changed with the boolean parameter `normHallstrom` in `vectorToHRT` and `calcAvHRT`.

[NTS1]: https://doi.org/10.1088/1361-6579/ab98b3 "Blesius et al., HRT assessment reviewed: a systematic review of heart rate turbulence methodology, 2020, Physiol. Meas. 41 08TR01"
[NTS2]: https://doi.org/10.1046/j.1540-8167.2004.03613.x "Cygankiewicz et al., Relationship between heart rate turbulence and heart rate, heart rate variability, and number of ventricular premature beats in coronary patients, 2004, J. Cardiovasc. Electrophysiol. 15 731–7"
[NTS3]: https://doi.org/10.1111/j.1542-474X.2009.00322.x "Chen, Impact of preceding ventricular premature beats on heart rate turbulence, 2009, Ann. Noninvas. Electro. 14 333–9"
[NTS4]: https://doi.org/10.1109/TBME.2004.828049 "Hallstrom et al., Structural relationships between measures based on heart beat intervals: potential for improved risk assessment, 2004, IEEE. Trans. Biomed. Eng. 51 1414–20"

## Reliability Check
The HRT parameter values pre se do not give any information about 1) how many VPCSs have been used to determine them and 2) how reliable the values are. However, two identical values are inherently different if one is calculated from a VPCS with a highly varying values and the other from a high amount of VPCS with hardly any variation. Still, HRT classification generally does not take this into account.

RHRT implements a reliability check to give the opportunity to only use HRT parameter values that are reliable to a desired extent. This check consists of a one-sided t-test (`t.test` of the stats package) off all separate values against the respective cut-off of the parameter. The resulting p-value implicates the possibility of the classification being true based on being the combination of average and variability of the parameter values and therefore the reliability of the averaged value.

These t-tests are being done automatically during `calcAvHRT` which is called by `vectorToHRT`. The default values of the cut-offs are 0 for `TO`, 2.5 for `TS` as well as `nTS` and 10 for `TT`.
`getResults` returns the results if reliable. However, it returns all results ignoring the reliability check via the boolean parameter `safe` and changes the p-value cut-off with `pmax` (0.05 per default).

Keep in mind that the parameter value cut-offs `coTO`, `coTS` and `coTT` are only used to compare the values and classify them. They are not related to the identically named parameters of `calcAvHRT` that are used for the t-tests.

## Calculation Order
The order in which the HRT parameters are calculated has an impact on the resulting values ([Chen 2011][CO1]). Though [Schmidt et al. 1999][CO2] proposed to first calculate an averaged tachogram and determine TS then and for TO to first assess it from the separate VPCSs and average the results afterwards, the order gets switched in some studies as reviewed in [Blesius et al. 2020][CO3]. Therefore, RHRT gives the opportunity to change the calculation order for TO and TS through the parameters `orTO` and `orTS` of `calcAvHRT`. By default the order is as suggested by Schmidt et al. Additionally with `av` you can switch between `mean` and `median` as averaging function.

[CO1]: https://doi.org/10.4081/hi.2011.e7 "Chen, Implications of turbulence slope variations in different approaches, 2011, Heart Int. 6 21–5"
[CO2]: https://doi.org/10.1016/S0140-6736(98)08428-1 "Schmidt et al., Heart-rate turbulence after ventricular premature beats as a predictor of mortality after acute myocardial infarction, 1999, The Lancet 353 1390–6"
[CO3]: https://doi.org/10.1088/1361-6579/ab98b3 "Blesius et al., HRT assessment reviewed: a systematic review of heart rate turbulence methodology, 2020, Physiol. Meas. 41 08TR01"

--------

Further information can be found in the other vignettes: [synopsis](synopsis.md), [objects & functions](objects_functions.md) and [example pipelines](examples.md).

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->
---
title: "RHRT: Example Pipelines"
author: "Valeria Blesius"
date: "2021-09-17"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{RHRT: Example Pipelines}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

--------

## Determining the HRT class of a person

The main focus of the package is to determine the HRT parameters or class of a person by a long-term ECG measurement. Load the data as a numeric vector and use `vectorToHRT` to find HRTs, then `plot` and `getResults` to check the HRT:


```r
library("RHRT")
hrtl <- vectorToHRT(testdataLong) # create the HRTList
plot(hrtl, main = "Zoomed in Tachogram") # plot the HRTs and check the variability
```

![](examples_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
getResults(hrtl) # get the averaged HRT parameters
```

```
## [1] "NR"
```

The results do not pass the reliability check so we get "NR" instead of an HRT class. The plot shows that firstly TO is near to zero and secondly there is a high variability in the VPCSs. We can go deeper into the data by checking the exact parameters (including TT as an additional hint to the person's status) and zooming out of the plot:


```r
round(
  getResults(hrtl, "full", TT = TRUE),
digits = 2) # get the parameters and p-values of the variability check
```

```
##    TO    TS    TT   pTO   pTS   pTT 
## -0.48 26.82  3.00  0.38  0.00  0.00
```

```r
plot(hrtl, cropped = FALSE, main = "Full Tachogram") # plot the full VPCSs
```

![](examples_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

As expected TO is not reliable with a p-value over 0.05. The VPCSs still seem to fluctuate a lot. We can can get a picture of the individual TO values by using `getHRTParams`:


```r
tos <- getHRTParams(hrtl, "TO")
tos
```

```
##  [1]  -3.08977601  -2.91673793  -3.86304382 -11.81622053   1.48341915
##  [6]  -2.04437316   2.02537861   6.54482521   0.04640218   1.24190434
## [11]  10.01034320  -3.41808350
```

```r
summary(tos)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -11.816  -3.172  -0.999  -0.483   1.619  10.010
```

```r
par(mar=c(0, 3, 0, 0))
boxplot(tos)
```

![](examples_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

These results can help to come to a well-founded decision on whether to classify the patient as HRT0/HRTA and trust the TO value or rather classify them conservatively as HRT1/HRTB. 

## Comparing HRT results with different methodological parameters

This is an example how the package can be used to analyse the HRT methodology. For instance, we can compare the difference in `TO` and `TS` values when the order of the calculation steps are switched.


```r
library("RHRT")
hrtl <- vectorToHRT(testdataLong)
getResults(hrtl, type = "parameter", safe = FALSE)
```

```
##         TO         TS 
## -0.4829969 26.8151839
```

```r
hrtl@avHRT <- calcAvHRT(hrtl, orTO = "avBefore", orTS = "avAfter")
getResults(hrtl, type = "parameter", safe = FALSE)
```

```
##         TO         TS 
## -0.6618164 35.2097614
```

--------

Further information can be found in the other vignettes: [synopsis](synopsis.md), [objects & functions](objects_functions.md) and [scientific background](background.md).

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# RHRT

RHRT provides tools for a Heart Rate Turbulence analysis of RR interval data. It can either find the ventricular premature complexes (VPCs) via a set of filter rules or can use annotation data to only check the beats with the correct annotation. VPC snippets are filtered for validity and HRT parameters calculated.

In addition to standard calculation methods the package allows to modify the filter and calculation parameters. It is therefore not only helpful to identify HRT classes of measurements for risk assessment but also for assessment of the methodology itself.

For more information please check the vignette of the package.
For more information about HRT have a look into the [original publication by Schmidt et al.](https://doi.org/10.1016/S0140-6736(98)08428-1) or our [review](https://doi.org/10.1088/1361-6579/ab98b3) with focus on the methodology.

## Installation

You can install the released version of RHRT from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RHRT")
```

And the development version from [GitHub](https://github.com/) with:

``` r
install.packages("devtools") # if not already done
devtools::install_github("VBlesius/RHRT/RHRT")
```
## Example

The general workflow of RHRT is to scan your interval data for HRT and check the results via HRT class and plot:

```{r example, fig.width=7, fig.height=4}
library(RHRT)
## scan your interval data and save the results as an HRTList
hrtl <- vectorToHRT(testdataLong)
## get the HRT class of your data
getResults(hrtl, type = "class")
## have a look at the data and the parameters
plot(hrtl)
```

## Data

Data to test the package can be found on [Physionet](https://physionet.org/). Via the [WFDB Toolkit](https://physionet.org/content/wfdb/10.6.2/) ECG data can be downloaded and/or converted, for example:

```{bash, eval = FALSE}
ann2rr -r chf2db/chf201 -a ecg -i s3 -w > ~/some/path/chf201.csv
```

Then load the data and use RHRT to find VPCSs:
``` {r, eval = FALSE}
chf201 <- read.table("~/some/path/chf201.csv")
hrtl <- RHRT::vectorToHRT(chf201[[1]]*1000, ann = chf201[[2]])
```

More example workflows can be found in the vignettes: [Quick start guide](vignettes/synopsis.md), [Objects & Functions](vignettes/objects_functions.md), [Scientific Background](vignettes/background.md) and [Example Pipelines](vignettes/examples.md).

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->---
title: "RHRT: Objects and Functions"
author: "Valeria Blesius"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{RHRT: Objects and Functions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This vignette gives a detailed overview of all objects and functions provided by RHRT.

--------

## HRT Object

### Slots

An HRT object saves the data of one VPCS and its HRT results. It consists of the following slots:

- **Intervals**
  - `preRRs`: preceding regular intervals, 5 per default
  - `couplRR`: CPI
  - `compRR`: CMI
  - `postRRs`: following regular intervals, 15 per default

- **HRT Parameters**
  - `TO`
  - `TS`
  - `TT`
  - `nTS`: normalised TS, for more Information see chapter *Normalisation of Turbulence Slope* in [Scientific Background](background.md)

- **Line coefficients**
  - `intercept`: The intercept of the TS regression line that is needed for plotting
  - `nintercept`: Analogously, the intercept of `nTS`

### Functions

- **`getRRs`** &nbsp;&nbsp; This function returns all intervals saved in the HRT object. This if helpful for i.e. calculating an averaged VPCS out of a set of HRT objects.

- **`plot`** &nbsp;&nbsp; Per default `plot` displays a zoomed in plot of the VPCS with highlighted TO and TS and a legend given the rounded HRT parameter values. In addition to the common parameters of graphics.plot the method accepts the following parameters:
  - `cropped`: switches between showing a zoomed in version and the full VPCS including CPI and CMI, the default is `TRUE`
  - `add`: adds the plot to the current one
  - `TT`: highlights TT and includes it to the legend
  - `paramsLegend`: switches between showing and hiding the parameter values in the legend, the default is `TRUE`
  - `colTO`, `colTS` and `colTT` determine the colour in which the different parameters are highlighted in the plot (red, blue and purple per default).

## HRTList Object

An HRTList object is created when `vectorToHRT` is called. All HRT parameters are calculated automatically in this process and can be called with the `getHRTParams`-methods.

### Slots

The HRTList object sums up all HRTs found in the given dataset. The slots are:

- `name`: The name of the object if given to `vectorToHRT`
- `IL`: the average length of all intervals in the given cleaned vector (see *vectorToHRT: Cleaning Input* for more information) 
- `pos`: the indices of the CPIs as found in the given vector
- `HRTs`: list of all HRT objects found
- `avHRT`: an avHRT object averaged from all HRTs
- `RMSSD`: the HRV parameter RMSSD calculated from all intervals in the given cleaned vector (see *vectorToHRT: Cleaning Input* for more information): RMSSD is the square root of the mean of the squared respective differences of the successive RR intervals, in R calculated as `sqrt(mean(diff(intervals)^2))`

### Functions

- **`calcAvHRT`** &nbsp;&nbsp; The function calculates the parameters of the averaged HRT. This is called automatically when using `vectorToHRT` with default parameters. If the avHRT should be calculated differently the following options are available:
  - `av`: The function with which averaging should be done: `mean` or `median`.
  - `orTO` and `orTS`: sets the order in which TO and TS should be calculated. With `avAfter` the parameter is assessed separately for every HRT and averaged afterwards, with `avBefore` the intervals of all VPCSs are averaged first and the parameter is assessed afterwards. The default is `avAfter` for `orTO` and `avBefore` for `orTS`.
  - `IL`: the average interval length that is needed to calculate nTS. The default value is automatically calculated from the whole cleaned vector (see *VectorToHRT: Cleaning Input* for more information) when calling `vectorToHRT` and saved in the `IL` slot of the HRTList object. 
  - `normIL`: the interval length to which the other parameters should be normalised. The default is 800 ms.
  - `normHallstrom`: Should nTS be normalised with the method by Hallstrom et al. or just based on the interval length? The default is `TRUE`.
  - `coTO`, `coTS` and `coTT`: The cut-off that should be used to calculate the reliability check for the different parameters. The default is `coTO` 0, `coTS` 2.5 and `coTT` 10.
- **`getResults`** &nbsp;&nbsp; This function returns either the HRT class or the parameter values. You can determine the output with
  - `type`: "class" returns the HRT class (system HRT0-2 or HRTA-C depending on `TT`), "parameter" returns the HRT parameters and "full" returns the HRT parameters and the p-values of the reliability check.
  - `TT`: Should TT be included in the return? The default is `TRUE`.
  - `nTS`: Switches between giving TS (default) and nTS.
  - `safe`: Per default `safe` is `TRUE` so only results that are reliable are returned. For not reliable results the function returns "NR" or, if `num` is `TRUE`, `NA`.
  - `pmax`: The cut-off of the p-value to determine reliability. Per default this is 0.05.
  - `num`: Forces the function to return numerics. Keep in mind that this is inapplicable when using `type` "class", in this case the function gives a warning and returns `NA`. With `type` "full" `num` is ignored, because in that case the result is already numeric.
  - `coTO`, `coTS` and `coTT`: Analogously to `calcAvHRT` the cut-off that should be used to determine whether the parameter values are normal. The default is `coTO` 0, `coTS` 2.5 and `coTT` 10. Be sure to give the same cut-offs to `calcAvHRT` and `getResults` if you don't use the result, otherwise the p-values won't match the results.
- **`getHRTParams`** &nbsp;&nbsp; Returns the values of the given slot of all HRT objects in the `HRTList.` This can be used to quickly list all separate HRT parameters of an `HRTList`. Although the function name focuses on the HRT parameters, it can return any other slot of the HRT objects.
- **`getPositions`** &nbsp;&nbsp; Returns the positions of the couplRRs which is identical to `HRTList@pos`.
- **`plot`** &nbsp;&nbsp; Analogously to the plot function of the HRT object HRTList objects can be plotted. This function plots the avHRT and adds the VPCSs of all HRTs as grey lines in the background.

## avHRT Object

An avHRT object is stored in an HRTList and inherits from the HRT object. It is averaged from the other HRTs in the HRTList automatically and can be recalculated with calcAvHRT.

### Slots

In addition to the HRT slots avHRT stores data about its calculation and the parameter validity:

- `av`, `orTO` and `orTS`: for more information see *HRTList: Functions (calcAvHRT)*
- `pTO`, `pTS`, `pTT` and `pnTS`: p-values from the reliability check, for more information see *Reliability Check* in [Scientific Background](background.md)
- `nRMSSD`: the RMSSD normalised to the a given heart rate, per default to 75 bpm

## vectorToHRT

This is the core function of the package. It finds VPCs, checks the respective VPCS for validity and saves the results in an HRTList object. Its parameters are

- `input`: RR interval data that should be searched for HRT. Data formatted as timestamps should be converted before using `vectorToHRT`.
- `annotations`: If no annotations are given `vectorToHRT` searches for matching patterns of interval lengths in the given vector regardless of any other information in respect to the type of beats. Therefore, the function could also save HRTs based on atrial premature complexes or other arrhythmia if the surrounding intervals match the filter rules (for more information see *Methods & Background: Filter Rules*). If annotations are given the function only checks the intervals marked to stem from ventricular beats. This leads to more accurate results and speeds up the runtime considerably. The annotations should match the beats *at the end* of the intervals. 
- `PVCAnn`: A character or string with which the VPCs are marked in the annotation vector. The default is "V".
- `normIL`: The interval length to which the other parameters should be normalised. The default is 800 ms.
- `normHallstrom`: Should nTS be normalised with the method by Hallstrom et al. or just based on the interval length? The default is `TRUE`.
- `numPreRRs` and `numPostRRs`: The number of regular intervals before and after the CPI and CMI, respectively, on which the filter rules are applied and from which TS and nTS are being calculated. The default is 5 and 15, respectively.
- `inputName`: You can give a name to the `HRTList` to easier organise your data. If no name is given, the slot is set with `NA`. 
- `minHRT`: This sets the minimal amount of HRTs that have to be found. Per default an `HRTList` is only created if the vector contains 5 or more HRTs.
- `cleaning`: To calculate `IL` and `RMSSD` the data is cleaned per default. (for more information see *VectorToHRT: Cleaning Input*).

### Cleaning Input
The `IL` and `RMSSD` can be highly biased through outliers. Since ECG data can include artefacts, especially at the end and the beginning, RHRT cleans the data before calculating these parameters. Intervals are removed if they

- are greater than 2000 or less than 300 ms or
- differ more than 20 % of their own value from the next interval.

--------

Further information can be found in the other vignettes: [synopsis](synopsis.md), [scientific background](background.md) and [example pipelines](examples.md).

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->---
title: "RHRT: Synopsis for the hasty (Quick-start guide)"
author: "Valeria Blesius"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{RHRT: Quick-start guide}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The RHRT package helps you assess **Heart Rate Turbulence** (HRT) in RR intervals and calculate **turbulence onset** (TO), **slope** (TS) and **timing** (TT). It can plot the tachograms and checks the results for reliability. The **ventricular premature beats** (VPCs) with **coupling** (CPI) and **compensatory interval** (CMI) can either be given with annotations or found on the basis of the filter rules as first described  by [Grimm et al. 2003](https://doi.org/10.1046/j.1542-474X.2003.08206.x). The type of average and order of calculation for all parameters can be set.

This vignette sums up the most common functions and parameters needed when using RHRT.

--------

## Loading package and data

```{r}
library("RHRT")
# testdataLong is a numeric vector of RR intervals in msec
data("testdataLong", package = "RHRT")
ints <- testdataLong
# testdataLong_Ann is a character vector of annotations corresponding to testdataLong
data("testdataLong_Ann", package = "RHRT")
ann <- testdataLong_Ann
```

## Checking interval data for HRTs

The **core function** of RHRT is `vectorToHRT` that finds valid VPCs in RR intervals and returns an `HRTList` object (see *HRTList object* in [Objects & Functions](objects_functions.md) for more information):

```{r}
hrtl <- vectorToHRT(ints) 
```

Every RR interval sequence that matches the needed interval lengths is considered to be a coupling and compensatory interval of a VPC, which can lead to wrong matches. If your data is annotated, you can provide the **annotation data** with the parameters `annotations` and `PVCAnn`.

```{r}
hrtl <- vectorToHRT(ints, annotations = ann, PVCAnn = "V")
```

Other parameters are:

* `numPreRRs` & `numPostRRs` are used to modify the **filter rules** to find HRTs (number of intervals before and after the VPC that have to match the filter criteria).
* `minHRT` is the **minimal number of HRTs** needed to calculate HRT / create a HRTList
* `normHallstrom` defines whether TS should be **normalised** with the method of Hallstrom et al. (see the chapter *Normalisation of Turbulence Slope* in the [scientific background](background.md) for more information). 

## Getting HRT parameters or class
```{r}
getResults(hrtl) # get the HRT class of the data
```

Per default `getResults` checks whether all needed HRT parameters can be calculated reliably. This is done via a t-test per parameter value (for more information see chapter *Reliability Check* in the [scientific background](background.md) vignette). If any of the parameter values is **not reliable** `getResults` returns NR (not reliable). 

```{r}
getResults(hrtl, safe = FALSE) # get the HRT class without safety check
```

In addition to the classification system HRT0-2 RHRT implements **HRTA-C** that is based on the three parameters TO, TS and TT. 

```{r}
getResults(hrtl, safe = FALSE, TT = TRUE) # include TT
```

With the parameter `type` you can choose between getting only the HRT **class**, all **parameter values** or the parameter values with the corresponding **p-values** (types "class", "parameter" or "full", respectively).

```{r}
getResults(hrtl, type = "parameter", TT = TRUE) # get the averaged HRT parameters
```

Other parameters are:

* `nTS`: the **normalised TS** is returned or used for classification instead of TS.
* `num`: forces the function to return **numerics** when using `type = parameter`. Depending on the results and your setting of `type` the `getResults` returns characters or numerics.
* `pmax`: changes the needed **significance level** for the parameters to be reliable.

## Plotting

```{r, fig.width=7, fig.height=4}
plot(hrtl, TT = TRUE) # plots the averaged VPCS and all underlying VPCSs in background
```

Per default the VPCS is **zoomed** in. If you want to also see the CPI and CMI use `cropped = FALSE`.

```{r, fig.width=7, fig.height=4}
plot(hrtl, cropped = FALSE) # shows also coupling and compensatory interval
```

--------

Further information can be found in the other vignettes about the [objects and functions](objects_functions.md), the [scientific background](background.md) or with [example pipelines](examples.md).

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->---
title: "RHRT: Example Pipelines"
author: "Valeria Blesius"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{RHRT: Example Pipelines}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

--------

## Determining the HRT class of a person

The main focus of the package is to determine the HRT parameters or class of a person by a long-term ECG measurement. Load the data as a numeric vector and use `vectorToHRT` to find HRTs, then `plot` and `getResults` to check the HRT:

```{r, fig.width=7, fig.height=4}
library("RHRT")
hrtl <- vectorToHRT(testdataLong) # create the HRTList
plot(hrtl, main = "Zoomed in Tachogram") # plot the HRTs and check the variability
getResults(hrtl) # get the averaged HRT parameters
```

The results do not pass the reliability check so we get "NR" instead of an HRT class. The plot shows that firstly TO is near to zero and secondly there is a high variability in the VPCSs. We can go deeper into the data by checking the exact parameters (including TT as an additional hint to the person's status) and zooming out of the plot:

```{r, fig.width=7, fig.height=4}
round(
  getResults(hrtl, "full", TT = TRUE),
digits = 2) # get the parameters and p-values of the variability check
plot(hrtl, cropped = FALSE, main = "Full Tachogram") # plot the full VPCSs
```

As expected TO is not reliable with a p-value over 0.05. The VPCSs still seem to fluctuate a lot. We can can get a picture of the individual TO values by using `getHRTParams`:

```{r}
tos <- getHRTParams(hrtl, "TO")
tos
summary(tos)
par(mar=c(0, 3, 0, 0))
boxplot(tos)
```

These results can help to come to a well-founded decision on whether to classify the patient as HRT0/HRTA and trust the TO value or rather classify them conservatively as HRT1/HRTB. 

## Comparing HRT results with different methodological parameters

This is an example how the package can be used to analyse the HRT methodology. For instance, we can compare the difference in `TO` and `TS` values when the order of the calculation steps are switched.

```{r}
library("RHRT")
hrtl <- vectorToHRT(testdataLong)
getResults(hrtl, type = "parameter", safe = FALSE)
hrtl@avHRT <- calcAvHRT(hrtl, orTO = "avBefore", orTS = "avAfter")
getResults(hrtl, type = "parameter", safe = FALSE)
```

--------

Further information can be found in the other vignettes: [synopsis](synopsis.md), [objects & functions](objects_functions.md) and [scientific background](background.md).

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->---
title: "RHRT: Scientific Background"
author: "Valeria Blesius"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{RHRT: Scientific Background}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->

This vignette gives descriptions and references regarding the scientific background of the package. RHRT includes the possibility for non-standard methodology in order to further analyse HRT and its optimal application. The respective sources and explanations of these variations can be found here. 

--------

## Filter Rules

To ensure snippets free of any bias and containing effective VPCs, the VPCSs are filtered based on their interval lengths. The first publication to mention filter rules was [Grimm et al.][FR1]. With little variations these are the criteria that are used in the package as possible VPCSs are only saved as HRT objects if they match the following criteria:

1) Filter rules for CPI and CMI:

    * CPI must have a maximal length of 80 % 
    * CMI must have a minimal length of 120 %

    Both intervals are compared to the **reference interval (RFI)**. This interval is calculated as the mean of the preceding intervals before the coupling interval.
<br><br>

2) Filter rules for regular intervals:

    * The length has to be between 300 ms and 2000 ms
    * They must not differ more than 20 % from RFI
    * or more than 200 ms from the preceding interval

    How many preceding and following intervals of CPI and CMI are checked is based on `numPreRRs` and `numPostRRs` of `vectorToHRT`. The default is 5 and 15, respectively.  If any of the intervals do not fit the rules, the complete set is neglected.

[FR1]: https://doi.org/10.1046/j.1542-474X.2003.08206.x "Grimm et al., Heart rate turbulence following ventricular premature beats in healthy controls, 2003, Ann. Noninvas. Electro. 8 127–31"

## Normalisation of Turbulence Slope

HRT is influenced by the heart rate. While there is no clear conclusion for TO, TS values clearly positively correlate with the RR interval length (reviewed in [Blesius et al. 2020][NTS1]). Therefore, RHRT calculates `nTS` that is normalised to a fixed interval length (800 ms per default) in addition to the common TS.

Beside the heart rate, TS is biased by the number of HRTs used to calculate it (reviewed in [Blesius et al. 2020][NTS1]). While physiological reasons were suggested for this phenomenon ([Cygankiewicz et al. 2004][NTS2] and [Chen 2009][NTS3]), [Hallstrom et al. 2004][NTS4] reasoned it to be a mathematically induced relation based on the number of VPCSs as well as the number of postRRs to determine TS. They proposed a method to normalise TS in which, firstly, TS is normalised to a HR of 75 bpm (which is 800 ms interval length). Here, it makes no mathematical difference whether TS is normalised or the intervals themselves before assessing TS. Secondly, the following formula is used:

&nbsp;&nbsp;&nbsp;&nbsp; nTS = TS - ( 0.02475 * (numPostRRs-2)^0.9449 * (RMSSD / √#VPCSs) )
    
RHRT uses this normalisation per default. This can be changed with the boolean parameter `normHallstrom` in `vectorToHRT` and `calcAvHRT`.

[NTS1]: https://doi.org/10.1088/1361-6579/ab98b3 "Blesius et al., HRT assessment reviewed: a systematic review of heart rate turbulence methodology, 2020, Physiol. Meas. 41 08TR01"
[NTS2]: https://doi.org/10.1046/j.1540-8167.2004.03613.x "Cygankiewicz et al., Relationship between heart rate turbulence and heart rate, heart rate variability, and number of ventricular premature beats in coronary patients, 2004, J. Cardiovasc. Electrophysiol. 15 731–7"
[NTS3]: https://doi.org/10.1111/j.1542-474X.2009.00322.x "Chen, Impact of preceding ventricular premature beats on heart rate turbulence, 2009, Ann. Noninvas. Electro. 14 333–9"
[NTS4]: https://doi.org/10.1109/TBME.2004.828049 "Hallstrom et al., Structural relationships between measures based on heart beat intervals: potential for improved risk assessment, 2004, IEEE. Trans. Biomed. Eng. 51 1414–20"

## Reliability Check
The HRT parameter values pre se do not give any information about 1) how many VPCSs have been used to determine them and 2) how reliable the values are. However, two identical values are inherently different if one is calculated from a VPCS with a highly varying values and the other from a high amount of VPCS with hardly any variation. Still, HRT classification generally does not take this into account.

RHRT implements a reliability check to give the opportunity to only use HRT parameter values that are reliable to a desired extent. This check consists of a one-sided t-test (`t.test` of the stats package) off all separate values against the respective cut-off of the parameter. The resulting p-value implicates the possibility of the classification being true based on being the combination of average and variability of the parameter values and therefore the reliability of the averaged value.

These t-tests are being done automatically during `calcAvHRT` which is called by `vectorToHRT`. The default values of the cut-offs are 0 for `TO`, 2.5 for `TS` as well as `nTS` and 10 for `TT`.
`getResults` returns the results if reliable. However, it returns all results ignoring the reliability check via the boolean parameter `safe` and changes the p-value cut-off with `pmax` (0.05 per default).

Keep in mind that the parameter value cut-offs `coTO`, `coTS` and `coTT` are only used to compare the values and classify them. They are not related to the identically named parameters of `calcAvHRT` that are used for the t-tests.

## Calculation Order
The order in which the HRT parameters are calculated has an impact on the resulting values ([Chen 2011][CO1]). Though [Schmidt et al. 1999][CO2] proposed to first calculate an averaged tachogram and determine TS then and for TO to first assess it from the separate VPCSs and average the results afterwards, the order gets switched in some studies as reviewed in [Blesius et al. 2020][CO3]. Therefore, RHRT gives the opportunity to change the calculation order for TO and TS through the parameters `orTO` and `orTS` of `calcAvHRT`. By default the order is as suggested by Schmidt et al. Additionally with `av` you can switch between `mean` and `median` as averaging function.

[CO1]: https://doi.org/10.4081/hi.2011.e7 "Chen, Implications of turbulence slope variations in different approaches, 2011, Heart Int. 6 21–5"
[CO2]: https://doi.org/10.1016/S0140-6736(98)08428-1 "Schmidt et al., Heart-rate turbulence after ventricular premature beats as a predictor of mortality after acute myocardial infarction, 1999, The Lancet 353 1390–6"
[CO3]: https://doi.org/10.1088/1361-6579/ab98b3 "Blesius et al., HRT assessment reviewed: a systematic review of heart rate turbulence methodology, 2020, Physiol. Meas. 41 08TR01"

--------

Further information can be found in the other vignettes: [synopsis](synopsis.md), [objects & functions](objects_functions.md) and [example pipelines](examples.md).

<!---
# Part of RHRT: R package to assess Heart Rate Turbulence from RR interval data 
# Copyright (C) 2021 Valeria Blesius

# RHRT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only.

# RHRT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RHRT.  If not, see <https://www.gnu.org/licenses/>.
-->% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conversion_to_HRT.R
\name{checkForHRT}
\alias{checkForHRT}
\title{Checks RR-intervals for HRT criteria and returns an HRT object}
\usage{
checkForHRT(intervals, numPreRRs = c_numPreRRs, numPostRRs = c_numPostRRs)
}
\arguments{
\item{intervals}{(Numeric vector) RR intervals in ms}

\item{numPreRRs}{(Numeric) Number of RRs before the coupling interval that are used for filtering}

\item{numPostRRs}{(Numeric) Number of RRs after the compensatory interval that are used for filtering}
}
\value{
(HRT) A single HRT object or NULL
}
\description{
Checks RR-intervals for HRT criteria and returns an HRT object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/avHRT.R
\docType{class}
\name{avHRT}
\alias{avHRT}
\alias{initialize,avHRT-method}
\title{S4 class to represent an avHRT object}
\usage{
\S4method{initialize}{avHRT}(
  .Object,
  av = mean,
  orTO = "avAfter",
  orTS = "avBefore",
  pTO = NA_real_,
  pTS = NA_real_,
  pTT = NA_real_,
  pnTS = NA_real_,
  nRMSSD = NA_real_,
  couplRR = NA_real_,
  compRR = NA_real_,
  preRRs = NA_real_,
  postRRs = NA_real_
)
}
\arguments{
\item{.Object}{The name of the class}

\item{av}{(Function) Type of averaging, either mean or median}

\item{orTO}{(Character) Order in which TO was calculated,
either "avAfter" (assessment of parameter and averaging)
or "avBefore" (averaging of the VPCSs and assessment of parameter)}

\item{orTS}{(Character) Order in which TS was calculated,
either "avAfter" (assessment of parameter and averaging)
or "avBefore" (averaging of the VPCSs and assessment of parameter)}

\item{pTO}{(Numeric) p-value of t-test checking the validity of TO}

\item{pTS}{(Numeric) p-value of t-test checking the validity of TS}

\item{pTT}{(Numeric) p-value of t-test checking the validity of TT}

\item{pnTS}{(Numeric) p-value of t-test checking the validity of normalised TS}

\item{nRMSSD}{(Numeric) RMSSD normalised to HR}

\item{couplRR}{(Numeric) Coupling interval}

\item{compRR}{(Numeric) Compensatory interval}

\item{preRRs}{(Numeric vector) Preceding intervals}

\item{postRRs}{(Numeric vector) Following intervals}
}
\value{
(avHRT) A new avHRT object
}
\description{
This class extends the HRT class. An avHRT is the average of an HRTList and
saves the way in which it was calculated.
}
\section{Slots}{

\describe{
\item{\code{av}}{(Function) Type of averaging, either mean or median}

\item{\code{orTO}}{(Character) Order in which TO was calculated,
either "avAfter" (assessment of parameter and averaging)
or "avBefore" (averaging of the VPCSs and assessment of parameter)}

\item{\code{orTS}}{(Character) Order in which TS was calculated,
either "avAfter" (assessment of parameter and averaging)
or "avBefore" (averaging of the VPCSs and assessment of parameter)}

\item{\code{pTO}}{(Numeric) p-value of t-test checking the validity of TO}

\item{\code{pTS}}{(Numeric) p-value of t-test checking the validity of TS}

\item{\code{pTT}}{(Numeric) p-value of t-test checking the validity of TT}

\item{\code{pnTS}}{(Numeric) p-value of t-test checking the validity of normalised TS}

\item{\code{nRMSSD}}{(Numeric) RMSSD normalised to HR}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRT.R
\docType{class}
\name{HRT}
\alias{HRT}
\alias{initialize,HRT-method}
\title{S4 class to represent an HRT object}
\usage{
\S4method{initialize}{HRT}(
  .Object,
  couplRR = NA_real_,
  compRR = NA_real_,
  preRRs = NA_real_,
  postRRs = NA_real_
)
}
\arguments{
\item{.Object}{(Character) The name of the class}

\item{couplRR}{(Numeric) Coupling interval}

\item{compRR}{(Numeric) Compensatory interval}

\item{preRRs}{(Numeric vector) Preceding intervals}

\item{postRRs}{(Numeric vector) Following intervals}
}
\value{
(HRT) A new HRT object
}
\description{
This class specifies an object to save the lengths of intervals surrounding a
premature ventricular beat. It saves the HRT parameters turbulence onset (TO),
slope (TS) and timing (TT) after calculation as well as the coefficients of an
ab-line used for the plot.
TS is saved after common calculation and after normalising.
}
\section{Slots}{

\describe{
\item{\code{couplRR}}{(Numeric) Coupling interval}

\item{\code{compRR}}{(Numeric) Compensatory interval}

\item{\code{preRRs}}{(Numeric vector) Preceding intervals}

\item{\code{postRRs}}{(Numeric vector) Following intervals}

\item{\code{TO}}{(Numeric) Turbulence onset}

\item{\code{TS}}{(Numeric) Turbulence slope}

\item{\code{TT}}{(Numeric) Turbulence timing}

\item{\code{intercept}}{(Numeric) Intercept of regression line of TS}

\item{\code{nTS}}{(Numeric) Normalised Turbulence slope}

\item{\code{nintercept}}{(Numeric) Intercept of regression line of nTS}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roll.R
\name{roll}
\alias{roll}
\title{Apply method on sliding window}
\usage{
roll(intervals, width, fun, ...)
}
\arguments{
\item{intervals}{vector}

\item{width}{window size}

\item{fun}{function to be applied}

\item{...}{additional arguments for FUN}
}
\value{
(list) List with return values of fun for each window
}
\description{
Applies a given function on a vector by rolling over it with a sliding window
mechanism.
}
\details{
This method was inspired by the function "wapply" by A. N. Spiess, University
Hospital Hamburg-Eppendorf (https://rmazing.wordpress.com/2013/04/23/wapply-a-faster-but-less-functional-rollapply-for-vector-setups/),
but adjusted for this package to speed it up.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conversion_to_HRT.R
\name{checkInput}
\alias{checkInput}
\title{Checks data input for compatibility}
\usage{
checkInput(input, numSnippet, label)
}
\arguments{
\item{input}{(Numeric vector) RR intervals in ms}

\item{numSnippet}{(Numeric) number of RRs in the the HRT snippet}

\item{label}{(Character) Name of the data given and formatted for output}
}
\value{
No return value, possibly throws errors/warnings
}
\description{
Checks data input for compatibility
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRTList.R
\name{getPositions}
\alias{getPositions}
\alias{getPositions,HRTList-method}
\title{Get positions of PVCs}
\usage{
getPositions(HRTListObj)

\S4method{getPositions}{HRTList}(HRTListObj)
}
\arguments{
\item{HRTListObj}{(HRTList object)}
}
\value{
No return value, possibly throws errors/warnings
}
\description{
Returns the positions of all ventricular premature complexes (VPCs) and
accordingly the coupling intervals that were found in the given vector
when the HRTList was created.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRT.R, R/HRTList.R
\name{checkValidity}
\alias{checkValidity}
\alias{checkValidity,HRT-method}
\alias{checkValidity,HRTList-method}
\title{Checks whether slots are set}
\usage{
checkValidity(x, ...)

\S4method{checkValidity}{HRT}(x)

\S4method{checkValidity}{HRTList}(x, av = FALSE, pos = FALSE)
}
\arguments{
\item{x}{HRTList}

\item{...}{Other parameters}

\item{av}{(Boolean) Should avHRT be checked?}

\item{pos}{(Boolean) Should pos be checked?}
}
\value{
No return value, possibly throws errors/warnings

No return value, possibly throws errors
}
\description{
Checks whether slots are set

Check for HRTList class
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRT.R
\name{calcTS}
\alias{calcTS}
\alias{calcTS,HRT-method}
\title{Calculate TS parameters}
\usage{
calcTS(HRTObj, normalising = FALSE, IL = c_normIL, normIL = c_normIL)

\S4method{calcTS}{HRT}(HRTObj, normalising = FALSE, IL = c_normIL, normIL = c_normIL)
}
\arguments{
\item{HRTObj}{(HRT) The HRT object, for which TS should be calculated}

\item{normalising}{(Boolean) Should the normalised TS be calculated?}

\item{IL}{(Numeric) The overall arithmetic mean of the interval length of the
measurement to normalise TS}

\item{normIL}{(Numeric) The interval length to which TS should be normalised}
}
\value{
(HRT) An HRT object with (re)calculated TS+intercept or nTS+nintercept
}
\description{
Calculates all TS parameters (TS itself, its index TT (turbulence timing)
and the intercept for the plot) and saves them in the corresponding slots.
Can also calculate normalised TS and intercept.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conversion_to_HRT.R
\name{checkAnnotations}
\alias{checkAnnotations}
\title{Checks annotations for compatibility}
\usage{
checkAnnotations(annotations, input, PVCAnn, label)
}
\arguments{
\item{annotations}{(Character vector) Annotations matching input}

\item{input}{(Numeric vector) RR intervals in ms}

\item{PVCAnn}{(Character) Character that marks a VPC in the annotations}

\item{label}{(Character) Name of the data given and formatted for output}
}
\value{
No return value, possibly throws errors/warnings
}
\description{
Checks annotations for compatibility
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRT.R
\name{calcHRTParams}
\alias{calcHRTParams}
\alias{calcHRTParams,HRT-method}
\title{Calculate HRT parameters}
\usage{
calcHRTParams(HRTObj, IL = c_normIL, normIL = c_normIL)

\S4method{calcHRTParams}{HRT}(HRTObj, IL = c_normIL, normIL = c_normIL)
}
\arguments{
\item{HRTObj}{(HRT) The HRT object of which the parameters should be calculated}

\item{IL}{(Numeric) The overall arithmetic mean of the interval length of the
measurement to normalise TS}

\item{normIL}{(Numeric) The interval length to which TS should be normalised}
}
\value{
(HRT) An HRT object with (re)calculated HRT parameters
}
\description{
Calculates all HRT parameters needed for an HRT object
and saves them in the corresponding slots.
}
\details{
This method is a wrapper for the methods calcTO and calcTS.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRTList.R
\name{plot,HRTList-method}
\alias{plot,HRTList-method}
\title{Plot an HRTList object}
\usage{
\S4method{plot}{HRTList}(
  x,
  cropped = TRUE,
  TT = FALSE,
  pch = 20,
  xlab = "# of RR interval",
  ylab = "length of RR interval (ms)",
  paramsLegend = TRUE,
  colTO = "#ec2023",
  colTS = "#006AFF",
  colTT = "#6800DE",
  ...
)
}
\arguments{
\item{x}{HRTList}

\item{cropped}{(Boolean) Should the plot be cut to focus on the HRT parameters?
To show all points use FALSE.}

\item{TT}{(Boolean) Should Turbulence timing be marked?}

\item{pch}{(Numeric) Plotting character, for other options see graphics::var}

\item{xlab}{(Character) Label for the x axis}

\item{ylab}{(Character) Label for the y axis}

\item{paramsLegend}{(Boolean) Should the parameter values of the HRT be plotted?}

\item{colTO}{(Character) Colour used to highlight TO}

\item{colTS}{(Character) Colour used to highlight TS}

\item{colTT}{(Character) Colour used to highlight TT}

\item{...}{Other arguments in tag = value form}
}
\value{
No return value
}
\description{
Plots RR-intervals saved in the HRT objects, especially the avHRT object,
and marks the HRT parameters.
}
\note{
Please note that some graphics parameters (par) cannot be modified,
 since they are needed to be set inside the function.
}
\examples{
# You need an HRTList
hrtl <- vectorToHRT(testdataLong, testdataLong_Ann)

# Plot your HRTList and zoom out
plot(hrtl, cropped = FALSE)

# Include TT and customise it
plot(hrtl, TT = TRUE, colTT = "green", pch = 7)

# Use standard graphics parameters
## Note: Some parameters are used inside the function and cannot be set
plot(hrtl, TT = TRUE, main = "Example plot", bty = "n", cex.lab = 1.2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Testdata-description.R
\docType{data}
\name{testdataLong_Ann}
\alias{testdataLong_Ann}
\title{Long term data annotations}
\format{
A vector of characters.
}
\usage{
testdataLong_Ann
}
\description{
Artificial dummy interval data: This dataset contains the annotations
matching testdataLong.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conversion_to_HRT.R
\name{vectorToHRT}
\alias{vectorToHRT}
\title{Convert a vector to HRTList}
\usage{
vectorToHRT(
  input,
  annotations = NULL,
  PVCAnn = "V",
  normIL = c_normIL,
  normHallstrom = TRUE,
  numPreRRs = c_numPreRRs,
  numPostRRs = c_numPostRRs,
  inputName = as.character(NA),
  minHRT = 5,
  cleaning = TRUE
)
}
\arguments{
\item{input}{(Numeric vector) RR intervals in ms}

\item{annotations}{(Character vector) Annotations matching input}

\item{PVCAnn}{(Character) Character that marks a VPC in the annotations}

\item{normIL}{(Numeric) The interval length to which TS should be normalised}

\item{normHallstrom}{(Boolean) Should the normalisation of Hallstrom be used?}

\item{numPreRRs}{(Numeric) Number of RRs before the coupling interval that are used for filtering}

\item{numPostRRs}{(Numeric) Number of RRs after the compensatory interval that are used for filtering}

\item{inputName}{(String) Name of the data}

\item{minHRT}{(Numeric) Minimal number of HRTs that are needed to create an HRTList object}

\item{cleaning}{(Boolean) Should the input be roughly cleaned from artefacts before calculating IL and RMSSD?}
}
\value{
(HRTList) An HRTList object
}
\description{
Scans for heart rate turbulence in a vector of RR-intervals and returns an
HRTList object including all found HRT objects.
The HRT criteria used were published by Schmidt et al.
(more information can be found in the vignette.)
}
\examples{
# You can use annotations to give the VPC indices
# Without annotation data RHRT will find VPCs based on common filtering criteria
vectorToHRT(testdataLong, annotations = testdataLong_Ann, PVCAnn = "V")

# Find HRTs with a broader range of sinus beats before and after the VPCs
vectorToHRT(testdataLong, inputName = "Dummy Measurement", numPreRRs = 10, numPostRRs = 20)
 
# Adjust the normalisation parameters
vectorToHRT(testdataLong, testdataLong_Ann, normHallstrom = FALSE, normIL = 900)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRTList.R
\name{getResults}
\alias{getResults}
\alias{getResults,HRTList-method}
\title{Get averaged HRT parameters}
\usage{
getResults(
  HRTListObj,
  type = "class",
  TT = FALSE,
  nTS = FALSE,
  safe = TRUE,
  pmax = 0.05,
  num = FALSE,
  coTO = COTO,
  coTS = COTS,
  coTT = COTT
)

\S4method{getResults}{HRTList}(
  HRTListObj,
  type = "class",
  TT = FALSE,
  nTS = FALSE,
  safe = TRUE,
  pmax = 0.05,
  num = FALSE,
  coTO = COTO,
  coTS = COTS,
  coTT = COTT
)
}
\arguments{
\item{HRTListObj}{HRTList object}

\item{type}{(String) Determining the amount of output: 'class' gives the HRT class, 'parameter' the parameter values and 'full' additionally the p-values describing parameter reliability}

\item{TT}{(Boolean) Should TT be given?}

\item{nTS}{(Boolean) Should the normalised TS (nTS) be given or used for the determination of the HRT class?}

\item{safe}{(Boolean) Should all values be given regardless of reliability checks? Note, that 'safe' is ignored when the type is 'full'.}

\item{pmax}{(Numeric) The significance level}

\item{num}{(Boolean) Should the results be numeric? This forces the results to stay numeric, but sets not reliable values as NA, if 'safe' is TRUE. Forced numeric values cannot be combined with type 'class'.}

\item{coTO}{(Numeric) Cut-off value for TO}

\item{coTS}{(Numeric) Cut-off value for TS and nTS}

\item{coTT}{(Numeric) Cut-off value for TT}
}
\value{
(Named vector, character or numeric) Either HRT classes, HRT parameter values and/or p-values
}
\description{
Returns the HRT parameters of the HRTList. Turbulence onset is calculated for
each HRT object and then averaged, turbulence slope is calculated via
averaging the intervals of all HRT objects to one HRT object and then
estimating the maximal slope.
}
\examples{
# You need an HRTList
hrtl <- vectorToHRT(testdataLong, testdataLong_Ann)

# Get the HRT classes of your HRTList
getResults(hrtl)
getResults(hrtl, TT = TRUE)

# Get the HRT parameter values of your HRTList
getResults(hrtl, type = "parameter", TT = TRUE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRTList.R
\name{getHRTParams}
\alias{getHRTParams}
\alias{getHRTParams,HRTList-method}
\title{Extracts all values of a special slot out of a HRTList}
\usage{
getHRTParams(HRTListObj, sl)

\S4method{getHRTParams}{HRTList}(HRTListObj, sl)
}
\arguments{
\item{HRTListObj}{HRTList object}

\item{sl}{(Character) Value of a slot saved by an HRT object}
}
\value{
(numeric vector or list) Vector or list of the numerics stored in the given slot
}
\description{
Extracts all values of the given slot in each HRT of the HRTList
and returns them in a list
}
\examples{
# You need an HRTList
hrtl <- vectorToHRT(testdataLong, testdataLong_Ann)

# Get all TOs of the HRTs in your HRTList
getHRTParams(hrtl, "TO")

# You can access all slots in the HRTs
getHRTParams(hrtl, "intercept")

# If you access slots that include more than one numeric, the function returns a list
getHRTParams(hrtl, "preRRs")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRTList.R
\docType{class}
\name{HRTList}
\alias{HRTList}
\alias{initialize,HRTList-method}
\title{S4 class to represent a list of HRT objects}
\usage{
\S4method{initialize}{HRTList}(
  .Object,
  name = NA_character_,
  IL = NA_real_,
  pos = NA_real_,
  HRTs = list(),
  avHRT = new("avHRT"),
  RMSSD = NA_real_
)
}
\arguments{
\item{.Object}{(Character) The name of the class}

\item{name}{(Character) Name of the vector if given}

\item{IL}{(Numeric) Arithmetic mean of the overall interval length of the vector}

\item{pos}{(Numeric vector) Positions of premature ventricular complexes in
given input}

\item{HRTs}{(List) All HRT objects}

\item{avHRT}{(avHRT object) The average of all HRTs}

\item{RMSSD}{(Numeric) Square root of the mean of the squared successive
differences between adjacent intervals of the whole measurement}
}
\value{
(HRTList) A new HRTList object
}
\description{
This class specifies an object to save all HRT objects of a given vector. It
also saves an averaged HRT for calculation of the averaged HRT parameters and
plotting of all HRTs in a single plot.
}
\section{Slots}{

\describe{
\item{\code{name}}{(Character) Name of the vector if given}

\item{\code{IL}}{(Numeric) Arithmetic mean of the overall interval length of the vector}

\item{\code{pos}}{(Numeric vector) Positions of premature ventricular complexes in
given input}

\item{\code{HRTs}}{(List) All HRT objects}

\item{\code{avHRT}}{(avHRT object) The average of all HRTs}

\item{\code{RMSSD}}{(Numeric) Square root of the mean of the squared successive
differences between adjacent intervals of the whole measurement}
}}

\note{
After using \code{vectorToHRT} all slots in the resulting HRTList
object are set. Please do not set them manually since many functions of the
HRTList class rely on valid values assigned to the needed slots.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Testdata-description.R
\docType{data}
\name{testdataLong}
\alias{testdataLong}
\title{Long term data}
\format{
A numeric vector.
}
\usage{
testdataLong
}
\description{
Artificial dummy interval data: This dataset represents a long-term
measurement and includes 15 VPCSs that fit the HRT filter rules.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conversion_to_HRT.R
\name{cleanInput}
\alias{cleanInput}
\title{Cleans data input for further checks or calculation}
\usage{
cleanInput(input)
}
\arguments{
\item{input}{(Numeric vector) RR intervals in ms}
}
\value{
(numeric vector) Input vector without possible bias
}
\description{
Cleans data input for further checks or calculation
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRT.R
\name{calcTO}
\alias{calcTO}
\alias{calcTO,HRT-method}
\title{Calculate TO parameters}
\usage{
calcTO(HRTObj)

\S4method{calcTO}{HRT}(HRTObj)
}
\arguments{
\item{HRTObj}{(HRT) The HRT object, for which TO should be calculated}
}
\value{
(HRT) An HRT object with (re)calculated TO
}
\description{
Calculates the TO parameters and saves it in the corresponding slot
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conversion_to_HRT.R
\name{getHRTs}
\alias{getHRTs}
\title{Finds HRTs}
\usage{
getHRTs(
  intervals,
  annotations = NULL,
  PVCAnn = "V",
  numPreRRs = c_numPreRRs,
  numPostRRs = c_numPostRRs,
  numSnippet
)
}
\arguments{
\item{intervals}{(Numeric vector) RR intervals in ms}

\item{annotations}{(Character vector) Annotations matching input}

\item{PVCAnn}{(Character) Character that marks a VPC in the annotations}

\item{numPreRRs}{(Numeric) Number of RRs before the coupling interval that are used for filtering}

\item{numPostRRs}{(Numeric) Number of RRs after the compensatory interval that are used for filtering}

\item{numSnippet}{(Numeric) Number of RRs in the HRT snippet}
}
\value{
(HRTList) HRTList with only pos and HRTs set
}
\description{
Scans for HRTs in the given vector and returns an HRTList object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRT.R
\name{getRRs}
\alias{getRRs}
\alias{getRRs,HRT-method}
\title{Returns the VPCS intervals in right order}
\usage{
getRRs(HRTObj)

\S4method{getRRs}{HRT}(HRTObj)
}
\arguments{
\item{HRTObj}{HRT}
}
\value{
(numeric vector) All VPCS intervals
}
\description{
Returns the VPCS intervals in right order
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRT.R
\name{plot,HRT-method}
\alias{plot,HRT-method}
\title{Plot an HRT object}
\usage{
\S4method{plot}{HRT}(
  x,
  cropped = TRUE,
  TT = FALSE,
  pch = 20,
  xlab = "# of RR interval",
  ylab = "length of RR interval (ms)",
  plotLegend = TRUE,
  paramsLegend = TRUE,
  pvalsLegend = FALSE,
  colTO = "#ec2023",
  colTS = "#006AFF",
  colTT = "#6800DE",
  add = FALSE,
  ...
)
}
\arguments{
\item{x}{(HRT) A HRT object}

\item{cropped}{(Boolean) Should the plot be cut to focus on the HRT parameters?
To show all points use FALSE.}

\item{TT}{(Boolean) Should Turbulence timing be marked?}

\item{pch}{(Numeric) Plotting character, for other options see graphics::var}

\item{xlab}{(Character) Label for the x axis}

\item{ylab}{(Character) Label for the y axis}

\item{plotLegend}{(Boolean) Should a legend be plotted?}

\item{paramsLegend}{(Boolean) Should the parameter values of the HRT be plotted?}

\item{pvalsLegend}{(Boolean) Should the p-values of the reliability check be plotted?}

\item{colTO}{(Character) Colour used to highlight TO}

\item{colTS}{(Character) Colour used to highlight TS}

\item{colTT}{(Character) Colour used to highlight TT}

\item{add}{(Boolean) Should the given HRT be added to a plot?}

\item{...}{Other arguments in tag = value form. See graphics::par for more information.}
}
\value{
No return value
}
\description{
Plots RR-intervals saved in the HRT object and marks
turbulence onset and turbulence slope.
}
\note{
Please note that some graphics parameters (par) cannot be modified,
 since they are needed to be set inside the function.
}
\examples{
# You need an HRT object
hrt <- vectorToHRT(testdataLong, testdataLong_Ann)@HRTs[[1]]

# Plot your HRT and zoom out
plot(hrt, cropped = FALSE)

# Include TT and customise it
plot(hrt, TT = TRUE, colTT = "green", pch = 7)

# Use standard graphics parameters
## Note: Some parameters are used inside the function and cannot be set
plot(hrt, TT = TRUE, main = "Example plot", bty = "n", cex.lab = 1.2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRTList.R
\name{calcAvHRT}
\alias{calcAvHRT}
\alias{calcAvHRT,HRTList-method}
\title{Calculate an avHRT object}
\usage{
calcAvHRT(
  HRTListObj,
  av = mean,
  orTO = "avAfter",
  orTS = "avBefore",
  IL = HRTListObj@IL,
  normIL = c_normIL,
  normHallstrom = TRUE,
  coTO = COTO,
  coTS = COTS,
  coTT = COTT
)

\S4method{calcAvHRT}{HRTList}(
  HRTListObj,
  av = mean,
  orTO = "avAfter",
  orTS = "avBefore",
  IL = HRTListObj@IL,
  normIL = c_normIL,
  normHallstrom = TRUE,
  coTO = COTO,
  coTS = COTS,
  coTT = COTT
)
}
\arguments{
\item{HRTListObj}{HRTList object}

\item{av}{(Function) Type of averaging the VPCSs, either mean or median}

\item{orTO}{(Character) Order in which TO was calculated,
either "avAfter" (assessment of parameter and averaging)
or "avBefore" (averaging of the VPCSs and assessment of parameter)}

\item{orTS}{(Character) Order in which TS was calculated,
either "avAfter" (assessment of parameter and averaging)
or "avBefore" (averaging of the VPCSs and assessment of parameter)}

\item{IL}{(Numeric) The overall arithmetic mean of the interval length of the
measurement to normalise TS}

\item{normIL}{(Numeric) The interval length to which TS should be normalised}

\item{normHallstrom}{(Boolean) Should the normalisation of Hallstrom be used?}

\item{coTO}{(Numeric) Cut-off value for TO}

\item{coTS}{(Numeric) Cut-off value for TS and nTS}

\item{coTT}{(Numeric) Cut-off value for TT}
}
\value{
(avHRT) The avHRT object of the given HRTList
}
\description{
For each index the average of the intervals across all HRTs in the HRTList
is calculated and the averaged HRT returned. The type of averaging, the order
of HRT parameter assessment and interval lengths for normalising TS can be passed.
}
\details{
To eliminate other RR variability TS is commonly assessed after averaging the
VPCSs. TO is commonly first calculated from the single VPCS and then
averaged. (See 'Heart Rate Turbulence: Standards of Measurement,
Physiological Interpretation, and Clinical Use, Axel Bauer et al.,
Journal of the American College of Cardiology, Volume 52, Issue 17,
Pages 1353-1365')
}
\examples{
# You need an HRTList
hrtl <- vectorToHRT(testdataLong, testdataLong_Ann)

# Recalculate the avHRT with different normalisation
calcAvHRT(hrtl, normIL = 1000, normHallstrom = FALSE)

# Recalculate the avHRT based on a different calculation order
calcAvHRT(hrtl, orTO = "avBefore", orTS = "avAfter")

# Set custom parameter cut-offs for the reliability check
## You should keep in mind to give the same cut-offs when calling getResults()
calcAvHRT(hrtl, coTO = 0.022, coTS = 1.42, coTT = 12)

}
