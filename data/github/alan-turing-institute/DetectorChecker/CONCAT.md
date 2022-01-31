# DetectorChecker <a><img src='logo_hex.png' align="right" height="139" /></a>

Master: [![Build Status](https://travis-ci.com/alan-turing-institute/DetectorChecker.svg?token=zxQwzfsqCyEouTqXAVUn&branch=master)](https://travis-ci.com/alan-turing-institute/DetectorChecker) Develop: [![Build Status](https://travis-ci.com/alan-turing-institute/DetectorChecker.svg?token=zxQwzfsqCyEouTqXAVUn&branch=develop)](https://travis-ci.com/alan-turing-institute/DetectorChecker)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![status](https://joss.theoj.org/papers/b6b67cd22d488e7bfa42b4074bc4eda8/status.svg)](https://joss.theoj.org/papers/b6b67cd22d488e7bfa42b4074bc4eda8) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4312870.svg)](https://doi.org/10.5281/zenodo.4312870)


Created by: [Julia Brettschneider](https://github.com/ejulia17) (original R code), [Tomas Lazauskas](https://github.com/tomaslaz) (R package engineering), [Oscar Giles](https://github.com/OscartGiles) (package development) and [Wilfrid Kendall](https://github.com/WilfridSKendall) (testing and editing).


## Overview

DetectorChecker is an R package to aid in the assessment of damage to CT scanners arising from exposure to high energy radiation.
While the target application concerns CT scanners, this package can also be used to analyze screen damage arising from other sources.


## Installation

To install from github you will need to have the [devtools](https://github.com/r-lib/devtools) package installed.

In R run one of the following, depending on whether you want to build the package Vignettes, removing the # if you do not have devtools installed:

```
# install.packages("devtools")
devtools::install_github("alan-turing-institute/DetectorChecker")
```


If you want to be able to view the `DetectorChecker_user-guide` Vignette you need to install with:

```
devtools::install_github("alan-turing-institute/DetectorChecker", 
     build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```

Installing with the vignettes may be slow (~10 min)


### Development version

```
# install.packages("devtools")
devtools::install_github("alan-turing-institute/DetectorChecker", ref = "develop")
```


## WebApp

The official release of the DetectorChecker WebApp is hosted at 
<https://detectorchecker.azurewebsites.net>.

<img src="https://raw.githubusercontent.com/alan-turing-institute/DetectorChecker/master/inst/img/DetectorChecker.png" width="500" align="center">

The source code for the WebApp implementation can be found on GitHub: 
<https://github.com/alan-turing-institute/DetectorCheckerWebApp>.


## Vignette

The user guide vignette provides detailed instructions for using the package and loading specific examples. Make sure you installed the package including vignette following the instructions above (see use of `build_vignettes = TRUE` in Section Installation) and then load the package followed by the vignette command:

```
library(detectorchecker)
vignette("DetectorChecker_user-guide", package = "detectorchecker")
```

## Manual
Detailed documentation is provided as a [pdf](docs/detectorchecker_1.0.8.pdf).
See also the [DetectorChecker_user-guide](vignettes/DetectorChecker_user-guide.html) (large html file, needs download).

## Examples
DetectorChecker includes a number of example datasets for five detector types:

1. Pilatus

2. PerkinElmer

3. PerkinElmer Refurbished

4. PerkinElmer Cropped

5. Excalibur

To load an example dataset, either call:

```
library(detectorchecker)

# Initiate a PerkinElmerFull detector object
detector <-  create_detector("PerkinElmerFull")

# Path of dataset
file_path <- system.file("extdata", "PerkinElmerFull",
                        "BadPixelMap_t1.bpm.xml", 
                        package = "detectorchecker")

# Load a pixel matrix into the detector object
detector <- load_pix_matrix(detector = detector, file_path = file_path)
```

or load one of the examples by calling:

```
library(detectorchecker)
data(PerkinElmerFull_exp_1)
```

which creates an appropriate detector module and loads an example pixel dataset.

For see the full list of example datasets call

```
data(package = "detectorchecker")
```

## Citation
If you use DetectorChecker in your work please cite our package.

BibTeX:

```
  @Misc{,
    title = {{DetectorChecker}: Assessment of damage to CT scanners},
    author = {Tomas Lazauskas and Julia Brettschneider and Oscar Giles and Wilfrid Kendall},
    url = {https://github.com/alan-turing-institute/DetectorChecker},
  }
```

## Getting help
If you found a bug or need support, please submit an issue [here](https://github.com/alan-turing-institute/DetectorChecker/issues/new).

## How to contribute
We welcome contributions! If you are willing to propose new features or have bug fixes to contribute, please submit a pull request [here](https://github.com/alan-turing-institute/DetectorChecker/pulls).
# Testing builds on multiple operating systems

We want to test the package on Mac OS, Linux and Windows. The following instructions details how to run tests on all these machines, assuming you are testing on a Mac OS machine without access to Windows or Linux.

It is **very important** to test the build on Windows before merging into the Master branch, as we are not using continuous integration for Windows at the moment.


## Testing on Mac

The simplest way to test the package is to use RStudio on that machine, open the package and then call

```devtools::check()``` 

in the console. This will run all the unit tests as well as a number of other checks on the package. 

Ideally we should have no errors, warnings or notes. Any of these will be a barrier to submitting on CRAN. 

If you have access to a Windows or Linux machine you can also do the above.




## Testing Windows and Linux without access to machine

### Linux

We are using Travis CI (https://travis-ci.com/alan-turing-institute/DetectorChecker). In the repository there is a file called `.travis.yml`. When a push is made to github travis will automatically test the package on a linux machine (running the Xenial version of the operating system). You can navigate to the link above to see whether the tests passed. If so, we are running successfully running on linux. 

### Windows

To test on windows we can use the devtools package. This works by bundling the package and uploading it to http://win-builder.r-project.org/

 Simply open the package in RStudio on the branch you wish to test and call one of the following in the console.

 1. To test on the latest development version of R:
 ```
 devtools::check_win_devel(pkg = ".", args = NULL, manual = TRUE, quiet = FALSE, ...)
 ```
2. To test on the current release version of R:
```
devtools::check_win_release(pkg = ".", args = NULL, manual = TRUE, quiet = FALSE, ...)
````
3. To test on the previous major release of R:
```
check_win_oldrelease(pkg = ".", args = NULL, manual = TRUE, quiet = FALSE, ...)
```


This will send an email to the package maintainer once the results are ready (takes ~20 min). The package maintainer is specified in the `DESCRIPTION` file in the main directory of the repository. In that file there is a section that looks like this:

```
Authors@R: c(
        person("Oscar", "Giles", email = "ogiles@turing.ac.uk", role = c("aut", "cre")),
        person("Tomas", "Lazauskas", email = "tlazauskas@turing.ac.uk", role = c("aut")),
        person("Wilfrid", "Kendall", email = "W.S.Kendall@warwick.ac.uk", role = c("aut")),
        person("Julia", "Brettschneider", email = "jabrettschneider@me.com", role = c("aut")))

```

To set yourself as the package maintainer, add `"cre"` to your role. 

For more information type ```??check_win_release``` in the RStudio console. 

# How to generate documentation files


We use roxygen2 (https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html) to produce all the code documentation. All code documentation is written in the package R files above the corresponding code (e.g. every function has its documentation written directly above it). 


## Updating the man folder

We then have a folder called `man` in the repository, which is generated by opening the package in RStudio and then type:
```
devtools::check_man()
```

## Building a PDF of the documentation

To build the documentation you need to have a working latex installation.

1. Open the DetectorChecker package in R studio (click file -> open project, then select the DetectorChecker folder).

2. Update the man folder. To do this run `devtools::check_man()` in the console.

3. Make sure R's current working directory is set to the DetectorChecker folder. Then in the console run:

```
devtools::build_manual()
```

This will create a pdf of the documentation. 

# Example files

This folder contains a set of example datasets, which can be used with the DetectorCheckerWebApp and are used internally by DetectorChecker.

## How to download

Please download the full repository by going to https://github.com/alan-turing-institute/DetectorChecker and clicking `clone or download`. Download as a zip file and then navigate to the `inst/extdata` folder.

You can then use these files with the DectorCheckerWebApp (https://github.com/alan-turing-institute/DetectorCheckerWebApp).

## Files by detector type

Examples are provided for six detector types. It is possible to load datasets as a single file, or split across multiple files. All examples, except Excalibur, use a single file. For Excalibur you will need to load all examples into the WebApp simultaneously:

1. PerkinElmer_Full
    - Two example datasets:
        - BadPixelMap_t1.bmp.xml
        - BadPixelMap_t2.bmp.xml

2. PerkinElmer_Refurbished
    - Two example datasets:
        - BadPixelMap_t1.bmp.xml
        - BadPixelMap_t2.bmp.xml

3. PerkinElmer_Cropped
    - Two example datasets:
        - BadPixelMap_t1.bmp.xml
        - BadPixelMap_t2.bmp.xml

4. Excalibur
    - One example dataset, split across 5 files (load them all in to the DetectorChecker WebApp):
        - pixelmask.fem1.hdf
        - pixelmask.fem2.hdf
        - pixelmask.fem3.hdf
        - pixelmask.fem4.hdf
        - pixelmask.fem5.hdf
        - pixelmask.fem6.hdf

5. Pilatus
    - One example dataset:
        - badpixel_mask.tif

6. user-defined
	- irregular
		- layout_par_irregular.txt - layout file
		- badpixelmap_irregular - example of a dead pixels file
		
	- photographic (photographic aspect ratio without any submodes and gaps)
   		- layout_par_photographic.txt - layout file
		- examples of dead pixels files:
  			- photographic.tif
  			- badpixelmap_photographic_hom.xml
  			- badpixelmap_photographic_inhom.xml
	
---
title: 'DetectorChecker: analyzing patterns of defects in detector screens'
author:  Julia A. Brettschneider, Oscar T. Giles, Wilfrid S. Kendall, Tomas Lazauskas
date: "30 June 2020"
authors:
- affiliation: 1, 2
  name: Julia A. Brettschneider
  orcid: 0000-0003-1763-466X
- affiliation: 2
  name: Oscar T. Giles
- affiliation: 1, 2
  name: Wilfrid S. Kendall
  orcid: 0000-0001-9799-3480
- affiliation: 2
  name: Tomas Lazauskas
output: pdf_document
bibliography: paper.bib
tags:
- R
- XCT
- bad pixel map
- defective pixels
- spatial statistics
affiliations:
- index: 1
  name: Department of Statistics, University of Warwick, Coventry, United Kingdom
- index: 2
  name: The Alan Turing Institute, London, United Kingdom
---


[_DetectorChecker_](https://github.com/alan-turing-institute/DetectorChecker)
refers to an R package and an associated web application environment
[_DetectorCheckerWebApp_](https://github.com/alan-turing-institute/DetectorCheckerWebApp),
intended to help
users who need to analyze spatial patterns of defects in images.
These images can be _panel-structured_, which is to say,
composed of sub-panels arranged in  an architecture which
the user can specify in the package or in the web application.
Primary beneficiaries are intended to be individuals responsible for
high-value digital detector screens used in X-ray computerised tomography (XCT),
where defects arise due to high radiation flux.
More generally the software can be used to analyse defects in other
panel-structured arrays, for example solar panels or very large display screens.
To maximize accessibility, and to avoid any issues arising from specific software environments,
we have created a web application which provides
the principal features of the software in standalone form.
The web application also affords the possibility of engaging with our team in the analysis of time-evolving defect patterns. To the best of our knowledge, this is the first and presently the only web application and R package facilitating spatial analysis of panel-structured images.


Digital detector screens are crucial high-value components of imaging systems used throughout
modern science, medicine, and engineering systems, particularly in XCT.
The US @FDA provides information for industry on X-ray imaging devices
and lists some common laboratory tests for
evaluation of general-use X-ray imaging devices.
It also notes the applicable standard
for each modality, thus forming the basis for a maintenance schedule.
Additionally a scheduled testing framework has been proposed by the Institute of Physics and Engineering in Medicine [@IPEM2].
In September 2019 the UK National Health Service
(NHS) announced a major investment of £200m to overcome outdated equipment,
noting that a significant proportion of CT, MRI and general X-ray equipment more than 10 years old [@UK-NHS].
Thus XCT system quality concerns are very topical.

Yaffe and Rowlands [-@YaffeRowlands-1997, especially section 3.8]
point out that XCT screen quality is linked to system performance.
[_DetectorChecker_](https://github.com/alan-turing-institute/DetectorChecker)
facilitates the inclusion of screen pixel assessment in a testing framework.
Note that screen replacement or refurbishment is expensive;
regular checks of screen pixels are needed (a) to quantify screen quality
and (b) to assess possible _special causes_ of defective pixels,
using Shewhart's [-@Shewhart-1939] terminology from classic quality control.
This is best done using spatial statistics,
both in order to determine the extent to which spatial patterns of defective pixels
can be accounted for by quantifiable independent random variation
and also by describing departures from spatial randomness in ways
which are suggestive of possible explanations (for example, stress due
to screen attachment, or failure at pixel level of data readout).
Theoretical spatial statistics methodology is crucial: foundations are discussed in @ChiuStoyanKendallMecke-2013
while
@BaddeleyRubakTurner-2015 describe the [_spatstat_](https://spatstat.org/) package, an implementation of spatial statistics methods in the
R statistical computing environment [@RFoundation-2019].
[_DetectorChecker_](https://github.com/alan-turing-institute/DetectorChecker)
[@tomas_lazauskas_2020_3662233] is an R package which adapts methods from [_spatstat_](https://spatstat.org/) to the case of panel-structured images,
and analyses point patterns arising either from
individual defects or from "clumps" of defects (determined in a manner specified by the user).
The associated [web application](https://detectorchecker.azurewebsites.net/)
[@tomas_lazauskas_2020_3662235]
is based on a self-contained R environment
[_DetectorCheckerWebApp_](https://github.com/alan-turing-institute/DetectorCheckerWebApp)
together
with a [_Shiny_](https://cran.r-project.org/web/packages/shiny/index.html)[@shiny] gui,
implemented and made available _via_ _Azure_ at <https://detectorchecker.azurewebsites.net/>.
The web application exposes the
basic functionality of the [_DetectorChecker_](https://github.com/alan-turing-institute/DetectorChecker) package without the need for users to install R.
In particular the web application  can be used
to define the geometry of the sub-panels of the detector screen (which is to say, the arrangement and size of the component sub-panels),
to upload the spatial arrangement of the defective pixels
(either
directly by means of "bad pixel maps" in XML format or inferred from test images in formats including TIFF),
and then to inspect the results using the facilities offered
by the package.
The software is freely available under MIT licence, accessible _via_ the Github repositories
in the above references.
To the best of our knowledge, there is no comparable package or web application
making methods of spatial statistics available for panel-structured image data of arbitrary structure architecture.


Defects are modelled as points in an image rectangle based on overall screen dimensions.
The pattern of defects can be modelled using the web application (workflow is
summarised in Figure \ref{fig:figure1}).

![Work flow for DetectorChecker web application. Feedback/skip paths illustrate various options: refocussing attention on subsets of the point pattern (isolated pixels, small clusters, linear clusters, \ldots); working through various graphical analyses; optionally emailing data to the DetectorChecker team; and statistically fitting a variety of models of damage intensity. \label{fig:figure1}](image/flowchartDCshort.pdf)

We now discuss selected steps of the workflow using
data derived from a Pilatus detector
screen and supplied to us by Diamond Lightsource, UK.

A. The user specifies the exact architecture
of the sub-panels of the panel-structured image.
This can be done either by using a drop-down menu to specify a predetermined option,
or by uploading a file giving the specific structure of sub-panels.
The data can then be uploaded.

B. Intensity maps can be produced _via_ kernel smoothing applied to the point pattern
(replacing each defect point by the corresponding translate of a fixed kernel function and then summing).
For example, the point pattern in Figure \ref{fig:figure2}(a) yields the intensity map given in Figure \ref{fig:figure2}(b).

![](image/fig2-a.pdf){ width=42% }\qquad\qquad  ![](image/fig2-b-trim.png){ width=50% }
\begin{figure}[!h]
\caption{Pilatus detector screen: (a) Example of point pattern of defects. (b) Intensity map resulting from point pattern of defects. The intensity map draws attention to the higher intensity of defects in the corners, which is born out by inspection of the point pattern. \label{fig:figure2}}
\end{figure}


C. Measuring possible departures from _complete spatial randomness_ (CSR).
CSR is what would be expected if the point pattern was in fact
    generated by a homogeneous Poisson
    point process of constant intensity $\lambda$,
can be assessed using visual inspection of graphs
of empirical estimates of $F$, $G$ and $K$ functions as described in @ChiuStoyanKendallMecke-2013. It is clear that the
point pattern of Figure \ref{fig:figure2}(a) is strongly inhomogeneous and therefore it is
not surprising that the corresponding graphical plots
indicate clear evidence of deviation from CSR:  

* The $F$ function or "empty space function"
computes the distribution of the nearest distance to a defect point from a typical location
chosen from the image rectangle uniformly at random (and independently of the point pattern).
If the point pattern did in fact satisfy CSR then
one could consider
an empirical estimate of the $F$ function to be a random perturbation of
the theoretical $F$ function under CSR, namely
$F_\text{pois}(r)=1-\exp(-\lambda \pi r^2)$.
Figure \ref{fig:figure3}(a)
graphs different
variants of $\hat{F}$, accounting in various ways for edge-effects.
Note the clear deviation of the $\hat{F}$ empirical estimates from what would be expected under CSR,
namely the theoretical $F_\text{pois}$.


* The $G$ function computes the distribution of nearest-neighbour distances between defect points; if the point pattern did in fact satisfy CSR then
one could view
an empirical estimate of the $G$ function as a random perturbation of
the theoretical $G$ function under CSR, namely
(by a conditional probability argument)
$G_\text{pois}(r)=1-\exp(-\lambda \pi r^2)$
(actually equal to $F_\text{pois}(r)$).
See Figure \ref{fig:figure3}(b), and note the clear deviation
from the theoretical $G_\text{pois}$
of the $\hat{G}$ empirical estimates (which again account in various ways for edge-effects),
hence again suggesting deviation from CSR.

![](image/fig3-a.pdf){ width=50% } ![](image/fig3-b.pdf){ width=50% }
\begin{figure}[!h]
\caption{Pilatus detector screen: (a) $F$ plot resulting from point pattern of defects. (b) $G$ plot resulting from point pattern of defects. The graphs contrast the theoretical curve arising from CSR (blue dotted curve) with several empirical curves involving different edge corrections. The different edge-corrected empirical curves agree with each other but indicate clear divergence from the CSR curve. Here and in the following figures the graphs correspond to output from the application: graph axes have not been harmonized. \label{fig:figure3}}
\end{figure}


* The $K$ function (Ripley's $K$ function) computes the mean number of defect points within a distance $r$ of a typical defect point, viewed as a function of $r$; if the point pattern did in fact satisfy CSR then
one could view
an empirical estimate of the $K$ function as a random perturbation of the theoretical $K$ function
under CSR, namely
$K_\text{pois}(r)=\pi r^2$. See Figure \ref{fig:figure4}(a), and note the deviation from $K_\text{pois}$
of the $\hat{K}$ empirical estimates (once more accounting for edge-effects in different ways), especially at short distances, once more suggesting deviation from CSR.
Note that for geometrical reasons the
$\hat{K}$ empirical estimates will exhibit substantially greater variation at large distances;
it is therefore appropriate to confine attention to the left-hand third of the $x$-axis.
The excess over the theoretical $K_\text{pois}$ at short distances, particularly for the estimate $\hat{K}_\text{iso}$, indicates that defects are more clustered than would be expected from CSR.


![](image/fig4-a.pdf){ width=50% } ![](image/fig4-b.pdf){ width=50% }
\begin{figure}[!h]
\caption{Pilatus detector screen: (a) $K$ plot resulting from point pattern of defects. (b) $K$ plot resulting from point pattern of defects, corrected for inhomogeneity.
For geometrical reasons it is appropriate to focus attention on small distances.
There is more variation between different edge-corrections of empirical curves
than for $F$ and $G$ curves. Empirical curves are closer at short distances (200 pixels or less) to the theoretical curve based on CSR (left panel, blue dotted curve) but still exhibit some discrepancy hinting at possibly
greater clustering relative to CSR.
However all curves agree closely for short distances in the right panel,
in which a correction has been made
for inhomogeneity (which has already been noted when considering the intensity map).
This suggests that an inhomogeneous Poisson process provides a good fit for the data.
\label{fig:figure4}}
\end{figure}


* Plots are also available which take account of inhomogeneity and compare these estimates to theoretical functions
computed for inhomogeneous Poisson point processes:
Figure \ref{fig:figure4}(b) gives an example of this in the case of the $K$ function. The plots of the $\hat{K}_\text{inhom}$
empirical inhomogeneity-adjusted estimates agree much more closely with
the theoretical $K^\text{pois}_\text{inhom}$ function
at short distances, supporting the hypothesis that the pattern of defects is what
might be expected to arise from an _inhomogeneous_ Poisson process of defects.



D. Finally the relationship of the defect points to sub-panel boundaries can be studied by means of various logistic regression options, which assess whether damage intensity appears to depend on distance from the centre of the image or horizontal or vertical distance from sub-panel edges. When this data set is modelled in terms of Euclidean distance from the centre, the web application reports substantial evidence for positive dependence of defect intensity on distance from the centre (see the highly significant coefficient for `as.vector(dist)` in the following web application output),
conforming with the visual impression given by Figure \ref{fig:figure2}(a),
and further explaining the nature of the spatial inhomogeneity indicated
by Figure \ref{fig:figure4}. In fact this positive dependence reflects
manufacturing details of
this particular screen design: Diamond reports that Pilatus detector screen panels are tested before installation, and better panels are placed in the centre of the structured display.

Here follows output from the web application
after performing the above logistic regression:

--------------------------------

```
Call:
glm(formula = as.vector(pix_matrix) ~ as.vector(dist),
family = binomial(link = logit))

Deviance Residuals:
    Min       1Q   Median       3Q      Max
-0.0262  -0.0211  -0.0193  -0.0171   4.3297  

Coefficients:
                  Estimate Std. Error z value Pr(>|z|)    
(Intercept)     -9.392e+00  9.628e-02 -97.545   <2e-16 ***
as.vector(dist)  8.021e-04  8.726e-05   9.192   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 22261  on 6224000  degrees of freedom
Residual deviance: 22173  on 6223999  degrees of freedom
AIC: 22177

Number of Fisher Scoring iterations: 11
```
--------------------------------



The web application also provides for further graphical options,
such as the study of direction from a typical defect point to
its nearest neighbour within the relevant sub-panel,
analysis at the level of "events" (appropriately defined grouping of clumps of defect pixels) rather than individual defect points,
and exclusion of regions of the image rectangle for which the defect intensity is clearly different
(this often arises in XCT, where corners of the image exhibit high defect intensity, presumably deriving from mechanical
stress due to supports of the screen).




An extended example of use of the R package, paralleled by corresponding use of the web application,
is available as a vignette in the Github repository [_DetectorChecker_](https://github.com/alan-turing-institute/DetectorChecker).

The R package and web application together offer significant
opportunities to address interesting and important challenges for the data analysis of defective pixel patterns.
The web application offers the possibility of uploading users' data to
a data repository, thus permitting the possibility of organizing cooperative
statistical investigations comparing patterns across different machines and
different modes of usage. In particular we envisage its use to collect
time sequences of images, to permit statistical investigation by the Warwick team
of deterioration over time, using latent Markov models
of the life and death of defective pixels which are currently being developed.
Such analysis requires sustained and regular monitoring of a diversity
of screens from various devices, together with recording of relevant metadata
such as detector usage.
Interested users are encouraged to make contact to discuss these possibilities,
which will permit evidence-based analysis
to support decisions on refurbishment and/or replacement
strategies.

# Acknowledgements

We gratefully acknowledge support from the UK EPSRC (grant EP/K031066/1)
and the Alan Turing Institute (under the EPSRC
grant EP/N510129/1) during this project.

We also wish to thank Martin O'Reilly (The Alan Turing Institute), Nicola Tartoni and Ian Horswell (Diamond Lightsource, UK) for guidance on detector types and sample data sets, and Tristan Lowe (Henry Moseley X-ray Imaging Facility, University of Manchester)
and Martin Turner (ITS Research IT, University of Manchester) for discussions and feedback.

# References
---
output:
  html_document:
    fig_caption: yes
  pdf_document:
    fig_caption: yes
title: "DetectorChecker Usage"
author: "Julia Brettschneider, Oscar Giles, Tomas Lazauskas, Wilfrid Kendall"
date: "`r Sys.Date()`"
vignette: >
  %\VignetteIndexEntry{DetectorChecker Usage}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}

references:
- id: Chiu2013
  title: Stochastic Geometry and its Applications
  author:
  - family: Chiu
    given: Sung Nok
  - family: Stoyan
    given: Dietrich
  - family: Kendall
    given: Wilfrid S.
  - family: Mecke
    given: Joseph 
  container-title: In Wiley Series in Probability and Statistics. 
  DOI: 10.1002/9781118658222
  publisher: Wiley
  type: book
  issued:
    year: 2013
---

---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# DetectorChecker

The DetectorChecker `R` package offers a quality check of digital X-ray detectors based on spatial analysis of damaged pixel distributions. It can also be used for analysing other types of defects on other types of instruments, as long as the data is similar. For example, instead of damaged pixels, the spatial arrangement of other types of dysfunctional pixels could be analysed. For an example of a different
type of instrument, consider CMOS sensors as used in digital cameras. 

Here is how to load the package. 
```{r}
library(detectorchecker)
```

## Loading and visualising detector layout and data

The elementary units (e.g. pixels) have to be arranged in a rectangular grid. The grid may be divided into rectangular modules, which may be separated by gaps. There are technical reasons for such gaps, and they are excluded from the analysis. 

### Check available detectors

Several predefined detector layouts are already included in the package. `Pilatus`, `Excalibur` and `PerkinElmerFull` (2000 by 2000 pixels) are common types of detectors built into professional X-ray machines. `PerkinElmerCropped1600` and `PerkinElmerRefurbished` are modified version of `PerkinElmerFull`; they were created by restricting the area. 
A list with all available detector layouts can be displayed.

```{r}
available_detectors
```

If your data has been created using one of the predefined layouts
then it suffices to upload the relevant detector layout.
Otherwise a user-defined layout can be specified by reading in an
appropriate text file as explained further below.

### Create detector object from available layouts

The first step is to create an `R` object with the layout matching your data set. 
Here are example invocations for `Pilatus` and `PerkinElmer_Full`.

```{r}
detector_pilatus <- create_detector("Pilatus")
detector_perkinfull <- create_detector("PerkinElmerFull")
```

A summary of the characteristics of the detector layout can be printed.

```{r}
summary(detector_pilatus)
```

You can check if an object is an instance of the `detector` with

```{r}
is.detector(detector_pilatus)
```

We will now explain the attributes of this object using `Pilatus` and `PerkinElmerFull` as examples. 
The first four attributes contain the name, a date (optional; use `NA` if there is no need to specify this),
the width (number of pixels horizontally) and the height (numbers of pixels vertically) of the detector. 

```{r}
paste("name:", detector_pilatus$name)
paste("date:", detector_pilatus$date)
paste("width:", detector_pilatus$detector_width)
paste("height", detector_pilatus$detector_height)
```

The grid layout is given by the following attributes The number of modules per row is given by the module_row_n attribute and per column in module_col_n.

```{r}
paste("module_col_n:", detector_pilatus$module_col_n)
paste("module_row_n:", detector_pilatus$module_row_n)
```

The sizes of the modules are stored in the next two attributes. `module_col_sizes` contains the widths (numbers of pixel columns in each module), while `module_row_sizes` contains the heights (numbers of pixel rows in each module). In the `Pilatus` detector all modules are of the same dimensions. In the `PerkinElmerFull` the first and last module of each row have fewer pixels.

```{r}
print("Pilatus: Cols and Rows")
detector_pilatus$module_col_sizes
detector_pilatus$module_row_sizes

print("Perkinfull: Cols and Rows")
detector_perkinfull$module_col_sizes
detector_perkinfull$module_row_size
```

The next two attributes contain the sizes of the gaps. The numbers of pixel columns per gap are in `gap_col_sizes`, while the numbers of pixels rows per gap are in `gap_row_sizes`
In the `Pilatus` detectors all horizontal gaps are of the same size, and all vertical gaps are of the same size. The various Perkin Elmer detectors and the `Excalibur` detector have no gaps at all. However, it is possible to create layouts with more irregular gap sizes (see the comments below on user-defined layouts).

```{r}
print("Pilatus: Gap Cols and Rows")
detector_pilatus$gap_col_sizes
detector_pilatus$gap_row_sizes

print("Perkinfull: Gap Cols and Rows")
detector_perkinfull$gap_col_sizes
detector_perkinfull$gap_row_sizes
```

The final attributes are calculated from the other slots; they contain the left and the right edge (column numbers) and the top and bottom edge (row numbers) of each module. If there are no gaps, such as in the three Perkin Elmer layouts, then the left edge of a module will be adjacent to the right edge of the previous one in the same row of the grid and the bottom edge of a module will be adjacent to the top edge of the previous one.

```{r}
print("Pilatus: Gap Cols and Rows")
detector_pilatus$module_edges_col
detector_pilatus$module_edges_row

print("Perkinful: Gap Cols and Rows")
detector_perkinfull$module_edges_col
detector_perkinfull$module_edges_row
```

The final attribute `detector_inconsistency` is a consistency check. If the sizes of the modules add up to the correct numbers, the status is 0. Otherwise, it is 1 and a suitable error message will be produced. The remaining slots will be filled when the detector damage information is read in.

```{r}
detector_pilatus$detector_inconsistency
```


### Create detector object from user-defined layouts

If your detector's layout is not contained in the list then you can create and read in the relevant layout parameters. You need to create a textfile that contains this information in a simple list setting all relevant parameters and store it in a directory of your choice.

We illustrate this by two examples. The first example is a simple layout for a detector with photographic aspect ratio with no subdivisions. We stored the lines below as a plain text file under the name  `layout_par_photographic.txt` in the subdirectory 
`extdata/user-defined/photographic/` of the `detectorchecker`package.

```{r}
detector_width = 600
detector_height = 400
module_col_n = 1
module_row_n = 1
module_col_sizes = c(600)
module_row_sizes = c(400)
gap_col_sizes = c()
gap_row_sizes = c()
module_edges_col = NA
module_edges_row = NA
```

The second example is an irregular layout with varying module dimensions and varying gap sizes. We stored the lines below as a plain text file under the name `layout_par_irregular.txt` in the subdirectory 
`extdata/user-defined/irregular/` of the `detectorchecker`package.

```{r}
detector_width = 1720
detector_height = 1060
module_col_n = 7
module_row_n = 5
module_col_sizes = c(100, 200, 300, 400, 300, 200, 100)
module_row_sizes = c(100, 200, 400, 200, 100)
gap_col_sizes = c(10,20,30,30,20,10)
gap_row_sizes = c(10,20,20,10)
module_edges_col = NA
module_edges_row = NA
```

To generate a detector object from user-defined parameters, we use a function reading the parameters from the text file as shown below for our examples. For your own user-defined detector you need to set the variable `file_path` to point to the directory where you stored the parameter file.

```{r}
file_path <-  system.file("extdata", "user-defined", 
                          "photographic", "layout_par_photographic.txt", 
                          package ="detectorchecker")
user_def_detector_photographic <- readin_detector(file_path)
file_path <-  system.file("extdata", "user-defined", 
                          "irregular", "layout_par_irregular.txt", 
                          package ="detectorchecker")
user_def_detector_irregular <- readin_detector(file_path)
```

### Visualise detector object

A visualisation of the layout shows modules and gaps as shown below for two of the pre-specified layouts. The default is to plot a title, but we can override this by setting the parameter caption to `FALSE`. The same rule applies for all plots produced using this package.

```{r}
plot(detector_pilatus)
plot(detector_pilatus, caption=FALSE)
plot(detector_perkinfull)
```

This also works for the two user-defined detectors.
```{r}
plot(user_def_detector_photographic)
plot(user_def_detector_irregular)
```


### Detector geometry

The class for detectors considered admissible for this package is illustrated by Figure \ref{fig:detector_layout}

![General layout with parameters\label{fig:detector_layout}](layoutDrawingPar.png)

### Upload files containing pixel damage locations

The next step is to upload a file containing information about damaged pixel locations. The format depends on what is provided for the detector type. 

#### `TIF` format

For `Pilatus`, the damaged pixels locations are usually provided as bad pixel masks in the form of `tif` files. A sample file with raw damaged pixel mask data is provided with the package and can be loaded by specifying the relevant folder location (here `extdata`) and file name (here `Pilatus/badpixel_mask.tif`). 

```{r}
detector_pilatus <- create_detector("Pilatus")  
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif", 
                          package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
  detector = detector_pilatus, file_path = file_path)
```

The spatial locations of the damaged pixels can be visualised.

```{r}
plot_pixels(detector = detector_pilatus)
```

Individual modules can be plotted separately. 
Here we plot the module located in the 4^th^ column and 5^th^ row of the module layout grid.

```{r}
plot_pixels(detector = detector_pilatus, 
            col = 4, row = 5)
```

For the user-defined layout used above we can also display damaged pixels. If the location are stored in the form of a `tif` file then this looks as follows.

```{r}
file_path <-  system.file("extdata", "user-defined", 
                          "photographic", "photographic.tif", 
                          package ="detectorchecker")
user_def_detector_photographic <- load_pix_matrix(
  detector = user_def_detector_photographic, file_path = file_path)
plot_pixels(detector = user_def_detector_photographic, caption = TRUE)
```

#### `XML` format

For the three Perkin Elmer layouts (`PerkinElmerFull`, `PerkinCropped1600`, and `PerkinElmerRefurbished`) the damaged pixels locations are usually stored as `xml` files. Sample files for each of the three Perkin Elmer layouts are supplied with the package and can be loaded by specifying the location and file name. 

```{r}
detector_perkinfull <- create_detector("PerkinElmerFull")
file_path <-  system.file("extdata", "PerkinElmerFull", 
                          "BadPixelMap_t1.bpm.xml",
                          package = "detectorchecker")
detector_perkinfull <- load_pix_matrix(
  detector = detector_perkinfull, file_path = file_path)
plot_pixels(detector = detector_perkinfull)
plot_pixels(detector = detector_perkinfull, 
            col = 11, row = 1)
```

To plot another of the `Perkin Elmer` layouts the detector type needs to be changed. For the refurbished detector layout this is done as follows:
To plot another of the Perkin Elmer layouts the detector type needs to be changed. For the refurbished detector layout this is done as follows:

```{r}
detector_perkinrefurb <- create_detector("PerkinElmerRefurbished")
file_path <-  system.file("extdata", "PerkinElmerRefurbished", 
                          "BadPixelMap_t1.bpm.xml", 
                          package = "detectorchecker")
detector_perkinrefurb <- load_pix_matrix(
  detector = detector_perkinrefurb, file_path = file_path)
plot_pixels(detector = detector_perkinrefurb)
```

To see the other sample data sets, change the file name, e.g.:

```{r}
detector_perkinrefurb <- create_detector("PerkinElmerRefurbished")
file_path <-  system.file("extdata", "PerkinElmerRefurbished", 
                          "BadPixelMap_t2.bpm.xml", 
                          package = "detectorchecker")
detector_perkinrefurb <- load_pix_matrix(
  detector = detector_perkinrefurb, file_path = file_path)
plot_pixels(detector = detector_perkinrefurb)
```

This also works for a user-defined detector if the damaged pixel locations are stored in `xml` format. (Dead pixels locations in the example file were generated by simulating a Poisson point pattern with inhomogeneous density.)

```{r}
file_path <-  system.file("extdata", "user-defined", "photographic",
                          "badpixelmap_photographic_inhom.xml", 
                          package ="detectorchecker")
user_def_detector_photographic <- load_pix_matrix(
  detector = user_def_detector_photographic, file_path = file_path)
plot_pixels(detector = user_def_detector_photographic, 
            caption = TRUE)
```

For the second user-defined example, we proceed similarly.

```{r}
file_path <-  system.file("extdata", "user-defined", 
                          "irregular", "badpixelmap_irregular.xml", 
                          package ="detectorchecker")
user_def_detector_irregular <- load_pix_matrix(
  detector = user_def_detector_irregular, file_path = file_path)
plot_pixels(detector = user_def_detector_irregular)
```
#### `HDF` format

For `Excalibur`, the damaged pixels locations are usually provided in a combination of 6 `hdf` files corresponding to the layout rows. One such collection of sample files is provided with the package.

```{r}
detector_exc <- create_detector("Excalibur")

file_path <-  c(
  system.file("extdata", "Excalibur", "pixelmask.fem1.hdf", 
              package ="detectorchecker"),
  system.file("extdata", "Excalibur", "pixelmask.fem2.hdf", 
              package ="detectorchecker"),
  system.file("extdata", "Excalibur", "pixelmask.fem3.hdf", 
              package ="detectorchecker"),
  system.file("extdata", "Excalibur", "pixelmask.fem4.hdf", 
              package ="detectorchecker"),
  system.file("extdata", "Excalibur", "pixelmask.fem5.hdf", 
              package ="detectorchecker"),
  system.file("extdata", "Excalibur", "pixelmask.fem6.hdf", 
              package ="detectorchecker")
  )
detector_exc <- load_pix_matrix(
  detector = detector_exc, file_path = file_path)
plot_pixels(detector = detector_exc)
```

As for the other layout types, individual modules can be plotted.

```{r}
plot_pixels(detector = detector_exc,
            col = 6, row = 4)
```


## Analysis of damaged pixels

### Density

Approximate damage intensity can be visualized by plotting a density. The parameter `adjust` must be positive and can be varied to tune the smoothness of the density estimation procedure.

```{r}
plot_pixels_density(detector = detector_pilatus)
plot_pixels_density(detector = detector_perkinfull)
plot_pixels_density(detector = user_def_detector_photographic)
```

The same plot can be carried out for individual modules. In the case of a density map, the smoothing chosen for the sub-panel is different from that chosen for the plot for the entire panel, so it will not simply be a detail of the larger image.

```{r}
plot_pixels_density(detector = detector_pilatus, 
                    row = 5, col = 1)
```

The parameter `adjust` influences the bandwidth used in the kernel based density estimation. The default is set at 0.5. Smaller values will emphasise the role of individual pixels, while larger values result in smoothing over larger areas.

```{r}
plot_pixels_density(detector = detector_pilatus, adjust=0.1)
plot_pixels_density(detector = detector_pilatus, adjust=2)
```


### Numerical summaries

A calculation of the total number of damaged pixels, the average number of damaged pixels per module, and a Chi-square test for its independence of the module can be performed. We show this here for `Pilatus`, preceded by a summary of the relevant layout characteristics.

```{r}
summary(detector_pilatus)
cat(dead_stats_summary(detector_pilatus))
```

A systematic display of the counts of damaged pixels per module can be displayed for the whole detector or for an individual module.

```{r}
detector_pilatus <- get_dead_stats(detector_pilatus)
plot_pixels_count(detector = detector_pilatus)
plot_pixels_count(detector = detector_pilatus, 
                  row = 5, col = 12)
```


### Arrows to nearest neighbours

For each damaged pixel, a nearest neighbour can be determined and all such relationships can be visualised by plotting directed arrows. Specifying coordinates in the layout grid plots only the corresponding modules. We illustrate this for one of the Perkin Elmer data examples.

```{r}
plot_pixels_arrows(detector = detector_perkinfull)
plot_pixels_arrows(detector = detector_perkinfull, 
                   row = 1, col = 1)
```


### Angles between nearest neighbours

Angles between neighbours can be summarised in a rose diagram.

```{r}
plot_pixels_angles(detector = detector_pilatus)
```

In some modules, the rose diagram is dominated by vertical and horizontal directions. However, in other modules this is not the case. 

```{r}
plot_pixels_angles(detector = detector_pilatus, 
                   row = 5, col = 11)
plot_pixels_angles(detector = detector_pilatus, 
                   row = 5, col = 4)
```

## Analysis of complete spatial randomness

A number of distance-based functions can be used to check whether
the point pattern created by the damaged pixels has 
features that would be expected from a pattern exhibiting complete spatial randomness (CSR) [@Chiu2013]: Ripley's $K$ function, the empty space function $F$ and the nearest-neighbour function $G$. Scales of the axes are chosen by `spatstat` default.


### F Function
F function or “empty space function” computes the distribution of the nearest distance to a defect
point from a typical location chosen uniformly from the image rectangle; if the point pattern did in fact
satisfy CSR then an empirical estimate of the F function could be viewed as a random perturbation
of the theoretical F function under CSR.

### G Function
The G function computes the distribution of nearest-neighbour distances between defect points; if the
point pattern did in fact satisfy CSR then an empirical estimate of the G function could be viewed as a
random perturbation of the theoretical G function under CSR.

### The K Function
The K function (Ripley’s K function) computes the mean number of defect points within a distance r
of a typical defect point, viewed as a function of r; if the point pattern did in fact satisfy CSR then an
empirical estimate of the K function could be viewed as a random perturbation of the theoretical K
function under CSR.

### $K$, $F$, $G$ function plots

The basic versions of these functions assume a homogeneous background density. 

```{r}
plot_pixels_kfg(detector = detector_pilatus, func = "K")
plot_pixels_kfg(detector = detector_pilatus, func = "F")
plot_pixels_kfg(detector = detector_pilatus, func = "G")
```

Individual modules can be selected as before. The functions consistently show a deviation from what would be expected from complete spatial randomness (which is graphed as `K_pois`). Note that there are several different empirical variants of each function, corresponding to different edge-correction techniques.

```{r, fig.show='hide'}
plot_pixels_kfg(detector = detector_pilatus, func = "K", 
                row = 5, col = 12)
plot_pixels_kfg(detector = detector_pilatus, func = "F", 
                row = 5, col = 12)
plot_pixels_kfg(detector = detector_pilatus, func = "G", 
                row = 5, col = 12)
```


### Inhomogeneous K, F, G plots

There are also versions of these functions that allow for an inhomogeneous background density (itself estimated from the data).

```{r}
plot_pixels_kfg(detector = detector_pilatus, func = "Kinhom")
plot_pixels_kfg(detector = detector_pilatus, func = "Finhom")
plot_pixels_kfg(detector = detector_pilatus, func = "Ginhom")
```

As is the case for their homogeneous counterparts, they can also be computed for individual modules.

```{r, fig.show='hide'}
plot_pixels_kfg(detector = detector_pilatus, func = "Kinhom", 
                row = 5, col = 12)
plot_pixels_kfg(detector = detector_pilatus, func = "Finhom", 
                row = 5, col = 12)
plot_pixels_kfg(detector = detector_pilatus, func = "Ginhom", 
                row = 5, col = 12)
```


## Analysis of damage events

Instead of interpreting every damaged pixel by itself, small clusters and lines can be interpreted as one damage event. Create a new object by specifying which kinds of small ensembles should be summarised as damage events (1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines, 6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines). To include all kinds of ensembles list them all. 

The differences between pixel and event level are demonstrated well e.g. in the `PerkinElmerFull` example data sets, because they possess lines of damaged pixels. 

```{r}
incl_event_list <- list(1,2,3,4,5,6,7,8)
detector_perkinfull_events <- detectorchecker::find_clumps(detector_perkinfull)
```

For the Perkin Elmer damaged pixel example files, there is an obvious difference between the pixel level (first plot below) and the event level (second plot below).

```{r}
plot_pixels(detector_perkinfull)
detectorchecker::plot_events(detector_perkinfull_events, 
                             incl_event_list = incl_event_list)
```

An individual module can be plotted by referring to the corresponding row and column.

```{r}
plot_pixels(detector_perkinfull, 
            row = 1, col = 11)
plot_events(detector_perkinfull_events, 
            incl_event_list = incl_event_list, 
            row = 1, col = 11)
```

Much of the analysis carried out on the pixel level can also be conducted on the event level.

### Numerical summaries

A graphical display of the counts of damage events per module can be displayed for the whole detector or a module.

```{r}
detectorchecker::plot_events_count(detector_perkinfull_events, 
                                   incl_event_list = incl_event_list)
```

### Density

For density plots the optional adjust parameter can be omitted (default=0.5) or changed to another positive value. Individual modules can be selected.

```{r}
detectorchecker::plot_events_density(detector_perkinfull_events, 
                                     incl_event_list = incl_event_list)
detectorchecker::plot_events_density(detector_perkinfull_events, 
                                     incl_event_list = incl_event_list, 
                                     row=1, col=16, adjust=1)
```

### Arrows between nearest neighbours 

For the complete detector or for individual modules:

```{r}
detectorchecker::plot_events_arrows(detector_perkinfull_events, 
                                    incl_event_list = incl_event_list)
detectorchecker::plot_events_arrows(detector_perkinfull_events, 
                                    incl_event_list = incl_event_list, row=1, col=16)
```

### Angles between nearest neighbours 

For the complete detector or for individual modules:

```{r}
detectorchecker::plot_events_angles(detector_perkinfull_events, 
                                    incl_event_list = incl_event_list)
detectorchecker::plot_events_angles(detector_perkinfull_events, 
                                    incl_event_list = incl_event_list, 
                                    row=1, col=16)
```

## Analysis of complete spatial randomness

This can be conducted in the same way as for the pixel level. It is worth comparing the results on the pixel level with those on the event level, because results may be different. We show below a number of ways to do this for the $K$ function. Similar approaches work for $F$ and $G$ functions.

```{r}
detectorchecker::plot_events_kfg(detector_perkinfull_events, func = "K",
                                 incl_event_list = incl_event_list)
detectorchecker::plot_events_kfg(detector_perkinfull_events, func = "K",
                                 incl_event_list = incl_event_list, row=1, col=16)
detectorchecker::plot_events_kfg(detector_perkinfull_events, func = "Kinhom", 
                                 incl_event_list = incl_event_list)
```

## Removing areas of high density pixel damage 

In some situations, the analysis may be dominated by an area of elevated damage. The investigation of complete spatial randomness then becomes uninformative. The area with elevated damage can be removed. 

```{r}
detector_perkinfull_modified <- remove_high_density_cluster(
  detector_perkinfull, min_pts = 30, eps_adjust = 0.05)
plot_pixels(detector = detector_perkinfull_modified)
``` 

The analysis for complete spatial randomness can be conducted on the remaining area.  
```{r}
plot_pixels_kfg(detector = detector_perkinfull_modified, func = "K")
plot_pixels_kfg(detector = detector_perkinfull_modified, func = "F")
plot_pixels_kfg(detector = detector_perkinfull_modified, func = "G")
```


## Models

This section concerns inferential statistical methods to analyse the state (damaged or intact) of the pixels. General linear models can include a spatial covariate reflecting the location of a pixel relative to the centre, to the edges or the corners of the detector. 
The spatial covariate can also express the location relative to its containing module (sub-panel).
This can be distance of the pixel to the closest vertical or horizontal edge of the containing module, or the minimum of both of these (i.e. the distance to the closest edge of the containg module).
Running times can be slow, but should not take longer than a minute
for detector types included in this package
when using a notebook computer of normal specification (as of 2020). 

We start with the spatial covariates that reflect the geometry of the pixels with respect to the detector as a whole. 
The following command fits a model that includes the euclidean distance of a pixel to the centre of the detector as covariate.
```{r}
glm_pixel_ctr_eucl(detector_pilatus)
```

The next command includes the maximum of the vertical and horizontal distances of a pixel to the centre of the detector (also called $L^\infty$ or L-infinity distance). 
```{r}
glm_pixel_ctr_linf(detector_pilatus)
```

It is also possible to include the distances of a pixel to the nearest corner of the detector as a covariate.

```{r}
glm_pixel_dist_corner(detector_pilatus)
```

Now we look at covariates that reflect the position of the pixel with respect to the containing module.
The function `glm_pixel_dist_edge_col()` fits a model that includes the distance of the relevant pixel
from the nearest vertical edge of the containing module.

```{r}
glm_pixel_dist_edge_col(detector_pilatus)
```

The function `glm_pixel_dist_edge_row()` does the same for rows, and `glm_pixel_dist_edge_min()` uses the smaller distance of the two. 

All these commands are also available on the event level by using the same name but 
substituting `events` for `pixel` in the command name.  
For example, the following command fits a model that includes
as covariate the euclidean distance of an event to the centre.

```{r}
glm_events_ctr_eucl(detector_pilatus)
```

For illustrative purposes it is also possible to visualise the spatial covariates. 
(This step would typically be missed out, since the plotting process is quite slow.)

The function `plot_pixel_ctr_eucl()` displays euclidean distance of pixels
to the centre of the detector using a colour scheme. 
The function `plot_pixel_ctr_linf()` displays the maximum of the vertical and horizontal distances of a pixel to the centre of the detector (also called $L^\infty$ or L-infinity distance). 
The next command displays euclidean distance of pixels to the nearest corner.

```{r}
plot_pixel_dist_corner(detector_pilatus)
```

The function `plot_pixel_dist_edge_col()` displays distance of pixels from the nearest vertical module edge, while `plot_pixel_dist_edge_row()` does the same for horizontal boundaries. 
In conclusion, the following command is based on 
the minimum of the distances displayed by `plot_pixel_dist_edge_col()`
and `plot_pixel_dist_edge_row()`;
in other words, it displays the distance of pixels to nearest module edge.

```{r}
plot_pixel_dist_edge(detector_pilatus)
```



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{Excalibur_Detector}
\alias{Excalibur_Detector}
\title{A S3 class to represent the Excalibur detector.}
\usage{
Excalibur_Detector()
}
\value{
Excalibur detector object
}
\description{
A S3 class to represent the Excalibur detector.
}
\examples{
detc <- Excalibur_Detector()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\docType{data}
\name{available_detectors}
\alias{available_detectors}
\title{A list of available preconfigured detectors. These can be created by \code{create_detector}}
\format{
An object of class \code{character} of length 5.
}
\usage{
available_detectors
}
\description{
A list of available preconfigured detectors. These can be created by \code{create_detector}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{plot_pixel_dist_corner}
\alias{plot_pixel_dist_corner}
\title{Calculates and plots pixel distances from corners}
\usage{
plot_pixel_dist_corner(detector, file_path = NA)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}
}
\description{
Calculates and plots pixel distances from corners
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_pixel_dist_corner(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{.detector_consist_check}
\alias{.detector_consist_check}
\title{Basic checks if parameters entered (slightly redundant on purpose) add up}
\usage{
.detector_consist_check(detector = NA)
}
\arguments{
\item{detector}{Detector object}
}
\description{
Basic checks if parameters entered (slightly redundant on purpose) add up
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{plot_pixel_dist_edge}
\alias{plot_pixel_dist_edge}
\title{Calculates and plots L-infinity distances from the module edges}
\usage{
plot_pixel_dist_edge(detector, file_path = NA)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}
}
\description{
Calculates and plots L-infinity distances from the module edges
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_pixel_dist_edge(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{.Default_Detector}
\alias{.Default_Detector}
\title{Detector module
A S3 class to represent a detector.}
\usage{
.Default_Detector(
  name = "Default",
  date = NA,
  detector_width = NA,
  detector_height = NA,
  module_col_n = NA,
  module_row_n = NA,
  module_col_sizes = NA,
  module_row_sizes = NA,
  gap_col_sizes = NA,
  gap_row_sizes = NA,
  module_edges_col = NA,
  module_edges_row = NA,
  detector_inconsistency = NA,
  pix_matrix = NA,
  pix_dead = NA,
  dead_stats = NA,
  pix_dead_modules = NA,
  clumps = NA,
  clumps_col = NA,
  clumps_row = NA
)
}
\arguments{
\item{name}{detector's name}

\item{date}{date}

\item{detector_width}{detector's width}

\item{detector_height}{detector's height}

\item{module_col_n}{number of columns in the grid of modules}

\item{module_row_n}{number of rows in the grid of modules}

\item{module_col_sizes}{vector with widths of the modules}

\item{module_row_sizes}{vector with heights of the modules}

\item{gap_col_sizes}{vector with widths of the gaps}

\item{gap_row_sizes}{vector with heights of the gaps}

\item{module_edges_col}{vector of columns that contain edges of modules}

\item{module_edges_row}{vector of rows that contain edges of modules}

\item{detector_inconsistency}{counts inconsistencies found in parameters entered}

\item{pix_matrix}{pixel matrix}

\item{pix_dead}{dead pixels coordinates}

\item{dead_stats}{dead pixel statistics}

\item{pix_dead_modules}{assigned module for each dead pixel}

\item{clumps}{clumps data (xyc_df data frame with pixels and their clump ID's, xyc_events data frame with clusters (clumps) and their clump ID's and centre coordinates)}

\item{clumps_col}{col number of the module on which analysis was performed}

\item{clumps_row}{row number of the module on which analysis was performed}
}
\value{
Detector object
}
\description{
Detector module
A S3 class to represent a detector.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{plot_pixels_density}
\alias{plot_pixels_density}
\title{A function to plot densities of dead pixels of detector or module}
\usage{
plot_pixels_density(
  detector,
  file_path = NA,
  adjust = 0.5,
  row = NA,
  col = NA,
  caption = TRUE,
  color = topo.colors(50)
)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}

\item{adjust}{Kernel density bandwidth}

\item{row}{Module row number}

\item{col}{Module column number}

\item{caption}{Flag to turn on/off figure caption}

\item{color}{a list of colors}
}
\description{
A function to plot densities of dead pixels of detector or module
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_pixels_density(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{plot_pixel_ctr_eucl}
\alias{plot_pixel_ctr_eucl}
\title{Calculates and plots pixel euclidean distance from the centre}
\usage{
plot_pixel_ctr_eucl(detector, file_path = NA)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}
}
\description{
Calculates and plots pixel euclidean distance from the centre
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_pixel_ctr_eucl(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{pixel_dist_ctr_eucl}
\alias{pixel_dist_ctr_eucl}
\title{Calculate euclidean distance from the center of a module for each pixel}
\usage{
pixel_dist_ctr_eucl(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Matrix of euclidean distances
}
\description{
Calculate euclidean distance from the center of a module for each pixel
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
pixel_dist_ctr_eucl(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{.matrix_from_xml}
\alias{.matrix_from_xml}
\title{Reads in xml file and returns a pixel matrix}
\usage{
.matrix_from_xml(detector, file_path)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Path to the xml file}
}
\value{
Data from an xml file
}
\description{
Reads in xml file and returns a pixel matrix
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{.xyc_pixels2events}
\alias{.xyc_pixels2events}
\title{Modifying clusters to events (consisting of 1 pixel representing the cluster)
Make into a point pattern of just events rather than pixels. Using xyc_ply object.
Collapse in one point using centres for clusters, but end points for lines, type dependend:
type 5 (closest to upper edge): ymin
type 6 (closest to lower edge): ymax
type 7 (closest to right edge): xmin
type 8 (closest to left edge):  xmax
This is inspired by Perkin Elmer Detector and be replaced by other choices if desired.}
\usage{
.xyc_pixels2events(xyc_ply)
}
\arguments{
\item{xyc_ply}{clums data frame}
}
\value{
events
}
\description{
Modifying clusters to events (consisting of 1 pixel representing the cluster)
Make into a point pattern of just events rather than pixels. Using xyc_ply object.
Collapse in one point using centres for clusters, but end points for lines, type dependend:
type 5 (closest to upper edge): ymin
type 6 (closest to lower edge): ymax
type 7 (closest to right edge): xmin
type 8 (closest to left edge):  xmax
This is inspired by Perkin Elmer Detector and be replaced by other choices if desired.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{plot_pixels_angles}
\alias{plot_pixels_angles}
\title{Plot nearest neighbour angles of dead pixels of detector or module}
\usage{
plot_pixels_angles(
  detector,
  file_path = NA,
  row = NA,
  col = NA,
  caption = TRUE
)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}

\item{row}{Module row number}

\item{col}{Module column number}

\item{caption}{Flag to turn on/off figure caption}
}
\description{
Plot nearest neighbour angles of dead pixels of detector or module
}
\examples{
detector_perkinfull <- create_detector("PerkinElmerFull")
file_path <-  system.file("extdata", "PerkinElmerFull",
                          "BadPixelMap_t1.bpm.xml",
                          package = "detectorchecker")
detector_perkinfull <- load_pix_matrix(
detector = detector_perkinfull, file_path = file_path)
plot_pixels_angles(detector_perkinfull)
plot_pixels_angles(detector_perkinfull, row = 1, col = 1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{.plot_kfg}
\alias{.plot_kfg}
\title{Plots K, F, G functions}
\usage{
.plot_kfg(ppp_obj, func, file_path = NA, caption = TRUE)
}
\arguments{
\item{ppp_obj}{ppp object}

\item{func}{Function name}

\item{file_path}{Output file path}

\item{caption}{Flag to turn on/off figure caption}
}
\description{
Plots K, F, G functions
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{.perform_glm}
\alias{.perform_glm}
\title{A simple wrapper around \code{glm()} with family = binomial(link = logit)}
\usage{
.perform_glm(symb_expr, family = binomial(link = logit))
}
\arguments{
\item{symb_expr}{symbolic description of the linear predictor}

\item{family}{a description of the error distribution}
}
\value{
Fitted model

glm_git fitted model
}
\description{
Calls glm(formula = symb_expr, family = family)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{.get_clump_event_ppp}
\alias{.get_clump_event_ppp}
\title{Creates ppp object for damaged detector events}
\usage{
.get_clump_event_ppp(
  detector,
  incl_event_list = NA,
  height = NULL,
  width = NULL
)
}
\arguments{
\item{detector}{Detector object}

\item{incl_event_list}{a list of events to be included}

\item{height}{Detector height}

\item{width}{Detector width}
}
\value{
ppp object for damaged detector events
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{.detector_edges}
\alias{.detector_edges}
\title{Defines the coordinates of detector's edges using module and gap sizes
Function is in 1d context to be applied to rows and cols separately.
Edges are inside the modules (first/last row/col of module).}
\usage{
.detector_edges(m, g)
}
\arguments{
\item{m}{vector of module sizes}

\item{g}{vectors of gap sizes}
}
\value{
Matrix with the information about the edges
}
\description{
Defines the coordinates of detector's edges using module and gap sizes
Function is in 1d context to be applied to rows and cols separately.
Edges are inside the modules (first/last row/col of module).
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{.mask_to_events}
\alias{.mask_to_events}
\title{Converts mask (dead pixels) to events}
\usage{
.mask_to_events(detector, dead_pix_mask, row = NA, col = NA)
}
\arguments{
\item{detector}{Detector object}

\item{dead_pix_mask}{Dead pixels mask}

\item{row}{Module row number}

\item{col}{Module col number}
}
\value{
list of pixels and events
}
\description{
Converts mask (dead pixels) to events
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{check_detector_avail}
\alias{check_detector_avail}
\title{Checks whether \code{detector_name} is preconfigured.}
\usage{
check_detector_avail(detector_name)
}
\arguments{
\item{detector_name}{The name of the detector}
}
\value{
Boolean
}
\description{
If TRUE can be created by \code{create_detector}
}
\examples{
check_detector_avail('Pilatus')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{.create_ppp_gaps_col}
\alias{.create_ppp_gaps_col}
\title{Creates ppp object of horizontal gaps}
\usage{
.create_ppp_gaps_col(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Point pattern dataset
}
\description{
Creates ppp object of horizontal gaps
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{readin_detector}
\alias{readin_detector}
\title{Reads in a user defined detector from a file}
\usage{
readin_detector(file_path)
}
\arguments{
\item{file_path}{A path to the user defined detector file}
}
\value{
Detector object
}
\description{
Reads in a user defined detector from a file
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{dead_stats_summary}
\alias{dead_stats_summary}
\title{Summary of damaged pixels for a given detector}
\usage{
dead_stats_summary(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
A string with damaged pixels overall statistics
}
\description{
Compute summary statistics for a detector object.
Ensure a damaged pixel matrix has been loaded with \code{load_pix_matrix}
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Calculate dead_stats_summary
dead_stats_summary(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{plot_pixel_dist_edge_col}
\alias{plot_pixel_dist_edge_col}
\title{Calculates and plots horizontal distances from the module edges}
\usage{
plot_pixel_dist_edge_col(detector, file_path = NA)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}
}
\description{
Calculates and plots horizontal distances from the module edges
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_pixel_dist_edge_col(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module.R
\name{which_module_idx}
\alias{which_module_idx}
\title{Function returns both col and row number of a dead pixel.}
\usage{
which_module_idx(x, y, module_edges_col, module_edges_row)
}
\arguments{
\item{x}{pixel x coordinate}

\item{y}{pixel y coordinate}

\item{module_edges_col}{vector of columns that contain edges of modules}

\item{module_edges_row}{vector of rows that contain edges of modules}
}
\value{
tmp
}
\description{
Function returns both col and row number of a dead pixel.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{.create_ppp_edges_row}
\alias{.create_ppp_edges_row}
\title{This is the create_ppp_edges_row creation function}
\usage{
.create_ppp_edges_row(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Point pattern dataset
}
\description{
This is the create_ppp_edges_row creation function
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{.ini_graphics}
\alias{.ini_graphics}
\title{Starts the graphics device driver for producing graphics with respect to a
chosen format}
\usage{
.ini_graphics(file_path)
}
\arguments{
\item{file_path}{Output path with an extension}
}
\description{
Starts the graphics device driver for producing graphics with respect to a
chosen format
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{.get_events_mask}
\alias{.get_events_mask}
\title{Generates events mask (a matrix with pixels as 0 and events as 1)}
\usage{
.get_events_mask(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
events mask
}
\description{
Generates events mask (a matrix with pixels as 0 and events as 1) indicating if a pixel is in an event
as calculated by \code{find_clumps(detc)}
}
\examples{
detc <- Excalibur_exp_1
detc_with_clumps <- find_clumps(detc)
.get_events_mask(detc_with_clumps)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{get_ppp_dead}
\alias{get_ppp_dead}
\title{Generates point pattern dataset (ppp) for the dead pixels}
\usage{
get_ppp_dead(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
ppp of dead pixels
}
\description{
Uses \code{spatstat::ppp} internally.
Creates an object of class "ppp" representing a point pattern dataset in the two-dimensional plane.
See \href{https://www.rdocumentation.org/packages/spatstat/versions/1.63-3/topics/ppp}{spatstat} docs for details.
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Create a point pattern dataset from the detector
dead_ppp <- get_ppp_dead(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{.matrix_from_tiff}
\alias{.matrix_from_tiff}
\title{I/O module
Reads in tiff file and returns a pixel matrix}
\usage{
.matrix_from_tiff(detector, file_path)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Path to the tiff file}
}
\value{
Pixel matrix with dead pixels flagged with 1
}
\description{
I/O module
Reads in tiff file and returns a pixel matrix
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{dist_edge_min}
\alias{dist_edge_min}
\title{Calculate L-infinity distance to module edge}
\usage{
dist_edge_min(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
distance matrix
}
\description{
Calculate L-infinity distance to module edge
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Calculate L-infinity distance to module edge
dist_edge_min(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{.plot_arrows}
\alias{.plot_arrows}
\title{Plots nearest neighbor oriented arrows}
\usage{
.plot_arrows(ppp_obj, caption, file_path = NA)
}
\arguments{
\item{ppp_obj}{spatstat ppp object}

\item{caption}{caption of the figure}

\item{file_path}{file path}
}
\description{
Plots nearest neighbor oriented arrows
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{glm_pixel_ctr_linf}
\alias{glm_pixel_ctr_linf}
\title{Predict dead pixels from the pixel's parallel maxima}
\usage{
glm_pixel_ctr_linf(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Fitted model
}
\description{
Fit a logistic regression model using \code{glm}.
Predicts dead pixels from the pixel's parallel maxima
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Fit logistic regression model
glm_pixel_ctr_linf(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{get_dead_pix_mask}
\alias{get_dead_pix_mask}
\title{Creates a mask matrix of dead pixels}
\usage{
get_dead_pix_mask(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
dead pixel mask
}
\description{
Converts the pix_dead attribute of a detector (NX2 list) to a matrix of 1 and 0 (1 for dead pixel)
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Calculate dead pixel mask
get_dead_pix_mask(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module.R
\name{.dist_edge}
\alias{.dist_edge}
\title{Function returns distance of a pixel to module edges.}
\usage{
.dist_edge(xy, module_edges)
}
\arguments{
\item{xy}{Coordinate of pixel}

\item{module_edges}{vector of edges of a module}
}
\value{
tmp Distance to edges
}
\description{
Function returns distance of a pixel to module edges.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{.assign_pixel_matrix}
\alias{.assign_pixel_matrix}
\title{Assign dead pixels to a detector}
\usage{
.assign_pixel_matrix(detector, pix_matrix)
}
\arguments{
\item{detector}{Detector object}

\item{pix_matrix}{A pixel matrix}
}
\description{
Assign dead pixels to a detector
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{plot_events_kfg}
\alias{plot_events_kfg}
\title{Plots K, F, G functions of a detector or module}
\usage{
plot_events_kfg(
  detector,
  func,
  file_path = NA,
  row = NA,
  col = NA,
  caption = TRUE,
  incl_event_list = NA
)
}
\arguments{
\item{detector}{Detector object}

\item{func}{Function name}

\item{file_path}{Output file path}

\item{row}{Module row number}

\item{col}{Module column number}

\item{caption}{Flag to turn on/off figure caption}

\item{incl_event_list}{a list of events to be included}
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
\examples{
detector_perkinfull <- create_detector("PerkinElmerFull")
file_path <-  system.file("extdata", "PerkinElmerFull",
                          "BadPixelMap_t1.bpm.xml",
                          package = "detectorchecker")
detector_perkinfull <- load_pix_matrix(
detector = detector_perkinfull, file_path = file_path)
detector_perkinfull_events = find_clumps(detector_perkinfull)
plot_events_kfg(detector_perkinfull_events, "K")
plot_events_kfg(detector_perkinfull_events, "F")
plot_events_kfg(detector_perkinfull_events, "G")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{plot_pixels}
\alias{plot_pixels}
\title{A function to plot detector with damaged pixels overlayed}
\usage{
plot_pixels(detector, col = NA, row = NA, file_path = NA, caption = TRUE)
}
\arguments{
\item{detector}{Detector object}

\item{col}{Module column number}

\item{row}{Module row number}

\item{file_path}{Output file path}

\item{caption}{Flag to turn on/off figure caption}
}
\description{
A function to plot detector with damaged pixels overlayed
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_pixels(detector_pilatus)
plot_pixels(detector_pilatus, row = 1, col = 1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{glm_events_dist_edge_min}
\alias{glm_events_dist_edge_min}
\title{Fits events distances to the nearest sub-panel edge using glm}
\usage{
glm_events_dist_edge_min(detector, incl_event_list = NA)
}
\arguments{
\item{detector}{Detector object}

\item{incl_event_list}{a list of events to be included}
}
\value{
Fitted model
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{dist_corner}
\alias{dist_corner}
\title{A function to calculate pixel distances from the closest corner}
\usage{
dist_corner(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Matrix containing pixel distances from closest corner
}
\description{
A function to calculate pixel distances from the closest corner
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Calculate distance from pixels to corners
dist_corner(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{.plot_density}
\alias{.plot_density}
\title{Plots module
Plots density}
\usage{
.plot_density(
  ppp_obj,
  caption,
  file_path = NA,
  adjust = 0.5,
  color = topo.colors(50)
)
}
\arguments{
\item{ppp_obj}{ppp object}

\item{caption}{caption of the figure}

\item{file_path}{file path}

\item{adjust}{Kernel density bandwidth}

\item{color}{a list of colors}
}
\description{
Plots module
Plots density
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{plot_module_events}
\alias{plot_module_events}
\title{Plots damaged detector module events}
\usage{
plot_module_events(
  detector,
  col,
  row,
  file_path = NA,
  caption = TRUE,
  incl_event_list = NA
)
}
\arguments{
\item{detector}{Detector object}

\item{col}{Module column number}

\item{row}{Module row number}

\item{file_path}{Output file path}

\item{caption}{Flag to turn on/off figure caption}

\item{incl_event_list}{a list of events to be included}
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
 #Plot module events
 plot_module_events(detector_pilatus, 1, 1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{Pilatus_Detector}
\alias{Pilatus_Detector}
\title{A S3 class to represent the PerkinElmerRefurbished detector.}
\usage{
Pilatus_Detector()
}
\value{
Pilatus detector object
}
\description{
A S3 class to represent the PerkinElmerRefurbished detector.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{find_clumps}
\alias{find_clumps}
\title{Locates and classifies clumps (Events) of damaged pixels in a detector.
Adds an events matrix attribute to a \code{detector}.}
\usage{
find_clumps(detector, row = NA, col = NA)
}
\arguments{
\item{detector}{Detector object}

\item{row}{Module row number}

\item{col}{Module column number}
}
\value{
Detector with events matrix attribute
}
\description{
Clumps are a collection of neighboring dead pixels.
Pixels are classified as a clumps if they share an edge.
A mathematical definition can be found \href{https://warwick.ac.uk/fac/sci/statistics/crism/research/17-02/17-02w.pdf}{here}.
See details for classifications types.
}
\details{
\tabular{ll}{
   Classification \tab Description \cr
   Singlton \tab Single connecting pixel. \cr
   Doublets \tab Two connecting pixels. \cr
   Triplets \tab Three connecting pixels. \cr
   Upper horizontal line \tab Horizontal line in upper half of detector. \cr
   Lower horizontal line \tab Horizontal line in lower half of detector. \cr
   Left vertical line \tab Vertical line on left side of detector. \cr
   Right vertical line \tab Vertical line on right side of detector. \cr
   Larger cluster \tab Any larger grouping of damaged pixels. \cr
}


Lines have a pixel width of 1 but my include up to 10\% stray adjacent pixels.
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
 #Find clumps
detector_pilatus_w_clumps <- find_clumps(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{plot_pixels_kfg}
\alias{plot_pixels_kfg}
\title{Plots K, F, G functions}
\usage{
plot_pixels_kfg(
  detector,
  func,
  file_path = NA,
  row = NA,
  col = NA,
  caption = TRUE
)
}
\arguments{
\item{detector}{Detector object}

\item{func}{Function name ("K', "F", or "G")}

\item{file_path}{Output file path}

\item{row}{module row number}

\item{col}{module column number}

\item{caption}{Flag to turn on/off figure caption}
}
\description{
Plots K, F, G functions
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_pixels_kfg(detector_pilatus, "K")
plot_pixels_kfg(detector_pilatus, "F")
plot_pixels_kfg(detector_pilatus, "G")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{dist_edge_col}
\alias{dist_edge_col}
\title{Calculate horizontal distance from each pixel to nearest module edge}
\usage{
dist_edge_col(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
distance matrix
}
\description{
Calculate horizontal distance from each pixel to nearest module edge
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Calculate horizontal distance from each pixel to nearest module edge
dist_edge_col(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{is.detector}
\alias{is.detector}
\title{Check that \code{x} is an S3 detector class}
\usage{
is.detector(x)
}
\arguments{
\item{x}{Any variable}
}
\value{
True if x is an instance of detector
}
\description{
Check that \code{x} is an S3 detector class
}
\examples{
x <- Excalibur_Detector()
is.detector(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{.classify_clump}
\alias{.classify_clump}
\title{Clasifies a clump}
\usage{
.classify_clump(detector, x, y)
}
\arguments{
\item{detector}{Detector object}

\item{x}{vector containing the x coordinates of a clump}

\item{y}{vector containing the y coordinates of a clump}
}
\value{
the class of a clump (1 - singleton, 2 - double, 3 - triplet,
4 - larger cluster, unless it actually has the shape of a line,
5 (6): vertical line where closest edge is the upper (lower) one,
7 (8): horizontal line where closest edge is the right (left) one)
}
\description{
Clasifies a clump
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{.clump_module}
\alias{.clump_module}
\title{Identifying modules for clumps}
\usage{
.clump_module(detector, rrc)
}
\arguments{
\item{detector}{Detector object}

\item{rrc}{raster clumps objects}
}
\value{
data frame of the modules relating the clump
}
\description{
Identifying modules for clumps
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{.norm_vec}
\alias{.norm_vec}
\title{Estimates the norm of a vector}
\usage{
.norm_vec(v)
}
\arguments{
\item{v}{vector}
}
\value{
norm of the vector v
}
\description{
Estimates the norm of a vector
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{pixel_dist_ctr_linf}
\alias{pixel_dist_ctr_linf}
\title{Calculate parallel maxima from the centre for each pixel}
\usage{
pixel_dist_ctr_linf(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Matrix of parallel maxima
}
\description{
Calculate parallel maxima from the centre for each pixel
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
pixel_dist_ctr_linf(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{glm_events_dist_edge_row}
\alias{glm_events_dist_edge_row}
\title{Fits events distances from the module edges by row using glm}
\usage{
glm_events_dist_edge_row(detector, incl_event_list = NA)
}
\arguments{
\item{detector}{Detector object}

\item{incl_event_list}{a list of events to be included}
}
\value{
Fitted model
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{glm_pixel_dist_edge_col}
\alias{glm_pixel_dist_edge_col}
\title{Predict dead pixels from pixel distances from the module edges by module column}
\usage{
glm_pixel_dist_edge_col(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Fitted model
}
\description{
Fit a logistic regression model using \code{glm}.
Predict dead pixels from pixel distances from the module edges by module column
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Fit logistic regression model
glm_pixel_dist_edge_col(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{remove_high_density_cluster}
\alias{remove_high_density_cluster}
\title{Remove high density cluster of dead pixels}
\usage{
remove_high_density_cluster(detector, min_pts = 30, eps_adjust = 0.05)
}
\arguments{
\item{detector}{Detector object}

\item{min_pts}{minimum points argument of dbscan function}

\item{eps_adjust}{adjust eps}
}
\value{
detector object with high density cluster of pixels removed
}
\description{
In some situations, the analysis may be dominated by an area of elevated damage.
The investigation of complete spatial randomness then becomes uninformative.
The area with elevated damage can be removed.
Dead statistics and clumps are recalculate if they were present in the Detector object.
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
# Find events
detector_pilatus_events <- find_clumps(detector_pilatus)
# Remove high density clusters
detector_pilatus_modified <- remove_high_density_cluster(detector_pilatus, min_pts = 30, eps_adjust = 0.05)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{glm_pixel_dist_edge_row}
\alias{glm_pixel_dist_edge_row}
\title{Predict dead pixels from pixel distances from the module edges by module row}
\usage{
glm_pixel_dist_edge_row(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Fitted model
}
\description{
Fit a logistic regression model using \code{glm}.
Predict dead pixels from pixel distances from the module edges by module row
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Fit logistic regression model
glm_pixel_dist_edge_row(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{.matrix_from_hdf}
\alias{.matrix_from_hdf}
\title{Reads in hdf file(s) and returns a pixel matrix}
\usage{
.matrix_from_hdf(detector, file_path)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{A list of paths to hdf files. Must be in the correct order.}
}
\value{
Data of a combined dataset from hdf files
}
\description{
Reads in hdf file(s) and returns a pixel matrix
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{plot_events_density}
\alias{plot_events_density}
\title{Plots density graph of events of a detector or module}
\usage{
plot_events_density(
  detector,
  file_path = NA,
  adjust = 0.5,
  row = NA,
  col = NA,
  caption = TRUE,
  incl_event_list = NA,
  color = topo.colors(50)
)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}

\item{adjust}{Kernel density bandwidth}

\item{row}{Module row number}

\item{col}{Module column number}

\item{caption}{Flag to turn on/off figure caption}

\item{incl_event_list}{a list of events to be included}

\item{color}{a list of colors}
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
\examples{
detector_perkinfull <- create_detector("PerkinElmerFull")
file_path <-  system.file("extdata", "PerkinElmerFull",
                          "BadPixelMap_t1.bpm.xml",
                          package = "detectorchecker")
detector_perkinfull <- load_pix_matrix(
detector = detector_perkinfull, file_path = file_path)
detector_perkinfull_events = find_clumps(detector_perkinfull)
plot_events_density(detector_perkinfull_events)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{plot_pixel_dist_edge_row}
\alias{plot_pixel_dist_edge_row}
\title{Calculates and plots vetical distances from the module edges}
\usage{
plot_pixel_dist_edge_row(detector, file_path = NA)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}
}
\description{
Calculates and plots vetical distances from the module edges
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_pixel_dist_edge_row(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{.plot_counts}
\alias{.plot_counts}
\title{Plots dead pixel counts}
\usage{
.plot_counts(module_count_arr, caption, file_path = NA)
}
\arguments{
\item{module_count_arr}{Counts per array}

\item{caption}{caption of the figure}

\item{file_path}{file path}
}
\description{
Plots dead pixel counts
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{plot.detector}
\alias{plot.detector}
\title{Plot detector}
\usage{
\method{plot}{detector}(detector, file_path = NA, caption = TRUE)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}

\item{caption}{Flag to turn on/off figure caption}
}
\description{
Plot detector
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{get_dead_stats}
\alias{get_dead_stats}
\title{Generate summary of damaged pixels and add as a dead_stats attribute to the detector object}
\usage{
get_dead_stats(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Detector object with dead_stats attribute
}
\description{
Generate summary of damaged pixels and add as a dead_stats attribute to the detector object
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Calculate dead_stats
detector_pilatus <- get_dead_stats(detector_pilatus)
# Print dead stats
print(detector_pilatus$dead_stats)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{plot_events_arrows}
\alias{plot_events_arrows}
\title{Plots arrows graph of events of a detector or module}
\usage{
plot_events_arrows(
  detector,
  file_path = NA,
  row = NA,
  col = NA,
  caption = TRUE,
  incl_event_list = NA
)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}

\item{row}{Module row number}

\item{col}{Module column number}

\item{caption}{Flag to turn on/off figure caption}

\item{incl_event_list}{a list of events to be included}
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{PerkinElmerFull_Detector}
\alias{PerkinElmerFull_Detector}
\title{A S3 class to represent the PerkinElmerFull detector.}
\usage{
PerkinElmerFull_Detector()
}
\value{
PerkinElmerFul detector object
}
\description{
A S3 class to represent the PerkinElmerFull detector.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{plot_pixels_arrows}
\alias{plot_pixels_arrows}
\title{A function to plot NN oriented arrows of dead pixels of detector or module}
\usage{
plot_pixels_arrows(
  detector,
  file_path = NA,
  row = NA,
  col = NA,
  caption = TRUE
)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}

\item{row}{Module row number}

\item{col}{Module column number}

\item{caption}{Flag to turn on/off figure caption}
}
\description{
A function to plot NN oriented arrows of dead pixels of detector or module
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_pixels_arrows(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.getmode}
\alias{.getmode}
\title{Returns the mode of a set of data}
\usage{
.getmode(v)
}
\arguments{
\item{v}{set of data}
}
\value{
uniqv the value of the mode
}
\description{
Returns the mode of a set of data
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{Dead_Stats}
\alias{Dead_Stats}
\title{A list to represent dead pixels statistics summary.
Added to a detector by \code{get_dead_stats}.}
\usage{
Dead_Stats(
  dead_n = NA,
  module_n = NA,
  module_count_arr = NA,
  module_count = NA,
  avg_dead_mod = NA,
  Chisq_s = NA,
  Chisq_df = NA,
  Chisq_p = NA
)
}
\arguments{
\item{dead_n}{Total number of damaged pixels:}

\item{module_n}{Total number of modules}

\item{module_count_arr}{Count of dead pixels in each quadrat}

\item{module_count}{Count of dead pixels in each quadrat}

\item{avg_dead_mod}{Average number of damaged pixels per module}

\item{Chisq_s}{The Chi-Squared test statistic value}

\item{Chisq_df}{Chi-Squared degrees of freedom}

\item{Chisq_p}{Chi-Squared p-value}
}
\value{
Dead_Stats list
}
\description{
A list to represent dead pixels statistics summary.
Added to a detector by \code{get_dead_stats}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{glm_events_dist_corner}
\alias{glm_events_dist_corner}
\title{Fits events distances to the nearest corner using glm}
\usage{
glm_events_dist_corner(detector, incl_event_list = NA)
}
\arguments{
\item{detector}{Detector object}

\item{incl_event_list}{a list of events to be included}
}
\value{
Fitted model
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{plot_events_angles}
\alias{plot_events_angles}
\title{Plots angles graph of events of a detector or module}
\usage{
plot_events_angles(
  detector,
  file_path = NA,
  row = NA,
  col = NA,
  caption = TRUE,
  incl_event_list = NA
)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}

\item{row}{Module row number}

\item{col}{Module column number}

\item{caption}{Flag to turn on/off figure caption}

\item{incl_event_list}{a list of events to be included}
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
\examples{
detector_perkinfull <- create_detector("PerkinElmerFull")
file_path <-  system.file("extdata", "PerkinElmerFull",
                          "BadPixelMap_t1.bpm.xml",
                          package = "detectorchecker")
detector_perkinfull <- load_pix_matrix(
detector = detector_perkinfull, file_path = file_path)
detector_perkinfull_events = find_clumps(detector_perkinfull)
plot_events_angles(detector_perkinfull_events)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.tr}
\alias{.tr}
\title{Utils module
Calculates the trace value of a square matrix}
\usage{
.tr(m)
}
\arguments{
\item{m}{A square matrix}
}
\value{
tr The trace value
}
\description{
Utils module
Calculates the trace value of a square matrix
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{.assign_module}
\alias{.assign_module}
\title{Pixel module
Function assign a module to each dead pixel}
\usage{
.assign_module(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
dead_modules
}
\description{
Pixel module
Function assign a module to each dead pixel
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.extract_detector_parameter}
\alias{.extract_detector_parameter}
\title{Checks whether a detector parameter is in the file string}
\usage{
.extract_detector_parameter(file_string, parameter)
}
\arguments{
\item{file_string}{String of a file context}

\item{parameter}{Detector parameter}
}
\value{
parameter value
}
\description{
Checks whether a detector parameter is in the file string
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{summary.detector}
\alias{summary.detector}
\title{Summary of a detector object}
\usage{
\method{summary}{detector}(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
String with the detector summary
}
\description{
Summary of a detector object
}
\examples{
detc <- create_detector("Pilatus")
summary(detc)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{.derive_detector}
\alias{.derive_detector}
\title{Deriving additional detector elements
Conditions additional elements of Detector object that are frequently used later
They are calculated from parameters defined in examples
Matrices that contains xy coordiantes of edges of modules
By definition, edges are part of modules (not part of gaps)
i.e. for each module two pairs: first/last col and first/last row.}
\usage{
.derive_detector(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Detector object
}
\description{
Deriving additional detector elements
Conditions additional elements of Detector object that are frequently used later
They are calculated from parameters defined in examples
Matrices that contains xy coordiantes of edges of modules
By definition, edges are part of modules (not part of gaps)
i.e. for each module two pairs: first/last col and first/last row.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{glm_pixel_dist_corner}
\alias{glm_pixel_dist_corner}
\title{Predict dead pixels from pixel distances to the nearest corner}
\usage{
glm_pixel_dist_corner(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Fitted model
}
\description{
Fit a logistic regression model using \code{glm}.
Predict dead pixels from pixel distances to the nearest corner
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Fit logistic regression model
glm_pixel_dist_corner(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{load_pix_matrix}
\alias{load_pix_matrix}
\title{A function to load pixel data and set as attribute on a detector}
\usage{
load_pix_matrix(detector, file_path)
}
\arguments{
\item{detector}{The name of the detector object to be used}

\item{file_path}{Path(s) to the file(s) containing dead pixel information}
}
\value{
Detector object
}
\description{
A function to load pixel data and set as attribute on a detector
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{glm_events_ctr_eucl}
\alias{glm_events_ctr_eucl}
\title{Fits events distances from the centre using glm}
\usage{
glm_events_ctr_eucl(detector, incl_event_list = NA)
}
\arguments{
\item{detector}{Detector object}

\item{incl_event_list}{a list of events to be included}
}
\value{
Fitted model
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{.create_ppp_gaps_row}
\alias{.create_ppp_gaps_row}
\title{Creates ppp object of vertical gaps}
\usage{
.create_ppp_gaps_row(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Point pattern dataset
}
\description{
Creates ppp object of vertical gaps
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{.dist_closest_edge}
\alias{.dist_closest_edge}
\title{A function to calculate closest distance to an edge for a pixel}
\usage{
.dist_closest_edge(x, size)
}
\arguments{
\item{x}{Coordinate of pixel}

\item{size}{Size of module}
}
\value{
distance to closest edge
}
\description{
A function to calculate closest distance to an edge for a pixel
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{plot_module_pixels}
\alias{plot_module_pixels}
\title{A function to plot detector module with damaged pixels}
\usage{
plot_module_pixels(detector, col, row, file_path = NA, caption = TRUE)
}
\arguments{
\item{detector}{Detector object}

\item{col}{Module column number}

\item{row}{Module row number}

\item{file_path}{Output file path}

\item{caption}{Flag to turn on/off figure caption}
}
\description{
A function to plot detector module with damaged pixels
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_module_pixels(detector_pilatus, 1, 1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{plot_pixel_ctr_linf}
\alias{plot_pixel_ctr_linf}
\title{Calculates and plots pixel parallel maxima from the centre}
\usage{
plot_pixel_ctr_linf(detector, file_path = NA)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}
}
\description{
Calculates and plots pixel parallel maxima from the centre
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
plot_pixel_ctr_linf(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{.dead_pix_coords}
\alias{.dead_pix_coords}
\title{Extracts a table of dead pixel coordinates from a pixel matrix}
\usage{
.dead_pix_coords(pix_matrix)
}
\arguments{
\item{pix_matrix}{pixel matrix with dead pixels flagged with 1}
}
\value{
Table containing dead pixel coordinates
}
\description{
Extracts a table of dead pixel coordinates from a pixel matrix
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{.check_clumps}
\alias{.check_clumps}
\title{Checks if correct clumps were found. If not, finds clumps}
\usage{
.check_clumps(detector, row = NA, col = NA)
}
\arguments{
\item{detector}{Detector object}

\item{row}{Module row number}

\item{col}{Module column number}
}
\value{
detector_events Detector object
}
\description{
Checks if correct clumps were found. If not, finds clumps
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{.plot_angles}
\alias{.plot_angles}
\title{Plots nearest neighbor angles}
\usage{
.plot_angles(ppp_obj, caption, file_path = NA)
}
\arguments{
\item{caption}{caption of the figure}

\item{file_path}{file path}

\item{ppp_object}{spatstat ppp object}
}
\description{
Plots nearest neighbor angles
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{get_events_matrix}
\alias{get_events_matrix}
\title{Generates events matrix for selected events}
\usage{
get_events_matrix(detector, incl_event_list = NA)
}
\arguments{
\item{detector}{Detector object}

\item{incl_event_list}{a list of events to be included}
}
\value{
Events matrix
}
\description{
Generates events mask (a matrix with pixels as 0 and events as 1) indicating if a pixel is in an event
as calculated by \code{find_clumps(detc)}
Note that the parameter incl_event_list is a list of integers from 1-8. For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Identify events
detector_pilatus <- find_clumps(detector_pilatus)
# Get events matrix
get_events_matrix(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{.extract_number}
\alias{.extract_number}
\title{Internal function to convert string values to numbers}
\usage{
.extract_number(s)
}
\arguments{
\item{s}{String expression?}
}
\value{
Numeric value
}
\description{
Internal function to convert string values to numbers
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{dist_edge_row}
\alias{dist_edge_row}
\title{Calculate vertical distance from each pixel to nearest module edge}
\usage{
dist_edge_row(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
distance matrix
}
\description{
Calculate vertical distance from each pixel to nearest module edge
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Calculate vertical distance from each pixel to nearest module edge
dist_edge_row(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module.R
\name{.check_select}
\alias{.check_select}
\title{Checks if the selected row and column are within the boundaries of the detector}
\usage{
.check_select(detector, row, col)
}
\arguments{
\item{detector}{Detector object}

\item{row}{module row}

\item{col}{module col}
}
\value{
Boolean
}
\description{
Checks if the selected row and column are within the boundaries of the detector
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{.get_clump_pixel_ppp}
\alias{.get_clump_pixel_ppp}
\title{Creates ppp object for damaged detector pixels}
\usage{
.get_clump_pixel_ppp(detector, incl_event_list = NA)
}
\arguments{
\item{detector}{Detector object}

\item{incl_event_list}{a list of events to be included}
}
\value{
ppp object for damaged detector pixels
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{plot_events}
\alias{plot_events}
\title{Plots damaged detector events}
\usage{
plot_events(
  detector,
  col = NA,
  row = NA,
  file_path = NA,
  caption = TRUE,
  incl_event_list = NA,
  plot_edges_gaps = TRUE
)
}
\arguments{
\item{detector}{Detector object}

\item{col}{Module column number}

\item{row}{Module row number}

\item{file_path}{Output file path}

\item{caption}{Flag to turn on/off figure caption}

\item{incl_event_list}{a list of events to be included}

\item{plot_edges_gaps}{Plot edges gaps}
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
 #Plot events
 plot_events(detector_pilatus)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{.xyc_ply_func}
\alias{.xyc_ply_func}
\title{Clasifies clumps with respect to xy coordinates.}
\usage{
.xyc_ply_func(detector, xyc_pixel_df)
}
\arguments{
\item{detector}{Detector object}

\item{xyc_pixel_df}{xyc_pixel_df}
}
\value{
data frame with clasification results
}
\description{
Clasifies clumps with respect to xy coordinates.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{plot_pixels_count}
\alias{plot_pixels_count}
\title{A function to plot detector with dead pixel counts per module}
\usage{
plot_pixels_count(detector, file_path = NA, row = NA, col = NA, caption = TRUE)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}

\item{row}{Module row number}

\item{col}{Module column number}

\item{caption}{Flag to turn on/off figure caption}
}
\description{
A function to plot detector with dead pixel counts per module
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(detector = detector_pilatus, file_path = file_path)
detector_pilatus_damage <- get_dead_stats(detector_pilatus)
plot_pixels_count(detector_pilatus_damage)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module.R
\name{which_module}
\alias{which_module}
\title{Module module
Returns row or column of a module that a dead pixel belongs to}
\usage{
which_module(coo, me)
}
\arguments{
\item{coo}{x or y coordinate of a dead pixel}

\item{me}{module edges}
}
\value{
row or column number
}
\description{
Module module
Returns row or column of a module that a dead pixel belongs to
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{create_detector}
\alias{create_detector}
\title{Create a Detector object}
\usage{
create_detector(detector_name)
}
\arguments{
\item{detector_name}{The name of the detector}
}
\value{
Detector S3 object
}
\description{
Create a Detector object.
If the \code{detector_name} is not available will raise an exception.
}
\examples{
detc <- create_detector(available_detectors[1])
detc <- create_detector("Pilatus")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{plot_events_count}
\alias{plot_events_count}
\title{Plots events count per detector or module}
\usage{
plot_events_count(
  detector,
  file_path = NA,
  row = NA,
  col = NA,
  caption = TRUE,
  incl_event_list = NA
)
}
\arguments{
\item{detector}{Detector object}

\item{file_path}{Output file path}

\item{row}{Module row number}

\item{col}{Module column number}

\item{caption}{Flag to turn on/off figure caption}

\item{incl_event_list}{a list of events to be included}
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
\examples{
detector_perkinfull <- create_detector("PerkinElmerFull")
file_path <-  system.file("extdata", "PerkinElmerFull",
                          "BadPixelMap_t1.bpm.xml",
                          package = "detectorchecker")
detector_perkinfull <- load_pix_matrix(
detector = detector_perkinfull, file_path = file_path)
detector_perkinfull_events = find_clumps(detector_perkinfull)
plot_events_count(detector_perkinfull_events)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{.get_ppp_dead_module}
\alias{.get_ppp_dead_module}
\title{Generates ppp for the dead pixels for a selected module}
\usage{
.get_ppp_dead_module(detector, row, col)
}
\arguments{
\item{detector}{Detector object}

\item{row}{module row number}

\item{col}{module column number}
}
\value{
ppp of dead pixels
}
\description{
Generates ppp for the dead pixels for a selected module
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{glm_pixel_dist_edge_min}
\alias{glm_pixel_dist_edge_min}
\title{Predict dead pixels from pixel distances to the nearest sub-panel edge}
\usage{
glm_pixel_dist_edge_min(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Fitted model
}
\description{
Fit a logistic regression model using \code{glm}.
Predict dead pixels from pixel distances to the nearest sub-panel edge
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Fit logistic regression model
glm_pixel_dist_edge_min(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{glm_events_dist_edge_col}
\alias{glm_events_dist_edge_col}
\title{Fits events distances from the module edges by column using glm}
\usage{
glm_events_dist_edge_col(detector, incl_event_list = NA)
}
\arguments{
\item{detector}{Detector object}

\item{incl_event_list}{a list of events to be included}
}
\value{
Fitted model
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{PerkinElmerCropped1600_Detector}
\alias{PerkinElmerCropped1600_Detector}
\title{A S3 class to represent the PerkinElmerCropped1600 detector.}
\usage{
PerkinElmerCropped1600_Detector()
}
\value{
PerkinElmerCropped1600 detector object
}
\description{
A S3 class to represent the PerkinElmerCropped1600 detector.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clump.R
\name{glm_events_ctr_linf}
\alias{glm_events_ctr_linf}
\title{Fits events parallel maxima from the centre using glm}
\usage{
glm_events_ctr_linf(detector, incl_event_list = NA)
}
\arguments{
\item{detector}{Detector object}

\item{incl_event_list}{a list of events to be included}
}
\value{
Fitted model
}
\description{
Note that the parameter incl_event_list is a list of integers from 1-8.  For all events this would be list(1,2,3,4,5,6,7,8).
The integer values declare which damage events to include:
1=singletons, 2=doublets, 3=triplets, 4=larger clusters, 5=upper horizontal lines,
6=lower horizontal lines, 7=left vertical lines, 8=right vertical lines
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{.get_detector_ppps}
\alias{.get_detector_ppps}
\title{Generate detector ppps for edges and gaps}
\usage{
.get_detector_ppps(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
a list of ppps for edges and gaps
}
\description{
Generate detector ppps for edges and gaps
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{.create_ppp_edges_col}
\alias{.create_ppp_edges_col}
\title{This is the ppp_edges_col creation function}
\usage{
.create_ppp_edges_col(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Point pattern dataset
}
\description{
This is the ppp_edges_col creation function
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{glm_pixel_ctr_eucl}
\alias{glm_pixel_ctr_eucl}
\title{Predict dead pixels from the pixel's euclidean distance from the detector center}
\usage{
glm_pixel_ctr_eucl(detector)
}
\arguments{
\item{detector}{Detector object}
}
\value{
Fitted model
}
\description{
Fit a logistic regression model using \code{glm}.
Predicts dead pixels from the pixel's euclidean distance from the detector center
}
\examples{
# Create a detector
detector_pilatus <- create_detector("Pilatus")
# Load a pixel matrix
file_path <-  system.file("extdata", "Pilatus", "badpixel_mask.tif",
                         package ="detectorchecker")
detector_pilatus <- load_pix_matrix(
 detector = detector_pilatus, file_path = file_path)
# Fit logistic regression model
glm_pixel_ctr_eucl(detector_pilatus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixel.R
\name{.plot_pixel}
\alias{.plot_pixel}
\title{Plots pixel distance analysis}
\usage{
.plot_pixel(data, width, height, file_path = NA)
}
\arguments{
\item{data}{Matrix containing pixel analysis data}

\item{width}{Plot width}

\item{height}{Plot height}

\item{file_path}{Output path with an extension}
}
\description{
Plots pixel distance analysis
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{PerkinElmerRefurbished_Detector}
\alias{PerkinElmerRefurbished_Detector}
\title{A S3 class to represent the PerkinElmerRefurbished detector.}
\usage{
PerkinElmerRefurbished_Detector()
}
\value{
PerkinElmerRefurbished detector object
}
\description{
A S3 class to represent the PerkinElmerRefurbished detector.
}
