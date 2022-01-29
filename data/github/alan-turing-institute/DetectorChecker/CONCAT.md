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
