[![Build Status](https://travis-ci.org/ropensci/colocr.svg?branch=master)](https://travis-ci.org/ropensci/colocr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/colocr?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/colocr)
[![codecov](https://codecov.io/gh/ropensci/colocr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/colocr)
[![Build Status](https://travis-ci.org/MahShaaban/colocr_app.svg?branch=master)](https://travis-ci.org/MahShaaban/colocr_app)
[![status](https://img.shields.io/badge/shinyapps.io-running-green.svg)](https://mahshaaban.shinyapps.io/colocr_app2/) 
[![](https://badges.ropensci.org/243_status.svg)](https://github.com/ropensci/onboarding/issues/243)

# colocr

An R package for conducting co-localization analysis.

## Overview

A few R packages are available for conducting image analysis, which is a very wide topic. As a result, some of us might feel at a loss when all they want to do is a simple co-localization calculations on a small number of microscopy images. This package provides a simple straight forward workflow for loading images, choosing regions of interest (ROIs) and calculating co-localization statistics. Included in the package, is a [shiny app](https://shiny.rstudio.com) that can be invoked locally to interactively select the regions of interest in a semi-automatic way. The package is based on the R package [`imager`](https://CRAN.R-project.org/package=imager).


## Installing `colocr`
`colocr` is available on [CRAN](https://cran.r-project.org) 
and can be installed using

```r
# install from cran
install.packages('colocr')
```

The package development version is available at [github](https://github.com/ropensci/colocr).

```r
# install from github
devtools::install_github('ropensci/colocr')
```

This package depends on `imager` which has some external dependencies. The instructions for installing `imager` can be found
[here](https://github.com/dahtah/imager).

## Getting started

To get started, load the required packages and the images. The images below
are from [DU145](https://en.wikipedia.org/wiki/DU145) cell line and were 
stained for two proteins; [RKIP](https://en.wikipedia.org/wiki/Raf_kinase_inhibitor_protein) and [LC3](https://en.wikipedia.org/wiki/MAP1LC3B).
Then, apply the appropriate parameters for choosing the regions of interest
using the `roi_select`. Finally, check the appropriateness of the 
parameters by highlighting the ROIs on the image.

```r
# load libraries
library(colocr)

# load images
fl <- system.file('extdata', 'Image0001_.jpg', package = 'colocr')
img <- image_load(fl)

# select ROI and show the results
par(mfrow = c(2,2), mar = rep(1, 4))

img %>%
  roi_select(threshold = 90) %>%
  roi_show()
```

The same can be achieved interactively using an accompanying **shiny** app.
To launch the app run.

```r
run_app()
```

The reset of the analysis depends on the particular kind of images. Now, `colocr`
implements two simple co-localization statistics; Pearson's Coefficient Correlation [(PCC)](https://www.ncbi.nlm.nih.gov/pubmed/20653013) and the Manders Overlap Coefficient [(MOC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3074624/).

To apply both measures of correlation, we first get the pixel intensities and call `roi_test` on the merge image.

```r
# calculate co-localization statistics
img %>%
  roi_select(threshold = 90) %>%
  roi_test(type = 'both')
```

The same analysis and more can be conducted using a web interface for the package available [here](https://mahshaaban.shinyapps.io/colocr_app2/)

## Acknowledgement

* The vignette images from [Lai Huyen Trang](https://www.researchgate.net/profile/Lai_Huyen_Trang)  
* The test examples from [Colocalization Benchmark Source (CBS)](https://www.colocalization-benchmark.com/index.html)  
* The implementation of the co-localization statistics from [Dunn et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3074624/)  

## Citation

```r
citation('colocr')
```

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# colocr 0.1.1

* Skip app testing on CRAN
* Add related article citation info

# colocr 0.1.0

* Packege added ROpenSci



# All contributions are most welcomed

I'd be glad to recieve any comments or ideas to help the package forward.

## Bugs

* To report a bug please use the [issue](https://github.com/ropensci/colocr/issues) page on github

## Code contributions

* Fork this repo to your github account
* Clone the repo to your machine and make changes
* Push to your account
* Submit a [pull](https://github.com/ropensci/colocr/pulls) request at this repo

## Email: [mahmoud.s.fahmy@students.kasralainy.egu.eg](mahmoud.s.fahmy@students.kasralainy.egu.eg)

## Thank you for contributing!
![version](https://img.shields.io/badge/version-v1.0-blue.svg)
[![Build Status](https://travis-ci.org/MahShaaban/colocr_app.svg?branch=master)](https://travis-ci.org/MahShaaban/colocr_app)
[![status](https://img.shields.io/badge/shinyapps.io-running-green.svg)](https://mahshaaban.shinyapps.io/colocr_app2/)  

# colocr_app

A shiny app for conducting co-localization analysis.
---
title: "Using colocr"
subtitle: "An R package for conducting co-localization analysis"
author: "Mahmoud Ahmed"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    self_contained: yes
vignette: >
  %\VignetteIndexEntry{Using colocr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = 'center',
  message = FALSE
)
```

## Overview

A few R packages are available for conducting image analysis, which is a very wide topic. As a result, some of us might feel at a loss when all they want to do is a simple co-localization calculations on a small number of microscopy images. This package provides a simple straight forward workflow for loading images, choosing regions of interest (ROIs) and calculating co-localization statistics. Included in the package, is a [shiny app](https://shiny.rstudio.com) that can be invoked locally to interactively select the regions of interest in a semi-automatic way. The package is based on the R package [`imager`](https://CRAN.R-project.org/package=imager).

## Getting started

To get started, load the required packages and the images. The images below
are from [DU145](https://en.wikipedia.org/wiki/DU145) cell line and were 
stained for two proteins; [RKIP](https://en.wikipedia.org/wiki/Raf_kinase_inhibitor_protein) and [LC3](https://en.wikipedia.org/wiki/MAP1LC3B).
Then, apply the appropriate parameters for choosing the regions of interest
using the `roi_select`. Finally, check the appropriateness of the 
parameters by highlighting the ROIs on the image.

```{r getting_started, fig.width=7, fig.height=7}
# load libraries
library(colocr)

# load images
fl <- system.file('extdata', 'Image0001_.jpg', package = 'colocr')
img <- image_load(fl)

# select ROI and show the results
par(mfrow = c(2,2), mar = rep(1, 4))

img %>%
  roi_select(threshold = 90) %>%
  roi_show()
```

The same can be achieved interactively using an accompanying **shiny** app.
To launch the app run.

```{r colocr_app, eval=FALSE}
colocr_app()
```

The rest of the analysis depends on the particular kind of images. Now, `colocr`
implements two simple co-localization statistics; Pearson's Coefficient Correlation [(PCC)](https://www.ncbi.nlm.nih.gov/pubmed/20653013) and the Manders Overlap Coefficient [(MOC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3074624/).

To apply both measures of correlation, we first get the pixel intensities and call `roi_test` on the merge image.

```{r call_roi_test}
# calculate co-localization statistics
img %>%
  roi_select(threshold = 90) %>%
  roi_test(type = 'both')
```

The same analysis and more can be conducted using a web interface for the package available [here](https://mahshaaban.shinyapps.io/colocr_app2/)

## Detailed Example

The following example uses images from the [DU145](https://en.wikipedia.org/wiki/DU145) prostate cancer cell line. In this experiment, the cell line was treated with probes for two proteins [RKIP](https://en.wikipedia.org/wiki/Raf_kinase_inhibitor_protein) and [LC3](https://en.wikipedia.org/wiki/MAP1LC3B). The aim of this experiment is to determine, how much of the two proteins are co-localized or co-distributed in the particular cell line.  

```{r show_images, fig.width=7, fig.height=2.5}
# load required libraries
library(colocr)

# load images
img <- image_load(system.file('extdata', 'Image0001_.jpg', package = 'colocr'))       # merge
img1 <- imager::channel(img, 1)  # red
img2 <- imager::channel(img, 2)  # green

# show images
par(mfrow = c(1,3), mar = rep(1,4))
plot(img, axes = FALSE, main = 'Merged')
plot(img1, axes = FALSE, main = 'RKIP')
plot(img2, axes = FALSE, main = 'LC3')
```

The `colocr` package provides a straight forward workflow for determining the amount of co-localization. This workflow consists of two steps only:

  * Choosing regions of interest (ROI)
  * Calculating the correlation between the pixel intensities

The first step can be achieved by calling `roi_select` on the image. In addition, `roi_show` can be used to visualize the regions that were selected to make sure they
match the expectations. Similarly, `roi_check` can be used to visualize the 
pixel intensities of the selected regions. The second step is calling `roi_test`
to calculate the co-localization statistics.

The calls to these functions can be piped using `%>%` to reduce the amount of typing
and make the code more readable.  

```{r workflow, echo=FALSE, out.width='300px'}
workflow_fig <- system.file('colocr_app', 'workflow.png', package = 'colocr')
knitr::include_graphics(workflow_fig)
```

### Choosing ROIs

The function `roi_select` relies on different algorithms from the `imager` package. However, using the functions to select the ROIs doesn't require any background knowledge in the workings of the algorithms and can be done through trying different parameters and choosing the most appropriate ones. Typically, one wants to select the 
regions of the image occupied by a cell or a group of cells. However, the package can
also be used to select certain areas/structures within the cell if they are distinct 
enough. By default, the largest contiguous region of the image is selected, more regions can be added using the argument `n`. The details of the other inputs are documented in the function help page `?roi_select`.
 
```{r call_roi_select}
# select regions of interest
img_rois <- img %>%
  roi_select(threshold = 90)
```

This function returns `cimg` object containing the original input image and an added attribute called `label` which indicates the selected regions. `label` is a vector of `integer`s; with 0 indicating the non-selected pixels and 1 for the selected regions. When the argument `n` is provided to `roi_select`, 1 is replaced by `integer` labels for each of the selected regions separately.

```{r output_str}
# class of the returned object
class(img); class(img_rois)

# name of added attribut
names(attributes(img)); names(attributes(img_rois))

# str of labels
label <- attr(img_rois, 'label')
str(label)

# unique labels 
unique(label)
```

Now, to make sure these parameters are appropriately encompassing the ROIs, call the `roi_show` to visualize side by side the original merge picture, a low resolution picture of the ROI and the images from the two different channels highlighted by the ROIs.

```{r call_roi_show, fig.width=7, fig.height=7}
# select ROI and show the results
par(mfrow = c(2,2), mar = rep(1, 4))

img_rois %>%
  roi_show()
```

Both the co-localization statistics implemented in this package quantify different aspects of the linear trend between the pixel intensities from the two channels of the image. Therefore, it is useful to visualize this trend and the distribution of the intensities to make sure the analysis is appropriate.

```{r call_coloc_show, fig.width=7, fig.height=3}
# show the scatter and density of the pixel values
par(mfrow=c(1,2), mar = c(4,4,1,1))

img_rois %>%
  roi_check()
```

Arguably, selecting the regions of interest is the most time consuming step in this kind of analysis. Usually, one has to do this selection by hand when using image analysis software such as [imageJ](https://imagej.nih.gov/ij/). This package only semi-automates this step, but still relies on the user's judgment on which parameters to use and whether or not the selected ROIs are appropriate. To make life easier, the package provides a simple shiny app to interactively determine these parameters and use it in the rest of the workflow. To launch the app run the following

```{r colocr_app2, eval=FALSE}
# run the shiny app
colocr_app()
```

And here is a screen shot from the app after applying the same parameters used previously.

```{r app_screenshot, echo=FALSE, out.width="600px"}
# show the screen shot of the shiny app
app_ss <- system.file('colocr_app/tests/one_image-expected', '001.png', package = 'colocr')
knitr::include_graphics(app_ss)
```

Although this app was designed to be invoked from within the package to help the users to choose the selection parameters interactively, it's a stand alone app and can run the same analysis described here. The app can be accessed from the web [here](https://mahshaaban.shinyapps.io/colocr_app2/)

### Calculating Correlation Statistics

The two different statistics implemented in this package are the PCC and SCC. The formal description and the rational for using them is detailed elsewhere. Invoking the test is a one function call on the selected regions of interest. 

```{r call_roi_test2}
# Calculate the co-localization statistics
tst <- img_rois %>%
  roi_test(type = 'both')
tst
```

`roi_test` returns a `data.frame` with a column for each of the desired statistics. When `n` is used in the selection of the regions of interest, a separate row is returned for each 
region.

```{r output_str2}
# str of the roi_test output
str(tst)
```

## Selecting a multiple regions in an image

In this example, we want to select the all cells as regions of interest 
and make sure no background is included in the selected region. To do that, we
use the different arguments in `roi_select`.

```{r roi_subset, fig.width=7, fig.height=7}
# load image
img2 <- image_load(system.file('extdata', 'Image0003_.jpg', package = 'colocr'))       # merge

# select ROI and show the results
par(mfrow = c(2,2), mar = rep(1, 4))

img2 %>%
  roi_select(threshold = 85,
             shrink = 10,
             clean = 10,
             n = 3) %>%
  roi_show()
```

## Analyzing a collection of images at once

To process a collection of images at once, we first make a `list` of the image objects
and pass it to the `roi_select`. Other arguments can be provided as a single value to
be applied to all images or specific values for each image.

```{r analyze_collection}
# make a list of images
fls <- c(system.file('extdata', 'Image0001_.jpg', package = 'colocr'),
         system.file('extdata', 'Image0003_.jpg', package = 'colocr'))
image_list <- image_load(fls)

# call roi_select on multiple images
image_list %>%
  roi_select(threshold = 90) %>%
  roi_test()
```

```{r analyze_collection2}
# make threshold input list
thresholds <- c(90, 95)

# call roi_select on multiple images and specific thresholds for each
image_list %>%
  roi_select(threshold = thresholds) %>%
  roi_test()
```

The same applies of the other two functions; `roi_show` and `roi_check`. When a `list`
of images is provided, they return same set of plots for each of the images.

## Reproducing the shiny app output

The interactive shiny app can be useful when trying the analysis for the first time.
Or when choosing the appropriate input values to select the regions of interest. 
To make things easier for the users, we enabled exporting the input and the output 
from the shiny app. These can be used to report an analysis directly, or to 
reproduce the results using R code.

This section describes how to use the exported input and output tables from the app
and using them to reproduce the same results. First, we start by loading `inputs_19.05.23_05.55.22.csv` which is the output of an analysis of two images 
performed by in the shiny app. Second, we load the `inputs_18.09.02_05.15.08.csv`
which contains the input parameters used in that analysis. Then, we call `roi_select`
and `roi_test` using these inputs. Finally, we compare the outputs form the app
and the script.

```{r read_output}
# show the output
stats <- read.csv(system.file('colocr_app', 'stats_19.05.23_05.55.21.csv', package = 'colocr'))

stats
```

```{r read_inputs}
# show the inputs
inputs <- read.csv(system.file('colocr_app', 'inputs_19.05.23_05.55.22.csv', package = 'colocr'), stringsAsFactors = FALSE)

inputs
```

```{r apply_inputs}
# read images
fls <- lapply(inputs$image, function(x) {
  system.file('extdata', x, package = 'colocr')
  })
imgs <- image_load(fls)

# use the app input to the roi_select function
rep_stats <- imgs %>%
  roi_select(threshold = inputs$threshold,
             shrink = inputs$shrink,
             grow = inputs$grow,
             fill = inputs$fill,
             clean = inputs$clean,
             tolerance = inputs$tolerance,
             n = inputs$roi_num) %>%
  roi_test(type = 'both')

rep_stats
```

```{r check_equal}
# check the app and the package output is equal
all.equal(round(stats$pcc, 2), round(c(rep_stats[[1]]$pcc, rep_stats[[2]]$pcc), 2))
all.equal(round(stats$moc, 2), round(c(rep_stats[[1]]$moc, rep_stats[[2]]$moc), 2))
```

## Description of the co-localization statistics

The following is a brief discussion of the theory and interpretations of different statistics that we used in this package as a measure of co-localization. For thorough and formal details, check this article by [Dunn et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3074624/).

### PCC

Pearson's correlation coefficient is the co-variance of the pixel intensity from the two channels. The mean of the intensities is subtracted from each pixel which makes the coefficient independent of the background level.

The PCC is calculated as follows

$$
PCC = \frac{\sum_i{(R_i-\bar R)\times(G_i-\bar G)}}{\sqrt{\sum_i{(R_i-\bar R)^2\times\sum_i(G_i-\bar G)^2}}}
$$

Where $R_i$ and the $G_i$ is the intensities of the red and green channels and the $\bar R$ and $\bar G$ are the average intensities.

The values of PCC are between 1 and -1 for perfect correlations in the positive and negative directions respectively and 0 means no correlation.

### MOC

On the other hand, the Manders Overlap Coefficient doesn't require subtraction of the mean. Therefore, the values are always between 0 and 1. Also, the MOC is independent from signal proportionality. 

$$
MOC = \frac{\sum_i{(R_i\times G_i)}}{\sqrt{\sum_i{R_i^2\times\sum_i G_i^2}}}
$$
Where $R_i$ and the $G_i$ is the intensities of the red and green channels.

## Colocalization Benchmark Source

> The Colocalization Benchmark Source (CBS) is a free collection of downloadable images to test and validate the degree of colocalization of markers in fluorescence microscopy studies. It consists of computer-simulated images with exactly known (pre-defined) values of colocalization ranging from 0% to 90%. They can be downloaded as sets as well as separately.

Here, we used two examples from the CBS to test the results of the analysis using `colocr`.

```{r example1_roi, fig.height=4, fig.width=4}
# load image
fl1 <- system.file('extdata', 'example1.png', package = 'colocr')
ex1 <- image_load(fl1)

# select and show regions of interest (based on the app)
par(mfrow=c(2,2), mar = rep(1, 4))

ex1 %>%
  roi_select(threshold = 90, shrink = 8, n = 50) %>%
  roi_show()
```

```{r example1_test}
# calculate co-localization statistics
# (expected pcc = 0.68, moc = 0.83)
ex1 %>%
  roi_select(threshold = 90, shrink = 8, n = 50) %>%
  roi_test(type = 'both') %>%
  colMeans()
```

```{r example3_roi, fig.height=4, fig.width=4}
# load image
fl3 <- system.file('extdata', 'example3.png', package = 'colocr')
ex3 <- image_load(fl3)

# select and show regioins of interest (based on the app)
par(mfrow=c(2,2), mar = rep(1, 4))

ex3 %>%
  roi_select(threshold = 90, shrink = 9, grow = 1, fill = 10, clean = 10, n = 20) %>%
  roi_show()
```

```{r example3_test}
# calculate co-localization statistics
# (expected pcc = 0.61, moc = 0.71)
ex3 %>%
  roi_select(threshold = 90, shrink = 9, grow = 1, fill = 10, clean = 10, n = 20) %>%
  roi_test() %>%
  colMeans()
```

## Acknowledgement

* The vignette images from [Lai Huyen Trang](https://www.researchgate.net/profile/Lai_Huyen_Trang)  
* The test examples from [Colocalization Benchmark Source (CBS)](https://www.colocalization-benchmark.com/index.html)  
* The implementation of the co-localization statistics from [Dunn et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3074624/)  
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{colocr_app}
\alias{colocr_app}
\title{Run the shiny App}
\usage{
colocr_app()
}
\description{
Run the shiny App
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{.intensity_get}
\alias{.intensity_get}
\title{Get pixel intensities}
\usage{
.intensity_get(img, ind = c(1, 2))
}
\arguments{
\item{img}{An object of class \code{\link[imager]{cimg}}}

\item{ind}{A \code{numeric} of length two for channel indexes}
}
\value{
A \code{list} of three items. The first two items are the values of
the pixel intensities of the channels indicated by \code{ind}. The third is
the labels of the individual regions of interest.
}
\description{
Get the pixel intensities of certain image channels
}
\examples{
# load image
fl <- system.file('extdata', 'Image0001_.jpg', package = 'colocr')
img <- image_load(fl)

# choose parameters
int <- roi_select(img, threshold = 90) \%>\%
  .intensity_get()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{.manders}
\alias{.manders}
\title{Calculate Marnders Overlap Coefficient}
\usage{
.manders(r, g)
}
\arguments{
\item{r}{A \code{numeric} vector}

\item{g}{A \code{numeric} vector}
}
\value{
A \code{numeric} of length one.
}
\description{
Calculates the manders overlap coefficient between two numeric vectors
}
\examples{
set.seed(123)
r <- rnorm(10)

set.seed(1234)
g <- rnorm(10)

.manders(r, g)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{roi_show}
\alias{roi_show}
\title{Show the selected regions of interest}
\usage{
roi_show(img, ind = c(1, 2))
}
\arguments{
\item{img}{A \code{\link[imager]{cimg}} object or a \code{list} of multiple
images such as the one returned from \code{\link{roi_select}}}

\item{ind}{A \code{numeric} object of length two. For the channel indexes.
or a \code{list} of similar vectors for each of \code{img} items.}
}
\description{
Show/highlight the selected regions of interest on different image channels
}
\details{
calling this function with \code{img} object which is returned from
\code{\link{roi_select}} returns four different plots. The original image, a
low resolution representation of the selected regions of interest and the
two channels indicated through \code{ind} highlighted.
}
\examples{
# load images
fl <- system.file('extdata', 'Image0001_.jpg', package = 'colocr')
img <- image_load(fl)

# choose and show ROI
oldpar <- par()
par(mfrow=c(2,2))

roi_select(img, threshold = 90) \%>\%
  roi_show()

par(oldpar)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{roi_test}
\alias{roi_test}
\title{Test Co-localization}
\usage{
roi_test(img, ind = c(1, 2), type = "pcc")
}
\arguments{
\item{img}{A \code{\link[imager]{cimg}} object or a \code{list} of multiple
images such as the one returned from \code{\link{roi_select}}}

\item{ind}{A \code{numeric} object of length two. For the channel indexes.
or a \code{list} of similar vectors for each of \code{img} items.}

\item{type}{A \code{character} vector of the desired co-localization
statistics. Default is 'pcc', other inputs are 'moc' or 'both'.}
}
\value{
A \code{data.frame} or a \code{list} of \code{data.frame}s.
}
\description{
Perform co-localization test statistics.
}
\details{
The co-localization stats requested in \code{type} is returned as
a column for each. When different labels are provided, the stats are
calculated for each label individually. When is \code{img} is a \code{list}
a \code{list} of such \code{data.frame}s is returned
}
\examples{
# load images
fl <- system.file('extdata', 'Image0001_.jpg', package = 'colocr')
img <- image_load(fl)

# choose roi and test colocalization
roi_select(img, threshold = 90) \%>\%
  roi_test()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{roi_select}
\alias{roi_select}
\title{Select regions of interest}
\usage{
roi_select(img, threshold, shrink = 5, grow = 5, fill = 5,
  clean = 5, tolerance = 0.1, n = 1)
}
\arguments{
\item{img}{An object of class \code{\link[imager]{cimg}} or a \code{list} of
multiple \code{\link[imager]{cimg}} items}

\item{threshold}{A \code{numeric} to be passed to
\code{\link[imager]{threshold}} or a \code{vector} of values for each
image in \code{img}}

\item{shrink}{A \code{numeric} to be passed to \code{\link[imager]{shrink}}
or a \code{vector} of values for each image in \code{img}}

\item{grow}{A \code{numeric} to be passed to \code{\link[imager]{grow}}
or a \code{vector} of values for each image in \code{img}}

\item{fill}{A \code{numeric} to be passed to \code{\link[imager]{fill}}
or a \code{vector} of values for each image in \code{img}}

\item{clean}{A \code{numeric} to be passed to \code{\link[imager]{clean}}
or a \code{vector} of values for each image in \code{img}}

\item{tolerance}{A \code{numeric} to be passed to \code{\link[imager]{label}}
or a \code{vector} of values for each image in \code{img}}

\item{n}{A \code{numeric} of the number of regions of interest
or a \code{vector} of values for each image in \code{img}}
}
\value{
A \code{\link[imager]{cimg}}. The original input \code{img} with an
additional attribute \code{label}. \code{label} is a \code{vector} of
\code{integer}s. The labels for the selected regions of interests starts
from 1 and 0 is ignored. When \code{img} is a list, a \code{list} is
returned.
}
\description{
Select regions of interest in an image using different morphological operations
}
\details{
The function applies several \code{\link{imager}} morphological
manipulations to select the regions of interest. These include
\code{\link[imager]{threshold}} which sets all values below certain cut to
0; \code{\link[imager]{shrink}}/\code{\link[imager]{grow}} for pixel set
dilation and erosion; \code{\link[imager]{fill}}/\code{\link[imager]{clean}}
for removing isolated regions and holes. When \code{n} is provided, the
individual regions (connected components) are selected where \code{tolerance}
is used to determine if two pixels belong to the same region.
}
\examples{
# load images
fl <- system.file('extdata', 'Image0001_.jpg', package = 'colocr')
img <- image_load(fl)

# choose ROI
newimg <- roi_select(img, threshold = 90)

# check the ROI labels
unique(attr(newimg, 'label'))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{.labels_add}
\alias{.labels_add}
\title{Label regions of interest}
\usage{
.labels_add(px, tolerance, n)
}
\arguments{
\item{px}{An object of class \code{\link[imager]{pixset}}}

\item{tolerance}{A \code{numeric} to be passed to \code{\link[imager]{label}}}

\item{n}{A \code{numeric}, the number of desired regions of interest}
}
\value{
An object of class \code{\link[imager]{cimg}}. The labels are coded
the values in the object starting from 1. The rest of the image is labeled 0.
}
\description{
Add labels to regions of interest in an image
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocr.R
\docType{package}
\name{colocr}
\alias{colocr}
\alias{colocr-package}
\title{\code{colocr}: Conduct Co-localization Analysis of Microscopy Images.}
\description{
Automate the co-localization analysis of fluorescence microscopy images.
Selecting regions of interest, extract pixel intensities from the image
channels and calculate different co-localization statistics.
}
\section{\code{colocr} functions}{

\code{\link{roi_select}}
\code{\link{roi_show}}
\code{\link{roi_check}}
\code{\link{roi_test}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{image_load}
\alias{image_load}
\title{Load images from files}
\usage{
image_load(image_file)
}
\arguments{
\item{image_file}{A \code{character} vector of one or more paths to image
files}
}
\value{
A \code{cimg} object or a \code{list} of \code{cimg} objects when
multiple files are passed to \code{image_file}.
}
\description{
A wrap around \code{\link[magick]{image_read}} and
\code{\link[imager]{magick2cimg}} to load one or more images from files.
}
\examples{
# load image
fl <- system.file('extdata', 'Image0001_.jpg', package = 'colocr')
img <- image_load(fl)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{.pearson}
\alias{.pearson}
\title{Calculate Pearson's Correlation Coefficient}
\usage{
.pearson(r, g)
}
\arguments{
\item{r}{A \code{numeric} vector}

\item{g}{A \code{numeric} vector}
}
\value{
A \code{numeric} of length one.
}
\description{
Calculates the Pearson's correlation coefficient between two numeric vectors
}
\examples{
set.seed(123)
r <- rnorm(10)

set.seed(1234)
g <- rnorm(10)

.pearson(r, g)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{roi_check}
\alias{roi_check}
\title{Show pixel intensities}
\usage{
roi_check(img, ind = c(1, 2))
}
\arguments{
\item{img}{A \code{\link[imager]{cimg}} object or a \code{list} of multiple
images such as the one returned from \code{\link{roi_select}}}

\item{ind}{A \code{numeric} object of length two. For the channel indexes.
or a \code{list} of similar vectors for each of \code{img} items.}
}
\description{
Show the pixel intensities of certain image channels
}
\details{
Calling this function returns two plots. The first is a scatter
plot of the pixel intensities from two channels. The second is the density
distribution of the intensities from the two channels.
}
\examples{
# load images
fl <- system.file('extdata', 'Image0001_.jpg', package = 'colocr')
img <- image_load(fl)

# choose ROI and show the pixel intensities
oldpar <- par()
par(mfrow = c(1, 2))

roi_select(img, threshold = 90) \%>\%
  roi_check()

par(oldpar)

}
