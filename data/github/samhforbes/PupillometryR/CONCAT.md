
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PupillometryR

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/PupillometryR)](https://CRAN.R-project.org/package=PupillometryR)
[![Travis build
status](https://app.travis-ci.com/samhforbes/PupillometryR.svg?branch=master)](https://app.travis-ci.com/samhforbes/PupillometryR)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02285/status.svg)](https://doi.org/10.21105/joss.02285)

<!-- badges: end -->

The goal of PupillometryR is to to pre-process and then analyze simple
pupil experiments in R.

## Installation

You can install the released version of PupillometryR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("PupillometryR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("samhforbes/PupillometryR")
```

## Setup

This package (and the example dataset) was designed in part, based on
Sylvain Sirois’ MATLAB tutorial, [which can be found
here](https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=).

The intention is an integrated pipeline for pupillometric experiments,
from data cleaning, pre-processing, various analysis techniques, and
visualising results.

To use all the functionality and plots that follow from the
PupillometryR pipeline, please start with *make\_pupillometryr\_data*,
e.g.:

``` r
library(PupillometryR)
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: ggplot2
#> Loading required package: rlang

data("pupil_data")

#Check that IDs are not numeric
pupil_data$ID <- as.character(pupil_data$ID)
#remove participant number 8, who had problematic data
pupil_data <- subset(pupil_data, ID != 8)
#blinks were registered as -1, so replace with NAs
pupil_data$LPupil[pupil_data$LPupil == -1] <- NA
pupil_data$RPupil[pupil_data$RPupil == -1] <- NA

Sdata <- make_pupillometryr_data(data = pupil_data,
                                 subject = ID,
                                 trial = Trial,
                                 time = Time,
                                 condition = Type)
```

All further functions associated with the package follow from there. For
example:

``` r
plot(Sdata, pupil = LPupil, group = 'condition')
#> Warning: Removed 3639 rows containing non-finite values (stat_summary).
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

To follow a detailed walkthrough, run:

``` r
vignette('PupillometryR')
```

or head to
[samforbes.me/PupillometryR](http://samforbes.me/PupillometryR/)

## Getting help

Please use the issues tab
(<https://github.com/samhforbes/PupillometryR/issues>) to file any bugs
or suggestions. For general pupillometry information, I recommend
[Sylvain’s
website](https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=),
[as well as Jackson and Sirois
(2009)](https://doi.org/10.1111/j.1467-7687.2008.00805.x). For reading
about using GAMs in pupillometry [this paper by van Rij et al. is
excellent](https://journals.sagepub.com/doi/10.1177/2331216519832483),
for general GAMs knowledge I recommend [this tutorial by Michael
Clark](https://m-clark.github.io/generalized-additive-models/case_for_gam.html)
as well as the mgcv documentation, and for general FDA information [this
website is
helpful](https://www.psych.mcgill.ca/misc/fda/resources.html), along
with the Ramsay and Silverman book (1997). Additionally, check out the
[raincloud plots paper by Allen et
al.](https://wellcomeopenresearch.org/articles/4-63#:~:text=In%20essence%2C%20raincloud%20plots%20combine,error%2C%20such%20as%20a%20boxplot.),
which is used for some of the in-built plotting in this package.

## Citation

Please cite the JOSS paper for this package if you use it: Forbes, S.
(2020). PupillometryR: An R package for preparing and analysing
pupillometry data. Journal of Open Source Software, 5(50), 2285.
<https://doi.org/10.21105/joss.02285>

## Acknowledgements

This package has had suggestions, encouragement, and help from a number
of people, but I wish to especially highlight Sylvain Sirois and Mihaela
Duta, whose input has been instrumental. I’d also like to thank Jacolien
van Rij for her input with the GAMs modelling portion of this tutorial,
and TJ Mahr for contributing to extending the use of GAMs in the
vignette.

## References

\[1\] Jackson, I., & Sirois, S. (2009). Infant cognition: Going full
factorial with pupil dilation. *Developmental Science*, 12(4), 670-679.
<https://doi.org/10.1111/j.1467-7687.2008.00805.x>

\[2\] Allen, M., Poggiali, D., Whitaker, K., Marshall, T. R., & Kievit,
R. (2019). Raincloud plots: a multi-platform tool for robust data
visualization. *Wellcome Open Research*, 4, 1-41.
<https://doi.org/10.12688/wellcomeopenres.15191.1>

\[3\] Ramsay, J.O., & Silverman, B.W. (1997). *Functional data
analysis*. New York: Springer-Verlag.

\[4\] van Rij, J., Hendriks, P., van Rijn, H., Baayen, R. H., & Wood, S.
N. (2019). Analyzing the time course of pupillometric data. *Trends in
Hearing*, 23, 233121651983248.
<https://doi.org/10.1177/2331216519832483>
# News

## Update Sept 2021 PupillometryR 0.0.4

* Updated to remove warnings from new version of dplyr
* More robust to missing data in a single eye

## Update June 2020 PupillometryR 0.0.3

* Documentation moved to vignettes
* Dependency on package *spectral* removed due to problems with CRAN

## Update May 2020 PupillometryR 0.0.2

* Minor patch to avoid breaks with new version of dplyr
* Warning errors removed from ggplot2 when plotting
---
title: 'PupillometryR: An R package for preparing and analysing pupillometry data'
tags:
  - R
  - Pupillometry
  - Eye-tracking
authors:
  - name: Samuel H. Forbes
    orcid: 0000-0003-1022-4676
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: School of Psychology, University of East Anglia
   index: 1
date: 1 April 2020 #update this
bibliography: paper.bib
---

# Summary

The study of pupil dilation, or pupillometry, has undergone a considerable revival in recent years.
Early research noted substantive changes in pupil dilation due to mental effort [@Hess1964; @Kahneman1969; @Kahneman1973]; and this method has seen a surge in usage following a rise in the popularity of automatic eye-tracking technology coupled with the search for more robust methods in psychological science [@vanderWel2018; @Winn2018]. Pupillometry, measured using automated eye-tracking technology, thus provides a methodology that is cheaper than most neuroimaging techniques, but more powerful than many behavioural measures, such as overall looking time measurements common to eye-tracking experiments.

Despite these considerable advantages when compared to other paradigms, the details of the actual methods and pipelines involved with pupillometry have been a topic for much debate [@Mathot2018; @vanRij2019].
These methodological differences are mostly attributable to two causes.
The first is that the nature of pupillometric studies require careful design and procedure.
Pupil dilation is highly sensitive to not just luminance, but also other properties of the experiment, including gaze position due to stimuli locations [@Gagl2011], as well as stimuli-specific considerations such as brightness or location on the screen, and participant-specific considerations, such as the stress levels and memory load.
In this way, the very nature of any pupillometric study demands careful planning on design, and thoughtful consideration about the details of the methods.
The second such cause of the variable methods in pupillometric studies stems from the analytical decisions made in the processing of the data.
Due to many researchers following either closed, in-house analytical pipelines, or open-source pipelines that require costly software, or particular equipment at some point in the analysis, there has been a lack of consensus on the analytical decisions to be made in processing.
These decisions include how to filter and clean the data [@Jackson2009], decisions to be made on baselining the data [@Mathot2018], and whether to take time-course or aggregate pupil size forward as a variable [@vanRij2019].
Thus `PupillometryR` aims to assist experimenters by providing a clear pipeline which is both free and open-source, and available for usage with most common brands of eye-tracker.

# Implementation

![An example of a raincloud plot of pupil data in PupillometryR.](Raincloud.pdf)

`PupillometryR` extends a suggested pipeline for implementation of pupillometric studies [@Jackson2009], making heavy use of the `signal` package [@signal2014] for pre-processing of the pupil data. The in-built plotting functions, designed for ease of use, rely on `ggplot2` [@Wickham2016] for the data visualisation, and include raincloud plots [@Allen2019] as an in-built data visualisation option. For analysis several options are given, including the use of Generalised Additive Models which uses the `mgcv` package [@Wood2017], and Functional Data Analysis, which uses functions from the `fda` package [@Ramsay2020].

![An example of a functional t-test in PupillometryR.](FDA.pdf)

Some comprehensive pupillometry pipelines exist already in MATLAB [@Hershman2019; @SiroisPupillometryWalkthrough], and in R [@Geller2020], and while each of these have their own merits, none offer the start-to-finish comprehensive pipeline in an open-source language compatible with most brands of eye-tracker, which includes in-built plotting functions, and flexibility of analysis style (time windows, FDA, GAMs) offered in `PupillometryR`. In addition, `PupillometryR` is available on CRAN, making it easy to download for R users, and subject to regular CRAN checks.

# Acknowledgements

This package benefitted greatly from discussions with and past work made available by Sylvain Sirois and advice from Mihaela Duta and Jacolien van Rij. I am also grateful to David Robinson for assistance with the implementation of raincloud plots in the package, and Micah Allen and the raincloud plots team for making the raincloud plot technique available.

# References
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

# PupillometryR

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/PupillometryR)](https://CRAN.R-project.org/package=PupillometryR)
[![Travis build status](https://app.travis-ci.com/samhforbes/PupillometryR.svg?branch=master)](https://app.travis-ci.com/samhforbes/PupillometryR)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02285/status.svg)](https://doi.org/10.21105/joss.02285)

<!-- badges: end -->

The goal of PupillometryR is to to pre-process and then analyze simple pupil
experiments in R.

## Installation

You can install the released version of PupillometryR from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("PupillometryR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("samhforbes/PupillometryR")
```
## Setup

This package (and the example dataset) was designed in part, based on Sylvain Sirois' MATLAB tutorial, [which can be found here](https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=).

The intention is an integrated pipeline for pupillometric experiments, from data cleaning, pre-processing, various analysis techniques, and visualising results.

To use all the functionality and plots that follow from the PupillometryR pipeline, please start with *make_pupillometryr_data*, e.g.:

```{r example}
library(PupillometryR)

data("pupil_data")

#Check that IDs are not numeric
pupil_data$ID <- as.character(pupil_data$ID)
#remove participant number 8, who had problematic data
pupil_data <- subset(pupil_data, ID != 8)
#blinks were registered as -1, so replace with NAs
pupil_data$LPupil[pupil_data$LPupil == -1] <- NA
pupil_data$RPupil[pupil_data$RPupil == -1] <- NA

Sdata <- make_pupillometryr_data(data = pupil_data,
                                 subject = ID,
                                 trial = Trial,
                                 time = Time,
                                 condition = Type)
```

All further functions associated with the package follow from there. 
For example:

```{r}
plot(Sdata, pupil = LPupil, group = 'condition')
```

To follow a detailed walkthrough, run:

```{r eval = F}
vignette('PupillometryR')
```

or head to [samforbes.me/PupillometryR](http://samforbes.me/PupillometryR/)

## Getting help

Please use the issues tab (https://github.com/samhforbes/PupillometryR/issues) to file any bugs or suggestions.
For general pupillometry information, I recommend [Sylvain's website](https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=), [as well as Jackson and Sirois (2009)](https://doi.org/10.1111/j.1467-7687.2008.00805.x). For reading about using GAMs in pupillometry [this paper by van Rij et al. is excellent](https://journals.sagepub.com/doi/10.1177/2331216519832483), for general GAMs knowledge I recommend [this tutorial by Michael Clark](https://m-clark.github.io/generalized-additive-models/case_for_gam.html) as well as the mgcv documentation, and for general FDA information [this website is helpful](https://www.psych.mcgill.ca/misc/fda/resources.html), along with the Ramsay and Silverman book (1997). Additionally, check out the [raincloud plots paper by Allen et al.](https://wellcomeopenresearch.org/articles/4-63#:~:text=In%20essence%2C%20raincloud%20plots%20combine,error%2C%20such%20as%20a%20boxplot.), which is used for some of the in-built plotting in this package.

## Citation

Please cite the JOSS paper for this package if you use it:
Forbes, S. (2020). PupillometryR: An R package for preparing and analysing pupillometry data. Journal of Open Source Software, 5(50), 2285. https://doi.org/10.21105/joss.02285

## Acknowledgements

This package has had suggestions, encouragement, and help from a number of people, but I wish to especially highlight Sylvain Sirois and Mihaela Duta, whose input has been instrumental. I'd also like to thank Jacolien van Rij for her input with the GAMs modelling portion of this tutorial, and TJ Mahr for contributing to extending the use of GAMs in the vignette.

## References

[1] Jackson, I., & Sirois, S. (2009). Infant cognition: Going full factorial with pupil dilation. *Developmental Science*, 12(4), 670-679. https://doi.org/10.1111/j.1467-7687.2008.00805.x

[2] Allen, M., Poggiali, D., Whitaker, K., Marshall, T. R., & Kievit, R. (2019). Raincloud plots: a multi-platform tool for robust data visualization. *Wellcome Open Research*, 4, 1-41.
https://doi.org/10.12688/wellcomeopenres.15191.1

[3] Ramsay, J.O., & Silverman, B.W. (1997). *Functional data analysis*. New York: Springer-Verlag.

[4] van Rij, J., Hendriks, P., van Rijn, H., Baayen, R. H., & Wood, S. N. (2019). Analyzing the time course of pupillometric data. *Trends in Hearing*, 23, 233121651983248. https://doi.org/10.1177/2331216519832483
---
title: "PupillometryR"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{PupillometryR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(PupillometryR)
```

This package is designed to make dealing with pupil data (perhaps more traditionally done in MATLAB) easier to wrangle in R. It makes heavy use of a few packages which should be acknowledged here, especially the (excellent) packages fda and signal.

As well as the above packages, it is very important to note that the type of analysis shown here has been available in MATLAB for a while, and there is an excellent tutorial on it, which I thoroughly recommend reading first, written by Sylvain Sirois, [here:](https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=).

It's worth making sure that your setup and experiment do facilitate the use of pupillometry - it may not be suited for all kinds of experiment.

## Getting started

We will first run through an example analysis with the data provided in the package, which, again, comes from Sylvain Sirois' tutorial on [his webpage](https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=).

The first thing I would recommend doing is having a close look at the pupil data. Eyetrackers have a couple of different ways of dealing with this, so it's important to know a few things:

1. What unit of measurement is being used
2. What value is given to missing data or blinks (typically ., -1, or *NA*)
3. What framerate were you recording at, and is this consistent.

The data here is an eyetracking experiment with hard and easy trials, performed on 9 participants. Participant 8 needs to be removed (more details on Sylvain's tutorial).

It's important here to make sure that the missing data get set to *NA* and not to a numeric value (such as -1) or deleted.

```{r}
data("pupil_data")

#Check that IDs are not numeric
pupil_data$ID <- as.character(pupil_data$ID)
#remove participant number 8, who had problematic data
pupil_data <- subset(pupil_data, ID != 8)
#blinks were registered as -1, so replace with NAs
pupil_data$LPupil[pupil_data$LPupil == -1] <- NA
pupil_data$RPupil[pupil_data$RPupil == -1] <- NA

```

The plotting is also a theme of this tutorial, so I will set a nice theme that makes the plots look pretty:

```{r}
library(ggplot2)
theme_set(theme_classic(base_size = 12))
```

First up, we need to put the data into pupillometryR format for further analysis

```{r}
Sdata <- make_pupillometryr_data(data = pupil_data,
                                 subject = ID,
                                 trial = Trial,
                                 time = Time,
                                 condition = Type)
```

In the current data, it is not a concern, but there may be a situation where certain timebins are missing from your data. This can be fixed here, and we will look at the raw data:

```{r}
new_data <- replace_missing_data(data = Sdata)

head(new_data)
```

Equally, if your data is not cut to the time windows that you are interested in, the subset_data function allows trimming. PupillometryR has some built-in plotting functions to allow you to look at certain data types. You simply need to specify a pupil column to display, and how you would like the data displayed in groups. The plots are ggplot items, so can be edited with themes, and arguments such as ylab(). Below we display it first by condition, then by subject:

```{r}
plot(new_data, pupil = LPupil, group = 'condition')

plot(new_data, pupil = LPupil, group = 'subject') 
```

## Smoothing and cleanup

PupillometryR offers a few smoothing options to make processing the data a little easier. We'll do the full set here.  A great reference for these is Sylvain's tutorial, and also Jackson and Sirois, 2009, Dev. Sci. First off, we can regress one pupil against the other to get some measure of smoothing.

```{r}
regressed_data <- regress_data(data = new_data,
                               pupil1 = RPupil,
                               pupil2 = LPupil)

```

Now that we have done that, we want the mean of the two pupil sizes, so let's see how that looks:

```{r}
mean_data <- calculate_mean_pupil_size(data = regressed_data, 
                                       pupil1 = RPupil, 
                                       pupil2 = LPupil)

plot(mean_data, pupil = mean_pupil, group = 'subject')
```

Now that we have a single pupil measure to work with, we can manipulate the data a little more easily. First thing we can do is to downsample the data. This is useful when we have large scale data, or when we have sampled at a really high rate, and we need to reduce it so we are measuring meaningful change. Here we have taken the median, but he mean could also be taken. We just need to specify the timebin size, in ms:

```{r}
mean_data <- downsample_time_data(data = mean_data,
                              pupil = mean_pupil,
                              timebin_size = 50,
                              option = 'median')
```

Now we need to clean up our data - let's first assess how much missing data there is:

```{r}
missing <- calculate_missing_data(mean_data, 
                                  mean_pupil)
head(missing)
```

We can see if we view the whole file that participant 6 has a fair amount of trials with a high missing data proportion. Now we need to clean this up. We have two parameters to do this - first is the proportion of *data* that we can accept missing in one trial before removing it from analysis. The second is what proportion of *trials* we can accept as missing before removing a participant for being unreliable. In this example, we will remove trials that have more than 75% of data missing, and we will remove participants that have more than 75% of trials removed.

```{r message = T}
mean_data2 <- clean_missing_data(mean_data,
                                 pupil = mean_pupil,
                                 trial_threshold = .75,
                                 subject_trial_threshold = .75)
```

Now we come to filtering the data. PupillometryR offers 3 filter types: A hanning window, a low-pass butterworth filter, and a median filter. The low-pass filter can be a little unstable at the beginning and end of each trial, so it's worth looking at your data to see if it's appropriate. Here we will use the median filter. The degree gives the size of the rolling window.

```{r}
filtered_data <- filter_data(data = mean_data2,
                             pupil = mean_pupil,
                             filter = 'median',
                             degree = 11)

plot(filtered_data, pupil = mean_pupil, group = 'subject')
```

The next step is to interpolate across blinks. Filtering before the interpolation generally allows more sensible interpolation, but the way this is done varies a bit on the data, and you will see examples without this. We can interpolate in this package either linear or cubic, but again, it's best to always check your data afterwards to make sure it looks the way you might expect. Here we opt for the linear interpolation:

```{r}
int_data <- interpolate_data(data = filtered_data,
                             pupil = mean_pupil,
                             type = 'linear')

plot(int_data, pupil = mean_pupil, group = 'subject')
```

Baselining the data is a powerful way of making sure we control for between-participant variance of average pupil size. If we are looking at analyses that are largely within-subject, as we do here, this may not be such an issue, but we will do this anyway. This function allows us to baseline to the mean pupil size within a time window. Here we are just taking the first 100ms of the trial. If your baseline period is just outside of your analysis window (which it often will be), you can use subset_data() to remove that after baselining.

```{r}
base_data <- baseline_data(data = int_data,
                           pupil = mean_pupil,
                           start = 0,
                           stop = 100)

plot(base_data, pupil = mean_pupil, group = 'subject')
```

## Window analyses

PupillometryR gives us a couple of options for window analysis. One is overall averages, the other is to break the data up into discrete time windows, and to analyse them. First we will opt for the overall averages.
We can plot these with any of boxplots, violin plots, or, since the new edition, Micah Allen-esque raincloud plots (Allen et al., 2018).

```{r}
window <- create_window_data(data = base_data,
                             pupil = mean_pupil)

plot(window, pupil = mean_pupil, windows = F, geom = 'boxplot')

head(window)
```

We could then simply analyse this with a t-test if we wished.

```{r}
t.test(mean_pupil ~ Type, paired = T, data = window)
```

Alternatively, we may wish to look at the data in chunks. Here we group the data in to 2000ms timebins for analysis (and we will opt for the raincloud plot in this instance):

```{r}
timeslots <- create_time_windows(data = base_data,
                                 pupil = mean_pupil,
                                 breaks = c(0, 2000, 4000, 6000, 8000, 10000))

plot(timeslots, pupil = mean_pupil, windows = T, geom = 'raincloud')

head(timeslots)
```

And again, we could analyse this with a linear model or an anova:

```{r}
lm(mean_pupil ~ Window * Type, data = timeslots)
```

## Modelling with Generalised Additive Models

Here we interfact with the mgcv package, an exceptionally powerful package for GAM data, by Simon Wood. I strongly encourage reading the vignette and checking out some of the great online tutorials (of which there are plenty; I quite like Michael Clark's one [here](https://m-clark.github.io/generalized-additive-models/case_for_gam.html)) before proceeding with these. 

We have to do a little bit of setting up of our variables (scaling and centering) before we continue. I need to make some variables numeric (the ones with an n on the end), and I am using the way trials are labelled to make this a numeric variable (this would probably be different for your data).

```{r}
library(mgcv)

base_data$IDn <- as.numeric(base_data$ID)
base_data$Typen <- ifelse(base_data$Type == 'Easy', .5, -.5)
base_data$Trialn <- as.numeric(substr(base_data$Trial, 5, 5))
base_data$Trialn <- ifelse(base_data$Typen == .5, base_data$Trialn, base_data$Trialn + 3)
base_data$ID <- as.factor(base_data$ID)
base_data$Trial <- as.factor(base_data$Trial)
```

Right, let's proceed with setting up a simple model. It's recommended for the amount of data points we might have for PupillometryR, bams might be a better option, but both gam() and bam() will work.

```{r}
m1 <- bam(mean_pupil ~ s(Time) +
            s(Time,  by = Typen),
          data = base_data,
          family = gaussian)

summary(m1)
```

We can use our default plotting function to see how it looks compared to the raw data, just by specifying the model= argument.

```{r}
plot(base_data, pupil = mean_pupil, group = 'condition', model = m1)
```

Of course there is the fact that we expect there to by some variation by trial, and that we should expect there to be differences for each participant. Our model also only accounts for a small amount of the variance. This model, therefore is no good. The way to check this is to assess model fit with the qqnorm, and to check the autocorrelation. We can do this with the help of the itsadug package:

```{r}
qqnorm(resid(m1))

itsadug::acf_resid(m1)
```

While the qqnorm looks to be almost passable, the autocorrelation in the second plot is very high. This is an important consideration in time-series data, and due consideration needs to be given to this. For a full discussion of how this issue affects time course data, and specifically pupil data, I highly recommend [Jacolien van Rij et al's paper here](https://journals.sagepub.com/doi/10.1177/2331216519832483). 

To reduce autocorrelation there are many methods we can try, as you will see from the above paper. I will stop short of repeating each of the steps taken in the excellent paper above, and jump straight away to a much more appropriate model for this data. First I will code in events (participant per trial). I will also create a second data frame to do this (model_data), so that the data we are working with doesn't lose its properties, and we can keep using the plotting functions.

```{r}
base_data$Event <- interaction(base_data$ID, base_data$Trial, drop = T)

model_data <- base_data
model_data <- itsadug::start_event(model_data,
                          column = 'Time', event = 'Event')

model_data <- droplevels(model_data[order(model_data$ID,
                                          model_data$Trial,
                                          model_data$Time),])
```

We now need to model this. We are setting an AR parameter, and allowing events to vary by time. You will see our deviance accounted for is now up around 96% - much better! The qqnorm is still far from perfect, and the commented-out model below m2 would do a bit better at this (again from van Rij et al) by using a scaled t distribution - but would take ages to run. 

```{r}
m2 <- bam(mean_pupil ~ Typen +
            s(Time,  by = Typen) +
            s(Time, Event, bs = 'fs', m = 1),
          data = base_data,
          family = gaussian,
          discrete = T,
          AR.start = model_data$start.event, rho = .6)

# m2 <- bam(mean_pupil ~ 
          #   s(Time,  by = Typen) +
          #   s(Time, Event, bs = 'fs', m = 1),
          # data = base_data,
          # family = scat,
          # discrete = T,
          # AR.start = model_data$start.event, rho = .6)

summary(m2)
qqnorm(resid(m2))
itsadug::acf_resid(m2)
plot(base_data, pupil = mean_pupil, group = 'condition', model = m2)
```

The summary from our second model indicates that there may be marginal evidence for this effect of condition. But how and when do they diverge???

(In fact, TJ Mahr was good enough to point out this elegant solution for this using GAM methods with the itsadug package, which I will [link to](https://gist.github.com/tjmahr/0d2b41ea1525205a99b19770fc916a90) rather than take credit for)

## Estimating divergences with functional data analysis

The above analyses may well suffice for what we have planned. However, sometimes it's useful for analysis to examine change over time, especially how and when two conditions diverge, and we can do that with Functional Data Analysis (FDA). This part of the package makes usage of the fda package. The complete guide really has been written in 1997 by Ramsay and Silverman, and there is a very helpful website on FDA [here](https://www.psych.mcgill.ca/misc/fda/resources.html). This package is currently only setup to use this analysis for two-condition experiments, but I hope to add options for functional ANOVA in the future.

To do this, first we want get the difference between the two conditions for each participant. By default this package wil take condition 2 - condition 1, so reorder the factors if required.

```{r}
differences <- create_difference_data(data = base_data,
                                      pupil = mean_pupil)

plot(differences, pupil = mean_pupil, geom = 'line')
```

We can now convert this to a functional data structure, made up of curves. To do this for this data we are going to make it up of cubics (order = 4) with 10 knots (basis = 10). The appropriate numbers here will depend on your dataset, and I strongly advise consulting Ramsay and Silverman's book, and the FDA website, as well as Sylvain's paper mentioned above. This interfaces with the fda package.

```{r}
spline_data <- create_functional_data(data = differences,
                                      pupil = mean_pupil,
                                      basis = 10,
                                      order = 4)


plot(spline_data, pupil = mean_pupil, geom = 'line', colour = 'blue')
```

That looks like it's done a pretty good job capturing the data. The advantage of this kind of analysis is that we can treat each curve as a function, and run a single functional t-test to work out during which window there are divergences. This package allows us to do that directly, and to observe the results.

```{r}
ft_data <- run_functional_t_test(data = spline_data,
                                 pupil = mean_pupil,
                                 alpha = 0.05)


plot(ft_data, show_divergence = T, colour = 'red', fill = 'orange')
```

If show_divergence is set to TRUE, the plot will highlight where the two conditions diverge at the alpha you set.

*NB* Remember the above discussion on autocorrelation in the GAMMs portion of this walkthrough? We are still dealing with time-series data, so this hasn't necessarily gone away. I am working on adding more powerful FDA techniques into this package to deal with these issues, so please watch this space.

## Acknowledgements

This package has had suggestions, encouragement, and help from a number of people, but I wish to especially highlight Sylvain Sirois and Mihaela Duta, whose input has been instrumental. I'd also like to thank Jacolien van Rij for her input with the GAMMs modelling portion of this tutorial.

## References

[1] Jackson, I., & Sirois, S. (2009). Infant cognition: Going full factorial with pupil dilation. *Developmental Science*, 12(4), 670-679. https://doi.org/10.1111/j.1467-7687.2008.00805.x

[2] Allen, M., Poggiali, D., Whitaker, K., Marshall, T. R., & Kievit, R. (2019). Raincloud plots: a multi-platform tool for robust data visualization. *Wellcome Open Research*, 4, 1-41.
https://doi.org/10.12688/wellcomeopenres.15191.1

[3] Ramsay, J.O., & Silverman, B.W. (1997). *Functional data analysis*. New York: Springer-Verlag.

[4] van Rij, J., Hendriks, P., van Rijn, H., Baayen, R. H., & Wood, S. N. (2019). Analyzing the time course of pupillometric data. Trends in Hearing, 23, 233121651983248. https://doi.org/10.1177/2331216519832483
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_smoothed_data.R
\name{regress_data}
\alias{regress_data}
\title{Regress one pupil against another for extra smoothing}
\usage{
regress_data(data, pupil1, pupil2)
}
\arguments{
\item{data}{a PupillometryR dataframe}

\item{pupil1}{Column name for first pupil data}

\item{pupil2}{Column name for second pupil data}
}
\value{
a PupillometryR dataframe with smoothed pupil values
}
\description{
regress_data runs a simple linear regression of pupil1 against pupil2 and the reverse.
This can help to account for small amount of bumpiness in the data.
The regression runs over each participant and each trial, per time.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
subject = ID,
trial = Trial,
time = Time,
condition = Type)
regressed_data <- regress_data(data = Sdata,
pupil1 = RPupil,
pupil2 = LPupil)
mean_data <- calculate_mean_pupil_size(data = regressed_data,
pupil1 = RPupil, pupil2 = LPupil)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset_data.R
\name{mean2}
\alias{mean2}
\title{Helper function mean2}
\usage{
mean2(x)
}
\arguments{
\item{x}{the object}
}
\description{
Somewhat useful function for ignoring NAs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_and_downsample.R
\name{downsample_time_data}
\alias{downsample_time_data}
\title{Downsample frequency to reduce number of samples and data size}
\usage{
downsample_time_data(data, pupil, timebin_size, option = c("mean", "median"))
}
\arguments{
\item{data}{your data of class PupillometryR}

\item{pupil}{a column name denoting pupil size}

\item{timebin_size}{the size of the new timebin you wish to use}

\item{option}{what should be calculated in each timebin - mean or median. Defaults to mean.}
}
\value{
A downsampled dataframe of class PupillometryR
}
\description{
This function is useful if you were sampling at a very high frequency (eg 500Hz)
causing the data size to be hard to manage, and high autocorrelation.
Careful decisions should be made about the time bin size and appropriateness
of this function, with respect to the data type.
}
\examples{
data(pupil_data)
Sdata <- make_pupillometryr_data(data = pupil_data,
subject = ID,
trial = Trial,
time = Time,
condition = Type)
new_data <- downsample_time_data(data = Sdata,
pupil = LPupil,
timebin_size = 50,
option = 'mean')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pupil_data.R
\name{plot.PupillometryR}
\alias{plot.PupillometryR}
\title{Pre-prepared plots of PupillometryR data}
\usage{
\method{plot}{PupillometryR}(
  x,
  pupil,
  group = c("none", "condition", "subject"),
  geom = c("point", "line", "pointrange"),
  model = NULL,
  ...
)
}
\arguments{
\item{x}{A PupillometryR dataframe}

\item{pupil}{Column name of pupil data to be plotted}

\item{group}{What to group the data by (none, condition, or subject)}

\item{geom}{Geom to pass to ggplot. Either point, line, or pointrange.}

\item{model}{Optional argument to plot agains a fitted model}

\item{...}{Ignored}
}
\value{
A ggplot object
}
\description{
The plot functions are designed to run with just data and pupil selections,
with some additional options for fun with plotting. This allows to see
raw data as points, grouped by either subject or condition.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
subject = ID,
trial = Trial,
time = Time,
condition = Type)
Sdata2 <- downsample_time_data(data = Sdata,
pupil = LPupil,
timebin_size = 100,
option = 'median')
p <- plot(Sdata2, pupil = LPupil, group = 'subject')
p

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_smoothed_data.R
\name{filter_data}
\alias{filter_data}
\title{Run a filter on the data to smooth it out.}
\usage{
filter_data(
  data,
  pupil,
  filter = c("median", "hanning", "lowpass"),
  degree = 11
)
}
\arguments{
\item{data}{a PupillometryR dataframe}

\item{pupil}{column name for pupil data}

\item{filter}{option for filtering the data}

\item{degree}{filter degree}
}
\value{
filtered pupil data
}
\description{
filter_data allows three different options for filtering, a butterworth lowpass filter, a hanning filter, or
a median filter. You can also set the degree of this filter; we recommend a default of 11.
This filters on one pupil, it can be re-run on a second pupil if needed. Lowpass makes use of the
butterworth filter and filtfilt from package signal, median makes use of runmed.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
subject = ID,
trial = Trial,
time = Time,
condition = Type)
mean_data <- calculate_mean_pupil_size(data = Sdata,
pupil1 = RPupil, pupil2 = LPupil)
filtered_data <- filter_data(data = mean_data,
pupil = mean_pupil,
filter = 'hanning',
degree = 11)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_and_downsample.R
\name{clean_missing_data}
\alias{clean_missing_data}
\title{Clean missing data above an acceptable threshold}
\usage{
clean_missing_data(
  data,
  pupil,
  trial_threshold = 1,
  subject_trial_threshold = 1
)
}
\arguments{
\item{data}{your data of class PupillometryR}

\item{pupil}{a column name denoting pupil size}

\item{trial_threshold}{a proportion of missing data over which a trial can be considered lost}

\item{subject_trial_threshold}{a proportion of missing trials over which a participant can be considered lost.}
}
\value{
A cleaned PupillometryR dataframe
}
\description{
This function can be used to remove trials and participants
who do not meet the threshold for a study. Note that there are two parameters for
cleaning, one to remove trials above a threshold,
the second to remove participants who drop more than a certain amount of trials.
}
\examples{
data(pupil_data)
Sdata <- make_pupillometryr_data(data = pupil_data,
subject = ID,
trial = Trial,
time = Time,
condition = Type)
new_data <- downsample_time_data(data = Sdata,
pupil = LPupil,
timebin_size = 50,
option = 'mean')
calculate_missing_data(data = new_data, pupil = LPupil)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset_data.R
\name{create_time_windows}
\alias{create_time_windows}
\title{Make PupillometryR dataframe into multiple time windows for easy analysis}
\usage{
create_time_windows(data, pupil, breaks)
}
\arguments{
\item{data}{a PupillometryR dataframe}

\item{pupil}{column name denoting pupil data to be used}

\item{breaks}{a vector or numbers indicating start times for each window}
}
\value{
a Pupil_window_data dataframe
}
\description{
This function creates a single collapsed data frame for easy analysis with an anova or model,
per condition.
By comparison create_window_data allows collapsing all into a single time window.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
regressed_data <- regress_data(data = Sdata, pupil1 = RPupil, pupil2 = LPupil)
mean_data <- calculate_mean_pupil_size(data = regressed_data,
pupil1 = RPupil, pupil2 = LPupil)
base_data <- baseline_data(data = mean_data, pupil = mean_pupil, start = 0, stop = 100)
time_window <- create_time_windows(data = base_data, pupil = mean_pupil,
breaks = c(1000, 2000))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pupil_data.R
\name{plot.Pupil_window_data}
\alias{plot.Pupil_window_data}
\title{Pre-prepared plots of PupillometryR data}
\usage{
\method{plot}{Pupil_window_data}(
  x,
  pupil,
  windows = c(FALSE, TRUE),
  geom = c("raincloud", "violin", "boxplot"),
  ...
)
}
\arguments{
\item{x}{A Pupil_window_data dataframe}

\item{pupil}{Column name of pupil data to be plotted}

\item{windows}{Whether you want to include time windows in the plot - logical}

\item{geom}{violin plots or boxplots. The newest version adds raincloud plots using Ben Marwick's flat violin plot.}

\item{...}{Ignored}
}
\value{
A ggplot object
}
\description{
The plot functions are designed to run with just data and pupil selections,
with some additional options for fun with plotting. To see these plots,
you must first use create_window_data.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
regressed_data <- regress_data(data = Sdata, pupil1 = RPupil, pupil2 = LPupil)
mean_data <- calculate_mean_pupil_size(data = regressed_data,
pupil1 = RPupil, pupil2 = LPupil)
base_data <- baseline_data(data = mean_data, pupil = mean_pupil, start = 0, stop = 100)
window <- create_window_data(data = base_data,pupil = mean_pupil)
p <-plot(window, pupil = mean_pupil, windows = FALSE, geom = 'boxplot')
p

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geom_flat_violin.R
\docType{data}
\name{GeomFlatViolin}
\alias{GeomFlatViolin}
\title{geom_flat_violin_HELPER2}
\description{
Borrowed from
\href{https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R}{Ben Marwick}.
Original author David Robinson.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pupil_data.R
\name{plot.Pupil_difference_data}
\alias{plot.Pupil_difference_data}
\title{Pre-prepared plots of PupillometryR data}
\usage{
\method{plot}{Pupil_difference_data}(x, pupil, geom = c("point", "line"), colour = "black", ...)
}
\arguments{
\item{x}{A Pupil_difference_data dataframe}

\item{pupil}{Column name of pupil data to be plotted}

\item{geom}{string indicating whether made of connected points or a line}

\item{colour}{string indicating colour of geom, passed to ggplot2}

\item{...}{Ignored}
}
\value{
A ggplot object
}
\description{
The plot functions are designed to run with just data and pupil selections,
with some additional options for fun with plotting. To see these plots,
you must first use create_difference_data.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
regressed_data <- regress_data(data = Sdata, pupil1 = RPupil, pupil2 = LPupil)
mean_data <- calculate_mean_pupil_size(data = regressed_data,
pupil1 = RPupil, pupil2 = LPupil)
base_data <- baseline_data(data = mean_data, pupil = mean_pupil, start = 0, stop = 100)
differences <- create_difference_data(data = base_data,
pupil = mean_pupil)
p <- plot(differences, pupil = mean_pupil, geom = 'line')
p
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geom_flat_violin.R
\name{geom_flat_violin}
\alias{geom_flat_violin}
\title{ggplot Flat Violin}
\usage{
geom_flat_violin(
  mapping = NULL,
  data = NULL,
  stat = "ydensity",
  position = "dodge",
  trim = TRUE,
  scale = "area",
  show.legend = NA,
  inherit.aes = TRUE,
  ...
)
}
\arguments{
\item{mapping}{A value}

\item{data}{A value}

\item{stat}{A value}

\item{position}{A value}

\item{trim}{A value}

\item{scale}{A value}

\item{show.legend}{A value}

\item{inherit.aes}{A value}

\item{...}{A value}
}
\description{
ggplot Flat Violin
}
\details{
Copy-pasted from https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R
somewhat hackish solution to:
https://twitter.com/EamonCaddigan/status/646759751242620928
based mostly on copy/pasting from ggplot2 geom_violin source:
https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r
The original seems to be: sourced from: https://gist.github.com/dgrtwo/eb7750e74997891d7c20,
Author is David Robinson.
A key internal function for the raincloud plots used as a plotting option in this package.
For information on raincloud plots see: Allen, M., Poggiali, D., Whitaker, K.,
Marshall, T. R., & Kievit, R. A. (2019). Raincloud plots: a multi-platform
tool for robust data visualization. Wellcome open research,
4, 63. doi:10.12688/wellcomeopenres.15191.1
}
\examples{
ggplot(diamonds, aes(cut, carat)) +
  geom_flat_violin() +
  coord_flip()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_and_downsample.R
\name{calculate_missing_data}
\alias{calculate_missing_data}
\title{Calculate the missing data amount}
\usage{
calculate_missing_data(data, pupil)
}
\arguments{
\item{data}{your data of class PupillometryR}

\item{pupil}{a column name denoting pupil size}
}
\value{
A summary table with number of missing samples in each trial
}
\description{
This function can be used to assess the amount of samples that have problematic
data from each trial, which helps assess cleaning parameters
}
\examples{
data(pupil_data)
Sdata <- make_pupillometryr_data(data = pupil_data,
subject = ID,
trial = Trial,
time = Time,
condition = Type)
new_data <- downsample_time_data(data = Sdata,
pupil = LPupil,
timebin_size = 50,
option = 'mean')
calculate_missing_data(data = new_data, pupil = LPupil)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functional_data.R
\name{create_functional_data}
\alias{create_functional_data}
\title{Makes a functional data with splines from a Pupil_difference_data dataframe.}
\usage{
create_functional_data(data, pupil, basis, order)
}
\arguments{
\item{data}{a Pupil_difference_data dataframe}

\item{pupil}{Column name indicating pupil data to fit}

\item{basis}{Integer specifying number of basis functions to create a b-spline basis}

\item{order}{Integer specifying order of b-splines (one higher than the degree)}
}
\value{
A Pupil_difference_data dataframe fitted with b-splines.
}
\description{
This function turns difference data into fitted splines in order to carry out functional data analysis.
Under the hood this passes basis and order to fda::Data2fd, and fda::create.bspline.basis, and is
mandatory before running run_functional_t_test. It is recommended to read the documentation for
package fda for further information.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
regressed_data <- regress_data(data = Sdata, pupil1 = RPupil, pupil2 = LPupil)
mean_data <- calculate_mean_pupil_size(data = regressed_data, pupil1 = RPupil, pupil2 = LPupil)
base_data <- baseline_data(data = mean_data, pupil = mean_pupil, start = 0, stop = 100)
differences <- create_difference_data(data = base_data, pupil = mean_pupil)
spline_data <- create_functional_data(data = differences, pupil = mean_pupil, basis = 10, order = 4)

}
\seealso{
fda package
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pupil_data.R
\docType{data}
\name{pupil_data}
\alias{pupil_data}
\title{Data collected in a pupillometry study by Sylvain Sirois}
\format{
A data frame with 28800 rows and 7 variables: \describe{
\item{ID}{Uniaue participant ID}
  \item{Trial}{Unique trial code (also unique for each participant)}
  \item{RPupil}{Right pupil size}
  \item{LPupil}{Left Pupil Size}
  \item{Timebin}{Ordered timebin within each trial}
  \item{Time}{Elapsed time within trial}
  \item{Type}{Hard or easy trial?} ... }
}
\source{
(https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=)
}
\usage{
pupil_data
}
\description{
Data from a simple study measuring pupil dilation as participants answer hard or easy maths problems.
Original data sourced and reformatted from Sylvain Sirois' Pupillometry tutorial available at https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=314&owa_no_fiche=3&owa_bottin=)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geom_flat_violin.R
\name{gfv_helper1}
\alias{gfv_helper1}
\alias{\%||\%}
\title{geom_flat_violin_HELPER1}
\description{
Borrowed from
https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R.
Original author David Robinson, from https://gist.github.com/dgrtwo/eb7750e74997891d7c20
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baseline_data.R
\name{baseline_data}
\alias{baseline_data}
\title{Baseline pupil data to the average pupil size within a window}
\usage{
baseline_data(data, pupil, start, stop)
}
\arguments{
\item{data}{a PupillometryR dataframe}

\item{pupil}{a column name denoting pupil data}

\item{start}{start time of baseline window}

\item{stop}{stop time of baseline window}
}
\value{
A PupillometryR dataframe, with baselined pupil
}
\description{
This function is for use with the PupillometryR package to baseline each participant's pupil size to the
mean pupil size within a window.
This may not be necessary if you are doing purely within-subject analyses, but it is convenient for
comparison across subjects, and makes results more uniform.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
regressed_data <- regress_data(data = Sdata, pupil1 = RPupil, pupil2 = LPupil)
mean_data <- calculate_mean_pupil_size(data = regressed_data,
pupil1 = RPupil, pupil2 = LPupil)
base_data <- baseline_data(data = mean_data, pupil = mean_pupil, start = 0, stop = 100)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pupil_data.R
\name{plot.Pupil_test_data}
\alias{plot.Pupil_test_data}
\title{Pre-prepared plots of PupillometryR data}
\usage{
\method{plot}{Pupil_test_data}(x, show_divergence = TRUE, colour = "black", fill = "grey", ...)
}
\arguments{
\item{x}{A Pupil_test_data dataframe}

\item{show_divergence}{logical indicating whether divergences are to be highlighted}

\item{colour}{string indicating colour of geom_line, passed to ggplot2}

\item{fill}{string indicating fill hue of divergence highlights, passed to ggplot2}

\item{...}{Ignored}
}
\value{
A ggplot object
}
\description{
The plot functions are designed to run with just data and pupil selections,
with some additional options for fun with plotting. To see these plots,
you must first use one of the run_functional tests.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
regressed_data <- regress_data(data = Sdata, pupil1 = RPupil, pupil2 = LPupil)
mean_data <- calculate_mean_pupil_size(data = regressed_data,
pupil1 = RPupil, pupil2 = LPupil)
base_data <- baseline_data(data = mean_data, pupil = mean_pupil, start = 0, stop = 100)
differences <- create_difference_data(data = base_data,
pupil = mean_pupil)
spline_data <- create_functional_data(data = differences, pupil = mean_pupil, basis = 10, order = 4)
ft_data <- run_functional_t_test(data = spline_data,
pupil = mean_pupil)
p <- plot(ft_data, show_divergence = TRUE, colour = 'red', fill = 'orange')
p
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_smoothed_data.R
\name{interpolate_data}
\alias{interpolate_data}
\title{Interpolate across the gaps in data}
\usage{
interpolate_data(data, pupil, type = c("linear", "cubic"))
}
\arguments{
\item{data}{a PupillometryR dataframe}

\item{pupil}{Column name for pupil data to be interpolated}

\item{type}{string indicating linear or cubic interpolation to be performed.}
}
\value{
interpolated pupillometry data
}
\description{
Once data is smoothed, it is important to deal with missing observations, such as blinks.
This allows simple interpolation over missing values, either linear, or cubic.
Depending on the analysis planed, this may not be a necessary option, but it is
strongly recommended for the functional analyses planned in this package.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
subject = ID,
trial = Trial,
time = Time,
condition = Type)
mean_data <- calculate_mean_pupil_size(data = Sdata,
pupil1 = RPupil, pupil2 = LPupil)
filtered_data <- filter_data(data = mean_data,
pupil = mean_pupil,
filter = 'hanning',
degree = 11)
int_data <- interpolate_data(data = filtered_data,
pupil = mean_pupil,
type = 'linear')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_mean_pupil_size.R
\name{calculate_mean_pupil_size}
\alias{calculate_mean_pupil_size}
\title{Calculate a mean size across two pupils over time}
\usage{
calculate_mean_pupil_size(data, pupil1, pupil2)
}
\arguments{
\item{data}{a PupillometryR dataframe}

\item{pupil1}{column name indicating pupil size}

\item{pupil2}{column name indicating pupil size}
}
\value{
A PupillometryR dataframe with a mean pupil column
}
\description{
This function is useful when you have left and right eye eyetracking data, and a mean of the two would be useful.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
regressed_data <- regress_data(data = Sdata, pupil1 = RPupil, pupil2 = LPupil)
mean_data <- calculate_mean_pupil_size(data = regressed_data, pupil1 = RPupil, pupil2 = LPupil)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functional_data.R
\name{run_functional_t_test}
\alias{run_functional_t_test}
\title{Run a functional t-test on a dataframe previously fitted with b-splines.}
\usage{
run_functional_t_test(data, pupil, alpha = 0.05)
}
\arguments{
\item{data}{a Pupil_difference_data fitted with b-splines}

\item{pupil}{column name indicating pupil data to test}

\item{alpha}{an alpha level to be set for the t-test}
}
\value{
A Pupil_test_data dataframe
}
\description{
This allows running of a functional t-test for a given alpha on pupil data that has been fitted with b-splines.
This is only appropriate for functional difference data, as it assumes we are dealing with condition A - condition B.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
regressed_data <- regress_data(data = Sdata, pupil1 = RPupil, pupil2 = LPupil)
mean_data <- calculate_mean_pupil_size(data = regressed_data, pupil1 = RPupil, pupil2 = LPupil)
base_data <- baseline_data(data = mean_data, pupil = mean_pupil, start = 0, stop = 100)
differences <- create_difference_data(data = base_data, pupil = mean_pupil)
spline_data <- create_functional_data(data = differences, pupil = mean_pupil, basis = 10, order = 4)
ft_data <- run_functional_t_test(data = spline_data, pupil = mean_pupil, alpha = 0.05)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_pupillometryR_data.R
\name{make_pupillometryr_data}
\alias{make_pupillometryr_data}
\title{Prepare data for pre-processing in PupillometryR}
\usage{
make_pupillometryr_data(data, subject, trial, time, condition, other)
}
\arguments{
\item{data}{a raw, long form dataframe organised by subject, trial, and time.
if your data is not long form, look at tidyr for examples of conversion.}

\item{subject}{column name indicating subject ID}

\item{trial}{column name indicating trial ID. This should be unique for participants}

\item{time}{column name indicating time column (should be numeric)}

\item{condition}{column name indicating experimental condition}

\item{other}{any other column you may wish to keep in the data frame for processing}
}
\value{
A dataframe ready to use in PupillometryR
}
\description{
This should be the first function you run as part of using PupillometryR.
This will make sure your data is in the right format for processing.
This package is designed to deal with data at it comes out of the eyetracker
in a long-form csv style format. Thus data input here would be a long
dataframe, wherein each row is a single frame collected by the eyetracker.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
subject = ID,
trial = Trial,
time = Time,
condition = Type)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functional_data.R
\name{create_difference_data}
\alias{create_difference_data}
\title{Create a difference data frame when dealing with a condition column with 2 levels}
\usage{
create_difference_data(data, pupil)
}
\arguments{
\item{data}{a PupillometryR dataframe}

\item{pupil}{column name for pupil data}
}
\value{
A Pupil_difference_data data frame
}
\description{
The difference data frame is used when creating a dataframe to do the functional t-test analysis.
This function would be the first step in that analysis, after doing the pre-processing.
It creates a frame where it treats the condition data as level2 - level1.
It will throw an error if there are more than two conditions.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
mean_data <- calculate_mean_pupil_size(data = Sdata,
pupil1 = RPupil, pupil2 = LPupil)
base_data <- baseline_data(data = mean_data, pupil = mean_pupil, start = 0, stop = 100)
differences <- create_difference_data(data = base_data, pupil = mean_pupil)
plot(differences, pupil = mean_pupil, geom = 'line')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset_data.R
\name{create_window_data}
\alias{create_window_data}
\title{Make PupillometryR dataframe into a single collapsed window for easy analysis}
\usage{
create_window_data(data, pupil)
}
\arguments{
\item{data}{a PupillometryR dataframe}

\item{pupil}{column name denoting pupil data to be used}
}
\value{
a Pupil_window_data dataframe
}
\description{
This function creates a single collapsed data frame for easy analysis with a t-test or anova,
per condition.
By comparison create_time_windows allows dividing it into multiple windows per time.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
regressed_data <- regress_data(data = Sdata, pupil1 = RPupil, pupil2 = LPupil)
mean_data <- calculate_mean_pupil_size(data = regressed_data,
pupil1 = RPupil, pupil2 = LPupil)
base_data <- baseline_data(data = mean_data, pupil = mean_pupil, start = 0, stop = 100)
window <- create_window_data(data = base_data, pupil = mean_pupil)
p <- plot(window, pupil = mean_pupil, windows = FALSE, geom = 'boxplot')
p
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset_data.R
\name{subset_data}
\alias{subset_data}
\title{Subset data to provide start and finish time windows}
\usage{
subset_data(data, start = NULL, stop = NULL, rezero = T, remove = T)
}
\arguments{
\item{data}{a PupillometryR dataframe}

\item{start}{a single number indicating start time of new dataframe}

\item{stop}{a single number indicating end time of new dataframe}

\item{rezero}{logical, whether time should start from zero}

\item{remove}{logical, remove observations outside of start and stop}
}
\value{
a subsetted PupillometryR dataframe
}
\description{
subset_data can be used on a PupillometryR dataframe to subset the time into relevant chunks.
This, ideally should be one of the first runctions run, before anything analytical.
Use this to indicate a start and stop time to create a new resized dataframe.
}
\examples{
Sdata <- make_pupillometryr_data(data = pupil_data,
                               subject = ID,
                               trial = Trial,
                               time = Time,
                               condition = Type)
subset_data(Sdata, start = 100, stop = 10000, rezero = TRUE, remove = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_smoothed_data.R
\name{replace_missing_data}
\alias{replace_missing_data}
\title{replaces missing observations if you have some degree of incomplete observations}
\usage{
replace_missing_data(data)
}
\arguments{
\item{data}{your data of class pupillometryR}
}
\value{
A time-stepped data frame
}
\description{
This is a useful function if you have a dataset where certain timepoints have been removed for whatever reason,
but you want continuous time data. This will make assumptions about trials being the same length though,
so may not be appropriate for all data types.
This should only be run after running make_pupillometry_data.
}
\examples{
data(pupil_data)
Sdata <- make_pupillometryr_data(data = pupil_data,
subject = ID,
trial = Trial,
time = Time,
condition = Type)
new_data <- replace_missing_data(data = Sdata)
}
