
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
