![](vignettes/GGIR-MASTERLOGO-RGB.png)

![GitHub Actions R-CMD-check](https://github.com/wadpac/GGIR/workflows/R-CMD-check-full/badge.svg)
[![codecov](https://codecov.io/gh/wadpac/GGIR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/wadpac/GGIR)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1051064.svg)](https://doi.org/10.5281/zenodo.1051064)
[![](https://cranlogs.r-pkg.org/badges/last-month/GGIR)](https://cran.r-project.org/package=GGIR)

## Getting started:
The package [vignette](https://CRAN.R-project.org/package=GGIR/vignettes/GGIR.html) and [this](https://youtu.be/S8YPTrYNWdU) short tutorial video provide an introduction to GGIR, including: How it can be installed, Key software features, and where to get help.

## Contribution guidelines:
We always welcome contributions to the package.
If you want to contribute to the development of GGIR, have a look at the [contribution guidelines](https://github.com/wadpac/GGIR/blob/master/CONTRIBUTING.md).

### Images usaged
The copyright of the GGIR logo as contained in the file vignettes/GGIR-MASTERLOGO-RGB.png lies with Accelting (Almere, The Netherlands), please contact v.vanhees@accelting.com to ask for permission to use this logo.

All other images in this repository are released under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

### Research notice

Please note that this repository is participating in a study into sustainability of open source projects. Data will be gathered about this repository for approximately the next 12 months, starting from June 2021.

Data collected will include number of contributors, number of PRs, time taken to close/merge these PRs, and issues closed.

For more information, please visit [the informational page](https://sustainable-open-science-and-software.github.io/) or download the [participant information sheet](https://sustainable-open-science-and-software.github.io/assets/PIS_sustainable_software.pdf).

# Version numbering

We use version encoding **A.B-C**:

- A increases with major changes that affect backward compatibility with previous releases like changes in function names, function arguments or file format.
- B increases with every CRAN release. We aim to avoid more than four CRAN releases per year.
- C increases with every GitHub release. We aim to avoid more than one GitHub release per month.

# GitHub releases

Before releasing, please make sure to check the following:

1. Create GitHub issue at least 1 weeks before the intended release to announce the release and indicate what will be in the release.
2. Make sure the change log `inst/NEWS.Rd` is up to date and that it says "GitHub-only-release date" rather than "release date"
3. Make sure the third (last) digit in the version number is incremented by one relative to the master branch and the date is the present date. This applies to the files `DESCRIPTION`, `CITATION.cff` (not the cff-version, but the version on line 56 of the .cff-file), `GGIR-package.Rd` and `NEWS.Rd` file. Use function `prepareNewRelease.R` in the root of GGIR to double check that version number and date are consistent between these files.
4. Update package contributor list if new people have contributed.
5. Run `R CMD check --as-cran` to make sure all tests and checks pass.

Note that GitHub releases require a release name. We typically choose a random name of a city or town in South America. Whatever you choose this should be an easy to read and remember word.

# CRAN releases

To do a CRAN release, follow the following steps:

1. Create GitHub issue at least 4 weeks before the intended CRAN release announcing the release and indicating what will be in the release and a to do list.
2. A CRAN release should not come with major changes that have not been covered by any of the GitHub-only releases.
3. When everything looks ready for the release, repeat the same process as for the GitHub release with a few differences:
    - In the change log it should now say "release date" rather than "GitHub-only-release date".
    - Second digit in the version number is incremented by 1 relative to the current CRAN version.
    - Check whether a new R version has been released or is coming up and make sure GGIR is also tested with that version.
    - Run in RStudio `devtools::check( manual = TRUE, remote = TRUE, incoming = TRUE)` which will help to check urls
4. Ask Vincent (GitHub tag: vincentvanhees) to submit the release to CRAN as it needs to come from my e-mail address.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the main project maintainer v.vanhees@accelting.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. [you have a question](#questions);
2. [you think you may have found a bug](#bugs) (including unexpected behavior);
3. [you want to make some kind of change to the code base](#changes-or-additions) (e.g. to fix a bug, to add a new feature, to update documentation);
4. [you want to make a new release of the code base](#new-release).

The sections below outline the steps in each case.

## Questions

1. use the search functionality [here](https://groups.google.com/g/RpackageGGIR) to see if someone already experienced the same issue;
2. if your search did not yield any relevant results, start a new conversation.

## Bugs

1. use the search functionality [here](https://github.com/wadpac/GGIR/issues) to see if someone already filed the same issue;
2. if your issue search did not yield any relevant results, make a new issue, and choose the Bug report type. This includes a checklist to make sure you provide enough information to the rest of the community to understand the cause and context of the problem.

## Changes or additions

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue. Choose the Feature request type, which includes a checklist of things to consider to get the discussion going;
2. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
3. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
4. make sure the existing tests still work by running the test suite from RStudio;
5. add your own tests (if necessary);
6. update or expand the documentation;
7. make sure the release notes in `inst/NEWS.Rd` are updated;
8. add your name to the contributors lists in the `DESCRIPTION` and `CITATION.cff` files;
9. push your feature branch to (your fork of) the GGIR repository on GitHub;
10. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/). The pull request template includes a checklist with the above items.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.

### Coding style

We loosely follow the [tidyverse style guide](https://style.tidyverse.org/), but do not enforce every rule strictly.
For instance, we prefer `=` instead of `<-` as the default assignment operator.
When in doubt about what style to use, don't hesitate to get in touch.

Some general guidelines that we try to adhere to:

- Use standard R as much as possible, keep dependencies to a minimum.
- Keep loops to a minimum.
- Don't make lines too long.

If you are a first time contributor, don't worry about coding style too much.
We will help you get things in shape.

## New release

GGIR follows the [release cycle process described in this document](RELEASE_CYCLE.md).<!-- Describe your PR here -->

<!-- Please, make sure the following items are checked -->
Checklist before merging:

- [ ] Existing tests still work (check by running the test suite, e.g. from RStudio).
- [ ] Added tests (if you added functionality) or fixed existing test (if you fixed a bug).
- [ ] Updated or expanded the documentation.
- [ ] Updated release notes in `inst/NEWS.Rd` with a user-readable summary. Please, include references to relevant issues or PR discussions.
- [ ] Added your name to the contributors lists in the `DESCRIPTION` and `CITATION.cff` files.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A short, clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior.

1. Sensor brand: '...'
2. Data format: '...'
3. Approximate recording duration '...' days
4. Are you using a sleep diary to guide the sleep detection: YES / NO
5. Copy of R command used: '....'
6. Have you tried processing your data based on GGIR's default argument values? Does the issue you report still appear? YES / NO
7.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem. Note that usually we are not only interested in see the error message in red, but all GGIR output to the console.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - GGIR Version [e.g. 2.2-0]

**Additional context**
Add any other context about the problem here.

**Before submitting**
- [ ] Have you tried the steps to reproduce? Do they include all relevant data and configuration? Does the issue you report still appear there?
- [ ] Have you tried this on the latest `master` branch from GitHub?
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
title: "Embedding external functions in GGIR"
author: "Vincent van Hees"
date: "January 27 2021"
output:
  pdf_document:
    toc: true
    toc_depth: 3
    number_sections: true
urlcolor: blue
vignette: >
  %\VignetteIndexEntry{Embedding external functions in GGIR}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
header-includes:
  - \usepackage{titling}
  - \pretitle{\begin{center}
    \includegraphics[]{GGIR-MASTERLOGO-RGB.png}\LARGE\\}
  - \posttitle{\end{center}}
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  python.reticulate = FALSE
)
```


# Introduction

If you like GGIR but want to use algorithms for raw data not included in GGIR then the external function embedding feature can be the solution. For example, you may want to pilot a new machine learned classifiction algorithm but you do not want to write all the data cleaning and aggregation steps needed for analysis of real life 'out of the lab' acceleormeter data.

**How it works:**

Internally GGIR loads the raw accelerometer data in memory blocks of about 24 hours. When the data is in memory, corrected for calibration error, and resampled to the sample rate required by your function, GGIR applies its own default algorithms as well as the external function provided by you (Python or R). The external function is expected to take as input: A three-column matrix with the acceleration data corresponding to the three acceleration axes, and an optional parameters argument which can be of any R format (character, list, vector, data.frame, etc). As output your external function is expected to produce a matrix or data.frame with one or multiple columns corresponding to the output of your external function.

# Example with external R function

In this example we will apply the function counts() from R package [activityCounts](https://CRAN.R-project.org/package=activityCounts) to the raw data, which produces an estimate of Actigraph counts per second.

## Write external function

Create file **calculateCounts.R** and insert the following code:

```{R,eval=FALSE}
calculateCounts = function(data=c(), parameters=c()) {
  # data: 3 column matrix with acc data
  # parameters: the sample rate of data
  library("activityCounts")
  if (ncol(data) == 4) data= data[,2:4]
  mycounts = counts(data=data, hertz=parameters, 
                    x_axis=1, y_axis=2, z_axis=3,
                    start_time = Sys.time())
  mycounts = mycounts[,2:4] #Note: do not provide timestamps to GGIR
  return(mycounts)
}
```

## Provide external function to GGIR

Create a new .R file for running the GGIR analysis, e.g. named `myscript.R`, and insert the following code.
Do not forget to update the filepath on the first line to point to your `calculateCounts.R` file.

```{R,eval=FALSE}
source("~/calculateCounts.R")
myfun =  list(FUN=calculateCounts,
              parameters= 30,
              expected_sample_rate= 30,
              expected_unit="g",
              colnames = c("countsX","countsY","countsZ"),
              outputres = 1,
              minlength = 1,
              outputtype="numeric",
              aggfunction = sum,
              timestamp=F,
              reporttype="scalar")
```

The above code creates object `myfun` of type `list` which is expected to come with the following elements:

- `FUN` A character string specifying the location of the external function you want to apply.
- `parameters` The parameters used by the function, which can be stored in any format (vector, matrix, list, data.frame). The user should make sure that the external function can handle this object.
- `expected_sample_rate` Expected sample rate, if the inputdata has a difference sample rate, then the data will be resampled.
- `expected_unit` Expected unit of the acceleration by external function: "mg", "g" or "ms2". If input data is different it will be converted.
- `colnames` Character vector with the names of the columns produced by the external function.
- `outputres` The resolution (seconds) of the output produced by the external function. Note, that this needs to be equal to or a multitude of the short epoch size of the g.part1 output (5 seconds) or the short epoch size should be a multitude of this resolution. In this way GGIR can aggregate or repeat the external function output to be used inside GGIR.
- `minlength` The minimum length (seconds) of input data needed, typically the window per which output is provided.
- `outputtype` Character to indicate the type of external function output. Set to "numeric" if data is stored in numbers (any numeric format), or "character" if it is a character string.
- `aggfunction` If the data needs to be aggregated to match the short epoch size of the g.part1 output (5 seconds) then this element specifies what
function should be used for the aggregation, e.g. mean, sum, median.
- `timestamp` Boolean to indicated whether timestamps  (seconds since 1-1-1970) should be passed on to the external 
    function as first columm of the data matrix..
- `reporttype` Character to indicate the type of reporting by GGIR: "scalar" if it should be averaged per day, "event" if it should be summed per day, or "type" if it is categorical variable that can only be aggregated per day by tabulating it.

Next, add a call to GGIR function g.shell.GGIR with `myfun` provided as one of its arguments:

```{R,eval=FALSE}
library(GGIR)
g.shell.GGIR(datadir="~/myaccelerometerdata",
             outputdir="~/myresults",
             mode=1:2,
             epochvalues2csv = TRUE,
             do.report=2,
             myfun=myfun) #<= this is where object myfun is provided to g.shell.GGIR
```

Please see the [\underline{general GGIR vignette}](https://cran.r-project.org/package=GGIR/vignettes/GGIR.html) for more information about function g.shell.GGIR.

# Example with external Python function

In this example we will use an external Python function to estimate the dominant signal frequency per acceleration axis. Note this can also be done in R, but it shows that even Python functions can be provided. 

## Write external function

Create **dominant_frequency.py** and insert the code shown below:

```{python, eval=FALSE}
import numpy

def dominant_frequency(x, sf):
  # x: vector with data values
  # sf: sample frequency
  fourier = numpy.fft.fft(x)
  frequencies = numpy.fft.fftfreq(len(x), 1/sf)
  magnitudes = abs(fourier[numpy.where(frequencies > 0)])
  peak_frequency = frequencies[numpy.argmax(magnitudes)]
  return peak_frequency
```

Create **dominant_frequency.R** that calls the python function and insert the following code:

```{R,eval=FALSE}
dominant_frequency = function(data=c(), parameters=c()) {
  # data: 3 column matrix with acc data
  # parameters: the sample rate of data
  source_python("dominant_frequency.py")
  sf=parameters
  N = nrow(data)
  ws = 5 # windowsize
  if (ncol(data) == 4) data= data[,2:4]
  data = data.frame(t= floor(seq(0,(N-1)/sf,by=1/sf)/ws),
                    x=data[,1], y=data[,2], z=data[,3])
  df = aggregate(data, by = list(data$t), 
                 FUN=function(x) {return(dominant_frequency(x,sf))})
  df = df[,-c(1:2)]
  return(df)
}
}
```

## Provide external function to GGIR

Create a new .R file for running the GGIR analysis, e.g. named myscript.R, and insert the following blocks of code.

Specification of Python environment to use, this can also be a conda environment or docker container (see documentation R package [\underline{reticulate}](https://rstudio.github.io/reticulate/) for further details). Make sure that that this Python environment has all the required dependencies for the external function, here we will only need `numpy`.

```{R,eval=FALSE}
  library("reticulate")
  use_virtualenv("~/myvenv", required = TRUE) # Local Python environment
  py_install("numpy", pip = TRUE)

```

Specify a `myfun` object as explained in the R example. Do not forget to update the filepath to the `"~/dominant_frequency.R"` file.

```{R,eval=FALSE}
source("~/dominant_frequency.R")
myfun =  list(FUN=dominant_frequency,
              parameters= 30,
              expected_sample_rate= 30,
              expected_unit="g",
              colnames = c("domfreqX", "domfreqY", "domfreqZ"),
              minlength = 5,
              outputres = 5,
              outputtype="numeric",
              aggfunction = median
              timestamp=F,
              reporttype="scalar")
```

Add a call to GGIR function g.shell.GGIR where `myfun` is provided as argument.
Note that, do.parallel is set to FALSE. Unfortunately Python embedding with R package reticulate and multi-threading with R package foreach as used in GGIR do not combine well.

```{R,eval=FALSE}
library(GGIR)
g.shell.GGIR(datadir="~/myaccelerometerdata",
             outputdir="~/myresults",
             mode=1:2,
             epochvalues2csv = TRUE,
             do.report=2,
             myfun=myfun,
             do.parallel = FALSE)
```


# Integration in GGIR output

## Part 1
The external function output is included in the time series produced by function GGIR function g.part1 and stored in an RData-file in `/output_nameofstudy/meta/basic`.
The resolution of these output in GGIR is set by g.shell.GGIR argument `windowsizes`, which is `c(5,900,3600)` by default. Here, the first element `5` specifies the short epoch size in seconds.
If the output of the external function is less then this resolution it will be aggregated with the function as specificied by aggfunction in the `myfun` object. In the count example we used the sum for this and for the dominant frequency example we used the median.

## Part 2
Next, in part2 GGIR aims to detect non-wear periods and imputes those. The impute time series can be found in the part 2 milestone data in folder: `/output_nameofstudy/meta/ms2.out`. If you want these to be directly stored in a csv file then set argument `epochvalues2csv = TRUE`.



# External functions released by GGIR collaborators:

- Wrist-based step detection algorithm: https://github.com/ShimmerEngineering/Verisense-Toolbox/tree/master/Verisense_step_algorithm
- Wrist-based sleep classification as described by Sundararajan et al. 2021 [link to paper](https://www.nature.com/articles/s41598-020-79217-x), and corresponding code is here: https://github.com/wadpac/SleepStageClassification/tree/master/ggir_ext---
title: "Day segment analyses with GGIR"
output:
   html_document:
    toc : true
    number_sections: true
    toc_depth: 3
    toc_float: true #by turning this on, the table of contents moves to the left of the page.
urlcolor: blue
vignette: >
  %\VignetteIndexEntry{Day segment analyses with GGIR}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


```{r, echo=FALSE, out.width = "60%", out.extra='style="border: 0; padding:20px"'}
knitr::include_graphics("GGIR-MASTERLOGO-RGB.png")
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Introduction

Is this specific person more active in the morning or in the afternoon? Are children more active during their work hours or their leisure time? Questions like these can be answered with GGIR but you first have to specify a few parameters.

The main input argument to be specified is `qwindow`, which can be used the following ways:

- To specify the clock hours in the day based on which the segmented day analyses should take place.
- To specify an activity log (diary) that should be used to guide the segmentation per individual and per day of the recording.

In the following sections I will discuss both scenarios.

# Clock hour-based segmentation

To perform clock hour segmentation, you will need to provide function g.shell.GGIR with
argument `qwindow` and assign it a numeric vector with the hours for the segmentation.
If the start and end of the day, are not explicitly provided in the 
vector GGIR will add them.  Please find below some example values for `qwindow`.
The number of values used by `qwindow` is unlimited, but be aware that some of the analyses are such as MX-metrics are
impossible for very small windows and will produce empty results.

qwindow value | Resulting segment(s) to be analysed
-------------------| -------------------
c(0,24) | midnight to following midnight (24 hours), the full day is the only segment.
c(8,24) | midnight-8:00 (8 hour segment), 8:00-midnight (16 hour segment), and midnight-midnight (24 hour segment).
c(6,11, 13, 17) | midnight-6:00 (6 hour segment), 6:00-11:00 (5 hour segment), 11:00-13:00 (2 hour segment), 13:00-17:00 (4 hour segment), 17:00-midnight (7 hour segment), and midnight-midnight (24 hour segment).
c(0:24) | 25 segments: 24 segments of 1 hour corresponding to each hour of the day, and midnight-midnight (24 hour segment).    

Day Saving Time (DST) is taken into account when identifying the start of the day, but
not when identifying the day segments. In other words, a 23 hour days is processed as the 24 hours after the first midnight.
This to ensure that segment length is identical across days of the week, which is needed to ease comparison of outcome variables across days.


# Segmentation guided by activity log

To perform activity-log based segmentation, you will need to provide function `g.shell.GGIR` with
argument `qwindow` and assign it the full path to your activity log in `.csv format`, e.g. `qwindow="C:/myactivitylog.csv"`.

The activity log is expected to be a .csv-file with the following structure:

ID | date | work | travelhome | home | date | work | travelhome | home 
-----| -------- | ------- | --------- | ----- | -------- | ------- | -------- | --------
1234 | 04-11-20 | 7:45:00 |  17:00:00 | 17:30 | 05-11-20 |  |  | 17:30
4567 | 24-11-20 | 7:45:00 |  17:00:00 | 17:30 | 25-11-20 | 7:45:00 |  17:00:00 | 17:30

**Rows:**
First row represents the column headers after which each row represents one accelerometer recording.

**ID-column:**
The first column is expected to hold the recording ID, which needs to match with the ID GGIR extracts from the accelerometer file. If unsure how to format the ID values, apply g.shell.GGIR to a sample of your accelerometer files using the default argument settings. The ID column in the generated part 2 .csv reports will show how the participant ID is extracted by GGIR. If no ID is extracted, see documentation for argument `idloc`, which helps you to specify the location of the participant in the file name or file header. If ID extraction fails the accelerometer files cannot be matched with the corresponding activity log entries.

**Date-column:**
The ID column is followed by a date column for the first log day. To ensure GGIR recognises this date correctly, specify argument `qwindow_dateformat`. The default format is `"\%d-\%m-\%Y"` as in 23-2-2021 to indicate the 23rd of February 2021. If your date is formatted as 2-23-21 then specify`"\%md-\%d-\%y"`. The column name of the date column needs to include the character combination "date" or "Date" or "DATE".

**Start-times:**
The date column is followed by one or multiple columns with start times for the activity types in that day format in hours:minutes:seconds. Do not provide dates in these cells. The header of the column will be used as label for each activity type. Insert a new date column before continuing with activity types for next day. Leave missing values empty. 

**Notes:**
- If an activity log was collected for some individuals then those will be processed with qwindow value c(0,24).
- Dates with no activity log data can be skipped, no need to have a column with the date followed by a column with the next date.
- The end time of one activity is assumed to be the start time of the next activity. We currently do not facilitate overlapping time segments.


# Examples

For more information about how to use the g.shell.GGIR function call
see explanation in the main [GGIR vignette](https://CRAN.R-project.org/package=GGIR/vignettes/GGIR.html).

## Clock-hour based segmentation:

```{R,eval=FALSE}
library("GGIR")
g.shell.GGIR(datadir = "/your/data/directory",
             outputdir = "/your/output/directory",
             mode = 1:2, # <= run GGIR parts 1 and 2
             do.report = 2, # <= generate csv-report for GGIR part 2
             qwindow = c(0, 6, 12, 18, 24))
```


## Activity log based segmentation:

```{R,eval=FALSE}
library("GGIR")
g.shell.GGIR(datadir = "/your/data/directory",
             outputdir = "/your/output/directory",
             mode = 1:2, # <= run GGIR parts 1 and 2
             do.report = 2, # <= generate csv-report for GGIR part 2
             qwindow = "/path/to/your/activity/log.csv")
```


After running this code GGIR creates an output folder in the output directory as specified with argument outputdir.
In the subfolder `results` you will then find three files:

- `part2_summary.csv` the recording level summary, with 1 row per recording and in recording level aggregates of day segments in columns. 
- `part2_daysummary.csv` the day level summary, with 1 row per day and day segment specific outcomes in columns.
- `part2_daysummary_longformat.csv` the day level summary in long format, such that each row represents one segment from one day in one recording.

In both `part2_summary.csv` and `part2_daysummary.csv` the column names tell you the day segment they correspond to. For example, column names ending with `_18-24hr` refer to the time segment 18:00-24:00. In `part2_daysummary_longformat.csv` the time segment is clarified via columns qwindow_timestamps and qwindow_name.


# Analyses performed per day segment

The analyses that GGIR per segment of the day, include:

**Acceleration distribution:**
Derived if argument `ilevels` is specified. You will find these under the variable names such as `[0,36)_ENMO_mg` which means time spent between 0 and 36 mg defined by acceleration metric ENMO.

**Number of valid hours of data:**
You will recognise these as `N_valid_hours_in_window` which tells yoyu the number of valid hours per time window, and `N_valid_hours` which is the number of valid hours per day.

**LXMX analysis:**
LXMX analysis, which stands for least and most active X hours of the segment. You will recognise these variable names
like `L5hr_ENMO_mg` which is the start time of the least active five hours defined by metric ENMO, and `L5_ENMO_mg` which is the average acceleration for those hours.

**Intensity gradient analysis:**
You will find these as variables that start with `ig_gradient_` See [description](https://cran.r-project.org/package=GGIR/vignettes/GGIR.html#41_Output_part_2) of GGIR part 2 output in the main GGIR vignette for further details.

**Time spent in Moderate or Vigorous Physical Activity (MVPA):**
You will find these as variables such as `MVPA_E5S_T201_ENMO` or `MVPA_E5S_B1M80%_T201_ENMO`. See [description](https://cran.r-project.org/package=GGIR/vignettes/GGIR.html#41_Output_part_2) of GGIR part 2 output in the main GGIR vignette for further details.


```{r, echo=FALSE, out.width = "60%", out.extra='style="border: 0; padding:20px"'}
knitr::include_graphics("GGIR-MASTERLOGO-RGB.png")
```

---
title: "Accelerometer data processing with GGIR"
output:
   html_document:
    toc : true
    number_sections: true
    toc_depth: 3
    toc_float: true #by turning this on, the table of contents moves to the left of the page.
urlcolor: blue
vignette: >
  %\VignetteIndexEntry{Accelerometer data processing with GGIR}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r, echo=FALSE, out.width = "60%", out.extra='style="border: 0; padding:20px"'}
knitr::include_graphics("GGIR-MASTERLOGO-RGB.png")
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Introduction

## What is GGIR?

[GGIR](https://CRAN.R-project.org/package=GGIR) is an R-package to process multi-day raw accelerometer data for physical activity and sleep research. The term _**raw**_ refers to data being expressed in m/s<sup>2</sup> or gravitational acceleration as opposed to the previous generation accelerometers which stored data in accelerometer brand specific units. The signal processing includes automatic calibration, detection of sustained abnormally high values, detection of non-wear and calculation of average magnitude of dynamic acceleration based on a variety of metrics. Next, GGIR uses this information to describe the data per recording, per day of measurement, and (optionally) per segment of a day of measurement, including estimates of physical activity, inactivity and sleep.
We published an overview paper of GGIR in 2019 [link](https://doi.org/10.1123/jmpb.2018-0063).

This vignette provides a general introduction on how to use GGIR and interpret the output, additionally you can find a [introduction video](https://youtu.be/RuFBCAqFJ2M) and a [mini-tutorial](https://youtu.be/S8YPTrYNWdU) on YouTube. If you want to use your own algorithms for raw data then GGIR facilitates this with it's external function embedding feature, documented in a separate vignette: [Embedding external functions in GGIR](https://CRAN.R-project.org/package=GGIR). GGIR is increasingly being used by research groups across the world. A non-exhaustive overview of academic publications related to GGIR can be found [here](https://github.com/wadpac/GGIR/wiki/Publication-list). R package GGIR would not have been possible without the support of the contributors listed in the author list at [GGIR](https://CRAN.R-project.org/package=GGIR), with specific code contributions over time since April 2016 (when GGIR development moved to GitHub) shown [here](https://github.com/wadpac/GGIR/graphs/contributors).

**Cite GGIR:**

When you use GGIR in publications do not forget to cite it properly as that makes your research more reproducible and it gives credit to it's developers. See paragraph on [Citing GGIR](#CitingGGIR) for details.


## Contributing, Support, and Keeping up to date

**How to contribute to the code?**

The development version of GGIR can be found on [github](https://github.com/wadpac/GGIR), which is also where you will find [guidance](https://github.com/wadpac/GGIR/blob/master/CONTRIBUTING.md) on how to contribute.


**How can I get service and support?**

GGIR is open source software and does not come with service or support guarantees. However, as user-community you can help each other via the [GGIR google group](https://groups.google.com/forum/#!forum/RpackageGGIR) or the [GitHub issue tracker](https://github.com/wadpac/GGIR/issues). Please use these public platform rather than private e-mails such that other users can learn from the conversations.

If you need dedicated support with the use of GGIR or need someone to adapt GGIR to your needs then Vincent van Hees is available as [independent consultant](https://www.accelting.com/).

**Change log**

Our log of main changes to GGIR over time can be found [here](https://cran.r-project.org/package=GGIR/news.html).

# Setting up your work environment

## Install R and RStudio

[Download and install R](https://cran.r-project.org/)

[Download and install RStudio](https://www.rstudio.com/products/rstudio/) (optional, but recommended)

Install GGIR with its dependencies from CRAN. You can do this with one command from the console command line:

```{R,eval=FALSE}
install.packages("GGIR", dependencies = TRUE)
```

Alternatively, to install the latest development version with the latest bug fixes use instead:

```{R,eval=FALSE}
install.packages("remotes")
remotes::install_github("wadpac/GGIR")
```


## Prepare folder structure

1. GGIR works with the following accelerometer brands and formats:
    - [GENEActiv](https://www.activinsights.com/) .bin and .csv
    - [Axivity](https://axivity.com/) AX3 and AX6 .wav, .csv and .cwa
    - [ActiGraph](https://actigraphcorp.com/) .csv and .gt3x (.gt3x only the newer format generated with firmware versions  above 2.5.0). Note for Actigraph users: If you want to work with .csv exports via the ActiLife then note that you have the option to export data with timestamps. Please do not do this as this causes memory issues for GGIR. To cope with the absence of timestamps GGIR will re-caculate timestamps from the sample frequency and the start time and date as presented in the file header.
    - [Movisens](https://www.movisens.com/en/) with data stored in folders.
    - Genea (an accelerometer that is not commercially available anymore, but which was used for some studies between 2007 and 2012) .bin and .csv
    - Any other accelerometer brand that generates csv output, such as [Shimmer](https://shimmersensing.com/), see function documentation for function read.myacc.csv.
2. All accelerometer data that needs to be analysed should be stored in one folder, or subfolders of that folder.
3. Give the folder an appropriate name, preferable with a reference to the study or project it is related to rather than just 'data', because the name of this folder will be used later on as an identifier of the dataset.

## GGIR shell function

GGIR comes with a large number of functions and optional settings (arguments) per functions.

To ease interacting with GGIR there is one central function, named `g.shell.GGIR`, to talk to all the other functions.

In this paragraph we will guide you through the main arguments relevant for 99% of research. First of all, it is important to understand that GGIR is structured in two ways. It has a computational structure of five parts which are applied sequentially to the data, and that `g.shell.GGIR` controls each of these parts:

- Part 1: Loads the data and stores derived features (aggregations) needed for the other parts. This is the time-consuming part. Once this is done, parts 2-5 can be run (or re-run with different parameters in parts 2-5) relatively quickly.
- Part 2: Data quality analyses and low-level description of signal features per day and per file. At this point a day is defined from midnight to midnight
- Part 3: Estimation of sustained inactivity and sleep periods, needed for input to Part 4 for sleep detection
- Part 4: Labels the sustained inactive periods detected in Part 3 as sleep, or daytime sustained inactivity, per night and per file
- Part 5: Derives sleep and physical activity characteristics by re-using information derived in part 2, 3 and 4. Total time in intensity categories, the number of bouts, time spent in bouts and average acceleration (overall activity) is calculated.

The reason why it split up in parts is that it avoids having the re-do all analysis if you only want to make a small change in the more downstream parts. The specific order and content of the parts has grown for historical and computational reasons.

Secondly, the function arguments which we will refer to as input parameters are structured thematically independently of the five parts they are used in:

- params_rawdata: parameters related to handling the raw data such as resampling or calibrating
- params_metrics: parameters related to aggregating the raw data to epoch level summary metrics
- params_sleep: parameters related to sleep detection
- params_physact: parameters related to physical activity
- params_247: parameters related to 24/7 behaviours that do not fall into the typical sleep or physical activity research category.
- params_output: parameters relating to how and whether output is stored.
- params_general: general parameters not covered by any of the above categories

This structure was introduced in GGIR version 2.5-6 to make the GGIR code and documentation easier to navigate.

To see the parameters in each parameter category and their default values do:
```{R,eval=FALSE}
library(GGIR)
print(load_params())
```

If you are only interested in one specific category like sleep:
```{R,eval=FALSE}
library(GGIR)
print(load_params()$params_sleep)
```


If you are only interested in parameter "HASIB.algo" from the sleep_params object:
```{R,eval=TRUE}
library(GGIR)
print(load_params()$params_sleep[["HASPT.algo"]])
```



Documentation for this parameter objects can be found in the ([GGIR pdf manual](https://CRAN.R-project.org/package=GGIR/GGIR.pdf)). All of these are accepted as argument to function `g.shell.GGIR`, because `g.shell.GGIR` is a shell around all GGIR functionality. However, the `params_` objects themselves can not be provided as input to `g.shell.GGIR`.


### Key general arguments

You will probably never need to think about most of the arguments listed above, because a lot of arguments are only included to facilitate methodological studies where researchers want to have control over every little detail.

The bare minimum input needed for `g.shell.GGIR` is:

```{R,eval=FALSE}
library(GGIR)
g.shell.GGIR(datadir="C:/mystudy/mydata",
             outputdir="D:/myresults")
```

Argument `datadir` allows you to specify where you have stored your accelerometer data and `outputdir` allows you to specify where you would like the output of the analyses to be stored. This cannot be equal to `datadir`. If you copy paste the above code to a new R script (file ending with .R) and Source it in R(Studio) then the dataset will be processed and the output will be stored in the specified output directory.

Below we have highlighted the key arguments you may want to be aware of. We are not giving a detailed explanation, please see the package manual for that.

- `mode` - which part of GGIR to run, GGIR is constructed in five parts.
- `overwrite` - whether to overwrite previously produced milestone output. Between each GGIR part, GGIR stores milestone output to ease re-running parts of the pipeline.
- `idloc` - tells GGIR where to find the participant ID (default: inside file header)
- `strategy` - informs GGIR how to consider the design of the experiment.
  - If `strategy` is set to value 1, then check out arguments `hrs.del.start` and `hrs.del.end`.
  - If `strategy` is set to value 3, then check out arguments `ndayswindow`.
- `maxdur` - maximum number of days you expect in a data file based on the study protocol.
- `desiredtz` - time zone of the experiment.
- `chunksize` - a way to tell GGIR to use less memory, which can be useful on machines with limited memory.
- `includedaycrit` - tell GGIR how many hours of valid data per day (midnight-midnight) is acceptable.
- `includenightcrit` - tell GGIR how many hours of a valid night (noon-noon) is acceptable.
- `qwindow` - argument to tell GGIR whether and how to segment the day for day-segment specific analysis.
- `mvpathreshold` and `boutcriter` - acceleration threshold and bout criteria used for calculating time spent in MVPA (only used in GGIR part2).
- `epochvalues2csv` - to export epoch level magnitude of acceleration to a csv files (in addition to already being stored as RData file)
- `dayborder` - to decide whether the edge of a day should be other than midnight.
- `iglevels` - argument related to intensity gradient method proposed by A. Rowlands.
- `do.report` - specify reports that need to be generated.
- `viewingwindow` and `visualreport` - to create a visual report, this only works when all five parts of GGIR have successfully run.

### Arguments from the individual parts and in which parameter object they are stored

The table below shows arguments in GGIR, the GGIR part they are used in, and the parameter object they belong too.
Note that a few parameters are not part of an parameter object.

Argument (parameter) | Used in GGIR part | Parameter object
-------------------- | ------------ | ----------------------
datadir | 1, 2, 4, 5 | not in parameter objects
f0 | 1, 2, 3, ,4 ,5 | not in parameter objects
f1 | 1, 2, 3, 4, 5 | not in parameter objects
windowsizes | 1, 5 | params_general
desiredtz | 1, 2, 3, 4, 5 | params_general
overwrite | 1, 2, 3, 4, 5 | params_general
do.parallel | 1, 2, 3, 5 | params_general
myfun | 1, 2, 3 | not in parameter objects
maxNcores | 1, 2, 3, 5 | params_general
outputdir | 1 | not in parameter objects
chunksize | 1 | params_rawdata
studyname | 1 | not in parameter objects
do.enmo | 1 | params_metrics
do.lfenmo | 1 | params_metrics
do.en | 1 | params_metrics
do.bfen | 1 | params_metrics
do.hfen | 1 | params_metrics
do.hfenplus | 1 | params_metrics
do.mad | 1 | params_metrics
do.anglex | 1 | params_metrics
do.angley | 1 | params_metrics
do.angle | 1 | params_metrics
do.enmoa | 1 | params_metrics
do.roll_med_acc_x | 1 | params_metrics
do.roll_med_acc_y | 1 | params_metrics
do.roll_med_acc_z | 1 | params_metrics
do.dev_roll_med_acc_x | 1 | params_metrics
do.dev_roll_med_acc_y | 1 | params_metrics
do.dev_roll_med_acc_z | 1 | params_metrics
do.lfen | 1 | params_metrics
do.lfx | 1 | params_metrics
do.lfy | 1 | params_metrics
do.lfz | 1 | params_metrics
do.hfx | 1 | params_metrics
do.hfy | 1 | params_metrics
do.hfz | 1 | params_metrics
do.bfx | 1 | params_metrics
do.bfy | 1 | params_metrics
do.bfz | 1 | params_metrics
do.zcx | 1 | params_metrics
do.zcy | 1 | params_metrics
do.zcz | 1 | params_metrics
lb | 1 | params_metrics
hb | 1 | params_metrics
n | 1 | params_metrics
do.cal | 1 | params_rawdata
spherecrit | 1 | params_rawdata
minloadcrit | 1 | params_rawdata
printsummary | 1 | params_rawdata
print.filename | 1 | params_general
backup.cal.coef | 1 | params_rawdata
*.rmc | 1 | params_raw
imputeTimegaps | 1 | params_rawdata
selectdaysfile | 1, 2 | params_cleaning
dayborder | 1, 2, 5 | params_general
dynrange | 1 | params_rawdata
configtz | 1 | params_general
minimumFileSizeMB | 1 | params_rawdata
interpolationType | 1 | params_rawdata
metadatadir | 2, 3, 4, 5 | not in parameter objects
minimum_MM_length.part5 | 5 | params_cleaning
strategy | 2, 5 | params_cleaning
hrs.del.start | 2, 5 | params_cleaning
hrs.del.end | 2, 5| params_cleaning
maxdur | 2, 5 | params_cleaning
includedaycrit | 2 | params_cleaning
L5M5window | 2 | params_247
M5L5res | 2, 5 | params_247
winhr | 2, 5| params_247
qwindow | 2 | params_247
qlevels | 2 | params_247
ilevels | 2 | params_247
mvpathreshold | 2 | params_phyact
boutcriter | 2 | params_phyact
ndayswindow | 2 | params_cleaning
idloc | 2, 4 | params_general
do.imp | 2 | params_cleaning
storefolderstructure | 2, 4, 5 | params_output
epochvalues2csv | 2 | params_output
mvpadur | 2 | params_phyact
window.summary.size | 2 | params_247
bout.metric | 2, 5| params_phyact
closedbout | 2 | params_phyact
IVIS_windowsize_minutes | 2 | params_247
IVIS_epochsize_seconds | 2 | params_247
IVIS.activity.metric | 2 | params_247
iglevels | 2, 5 | params_247
TimeSegments2ZeroFile | 2 | params_cleaning
qM5L5 | 2 | params_247
MX.ig.min.dur | 2 | params_247
qwindow_dateformat | 2 | params_247
anglethreshold | 3 | params_sleep
timethreshold | 3 | params_sleep
acc.metric | 3, 5 | params_general
ignorenonwear | 3 | params_sleep
constrain2range | 3 | params_sleep
do.part3.pdf | 3 | params_output
sensor.location | 3, 4 | params_general
HASPT.algo | 3 | params_sleep
HASIB.algo | 3 | params_sleep
Sadeh_axis | 3 | params_sleep
longitudinal_axis | 3 | params_sleep
HASPT.ignore.invalid | 3 | params_sleep
loglocation | 4, 5 | params_sleep
colid | 4 | params_sleep
coln1 | 4 | params_sleep
nnights | 4 | params_sleep
sleeplogidnum | 4, 5 | params_sleep
do.visual | 4 | params_output
outliers.only | 4 | params_output
excludefirstlast | 4 | params_cleaning
criterror | 4 | params_output
includenightcrit | 4 | params_cleaning
relyonguider | 4 | params_sleep
relyonsleeplog | 4 | not in parameter objects
def.noc.sleep | 4 | params_sleep
data_cleaning_file | 4, 5 | params_cleaning
excludefirst.part4 | 4 | params_cleaning
excludelast.part4 | 4 | params_cleaning
sleeplogsep | 4 | params_cleaning
sleepwindowType | 4 | params_cleaning
excludefirstlast.part5 | 5 | params_cleaning
boutcriter.mvpa | 5 | params_phyact
boutcriter.in | 5 | params_phyact
boutcriter.lig | 5 | params_phyact
threshold.lig | 5 | params_phyact
threshold.mod | 5 | params_phyact
threshold.vig | 5 | params_phyact
timewindow | 5 | params_output
boutdur.mvpa | 5 | params_phyact
boutdur.in | 5 | params_phyact
boutdur.lig | 5 | params_phyact
save_ms5rawlevels | 5 | params_output
part5_agg2_60seconds | 5 | params_general
save_ms5raw_format | 5 | params_output
save_ms5raw_without_invalid | 5 | params_output
includedaycrit.part5 | 5 | params_cleaning
frag.metrics | 5 | params_phyact
LUXthresholds | 5 | params_247
LUX_cal_constant | 5 | params_247
LUX_cal_exponent | 5 | params_247
LUX_day_segments | 5 | params_247
do.sibreport | 5 | params_output

### Key arguments related to sleep analysis

For an explanation on how sleep is detected and the specific role of the various function arguments see section [Sleep analysis](#Sleep_analysis).

- Arguments related to configuring the sleep detection algorithm: `anglethreshold`, `timethreshold`, `HASPT.algo`, `HASIB.algo`, `Sadeh_axis`, and `HASPT.ignore.invalid`.
- `ignorenonwear` if set to TRUE then ignore detected monitor non-wear periods in the detection of sustained inactivity bouts to avoid confusion between monitor non-wear time.
- If you want to create a visualisation of how sleep period time and sustained inactivity bouts match throughout a day then consider arguments `do.visual`, `outliers.only`, and `criterror`.
- If you want to exclude the first and last night from the sleep analysis then used `excludefirstlast`.
- `def.noc.sleep` specifies how the sleep period time window should be estimated if no sleeplog is used.
- `includenightcrit` Minimum number of valid hours per night (24 hour window between noon and noon or 6pm-6pm).
- `data_cleaning_file` to ginore specific nights for specific individuals, see also section [Data cleaning file](#Data_cleaning_file).
- If you want the sleep analysis to be guided by a sleeplog (diary) then check out arguments `loglocation` which specifies the location of the spreadsheet (csv) with sleep log information. Further, use arguments `colid`, `coln1`, `nnights`, and `sleeplogsep` to specify the details of your sleep log structure. 

GGIR facilitates two possible sleeplog file structures:

#### Basic sleep log

Example of a basic sleeplog:

ID | onset_N1 | wakeup_N1 | onset_N2 | wakeup_N2 | onset_N3 | wakeup_N3 | onset_N4 | ...
---|---|---|---|---|---|---|---|---
123 | 22:00:00 | 09:00:00 | 21:00:00 | 08:00:00 | 23:45:00 | 06:30:00 | 00:00:00 | ...
345 | 21:55:00 | 08:47:00 |  |  | 23:45:00 | 06:30:00 | 00:00:00 | ...

- One column for participant id, this does not have to be the first column. Specify which column it is with  argument `colid`.
- Alternatingly one column for onset time and one column for waking time. Specify which column is the column for the first night by argument `coln1`, in the above example `coln1=2`.
- Timestamps are to be stored without date as in hh:mm:ss with hour values ranging between 0 and 23 (not 24). If onset corresponds to lights out or intention to fall asleep, then it specify `sleepwindowType = "TimeInBed"`.
- There can be multiple sleeplogs in the same spreadsheet. Each row representing a single recording. 
- First row: The first row of the spreadsheet needs to be filled with column names. For the basic sleep log format it does not matter what these column names are.


#### Advanced sleep log

Example of an advanced sleeplog for two recordings:

ID | D1_date | D1_wakeup | D1_inbed | D1_nap_start | D1_nap_end | D1_nonwear1_off | D1_nonwear1_on | D2_date | ...
---|---|---|---|---|---|---|---|---|---
123 | 30/03/2015 |  09:00:00 | 22:00:00 | 11:15:00 | 11:45:00 | 13:35:00 | 14:10:00 | 31/03/2015 | ...
567 | 20/04/2015 |  08:30:00 | 23:15:00 |          |          |          |          | 21/04/2015 | ...

Relative to the basic sleeplog format the advanced sleep log format comes with the following changes:

- Recording are stored in rows, while all information per days are stored in columns.
- Information per day is preceded by one columns that holds the calendar date. GGIR has been designed to recognise and handle any date format but assumes that you have used the same date format consistently through the sleeplog.
- Per calendar date there is a column for wakeup time and followed by a column for onset or in-bed time. Note that this is different from the basic sleep log, where wakeup time follows the column for onset or in-bed time. So, the advanced sleep log is calendar date oriented. However, if the sleep onset time is at 2am, you should still fill in the 02:00:00, even thought it is the 02:00:00 of the next calendar date.
- You can add columns relating to napping time and nonwear time. These are not used for the sleep analysis in g.part3 and g.part4, but used in g.part5 to facilitate napping analysis, see argument `do.sibreport`. Multiple naps and multiple nonwear periods can be entered per day.
- Leave cells for missing values blank.
- Column names are critical for the advanced sleeplog format: Date columns are recognised by GGIR as any column name with the word "date" in it. The advanced sleep log format is recognised by GGIR by looking out for the occurrence of at least two columnnames with the word "date" in their name. Wakeup times are recognised with the words "wakeup" in the column name. Sleeponset, to bed or in bed times are recognised as the columns with any of the following character combinations in their name: "onset", "inbed", "tobed", or "lightsout". Napping times are recognised by columns with the word "nap" in their name. Nonwear times are recognised by columns with the word "nonwear" in their name.


### Key arguments related to time use analysis

For an explanation on how time use analysis is performed see section [Waking-waking or 24 hour time-use analysis](#Waking-waking_or_24_hour).

- `excludefirstlast.part5` - whether to ignore the last and first day.
- `includedaycrit.part5` - tell GGIR what fraction of the waking hours in a day (value below 1) is acceptable.
- `minimum_MM_length.part5` - tell GGIR what the minimum length (hours) should be of the MM window in part 5.
- `bout.metric` - choose metric for calculating bouts (we recommend using the default setting).
- Configure thresholds for acceleration levels (some may want to interpret this as intensity levels): `threshold.lig`, `threshold.mod`, `threshold.vig`.
- Configure what fraction of a bout needs to meet the threshold (cut-point) crtieria `boutcriter.in`, `boutcriter.lig`, `boutcriter.mvpa`.  Note that bout.metric and the boutcriter arguments are complementary.
When bout.metric = 6 combined with boutcriter.mvpa=0.8 means that an MVPA bout can have interruptions (i.e., the time out of MVPA intensity) that meet the following criteria:
  (1) A single interruption can last < 1 min
  (2) Repeated interruptions are allowed provided that their total time does not exceed 20% of the bout duration
  (3) The time spent in the interruptions is included in the duration of the MVPA bout. For example: A 25-minute bout can have two 1 minute interruption, but not a single 2-minute interruption. Here, the full 25 minutes would count towards the duration of the MVPA bout.
- `timewindow` to specify whether days should be defined from midnight to midnight `"MM"`, from waking-up to waking-up `"WW"`, or both `c("MM","WW")`.
- Configure durations of bouts: `boutdur.mvpa`, `boutdur.in`, and `boutdur.lig`. Note that this can be a vector of multiple values indicating the minimum and maximum duration of subsequent bout types, e.g. 1-5 minutes MVPA, 5-10 minutes MVPA, and longer than 10 minutes MVPA.

### Published cut-points and how to use them

Cut-points to estimate time spent in acceleration levels that are roughly liked to levels of energy metabolism have been proposed by:

- [Esliger et al 2011](https://journals.lww.com/acsm-msse/Fulltext/2011/06000/Validation_of_the_GENEA_Accelerometer.22.aspx): wrist and waist in adults.
- [Phillips et al 2013](https://www.jsams.org/article/S1440-2440(12)00112-0/fulltext): wrist and hip in 8-14 year olds.
- [Schaefer et al 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3960318/): wrist in 6-11 year old children.
- [Hildebrand et al 2014](https://journals.lww.com/acsm-msse/Fulltext/2014/09000/Age_Group_Comparability_of_Raw_Accelerometer.17.aspx) and [2016](https://onlinelibrary.wiley.com/doi/abs/10.1111/sms.12795): wrist and hip in - [Vaha-Ypya et al 2015](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0134813): hip in adults.
- [Roscoe et al 2017](https://link.springer.com/article/10.1007/s00431-017-2948-2): wrist in 4-5 year old pre-school children.
7-11 and 21-61 years old.
- [Sanders et al. 2018](https://doi.org/10.1080/02640414.2018.1555904): wrist and hip in older adults (mean: 69.6 years)
- [Migueles et al  2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8150960/): wrist and hip in older adults aged above 70 (mean: 78.7) years. 
- [Fraysse et al. 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7843957/): wrist in older adults aged above 70 (mean: 77) years.
- [Dibben et al. 2020](https://bmcsportsscimedrehabil.biomedcentral.com/articles/10.1186/s13102-020-00196-7): wrist- and hip in adults with heart failure (mean age 71)
- If you are aware of any additional publications then let us know.

**Acceleration metric not available in GGIR?**
Some of the above publications make use of acceleration metrics that sum their values per epoch rather than average them per epoch like GGIR does. However, no need to worry: Those metrics are effectively the same as the metrics included in GGIR but sum rather than average their raw values per epoch. So, to use their cut-point we need to multiply the proposed cut-point by the sample frequency used in the study that proposed it. For each of the studies this is detailed below. Note that GGIR intentionally does not sum values per epoch because that approach makes the cut-point sample frequency dependent, which complicates comparisons and harmonisation of literature. The explained variance and accuracy remains identical because we are only multiplying with a constant.

**Warning about validity in new studies**
In all of the studies above, excluding Hildebrand et al. 2016, no effort was made to calibrate the acceleration sensors relative to gravitational acceleration prior to cut-point development. Theoretically this can be expected to cause a bias in the cut-point estimates proportional to the calibration error in each device, especially for cut-points based on acceleration metrics which rely on the assumption of accurate calibration such as metrics: ENMO, EN, ENMOa, and by that also metric SVMgs used by studies such as Esliger 2011, Phillips 2013, and Dibben 2020.

Studies done with ActiGraph devices when configured with 'idle sleep mode' turned on, will have zero-strings in all three axes during periods of no movement. Studies do not clarify how these zeros strings are accounted for. The insertion of zero strings is problematic as raw data accelerometers should always measure the gravitational component when not moving. This mistake by the ActiGraph company directly affect metrics that rely on the presence of a gravitational component such as ENMO, EN, ENMOAa, and SVMgs. Further, also other metrics may be affected as the sudden disappearance of gravitational acceleration will result in a spike at the start and end of the no movement period.

Several studies that aimed to independently evaluate cut-point methods seem to be failing in recognise these challenges. Therefore, it is best to interpret cut-points with caution and just as a crude boundary between intensity levels. Future methodological studies around cut-points are advised to account for accelerometer calibration error and the problematic zero strings in raw ActiGraph data.

**Esliger 2011, Phillips 2013, Fraysse 2020, Dibben2020:**

-	In GGIR use metric ENMOa instead of ENMO with arguments `do.enmoa = TRUE`, `do.enmo = FALSE`, and `acc.metric=”ENMOa”`.
- `threshold.lig = (LightCutPointFromPaper_in_gmins/(sampleRateInStudy*60)) * 1000`
- `threshold.mod = (ModerateCutPointFromPaper_in_gmins/(sampleRateInStudy*60)) * 1000`
- `threshold.vig = (VigorousCutPointFromPaper_in_gmins/(sampleRateInStudy*60)) * 1000`
- `mvpathreshold = (ModerateCutPointFromPaper_in_gmins/(sampleRateInStudy*60)) * 1000`
- The value for `sampleRateInStudy` was 80 for Esliger and Phillips and 100 for Fraysse.
- [Dibben et al. 2020](https://bmcsportsscimedrehabil.biomedcentral.com/articles/10.1186/s13102-020-00196-7) does not report sample frequency by which the proposed cut-points cannot be replicated.
-	In the part 2 results you will need the MVPA estimates that are related to ENMOa, not ENMO.
-	In the part 5 results everything will be based on the new cut-points.

**Roscoe 2017:**

-	In GGIR use metric ENMOa instead of ENMO with arguments `do.enmoa = TRUE`, `do.enmo = FALSE`, and `acc.metric=”ENMOa”`.
- `threshold.lig = (LightCutPointFromPaper_in_gsecs/85.7) * 1000`
- `threshold.mod = (ModerateCutPointFromPaper_in_gsecs/85.7) * 1000`
- `threshold.vig = (VigorousCutPointFromPaper_in_gsecs/85.7) * 1000`
- `mvpathreshold = (ModerateCutPointFromPaper_in_gsecs/85.7) * 1000`
-	In the part 2 results you will need the MVPA estimates that are related to ENMOa, not ENMO.
-	In the part 5 results everything will be based on the new cut-points.

**Schaeffer 2014:**

-	In GGIR use metric EN instead of ENMO with arguments `do.en = TRUE`, `do.enmo = FALSE`, and `acc.metric=”EN”`.
-	Specify Schaeffer cut-points as:
- `threshold.lig = (LightCutPointFromPaper/75) * 1000`
- `threshold.mod = (ModerateCutPointFromPaper/75) * 1000`
- `threshold.vig = (VigorousCutPointFromPaper/75) * 1000`
- `mvpathreshold = (ModerateCutPointFromPaper/75) * 1000`
-	In the part2 results you will need the MVPA estimates that are related to EN, not ENMO.
-	In the part 5 results everything will be based on the new cut-points.

**Vaha-Ypya et al 2015:**

-	Use default setting do.mad= TRUE, acc.metric=”MAD”
-	Use the cut-points as provided by Vaha-Ypya directly. No need for scaling.

**Hildebrand 2014, Hildebrand 2016, Migueles 2021. Sanders 2018:**

-	Use default setting do.enmo= TRUE, acc.metric=”ENMO”
-	Use the cut-points as provided by Hildebrand directly. No need for scaling.


### Example call

If you consider all the arguments above you me may end up with a call to `g.shell.GGIR` that could look as follows.

```{R,eval=FALSE}
library(GGIR)
g.shell.GGIR(
             mode=c(1,2,3,4,5),
             datadir="C:/mystudy/mydata",
             outputdir="D:/myresults",
             do.report=c(2,4,5),
             #=====================
             # Part 2
             #=====================
             strategy = 1,
             hrs.del.start = 0,          hrs.del.end = 0,
             maxdur = 9,                 includedaycrit = 16,
             qwindow=c(0,24),
             mvpathreshold =c(100),
             bout.metric = 4,
             excludefirstlast = FALSE,
             includenightcrit = 16,
             #=====================
             # Part 3 + 4
             #=====================
             def.noc.sleep = 1,
             outliers.only = TRUE,
             criterror = 4,
             do.visual = TRUE,
             #=====================
             # Part 5
             #=====================
             threshold.lig = c(30), threshold.mod = c(100),  threshold.vig = c(400),
             boutcriter = 0.8,      boutcriter.in = 0.9,     boutcriter.lig = 0.8,
             boutcriter.mvpa = 0.8, boutdur.in = c(1,10,30), boutdur.lig = c(1,10),
             boutdur.mvpa = c(1),
             includedaycrit.part5 = 2/3,
             #=====================
             # Visual report
             #=====================
             timewindow = c("WW"),
             visualreport=TRUE)
```

Once you have used `g.shell.GGIR` and the output directory (outputdir) will be filled with milestone data and results.

### Configuration file

Function `g.shell.GGIR` stores all the explicitly entered argument values and default values for the argument that are not explicitly provided in a csv-file named config.csv stored in the root of the output folder. The config.csv file is accepted as input to `g.shell.GGIR`
with argument `configfile` to replace the specification of all the arguments, except `datadir` and `outputdir`, see example below.

```{R,eval=FALSE}
library(GGIR)
g.shell.GGIR(datadir="C:/mystudy/mydata",
             outputdir="D:/myresults", configfile = "D:/myconfigfiles/config.csv")
```

The practical value of this is that it eases the replication of analysis, because instead of having to share you R script, sharing your config.csv file will be sufficient. Further, the config.csv file contribute to the reproducibility of your data analysis.

Note 1: When combining a configuration file with explicitly provided argument values, the explicitly provided argument values will overrule the argument values in the configuration file.
Note 2: The config.csv file in the root of the output folder will be overwritten every time you use `g.shell.GGIR`. So, if you would like to add annotations in the file, e.g. in the fourth column, then you will need to store it somewhere outside the output folder and explicitly point to it with `configfile` argument.

# Time for action: How to run your analysis?

## From the console

You can use

`source("pathtoscript/myshellscript.R")`

or use the Source button in RStudio if you use RStudio.

## In a cluster

GGIR by default support multi-thread processing, which can be turned off by seting argument `do.parallel = FALSE`. If this is still not fast enough then I advise using a GGIR on a computing cluster. The way I did it on a Sun Grid Engine cluster is shown below, please note that some of these commands are specific to the computing cluster you are working on. Also, you may actually want to use an R package like clustermq or snowfall, which avoids having to write bash script. Please consult your local cluster specialist to tailor this to your situation. In my case, I had three files for the SGE setting:

**submit.sh**

```{bash,eval=FALSE}
for i in {1..707}; do
    n=1
    s=$(($(($n * $[$i-1]))+1))
    e=$(($i * $n))
    qsub /home/nvhv/WORKING_DATA/bashscripts/run-mainscript.sh $s $e
done
```

**run-mainscript.sh**

```{bash,eval=FALSE}
#! /bin/bash
#$ -cwd -V
#$ -l h_vmem=12G
/usr/bin/R --vanilla --args f0=$1 f1=$2 < /home/nvhv/WORKING_DATA/test/myshellscript.R
```

**myshellscript.R**

```{R,eval=FALSE}
options(echo=TRUE)
args = commandArgs(TRUE)
if(length(args) > 0) {
  for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
  }
}
g.shell.GGIR(f0=f0,f1=f1,...)
```

You will need to update the `...` in the last line with the arguments you used for `g.shell.GGIR`. Note that `f0=f0,f1=f1` is essential for this to work. The values of `f0` and `f1` are passed on from the bash script.

Once this is all setup you will need to call `bash submit.sh`
from the command line.

**Important Note:**

Please make sure that you process one GGIR part at the same time on a cluster, because each part assumes that preceding parts have been ran. You can make sure of this by always specifying argument `mode` to a single part of GGIR. Once the analysis stops update argument `mode` to the next part until all parts are done. The speed of the parallel processing is obviously dependent on the capacity of your computing cluster and the size of your dataset.


# Inspecting the results

GGIR generates the following types of output.
-	csv-spreadsheets with all the variables you need for physical activity, sleep and circadian rhythm research
-	Pdfs with on each page a low resolution plot of the data per file and quality indicators
- R objects with milestone data
-	Pdfs with a visual summary of the physical activity and sleep patterns as identified (see example below)

```{r, out.width = "700px",echo=FALSE}
knitr::include_graphics("reportexample.jpg")
```

## Output part 2 {.tabset}

Part 2 generates the following output:

- part2_summary.csv: Person level summary (see below)
- part2_daysummary.csv: Day level summary (see below)
- QC/data_quality_report.csv: Overview of calibration results and whether or not a file was corrupt or too short to be processed,
- QC/plots to check data quality 1.pdf: A pdf with visualisation of the acceleration time series in 15 minute resolution and with invalid data segments highlighted in colours (yellow: non-wear based on standard deviation threshold, brown: non-wear after extra filtering step (introduced in 2013), and purple: clipping)

### Person level summary

(Part of) variable name   |   Description
------------------------ | ------------------------------------------------------
ID | Participant id
device_sn | Device serial number
bodylocation | Body location extracted from file header
filename | Name of the data file
start_time | Timestamp when recording started
startday |	Day of the week on which recording started
samplefreq |	Sample frequency (Hz)
device | Accelerometer brand, e.g. GENEACtiv
clipping_score | The Clipping score: Fraction of 15 minute windows per file for which the acceleration in one of the three axis was close to the maximum for at least 80% of the time. This should be 0.
meas_dur_dys} | Measurement duration (days)
complete_24hcycle | Completeness score: Fraction of 15 minute windows per 24 hours for which no valid data is available at any day of the measurement.
meas_dur_def_proto_day | measurement duration according to protocol (days):	Measurement duration (days) minus the hours that are ignored at the beginning and end of the measurement motivated by protocol design
wear_dur_def_proto_day | wear duration duration according to protocol (days): So, if the protocol was seven days of measurement, then wearing the accelerometer for 8 days and recording data for 8 days will still make that the wear duration is 7 days
calib_err | Calibration error (static estimate)	Estimated based on all ‘non-movement’ periods in the measurement after applying the autocalibration.
calib_status | Calibration status: Summary statement about the status of the calibration error minimisation
ENMO_fullRecordingMean | ENMO is the main summary measure of acceleration. The value presented is the average ENMO over all the available data normalised per 24-hour cycles (diurnal balanced), with invalid data imputed by the average at similar time points on different days of the week. In addition to ENMO it is possible to extract other acceleration metrics (i.e. BFEN, HFEN, HFENplus). We emphasize that it is calculated over the full recording because the alternative is that a variable is only calculated overmeasurement days with sufficient valid hours of data.
ENMO | (only available if set to true in part1.R)	ENMO is the main summary measure of acceleration. The value presented is the average ENMO over all the available data normalised per 24 hour cycles, with invalid data imputed by the average at similar timepoints on different days of the week. In addition to ENMO it is possible to extract other acceleration metrics in part1.R (i.e. BFEN, HFEN, HFENplus)	See also [van Hees PLoSONE April 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061691) for a detailed description and comparison of these techniques.
pX_A_mg_0-24h_fullRecording |	This variable represents the Xth percentile in the distribution of short epoch metric value A of the average day. The average day may not be ideal for describing the distribution. Therefore, the code also extracts the following variable.
AD_pX_A_mg_0-24h |	This variable represents the Xth percentile in the distribution of short epoch metric value A per day averaged across all days.
L5_A_mg_0-24 | Average of metric A during the least active five* hours in the day that is the lowest rolling average value of metric A. (* window size is modifiable by argument `winhr`)
M5_A_mg_0-24 |  Average of metric A during the most active five* hours in the day that is the lowest rolling average value of metric A. (* window size is modifiable by argument `winhr`)
L5hr_A_mg_0-24 | Starting time in hours and fractions of hours of L5_A_mg_0-24
M5hr_A_mg_0-24 | Starting time in hours and fractions of hours of M5_A_mg_0-24
ig_gradient_ENMO_0-24hr_fullRecording | Intensity gradient calculated over the full recording.
1to6am_ENMO_mg | Average metric value ENMO between 1am and 6am
N valid WEdays | Number of valid weekend days
N valid WKdays | Number of valid week days
IS_interdailystability | inter daily stability. The movement count that is derived for this was an attempt to follow the original approach by Eus J. W. Van Someren (Chronobiology International. 1999. Volume 16, issue 4).
IV_intradailyvariability | intra daily variability. In contrast to the original paper, we ignore the epoch transitions between the end of a day and the beginning of the next day for the numerator of the equation, this to make it a true measure of intradaily variability. Same note as for IS: The movement count that is derived for this was an attempt to follow the original approach.
IVIS_windowsize_minutes | Sizes of the windows based on which IV and IS are calculated (note that this is modifiable)
IVIS_epochsize_seconds | Argument has been deprecated
`AD_` | All days (plain average of all available days, no weighting). The variable was calculated per day and then averaged over all the available days
`WE_` | Weekend days (plain average of all available days, no weighting). The variable was calculated per day and then averaged over weekend days only
`WD_` | Week days (plain average of all available days, no weighting). The variable was calculated per day and then averaged over week days only
`WWE_` | Weekend days (weighted average)	The variable was calculated per day and then averaged over weekend days. Double weekend days are averaged. This is only relevant for experiments that last for more than seven days.
`WWD_` | Week days (weighted average)	The variable was calculated per day and then averaged over week days. Double weekend days were averaged. This is only relevant for experiments that last for more than seven days)
WWD_MVPA_E5S_T100_ENMO | Time spent in moderate-to-vigorous based on 5 second epoch size and an ENMO metric threshold of 100
`WWE_MVPA_E5S_B1M80%_T100_ENMO` | Time spent in moderate-to-vigorous based on 5 second epoch size and an ENMO metric threshold of 100 based on a bout criteria of 100
`WE_[100,150)_mg_0-24h_ENMO` | Time spent between (and including) 100 mg and 150 (excluding 150 itself) between 0 and 24 hours (the full day) using metric ENMO data exclusion strategy (value=1, ignore specific hours; value=2, ignore all data before the first midnight and after the last midnight) |	A log of decision made in part2.R
`_MVPA_E5S_B1M80_T100` | MVPA calculated based on 5 second epoch setting bout duration 1 Minute and inclusion criterion of more than 80 percent. This is only done for metric ENMO at the moment, and only if mvpa threshold is not left blank
`_ENMO_mg` | ENMO or other metric was first calculated per day and then average according to AD, WD, WWE, WWD
data exclusion strategy | A log of the decision made when calling  g.impute: value=1 mean ignore specific hours; value=2 mean ignore all data before the first midnight and after the last midnight
n hours ignored at start of meas (if strategy=1) | number of hours ignored at the start of the measurement (if strategy = 1)	A log of decision made in part2.R
n hours ignored at end of meas (if strategy=1) | number of hours ignored at the end of the measurement (if strategy = 1).	A log of decision made in part2.R
n hours ignored at end of meas (if strategy=1) | number of days of measurement after which all data is ignored (if strategy = 1)	A log of decision made in part2.R
epoch size to which acceleration was averaged (seconds) | A log of decision made in part1.R
pdffilenumb |	Indicator of in which pdf-file the plot was stored
pdfpagecount |	Indicator of in which pdf-page the plot was stored


### Day level summary

This is a non-exhaustive list, because most concepts have been explained in summary.csv

(Part of) variable name	| Description
------------------------ | ---------------------------------------------
ID | Participant id
filename | Name of the data file
calender_date | Timestamp and date on which measurement started
bodylocation | Location of the accelerometer as extracted from file header
N valid hours | Number of hours with valid data in the day
N hours	| Number of hours of measurement in a day, which typically is 24, unless it is a day on which the clock changes (DST) resulting in 23 or 25 hours. The value can be less than 23 if the measurement started or ended this day
weekday | Name of weekday
measurement | Day of measurement	Day number relative to start of the measurement
L5hr_ENMO_mg_0-24h | Hour on which L5 starts for these 24 hours (defined with metric ENMO)
L5_ENMO_mg_0-24h | Average acceleration for L5 (defined with metric ENMO)
`[A,B)_mg_0-24h_ENMO` |	Time spent in minutes between (and including) acceleration value A in mg and (excluding) acceleration value B in mg based on metric ENMO
ig_gradient_ENMO_0-24hr | Gradient from intensity gradient analysis proposed by [Rowlands et al. 2018](https://journals.lww.com/acsm-msse/Fulltext/2018/06000/Beyond_Cut_Points__Accelerometer_Metrics_that.25.aspx) based on metric ENMO for the time segment 0 to 24 hours
ig_intercept_ENMO_0-24hr | Intercept from intensity gradient analysis proposed by [Rowlands et al. 2018](https://journals.lww.com/acsm-msse/Fulltext/2018/06000/Beyond_Cut_Points__Accelerometer_Metrics_that.25.aspx)  based on metric ENMO for the time segment 0 to 24 hours
ig_rsquared_ENMO_0-24hr | r squared from intensity gradient analysis proposed by [Rowlands et al. 2018](https://journals.lww.com/acsm-msse/Fulltext/2018/06000/Beyond_Cut_Points__Accelerometer_Metrics_that.25.aspx)  based on metric ENMO for the time segment 0 to 24 hours


## Output part 4 {.tabset}

Part 4 generates the following output:


### Night level summaries

- part4_nightsummary_sleep_cleaned.csv
- QC/part4_nightsummary_sleep_full.csv

The csv. files contain the variables as shown below.

(Part of) variable name |   Description
------------------------ | ------------------------------------------------------
ID | Participant ID extracted from file
night | Number of the night in the recording
sleeponset | Detected onset of sleep expressed as hours since the midnight of the previous night.
wakeup | Detected waking time (after sleep period) expressed as hours since the midnight of the previous night.
SptDuration | Difference between onset and waking time.
sleepparam | Definition of sustained inactivity by accelerometer.
guider | guider used as discussed in paragraph [Sleep analysis](#Sleep_analysis).
guider_onset | Start of Sleep Period Time window derived from the guider.
guider_wake | End of Sleep Period Time window derived guider.
guider_SptDuration | Time SPT duration derived from guider_wake and guider_onset.
error_onset |  Difference between sleeponset and guider_onset
error_wake | Difference between wakeup and guider_wake
fraction_night_invalid | Fraction of the night (noon-noon or 6pm-6pm) for which the data was invalid, e.g. monitor not worn or no accelerometer measurement started/ended within the night.
SleepDurationInSpt | Total sleep duration, which equals the accumulated nocturnal sustained inactivity bouts within the Sleep Period Time.
duration_sib_wakinghours | Accumulated sustained inactivity bouts during the day. These are the periods we would label during the night as sleep, but during the day they form a subclass of inactivity, which may represent day time sleep or wakefulness while being motionless for a sustained period of time number_sib_sleepperiod} Number of nocturnal sleep periods, with nocturnal referring to the Sleep Period Time window.
duration_sib_wakinghours_atleast15min | Same as duration_sib_wakinghours, but limited to SIBs that last at least 15 minutes.
number_sib_wakinghours | Number of sustained inactivity bouts during the day, with day referring to the time outside the Sleep Period Time window.
sleeponset_ts | sleeponset formatted as a timestamp
wakeup_ts | wakeup formatted as a timestamp
guider_onset_ts | guider_onset formatted as a timestamp
guider_wake_ts | guider_wake formatted as a timestamp
page | pdf page on which the visualisation can be found
daysleeper | If 0 then the person is a nightsleeper (sleep period did not overlap with noon) if value=1 then the person is a daysleeper (sleep period did overlap with noon)
weekday | Day of the week on which the night started
calendardate | Calendar date on which the night started in day/month/year format.
filename | Name of the accelerometer file
cleaningcode | see paragraph [Cleaningcode](#Cleaningcode).
sleeplog_used | Whether a sleep log was used (TRUE/FALSE)
acc_available | Whether accelerometer data was available (TRUE/FALSE).
WASO | Wake After Sleep Onset: SptDuration - SleepDurationInSpt
SptDuration | Sleep Period Time window duration: wakeup - sleeponset
error_onset | Difference between sleeponset and guider_onset (this variable is only available in the full report as stored in the QC folder)
error_wake | Difference between wakeup and guider_wake (this variable is only available in the full report as stored in the QC folder)
SleepRegularityIndex | The Sleep Regularity Index as proposed by [Phillips et al. 2017](https://www.nature.com/articles/s41598-017-03171-4), but calculated per day-pair to 
enable user to study patterns across days.
SriFractionValid | Fraction of the 24 hour period that was valid in both current as well
as in matching timestamps for the next calendar day. See GGIR function manual for details.


#### Non-default variables in part 4 csv report

These additional are only stored if you used a sleeplog that captures time in bed, 
or when using guider HorAngle for hip-worn accelerometer data. If either of these 
applies set argument `sleepwindowType` to "TimeInBed".

(Part of) variable name |   Description
------------------------ | ------------------------------------------------------
guider_guider_inbedStart | Time of getting in bed
guider_guider_inbedEnd | Time of getting out of bed
guider_inbedDuration | Time in Bed: guider_inbedEnd - guider_inbedStart
sleepefficiency | Sleep efficiency, calculated as: SleepDurationInSpt / guider_inbedDuration
sleeplatency | Sleep latency, calculated as: sleeponset - guider_inbedStart


### Person level summaries

- part4_summary_sleep_cleaned.csv
- QC/part4_summary_sleep_full.csv

In the person level report the variables are derived from the variables in the night
level summary. Minor extensions to the variable names explain how variables are aggregated across the days.
Please find below extra clarification on a few of the variable names for which the meaning may not be obvious:

(Part of) variable name | Description
------------------------ | ------------------------------------------------------
`_mn` | mean across days
`_sd` | standard deviation across days
`_AD` | All days
`_WE` | Weekend days
`_WD` | Week days
sleeplog_used | Whether a sleeplog was available (TRUE) or not (FALSE)
sleep_efficiency | Accelerometer detrive sleep efficiency within the sleep period time calculated as the ratio between acc_SleepDurationInSpt and acc_SptDuration (denominator). Only available at person level, because at night level the user can calculate this from existing variables.
n_nights_acc | Number of nights of accelerometer data
n_nights_sleeplog | Number of nights of sleeplog data.
n_WE_nights_complete | Number of weekend nights complete which means both accelerometer and estimate from guider.
n_WD_nights_complete | Number of weekday nights complete which means both accelerometer and estimate from guider.
n_WEnights_daysleeper | Number of weekend nights on which the person slept until after noon.
n_WDnights_daysleeper | Number of weekday nights on which the person slept until after noon.
duration_sib_wakinghour | Total duration of sustained inactivity bouts during the waking hours.
number_sib_wakinghours | Number of sustained inactivity bouts during the waking hours.
average_dur_sib_wakinghours | Average duration of the sustained inactivity bouts during the day (outside the sleep period duration). Calculated as duration_sib_wakinghour divided by number_sib_wakinghours per day, after which the mean and standard deviation are calculated across days.

### visualisation_sleep.pdf
Visualisation to support data quality checks:
- visualisation_sleep.pdf (optional)

When input argument `do.visual` is set to TRUE GGIR can show the following visual comparison between the time window of being asleep (or in bed) according to the sleeplog and the detected sustained inactivity bouts according to the accelerometer data. This visualisation is stored in the results folder as `visualisation_sleep.pdf`.

Explanation of the image: Each line represents one night. Colours are used to distinguish definitions of sustained inactivity bouts (2 definitions in this case) and to indicate existence or absence of overlap with the sleeplog. When argument `outliers.only` is set to FALSE it will visualise all available nights in the dataset. If `outliers.only` is set to TRUE it will visualise only nights with a difference in onset or waking time between sleeplog and sustained inactivity bouts larger than the value of argument `criterror`.

This visualisation with outliers.only set to TRUE and critererror set to 4 was very powerful to identify entry errors in sleeplog data in van Hees et al PLoSONE 2015. We had over 25 thousand nights of data, and this visualisation allowed us to quickly zoom in on the most problematic nights to investigate possible mistakes in GGIR or mistakes in data entry.

```{r, out.width = "700px",echo=FALSE}
knitr::include_graphics("example_dovisual.jpg")
```

## Output part 5 {.tabset}

The output of part 5 is dependent on the parameter configuration, it will generate as many output files as there are unique combination of the three thresholds provide. For example, the output could be:

For example, the following files will be generated if the threshold configuration was 30 for light activity, 100 for moderate and 400 for vigorous activity:
- part5_daysummary_MM_L30M100V400_T5A5.csv
- part5_daysummary_WW_L30M100V400_T5A5.csv
- part5_personsummary_MM_L30M100V400_T5A5.csv
- part5_personsummary_WW_L30M100V400_T5A5.csv
- file summary reports/Report_nameofdatafile.pdf

### Day level summary

(Term in) variable name | Description
------------------------ | ---------------------------------------------
sleeponset | onset of sleep expressed in hours since the  midnight in the night preceding the night of interest, e.g. 26 is 2am.
wakeup | waking up time express in the same way as sleeponset.
sleeponset_ts | onset of sleep expressed as a timestamp hours:minutes:seconds
daysleeper | if 0 then the person woke up before noon, if 1 then the person woke up after noon
cleaningcode | See paragraph [Cleaningcode](#Cleaningcode).
dur_day_spt_min | Total length of daytime waking hours and spt combined (typically 24 hours for MM report).
`dur_` | duration of a behavioral class that will be specified int he rest of the variable name
`ACC_` | (average) acceleration according to default metric specific by acc.metric
`_spt_wake_` | Wakefulness within the Sleep period time window.
`_spt_sleep_` | Sleep within the Sleep period time window.
`_IN_` | Inactivity
`_LIG_` | Light activity
`_MOD_` | Moderate activity
`_VIG_` | Vigorous activity
`_MVPA_` | Moderate or Vigorous activity
`_unbt_` | Unbouted
`_bts_` | Bouts (also known as sojourns), which are segments that for which the acceleration is within a specified range for a specified fraction of the time.
`_bts_1_10_` | Bouts lasting at least 1 minute and less than 10 minutes (1 and 9.99 minutes are included, but 10 minutes is not).
Nblock | number of blocks of a certain behavioral class, not these are not bouts but a count of the number of times the behavioral class occurs without interruptions.
WW | in filename refers to analyses based on the timewindow from waking to waking up
MM | in filename refers to analyses done on windows between midnight and midnight
calendar_date | calendar date on which the window started in day/month/year format. So, for WW window this could mean that you have two windows starting on the same date.
weekday | weekday on which the window started. So, for WW window this could mean that you have two windows starting on the weekday.
`_total_IN` | total time spent in inactivity (no distinction between bouted or unbouted behavior, this is a simple count of the number of epochs that meet the threshold criteria.
`_total_LIG` | total time spent in light activity.
nonwear_perc_day | Non-wear percentage during the waking hours of this day.
nonwear_perc_spt | Non-wear percentage during the spt hours of this day.
nonwear_perc_day_spt | Non-wear percentage during the whole day, including waking and spt.
dur_day_min | Duration of waking hours within this day window
dur_spt_min | Duration of Sleep Period Time within this day window.
dur_day_spt_min | Duration this day window, including both waking hours and SPT.
sleep_efficiency | sleep_efficiency in part 5 is not the same as in part 4, but calculated as the percentage of sleep within the sleep period time window. The conventional approach is the approach used in part 4.
L5TIME	| Timing of least active 5hrs, expressed as timestamp in the day
M5TIME | Timing of most active 5hrs
L5TIME_num, M5TIME_num | Timing of least/most active 5hrs, expressed as hours in the day. Note that L5/M5 timing variables are difficult to average across days because 23:00 and 1:00 would average to noon and not to midnight. So, caution is needed when interpreting person averages.
L5VALUE	| Acceleration value for least active 5hrs
M5VALUE	| Acceleration value for most active 5hrs
`FRAG_` | All variables related to behavioural fragmentation analysis
`TP_` | Transition probability
PA2IN | Physical activity fragments followed by inactivity fragments
IN2PA | Physical inactivity fragments followed by activity fragments
Nfrag | Number of fragments
IN2LIPA | Inactivity fragments followed by LIPA
IN2MVPA | Inactivity fragments followed by MVPA
`mean_dur` | mean duration of a fragment category
`Gini_dur` | Gini index
`CoV_dur` | Coefficient of Variation
alpha | Power law exponent
`x0.5` | Derived from power law exponent alpha, see [Chastin et al. 2010](https://www.sciencedirect.com/science/article/abs/pii/S096663620900602X?via%3Dihub)
`W0.5` | Derived from power law exponent alpha, see [Chastin et al. 2010](https://www.sciencedirect.com/science/article/abs/pii/S096663620900602X?via%3Dihub)
nap_count | Total number of naps, only calculated when argument do.sibreport = TRUE, currently optimised for 3.5-year olds. See function documentation for function `g.part5.classifyNaps`.
nap_totalduration | Total nap duration, only calculated when argument do.sibreport = TRUE, currently optimised for 3.5-year old. See function documentation for function `g.part5.classifyNaps`.

**Special note if you are working on compositional data analysis:**

The duration of all `dur_` variables that have `_total_` in their name should add up to the total length of the waking hours in a day.
Similarly, the duration of all other `dur_` variables excluding the variables `_total_` in their name and excluding the variable with `dur_day_min`, `dur_spt_min`, and `dur_day_spt_min` should also add up to the length of the full day.

**Motivation for default boutcriter.in = 0.9:**

The idea is that if you allow for bouts of 30 minutes it would
not make sense to allow for breaks of 20 percent (6 minutes!) this is why I
used a more stringent criteria for the highest category. Please note that
you can change these criteria via arguments `boutcriter.mvpa`, `boutcriter.in`,
and `boutcriter.lig`.


### Person level summary

Most variables in the person level summary are derived from the day level summary, but extended with `_pla` to indicate that the variable was calculated as the plain average across all valid days. Variables extended with `_wei` represent the weighted average of across all days where weekend days always weighted 2/5 relative to the contribution of week days.

Variable name	| Description
------------------------ | ---------------------------------------------
Nvaliddays | Total number of valid days.
Nvaliddays_WD | Number of valid week days.
Nvaliddays_WE | Number of valid weekend days, where the days that start on Saturday or Sunday are considered weekend.
NcleaningcodeX | Number of days that had cleaning code X for the corresponding sleep analysis in part 4. In case of MM analysis this refers to the night at the end of the day.
Nvaliddays_AL10F_WD	| Number of valid week days with at least 10 fragments (5 inactivity or 5 inactive)
Nvaliddays_AL10F_WE	| Number of valid weekend days with at least 10 fragments (5 inactivity or 5 inactive)
`_wei` | weighted average of weekend and week days, using a 2/5 ratio, see above.
`_pla` | plain average of all days, see above

# Motivation and clarification

In this chapter we will try to collect motivations and clarification behind GGIR which may not have been clear from the existing publications.

## Reproducibilty of GGIR analyses

Some tips to increase reproducibility of your findings:

1. When you publish your findings, please remember to add the GGIR package version number. All of GGIR are archived by CRAN and available from the archive section on the package [website](https://CRAN.R-project.org/package=GGIR). GGIR has evolved over the years. To get a better understanding of how versions differ you should check the NEWS sections from the package [website](https://CRAN.R-project.org/package=GGIR)
2. Report how you configured the accelerometer
3. Report the study protocol and wear instructions given to the participants
4. Report GGIR version
5. Report how GGIR was used: Share the config.csv file or your R script
6. Report how you post-processed / cleaned GGIR output
7. Report how reported outcomes relate to the specific variable names in GGIR


## Auto-calibration

An acceleration sensor works on the principle that acceleration is captured mechanically and converted into an electrical signal. The relationship between the electrical signal and the acceleration is usually assumed to be linear, involving an offset and a gain factor. We shall refer to the establishment of the offset and gain factor as the sensor calibration procedure. Accelerometers are usually calibrated as part of the manufacturing process under non-movement conditions using the local gravitational acceleration as a reference. The manufacturer calibration can later be evaluated by holding each sensor axis parallel (up and down) or perpendicular to the direction of gravity; readings for each axis should be ±1 and 0 g, respectively. However, this procedure can be cumbersome in studies with a high throughput. Furthermore, such a calibration check will not be possible for data that have been collected in the past and for which the corresponding accelerometer device does not exist anymore. Techniques have been proposed that can check and correct for calibration error based on the collected triaxial accelerometer data in the participant's daily life without additional experiments, referred to as autocalibration. The general principle of these techniques is that a recording of acceleration is screened for nonmovement periods. Next, the moving average over the nonmovement periods is taken from each of the three orthogonal sensor axes and used to generate a three-dimensional ellipsoid representation that should ideally be a sphere with radius 1 g. Here, deviations between the radius of the three-dimensional ellipsoid and 1 g (ideal calibration) can then be used to derive correction factors for sensor axis-specific calibration error. This auto-calibration performed by GGIR uses this technique and a more detailed description and demonstration can be found in the published paper.

Reference:

- van Hees VT, Fang Z, Langford J, Assah F, Mohammad A, da Silva IC, Trenell MI, White T, Wareham NJ, Brage S. Autocalibration of accelerometer data for free-living physical activity assessment using local gravity and temperature: an evaluation on four continents. J Appl Physiol (1985). 2014 Oct 1;117(7):738-44.
PMID: 25103964 [link](https://journals.physiology.org/doi/full/10.1152/japplphysiol.00421.2014)

Key decisions to be made:

1.	Whether to apply auto-calibration or not (default and recommended setting is YES). You  can turn this off by setting argument `do.call` in `g.shell.GGIR` to `do.call=FALSE`.
2.	Other variables are probably best left in their default setting

Key output variables:

1. Variable value `cal.error.end` as stored in data_quality_report.csv or variable value calib_err in summary.csv. These should be less than 0.01 g (10mg).

## Non-wear detection
Accelerometer non-wear time is estimated on the basis of the standard deviation and the value range of the raw data from each accelerometer axis. Classification is done per 15 minute (or ws2) block and based on the characteristics of the 60 minute (or ws) window centred at these 15 minutes. A block is classified as non-wear time if the standard deviation of the 60 minute window is less than 13.0 mg ($1 mg = 0.00981 m·s^−2$) and the value range of the 60 minute window is less than 50 mg, for at least two out of the three accelerometer axes. The procedure for non-wear detection was modified in comparison to the procedure as applied in the 2011 PLoSONE publication [link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022922). Instead of 30-minute time windows 60-minute time windows were used to decrease the chance of accidently detecting short sedentary periods as non-wear time. The windows were overlapping (15 minute steps, window overlap of 45 minutes), which was done to improve the accuracy of detecting the boundaries of non-wear time as opposed to non-overlapping time windows. Inspection of unpublished data on non-wear classification by the algorithm as described in our published work indicated that the algorithm does not cope well with periods of monitor transportation per post. Here, long periods of non-wear are briefly interrupted by periods of movement, which are normally interpreted as monitor wear. Therefore, the algorithm was expanded with an additional stage in which the plausibility of “wear-periods” in-between non-wear periods is tested. Short periods of detected wear-time in-between longer periods of detected non-wear were classified as non-wear time based on the duration and the proportion of the duration relative to the bordering periods of detected non-wear-periods. The following criteria were derived from visual observation of various datasets using knowledge about study protocols. All detected wear-periods of less than six hours and less than 30% of the combined duration of their bordering non-wear periods were classified as non-wear. Additionally, all wear-periods of less than three hours and which formed less than 80% of their bordering non-wear periods were classified as non-wear. The motivation for selecting a relatively high criteria (< 30%) in combination with a long period (6hrs) and a low criteria (< 80%) in combination with a short period (3 hrs) was that long period are more likely to be actually related to monitor wear time. A visual model was created, see Figure 1. Here, units of time are presented in squares and marked grey if detected as non-wear time. Period C is detected as wear-time and borders to non-wear periods B and D, see Figure 1. If the length of C is less than six hours and C divided by the sum of B and D is less than 0.3 then the first criteria is met and block C is turned into a non-wear period.

```{r, out.width = "400px",echo=FALSE}
knitr::include_graphics("nonwearimage.jpg")
```

By visual inspection of >100 traces from a large observational study it turned out that applying this stage in three iterative stages allowed for improved classification of periods characterised by intermittent periods of non-wear and apparent wear. Further, an additional rule was introduced for the final 24 hours of each measurement. The final 24 hours are often considered the period in which the accelerometer is potentially taken off but moved because of transportation, e.g. by the mail service. All wear-periods in the final 24 hrs of each measurement shorter than three hours and preceded by at least one hour of non-wear time were classified as non-wear. Finally, if the measurement starts or ends with a period of less than three hours of wear followed by non-wear (any length) then this period of wear is classified as non-wear. These additional criteria for screening the beginning and end of the accelerometer file reflect the likelihood of movements that are involved when starting the accelerometer or downloading the data from the accelerometer.

Reference:

- van Hees VT, Gorzelniak L, Dean León EC, Eder M, Pias M, Taherian S, Ekelund U, Renström F, Franks PW, Horsch A, Brage S. Separating movement and gravity components in an acceleration signal and implications for the assessment of human daily physical activity. PLoS One. 2013 Apr 23 [link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061691).

Key decisions to be made:

1.	Size of windows
2.	Whether to utilize the non-wear detection

Key output variables:

1.	Raw classification
2.	Non-wear duration
3.	Non-wear duration taking into account the protocol

## Clipping score

The acceleration signal was screened for ‘clipping’. If more than 50% of the data points in a 15 minute time window are higher than 7.5g (close to the maximal dynamic range of this sensor) the corresponding time period is considered as potentially corrupt data, which may be explained by the sensor getting stuck at its extreme value.

Reference:

- van Hees VT, Gorzelniak L, Dean León EC, Eder M, Pias M, Taherian S, Ekelund U, Renström F, Franks PW, Horsch A, Brage S. Separating movement and gravity components in an acceleration signal and implications for the assessment of human daily physical activity. PLoS One. 2013 Apr 23 [link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061691)

## Why collapse information to epoch level?

Although many data points are collected we decide to only work with aggregated values (e.g. 1 or 5 second epochs) for the following reasons:

1.	Accelerometers are often used to describe patterns in metabolic energy expenditure. Metabolic energy expenditure is typically defined per breath or per minute (indirect calorimetry), per day (room calorimeter), or per multiple days (doubly labelled water method). In order to validate our methods against these reference standards we need to work with a similar time resolution.

2.	Collapsing the data to epoch summary measures helps to standardise for differences in sample frequency between studies.

3.	There is little evidence that the raw data is an accurate representation of body acceleration. All scientific evidence on the validity of accelerometer data has so far been based on epoch averages.

4.	Collapsing the data to epoch summary measures may help to average out different noise levels and make sensor brands more comparable.


## Sleep analysis {#Sleep_analysis}

In GGIR sleep analysis has been implemented in part 3 and 4. Sleep analysis comes at two levels: The identification of the main Sleep Period Time (SPT) window or the time in bed window (TIB), and the discrimination of sleep and wakefulness periods. The term sleep is somewhat controversial in the context of accelerometry, because accelerometer only capture lack of movement. To acknowledge this challenge GGIR refers to these classified sleep periods as **sustained inactivity bouts (abbreviated as SIB)**. 

Current, GGIR offers the user the choice to identify SIB period using any of the following algorithms:

- **vanHees2015:** Heuristic algorithm proposed in 2015 [link](https://doi.org/10.1371/journal.pone.0142533) which looks for periods of time where the z-angle does not change by more than 5 degrees for at least 5 minutes. This in contrast to conventional sleep detection algorithms such as Sadeh and Galland which rely on acceleration to estimate sleep. The vanHees2015 algorithm is the default.
- **Sadeh1994:** The algorithm proposed by Sadeh et al. [link](https://doi.org/10.1093/sleep/17.3.201). 
To use the GGIR implementation of the zero-crossing counts and Sadeh algorithm, specify argument `HASIB.algo = “Sadeh1994"` and argument `Sadeh_axis = "Y"` to indicate that the algorithm should use the Y-axis of the sensor.
- **Galland2012:** The count-scaled algorithm proposed by Galland et al. [link](http://dx.doi.org/10.1016/j.sleep.2012.01.018). To use our implementation of the Galland2012 algorithm specify argument `HASIB.algo = “Galland2012"`. Further, set `Sadeh_axis = "Y"` to specify that the algorithm should use the Y-axis.

**Notes on the replication of the movement counts needed for the Sadeh and Galland algorithms:**

The implementation of the zero-crossing count in GGIR is not an exact copy of the original approach as used in the AMA-32 Motionlogger Actigraph by Ambulatory-monitoring Inc. ("AMI") that was used in the studies by Sadeh and colleagues in late 1980s and 1990s. No complete publicly accessible description of that approach exists. From personal correspondence with AMI, I learnt that the technique has been kept proprietary and has never been shared with or sold to other actigraphy manufacturers (time of correspondence October 2021). Therefore, if you would like to replicate the exact zero-crossing counts calculation used by Sadeh and colleague's consider using AMI's actigraph device. However, if you prioritise openness over methodological consistency with the original studies by Sadeh and colleagues then you may want to consider any of the open source techniques in this library.

Missing information about the calculation includes: (1) Sadeh specified that calculations were done on data from the Y-axis but the direction of the Y-axis was not clarified. Therefore, it is unclear whether the Y-axis at the time corresponds to the Y-axis of modern sensors, (2) Properties of the frequency filter are missing like the filter order and more generally it is unclear how to simulate original acceleration sensor behaviour with modern sensor data, and (3) Sensitivity of the sensor, we are now guessing that the Motionlogger had a sensitivity of 0.01 _g_ but without direct proof.

The method proposed by Galland and colleagues in 2012 was designed for counts captured with the Actical device (Mini Mitter Co, Inc Bend OR). Based on the correspondence with AMI we can conclude that Actical counts are not identical to AMI's actigraph counts. Further, a publicly accessible complete description of the Actical calculation does not exist. Therefore, we can also conclude that methodological consistency cannot be guaranteed for  Actical counts.


### Guiders

SIBs (explained above) can occur anytime in the day. In order to differentiate SIBs that correspond to daytime rest/naps from SIBs that correspond to the main Sleep Period Time window (abbreviated as SPT), a guiding method referred as **guider** is used. All SIBs that overlap with the window defined by guider are considered as sleep within the SPT window. The start of the first SIB identified as sleep period and the end of the last SIB identified as sleep period define the beginning and the end of the SPT window. In this way the classification relies on the accelerometer for detecting the timing of sleep onset and waking up time, but the guider tells it in what part of the day it should look, as SPT window will be defined only if SIB is detected during the guider specified window.

If a guider reflects the Time in Bed the interpretation of the Sleep Period Time, Sleep onset time and Wakeup time remains unchanged. However, we can then also assess sleep latency and and sleep efficiency, which will be included in the report.

The guiding method as introduced above can be one of the following methods:

-	**Guider = sleep log**: As presented in [before mentioned 2015 article](https://doi.org/10.1371/journal.pone.0142533).See section [on sleep analysis related arguments][Key arguments related to sleep analysis] for a discussion fo sleep log data formats. Specify argument `sleepwindowType` to clarify whether the sleeplog capture "TimeInBed" or "SPT". If it is set to "TimeInBed", GGIR will automatically expand its part 4 analyses with sleep latency and sleep efficiency assessment.
-	**Guider = HDCZA**: As presented in [our 2018 article](https://dx.doi.org/10.1038/s41598-018-31266-z). The HDCZA algorithm does not require access to a sleep log, and is designed for studies where no sleep log is available. The time segment over which the HDCZA is derived are by default from noon to noon. However, if the HDCZA ends between 11am and noon then it will be applied again but to a 6pm-6pm window.
-	**Guider = L5+/-12**: As presented in [our 2018 article](https://dx.doi.org/10.1038/s41598-018-31266-z). Twelve hour window centred around the least active 5 hours of the day.
- **Guider = setwindow**: Window times are specified by user, constant at specific clock times with argument `def.noc.sleep`.
- **Guider = HorAngle**: Only used if argument `sensor.location="hip"`, because this will trigger the identification of the longitudinal axis based on 24-hour lagged correlation. You can also force GGIR to use a specific axis as longitudinal axis with argument `longitudinal_axis`. Next, it identifies when the horizontal axis is between -45 and 45 degrees and considers this a horizontal posture. Next, this is used to identify the largest time in bed period, by only considering horizontal time segments of at least 30 minutes, and then looking for longest horizontal period in the day where gaps of less than 60 minutes are ignored. Therefore, it is partially similar to the HDCZA algorithm. When "HorAngle" is used, `sleepwindowType` is automatically set to "TimeInBed".

For all guiders other than "HorAngle" and "sleep log" argument `sleepwindowType` is automatically switched to "SPT", such that no attempt is made to estimate sleep latency or sleep efficiency.

GGIR uses by default the sleep log, if the sleep log is not available it falls back on the HDCZA algorithm (or HorAngle if `sensor.location="hip"`). If HDCZA is not successful GGIR will falls back on the L5+/-12 definition, and if this is not available it will use the setwindow. The user can specify the priority with argument `def.noc.sleep`. So, when we refer to guider then we refer to one of these five methods.

### Daysleepers (nights workers)

If the guider indicates that the person woke up after noon, the sleep analysis in part 4 is performed again on a window from 6pm-6pm. In this way our method is sensitive to people who have their main sleep period starting before noon and ending after noon, referred as daysleeper=1 in daysummary.csv file, which you can interpret as night workers. Note that the L5+/-12 algorithm is not configured to identify daysleepers, it will only consider the noon-noon time window.

### Cleaningcode {#Cleaningcode}

To monitor possible problems with the sleep assessment, the variable **cleaningcode** is recorded for each night. Cleaningcode per night (noon-noon or 6pm-6pm as described above) can have one of the following values:

-	0: no problem, sleep log available and SPT is identified;
-	1: sleep log not available, thus HDCZA is used and SPT is identified,
-	2: not enough valid accelerometer data based on the non-wear and clipping detection from part summarised over the present night where the argument `includenightcrit` indicates the minimum number of hours of valid data needed within those 24 hours.
-	3: no accelerometer data available,
-	4: there were no nights to be analysed for this person,
-	5: SPT estimated based on guider only, because either no SIB was found during the entire guider window, which complicates defining the start and end of the SPT, or the user specified the ID number of the recording and the night number in the [data_cleaning_file](#Data_cleaning_file) to tell GGIR to rely on the guider and not rely on the accelerometer data for this particular night
-	6: no sleep log available and HDCZA also failed for this specific night then use average of HDCZA estimates from other nights in the recording as guider for this night. If HDCZA estimates are not available during the entire recording then use L5+/-12 estimate for this night. The last scenario seems highly unlikely in a recording where the accelerometer was worn for at least one day.

### Difference between cleaned and full output

All the information for each night is stored in the `results/QC` folder allowing tracing of the data analysis and night selection. The cleaned results stored in the results folder. In part 4 a night is excluded from the ‘cleaned’ results based on the following criteria:

-	If the study proposed a sleep log to the individuals, then nights are excluded for which the sleep log was not used as a guider (i.o.w. nights with cleaningcode not equal to 0 or variable sleep log used equals FALSE).
-	If the study did not propose a sleep log to the individuals, then all nights are removed with cleaningcode higher than 1.

Be aware that if using the full output and working with wrist accelerometer data, then missing entries in a sleep log that asks for Time in Bed will be replaced by HDCZA estimates of SPT. Therefore, extra caution should be taken when working with the full output.

Notice that part 4 is focused on sleep research, by which the cleaned reported is the way it is. In the next section we will discuss the analysis done by part 5. There, the choice of guider may be considered less important, by which we use different criteria for including nights. So, you may see that a night that is excluded from the cleaned results in part 4 still appears in the cleaned results for part 5.

### Data cleaning file {#Data_cleaning_file}

The package allows some adjustments to be made after data quality check. The `data_cleaning_file` argument allows you to specify individuals and nights for whom part4 should entirely rely on the guider (for example if we decide to use sleep log only information).
The first column of this file should have header `ID` and there should be a column `relyonguider_part4` to specify the night.
The `data_cleaning_file` allows you to tell GGIR which person(s) and night(s) should be omitted in part 4. The the night numbers to be excluded should be listed in a column `night_part4` as header.

## Waking-waking or 24 hour time-use analysis {#Waking-waking_or_24_hour}

In part 5 the sleep estimates from part 4 are used to describe 24-hour time use. Part 5 allows you to do this in two ways: Literally 24 hours which start and end a calendar day (default midnight, but modifiable with argument `dayborder`) or from waking up to waking up. In GGIR we refer to the former as **MM** windows and to the latter as **WW** windows. The onset and waking times are guided by the estimates from part 4, but if they are missing part 5 will attempt to retrieve the estimate from the guider method, because even if the accelerometer was not worn during the night, or a sleep log is missing in a study where sleep log was proposed to the participants, estimates from a sleep log or HDCZA can still be considered a reasonable estimate of the SPT window in the context of 24-hour time use analysis.

If WW is used in combination with ignoring the first and last midnight, argument `excludefirstlast`, then the first wake-up time (on the second recording day) needs to be extracted for the first WW day. This is done with the guider method. This also means that the last WW window ends on the before last morning of the recording.

### Time series output files

If you want to inspect the time series corresponding to these windows then see argument `save_ms5rawlevels`, which allows you to export the time series including behavioral classes and non-wear information to csv files. The behavioral classes are included as numbers, the legend for these classes is stored as a separate legend file in the meta/ms5.outraw folder named "behavioralcodes2020-04-26.csv" where the date will correspond to the date of analysis.
A distinction is made between the full results stored in the `results/QC` folder and the cleaned results stored in the results folder.

### Day inclusion criteria

The full part 5 output is stored in the `results/QC` folder. The default inclusion criteria for days in the cleaned output from part 5 (stored in the `results` folder) are:

-	For both MM and WW defined days: The valid (sensor worn) time fraction of the day needs to be above the fraction specified with argument `includedaycrit.part5` (default 2/3).
-	For MM defined days only: The length of the day needs to be at least the number of hours as specified by `minimum_MM_length.part5` (default 23). Note that if your experiment started and ended in the middle of the day then this default setting will exclude those incomplete first and last days. If you think including these days is still meaningful for your work then adjust the argument `minimum_MM_length.part5`.

**Important notes:**

- No criteria is set for the amount of valid data during the SPT window, because all we are interested in part 5 is knowing the borders of the night and we trust that this was sufficiently estimated by part 4. If you disagree then please notice that all the days are included in the full report available in `results/QC` folder.
- This means that argument `includenightcrit` as used for part 4 is not used in part 5.

The `data_cleaning_file` argument discussed in [Data_cleaning_file](#Data_cleaning_file) also allows you to tell GGIR which person(s) and day(s) should be omitted in part 5. The the day numbers to be excluded should be listed in a column `day_part5` as header.

### Fragmentation metrics

When setting input argument as `frag.metrics="all"` GGIR part 5 will perform daytime behavioural fragmentation analysis. Do this in combination with argument `part5_agg2_60seconds=TRUE` as that will aggregate the time series to 1 minute resolution as is common in behavioural fragmentation literature.

In GGIR, a fragment is a defined as a sequence of epochs that belong to one of the four categories:

1. Inactivity
2. Light Physical Activity (LIPA)
3. Moderate or Vigorous Physical Acitivty (MVPA)
4. Physical activity (can be either LIPA or MVPA)

Each of these categories represents the combination of bouted and unbouted time in the respective categories. Inactivity and physical activity add up to a full day, as well as inactivity, LIPA and MVPA.
The fragmentation metrics are applied in function `g.fragmentation`.

**Literature about these metrics:**

- Coefficient of Variance (`CoV`) is calculated according to [Blikman et al. 2014](https://www.archives-pmr.org/article/S0003-9993(14)01063-6/fulltext).
- Transition probability (`TP`) from Inactivity (IN) to Physical activity (IN2PA) and from Physical activity to inactivity (PA2IN) are calculated as 1 divided by the mean fragment duration. The transition probability from Inactivity to LIPA and MVPA are calculated as: (Total duration in IN followed by LIPA or MVPA, respectively, divided by total duration in IN) divided by the average duration in IN.
- Gini index is calculated with function `Gini` from the `ineq` R package, and with it's argument `corr` set to TRUE.
- Power law exponent metrics: Alpha, x0.5, and W0.5 are calculated according to [Chastin et al. 2010](https://www.sciencedirect.com/science/article/abs/pii/S096663620900602X?via%3Dihub).
- Number of fragment per minutes (`NFragPM`) is calculated identical to metric `fragmentation` in [Chastin et al. 2012](https://academic.oup.com/ageing/article/41/1/111/46538), but it is renamed here to be a more specific reflection of the calculation. The term
`fragmentation` appears too generic given that all fragmentation metrics inform us about fragmentation. Please not that this is effectively the same metric as the transition probability, because total number divided by total sum in duration equals 1 divided by average duration. This is just different terminology for the same construct.

**Conditions for calculation and value when condition is not met:**

- Metrics `Gini` and `CoV` are only calculated if there are at least 10 fragments (e.g. 5 inactive and 5 active). If this condition is not met the metric value will be set to missing.
- Metrics related to power law exponent alpha are also only calculated when there are at least 10 fragments, but with the additional condition that the standard deviation in fragment duration is not zero. If these conditions are not met the metric value will be set to missing.
- Other metrics related to binary fragmentation (`mean_dur_PA` and `mean_dur_IN`), are calculated when there are at least 2 fragments (1 inactive, 1 active). If this condition is not met the value will is set to zero.
- Metrics related to `TP` are calculated if: There is at least 1 inactivity fragment AND (1 LIPA OR 1 MVPA fragment). If this condition is not met the `TP` metric value is set to zero.

To keep an overview of which recording days met the criteria for non-zero standard deviation and at least ten fragments, GGIR part5 stores variable `Nvaliddays_AL10F` at person level (=Number of valid days with at least 10 fragments), and `SD_dur` (=standard deviation of fragment durations) at day level as well as aggregated per person.

**Difference between fragments and blocks:**

Elsewhere in the part5 we use the term `block`. A `block` is a sequence of epochs that belong to the same behavioural class. This may sound similar to the definition of a fragment, but for blocks we distinguish every behavioural class, which includes the subcategories such as bouted and unbouted behaviour. This means that variables `Nblock_day_total_IN` and `Nblock_day_total_LIG` are identical to `Nfrag_IN_day` and `Nfrag_LIPA_day`, respectively.


**Differences with R package ActFrag:**

The fragmentation functionality is loosely inspired on the great work done by
Dr. Junrui Di and colleages in R package ActFrag, as described in
[Junrui Di et al. 2017](https://www.biorxiv.org/content/10.1101/182337v1).

However, we made a couple of a different decisions that may affect comparability:

- GGIR derives fragmentation metrics for day. This allows us to avoid the issue of
dealing with appending days. Further, it allows us to test for behavioural differences
between days of the week. It is well known that human behaviour can be driven by
weekly rhythms, e.g. work days versus weekend. Estimating fragmentation per day of
the week allows us to study and account for such possible variation. As with all other
GGIR variables we also report recording level aggregates of the daily estimates.
- Power law alpha exponent metrics were calculated according to [Chastin et al. 2010](https://www.sciencedirect.com/science/article/abs/pii/S096663620900602X?via%3Dihub) using
the theoretical minimum fragment duration instead of the observed minimum fragment
duration.

## Why use data metric ENMO as default?

GGIR offers a range of acceleration metrics to choose from, but only one metric can be the default. Acceleration metric ENMO (Euclidean Norm Minus One with negative values rounded to zero) has been the default metric in GGIR. In 2013 we wrote a paper in which we investigated different ways of summarising the raw acceleration data. In short, different metrics exist and there is very little literature to support the superiority of any metric at the time. As long as different studies use different metrics their findings will not be comparable. Therefore, the choice for metric ENMO is partially pragmatic. GGIR uses ENMO as default because:

1.	ENMO has demonstrated value in describing variance in energy expenditure, correlated with questionnaire data, able to describe patterns in physical activity
2.	ENMO is easy to describe mathematically and by that improves reproducibility across studies and software tools
3.	ENMO attempts to quantify the actual biomechanical acceleration in universal units.
4.  The 2013 paper showed that when ENMO is used in combination with auto-calibration it has similar validity to filter-based metrics like HFEN and BFEN, which are conceptually similar to metrics proposed later such as MIMSunit, MAD, AI0.
5.  Studies who have criticised ENMO consistently failed to apply auto-calibration, or attempted to apply auto-calibration in a lab setting, ignoring the fact that the auto-calibration is not designed for short lasting lab settings. It needs free-living data to work properly. Further, studies are often not clear about how the problematic zero imputation during the idle sleep mode in ActiGraph devices is dealt with. See also paragraph: [Published cut-points and how to use them]([#Published cut-points and how to use them).

See also [this](https://medium.com/@vincentvanhees/ten-misunderstandings-surrounding-information-extraction-from-wearable-accelerometer-data-a4f767a865b6) blog post on this topic.

## What does GGIR stand for?

I wanted a short name and not to spend too much time finding it. At the time I was primarily working with GENEActiv and GENEA data In R, and that's how the name GGIR was born: Short, easy to remember, and as acronym sufficiently vague to not be tight up with a specific functionality. However, later the functionality expanded to other sensor brands, so the abbreviation has lost its functional meaning.

# Other Resources
- The [GGIR package manual](https://CRAN.R-project.org/package=GGIR) provides documentation on individual functions.
- For general questions about how to use GGIR join our [google group/mailing list](https://groups.google.com/forum/#!forum/rpackageggir).
- For bug reports please post them [here](https://github.com/wadpac/GGIR/issues).

# Citing GGIR {#CitingGGIR}
A correct citation of research software is important to make your research reproducible and to acknowledge the effort that goes into the development of open-source software.

To do so, please report the GGIR version you used in the text. Additionally, please also cite:

1. Migueles JH, Rowlands AV, et al. GGIR: A Research Community–Driven Open Source R Package for Generating Physical Activity and Sleep Outcomes From Multi-Day Raw Accelerometer Data. Journal for the Measurement of Physical Behaviour. 2(3) 2019. doi: 10.1123/jmpb.2018-0063.

If your work depends on the quantification of **physical activity** then also cite:

2. van Hees VT, Gorzelniak L, et al. Separating Movement and Gravity Components in an Acceleration Signal and Implications for the Assessment of Human Daily Physical Activity. PLoS ONE 8(4) 2013. [link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061691)
3. Sabia S, van Hees VT, Shipley MJ, Trenell MI, Hagger-Johnson G, Elbaz A, Kivimaki M, Singh-Manoux A.
Association between questionnaire- and accelerometer-assessed physical activity: the role of sociodemographic factors. Am J Epidemiol. 2014 Mar 15;179(6):781-90. doi: 10.1093/aje/kwt330. Epub 2014 Feb 4.
PMID: 24500862 [link](https://pubmed.ncbi.nlm.nih.gov/24500862/)

If you used the **auto-calibration functionality** then also cite:

4. van Hees VT, Fang Z, et al. Auto-calibration of accelerometer data for free-living physical activity assessment using local gravity and temperature: an evaluation on four continents. J Appl Physiol 2014. [link](https://journals.physiology.org/doi/full/10.1152/japplphysiol.00421.2014)

If you used the **sleep detection** then also cite:

5. van Hees VT, Sabia S, et al. A novel, open access method to assess sleep duration using a wrist-worn accelerometer, PLoS ONE, 2015 [link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0142533)

If you used the **sleep detection without relying on sleep diary** then also cite:

6. van Hees VT, Sabia S, et al. Estimating sleep parameters using an accelerometer without sleep diary. Scientific Reports 2018. doi: 10.1038/s41598-018-31266-z. [link](https://www.nature.com/articles/s41598-018-31266-z)


## Copyright for GGIR logo
The copyright of the GGIR logo lies with Accelting (Almere, The Netherlands), please contact v.vanhees@acceleting.com to ask for permission to use this logo.


```{r, echo=FALSE, out.width = "60%", out.extra='style="border: 0; padding:20px"'}
knitr::include_graphics("GGIR-MASTERLOGO-RGB.png")
```
\name{NEWS}
\title{News for Package \pkg{GGIR}}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\section{Changes in version 2.5-7 (release date:25-01-2022)}{
\itemize{
  \item Fixing issue introduced in 2.5-6 release: maxNcores was not defined in part 1, params_output not passed on to part 3, and params objects were unneccesarily repeatedly checked
  }
}
\section{Changes in version 2.5-6 (GitHub-only-release date:18-01-2022)}{
\itemize{
  \item Internal re-structuring: Arguments are now internally passed via parameter objects to 
    make code and code documentation easier to navigate. This new structure preserve backward 
    compatibility with the user's older R scripts.
  }
}
\section{Changes in version 2.5-5 (GitHub-only-release date:18-01-2022)}{
\itemize{
  \item Fixed bug #484 affecting part 5 report not being able to generated because non-matching columns in milestone data
  \item Part 2: BrondCounts and ZeroCrosssing Counts no long auto-scaled by 1000 like all the other metrics.
  \item Part 3: Fixed warning when SRI is calculated on DST day with 23 hours
  \item Part 5: sibreport and timeseries now stored per sib definition.
}
}
\section{Changes in version 2.5-4 (GitHub-only-release date:11-01-2022)}{
  \itemize{
    \item Part 1: Now able to process modern .gt3x data with thanks to R package read.gt3x.
    \item Part 1: Timegaps and zeros across all three axes in ActiGraph data (.csv and .gt3x) now automatically imputed.
    \item Adding BrondCounts as optional acceleration metric and dependency on activityCounts, and enabling Sadeh and Galland algorithms to use it.
  }
}
\section{Changes in version 2.5-3 (GitHub-only-release date:23-12-2021)}{
  \itemize{
    \item Part 1: Depricating function g.metric as its functionality has been taken over by g.applymetrics.
    \item Part 2: Adding warning when ID cannot be extracted from file based on argument idloc.
    \item Part 4: Empty or incomplete sleeplog rows now better ignored
    \item Part 4: No longer an R warning is given when ID is missing in sleeplog, because this is common and user can already see it in csv reports.
  }
}
\section{Changes in version 2.5-2 (GitHub-only-release date:09-12-2021)}{
  \itemize{
    \item Part 3: Now makes sure that HASPT is skipped when user configures def.noc.sleep as a set time window. This
    feature that probaly few people use nowadays broke with the 2.5-0 release.
    \item g.shell.GGIR: Now gives warning when user supplies double arguments.
    \item Part 4: Now warns when none of the IDs in the sleeplog could be matched with accelerometer data.
    \item Part 1: read.myacc.csv fix bug with argument rmc.check4timegaps
    \item Part 3: Fix #472 SRI calculation not possible when complete absence of sleep in recording
    \item Part 5: Experimental nap detection added to report and time series, currently only model for 3.5 year olds available.
  }
}
\section{Changes in version 2.5-1 (GitHub-only-release date:03-11-2021)}{
  \itemize{
    \item g.shell.GGIR: Removing forced assignment of sleepwindowType argument to "SPT"
    \item g.shell.GGIR: Reviving setwindow option def.noc.sleep, which broke
  }
}
\section{Changes in version 2.5-0 (release date:25-10-2021)}{
\itemize{
  \item Minor updates to function documentation and vignettes
  \item visualreport: Month now displayed by name rather than number to avoid confusion about format.
}
}
\section{Changes in version 2.4-3 (GitHub-only-release date:20-10-2021)}{
  \itemize{
    \item Vignette expanded with documentation on non-default variables in part 4 csv report
    \item Vignette added with a tutorial on how to perform segmented day analysis
    \item Part 1: The sensor fusion functionality introduced in version 2.2-2, as a bit
    of an experimental development, has now been removed from GGIR as added value turned out
    to be limited.
    \item Part 2: Code for ID extraction tidied up, and new idloc argument options added.
    \item Part 2: Activity log can now also have empty days.
    \item Part 2: Now also exported in long format if qwindow has length longer than 2.
    \item Part 3 and 4 expanded with Sleep Regularity Index.
    \item Part 4 number of awakenings is now also in the output.
    \item Part 5 LUX per segment calculation bug fixes.
    \item visualreport code modified such that it visualises any day with at least 30 minutes of data
    unless sleep could not be estimated from the corresponding night.
  }
}
\section{Changes in version 2.4-2 (GitHub-only-release date:14-09-2021)}{
  \itemize{
    \item Part 1 now able to derive zero-crossing counts needed for Sadeh algorithm
    \item Part 2 function g.conv.actlog to use activity log to guide qwindow now able to tailor dateformat
    \item Part 3 now option to tailor SPT detection for hip data by using horizontal
      longitudinal angle as indicator of lying
    \item Part 3 code restructured to estimate SPT and SIB in separated function (HASIB and HASPT)
    such that g.sib.det is more readable.
    \item Part 3 now option to select SIB algorithm, including: vanHees2015, Sadeh1994,
    and Galland2012
    \item Part 4 now able to handle more advance sleeplog format that also contain nap
    and nonwear entries per date
    \item Part 4 report now includes guider estimates also in 'cleaned' report
    \item Part 4 WASO and NOA (Number of awakenings) added to output
    \item Part 4 Able to handle time in bed diaries/algorithms, and generates sleep
    latency and sleep efficiency estimate if available
    \item Part 5 now able to export sustained inactivity bouts and self-reported
    naps and nonwear to csv report to aid nap analysis
    \item Part 5 numeric timing of MX-metrics changed to be hour of the day
    \item Part 5 Lux per segment variables now turn LUX during SPT to zero
    \item Part 4 legend added to visualisation
    \item Part 5 bug fixed that caused part5 single row output to be stored as a single column
  }
}
\section{Changes in version 2.4-1 (GitHub-only-release date:12-07-2021)}{
  \itemize{
    \item Part 1 fixed bug related to g.getmeta reading movisens files
    \item Part 1 resampling function expanded with option to choose interpolation technique
    \item Part 3 fixed bug introduced in 2.4-0 making it impossible to process long
    recordings with many sustained inactivity bouts
  }
}
\section{Changes in version 2.4-0 (release date:03-06-2021)}{
  \itemize{
    \item Part 5 LUX variables per segment of the day improved
    \item Removing calls to closeAllConnections as requested by CRAN
  }
}
\section{Changes in version 2.3-3 (GitHub-only-release date:20-05-2021)}{
  \itemize{
    \item Part 5 fragmentation analysis added
    \item Part 5 report, added argument week_weekend_aggregate.part5
    \item Part 5 intensity gradient now also extracted here for waking hours
    \item Part 5 LUX variables added (beta-implementation)
    \item Part 2 and 5, Bout detection algorithm improved, see bout.metric 6
    \item Parts 1, 2, 3 and 5: Added argument maxNcores to control maximum number
    of cores to be used in parallel processing mode
    \item Part 1 fixed bug related with ENMOa calculation which made the metashort data frame to be empty
    \item Part 2 fixed bug to report generation in part 2
    \item Fixed bug related to visualization report when a file has not been processed in part 4
  }
}
\section{Changes in version 2.3-2 (GitHub-only-release date:21-04-2021)}{
  \itemize{
    \item Part 1 when using ad-hoc data format: now able to handle timestamp column.
    \item Part 2 report when using ad-hoc data format with header: now able to extract
    and store serial number in part2 report.
    \item Test file generation for unit tests improved. Thanks Lena Kushleyeva.
    \item Part 2 estimation of longitudinal axis and related error handling improved
  }
}
\section{Changes in version 2.3-1 (GitHub-only-release date:26-03-2021)}{
  \itemize{
    \item Part 3 bug fixed in usage of argument dayborder when not 0.
    Thanks Ruben Brondeel for spotting this.
    \item Part 3 bug fixed in recognition fraction night invalid for first night.
    Thanks Ruben Brondeel for spotting.
    \item Part 1 fixed warning related to closing connection to bin files. Thanks John Muschelli.
    \item Part 2 fixed warning that deprecated IVIS argument does not exist
  }
}
\section{Changes in version 2.3-0 (release date:17-02-2021)}{
  \itemize{
    \item Part 5 bug fixed for g.part5.fixmissingnight when night is also missing in sleeplog
    \item Part 1 Documentation for sensor fusion functionality improved
    \item Part 2 csv report generation speeded up
  }
}
\section{Changes in version 2.2-2 (GitHub-only-release date:13-02-2021)}{
  \itemize{
    \item Part 5 bug fix in waking up time for first night, in WW mode and HDCZA guider
    \item Various minor function expansions to documentation and vignette
    \item Part 2 Earlier draft implementation of IV and IS metrics now better tested and bugs fixed.
    \item Part 2 Default value for argument IVIS.activity.metric changed to 2.
    \item Migrate codecoverage testing to Github Actions
    \item Part 1 Implemented initial sensor fusion functionality for AX6 cwa data
  }
}
\section{Changes in version 2.2-1 (GitHub-only-release date:12-01-2021)}{
  \itemize{
    \item Part 3 Fix bug introduced in 2.2-0 for recordings without sustained inactivity bouts
    \item Transitioned from Travis+Appveyor CI to GitHub Actions
    \item Part 5 fixed issue with expected max number of files to process
    \item Enabled console display of GGIR version number when running g.shell.shell under Windows
    \item Embedded new GGIR logo and update layout of vignette
  }
}
\section{Changes in version 2.2-0 (release date:22-11-2020)}{
  \itemize{
    \item Part 5 improved midnights identification when the recording starts at 0am
    \item Part 5 fixed days misclassification for MM output when recordings starts with some non-wear days
    \item Documentation added to vignette for L5TIME_num and M5TIME_num
  }
}
\section{Changes in version 2.1-3 (GitHub-only-release date:20-10-2020)}{
  \itemize{
    \item Depricating sessionInfo storage as it caused problems for large scale parallel processing
    \item Part 1 fix to 1 minute time shift when recording starts at full minute.
    \item Improved compatibility with older .cwa file formats.
  }
}
\section{Changes in version 2.1-2 (GitHub-only-release date:25-09-2020)}{
  \itemize{
    \item Part 5 bug fixed with day name and date allocation for daysleepers for MM report
    \item Fix bug part 4 and 5 with missing nights which are not accounted for when assessing max night number.
    \item Fix bug part 3 introduced in version 2.1-0 re. SPT detection for daysleepers.
    \item Fix redefinition of the time windows (with "MM") in part 5 to include all recording days when the measurement starts with multiple non-wear days.
    \item GGIR version now displayed in console when running GGIR
    \item Bout metric 5 added to enable not allowing for gaps in bouts.
  }
}

\section{Changes in version 2.1-1 (GitHub-only-release date:07-09-2020)}{
  \itemize{
    \item Part2 argument qwindow now able to handle person specific activity diary
    \item Intensity gradient now also extracted from long MX windows, see MX.ig.min.dur
  }
}
\section{Changes in version 2.1-0 (release date:27-08-2020)}{
  \itemize{
    \item Part 2 now able to be guided by activity log (see qwindow)
    \item Part 2 now stores output also in long format if day is segmented
    \item Actigraph serial number now extracted from the csv fileheader
    \item Actigraph serial number now used to assign correct dynamic range per
    Actigraph generation (6 and 8 g respectively). Previously 8 assumed.
    \item Part 2 indicator added of which axis is most correlated with itself
    over 24 hours, for hip data this may be indicator of vertical axis orientation.
    \item Part 1 bug fixed in timestamp generation that caused measurements
    that start exactly at the hour and 0 minutes to have a 1 minute offset in time.
    \item Part 1 to reduce data storage size epoch level metric values are now
    rounded to 4 decimal places.
    \item Part 4 report, issue fixed with double night entries.
    \item Bug fixed with Actigraph starttime recognition when machine language
    is not English and months are expressed in b or hh format.
    \item Bug fixed in accounting for daysleeper in first night ofr recording.
  }
}
\section{Changes in version 2.0-2 (GitHub-only-release date:18-06-2020)}{
  \itemize{
    \item Visual report, fixed argument viewingwindow for g.plot5 (visualreport)
    \item Part 1 rmc.myacc.csv now raises warning when timestamps are not recognised
    \item Part 1 fixed rmc.dynamic_range not correctly passed on in 2.0-1
    \item Part 5 sub-function round_seconds_to_ws3new now handles missing values
  }
}
\section{Changes in version 2.0-1 (GitHub-only-release date:04-06-2020)}{
  \itemize{
    \item Now able to read gzipped csv files
    \item Part 2 Fixed midnight selection being off by one sample
    \item Part 1 tidied up metric calculation and added a few metrics
    \item Part 2 Added strategy 4 (ignore data before first midnight)
    \item Data from movisens accelerometers can be now processe
    \item Part 5 timestamps in timeseries output in RData format (new since 2.0)
    now correctly accounts for timezone
    \item Part 4 now has option to only exclude the first or only the last night
  }
}
\section{Changes in version 2.0-0 (release date:30-04-2020)}{
  \itemize{
    \item Now R version 4.0 compatible
    \item Part 1 Clipping detection expanded: If any value in block more than 150 percent
    dynamic range then ignore entire block.
    \item Part 2 report now able to handle changing variable count due to missing data
    \item Part 2 L5M5 better able to handle small qwindow intervals
    \item Part 3 HDCZA algorithm expanded to be able to detect daysleepers
    \item Part 3 various improvements to qc plots.
    \item Part 5 now also stores full and cleaned output
    \item Part 5 now better handles missing days in part 4 output.
    \item Part 5 behavioral class SIB removed from daytime
    \item Part 5 time series export more user-friendly.
    \item Part 5 function code split up in 7 new functions.
    \item Part 4 + 5 argument data_cleaning_file added.
    \item Part 4 + 5 output variable names improved and documented in vignette
    \item Numunpack function moved back to c++
    \item Various updates to visualreport (plot5 function)
    \item External function embedding feature added
    \item We now consistently refer to ID (not id) and calendar_date, spelling was
    inconsistent.
    \item Vignette now has documentation on sleep and time-use analysis.
  }
}
\section{Changes in version 1.11-1 (GitHub-only-release date:28-02-2020)}{
  \itemize{
    \item Metric ENMOa now facilitated for MVPA calculation
    \item Bug fixed with part4 daysleeper handling
    \item Part5 WW window calculation improved, first day now uses sleeplog or
    HDCZA algorithm estimate, and last day is ignored if no sleep estimates are
    available. This also affects csv exports by argument save_ms5rawlevels.
    \item Added explanation to vignette on how to use published cut-points.
    \item Axivity AX6 (acc + gyro) in cwa format now supported for file reading,
    actually using the gyro data for feature calculation is future work.
    In the mean time, gyro signal will be ignored by the rest of GGIR.
    \item Axivity AX3 acc data in cwa format can now also be read if dynamic
    range is not 8g. Previously this was not possible.
    \item Fixed bug related to visual report generation when qwindow is set to
    non-default value.
    \item Added way to handle Actigraph files which start with several days of
    zeros which complicates the auto-calibration.
    \item Default for desiredtz to timezone of machine.
  }
}
\section{Changes in version 1.11-0 (release date:03-12-2019)}{
\itemize{
  \item Fixed bug that emerged in previous version with GENEActiv .bin data not
  being processed by g.getmeta for some files.
  \item read.myacc.csv is now able to resample data with irregular sample rates
  and handle timestamps in character format.
  \item function resample can now handle any matrix size, previous only 3 columns.
  \item Fixed bug when using multiple non-angle metrics in part1 and trying to calculate
  1to6am average metric value in part 2.
  \item Expanded Actigraph date format recognition ability.
  \item visualisationreport (function g.plot5) enhanced with colour coding for activity
  classes.
  \item Fixed bug in sleep period time recognition for first day of measurement.
}
}
\section{Changes in version 1.10-7 (release date:06-10-2019)}{
\itemize{
  \item Fixed functionality to supply calibration coefficients file to backup.cal.coef.
  \item Fix OSx flavor not being released on CRAN in previous version.
  \item Upgrades to foreach loop to ease package maintenance
  \item Documentation part4 expanded to clarify difference between full and cleaned report.
  \item Non-wear detection now possible at 1 minute resolution, previously 5 minute.
  \item Function read.myacc.csv now able to utilize 3rd party wear detection.
}
}
\section{Changes in version 1.10-5 (release date:13-09-2019)}{
\itemize{
  \item Device serial number recognition in Axivity cwa files fixed
  \item New GitHub release, because previous version did not install.
}
}

\section{Changes in version 1.10-4 (release date:08-09-2019)}{
\itemize{
  \item Fixed bug introduced in version 1.10-1 in the conversion from numeric to character sleep times
  \item Dependencies of dependencies removed from the DESCRIPTION file
  \item Fixed 1to6am variables, which was wrong in version 1.10-1
  \item Added functionality to handle accelerometer data from any accelerometer brand
  stored in csv files via read.myacc.csv. Pass on the arguments of this function to
  g.shell.GGIR to use this functionality.
  \item Added reference to new GGIR paper to the documentation
}
}
\section{Changes in version 1.10-1 (release date:23-08-2019)}{
\itemize{
  \item Configuration file option now added to g.shell.GGIR and documented in vignette
  \item metric lfen added (low-pass filter signels followed by Euclidean norm)
  \item issues fixed with passing on of hb and lb arguments
  \item argument backup.cal.coef can now also handle data_quality_reports.csv files
  \item part 1 now automatically uses previously generated calibration coefficients if the datafile
  was previously processed, see documentation g.part1 for further details.
  \item Enabled multiple values for argument winhr, by which part2 can now calculated
  for example L3M3, L5M5, L6M6, L10M10 all in one go. Further, option added (qM5L5) to
  extract percentiles (quantiles) from the value distribution corresponding to these windows.
  \item Moved IVIS calculation to seperate function, and split up function g.analyse.
  \item Now possible to specify time windows that need to be ignored for imputation, see TimeSegments2ZeroFile.
  \item Default value for argument mode changed from mode = c(1,2) to mode = 1:5, to perform all the parts.
  \item Checks added for user write and read access permission, and subsequent warnings given..
  \item Parts 1, 2, 3, and 5 can now use multi CPU cores which speeds up the processing.
  \item Argument minimumFileSizeMB added to g.part1 to aid filtering out too small data files.
  }
  }
\section{Changes in version 1.9-2 (release date:03-07-2019)}{
\itemize{
  \item Added functionality to work with studies where accelerometer is configured in one
  timezone and used in other timezone. Only functional for AX3 cwa data at the moment.
  See argument 'configtz'.
  \item Sleep estimation is now skipped if a day only has one sustained inactivity bout
  \item Arguments ignorenonwear default value changed to TRUE and def.noc.sleep default changed
  to 1 in line with literature.
  \item Fixed AX3 csv format starttimestamp recognition
  \item Part5 csv export now also includes class labels (previously only class numbers).
  }
  }
\section{Changes in version 1.9-1 (release date:08-05-2019)}{
\itemize{
    \item Fixed part5 output midnight-midnight window when monitor not worn during first days.
    \item Fixed assumption that when using argument idloc=2 the ID has a letter at the end, and automatically
    removes the last value in the index. The code now first checks for this assumption.
    \item Update vignette with a more elaborate explanation of the optional arguments to g.shell.GGIR.
  }
}
\section{Changes in version 1.9-0 (release date:14-03-2019)}{
\itemize{
  \item functionality storefilestructure should now store filestructure in output of part 2, 4 and 5.
  \item filelocationkey.csv that was previously written by storefilestructure was redundant and removed.
  \item sessioninfo storage improved.
  \item Fixed bug that caused part2 to provide incorrect window specific estimates on first day.
  of measurement if day is incomplete (not 23, 24 or 25 hours).
  \item Calibrate function now better able to handle files with more than a week of data, where
  auto-calibration struggles to find enough sphere data in the first week.
  \item Fixed part5 output midnight-midnight window when monitor not worn during first days.
  }
}
\section{Changes in version 1.8-1 (release date:11-01-2019)}{
\itemize{
  \item Part4 handling of clocktime 9am corrected in addition to fix from version 1.7-1.
  \item Fixed bug for Actigraph csv header recognition when column 2 and 3 are NA (prevented processing before)
  \item Fixed bug in g.report.part5 in the calculation of the total number of valid days per person.
  \item Fixed bug that caused part5 to struggle with timezones west of greenwhich time.
  \item desiredtz added as explicit argument to g.inspectfile, g.cwaread, and g.dotorcomma.
  \item Fixed bug in pageindexing in g.readaccfile when machine runs out of memory.
  \item Fixed bug in pageindexing in Actigraph csv files (10 rows in raw data skipped every block (day)) of data.
  \item Added more informative warning message in g.report.part2 if file cannot be read.
  \item Fixed bug in scenario when person is daysleeper and wakingup time occurs before noon, and 12 hours too early.
  \item Fixed bug re. storefolderstructure=TRUE causing 2 variables to drop in g.report.part4 if storefolderstructure=TRUE.
  \item g.intensitygradient enabled to handle absense of data.
  \item tidied up some of the redundant or even confusiong information printer to the console
  }
}
\section{Changes in version 1.7-1 (release date:25-11-2018)}{
\itemize{
  \item Fixed order of Nbouts output in g.part5 was not consistent with bout duration.
  \item Functionality added to read AX3 Axivity csv files that have the following
  characteristics: Raw data in g-units, not resampled, and with timestamps stored in
  the first column.
  \item Fixed bug that caused the date of the last measurement day in part5 output
  to be one day ahead if argument dayborder=12.
  \item Part5 struggled to process measurements with more than 40 days, now fixed.
  \item g.getstarttime can now also handle dates in Actigraph csv-file headers that are
  dot separated, e.g. 20.05.2016, before it only handled 20-05-2016 and 20/05/2016.
  \item Part4 handling of clocktime 9am corrected.
  \item Filename identification in part5 improved when saving raw level data.
  \item Fixed ability to read wav file header when shorter than expected.
  \item Fixed issue with CWA read functionality causing some files not to be completely read
  }
}
\section{Changes in version 1.6-7 (release date:23-9-2018)}{
\itemize{
  \item Link with Zenodo for doi generation removed.
  \item Broken url fixed in vignette.
  }
}
\section{Changes in version 1.6-1 (release date:21-9-2018)}{
\itemize{
  \item Report part 4: Count of available nights with accelerometer data fixed
  \item Report part 4: NA and NaN values replaced by empty cells like in other reports
  \item Intensity gradient analysis added to part2 output according to Rowlands et al.
  MSSE 2018, doi: 10.1249/MSS.0000000000001561
  \item Documentation on part4 output variables improved.
  \item Providing incorrect value of sleeplogidnum in part4 should now provide a more
  informative error message
  \item Added functionality to handle timestamps that start with the year.
  }
}
\section{Changes in version 1.6-0 (release date:29-7-2018)}{
\itemize{
  \item Fixed timezone dependency of g.analyse (affected only order of columns), g.part5
  (affected time detection), and consequently test_chainof5parts.R
  \item Read functionality for Actigraph csv files speeded up by replacing
  read.csv by data.table::fread
  \item Argument qwindow is now able to handle input vectors longer than 2
  and will derive all part2 variables for each time window that can be defined
  from the values of qwindow, e.g. value =c(0,8,24) will analyse the windows:
  0-8, 8-24 and 0-24 hour.
  \item Argument L5M5window depricated because argument qwindow now defines the
  window over which L5M5 analysis is performed.
  \item Argument winhr is now reflected in g.analyse/g.part2 output variable
  names, previously this variable name was hardcoded as L5M5, even if winhr was not 5.
  \item Part2 output variable names updated to be more consistent in
  structure and more explicit about the timewindow over which they
  are calculated. The variables that were calculated over the full
  recording (using diurnal normalization) now have the extension
  "fullRecording", this in contrast to the variables that are only
  calculated from measurement days with 'enough' valid data.
  \item Fixed calculation of N valid WEdays and Nvalid WKdays in part2 that was wrong
  since version 1.5-21. It counted all days and did not exclude days with insufficient
  amount of data.
  \item Fixed warning message in test of g.part5.
  \item output variable acc_timeinbed renamed to acc_SptDuration to avoid confusion
  with terminology used in supporting papers, which are about SPT (sleep period time)
  detection and not about time in bed detection.
  \item output variable acc_dur_noc renamed to acc_SleepDurationInSpt to improve clarity
  of variable name relative to acc_SptDuration
  }
}
\section{Changes in version 1.5-24 (release date:9-7-2018)}{
\itemize{
  \item Variable names ENMO accidentatly disappeared from g.analyse output
  in 1.5-23, this has now been reversed.
  \item Fixed issue in g.part5 in handling the last day of measurement when using
  'MM' windows and dayborder not equal to 0 (midnight) sometimes resulting in
  the last day longer than 25 hours.
  \item Fixed error message in g.part5 that occurs when calculating L5M5 when day
  is shorter than L5M5 time window.
  }
}
\section{Changes in version 1.5-23 (release date:4-7-2018)}{
\itemize{
  \item Unit tests speeded up by using a smaller test file.
  }
}
\section{Changes in version 1.5-22 (release date:3-7-2018)}{
\itemize{
  \item g.part4 bug fixed that was introduced in version 1.5-21 regarding the handling
  of daysleepers (people who wake up after noon) causing sleep estimates to be zero.
  }
}
\section{Changes in version 1.5-21 (release date:22-4-2018)}{
\itemize{
  \item g.part5 is now able to generate summary from all measurment days (thanks Jairo).
  \item MM results in g.part5 now correctly stores onset and waking for single date per row.
  \item Various checks added to g.part4 to ensure measurements in all shapes and size can be
  processed.
  \item minor improvements to documentation g.part5
  \item size of example data and vignette images reduced to reduce package size
  }
}
\section{Changes in version 1.5-18 (release date:18-3-2018)}{
\itemize{
  \item Improved handling of day saving time in g.part3 and g.part4 within the two nights
  in a year when the clock changes (DST was already correctly handled in the rest
  of the year)
  \item Algorithm HDCZA is now the default algorithm to use, if a sleeplog
  is used and an entry is missing for a particular night.
  \item pdf generation in part 3 is now optional (argument do.part3.pdf), this may
  be useful for slightly speeding up data processing as it takes a few second to generate
  the plot.
  \item test added for g.getbout function.
  \item bug fixed in timestamp recognition for object timebb in function g.part5
  }
}

\section{Changes in version 1.5-17 (release date:19-2-2018)}{
\itemize{
  \item SPT-window detection now updated with a constrained threshold to make it more robust
  against between accelerometer brand differences. This is the approach used for our PSG in
  <https://www.biorxiv.org/content/early/2018/02/01/257972>
  }
}

\section{Changes in version 1.5-16 (release date:17-1-2018)}{
\itemize{
  \item cwaread issue #57 fixed
  \item SPT-window detection included (work in progress)
  \item cpp code fixed which did not compile anymore
  \item a machine specific function test removed following feedback from CRAN maintainers
  }
}

\section{Changes in version 1.5-12 (release date:8-8-2017)}{
\itemize{
  \item Fixed bug introduced in 1.5-1: large window size of 3600 seconds was accidentally
  hardcoded when the g.readaccfile function was added to GGIR in version 1.5-1
  \item Codecov testing added and badge added to the README file
  \item Functions added to create dummy accelerometer file (csv) and dummy sleeplog (csv),
	needed for testing
	\item Bug fixed in g.wavread to recognize .wav extension file header for files with alternative
	header size
  \item Default IVIS_epochsize_seconds parameter updated from 30 to 3600
  \item g.part1 messages on the consolo are now condensed printing only the number of the
	blocks loaded separated by spaces rather than new lines
  \item Split g.part5 function into multiple smaller functions
  \item Replace hard-coded "Europe/London" in g.part5 by desiredtz, to make the function
	work for users outside the UK
  \item Data frame output from g.part5 is now tidied up by removing empty rows and columns generated
  \item Calculation of mean amplitude deviation (MAD) is now implemented in g.part1 by the argument do.mad
  \item Percentiles and levels in g.part2 are now calculated for all the acceleration metrics selected
  \item g.part3, g.part4 and g.part5 are now independent on metric ENMO to work, argument acc.metric allows
	the user to select which metric to use to describe behavior
  \item Argument dayborder is now included in g.part5 to consider the whole measurement in case the protocol
	starts after midnight
	\item Jairo Migueles added to list of contributors
  }
}

\section{Changes in version 1.5-10 (release date:12-7-2017)}{
\itemize{
  \item Date format recognition improved for Actigraph csv files
  }
}

\section{Changes in version 1.5-9 (release date:21-5-2017)}{
\itemize{
  \item g.inspectfile now also functional with cwa data
  \item option added to enforce dynamic range with argument dynrange
  \item vignette expanded
}
}
\section{Changes in version 1.5-7 (release date:9-5-2017)}{
\itemize{
  \item vignette expanded
  \item Bugs fixes in relation to new cwa-read functionality c++ routine registration
  \item Documentation added for all underlying functions
}
}
\section{Changes in version 1.5-3 (release date:29-4-2017)}{
\itemize{
  \item vignette added
  \item Bugs fixes in relation to new cwa-read functionality
  \item Bugs fixes in correct number of days recognition in part5
}
}
\section{Changes in version 1.5-1 (release date:23-4-2017)}{
\itemize{
  \item Removed teLindert2013 metric, because it was not used and not verified
  \item Split g.getmeta function into multiple smaller functions
  \item Added IS and IV variables to g.analyse (still in explorative version)
  \item bug fixed related to wav file read errors
  \item function g.cwaread added (credits to E. Mirkes) for reading
  Axivity .cwa-format data. g.shellGGIR will automatically use this function
  when input data has .cwa extension
  \item bug fixed related to for GENEactiv starttime recognition which was
  introduced in version 1.2-11
  \item g.part5 documentation added on its output
}
}
\section{Changes in version 1.4 (release date:22-1-2017)}{
\itemize{
  \item bug fixed in functionlity to process only specific days in measurement
  (credits to J Heywood)
  \item bug fixed in midnight recognition in g.part5
  \item improvement to handling of measurements that start a few minutes before
  midnight (credits to E Mirkes)
  \item bug fixed related processing files shorter than 1 day, introduced in
  previous version
  \item documentationa added for Axivity wav-format data
  \item start made with implementing code testing functionality using testthat and covr
  \item documentation improved for argument def.noc.sleep in function g.part4
}
}
\section{Changes in version 1.3-2 (release date:28-11-2016)}{
\itemize{
  \item g.part5 added. g.part5 merges the output from the first four parts
  \item Functionality added to read Axivity wav-format files with acceleration in first
  three channels. No documentation added yet until I have more confirmation that it works
  well
}
}
\section{Changes in version 1.2-11 (release date:31-8-2016)}{
\itemize{
  \item Bug fixed related to MVPA variable. The bug was a result of
  the updates in version 1.2-10
}
}
\section{Changes in version 1.2-10 (release date: 28-8-2016)}{
\itemize{
  \item Changed function argument 'mvpa.2014' into 'bout.metric' across
  the package in preparation for a central defintion of bouts for future
  GGIR version which will not only provide bout calculations for MVPA
  but also for inactivity. Further, I added function g.getbout to improve
  transparency about the bout calculations
  \item Updated documentation for function g.analyse to clarify different
  bout metric definitions
  \item Improvements to the functionality to only process specific days from
  a long accelerometer file using the argument selectdaysfile
  \item Timestamps are now in ISO 8601 format. I have updated the code such that it
  can still handle old timestamp format, but newly processed files will produce
  timestamps in the ISO 8601 format. The practical difference is that it will
  include a numeric timezone indicator.
  \item Bugs fixed in data selection in g.getmeta function. In the old code it tended
  to drop the last 30-45 minutes of a file.
  \item Added more optional features to be generated by g.getmeta, including
  rolling medians of the acceleration signals.
}
}
\section{Changes in version 1.2-8 (release date: 24-5-2016)}{
\itemize{
  \item Updated documentation for function g.analyse to clarify different
  mvpa bout definitions
  \item mvpa.2014 = TRUE turned back on again (was disable in last version)
}
}
\section{Changes in version 1.2-7 (release date: 12-5-2016)}{
\itemize{
  \item mvpa.2014 TRUE/FALSE was swapped, FALSE is now the default
  \item mvpa.2014 = TRUE disabled
}
}
\section{Changes in version 1.2-6 (release date: 10-5-2016)}{
\itemize{
  \item Modified warning message in relation to the change in MVPA bout defintion
}
}
\section{Changes in version 1.2-5 (release date: 8-5-2016)}{
\itemize{
  \item Accelerometer non-wear time is now also reported in the output of part 4
  \item Part 1 is now able to only process specific days of a measurement via argument
  selectdaysfile. This is useful when measurement lasts for a week and the participant
  is instructed to only wear the accelerometer on one or two specific days.
  \item Argument mvpa.2014 and closedbout added to function g.analyse. The calculation of
  MVPA (moderate and vigorous physical activity) has been available since 2014.
  This calculation has been improved, but the user has the option to continue using
  the old calculation.
}
}
\section{Changes in version 1.2-2 (release date: 7-1-2016)}{
\itemize{
  \item Bug fixed in the loading of data files with (very) large amounts of data
  \item Bug fixed in starttime allocation for measurements starting in the 15 minutes before midnight
}
}
\section{Changes in version 1.2-1 (release date: 9-12-2015)}{
\itemize{
  \item Literature reference for sleep detection updated
  \item Argument backup.cal.coef now with improved feedback if something goes wrong
  \item Report generation for part 4 much faster now
  \item Bug fixed in part 4 in assignment of dayname when a person sleeps during the day
  \item g.shell.GGIR now capable of handling minimal input argument specifications
  \item Console output from part 3 and 4 more compact now
}
}
\section{Changes in version 1.2-0 (release date: 27-10-2015)}{
\itemize{
  \item Package expanded with functions for detecting sleep (g.part3 and g.part4)
}
}
\section{Changes in version 1.1-5 (release date: 11-05-2015)}{
\itemize{
  \item Additional bugs fixs related to dealing with csv-format data from the Actigraph accelerometer brand
  \item g.part2 now also stores its output as milestone data just like g.part1. This to facilitate parallel processing of large  amounts of data files on clusters. \item The orginal report generation functionality in g.part2 has now been moved to shell function g.shell.GGIR because part3, 4 and 5 which are scheduled for later
  this year will combine milestone data from multiple analysis parts. It therefore, is more logical to control all report generation from the top level in the function hierarchy (g.shell.GGIR).
  \item g.part1 now comes with the option to provide backup calibration coefficients in case auto-calibration is unsuccessful. This function is primarily designed for studies with  short lasting experiments.
  \item g.part2 now has the option to export epoch values to a csv-file. Note that these same epoch values are also stored in the .RData milestone file from part2.
  The export option is mainly to ease access to epoch level data outside the R environment.
  \item g.shell.GGIR now offers the option to overwrite previously generated
  milestone data with variable 'overwrite'. The default setting (FALSE) is still to skip previously analysed files, which is intended to avoid having to do the same analyses twice after an interruption. However, overwriting previously generated milestone data could be useful when modifications are being made to the input arguments.
  \item g.shell.GGIR now offers the option to record the folderstructure in which
  an accelerometer file is located, especially useful for studies where accelerometer files are stored hierarchally in line with the study design.
}
}
\section{Changes in version 1.1-4 (release date: 06-11-2014)}{
\itemize{
  \item Additional bugs fixs related to recognising data format in Actigraph data
  \item Angle variables are now extracted based on 5 second rolling median as
  opposed to 501 sample rolling median
}
}
\section{Changes in version 1.1-3 (release date: 21-10-2014)}{
\itemize{
  \item Package expanded with functions: g.part1, g.part2, and g.shell.GGIR.
  These shell functions should help movement scientists
   to utilize the package without too much prior knowledge about R
  \item Additional bugs fixs related to recognising data format in Actigraph data
  \item Package expanded with axis angle metrics
  \item Package expanded with metric for replicating teLindert2013 counts, see
  \link{g.getmeta}
  \item Package expanded with metric ENMOa, see function tutorial g.getmeta
  }
}
\section{Changes in version 1.0-6 (release date: 1-9-2014)}{
\itemize{
  \item Bug fixed related to recognising date format in csv-file header from
  Actigraph accelerometer brand
  \item Literature reference added to function g.calibrate
  \item function g.getmeta expanded with argument 'chunksize'
  }
}
\section{Changes in version 1.0-4 (release date: 29-4-2014)}{
\itemize{
  \item Implemented functionality for csv-fromat data from GENEActiv and Actigraph.
  It seems to work for the test files I have, more testing may be necessary
  \item Cleaned up some of the NaN and NA output to aim for consistent annotation
  of missing data
  }
}
\section{Changes in version 1.0-3 (release date: 27-3-2014)}{
\itemize{
  \item Fixed bug in modified \link{g.analyse} in version 1.0-2
  }
}
\section{Changes in version 1.0-2 (release date: 14-3-2014)}{
\itemize{
  \item Expanded \link{g.analyse} with estimates of time spent in moderate and
  vigorous activity (a construct popular among physical activity researchers)
  \item Re-named a number of variables in the output from \link{g.analyse} to be
  more friendly for re-use in stata or sas. The majority of variable names are now
  shorter and do not include spaces, dots or commas
  }
}
\section{Changes in version 1.0-1 (release date: 29-1-2014)}{
\itemize{
  \item Fixed Linux-Windows sensitivty in \link{g.getmeta}. Certain damaged files
  can only be read with mmap.load set to FALSE in package GENEAread. Function
  \link{g.getmeta} in GGIR catches this problem and turns mmap.load to FALSE if
  necessary. This catch worked  well under Linux, but not for R in Windows.
  I have now fixed this
  }
}
\section{Changes in version 1.0-00 (release date: 21-12-2013)}{
\itemize{
	\item Added examples
  \item Expanded documentation for function \link{g.analyse}
	\item Fixed bug in extraction of starttime that caused the starttime to be
  truncated to  00:00:00 in a fraction of measurements.
	\item Fixed bug in extract of temperature in function \link{g.calibrate}
  \item Deleted three explorative variables that were only extracted in
  \link{g.analyse} if argument doangle was set to TRUE in function \link{g.getmeta}.
  A number of bugs and the lack of referable journal publications made me decide
  to remove these variables while working on them. I intend on re-releasing these
  variables during the course of 2014. Please contact me if this causes you problems
}
}
\section{Changes in version 0.6-17 (release date: 8-8-2013)}{
\itemize{
  \item This was the first version
}
}
\name{g.part5.wakesleepwindows}
\alias{g.part5.wakesleepwindows}
\title{
  Label wake and sleepperiod window
}
\description{
  Not intended for direct use by GGIR users.
  Label wake and sleepperiod window as part of \link{g.part5}.
}
\usage{
  g.part5.wakesleepwindows(ts, summarysleep_tmp2, desiredtz,
  nightsi, sleeplog, ws3new, Nts, ID, Nepochsinhour)
}
\arguments{

  \item{ts}{
    data.frame with time series
  }
  \item{summarysleep_tmp2}{
    cleaned output from part 4
  }
  \item{desiredtz}{
  }
  \item{nightsi}{
  }
  \item{sleeplog}{
  }
  \item{ws3new}{
  }
  \item{Nts}{
  }
  \item{ID}{
  }
  \item{Nepochsinhour}{
  }

}
\value{
  Object ts
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.analyse}
\alias{g.analyse}
\title{
  Function to analsyse meta-data generated by \link{g.getmeta} and \link{g.impute}
}
\description{
  Analyses the output from other functions within the packages to generate a basic
  descriptive summary for each accelerometer data file. Analyses include: Average 
  acceleration per day, per measurement, L5M5 analyses (assessment of the five
  hours with lowest acceleration and with highest acceleration). Further, the 
  traditionally popular variable MVPA is automatically extracted in six variants:
  without bout criteria in combination with epoch = epoch length as defined in
  g.getmeta (first value of the input argument windowsizes), 1 minute,
  and 5 minutes, and for bout durations 1 minute, 5 minutes or 10 minutes in
  combination with the epoch length as defined in g.getmeta.
}
\usage{
  g.analyse(I, C, M, IMP, params_247 = c(), params_phyact = c(),
            quantiletype = 7, includedaycrit = 16, 
            idloc = 1, snloc = 1, selectdaysfile=c(), 
            dayborder=0,  desiredtz = "", myfun=c(), ...)
}
\arguments{
  \item{I}{
    the output from function \link{g.inspectfile}
  }
  \item{C}{
    the output from function \link{g.calibrate}
  }
  \item{M}{
    the output from function \link{g.getmeta}
  }
  \item{IMP}{
    the output from function \link{g.impute}
  }
  \item{params_247}{
    See \link{g.part2}
  }
  \item{params_phyact}{
    See \link{g.part2}
  }
  \item{quantiletype}{
    type of quantile function to use (default recommended). For details, see 
    quantile function in STATS package
  }
  \item{includedaycrit}{
   See \link{g.part1}
  }
  \item{idloc}{
   See \link{g.part1}
  }
  \item{snloc}{
    If value = 1 (default) the code assumes that device serial number is stored in
    the obvious header field. If value = 2 the code uses the character string between
    the first and second character '_' in the filename as the serial number
  }
  \item{selectdaysfile}{
    See \link{g.part1}
  }
  \item{dayborder}{
    See \link{g.part1}
  }
  \item{desiredtz}{
    See \link{g.part1}
  }
  \item{myfun}{
    External function object to be applied to raw data, see \link{g.getmeta}.
  }
  \item{...}{
     Any argument used in the previous version of g.analyse, which will now
     be used to overrule the arguments specified with the parameter objects.
  }
}
\value{
  g.analyse generated two data,franeL
  \item{\code{summary}}{summary for the file that was analysed}
  \item{\code{daysummary}}{summary per day for the file that was analysed}
  These data.frames are used by function g.report.part2 to generate csv reports.
  An exaplantion of all the columns in the data.frame and subsequent csv reports
  can be found in the package vignette (Output part 2).
}
\examples{
  data(data.getmeta)
  data(data.inspectfile)
  data(data.calibrate)
  \dontrun{
    #inspect file:
    I = g.inspectfile(datafile)
    
    #autocalibration:
    C = g.calibrate(datafile) 
    
    #get meta-data:
    M = g.getmeta(datafile, desiredtz = "Europe/London", 
    windowsizes = c(5, 900, 3600),
    daylimit = FALSE, offset = c(0, 0, 0), 
    scale = c(1, 1, 1), tempoffset = c(0, 0, 0))
  }
  #impute meta-data:
  IMP = g.impute(M = data.getmeta, I = data.inspectfile)
  
  #analyse and produce summary:
  A = g.analyse(I = data.inspectfile, C = data.calibrate,
  M = data.getmeta, IMP)
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{create_test_acc_csv}
\alias{create_test_acc_csv}
\title{
Creates csv data file for testing purposes
}
\description{
Creates file in the Actigraph csv data format with dummy data
that can be used for testing. The file includes accelerometer data 
with bouts of higher acceleration, variations non-movement periods
in a range of accelerometer positions to allow for testing the
auto-calibration functionality.}
\usage{
  create_test_acc_csv(sf=3,Nmin=2000,storagelocation=c())
}
\arguments{
  \item{sf}{
    Sample frequency in Hertz, the default here is low to minimize file size
  }
  \item{Nmin}{
    Number of minutes (minimum is 2000)
  }
  \item{storagelocation}{
    Location where the test file named testfile.csv will be stored
    If no value is provided then the function uses the current 
    working directory
  }
}
\value{
 The function does not produce any output values. Only the file is
 stored
}
  
\examples{
  \dontrun{
    create_test_acc_csv()
  }
}
\name{g.part5.savetimeseries}
\alias{g.part5.savetimeseries}
\title{
  Saves oart 5 time series to csv files
}
\description{
  Not intended for direct use by GGIR users.
  Saves oart 5 time series to csv files as part of \link{g.part5}.
}
\usage{
  g.part5.savetimeseries(ts, LEVELS, desiredtz, rawlevels_fname,
  save_ms5raw_format="csv",
  save_ms5raw_without_invalid=TRUE,
  DaCleanFile=c(), includedaycrit.part5=2/3, ID=c())
}
\arguments{
  \item{ts}{
  }
  \item{LEVELS}{
  }
  \item{desiredtz}{
     See \link{g.getmeta}.
  }
  \item{rawlevels_fname}{
  }
  \item{save_ms5raw_format}{
  See \link{g.part5}
  }
  \item{save_ms5raw_without_invalid}{
  See \link{g.part5}
  }
  \item{DaCleanFile}{
  Content of data_cleaning_file as documented in \link{g.report.part5}.
  Only used in this function if save_ms5rawlevels is TRUE,  and it 
  only affects the time 
  series files stored.
  }
  \item{includedaycrit.part5}{
  See \link{g.report.part5}. Only used in this function if
  save_ms5rawlevels is TRUE,  and it only affects the time 
  series files stored.
  }
  \item{ID}{
  If data_cleaning_file is used then this argument specifies
  which participant ID the data correspond with.
  }
}
\value{
  Function does not provide output, it only prepare data for saving
  and saves it to a file.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.part5.handle_lux_extremes}
\alias{g.part5.handle_lux_extremes}
\title{
  Check lux values for extremes and imputes or removes them
}
\description{
  Extreme values are imputed by mean of neightbours if they occur isolated or
  in a sequence of two, and removed if they occure in a sequence of 3 or longer.
}
\usage{
  g.part5.handle_lux_extremes(lux)
}
\arguments{
  \item{lux}{
  Vector with lux values
  }
}
\value{
  List of imputed lux values and a vector with matching length named
  correction_log indicating which timestamps where imputed (value=1),
  replaced by NA (value=2) or untouched (value=0).
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.analyse.perfile}
\alias{g.analyse.perfile}
\title{
  Function supports \link{g.analyse}. Not intended for direct use by user.
}
\description{
  Generates recording specific analyses and fills corresponding
  output matrix, \link{g.analyse}.
}
\usage{
g.analyse.perfile(ID, fname, deviceSerialNumber,
  BodyLocation, startt, I, LC2, LD, dcomplscore,
  LMp, LWp, C, lookat, AveAccAve24hr, 
  colnames_to_lookat, QUAN,  ML5AD,
  ML5AD_names, igfullr, igfullr_names,
  daysummary, ds_names, includedaycrit, strategy, hrs.del.start,
  hrs.del.end, maxdur, windowsizes, idloc, snloc, wdayname, doquan,
  qlevels_names, doiglevels, tooshort, InterdailyStability, 
  IntradailyVariability,
  IVIS_windowsize_minutes, qwindow, longitudinal_axis_id)
}

\arguments{
\item{ID}{Person Identification number, this can be numeric or character}
\item{fname}{see \link{g.analyse.perday}}
\item{deviceSerialNumber}{As produced by \link{g.extractheadervars}}
\item{BodyLocation}{as produced by \link{g.extractheadervars}}
\item{startt}{First timestamp in metalong}
\item{I}{output \link{g.inspectfile}}
\item{LC2}{see \link{g.impute}}
\item{LD}{length data in minutes}
\item{dcomplscore}{see \link{g.impute}}
\item{LMp}{length measurement based on study protocol (minutes)}
\item{LWp}{length of sensor worn based on study protocol (minutes)}
\item{C}{output \link{g.calibrate}}
\item{lookat}{indices of metashort column to analyse}
\item{AveAccAve24hr}{Average acceleration in an average 24 hour cycle}
\item{colnames_to_lookat}{Names of columns to look at, corresponding 
to argurment lookat}
\item{QUAN}{Results quantile analysis on the average day produced by \link{g.analyse.avday}}
\item{ML5AD}{Results ML5 analyses on the average day produced by \link{g.analyse.avday}}
\item{ML5AD_names}{Columns names corresponding to ML5AD}
\item{igfullr}{Results intensity gradient (ig) analysis on
the average day produced by \link{g.analyse.avday}}
\item{igfullr_names}{Columns names corresponding to igfullr}
\item{daysummary}{object produced by \link{g.analyse.perday}}
\item{ds_names}{column names corresponding to daysummary}
\item{includedaycrit}{see \link{g.analyse}}
\item{strategy}{see \link{g.analyse}}
\item{hrs.del.start}{see \link{g.analyse}}
\item{hrs.del.end}{see \link{g.analyse}}
\item{maxdur}{see \link{g.analyse}}
\item{windowsizes}{see \link{g.getmeta}}
\item{idloc}{see \link{g.analyse}}
\item{snloc}{see \link{g.analyse}}
\item{wdayname}{character with weekdayname}
\item{doquan}{Boolean whether quantile analysis should be done}
\item{qlevels_names}{object produced by \link{g.analyse.avday}}
\item{doiglevels}{Boolean to indicate whether iglevels should be calculated}
\item{tooshort}{0 (file not too short) or 1 (file too short)}
\item{InterdailyStability}{see \link{g.IVIS}}
\item{IntradailyVariability}{see \link{g.IVIS}}
\item{IVIS_windowsize_minutes}{see \link{g.IVIS}}
\item{qwindow}{see \link{g.analyse}}
\item{longitudinal_axis_id}{Index of axis for which the angle correlates most 
  strongly across 24 hoursas calculated inside g.analyse. For hip worn 
  accelerometer this helps to check which axis was the veritcal axis. 
  The estimate may not be informative for other attachment locations.}

}
\value{
\item{\code{filesummary}}{summary for the file that was analysed}
\item{\code{daysummary}}{Summary per day for the file that was analysed}
}

\author{
Vincent T van Hees <v.vanhees@accelting.com>
}\name{read.gt3x_ggir}
\alias{read.gt3x_ggir}
\title{Read GT3X}
\usage{
read.gt3x_ggir(
  path,
  verbose = FALSE,
  asDataFrame = FALSE,
  imputeZeroes = FALSE,
  flag_idle_sleep = FALSE,
  cleanup = FALSE,
  ...,
  add_light = FALSE
)
}
\arguments{
\item{path}{Path to gt3x folder}

\item{verbose}{print diagnostic messages}

\item{asDataFrame}{convert to an \code{activity_df}, see
\code{as.data.frame.activity}}

\item{imputeZeroes}{Impute zeros in case there are missingness?
  Default is FALSE, in which case
  the time series will be incomplete in case there is missingness.}

\item{flag_idle_sleep}{flag idle sleep mode.  If \code{imputeZeroes = TRUE},
this finds where all 3 axes are zero.}

\item{cleanup}{should any unzipped files be deleted?}

\item{...}{additional arguments to pass to \code{parseGT3X} C++ code}

\item{add_light}{add light data to the \code{data.frame} if data exists in the
GT3X}
}
\value{
  A numeric matrix with 3 columns (X, Y, Z) and the following
  attributes:
  \itemize{
    \item \code{start_time} :  Start time from info file in \code{POSIXct} format.
    \item \code{subject_name} : Subject name from info file
    \item \code{time_zone} : Time zone from info file
    \item \code{missingness} : Named integer vector. Names are \code{POSIXct}
    timestamps and values are the number of missing values.
  }
}
\description{
  Read activity samples from a GT3X file as a matrix.
  Please note that all timestamps are in local time (of the device)
  even though they are represented as \code{POSIXct} with GMT timezone.
  
  The code in this function is a modified version of the read.gt3x in that it aids batch-loading of modern gt3x files. A pull request has been made to feed these enhancements back into the original code base
  https://github.com/THLfi/read.gt3x/pull/40. If and when merged we intend to deprecate the GGIR version of the code and make a direct dependency.
}
\note{
  The timestamps in the .gt3x data format are saved in .NET format, which is
  nanoseconds in local time since 0001-01-01.
  This is a bit tricky to parse into an R datetime format. DateTimes are
  therefore represented as \code{POSIXct} format with the
  'GMT' timezone attribute, which is false; the datetime actually
  represents local time.
}
\name{CalcSleepRegularityIndex}
\alias{CalcSleepRegularityIndex}
\title{
  Calculates Sleep Regularity Index
}
\description{
  Calculates Sleep Regularity Index per day pair proposed by Phillips and 
  colleagues in 2017 expanded with day-pair level estimates.
}
\usage{
  CalcSleepRegularityIndex(data = c(), epochsize = c(), desiredtz= c())
}
\arguments{
  \item{data}{
    Data.frame produced by function \link{g.sib.det}.
  }
  \item{epochsize}{
    Numeric value of epoch size in seconds.
  }
  \item{desiredtz}{
    Character with timezone database name, see also \link{g.getmeta}
  }
}
\value{
  Data.frame with columns: day (day number); Sleep Regularity Index, which by definition must lie in the 
  range -100 (reversed regularity), to 0 (random pattern), to 100 (perfect regularity);
  weekday (e.g. Wednesday); frac_valid, number between 0 and 1 indicating the fraction of the 24 hour
  period for which valid data was available in both the current and the next day, and; date.
}
\details{
  Calculates Sleep Regularity Index per day pair. Absense of missing data is not used 
  as a criteria for calculation. Instead the code asses the fraction of the time for which
  matching valid data points were found in both days. Later in g.part4 this fraction is used
  to include or exclude days based on the excludenightcrit criteria it also uses for the other
  sleep variables. In g.report.part4 these day-level SRI values are stored, but also
  aggregated across all recording days, all weekend days, and all weekend days, respectively. 
  Therefore, this function is broader in functionality than the algorithm proposed by Phillips and 
  colleagues in 2017.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\references{
\itemize{
  \item Andrew J. K. Phillips, William M. Clerx, et al. Irregular sleep/wake patterns 
  are associated with poorer academic performance and delayed circadian and 
  sleep/wake timing. Scientific Reports. 2017 June 12
}
}\name{g.getidfromheaderobject}
\alias{g.getidfromheaderobject}
\title{
  Extracts participant identifier from header object
}
\description{
  Extracts participant identifier from header object, if
  it can not be found then the filename is used as identifier.
  Function is not intended for direct interaction by package end user
  
}
\usage{
  g.getidfromheaderobject(filename,header,dformat,mon)	
}
\arguments{
  \item{filename}{
    File name
  }
  \item{header}{
    header object as extracted with \link{g.inspectfile}
  }
  \item{dformat}{
    Data format code, same as for \link{g.dotorcomma}
  }
  \item{mon}{
    Monitor code, same as for \link{g.dotorcomma}
  }
}
\value{
Participant identifier as character 
}
\examples{
\dontrun{
  g.getidfromheaderobject(filename="C:/myfile.bin",header,dformat=2,mon=2)	
}
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.part2}
\alias{g.part2}
\title{
  function to analyse and summarize pre-processed output from \link{g.part1}
}
\description{
  Loads the output from \link{g.part1} and then applies \link{g.impute} and
  \link{g.analyse}, after which the output is converted to .RData-format
  which will be used by \link{g.shell.GGIR} to generate reports.
  The variables in these reports are the same variables as described in
  \link{g.analyse}.
}
\usage{
  g.part2(datadir=c(), metadatadir=c(), f0=c(), f1=c(), 
        myfun=c(), params_cleaning = c(), params_247 = c(),
        params_phyact = c(), params_output = c(), params_general = c(), ...)
}
\arguments{
  \item{datadir}{
    Directory where the accelerometer files are stored or
    list, e.g. "C:/mydata" of accelerometer filenames
    and directories, e.g.
    c("C:/mydata/myfile1.bin", "C:/mydata/myfile2.bin").
  }
  \item{metadatadir}{
    Directory where the output from \link{g.part1} was stored
  }
  \item{f0}{
    File index to start with (default = 1). Index refers to the filenames sorted
    in increasing order
  }
  \item{f1}{
    File index to finish with (defaults to number of files available)
  }
  \item{myfun}{
    External function object to be applied to raw data.
    See details \link{applyExtFunction}.
  }
  \item{params_cleaning}{
    See \link{g.part1}
  }
  \item{params_247}{
    See details
  }
  \item{params_phyact}{
    See details
 }
  \item{params_output}{
    See details
  }
  \item{params_general}{
    See \link{g.part1}
  }
  \item{...}{
    To enable compatibility with R scripts written for older GGIR versions,
    the user can also provide parameters listed in the params_ objects as direct argument.
  }
}
\value{
  The function provides no values, it only ensures that other functions are called
  and that their output is stored in the folder structure as created with \link{g.part1}.
}
\details{
  GGIR comes with many processing parameters, which have been thematically grouped in
  parameter objects (R list). By running print(load_params()) you can
  see the default values of all the parameter objects. When g.part 2 is used via \link{g.shell.GGIR}
  you have the option to specifiy a configuration file, which will overrule the default
  parameter values. Further, as user you can set parameter values as input argument to both g.part2
  and \link{g.shell.GGIR}. Directly specified argument overrule the configuration file and default values.
  
  The parameter objects used by GGIR part 2 (g.part2) that are no already discussed in
  \link{g.part1} are:
  
  \subsection{params_output}{
    A list of parameters used to specify whether and how GGIR stores its output at various stages of the
    process.
    \describe{
      \item{storefolderstructure}{Boolean. Store folder structure of the accelerometer data.}
      \item{do.part3.pdf}{Boolean. In g.part3: Whether to generate a pdf for part 3 (default is TRUE).}
      \item{timewindow}{In g.part5: Timewindow over which summary statistics are derived.
        Value can be "MM" (midnight to midnight), "WW" (waking time to waking time), 
        or both c("MM","WW").}
      \item{save_ms5rawlevels}{Boolean, whether to save the time series 
        classification (levels) as a csv files.}
      \item{save_ms5raw_format}{Character string to specify how data should
        be stored: either "csv" (default) or "RData". Only used if save_ms5rawlevels=TRUE.}
      \item{save_ms5raw_without_invalid}{Boolean to indicate whether to remove
        invalid days from the time series output files. Only used if save_ms5rawlevels=TRUE.}
      \item{epochvalues2csv}{Boolean. If TRUE then epoch values are exported to a
        CSV spreadsheet. Here, non-wear time is imputed where possible (default = FALSE).}
      \item{do.sibreport}{Boolean. Applied in g.part5. Boolean to indicate whether
        to generate report for the sustained inactivity bouts (sib).}
      \item{do.visual}{Boolean. If g.part4 is run with do.visual == TRUE then
        the function will generate a pdf with a visual representation of the
        overlap between the sleeplog entries and the accelerometer detections.
        This can be used to visualy verify that the sleeplog entries do
        not come with obvious mistakes.}
      \item{outliers.only}{Boolean. Relevant for do.visual == TRUE. Outliers.only == FALSE
        will visualise all available nights in the data. Outliers.only == TRUE will visualise
        only for nights with a difference in onset or waking time
        larger than the variable of argument criterror.}
      \item{criterror}{Numeric. Relevant for do.visual == TRUE and outliers.only == TRUE.
        criterror specifies the number of minimum number of hours difference
        between sleep log and  accelerometer estimate for the night to be
        included in the visualisation.}
      \item{visualreport}{Boolean. If TRUE then generate visual report based on combined output from part 2
        and 4. This is in beta-version at the moment.}
      \item{viewingwindow}{Numeric. Centre the day as displayed around noon (value = 1) or
        around midnight (value = 2).}
      \item{week_weekend_aggregate.part5}{Boolean, see \link{g.report.part5}}
      \item{dofirstpage}{Boolean, see \link{g.plot5}}
      \item{timewindow}{Timewindow over which summary statistics are derived.
        Value can be "MM" (midnight to midnight), "WW" (waking time to waking time),
        or both c("MM","WW").}
    }
  }
  \subsection{params_phyact}{
    A list of parameters releated to physical activity as used in GGIRpart2 and GGIRpart5.
    \describe{
      \item{threshold.lig}{Numeric. In g.part5: Threshold for light physical activity to
        separate inactivity from light. Value can be one number or an array of multiple
        numbers, e.g. threshold.lig =c(30,40). If multiple numbers are entered then
        analysis will be repliced for each combination of threshold values. Threshold is
        applied to the first metric in the milestone data, so if you have only specified
        do.ENMO == TRUE then it will be applied to ENMO.}
      \item{threshold.mod}{Numeric. In g.part5: Threshold for moderate physical activity 
        to separate light from moderate. Value can be one number or an array of 
        multiple numbers, e.g. threshold.mod =c(100,110).
        If multiple numbers are entered then analysis will be repliced for each
        ombination of threshold values. Threshold is applied to the first metric in the
        milestone data, so if you have only specified do.ENMO == TRUE then it will be
        applied to ENMO.}
      \item{threshold.vig}{Numeric. In g.part5: Threshold for vigorous physical activity 
        to separate moderate from vigorous. Value can be one number or an array of 
        multiple numbers, e.g. threshold.mod =c(400,500). If multiple numbers are
        entered then analysis will be repliced for each combination of threshold values.
        Threshold is applied to the first metric in the milestone data, so if you
        have only specified do.ENMO == TRUE then it will be applied to ENMO.}
      \item{closedbout}{Boolean, see \link{g.getbout}}
      \item{frag.metrics}{Character, see \link{g.fragmentation}}
      \item{mvpathreshold}{Numeric, Acceleration threshold for MVPA estimation in GGIR part2.
        This can be a single number or an array of numbers, e.g. c(100,120). In the later case 
        the code will estimate MVPA seperately for each threshold. If this variable 
        is left blank c() then MVPA is not estimated}
      \item{boutcriter}{Numeric, The variable boutcriter is a number between
        0 and 1 and defines what fraction of a bout needs to be above the 
        mvpathreshold, only used in GGIR part 2}
      \item{mvpadur}{Numeric, default = c(1,5,10). Three bout duration for which 
        MVPA will be calculated. Only used in GGIR part 2}
      \item{bout.metric}{Numeric, Specify a metric for bout detection. A description of these
        bout metrics can be found in the new function \link{g.getbout}} 
      \item{boutdur.mvpa}{Numeric, see Durations of mvpa bouts in minutes to be extracted.
        The default values is c(1,5,10) and will start with the identification of
        10 minute bouts, followed by 5 minute bouts in the rest of the data, and followed
        by 1 minute bouts in the rest of the data.}
      \item{boutdur.in}{Numeric, see  Durations of inactivty bouts in minutes to be
        extracted. Inactivity bouts are detected in the segments of the data which
        were not labelled as sleep or MVPA bouts. The default duration values
        is c(10,20,30), this will start with the identification of 30 minute bouts,
        followed by 20 minute bouts in the rest of the data, and followed by 10 minute
        bouts in the rest of the data.}
      \item{boutdur.lig}{Numeric, see  Durations of light activty bouts in minutes
        to be extracted. Light activity bouts are detected in the segments of the data
        which were not labelled as sleep, MVPA, or inactivity bouts. The default
        duration values is c(1,5,10), this will start with the identification of 
        10 minute bouts, followed by 5 minute bouts in the rest of the data, and followed
        by 1 minute bouts in the rest of the data.}
      \item{boutcriter.in}{Numeric, see  A number between 0 and 1 and defines what fraction
        of a bout needs to be below the light threshold}
      \item{boutcriter.lig}{Numeric, see A number between 0 and 1 and defines what
        fraction of a bout needs to be between the light and moderage threshold}
      \item{boutcriter.mvpa}{Numeric, see A number between 0 and 1 and defines
        what fraction of a bout needs to be above the mvpathreshold}
    }
  }
  \subsection{params_247}{
    A list of parameters releated to description of 24/7 behaviours that do not fall
    under conventional physical activity or sleep outcomes, these parameters are used
    in GGIRpart2 and GGIRpart5: 
    \describe{
      \item{qwindow}{Numeric or character, To specify windows over which all
        variables are calculated, e.g. acceleration distirbution, number of valid
        hours, LXMX analysis, MVPA. If value = c(0,24), which is the default, 
        all variables will only be calculated over the full 24 hours in a day, If
        value =c(8,24) variables will be calculated over the window 0-8, 8-24 and 0-24.
        All days in the recording will be segmented based on these values. If you want
        to use a day specific segmentation then you can set qwindow to be
        the full path to activity diary file. See documentation \link{g.conv.actlog}
        for details.}
      \item{qwindow_dateformat}{Numeric, see \link{g.conv.actlog}}
      \item{M5L5res}{Numeric, resolution of L5 and M5 analysis in minutes 
        (default: 10 minutes)}
      \item{winhr}{Numeric, Vector of window size(s) (unit: hours) of L5 and
        M5 analysis (dedault = 5 hours)}
      \item{qlevels}{Numeric, array of percentiles for which value needs to be extracted.
        These need to be expressed as a fraction of 1, e.g. c(0.1, 0.5, 0.75). 
        There is no limit to the number of percentiles. If left empty then percentiles will not be 
        extracted. Distribution will be derived from short epoch metric data.}
      \item{ilevels}{Numeric, Levels for acceleration value frequency
        distribution in mg, e.g. c(0,100,200). There is no limit to the number of levels.}
      \item{window.summary.size}{Numeric, Functionality designed for the London Centre
        of Longidutinal studies. Size in minutes of the summary window}
      \item{iglevels}{Numeric, Levels for acceleration value frequency distribution
        in mg used for intensity gradient calculation (according to the method by 
        Rowlands 2018). By default this is argument is empty and the intensity gradient
        calculation is not done. The user can either provide a single value (any) to 
        make the intensity gradient use the bins c(seq(0,4000,by=25),8000) or the 
        user could specify their own distribution. There is no constriction to the 
        number of levels.}
      \item{IVIS_windowsize_minutes}{Numeric, see \link{g.IVIS}}
      \item{IVIS_epochsize_seconds}{depricated. Numeric, see \link{g.IVIS}}
      \item{IVIS.activity.metric}{Numeric, see \link{g.IVIS}}
      \item{qM5L5}{Numeric, see \link{g.getM5L5}}
      \item{MX.ig.min.dur}{Numeric, see \link{g.getM5L5}}
      \item{LUXthresholds}{Numeric. Vector with numeric sequece corresponding to
        the thresholds used to calculated time spent in LUX ranges.}
      \item{LUX_cal_constant}{Numeric, if both LUX_cal_constant and LUX_cal_exponent are 
        provided LUX LUX values are converted based on formula y = constant * exp(x * exponent)}
      \item{LUX_cal_exponent}{Numeric, if both LUX_cal_constant and LUX_cal_exponent are provided LUX
        LUX values are converted based on formula y = constant * exp(x * exponent)}
      \item{LUX_day_segments}{Numeric vector with hours at which the day should be segmented for
        the LUX analysis.}
      \item{L5M5window}{Argument depricated after version 1.5-24. 
        This argument used to define the start and end time, in 24 hour clock hours,
        over which L5M5 needs to be calculated. Now this is done with argument qwindow}
    }
  }
}
\examples{
  \dontrun{
    metadatadir = "C:/myresults/output_mystudy"
    g.part2(metadatadir)
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\references{
  \itemize{
    \item van Hees VT, Gorzelniak L, Dean Leon EC, Eder M, Pias M, et al. (2013) Separating
    Movement and Gravity Components in an Acceleration Signal and Implications for the
    Assessment of Human Daily Physical Activity. PLoS ONE 8(4): e61691.
    doi:10.1371/journal.pone.0061691
    \item van Hees VT, Fang Z, Langford J, Assah F, Mohammad A, da Silva IC, Trenell MI,
    White T, Wareham NJ, Brage S. Auto-calibration of accelerometer data for
    free-living physical activity assessment using local gravity and temperature:
    an evaluation on four continents. J Appl Physiol (1985). 2014 Aug 7
  }
}\name{g.calibrate}
\alias{g.calibrate}
\title{
  function to estimate calibration error and make recommendation for addressing it
}
\description{
  Function starts by identifying ten second windows of non-movement. Next, the
  average acceleration per axis per window is used to estimate calibration error 
  (offset and scaling) per axis. The function provides recommended correction factors
  to address the calibration error and a summary of the callibration procedure.
}
\usage{
  g.calibrate(datafile, params_rawdata = c(), params_general = c(),
              params_cleaning = c(), ...)
}
\arguments{
  \item{datafile}{
    Name of accelerometer file
  }
  \item{params_rawdata}{
    See \link{g.part1}
  }
  \item{params_general}{
    See \link{g.part1}
  }
  \item{params_cleaning}{
    See \link{g.part1}
  }
  \item{...}{
     Any argument used in the previous version of g.calibrate, which will now
     be used to overrule the arguments specified with the parameter objects.
  }  
}
\value{
 \item{\code{scale}}{scaling correction values, e.g. c(1,1,1) }
  \item{\code{offset}}{offset correction values, e.g. c(0,0,0)}
  \item{\code{tempoffset}}{correction values related to temperature, e.g. c(0,0,0)}
  \item{\code{cal.error.start}}{absolute difference between Euclidean norm during all
  non-movement windows and 1 g before autocalibration}
  \item{\code{cal.error.end}}{absolute difference between Euclidean norm during all 
  non-movement windows and 1 g after autocalibration}
  \item{\code{spheredata}}{average, standard deviation, Euclidean norm and temperature
  (if available) for all ten second non-movement windows as used for the
  autocalibration procedure}
  \item{\code{npoints}}{number of 10 second no-movement windows used to populate
  the sphere}
  \item{\code{nhoursused}}{number of hours of measurement data scanned to find
    the ten second time windows with no movement}
  \item{\code{meantempcal}}{mean temperature corresponding to the data as used for 
    autocalibration. Only applies to data collected with GENEActiv monitor. 
  }
}
\examples{
  \dontrun{
    datafile = "C:/myfolder/testfile.bin"
    
    #Apply autocalibration:
    C = g.calibrate(datafile)
    print(C$scale)
    print(C$offset)
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
  Zhou Fang
}
\references{
  \itemize{
    \item van Hees VT, Fang Z, Langford J, Assah F, Mohammad A, da Silva IC, Trenell MI, 
    White T, Wareham NJ, Brage S. Auto-calibration of accelerometer data for
    free-living physical activity assessment using local gravity and temperature: 
    an evaluation on four continents. J Appl Physiol (1985). 2014 Aug 7
  }
}\name{g.analyse.avy}
\alias{g.analyse.avday}
\title{
  Function supports \link{g.analyse}. Not intended for direct use by user.
}
\description{
  Generatess average day analyses and fills corresponding output
  matrix, \link{g.analyse}.
}
\usage{
  g.analyse.avday(doquan, averageday, M, IMP, t_TWDI, quantiletype,
                   ws3, doiglevels, firstmidnighti, ws2, midnightsi,
                   params_247 = c(), ...)
}
\arguments{
  \item{doquan}{Boolean whether quantile analysis should be done}
  \item{averageday}{ As produced by \link{g.impute}}
  \item{M}{ As produced by \link{g.getmeta}}
  \item{IMP}{ As produced by \link{g.impute}}
  \item{t_TWDI}{ Same as qwindow as described in \link{g.analyse}}
  \item{quantiletype}{see \link{g.analyse}}
  \item{ws3}{ Epoch size in seconds}
  \item{doiglevels}{Boolean to indicate whether iglevels should be calculated}
  \item{firstmidnighti}{see \link{g.detecmidnight}}
  \item{ws2}{see \link{g.weardec}}
  \item{midnightsi}{see \link{g.detecmidnight}}
  \item{params_247}{
    See \link{g.part2}
  }
  \item{...}{
   Any argument used in the previous version of g.analyse.avday, which will now
   be used to overrule the arguments specified with the parameter objects.
  }
}
\value{
  \item{\code{InterdailyStability}}{}
  \item{\code{IntradailyVariability}}{}
  \item{\code{igfullr_names}}{}
  \item{\code{igfullr}}{}
  \item{\code{QUAN}}{}
  \item{\code{qlevels_names}}{}
  \item{\code{ML5AD}}{}
  \item{\code{ML5AD_names}}{}
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.plot}
\alias{g.plot}
\title{
function to generate a plot for quality check purposes
}
\description{
Function takes meta-data as generated by \link{g.getmeta} and \link{g.impute}
to create a visual representation of imputed time periods
}
\usage{
g.plot(IMP, M, I, durplot)
}
\arguments{
  \item{IMP}{
  output from \link{g.impute}
}
  \item{M}{
  output from \link{g.getmeta}
}
  \item{I}{
  output from \link{g.inspectfile}
}
  \item{durplot}{
number of days to plot
}
}
\value{
function only produces a plot, no values
}
\examples{
\dontrun{
  #inspect file:
  I = g.inspectfile(datafile)
  
  #autocalibration:
  C = g.calibrate(datafile) 
  
  #get meta-data:
  M = g.getmeta(datafile)
}
data(data.getmeta)
data(data.inspectfile)

#impute meta-data:
IMP = g.impute(M = data.getmeta, I = data.inspectfile, strategy = 1,
hrs.del.start = 0, hrs.del.end = 0, maxdur = 0)

#plot data
g.plot(IMP, M = data.getmeta, I = data.inspectfile, durplot=4)
}

\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.imputeTimegaps}
\alias{g.imputeTimegaps}
\title{
  Impute gaps in three axis raw accelerometer data
}
\description{
  Removes all sample with a zero in each of the three axes, and then imputes time 
  gaps by the last recorded value per axis normalised to 1 _g_
}
\usage{
  g.imputeTimegaps(x, xyzCol, timeCol, sf, k=0.25)
}
\arguments{
  \item{x}{
    Data.frame with raw accelerometer data, and a timestamp column with millisecond resolution.
  }
  \item{xyzCol}{
    Columnnames or numbers for the x, y and z column
  }
  \item{timeCol}{
    Column name or number for the timestamp column
  }
  \item{sf}{
    Sample frequency in Hertz
  }
   \item{k}{
    Minimum time gap length to be imputed
  }
}
\value{
  Data.frame based on input x with timegaps inputed

}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.inspectfile}
\alias{g.inspectfile}
\title{
function to inspect accelerometer file for brand, sample frequency and header
}
\description{
Inspects accelerometer file for key information, including: monitor brand, sample frequency and file header
}
\usage{
g.inspectfile(datafile, desiredtz = "", params_rawdata = c(),
                         configtz = c(), ...)
}
\arguments{
  \item{datafile}{
  	name of data file
  }
  \item{desiredtz}{
  	Desired timezone, see documentation \link{g.getmeta}
  }
  \item{params_rawdata}{
    See \link{g.part1}
  }
  \item{configtz}{
    ...
  }
  \item{...}{
   Any argument used in the previous version of g.getmeta, which will now
    be used to overrule the arguments specified with the parameter objects.
  }
}
  
\value{
 \item{header}{fileheader}
  \item{monn}{monitor name (genea, geneactive)}
  \item{monc}{monitor brand code (0 - ad-hoc file format, 1 = genea (non-commercial),
      2 = GENEActive, 3 = actigraph, 4 = Axivity (AX3, AX6), 5 = Movisense, 6 = Verisense)}
   \item{dformn}{data format name, e.g bin, csv, cwa, gt3x}
  \item{dformc}{data format code (1 = .bin, 2 = .csv, 3 = .wav, 4 = .cwa, 5 = ad-hoc .csv, 6 = .gt3x)}
  \item{sf}{samplefrequency in Hertz}
  \item{filename}{filename}
}

\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{get_nw_clip_block_params}
\alias{get_nw_clip_block_params}
\title{
  Set monitor brand specific parameters
}
\description{
  Set monitor brand specific thresholds for non-wear detection, clipping
  etection, and blocksizes to be loaded.
  Not designed for direct use by user.
}
\usage{
  get_nw_clip_block_params(chunksize, dynrange, monc, rmc.noise=c(),
  sf, dformat, rmc.dynamic_range)
}
\arguments{
  \item{chunksize}{
    See \link{g.calibrate}
  }
  \item{dynrange}{
    See \link{g.getmeta}
  }
  \item{monc}{
    See \link{g.inspectfile}
  }
  \item{rmc.noise}{
    Noise level of acceleration signal in _g_-units, used when working ad-hoc data formats via \link{read.myacc.csv}.
  }
   \item{sf}{
    Numeric, sample frequency in Hertz
  }
  \item{dformat}{
    See \link{g.dotorcomma}
  }
  \item{rmc.dynamic_range}{
    Optional, please see \link{read.myacc.csv}
  }
  
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.part5.onsetwaketiming}
\alias{g.part5.onsetwaketiming}
\title{
  Identify wake and sleepperiod window timing
}
\description{
  Not intended for direct use by GGIR users.
  Labels timing of wakeing up and sleep onset as part of \link{g.part5}.
}
\usage{
  g.part5.onsetwaketiming(qqq, ts, min, sec, hour, timewindowi, skiponset, skipwake)
}
\arguments{
  \item{qqq}{
  }
  \item{ts}{
  }
  \item{min}{
  }
  \item{sec}{
  }
  \item{hour}{
  }
  \item{timewindowi}{
  }
  \item{skiponset}{
  }
  \item{skipwake}{
  }

}
\value{
  A list with objects: wake, onset, wakei, onseti, skiponset, and skipwake.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{getStartEnd}
\alias{getStartEnd}
\title{
Generate start and end time of a day
}
\description{
Generate start and end time of a day when working with argument selectdaysfile
in \link{g.part1}. The user provides a date and a start hour
which is used to generate the timestamps of the start hour minutes 5 minutes and the
start hour plus 24 hours. Function not designed for direct use by package user.
}
\usage{
getStartEnd(d, startHour, outputFormat = "\%d/\%m/\%Y \%H:\%M:\%S",
  tz = "Europe/London")
}
\arguments{
  \item{d}{
  character with date (without time) format
  }
   \item{startHour}{
  Hour that analysis starts at
  }
  \item{outputFormat}{
  Characterstring indicating outputFormat
  }
  \item{tz}{
  Same as desiredtz in \link{g.part1}
  }
}
\value{
Data.frame with two columns: a start time five minutes before startHour on day d
  and an endtime 24 hours after startHour
}

\examples{
startandendtime = getStartEnd(d="20/5/2017", startHour=4)
}
\author{
  Joe Heywood <j.heywood@ucl.ac.uk>
}\name{g.create.sp.mat}
\alias{g.create.sp.mat}
\title{
Converts sleep period information. Not intended for direct use
}
\description{
Function to convert data into sleep period matrix part of g.part4.R.
Not intended for direct use by package user
}
\usage{
g.create.sp.mat(nsp,spo,sleepdet.t,daysleep=FALSE)

}
\arguments{
  \item{nsp}{
  Integer indicating the number of sleep periods
  } 
 \item{spo}{
  Empty matrix with overview of sleep periods, 5 columns and
  as along as nps
  }
 \item{sleepdet.t}{
   Part of detected sleep from g.sib.det for one night and one 
   sleep definition
  }
 \item{daysleep}{
  Boolean to indicator whether this person woke up
  after noon (daysleeper)
  }
 
}
\value{
 \itemize{
\item spo matrix with start and end of each sleep period
\item calendardate date corresponding to the day on which the night started
\item item wdayname weekdayname
}
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.analyse.perday}
\alias{g.analyse.perday}
\title{
  Function supports \link{g.analyse}. Not intended for direct use by user.
}
\description{
  Generates day specific analyses and fills corresponding
  output matrix, \link{g.analyse}.
}
\usage{
g.analyse.perday(selectdaysfile, ndays, firstmidnighti, time, nfeatures, 
                midnightsi, metashort, averageday,
                doiglevels, nfulldays,lastmidnight, ws3, ws2, qcheck,
                fname, idloc, BodyLocation, wdayname, tooshort, includedaycrit,
                doquan, quantiletype, doilevels, domvpa,
                mvpanames, wdaycode, ID,
                deviceSerialNumber, ExtFunColsi, myfun, desiredtz = "",
                params_247 = c(), params_phyact = c(),
                ...)
}
\arguments{
  \item{selectdaysfile}{see \link{g.analyse}}
  \item{ndays}{Number of days in file} 
  \item{firstmidnighti}{see \link{g.detecmidnight}}
  \item{time}{timestamp column from metalong converted to character}
  \item{nfeatures}{estimate of number of variables that need to be stored in
  the output matrix}
  \item{midnightsi}{see \link{g.detecmidnight}}
  \item{metashort}{see \link{g.impute}}
  \item{averageday}{As produced by \link{g.impute}}
  \item{doiglevels}{Boolean to indicate whether iglevels should be calculated}
  \item{nfulldays}{Number of days between the first and last midnight in the recording}
  \item{lastmidnight}{see \link{g.detecmidnight}}
  \item{ws3}{Epoch size in seconds}
  \item{ws2}{see \link{g.weardec}}
  \item{qcheck}{vector with zeros and ones for each epoch, respenting the
    quality check derived with g.impute}
  \item{fname}{RData filename produced by g.part1}
  \item{idloc}{see \link{g.analyse}}
  \item{BodyLocation}{as produced by \link{g.extractheadervars}} 
  \item{wdayname}{character with weekdayname}
  \item{tooshort}{0 (file not too short) or 1 (file too short)} 
  \item{includedaycrit}{see \link{g.analyse}} 
  \item{doquan}{Boolean whether quantile analysis should be done}
  \item{quantiletype}{see \link{g.analyse}} 
  \item{doilevels}{Boolean whether to generate ilevels, see \link{g.analyse}} 
  \item{domvpa}{Boolean whether to do mvpa analysis}
  \item{mvpanames}{Matrix with 6 columns and 1 row holding the names for the six 
    mvpa variables}
  \item{wdaycode}{Equal to M$wday as produced by \link{g.getmeta}}
  \item{ID}{Person Identification number, this can be numeric or character}
  \item{deviceSerialNumber}{As produced by \link{g.extractheadervars}}
  \item{ExtFunColsi}{column index of metashort where metric is stored}
  \item{myfun}{External function object to be applied to raw data, see \link{g.getmeta}.}
  \item{desiredtz}{see \link{g.part1}}
  \item{params_247}{
    See \link{g.part2}
  }
  \item{params_phyact}{
    See \link{g.part2}
  }
  \item{...}{
   Any argument used in the previous version of g.analyse.perday, which will now
   be used to overrule the arguments specified with the parameter objects.
  }
}
\value{
  \item{\code{daysummary}}{Summary per day for the file that was analysed}
  \item{\code{ds_names}}{Variable names in daysummary}
  \item{\code{windowsummary}}{Window summary, only used when
  selectdayfile is specified}
  \item{\code{ws_names}}{Variable names in windowsummary}
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.report.part2}
\alias{g.report.part2}
\title{
Generate report from milestone data produced by \link{g.part2}
}
\description{
Creates report from milestone data produced by \link{g.part2}. Not intended
  for direct use by package user
}
\usage{
g.report.part2(metadatadir=c(), f0=c(), f1=c(), maxdur = 0, 
selectdaysfile=c(), store.long=FALSE)
}
\arguments{
  \item{metadatadir}{
  see \link{g.part2}
  }
  \item{f0}{
  see \link{g.part2}
  }
  \item{f1}{
  see \link{g.part2}
  }
  \item{maxdur}{
  see \link{g.part2}
  }
  \item{selectdaysfile}{
  see \link{g.part2}
  }
  \item{store.long}{
  Booelean to indicate whether output should stored in long
  format in addition to default wide format. Automatically turned
  to TRUE if using day segmentation with qwindow.
  }
}
\value{
Function does not produce data, but only writes reports
in csv format and visual reports in pdf format
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{numUnpack}
\alias{numUnpack}
\docType{package}
\title{
  Simple function using Rcpp
}
\description{
  Simple function using Rcpp
}
\usage{
  numUnpack(pack)	
}
\arguments{
  \item{pack}{
    vector of integer
  }
}
\examples{
  \dontrun{
    numUnpack()
  }
}\name{g.weardec}
\alias{g.weardec}
\title{
 Detects whether accelerometer is worn
}
\description{
Uses the object produced by \link{g.part1} to assess
whether the accelerometer was worn
}
\usage{
g.weardec(M,wearthreshold,ws2)
}
\arguments{
  \item{M}{
  Object produced by \link{g.getmeta}
  }
   \item{wearthreshold}{
  Number of axis that at least need to meet the non-wear criteria
  }
  \item{ws2}{
  Large windowsize used in seconds to apply non-wear detection
  Small window size not needed, because this is inherent to the object M
  }
  
}
\value{
\itemize{
\item \code{r1} Participant id extracted from file
\item \code{r2} Night number
\item \code{r3} Detected onset of sleep expressed as hours 
since the previous midnight
\item \code{LC} fraction of 15 minute windows with more than 5 percent
clipping
\item \code{LC2} fraction of 15 minute windows with more than 80
percent clipping
}

}
\examples{
data(data.getmeta)
output = g.weardec(M=data.getmeta,wearthreshold=2,ws2=3600)
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{HASPT}
\alias{HASPT}
\title{
  Heuristic Algorithms estimating SPT window.
}
\description{
  As used in function \link{g.sib.det}. Function is not intended for direct use
  by GGIR user.
}
\usage{
HASPT(angle, perc = 10, spt_threshold = 15, sptblocksize = 30, 
      spt_max_gap = 60, ws3 = 5, constrain2range = FALSE,
      HASPT.algo="HDCZA", invalid, HASPT.ignore.invalid=FALSE)
}
\arguments{
  \item{angle}{
    Vector of epoch level estimates of angle
  }
  \item{perc}{
    Number to indicate percentage threshold (default 10 corresponds to 2018
    paper)
  }
  \item{spt_threshold}{
    Numeric threshold used in HASPT algorithm (default 15 corresponds to 
    2018 paper)
  }
  \item{sptblocksize}{
    Number to indicate minimum SPT block size (minutes)
  }
  \item{spt_max_gap}{
    Number to indicate maximum gap (minutes) in SPT window blocks.
  }
  \item{ws3}{
    Number representing epoch length in seconds
  }
  \item{constrain2range}{
    Bolean to indicate whether threshold should be constrained to a range
  }
  \item{HASPT.algo}{
    Character to indicate what algortihm should be used. Default "HDCZA" is 
    Heuristic algorithm looking at Distribution of Change in Z-Angle as
    described in van Hees et al. 2018. Other options included:
    "HorAngle", which is based on HDCZA but replaces non-movement detection of 
    the HDCZA algorithm by looking for time segments where the angle of the 
    longitudinal sensor axis has an angle relative to the horizontal plane
    between -45 and +45 degrees.
  }
  \item{invalid}{
    Integer vector with per epoch an indicator of valid(=0) or invalid(=1) epoch.
  }
  \item{HASPT.ignore.invalid}{
    Boolean to indicate whether invalid time segments should be ignored
  }
}
\value{
  List with start and end times of the SPT window and the threshold as used.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.applymetrics}
\alias{g.applymetrics}
\title{
  Extract metrics from acceleration signals
}
\description{
  Function to extract metrics from acceleration signal. Not intended
  for direct use by user
}
\usage{
  g.applymetrics(data, sf, ws3, metrics2do,
                          n = 4, lb = 0.2, hb = 15)
}
\arguments{
  \item{data}{
    Three column matrix with x, y, and z acceleration data
  }
  \item{n}{
    filter order, only needed if a metric is selected
    that involves a frequency filter
  }
  \item{sf}{
     sample frequency
  }
   \item{ws3}{
     Epoch size in seconds
  }
  \item{metrics2do}{
    Dataframe with Boolean indicator for all metrics whether
    they should be extracted or not. For instance,
    metrics2do$do.bfen = TRUE, indicates that the bfen metric
    should be extracted
  }
  \item{lb}{
    Lower boundery of cut-off frequencies
  }
  \item{hb}{
    Higher boundery of cut-off frequencies
  }
  
}
\value{
Dataframe with metric values in columns average per epoch (ws3)
}

\examples{
  Gx = runif(n=10000,min=0,max=2)
  Gy = runif(n=10000,min=1,max=3)
  Gz = runif(n=10000,min=0,max=2)
  data = cbind(Gx, Gy, Gz)
  metrics2do = data.frame(do.bfen=TRUE,do.enmo=TRUE,do.lfenmo=FALSE,
  do.en=FALSE,do.hfen=FALSE,do.hfenplus=FALSE,do.mad=FALSE,do.anglex=FALSE,
  do.angley=FALSE,do.anglez=FALSE,do.roll_med_acc_x=FALSE,
  do.roll_med_acc_y=FALSE,do.roll_med_acc_z=FALSE,
  do.dev_roll_med_acc_x=FALSE,do.dev_roll_med_acc_y=FALSE,
  do.dev_roll_med_acc_z=FALSE,do.enmoa=FALSE,
  do.lfx=FALSE, do.lfy=FALSE, do.lfz=FALSE, 
  do.hfx=FALSE, do.hfy=FALSE, do.hfz=FALSE, 
  do.bfx=FALSE, do.bfy=FALSE, do.bfz=FALSE,
  do.zcx=FALSE, do.zcy=FALSE, do.zcz=FALSE, do.brondcounts=FALSE)
  
  extractedmetrics = g.applymetrics(data,n=4,sf=40,ws3=5,metrics2do)
}

\author{
Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.conv.actlog}
\alias{g.conv.actlog}
\title{
  Function to read activity log and make it useful for the rest of GGIR.
}
\description{
 Function to read activity log and convert it into data.frame
 that has for each ID and date a different qwindow vector
  
}
\usage{
  g.conv.actlog(qwindow, qwindow_dateformat="\%d-\%m-\%Y")	
}
\arguments{
  \item{qwindow}{
    Path to csv file with activity log. Expected format of the activity diary is:
    First column headers followed by one row per recording, first column is recording ID,
    which needs to match with the ID GGIR extracts from the accelerometer file.
    Followed by date column in format "23-04-2017", where date format is specified by 
    argument qwindow_dateformat (below). Use the character combination date, Date or 
    DATE in the column name. This is followed by 
    one or multiple columns with start times for the activity types in that day format in 
    hours:minutes:seconds. The header of the column will be used as label for each activity
    type. Insert a new date column before continuing with activity types for next day.
    Leave missing values empty. If an activitylog is used then individuals who do 
    not appear in the activitylog will still be processed with value c(0,24).
    Dates with no activiy log data can be skipped, no need to have a column with the 
    date followed by a column with the next date.
  }
  \item{qwindow_dateformat}{
    Character specifying the date format used in the activity log.
  }
}
\value{
Data.frame with column ID, date and qwindow, where each
qwindow value is a qwindow vector
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}
\name{resample}
\alias{resample}
\docType{package}
\title{
  Simple function using Rcpp
}
\description{
  Simple function using Rcpp
}
\usage{
  resample(raw, rawTime, time, stop, type=1)	
}
\arguments{
  \item{raw}{
    stop-by-3 matrix with raw values of x, y and z.
  }
  \item{rawTime}{
    vector with stop elements of raw time.
  }
  \item{time}{
    array with required time points.
  }
  \item{stop}{
    Number of rows in raw
  }
  \item{type}{
    integer to indicate type of interpolation, 1=linear, 2=nearest neighbour
  }
}
\examples{
  \dontrun{
    resample()
  }
}
\name{g.part5}
\alias{g.part5}
\title{
  Merge output from physical activity and sleep analysis into one report
}
\description{
  Function to merge the output from g.part2 and g.part4 into one report enhanced with
  profiling of sleep and physical activity stratified across intensity levels and
  based on bouted periods as well as non-bouted periods.
}
\usage{
g.part5(datadir = c(), metadatadir = c(), f0=c(), f1=c(),
                   params_sleep = c(), params_metrics = c(),
                   params_247 = c(), params_phyact = c(), 
                   params_cleaning = c(), params_output = c(),
                   params_general = c(), ...)
}
\arguments{
  \item{datadir}{
    Directory where the accelerometer files are stored or list of accelerometer
    filenames and directories
  }
  \item{metadatadir}{
    Directory that holds a folders 'meta' and inside this a folder 'basic' which
    contains the milestone data produced by g.part1. The folderstructure
    is normally created by g.part1 and g.shell.GGIR will recognise what the value
    of metadatadir is.
  }
  \item{f0}{
    File index to start with (default = 1). Index refers to the filenames sorted
    in increasing order
  }
  \item{f1}{
    File index to finish with (defaults to number of files available)
  }
  \item{params_sleep}{
    List of parameters used for sleep analysis (GGIR part 3, 4, and 5): see documentation \link{g.part3}.
  }
  \item{params_metrics}{
    List of parameters used for metrics extraction (GGIR part 1): see documentation \link{g.part1}.
  }
  \item{params_247}{
    List, see \link{g.part2}
  }
  \item{params_phyact}{
    List, see \link{g.part2}
  }
  \item{params_cleaning}{
    List, see \link{g.part1}
  }
  \item{params_output}{
    List, see \link{g.part2}
  }
  \item{params_general}{
    List, see \link{g.part1}
  }
  \item{...}{
    To enable compatibility with R scripts written for older GGIR versions,
    the user can also provide parameters listed in the params_ objects as direct argument.
  }
}
\value{
  The function does not produce values but generates an RData file
  in the milestone subfolder ms5.out which incudes a dataframe
  named \code{output}. This dataframe is used in g.report.part5 to create
  two reports one per day and one per person. See package vignette
  paragraph "Output part 5" for description of all the variables.
}
\examples{
  \dontrun{
    metadatadir = "C:/myfolder/meta"
    g.part5(metadatadir=metadatadir)
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{parseGT3Xggir}
\alias{parseGT3Xggir}
\title{Parse activity samples from a GT3X file with batch loading}
\usage{
  parseGT3Xggir(
    filename,
    max_samples,
    scale_factor,
    sample_rate,
    start_time,
    batch_begin = 0L,
    batch_end = 0L,
    verbose = FALSE,
    debug = FALSE,
    impute_zeroes = FALSE
  )
}
\arguments{
\item{filename}{(char*) path to a log.bin file inside the unzipped gt3x folder, which contains the activity samples}

  \item{max_samples}{Maximum number of rows to parse. The returned matrix will always contain this number of rows, having zeroes if
  not data is found.}

\item{scale_factor}{Scale factor for the activity samples.}

\item{sample_rate}{sampling rate for activity samples.}

\item{start_time}{starting time of the sample recording.}

\item{batch_begin}{first second in time relative to start of raw non-imputed recording to include in this batch}

\item{batch_end}{last second in time relative to start of raw non-imputed recording to include in this batch}

\item{verbose}{Print the parameters from the log.bin file and other messages?}

\item{debug}{Print information for every activity second}

\item{impute_zeroes}{Impute zeros in case there are missingness?}
}
\value{
  Returns a matrix with max_samples rows and 3 columns with the acceleration samples.
  The matrix has attributes
  "time_index", "missingness", "start_time_log", "sample_rate", "impute_zeroes".
}
\description{
  Parse activity samples from a GT3X file.
  The code in this function is a modified version of the read.gt3x in that it aids batch-loading of modern gt3x files. A pull request has been made to feed these enhancements back into the original code base
  https://github.com/THLfi/read.gt3x/pull/40. If and when merged we intend to deprecate the GGIR version of the code and make a direct dependency.
}
\name{updateBlocksize}
\alias{updateBlocksize}
\title{
Update blocksize of data to be read depending on available memory.
}
\description{
Function queries available memory to either lower or increase the blocksize
used by function \link{g.readaccfile}
}
\usage{
updateBlocksize(blocksize, bsc_qc)
}
\arguments{
  \item{blocksize}{
Number of filepages (binary data) or rows (other dataformats).
}
\item{bsc_qc}{
  Data.frame with columns time (timestamp from Sys.time) and size (memory size).
  This is used for housekeeping in \link{g.calibrate} and \link{g.getmeta}
}
}
\value{
List with blocksize and bsc_qc, same format as input, although bsc_qc has one new
row.
}\name{g.wavread}
\alias{g.wavread}
\title{
function to read .wav files as produced by the accelerometer named 'Axivity'

}
\description{
For reading the wav accelerometer data as collected with an Axivity accelerometer 
}
\usage{
g.wavread(wavfile, start = 1, end = 100,units="minutes")
}
\arguments{
  \item{wavfile}{
filename (required)
}
  \item{start}{
 start point for reading data, see also units
}
  \item{end}{
end point for reading data, see also units
}
\item{units}{
units used for defining start and end
}
}
\details{
If only \code{start} is defined then \code{g.binread} will read all data beyond
\code{start} until the end of the file is reached

}
\value{
 \item{\code{rawxyz}}{matrix with raw x, y, and, z acceleration values}
  \item{\code{header}}{file header}
  \item{\code{timestamps}}{local timestamps for \code{rawxyz}}
 }
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.getM5L5}
\alias{g.getM5L5}
\title{
Extract M5 and L5 from time series
}
\description{
Extract M5 and L5 from time series, function used by \link{g.analyse} and
not intended for direct use by package user. Please see \link{g.analyse}
for further clarification on functionalities
}
\usage{
g.getM5L5(varnum,ws3,t0_LFMF,t1_LFMF,M5L5res,winhr,qM5L5=c(), 
iglevels=c(), MX.ig.min.dur=10)	
}
\arguments{
  \item{varnum}{
  Numeric vector of epoch values
}
  \item{ws3}{
  Small epoch size in seconds

}
\item{t0_LFMF}{
Start hour of the day for the M5L5 analyses, e.g. 0 for midnight
}
\item{t1_LFMF}{
End hour of the day for the M5L5 analyses, e.g. 24 for midnight
}
\item{M5L5res}{
Resolution of hte M5L5 analyses in minutes 

}
\item{winhr}{
windowsize of M5L5 analyses, e.g. 5 hours
}
\item{qM5L5}{
Percentiles (quantiles) to be calculated over L5 and M5 window.
}
\item{iglevels}{
 See  \link{g.analyse}. If provided then the intensity gradient will be calculated
 for all MX windows larger or equal than argument MX.ig.min.dur
}
\item{MX.ig.min.dur}{
  Minimum MX duration needed in order for intensity gradient to be calculated
}

}
\value{
\itemize{
\item DAYL5HOUR = Starting time in hours of L5
\item DAYL5VALUE = average acceleration during L5
\item DAYM5HOUR = Starting time in hours of M5
\item DAYM5VALUE = average acceleration during M5
\item V5NIGHT = average acceleration between 1am and 6am
}
}
\examples{
data(data.getmeta)
g.getM5L5 = function(varnum=data.getmeta,ws3=5,t0_LFMF=0,
t1_LFMF=24,M5L5res=10,winhr=5)
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}\name{data.calibrate}
\alias{data.calibrate}
\docType{data}
\title{
Example output from g.calibrate
}
\description{
data.calibrate is example output from \link{g.calibrate}
}
\usage{data(data.calibrate)}
\format{
  The format is:
 chr "data.calibrate"
}
\source{
The data was collected on one individual for testing purposes
}
\examples{
data(data.calibrate)
}
\keyword{datasets}
\name{g.report.part4}
\alias{g.report.part4}
\title{
Generate report from milestone data produced by \link{g.part4}
}
\description{
Creates report from milestone data produced by \link{g.part4}. Not intended
  for direct use by package user
}
\usage{
g.report.part4(datadir=c(),metadatadir=c(),loglocation = c(),f0=c(),
  f1=c(),storefolderstructure=TRUE, data_cleaning_file=c(), sleepwindowType = "SPT")
}
\arguments{
  \item{datadir}{
    see \link{g.part4}
  }
  \item{metadatadir}{
    see \link{g.part4}
  }
  \item{loglocation}{
    see \link{g.part4}
  }
  \item{f0}{
    see \link{g.part4}
  }
  \item{f1}{
    see \link{g.part4}
  }
  \item{storefolderstructure}{
    see \link{g.part4}
  }
  \item{data_cleaning_file}{
    see \link{g.part4}
  }
  \item{sleepwindowType}{
    see \link{g.part4}
  }
}
\value{
Function does not produce data, but only writes reports
in csv format and a visual report in pdf.

The following files are stored in the root of the results folder:
part4_nightsummary_sleep_cleaned.csv
part4_summary_sleep_cleaned.csv

The following files are stored in the folder results/QC:
part4_nightsummary_sleep_full.csv
part4_summary_sleep_full.csv

If a sleeplog is used *_full.csv as stored in the QC folder includes estimates 
for all nights in the data, and *_cleaned.csv in the results folder includes
estimates for all nights in the data excluding the nights that did not had a
sleeplog entry or had no valid accelerometer data.

If a sleep log is not used then * _cleaned.csv includes the nights that are
in *_full.csv excluding the nights with insufficient data.

If you have a study where the sleeplog was available for a subset of the
participants, but you want to include all individuals in your analysis, 
then use the *_full.csv output and clean the night level data yourself
by excluding  rows with cleaningcode > 1 which are the cases where no
or invalid accelerometer data was present.

The above means that for studies with missing sleeplog entries for some
individuals and some nights using the *_full.csv output and excluding
rows (nights) with cleaningcode > 1 will lead to the same as
* _cleaned.csv plus sleep estimates for the nights with missing
sleeplog, providing that there was enough accelerometer data for 
those nights.

In other words, *_cleaned.csv is perfect if you only want to rely on 
nights with a sleeplog or if you do not use a sleeplog at all. For all
other scenarios We advise using the *_full.csv report and to clean it 
yourself.

See package vignette sections "Sleep analysis" and "Output part 4"
for a more elaborative description of the sleep analysis and reporting.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.sibreport}
\alias{g.sibreport}
\title{
  Generate sustiained inactivty bouts report
}
\description{
  Generate sustained inactivity bout report. Function not intended
  for direct use by package user
}
\usage{
  g.sibreport(ts, ID, epochlength, logs_diaries=c(), desiredtz="")
}
\arguments{
  \item{ts}{
    Data frame with time series as created inside function \link{g.part5}
  }
  \item{ID}{
   Recording identifier (character or numeric)
  }
  \item{epochlength}{
    Numeric to indicate epoch length in seconds in the ts object
  }
  \item{logs_diaries}{
  Object produced by \link{g.loadlog} function
  }
  \item{desiredtz}{
    See \link{g.getmeta}
  }
}
\value{
  Dataframe with one row per sustained inactivity bout and corresponding 
  properties stored in the data.frame columns.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{POSIXtime2iso8601}
\alias{POSIXtime2iso8601}
\title{
Convert POSIX to iso8601 timestamp
}
\description{
To avoid ambiguities when sharing and comparing timestamps. All timestamps
are expressed in iso8601 format: https://en.wikipedia.org/wiki/ISO_8601
}
\usage{
POSIXtime2iso8601(x,tz)	
}
\arguments{
  \item{x}{
Vector of timestamps in POSIX format
}
  \item{tz}{
 Timezone of data collection, e.g. "Europe/London".
 See
 https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
 for full list
}
}
  
\examples{
\dontrun{
x ="2017-05-07 13:15:17 CEST"
tz = "Europe/Amsterdam"
x_converted = POSIXtime2iso8601(x,tz)
}
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.sib.plot}
\alias{g.sib.plot}
\title{
  Create plot of sustained inactivity bouts
}
\description{
  Function create plot of sustained inactivity bouts for quality
  check purposes as part of \link{g.part3}. Not intended for direct use by package user
}
\usage{
  g.sib.plot(SLE, M, I, plottitle, nightsperpage=7, desiredtz="")
}
\arguments{
  \item{SLE}{
    Output from \link{g.sib.det}
  }
  \item{M}{
    Output from \link{g.getmeta}
  }
  \item{I}{
    Output from \link{g.inspectfile}
  }
  \item{plottitle}{
    Title to be used in the plot
  }
  \item{nightsperpage}{
     Number of nights to show per page
  }
  \item{desiredtz}{
    See \link{g.part3}
  }
  
}
\value{
Function has no output other than the plot
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.IVIS}
\alias{g.IVIS}
\title{
  Calculates IV and IS
}
\description{
  To extract interdaily stability and interdaily variability as originally proposed by
  van Someren.
}
\usage{
  g.IVIS(Xi, epochsizesecondsXi = 5, IVIS_epochsize_seconds = c(),
  IVIS_windowsize_minutes = 60, IVIS.activity.metric = 1)
}
\arguments{
  \item{Xi}{
  Vector with acceleration values, e.g. ENMO metric.
  }
  \item{epochsizesecondsXi}{
  Epoch size of the values in Xi expressed in seconds.
  }
  \item{IVIS_epochsize_seconds}{ 
   This argument has been depricated.
  }
  \item{IVIS_windowsize_minutes}{
    Window size of the Intradaily Variability (IV) and Interdaily
  Stability (IS) metrics in minutes, needs to be able to add up to 24 hours.
  }
  \item{IVIS.activity.metric}{
  Metric used for activity calculation.
  Value = 1, uses continuous scaled acceleration.
  Value = 2, tries to collapse acce;eration into a binary score of rest
  versus active to try to similate the original approach.
  }
}
\value{
  \item{InterdailyStability}{}
  \item{IntradailyVariability}{}

}
\examples{
  Xi = abs(rnorm(n = 10000,mean = 0.2))
  IVISvariables = g.IVIS(Xi=Xi)
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\references{
\itemize{
  \item Eus J. W. Van Someren, Dick F. Swaab, Christopher C. Colenda, Wayne Cohen, W. Vaughn McCall & Peter B. Rosenquist.
  Bright Light Therapy: Improved Sensitivity to Its Effects on Rest-Activity Rhythms in Alzheimer Patients by Application of Nonparametric Methods/
  Chronobiology International. 1999. Volume 16, issue 4.
}
}
\name{g.loadlog}
\alias{g.loadlog}
\title{
  Load and clean sleeplog information
}
\description{
  Loads sleeplog from a csv input file and applies sanity checks
  before storing the output in a dataframe
}
\usage{
  g.loadlog(loglocation=c(),coln1=c(),colid=c(),nnights=c(),
    sleeplogidnum=TRUE, sleeplogsep=",", meta.sleep.folder = c(),
  desiredtz="")
}
\arguments{
  \item{loglocation}{
    Location of the spreadsheet (csv) with sleep log information. 
    See package vignette for explanation on expected format
  }
  \item{coln1}{
    Column number in the sleep log spreadsheet where the onset of the first
  night starts
  }
  \item{colid}{
      Column number in the sleep log spreadsheet in which the participant
  ID code is stored (default = 1)

  }
  \item{nnights}{
    Number of nights for which sleep log information should be available.
  It assumes that this is constant within a study. If sleep log information
  is missing for certain nights then leave these blank
  
  }
  \item{sleeplogidnum}{
    Should the participant identifier as stored in the sleeplog be
    interpretted as a number (TRUE=default) or a character (FALSE)?
  }
  \item{sleeplogsep}{
    Value used as sep argument for reading sleeplog csv file.
  }
  \item{meta.sleep.folder}{
    Path to part3 milestone data, only specify if sleeplog is in advanced format.
  }
  \item{desiredtz}{
    See \link{g.part4}
  }
}
\value{
  Data frame with sleeplog, which can be either in basic format or in advanced 
  format. See GGIR package vignette for discussion of these two formats.
}
\examples{
\dontrun{
  sleeplog = g.loadlog(loglocation="C:/mysleeplog.csv",coln1=2,
  colid=1,nnights=5,sleeplogidnum=TRUE)
}
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.part5.definedays}
\alias{g.part5.definedays}
\title{
  Fix missing night in part 4 output
}
\description{
  Not intended for direct use by GGIR users.
  Defines when day windows start and end as part of \link{g.part5}.
}
\usage{
  g.part5.definedays(nightsi, wi, indjump, nightsi_bu, 
                              ws3new, qqq_backup=c(), ts, Nts, timewindowi, Nwindows)
}
\arguments{
  \item{nightsi}{
  }
  \item{wi}{
  }
  \item{indjump}{
  }
  \item{nightsi_bu}{
  }
  \item{ws3new}{
  }
  \item{qqq_backup}{
  }
  \item{ts}{
  }
  \item{Nts}{
  }
  \item{timewindowi}{
  }
  \item{Nwindows}{
  }

}
\value{
  List of qqq and qqq_backup 
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.part5.classifyNaps}
\alias{g.part5.classifyNaps}
\title{
  Classify Naps from identified sustained inactivty bouts
}
\description{
  Classify Naps from identified sustained inactivty bouts, based on model
  that was originally trained with hip-worn accelerometer data in 3-3.5 year olds.
  Assume that metric ENMO is used and HASIB.algo is set to vanHees2015.
}
\usage{
  g.part5.classifyNaps(sibreport = c(), desiredtz = "", 
        possible_nap_window = c(9, 18),
        possible_nap_dur = c(15, 240),
        nap_model = "hip3yr", HASIB.algo = "vanHees2015")
}
\arguments{
  \item{sibreport}{
    Object generated by \link{g.sibreport}
  }
  \item{desiredtz}{
    See \link{g.getmeta}.
  }
  \item{possible_nap_window}{
    Numeric vector of length two with range in clock hours during which naps are
    assumed to take place.
  }
  \item{possible_nap_dur}{
   Numeric vector of length two with range in duration (minutes) of a nap.
  }
  \item{nap_model}{
    Character to specify classification model. Currently the only option is "hip3yr", which
    corresponds to a model trained with hip data in 3-3.5 olds trained with parent diary data.
  }
  \item{HASIB.algo}{
    See \link{g.part3}.
  }
}
\value{
  Data.frame with classified naps and newly detected non-wear periods.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{getStartEndNumeric}
\alias{getStartEndNumeric}
\title{
Generate start and end page of a day
}
\description{
Generate start and end page of a day when working with argument selectdaysfile
in \link{g.part1}. The user provides a date and a start hour
which is used to generate the pages of the start hour minutes 5 minutes and the
start hour plus 24 hours. Function not designed for direct use by package user.
}
\usage{
  getStartEndNumeric(d, hhr, startHour = 4)
}
\arguments{
  \item{d}{
    Character with date (without time) format
  }
  \item{hhr}{
    GENEActiv::header.info(f) output
  }
  \item{startHour}{
    Hour that analysis starts at
  }
}
\value{
  Data.frame with two columns: a start page five minutes before startHour
  on day d and an end page 24 hours after startHour
}

\examples{
\dontrun{
hhr = GENEActiv::header.info("C:/myfile.bin")
mystartandendpage = getStartEndNumeric(d="20/5/2017", hhr, startHour = 4)
}
}
\author{
  Joe Heywood <j.heywood@ucl.ac.uk>
}\name{g.convert.part2.long}
\alias{g.convert.part2.long}
\title{
  Convert part 2 report to long format
}
\description{
  Not for direct access by used. This function is used 
  inside g.report.part2 and convert2 part 2 report to long 
  ormat if there
  are multiple segments per day
}
\usage{
g.convert.part2.long(daySUMMARY)
}
\arguments{
  \item{daySUMMARY}{
  Object available inside g.report.part2
  }
}
\value{
  Data.frame with long format version of daySUMMARY
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{applyExtFunction}
\alias{applyExtFunction}
\title{
  Apply external function to acceleration data.
}
\description{
  Applies external function to the raw acceleration data within GGIR.
  This makes it easier for new algorithms developed to be pilotted
  on accelerometer data while taking advantage of the existing comprehensive GGIR
  data management and analysis infrastructure.
  This function is not for direct interaction by user, please supply object 
  \code{myfun} to \link{g.shell.GGIR} or \link{g.part1}. Object \code{myfun}
  is a list as detailed below.
}
\usage{
  applyExtFunction(data, myfun, sf, ws3, interpolationType=1)
  
}
\arguments{
  \item{data}{
    Data data.frame as present internally in \link{g.getmeta}. It has at least
    four columns of which the first is the timestamp followed by the x, y,
    and z acceleration.
  }
  \item{myfun}{
    See details, in short: myfun is a list object that holds the external function
    to be applied to the data and various parameters to aid in the process.
  }
  \item{sf}{
    Sample frequency (Hertz) of the data object
  }
  \item{ws3}{
    Short epoch size (first value of windowsizes in \link{g.getmeta}).
  }
  \item{interpolationType}{
    See \link{g.getmeta}
  }
}
\value{
  The output of the external algorithm aggregated or repeated to fit the
  short epoch length of GGIR. Therefore, the short epoch length of GGIR
  should be a multitude of the resolution of the external function output,
  or visa versa.
}
\details{
    See package vignette for detailed tutorial with examples
    on how to use the function embedding:
    https://cran.r-project.org/web/package=GGIR/vignettes/applyExtFunction.pdf
    Function applyExtFunction
    is typically not used by the GGIR user directly.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.getstarttime}
\alias{g.getstarttime}
\title{
  Extract start time of a measurement
}
\description{
  Extract start time of a measurement. GGIR calculates all timestamps by
  using the first timestamp and sample frequency. Not intended
  for direct use by package user
}
\usage{
  g.getstarttime(datafile, P, header, mon, dformat, desiredtz,
  selectdaysfile)
}
\arguments{
  \item{datafile}{
   Full path to data file
  }
  \item{P}{
    Object extracted with \link{g.readaccfile}
  }
  \item{header}{
    File header extracted with \link{g.inspectfile}
  }
  \item{mon}{
    Same as in \link{g.dotorcomma}
  }
  \item{dformat}{
    Same as in \link{g.dotorcomma}
  }
  \item{desiredtz}{
    Same as in \link{g.dotorcomma}
  }
  \item{selectdaysfile}{
    See \link{g.part1}
  }
}
\value{
  The starttime
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.binread}
\alias{g.binread}
\title{
function to read binary files as produced by the accelerometer named 'Genea',
not to be confused with the 'GENEActiv' (see package GENEAread for this)
}
\description{
For reading the binary data as collected with a Genea accelerometer 
(Unilever Discover, UK). For reading GENEActive binary data, see package GENEAread.
}
\usage{
g.binread(binfile, start = 0, end = 0)
}
\arguments{
  \item{binfile}{
filename (required)
}
  \item{start}{
 start point for reading data, this can either be a timestamp
 "year-month-day hr:min:sec" or a page number (optional)
}
  \item{end}{
end point for reading data, this can either be a timestamp 
"year-month-day hr:min:sec" or a page number (optional)
}
}
\details{
If only \code{start} is defined then \code{g.binread} will read all data beyond
\code{start} until the end of the file is reached

}
\value{
 \item{\code{rawxyz}}{matrix with raw x, y, and, z acceleration values}
  \item{\code{header}}{file header}
  \item{\code{timestamps1}}{timestamps for \code{rawxyz} in seconds since 1970-01-01 00:00}
  \item{\code{timestamps2}}{timestamps for \code{rawxyz} in day time format }
  \item{\code{batt.voltage}}{matrix with battery voltage and corresponding timestamps}
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
Jing Hua Zhao <jinghua.zhao@mrc-epid.cam.ac.uk>
}\name{g.dotorcomma}
\alias{g.dotorcomma}
\title{
Assesses whether decimals in fileheader are stored with comma
or dot separated decimals
}
\description{
The function is used by \link{g.readaccfile} to assess how numeric data
should be interpretted
}
\usage{
g.dotorcomma(inputfile,dformat,mon, desiredtz = "", ...)
}
\arguments{
  \item{inputfile}{
  full path to inputfile
  }
   \item{dformat}{
  Data format code: 1=.bin, 2=.csv, 3=.wav, 4=.cwa, 5=.csv for ad-hoc monitor 
  brand
  }
  \item{mon}{
  Monitor code (accelorometer brand): 0=undefined, 1=GENEA, 2=GENEActiv,
  3=Actigraph, 4=Axivity, 5=Movisense, 6=Verisense
  }
  \item{desiredtz}{
  Desired timezone, see documentation \link{g.getmeta}
  }
  \item{...}{
   Any input arguments needed for function \link{read.myacc.csv} if you
   are working with a non-standard csv formatted files.
  }
  
}
\value{
Character object showing how decimals are separated
}

\examples{
\dontrun{
decn = g.dotorcomma(inputfile="C:/myfile.bin",dformat=1,mon=2)
}
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.readtemp_movisens}
\alias{g.readtemp_movisens}
\title{
Reads the temperature from movisens files.
}
\description{
  Reads the temperature from movisens files, resamples it and adds 
  it to the matrix where accelerations are stored
}
\usage{
g.readtemp_movisens(datafile, desiredtz = "", from = c(), to = c(),
  interpolationType=1)
}
\arguments{
  \item{datafile}{
    Full path to the folder where the movisens bin files are stored. Note that 
    movisens store a set of bin file in one folder per recording. GGIR will read 
    the pertinent bin file to access to the temperature data.
  }
  \item{desiredtz}{
    See \link{g.getmeta}
  }
  \item{from}{
    Origin point to derive the temperature from movisens files (automatically 
    calculated by GGIR)
  }
  \item{to}{
    End point to derive the temperature from movisens files (automatically 
    calculated by GGIR)
  }
  \item{interpolationType}{
    Integer, see \link{g.getmeta}
  }
}
\value{
  Data matrix with the temperature values resampled at 64 Hz.
}
\examples{
\dontrun{
  P = g.readtemp_movisens(datafile, desiredtz = "", from = c(), to = c())
}
}
\name{g.abr.day.names}
\alias{g.abr.day.names}
\title{
  Abbreviates daynames to numbers, needed for report generation in
  \link{g.plot5}
}
\description{
  Abbreviates daynames Monday becomes MON and Sunday becomes SUN
}
\usage{
  g.abr.day.names(daynames)	
}
\arguments{
  \item{daynames}{
    Vector of daynames in character format
  }
}

\examples{
  daynames = c("Monday","Friday")
  daynames_converted = g.abr.day.names(daynames)
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.cwaread}
\alias{g.cwaread}
\title{
Function to read .cwa-format files as produced by the accelerometer named 'Axivity'
}
\description{
For reading .cwa-format data, if you have .wav format data then see function \link{g.wavread}
}
\usage{
g.cwaread(fileName, start = 0, end = 0, progressBar = FALSE, 
  desiredtz = "", configtz = c(), interpolationType=1)
}
\arguments{
  \item{fileName}{
    filename (required)
  }
  \item{start}{
    start point for reading data, this can either be a timestamp
    "year-month-day hr:min:sec" or a page number (optional)
  }
  \item{end}{
    end point for reading data, this can either be a timestamp 
    "year-month-day hr:min:sec" or a page number (optional)
  }
  \item{progressBar}{
    Is trigger to switch on/off the text progress bar. If progressBar
    is TRUE then the function displays the progress bar but it works
    slightly slower
  }
  \item{desiredtz}{
    Desired timezone, see documentation \link{g.getmeta}
  }
  \item{configtz}{
    Only functional for AX3 cwa data at the moment. Timezone in which the accelerometer
    was configured. Only use this argument if the timezone of configuration and
    timezone in which recording took place are different.
  }
  \item{interpolationType}{
    See \link{g.getmeta}
  }
}
\value{
 \item{\code{data}}{dataframe with timestamp, raw x, -y, and, -z acceleration values,
  temperature, battery and light}
  \item{\code{header}}{file header}
}
\author{
  Evgeny Mirkes <em322@leicester.ac.uk>
  Vincent van Hees <v.vanhees@accelting.com>
}\name{data.inspectfile}
\alias{data.inspectfile}
\docType{data}
\title{
Example output from g.inspectfile
}
\description{
data.inspectfile is example output from \link{g.inspectfile}
}
\usage{data(data.inspectfile)}
\format{
  The format is:
 chr "data.inspectfile"
}
\source{
The data was collected on one individual for testing purposes
}
\examples{
data(data.inspectfile)
}
\keyword{datasets}
\name{g.part4_extractid}
\alias{g.part4_extractid}
\title{
 Extracts ID from filename and finds matching rows in sleeplog
}
\description{
  Extracts ID from filename and finds matching rows in sleeplog. Function not
  designed for direct use by GGIR users.
}
\usage{
g.part4_extractid(idloc, fname, dolog, sleeplogidnum, sleeplog)
}
\arguments{
  \item{idloc}{
   See \link{g.part4}
  }
  \item{fname}{
   Full patth to filename
  }
  \item{dolog}{
    Boolean to indicate whether to rely on a sleeplog
  }
   \item{sleeplogidnum}{
  Should the participant identifier as stored in the sleeplog be
  interpretted as a number (TRUE=default) or a character (FALSE)?
  }
  \item{sleeplog}{
   Sleeplog data.frame passed on from g.part4
  }
}
\value{
 List with accid the ID and matching_indices_sleeplog a vector with matching row
 indices in the sleeplog
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{check_myfun}
\alias{check_myfun}
\title{
  Checks myfun object before it is passed to applyExtfunction
}
\description{
  Checks that object myfun is a list and check the elements of the
  list for: that element names are as expected, that value of each
  element is of the expected type and length.
}
\usage{
  check_myfun(myfun, windowsizes)
  
}
\arguments{
  \item{myfun}{
    See \link{applyExtFunction}
  }
  \item{windowsizes}{
    See \link{g.getmeta}).
  }
}
\value{
  0 if all checkes passed, 1 if one or more checks did not pass. Error
  message are printed to the console with feedback on which checks did 
  not pass.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.part1}
\alias{g.part1}
\title{
  function to load and pre-process acceleration files
}
\description{
  Calls function \link{g.getmeta} and \link{g.calibrate}, and converts the
  output to .RData-format which will be the input for \link{g.part2}. Here,
  the function generates a folder structure to keep track of various output files.
  The reason why these \link{g.part1} and \link{g.part2} are not merged as one
  generic shell function is because g.part1 takes much
  longer to and involves only minor decisions of interest to the movement scientist.
  Function g.part2 on the other hand is relatively fast and comes with all the
  decisions that directly impact on the variables that are of interest to the
  movement scientist. Therefore, the user may want to run g.part1 overnight
  or on a computing cluster, while g.part2 can then be the main playing ground
  for the movement scientist. Function \link{g.shell.GGIR} provides the main shell
  that allows for operating g.part1 and g.part2.
}
\usage{
  g.part1(datadir = c(), outputdir = c(), f0 = 1, f1 = c(),
          studyname = c(), myfun = c(), params_metrics = c(), params_rawdata = c(),
          params_cleaning = c(), params_general = c(), ...)
}

\arguments{
  \item{datadir}{
    Directory where the accelerometer files are stored or list of accelerometer
    filenames and directories
  }
  \item{outputdir}{
    Directory where the output needs to be stored. Note that this function will
    attempt to create folders in this directory and uses those folder to organise
    output
  }
  \item{f0}{
    File index to start with (default = 1). Index refers to the filenames sorted
    in increasing order
  }
  \item{f1}{
    File index to finish with (defaults to number of files available)
  }
  \item{studyname}{
    If the datadir is a folder then the study will be given the name of the
    data directory. If datadir is a list of filenames then the studyname will be used
    as name for the analysis
  }
  \item{myfun}{
    External function object to be applied to raw data.
    See details \link{applyExtFunction}.
  }
  \item{params_metrics}{
   See details
  }
  \item{params_rawdata}{
   See details
  }
  \item{params_cleaning}{
   See details
  }
  \item{params_general}{
   See details.
  }
  \item{...}{
    Any input arguments needed for function \link{read.myacc.csv} if you
    are working with a non-standard csv formatted files. To enable compatibility with R 
    scripts written for older GGIR versions, the user can also provide parameters 
    listed in the params_ objects as direct argument.
  }
}
\details{
  GGIR comes with many processing parameters, which have been thematically grouped in
  parameter objects (R list). By running print(load_params()) you can
  see the default values of all the parameter objects. When g.part 1 is used via \link{g.shell.GGIR}
  you have the option to specifiy a configuration file, which will overrule the default
  parameter values. Further, as user you can set parameter values as input argument to both g.part1
  and \link{g.shell.GGIR}. Directly specified argument overrule the configuration file and default values.
  
  GGIR part 1 (g.part1) takes the following parameter objects as input:
  
  \subsection{params_metrics}{
  A list of parameters used to specify the signal metrics that need to be extract in GGIR part 1.
  \describe{
      \item{do.anglex}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.angley}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.anglez}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.zcx}{Boolean, if TRUE calculate metric zero-crossing count for x-axis.
        For computation specifics see source code of function \link{g.applymetrics}}
      \item{do.zcy}{Boolean, if TRUE calculate metric zero-crossing count for y-axis.
        For computation specifics see source code of function \link{g.applymetrics}}
      \item{do.zcz}{Boolean, if TRUE calculate metric zero-crossing count for z-axis.
        For computation specifics see source code of function \link{g.applymetrics}}
      \item{do.enmo}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.lfenmo}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.en}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.mad}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.enmoa}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.roll_med_acc_x}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.roll_med_acc_y}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.roll_med_acc_z}{Boolean, if TRUE calculate metric.
        For computation specifics see source code of function \link{g.applymetrics}}
      \item{do.dev_roll_med_acc_x}{Boolean, if TRUE calculate metric.
        For computation specifics see source code of function \link{g.applymetrics}}
      \item{do.dev_roll_med_acc_y}{Boolean, if TRUE calculate metric.
        For computation specifics see source code of function \link{g.applymetrics}}
      \item{do.dev_roll_med_acc_z}{Boolean, if TRUE calculate metric.
        For computation specifics see source code of function \link{g.applymetrics}}
      \item{do.bfen}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.hfen}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.hfenplus}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.lfen}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.lfx}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.lfy}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.lfz}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.hfx}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.hfy}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.hfz}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.bfx}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.bfy}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.bfz}{Boolean, if TRUE calculate metric. For computation
        specifics see source code of function \link{g.applymetrics}}
      \item{do.brondcounts}{Boolean, if TRUE calculate metric via R
        package activityCounts. We call them BrondCounts because
        there are large number of acitivty counts in the physical activity and
        sleep research field. By calling them Brond Counts we clarify that 
        these are the counts proposed by Jan Brond and implemented in R by
        Ruben Brondeel. The Brond Counts are intended to be an imitation of
        one the counts produced by one of the closed source ActiLife software
        by ActiGraph.}
      \item{lb}{Numeric, lower boundary of the frequency filter (in Hertz) as
        used in the filter-based metrics.}
      \item{hb}{Numeric, higher boundary of the frequency filter (in Hertz) as
        used in the filter-based metrics.}
      \item{n}{Numeric, order of the frequency filter as used in a variety of metrics.}
   }
  }
  
  \subsection{params_rawdata}{
    A list of parameters used to related to reading and pre-processing 
    raw data, excluding parameters related to metrics as those are in
    the params_metrics object.
    \describe{
      \item{backup.cal.coef}{Character. Default value is "retrieve".
        Option to use backed-up calibration coefficient instead of
        deriving the calibration coefficients when analysing the same file twice.
        Argument backup.cal.coef has two usecase. Use case 1: If the auto-calibration
        fails then the user has the option to provide back-up
        calibration coefficients via this argument. The value of the argument needs to
        be the name and directory of a csv-spreadsheet with the following column names
        and subsequent values: 'filename' with the names of accelerometer files on which
        the calibration coefficients need to be applied in case auto-calibration fails;
        'scale.x', 'scale.y', and 'scale.z' with the scaling coefficients; 'offset.x',
        'offset.y', and 'offset.z' with the offset coefficients, and;
        'temperature.offset.x', 'temperature.offset.y', and 'temperature.offset.z'
        with the temperature offset coefficients. This can be useful for analysing
        short lasting laboratory experiments with insufficient sphere data to perform
        the auto-calibration, but for which calibration coefficients can be derived
        in an alternative way.  It is the users responsibility to compile the
        csv-spreadsheet. Instead of building this file the user can also
        Use case 2: The user wants to avoid performing the auto-calibration repeatedly
        on the same file. If backup.cal.coef value is set to "retrieve" (default) then
        GGIR will look out for the  data_quality_report.csv  file in the outputfolder
        QC, which holds the previously generated calibration coefficients. If you
        do not want this happen, then deleted the data_quality_report.csv from the
        QC folder or set it to value "redo".}
      \item{minimumFileSizeMB}{Numeric. Minimum File size in MB required to enter processing,
        default 2MB. This argument can help
        to avoid having short uninformative files to enter the analyses. Given that a typical accelerometer
        collects several MBs per hour, the default setting should only skip the very tiny files.}
      \item{do.cal}{Boolean. Whether to apply auto-calibration or not by \link{g.calibrate}. Default and
        recommended setting is TRUE.}
      \item{imputeTimegaps}{Boolean to indicate whether timegaps larger than 1 sample should be imputed.
        Currently onlly used for .gt3x data and ActiGraph .csv format, where timegaps can be expected as a result of
        Actigraph's idle sleep.mode configuration that is turned on in some studies.}
      \item{spherecrit}{The minimum required acceleration value (in g) on both sides of 0 g
        for each axis. Used to judge whether the sphere is sufficiently populated}
      \item{minloadcrit}{The minimum number of hours the code needs to read for the
        autocalibration procedure to be effective (only sensitive to multitudes of 12 hrs, 
        other values will be ceiled). After loading these hours only extra data is loaded 
        if calibration error has not been reduced to under 0.01 g.}
      \item{printsummary}{Boolean. If TRUE will print a summary when done}
      \item{chunksize}{Numeric. Value between 0.2 and 1 to specificy the size of chunks to be 
        loaded as a fraction of a 12 hour period, e.g. 0.5 equals 6 hour chunks.
        The default is 1 (12 hrs). For machines with less than 4Gb of RAM memory a value
        below 1 is recommended.}
      \item{dynrange}{Numeric, provide dynamic range for accelerometer data to
        overwrite hardcoded 6 g for GENEA and 8 g for other brands}
      \item{interpolationType}{Integer to indicate type of interpolation to be used
        when resampling time series (mainly relevant for Axivity sensors),
        1=linear, 2=nearest neighbour}
      \item{all arguments that start with "rmc.".}{see function \link{read.myacc.csv}}
    }
  }   
  \subsection{params_cleaning}{
    A list of parameters used across all GGIR parts releated to masking or 
    imputing data, abbreviated as 'cleaning'.
    \describe{
      \item{do.imp}{Boolean. Whether to impute missing values (e.g. suspected of monitor non-wear) or not
        by \link{g.impute} in GGIR part2. Default and recommended setting is TRUE}
      \item{TimeSegments2ZeroFile}{Character. Path to csv-file holding the data.frame used for argument
        TimeSegments2Zero in function \link{g.impute}}
      \item{data_cleaning_file}{Character. Optional path to a csv file you create that holds four
        columns: ID, day_part5, relyonguider_part4, and night_part4. ID should hold the participant ID.
        Columns day_part5 and night_part4 allow you to specify which day(s) and
        night(s) need to be excluded from part 5 and 4, respectively. So, this will be done regardless
        of whether the rest of GGIR thinks those day(s)/night(s)
        are valid. Column relyonguider_part4 allows you to specify for which nights
        part 4 should fully rely on the guider. See also package vignette.}
      \item{excludefirstlast.part5}{Boolean. If TRUE then the first and last window 
      (waking-waking or midnight-midnight) are ignored in part 5.}
      \item{excludefirstlast}{Boolean. If TRUE then the first and last night of the measurement are
        ignored for the sleep assessment (part 4).}
      \item{excludefirst.part4}{Boolean. If TRUE then the first night of the measurement are
        ignored for the sleep assessment (part 4.}
      \item{excludelast.part4}{Boolean. If TRUE then the last night of the measurement are
        ignored for the sleep assessment.}
      \item{includenightcrit}{Numeric. Minimum number of valid hours per night (24 hour window between
        noon and noon), used for sleep assessment (part 4).}
      \item{minimum_MM_length.part5}{Numeric. Minimum length in hours of a MM day to be included
      in the cleaned part 5 results.}
      \item{selectdaysfile}{Character, Functionality designed for the London Centre
        of Longidutinal studies. Csv file holding the relation between device 
        serial numbers and measurement days of interest.}
      \item{strategy}{Numeric, how to deal with knowledge about study protocol.
        value = 1 means select data based on \code{hrs.del.start}, \code{hrs.del.end}, 
        and \code{maxdur}. Value = 2 makes that only the data between the first
        midnight and the last midnight is used for imputation. Value = 3 only selects
        the most active X days in the file where X is specified by argument \code{ndayswindow}. 
        Value = 4 to only use the data after the first midnight. Used in GGIR part 2}
      \item{hrs.del.start}{Numeric, how many HOURS after start of experiment did wearing
        of monitor start? Used in GGIR part 2}
      \item{hrs.del.end}{Numeric, how many HOURS before the end of the experiment did 
        wearing of monitor definitely end? Used in GGIR part 2}
      \item{maxdur}{Numeric, How many DAYS after start of experiment did experiment
        definitely stop? (set to zero if unknown = default). Used in GGIR part2}
      \item{ndayswindow}{Numeric,  If \code{strategy} is set to 3 then this is the 
        size of the window as a number of days. Used in GGIR part2}
      \item{includedaycrit.part5}{Numeric. see \link{g.report.part5}}
      \item{includedaycrit}{Numeric, minimum required number of valid hours
        in day specific analysis (NOTE: there is no minimum required number of
        hours per day in the summary of an entire measurement, every available 
        hour is used to make the best possible inference on average metric value
        per average day)}
    }
  }
   \subsection{params_general}{
    A list of parameters used across all GGIR parts that do not fall in any of the other
    categories.
    \describe{
      \item{overwrite}{Boolean. Do you want to overwrite analysis for which milestone data exists?
        If overwrite=FALSE then milestone data from a previous analysis will
        be used if available and visual reports will not be created again.}
      \item{selectdaysfile}{Character. Do not use, this is legacy code for one specific data study.
        Character pointing at a csv file holding the relationship between device serial
        numbers (first column) and measurement dates of interest
        (second and third column). The date format should be dd/mm/yyyy. And the first row
        if the csv file is assumed to have a character variable names, e.g. "serialnumber"
        "Day1" and "Day2" respectively. Raw data will be extracted and stored in the output
        directory in a new subfolder named 'raw'.}
      \item{dayborder}{Numeric. Hour at which days start and end (default = 0), 
      value = 4 would mean 4 am}
      \item{do.parallel}{Boolean. whether to use multi-core processing
        (only works if at least 4 CPU cores are available).}
      \item{maxNcores}{Numeric. Maximum number of cores to use when argument do.parallel is set to true.
        GGIR by default uses the maximum number of available cores, but this argument
        allows you to set a lower maximum.}
      \item{acc.metric}{Boolean. Which one of the metrics do you want to consider to analyze L5.
        The metric of interest need to be calculated in M.}
      \item{part5_agg2_60seconds}{Boolean. Wether to use aggregate epochs to 60 seconds
        as part of the part 5 analysis.}
      \item{print.filename}{Boolean. Whether to print the filename before before analysing
        it (default is FALSE). Printing the filename can be useful to investigate
        problems (e.g. to verify that which file is being read).}
      \item{desiredtz}{Character, desired timezone: see also http://en.wikipedia.org/wiki/Zone.tab}
      \item{configtz}{Character, Only functional for AX3 cwa data at the moment. 
        Timezone in which the accelerometer was configured. Only use this argument
        if the timezone of configuration and timezone in which recording took
        place are different.}
      \item{sensor.location}{Character, see \link{g.sib.det}}
      \item{acc.metric}{Character, see \link{g.sib.det}}
      \item{windowsizes}{Numeric vector, three values to indicate the lengths of the 
        windows as in c(window1,window2,window3): window1 is the short epoch length
        in seconds and by default 5 this is the time window over which acceleration and
        angle metrics are calculated, window2 is the long epoch length in seconds 
        for which non-wear and signal clipping are defined, default 900. However, 
        window3 is the window length of data used for non-wear detection and by default
        3600 seconds. So, when window3 is larger than window2 we use overlapping windows,
        while if window2 equals window3 non-wear periods are assessed by non-overlapping
        windows. Window2 is expected to be a multitude of 60 seconds.}
      \item{idloc}{Numeric. If idloc = 1 (default) the code assumes that ID
        number is stored in the obvious header field. Note that for ActiGraph data
        the ID never stored in the file header.  For value set to 2, 5, 6, and 7, GGIR
        looks at the filename and extracts the character string preceding the first 
        occurance of a '_', ' ' (space), '.' (dot), and '-', respecitvely. You may have 
        noticed that idloc 3 and 4 are skipped, they were used for one study in 2012,
        and not actively maintained anymore, but because it is legacy code not omitted.}
    }
  }
}
\value{
 The function provides no values, it only ensures that the output from other
 functions is stored in .RData(one file per accelerometer file) in folder structure
}
\examples{
  \dontrun{
    datafile = "C:/myfolder/mydata"
    outputdir = "C:/myresults"
    g.part1(datadir,outputdir)
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\references{
  \itemize{
    \item van Hees VT, Gorzelniak L, Dean Leon EC, Eder M, Pias M, et al. (2013) Separating
      Movement and Gravity Components in an Acceleration Signal and Implications for the
      Assessment of Human Daily Physical Activity. PLoS ONE 8(4): e61691.
      doi:10.1371/journal.pone.0061691
    \item van Hees VT, Fang Z, Langford J, Assah F, Mohammad A, da Silva IC, Trenell MI,
      White T, Wareham NJ, Brage S. Auto-calibration of accelerometer data for
      free-living physical activity assessment using local gravity and temperature:
      an evaluation on four continents. J Appl Physiol (1985). 2014 Aug 7
    \item Aittasalo M, Vaha-Ypya H, Vasankari T, Husu P, Jussila AM, and Sievanen H. Mean
      amplitude deviation calculated from raw acceleration data: a novel method for
      classifying the intensity of adolescents physical activity irrespective of accelerometer
      brand. BMC Sports Science, Medicine and Rehabilitation (2015).
  }
}
\name{g.impute}
\alias{g.impute}
\title{
  Function to identify invalid periods in the meta-data as generated by \link{g.getmeta}
  and to impute these invalid periods with the average of similar timepoints on other
  days of the measurement
}
\description{
  Functions takes the output from \link{g.getmeta} and information about the study
  protocol to label impute invalid time segments in the data.
}
\usage{
  g.impute(M, I, params_cleaning = c(),
  desiredtz="", dayborder= 0, TimeSegments2Zero =c(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
    output from \link{g.getmeta}
  }
  \item{I}{
    output from \link{g.inspectfile}
  }
  \item{params_cleaning}{
    See \link{g.part1}
  }
  \item{desiredtz}{
    See \link{g.part1}
  }
  \item{dayborder}{
    See \link{g.part1}
  }
  \item{TimeSegments2Zero}{
    Optional data.frame to specify which time segments need to be ignored for the imputation,
    and acceleration metrics to be imputed by zeros. The data.frame is expected
    to contain two columns named windowstart and windowend, with the start- and end
    time of the time segment in POSIXlt class.
  }
  \item{...}{
     Any argument used in the previous version of g.impute, which will now
     be used to overrule the arguments specified with the parameter objects.
  }
}
\value{
  \item{metashort}{imputed short epoch variables}
  \item{rout}{matrix to clarify when data was imputed for each long epoch time window
  and the reason for imputation. Value = 1 indicates imputation. 
  Columns 1 = monitor non wear, column 2 = clipping, column 3 = additional nonwear,
  column 4 = protocol based exclusion and column5 = sum of column 1,2,3 and 4. }
  \item{averageday}{matrix with n columns for n metrics values and m rows for
  m short epoch time windows in an average 24 hours period}
}
\examples{
  \dontrun{
    #inspect file:
    I = g.inspectfile(datafile)
    #autocalibration:
    C = g.calibrate(datafile) 
    #get meta-data:
    M = g.getmeta(datafile)
  }
  data(data.getmeta)
  data(data.inspectfile)
  #impute meta-data:
  IMP = g.impute(M=data.getmeta, I=data.inspectfile)
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.extractheadervars}
\alias{g.extractheadervars}
\title{
Extracts header variables from header object
}
\description{
Function is not intended for direct interaction by package end user
}
\usage{
g.extractheadervars(I)
}
\arguments{
  \item{I}{
  Object produced by \link{g.inspectfile}
  }
}
\value{
\itemize{
\item ID = participant identifier
\item iid = investigator identifier
\item HN = handedness
\item BodyLocation = Attachement location of the sensor
\item SX = sex
\item deviceSerialNumber = serial number
}
}
\examples{
data(data.inspectfile)
headervars = g.extractheadervars(I=data.inspectfile)
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.report.part5}
\alias{g.report.part5}
\title{
Generate report from milestone data produced by \link{g.part5}
}
\description{
Creates report from milestone data produced by \link{g.part5}. Not intended
  for direct use by package user
}
\usage{
g.report.part5(metadatadir=c(),f0=c(),f1=c(),loglocation=c(),
                          includenightcrit=c(),includedaycrit=c(),
                          data_cleaning_file=c(),
                          includedaycrit.part5=2/3,
                          minimum_MM_length.part5=23,
                          week_weekend_aggregate.part5=FALSE,
                          LUX_day_segments=c())
}
\arguments{
  \item{metadatadir}{
  see \link{g.part5}
  }
  \item{f0}{
  see \link{g.part5}
  }
  \item{f1}{
  see \link{g.part5}
  }
  \item{loglocation}{
  see \link{g.part4}
  }
  \item{includenightcrit}{
    Despricated as of version 2.0, not used anymore in part 5 report
  }
  \item{includedaycrit}{
    Despricated as of version 2.0, not used anymore in part 5 report
  }
  \item{data_cleaning_file}{
  see \link{g.part4}
  }
  \item{includedaycrit.part5}{
    Inclusion criteria for number of valid hours, either
    as expressed as a ratio of 1 or as the number of hours in a 24 hour day.
  }
  \item{minimum_MM_length.part5}{
    Minimum length in hours of a MM day to be included in the cleaned
    part 5 results.
  }
   \item{week_weekend_aggregate.part5}{
   Boolean to indicate whether week and weekend-days aggregates
   should be stored. This is turned off by default as it generates a 
   large number of extra columns in the output report.
   }
   \item{LUX_day_segments}{
   see \link{g.part5}
   }
}
\value{
Function does not produce data, but only writes reports
in csv format

The following files are stored in the root of the results folder:
part5_daysummary_*
part5_personsummary_*

The following files are stored in the folder results/QC:
part5_daysummary_full_*

See package vignette paragraph "Waking-waking or 24 hour time-use analysis"
and "Output part 5" for a more elaborative description of
the full day time-use and analysis and reporting.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{create_test_sleeplog_csv}
\alias{create_test_sleeplog_csv}
\title{
Creates csv sleeplog file for testing purposes
}
\description{
Creates sleeplog file in the format as expected by g.part4
with dummy data (23:00 onset, 07:00 waking time for every night).
}
\usage{
  create_test_sleeplog_csv(Nnights=7,storagelocation=c(), advanced=FALSE)
}
\arguments{
  \item{Nnights}{
    Number of nights (minimum is 1)
  }
  \item{storagelocation}{
    Location where the test file named testfile.csv will be stored
    If no value is provided then the function uses the current 
    working directory
  }
  \item{advanced}{
    Boolean to indicate whether to create an advanced sleeplog that also includes
    logs of nap times and nonwear
  }
}
\value{
 The function does not produce any output values. Only the file is
 stored
}
  
\examples{
  \dontrun{
    create_test_sleeplog_csv()
  }
}
\name{identify_levels}
\alias{identify_levels}
\title{
  Identifies levels of behaviour for g.part5 function.
}
\description{
  Identifies levels of behaviour from acceleratioon
  and sustained inactivity sibdetection (using angles). Function not
  intended for direct use by package user.
}
\usage{
  identify_levels(ts, TRLi,TRMi,TRVi,
                  ws3, params_phyact, ...)
}
\arguments{
  \item{ts}{
    Data.frame with time series genrated in .gpart5
  }
  \item{TRLi}{
    Numeric acceleration threshold light
  }
  \item{TRMi}{
    Numeric acceleration threshold moderate
  }
  \item{TRVi}{
    Numeric acceleration threshold vigorous
  }
  \item{ws3}{
    Numeric size of epoch in seconds
  }
  \item{params_phyact}{
    See \link{g.part2}
  }
  \item{...}{
     Any argument used in the previous version of identify_level, which will now
     be used to overrule the arguments specified with the parameter objects.
  }
}
\value{
  List with items:
  item{LEVELS}{}
  item{OLEVELS}{}
  item{Lnames}{}
  item{bc.mvpa}{}
  item{bc.lig}{}
  item{bc.in}{}
  item{ts}{}
}
\examples{
  \dontrun{
    levels = identify_levels(TRLi,TRMi,TRVi,
                               boutdur.mvpa,boutcriter.mvpa,
                               boutdur.lig,boutcriter.lig,
                               boutdur.in,boutcriter.in,
                               ws3,bout.metric)
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{datadir2fnames}
\alias{datadir2fnames}
\title{
Generates vector of file names out of datadir input
argument
}
\description{
Uses input argument datadir from \link{g.part1} and
the output from \link{isfilelist} to generate vector of filenames
}
\usage{
datadir2fnames(datadir,filelist)
}
\arguments{
  \item{datadir}{
  See \link{g.part1}
  }
   \item{filelist}{
  Produced by \link{isfilelist}
  }
}
\value{
Character vector of filenames
}

\examples{
\dontrun{
datadir2fnames(datadir = "C:/mydatafolder",filelist=TRUE)
}
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{get_starttime_weekday_meantemp_truncdata }
\alias{get_starttime_weekday_meantemp_truncdata}
\title{
  Get starttime (adjusted), weekday, mean temp, and adjust data accordingly.
}
\description{
  Function not intended for direct use by user.
  Used inside \link{g.getmeta} as an intermediate step between
  loading the raw data and calibrating it. This step includes extracting
  the starttime and adjusting it to nearest integer number of long epoch window
  lengths in an hour, truncating the data accordingly, extracting the 
  corresponding weekday and mean temperature (if temperature is available).
}
\usage{
  get_starttime_weekday_meantemp_truncdata(temp.available, monc, 
  dformat, data, selectdaysfile, P, header, desiredtz, sf, i,
  datafile, ws2, starttime, wday, weekdays, wdayname)
}
\arguments{
  \item{temp.available}{
    Boolean whether temperate is available.
  }
  \item{monc}{
    See \link{g.inspectfile}
  }
  \item{dformat}{
    See \link{g.dotorcomma}
  }
  \item{data}{
    Data part of \link{g.readaccfile} output
  }
  \item{selectdaysfile}{
    See \link{g.dotorcomma}
  }
  \item{P}{
    data loaded from accelerometer file with \link{g.readaccfile}
  }
  \item{header}{
    Header part of \link{g.readaccfile} output
  }
  \item{desiredtz}{
    See \link{g.getmeta}
  }
  \item{sf}{
    Numeric, sample frequency in Hertz
  }
  \item{i}{
    Integer index of passed on from \link{g.getmeta}
    to indicate what data block is being read.
  }
  \item{datafile}{
    See \link{g.getmeta}
  }
  \item{ws2}{
    Long epoch length
  }
  \item{starttime}{
    Once calculate it is remembered and fed into this function again,
    such that it does not have to be recalulated.
  }
  \item{wday}{
    Once calculate it is remembered and fed into this function again,
    such that it does not have to be recalulated.
  }
  \item{weekdays}{
    Once calculate it is remembered and fed into this function again,
    such that it does not have to be recalulated.
  }
  \item{wdayname}{
    Once calculate it is remembered and fed into this function again,
    such that it does not have to be recalulated.
  }
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.intensitygradient}
\alias{g.intensitygradient}
\title{
Intensity gradient calculation
}
\description{
Calculates the intensity gradient based on Rowlands et al. 2018.
The function assumes that the user has already calculated the value
distribution.
}
\usage{
g.intensitygradient(x,y)
}
\arguments{
  \item{x}{
  Numeric vector of mid-points of the bins (mg)
  }
  \item{y}{
  Numeric vector of time spent in bins (minutes)
  }
}
\value{
 \item{y_intercept}{y-intercept of a linear regression line in log-log space}
 \item{gradient}{Beta coefficient of a linear regression line in log-log space}
 \item{rsquared}{R squared of x and y values in log-log space}
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}
\references{
Rowlands A, Edwardson CL, et al. (2018) Beyond Cut Points: Accelerometer Metrics 
that Capture the Physical Activity Profile. MSSE 50(6):1. 
doi:10.1249/MSS.0000000000001561
}
\name{data.getmeta}
\alias{data.getmeta}
\docType{data}
\title{
Example output from g.getmeta
}
\description{
data.getmeta is example output from \link{g.getmeta}
}
\usage{data(data.getmeta)}
\format{
  The format is:
 chr "data.getmeta"
}
\source{
The data was collected on one individual for testing purposes
}
\examples{
data(data.getmeta)
}
\keyword{datasets}
\name{check_params}
\alias{check_params}
\title{
  Check default parameters
}
\description{
  Checks parameter objects for class and logical combinations.
  Called from \link{extract_params}. Not intended for direct use by GGIR users.
}
\usage{
  check_params(params_sleep = c(), params_metrics = c(),
                        params_rawdata = c(), params_247 = c(),
                        params_phyact = c(), params_cleaning = c(),
                         params_output = c(), params_general = c())
}
\arguments{
 \item{params_sleep}{
    List with sleep parameters
  }
  \item{params_metrics}{
    List with parameters related to metrics
  }
  \item{params_rawdata}{
    List with parameters related to raw data reading and processing
  }
  \item{params_247}{
    List with parameters related to 24/7 behavioural analysis, which includes anything
    that does not fit with physical activity or sleep research
  }
  \item{params_phyact}{
    List with parameters related to physical activity analysis
  }
  \item{params_cleaning}{
    List with parameters related to cleaning the time series, including masking and imputation
  }
  \item{params_output}{
    List with parameters related to how GGIR stores its output
  }
  \item{params_general}{
    List with parameters related to general topics
  }
}
\value{
  Lists of updated parameter objects 
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.downsample}
\alias{g.downsample}
\title{
  Downsample a vector of numeric values at three time resolutions
}
\description{
  Downsamples a vector of numeric values at three time resolutions:
  1 seconds, ws3 seconds, and ws2 second. Function is not intended
  for direct interaction by package end user
  
}
\usage{
  g.downsample(sig,fs,ws3,ws2)	
}
\arguments{
  \item{sig}{
    Vector of numeric values
  }
  \item{fs}{
    Sample frequency
  }
  \item{ws3}{
    ws3 epoch size, e.g. 5 seconds
  }
  \item{ws2}{
    ws2 epoch size, e.g. 90 seconds
  }
}
\value{
List with three object: var1, var2, and var3
corresponding to downsample time series at
1 seconds, ws2 seconds, and ws3 seconds resoluton, respectively
}
\examples{
  sig = runif(n=10000,min=1,max=10)
 downsampled_sig = g.downsample(sig,fs=20,ws3=5,ws2=15)
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.part3}
\alias{g.part3}
\title{
  Detection of sustained inactivity periods as needed for sleep detection
  in g.part4.
}
\description{
  Function called by g.shell.GGIR. It estimates the sustained inactivity
  periods in each day, which are used as input for g.part4 which then labels
  them as nocturnal sleep or day time sustained inactivity periods.
  Typical users should work with function g.shell.GGIR only.
}
\usage{
g.part3(metadatadir=c(), f0, f1, myfun=c(), 
  params_sleep = c(), params_metrics = c(), params_output = c(), params_general = c(),
  ...)
}
\arguments{
  \item{metadatadir}{
    Directory that holds a folder 'meta' and inside this a folder 'basic'
    which contains the milestone data produced by \link{g.part1}. The folderstructure
    is normally created by g.part1 and g.shell.GGIR will recognise what
    the value of metadatadir is.
  }
  \item{f0}{
    File index to start with (default = 1). Index refers to the filenames
    sorted in increasing order
  }
  \item{f1}{
    File index to finish with (defaults to number of files available)
  }
  \item{myfun}{
    External function object to be applied to raw data.
    See details \link{applyExtFunction}.
  }
  \item{params_sleep}{
    See details
  }
  \item{params_metrics}{
    List, see \link{g.part1}.
  }
  \item{params_output}{
    List, see \link{g.part1}
  }
  \item{params_general}{
    List, see \link{g.part1}
  }
  \item{...}{
    To enable compatibility with R scripts written for older GGIR versions,
    the user can also provide parameters listed in the params_ objects as direct argument.
  }
}
\details{
  GGIR comes with many processing parameters, which have been thematically grouped in
  parameter objects (R list). By running print(load_params()) you can
  see the default values of all the parameter objects. When g.part 3 is used via \link{g.shell.GGIR}
  you have the option to specifiy a configuration file, which will overrule the default
  parameter values. Further, as user you can set parameter values as input argument to both g.part3
  and \link{g.shell.GGIR}. Directly specified argument overrule the configuration file and default values.
  
  The parameter objects used by GGIR part 3 (g.part3) that are no already discussed in
  \link{g.part1} or \link{g.part2} are:
  
  \subsection{params_sleep}{
    A list of parameters used to configure the sleep analysis as performend in
    GGIR part 3 and 4.
    \describe{
      \item{relyonguider}{Boolean. If TRUE then sleep onset and waking time are defined based on
        timestamps derived from the guider. If participants were instructed NOT to wear the accelerometer
        during waking hours then set to TRUE, in all other scenarios
        set to FALSE (default).}
      \item{relyonsleeplog}{Do not use, now replaced by argument relyonguider.
        Values provided to argument relyonsleeplog will be passed on to 
        argument relyonguider to not preserve functionality of old R scripts.}
      \item{def.noc.sleep}{Numeric. The time window during which sustained
        inactivity will be assumed to represent sleep, e.g. def.noc.sleep=c(21,9).
        This is only used if no sleep log entry is available. If def.noc.sleep is
        left blank 'def.noc.sleep=c()' then the 12 hour window centred
        at the least active 5 hours of the 24 hour period will be used
        instead. Here, L5 is hardcoded and will not change by changing
        argument winhr in function \link{g.part2}. If def.noc.sleep is filled
        with a single integer, e.g. def.noc.sleep=c(1) then the window
        will be detected with based on built in algorithms.
        See argument HASPT.algo from \link{HASPT} for specifying which of those
        algorithms to use.}
      \item{sleepwindowType}{Character to indicate type of sleeplog, default "SPT".
        Set to "TimeInBed" if sleep log recorded time in bed to enable calculation
        of sleep latency and sleep efficiency.}
      \item{nnights}{Number of nights for which sleep log information should be available. 
        It assumes that this is constant within a study. If sleep log information
        is missing for certain nights then leave these blank.}
      \item{loglocation}{Character. Path to csv file with sleep log
        information. See package vignette for how to format this file.}
      \item{colid}{Numeric. Column number in the sleep log spreadsheet in which
          the participant ID code is stored (default = 1)}
      \item{coln1}{Numeric. Column number in the sleep log spreadsheet where 
          the onset of the first night starts}
      \item{sleeplogidnum}{Boolean. Should the participant identifier as stored in
          the sleeplog be interpretted as a number (TRUE=default) or character (FALSE)?}
      \item{ignorenonwear}{Boolean, see \link{g.sib.sum}}
      \item{constrain2range}{Boolean,  Whether or not to constrain the range of
        threshold used in the diary free sleep period time window detection}
      \item{HASPT.ignore.invalid}{Boolean, see \link{HASPT}}
      \item{HASPT.algo}{Character, character to indicate what heuristic algorithm
        to use for detecting the SPT window, see \link{HASPT}}
      \item{HASIB.algo}{Character, character to indicate what heuristic algorithm to use 
        for detecting the sustained inactivity bouts (SIB), see \link{HASIB}}
      \item{Sadeh_axis}{Character, character to indicate which axis to use for
        the Sadeh1994 algorithm, and  other algortihms that related on count-based
        Actigraphy such as Galland2012.}
      \item{sleeplogsep}{Character, \link{g.loadlog}}
      \item{nap_model}{Character, see \link{g.part5.classifyNaps}}
      \item{longitudinal_axis}{Integer to indicate which axis is the longitudinal axis. 
        If not provided function will estimate longitudinal axis. Only used when
        sensor.location="hip" or HASPT.algo="HorAngle".} 
      \item{anglethreshold}{Numeric, Angle threshold (degrees) for sustained 
        inactivity periods detection, default = 5}
      \item{timethreshold}{Numeric, time threshold (minutes) for sustained
        inactivity periods detection, default = 5. This can be specified
        as multiple thresholds, each of which will be implemented. 
        For example, timethreshold = c(5,10)}
      \item{possible_nap_window}{Numeric, see \link{g.part5.classifyNaps}}
      \item{possible_nap_dur}{Numeric, see \link{g.part5.classifyNaps}}
    }
  }
}
\value{
  The function provides no values, it only ensures that other functions
  are called and that their output is stored in .RData files.
  \cr
  \itemize{
    \item \code{night} nightnumber
    \item \code{definition} definition of sustained inactivity. For example,
    T10A5 refers to 10 minute window and a 5 degree angle (see paper for
    further explaination).
    \item \code{start.time.day} timestamp when the day started
    \item \code{nsib.periods} number of sustained inactivity bouts
    \item \code{tot.sib.dur.hrs} total duration of all sustained inactivity bouts
    \item \code{fraction.night.invalid} fraction of the night for which
    accelerometer data was invalid, e.g. monitor not worn
    \item \code{sib.period} number of sustained inactivity period
    \item \code{sib.onset.time} onset time of sustained inactivity period
    \item \code{sib.end.time} end time of sustained inactivity period
  }
}
\examples{
  \dontrun{
    metadatadir = "C:/myfolder/meta" # assumes that there is a subfolder in
    # metadatadir named 'basic' containing the output from g.part1
    g.part3(metadatadir=metadatadir, anglethreshold=5,
    timethreshold=5, overwrite=FALSE)
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\references{
  \itemize{
    \item van Hees VT, Sabia S, et al. (2015) A novel, open access method to
    assess sleep duration using a wrist-worn accelerometer, PLoS ONE, November 2015
    \item van Hees VT, Sabia S, et al. (2018) Estimating sleep parameters
    using an accelerometer without sleep diary. Scientific Reports.
  }
}
\name{g.part5.lux_persegment}
\alias{g.part5.lux_persegment}
\title{
  Extract key lux variables per segment of the data.
}
\description{
  Extracts per segment of the day: mean lux, time above 1000 lux, time awake,
  and time LUX imputed. Function not intended
  for direct use by package user.
}
\usage{
  g.part5.lux_persegment(ts, sse, LUX_day_segments, ws3new)
}
\arguments{
  \item{ts}{
  }
  \item{sse}{
  }
  \item{LUX_day_segments}{
  }
  \item{ws3new}{
  }
}
\value{
 List with values (vector) of the derived variables and corresponding names (vector).
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{is.ISO8601}
\alias{is.ISO8601}
\title{
Check whether character timestamp is in iso8601 format.
}
\description{
Checks whether timestamp stored in character format is in ISO8601 format or not}
\usage{
is.ISO8601(x)	
}
\arguments{
  \item{x}{
Timestamps in character format either in ISO8601 or as "yyyy-mm-dd hh:mm:ss".
}
}
\examples{
x ="1980-1-1 18:00:00"
is.ISO8601(x)
}
\name{extract_params}
\alias{extract_params}
\title{
  Extract parameters from input and add them to params
}
\description{
  Extracts parameters separately provided by input and adds them to the params objects.
  Not intended for direct use by GGIR users.
}
\usage{
  extract_params(params_sleep = c(), params_metrics = c(),
                 params_rawdata = c(), params_247 = c(),
                 params_phyact = c(), params_cleaning = c(),
                 params_output = c(), params_general = c(), input = c(),
                 configfile_csv = c())
}
\arguments{
  \item{params_sleep}{
    List with sleep parameters
  }
  \item{params_metrics}{
    List with parameters related to metrics
  }
  \item{params_rawdata}{
    List with parameters related to raw data reading and processing
  }
  \item{params_247}{
    List with parameters related to 24/7 behavioural analysis, which includes anything
    that does not fit with physical activity or sleep research
  }
  \item{params_phyact}{
    List with parameters related to physical activity analysis
  }
  \item{params_cleaning}{
    List with parameters related to cleaning the time series, including masking and imputation
  }
  \item{params_output}{
    List with parameters related to how GGIR stores its output
  }
  \item{params_general}{
    List with parameters related to general topics
  }
  \item{input}{
    All objects provided by users
  }
  \item{configfile_csv}{
    Csv configuration file
  }
}
\value{
  Lists of updated parameter objects 
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{iso8601chartime2POSIX}
\alias{iso8601chartime2POSIX}
\title{
Convert iso8601 timestamps to POSIX timestamp
}
\description{
To avoid ambiguities when sharing and comparing timestamps. All timestamps
are expressed in iso8601 format: https://en.wikipedia.org/wiki/ISO_8601
However, to generate plots in R we need to convert them back to POSIX
}
\usage{
iso8601chartime2POSIX(x,tz)
}
\arguments{
  \item{x}{
Vector of timestamps in iso8601 in character format
}
  \item{tz}{
 Timezone of data collection, e.g. "Europe/London". 
 See List_of_tz_database_time_zones on Wikipedia
 for full list.
}
}
\examples{
x ="2017-05-07T13:00:00+0200"
tz = "Europe/Amsterdam"
x_converted = iso8601chartime2POSIX(x,tz)
}\name{GGIR-package}
\alias{GGIR-package}
\alias{GGIR}
\docType{package}
\title{
  A package to process multi-day raw accelerometer data
}
\description{
  Disclaimer: If you are a new GGIR user then please see
  \href{https://cran.r-project.org/package=GGIR/vignettes/GGIR.html}{package vignette}
   for an introduction to GGIR.\cr
  \cr
  This document is primarily aimed at documenting the functions and their input arguments.\cr
  \cr
  Please note that there is google discussion group for this package (link below).\cr
  \cr
  You can thank us for sharing the code in this package and for developing
  it as a generic purpose tool by citing the package name and by
  citing the supporting publications (e.g. Migueles et al. 2019) in your publications.
}
\details{
  \tabular{ll}{
  Package: \tab GGIR\cr
  Type: \tab Package\cr
  Version: \tab 2.5-6\cr
  Date: \tab 2022-01-18\cr
  License: \tab LGPL (>= 2.0, < 3)\cr
  Discussion group: \tab https://groups.google.com/forum/#!forum/rpackageggir\cr
  }
}
\examples{
  \dontrun{
    #inspect file:
    I = g.inspectfile(datafile)

    #autocalibration:
    C = g.calibrate(datafile)

    #get meta-data:
    M = g.getmeta(datafile)
  }
  data(data.getmeta)
  data(data.inspectfile)
  data(data.calibrate)

  #impute meta-data:
  IMP = g.impute(M = data.getmeta, I = data.inspectfile)
  #analyse and produce summary:
  A = g.analyse(I = data.inspectfile, C = data.calibrate, M = data.getmeta, IMP)
  #plot data
  g.plot(IMP, M = data.getmeta, I = data.inspectfile, durplot=4)
}
\author{
  \itemize{
    \item Vincent T van Hees <v.vanhees@accelting.com> main creator and developer
    \item Zhou Fang developed calibration algorithm used in function \link{g.calibrate}
    \item Jing Hua Zhao <jinghua.zhao@mrc-epid.cam.ac.uk> co-developed function \link{g.binread}
    \item Joe Heywood helped develop the functionality to process specific recording days
    \item Evgeny Mirkes created function \link{g.cwaread}
    \item Severine Sabia, Mathilde Chen, and Manasa Yerramalla extensively tested and provided feedback on various functions
    \item Joan Capdevila Pujol helped to improve various functions
    \item Jairo H Migueles <jairohm@ugr.es> helped to improve various functions
    \item Dan Jackson helped with unpack function for AX3 data.
    \item Matthew R Patterson helped with enhancing the visual report.
    \item Lena Kushleyeva helped fix bug in sleep detection.
  }
}
\references{
  \itemize{
    \item Migueles JH, Rowlands AV, et al. GGIR: A Research Community-Driven Open Source
    R Package for Generating Physical Activity and Sleep Outcomes From Multi-Day Raw
    Accelerometer Data. Journal for the Measurement of Physical Behaviour. 2(3) 2019.
    doi:10.1123/jmpb.2018-0063.
    \item van Hees VT, Gorzelniak L, Dean Leon EC, Eder M, Pias M, et al. (2013) Separating
    Movement and Gravity Components in an Acceleration Signal and Implications for the
    Assessment of Human Daily Physical Activity. PLoS ONE 8(4): e61691.
    doi:10.1371/journal.pone.0061691
    \item van Hees VT, Fang Z, Langford J, Assah F, Mohammad A, da Silva IC, Trenell MI,
    White T, Wareham NJ, Brage S. Auto-calibration of accelerometer data for
    free-living physical activity assessment using local gravity and temperature:
    an evaluation on four continents. J Appl Physiol (1985). 2014 Aug 7
    \item van Hees VT, Sabia S, et al. (2015) A novel, open access method to
    assess sleep duration using a wrist-worn accelerometer, PLoS ONE, November 2015
  }
}
\name{getFirstTimestamp}
\alias{getFirstTimestamp}
\title{
  Extract first timestamp from GENEActiv file
}
\description{
  Extract first timestamp from GENEActiv file, only used
  when using the selectdaysfile argument. Function not designed for
  direct use by package user.
}
\usage{
getFirstTimestamp(f, p1)
}
\arguments{
  \item{f}{
    GENEActiv filename
  }
  \item{p1}{
    First value of timestamps object
  }
}
\value{
  POSIX object withstarttime
}
\author{
  Joe Heywood <j.heywood@ucl.ac.uk>
}\name{g.part5.addsib}
\alias{g.part5.addsib}
\title{
  Adds the sustained inactivity bout to the ts series.
}
\description{
  Not intended for direct use by GGIR users.
  Adds the sustained inactivity bout to the ts series
  as part of \link{g.part5}.
}
\usage{
  g.part5.addsib(ts,ws3, Nts, S2, desiredtz, j, nightsi)
}
\arguments{
  \item{ts}{
  }
  \item{ws3}{
  }
  \item{Nts}{
  }
  \item{S2}{
  }
  \item{desiredtz}{
  }
  \item{j}{
  }
  \item{nightsi}{
  }
}
\value{
  Data.frame ts
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.fragmentation}
\alias{g.fragmentation}
\title{
Fragmentation metrics from time series.
}
\description{
The function is used by \link{g.part5} to derive time series 
fragmentation metrics. The function assumes that NA values and nonwear time is 
accounted for before the data enters the function.
}
\usage{
  g.fragmentation(frag.metrics = c("mean", "TP", "Gini", "power",
        "CoV", "NFragPM", "all"), LEVELS = c(), Lnames=c(), xmin=1)
}
\arguments{
   \item{frag.metrics}{
   Character with fragmentation metric to exract. Can be "mean", "TP", "Gini", 
   "power", or "CoV", "NFragPM", or all the above metrics with "all". See details.
  }
  \item{LEVELS}{
    Numeric vector of behavioural level classes derived with \link{identify_levels}
  }
  \item{Lnames}{
    Character vector with names of classes used in LEVELS, see details.
  }
  \item{xmin}{
    Numeric scalar to indicate the minimum recordable fragment length. In \link{g.part5}
    this is derived from the epoch length.
  }
}
\value{
List with Character object showing how decimals are separated
\item{TP_PA2IN}{Transition probability physical activity to inactivity}
\item{TP_IN2PA}{Transition probability physical inactivity to activity}
\item{Nfrag_IN2LIPA}{Number of inacitivty fragments succeeded by LIPA
(light physical activity)}
\item{TP_IN2LIPA}{Transition probability physical inactivity to LIPA}
\item{Nfrag_IN2MVPA}{Number of inacitivty fragments succeeded by MVPA
(moderate or vigorous physical activity)}
\item{TP_IN2MVPA}{Transition probability physical inactivity to MVPA}
\item{Nfrag_MVPA}{Number of MVPA fragments}
\item{Nfrag_LIPA}{Number of LIPA fragments}
\item{mean_dur_MVPA}{mean MVPA fragment duration}
\item{mean_dur_LIPA}{mean LIPA fragment duration}
\item{Nfrag_IN}{Number of inactivity fragments}
\item{Nfrag_PA}{Number of activity fragments}
\item{mean_dur_IN}{mean duration inactivity fragments}
\item{mean_dur_PA}{mean duration activity fragments}
\item{Gini_dur_IN}{Gini index corresponding to inactivity fragment durations}
\item{Gini_dur_PA}{Gini index corresponding to activity fragment durations}
\item{CoV_dur_IN}{Coefficient of Variance corresponding to inactivity fragment durations}
\item{CoV_dur_PA}{Coefficient of Variance corresponding to activity fragment durations}
\item{alpha_dur_IN}{Alpha of the fitted power distribution through inactivity fragment durations}
\item{alpha_dur_PA}{Alpha of the fitted power distribution through activity fragment durations}
\item{x0.5_dur_IN}{x0.5 corresponding to alpha_dur_IN}
\item{x0.5_dur_PA}{x0.5 corresponding to alpha_dur_PA}
\item{W0.5_dur_IN}{W0.5 corresponding to alpha_dur_IN}
\item{W0.5_dur_PA}{W0.5 corresponding to alpha_dur_PA}
\item{NFragPM_IN}{Number of IN fragments per minutes in IN}
\item{NFragPM_PA}{Number of PA fragments per minutes in PA}
\item{SD_dur_IN}{Standard deviation in the duration of inactivity fragments}
\item{SD_dur_PA}{Standard deviation in the duration of physical activity fragments}

}
\details{
See package vignette for description of fragmentation metrics.
In short, abbreviation "TP" refers to transition probality metrics,
abbreviation "CoV" refers to Coefficient of Variance, and
metric "NFragPM" refers to the Number of fragments per minute.

Regarding the Lnames argument. The class names included in this are categorised
as follows:
\itemize{
  \item{Inactive - if name includes the character strings "day_IN_unbt" or "day_IN_bts".}
  \item{LIPA - If name includes the character strings "day_LIG_unbt" or "day_LIG_bts".}
  \item{MVPA - If name includes the character strings "day_MOD_unbt", "day_VIG_unbt", or "day_MVPA_bts"}
}
}
\examples{
\dontrun{
    x = c(6, 5, 6, 7, 6, 6, 7, 6, 6, 5, 6, 6, 6, 5, 7, 6, 6, 5, 5, 5, 6, 7, 6,
        6, 6, 6, 7, 6, 5, 5, 5, 5, 5, 6, 6, 6, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        7, 7, 6, 5, 6, 5, 6, 5, rep(12, 11), 5, 6, 6, 6, 5, 6, rep(9, 14), 6,
        5, 7, 7, 6, 7, 7, 7, 6, 6, 6, 5, 6, 5, 5, 5, 6, 5, 5, 5, 5, 5, 5, 5)
  Lnames = c("spt_sleep", "spt_wake_IN", "spt_wake_LIG", "spt_wake_MOD",
            "spt_wake_VIG", "day_IN_unbt", "day_LIG_unbt", "day_MOD_unbt",
            "day_VIG_unbt", "day_MVPA_bts_10", "day_IN_bts_30",
             "day_IN_bts_10_30", "day_LIG_bts_10")
  out = g.fragmentation(frag.metrics = "all",
                        LEVELS = x,
                        Lnames=Lnames)}
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{ismovisens}
\alias{ismovisens}
\title{
Checks whether the files to process are collected with movisens accelerometers.
}
\description{
Checks whether the files in the datadir folder are files collected
with movisens accelerometers. Note that movisens data are stored in one folder per recording that includes multiple bin-files
(instead of one file per recording as usual in other accelerometer brands). Therefore, datadir indicates
the directory where all the recording folders are stored, then, GGIR reads the pertinent bin files from
every folder.}
\usage{
ismovisens(data)	
}
\arguments{
  \item{data}{
Full path to the recording folder (with the bin files inside) or the datadir (where all the recording folders are stored).
}
}
\value{
Boolean whether it is a movisens file (TRUE) or not (FALSE)
}
  
\examples{
\dontrun{
is.mv = ismovisens(data)
}
}
\name{g.createcoordinates}
\alias{g.createcoordinates}
\title{
  Create coordinates for \link{g.plot}
}
\description{
  Function creates the coordinates for the blocks \link{g.plot}
  Function not designed for direct use by package user.
}
\usage{
g.createcoordinates(r,timeline)
}
\arguments{
  \item{r}{
    Vector of zeros and ones reflecting the moments in time when there
    should be a block (1)
  }
  \item{timeline}{
    Vector of time indicators, this can be numbers or actual timestamps
    The length of timeline needs to match the length of argument r
  }
}
\value{
  List with two objects: x0 with all the coordinates correspoding to the
  start of each blocks on the timelines and x1 with all the coordinates
  corresponding to the end of each block on the timeline
}
\author{
  Vincent van Hees <v.vanhees@accelting.com>
}\name{g.readaccfile}
\alias{g.readaccfile}
\title{
  Generic functiont to read large blocks of accelerometer data
}
\description{
  The function is used by \link{g.getmeta} and \link{g.calibrate} to
  read large blocks of the accelerometer file, which are processed and
  then deleted from memory. This is needed for memory management.
}
\usage{
  g.readaccfile(filename, blocksize, blocknumber, selectdaysfile = c(), filequality,
                           decn, ws, PreviousEndPage = 1, inspectfileobject = c(),
                           params_rawdata = c(), params_general = c(), ...)
}
\arguments{
  \item{filename}{
    filename
  }
  \item{blocksize}{
    Size of blocks (in file pages) to be read
  }
  \item{blocknumber}{
    Block number relative to start of file
  }
  \item{selectdaysfile}{
    See documentation \link{g.getmeta}
  }
  \item{filequality}{
    Single row dataframe with columns: filetooshort, filecorrupt,
    and filedoesnotholdday. All with the value TRUE or FALSE
  }
  \item{decn}{
    Character with a dot or a comma, used for interpretting
    samplefrequency in the file header. decn is derived with
    \link{g.dotorcomma}
  }
  \item{ws}{
    Larger windowsize for non-detection, see documentation \link{g.part2}
  }
  \item{PreviousEndPage}{
    Page number on which previous block ended (automatically assigned within
    g.getmeta and g.calibrate).
  }
  \item{inspectfileobject}{
    Output from the function \link{g.inspectfile}.
  }
  \item{params_rawdata}{
    See \link{g.part1}
  }
  \item{params_general}{
    See \link{g.part1}
  }
  \item{...}{
    Any input arguments needed for function \link{read.myacc.csv} if you
    are working with a non-standard csv formatted files. Furter,
    any argument used in the previous version of g.readaccfile, which will now
    be used to overrule the arguments specified with the parameter objects.
  }
}
\value{
  \itemize{
    \item \code{P} Block object extracted from file with format specific to
    accelerometer brand
    \item \code{filequality} Same as in function arguments
    \item \code{switchoffLD} Boolean to indicate whether it is worth
    continueing to read the next block of data or not
    \item \code{endpage} Page number on which blocked ends, this will be
    used as input for argument PreviousEndPage when reading the next block.
  }
}
\examples{
  \dontrun{
    filequality = data.frame(filetooshort = FALSE, filecorrupt = FALSE,
    filedoesnotholdday = FALSE)
    output = g.readaccfile(filename = "C:/myfile.bin", 
    blocksize = 20000, blocknumber = 1,
    selectdaysfile = c(), filequality = filequality,
    decn = ".", dayborder = 0, PreviousEndPage = c()) 
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.sib.det}
\alias{g.sib.det}
\title{
  sustiained inactivty bouts detection
}
\description{
  Detects sustiained inactivty bouts. Function not intended
  for direct use by package user
}
\usage{
  g.sib.det(M, IMP, I, twd = c(-12, 12),
             acc.metric = "ENMO", desiredtz = "",
             myfun=c(), sensor.location = "wrist", params_sleep = c(), ...)
}
\arguments{
  \item{M}{
    Object produced by \link{g.getmeta}
  }
  \item{IMP}{
    Object produced by \link{g.impute}
  }
  \item{I}{
    Object produced by \link{g.inspectfile}
  }
  \item{twd}{
    Vector of length 2, indicating the time window to consider
    as hours relative to midnight.
  }
    \item{acc.metric}{
    Which one of the metrics do you want to consider to analyze L5. 
    The metric of interest need to be calculated in
    M (see \link{g.part1})
  }
  \item{desiredtz}{
    See \link{g.part3}
  }
  \item{myfun}{
    External function object to be applied to raw data.
    See details \link{applyExtFunction}.
  }
  \item{sensor.location}{
    Character to indicate sensor location, default is wrist.
    If it is hip HDCZA algorithm also requires longitudinal axis of sensor to be
    between -45 and +45 degrees.
  }
   \item{params_sleep}{
    See \link{g.part3}
  }
  \item{...}{
     Any argument used in the previous version of g.sib.det, which will now
     be used to overrule the arguments specified with the parameter objects.
  }
}
\value{
  \itemize{
    \item output = Dataframe for every epoch a classification
    \item detection.failed = Boolean whether detection failed
    \item L5list = L5 for every day (defined from noon to noon)
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{chartime2iso8601}
\alias{chartime2iso8601}
\title{
Convert character timestamps to iso8601 timestamp
}
\description{
To avoid ambiguities when sharing and comparing timestamps. All timestamps
are expressed in iso8601 format: https://en.wikipedia.org/wiki/ISO_8601
}
\usage{
chartime2iso8601(x,tz)	
}
\arguments{
  \item{x}{
Vector of timestamps in character format:
year-month-date and optional followed by hour:minute:second 
For example, "1980-01-01 18:00:00"
}
  \item{tz}{
 Timezone of data collection, e.g. "Europe/London".
 See https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
 for full list
}
}
\examples{
x ="1980-1-1 18:00:00"
tz = "Europe/Amsterdam"
x_converted = chartime2iso8601(x,tz)
}
\name{HASIB}
\alias{HASIB}
\title{
  Heuristic algorithms for sustiained inactivty bouts detection
}
\description{
  Apply heuristic algorithms for sustiained inactivty bouts detection.
  Function not intended for direct use by package user
}
\usage{
  HASIB(HASIB.algo = "vanHees2015", timethreshold=c(), anglethreshold=c(), 
                 time=c(), anglez=c(), ws3=c(), zeroCrossingCount=c(), BrondCount = c())
}
\arguments{
  \item{HASIB.algo}{
    Character to indicator which sib algorithm should be used.
    Default value: "vanHees2015". Other options: "Sadeh1994", "Galland2012"
  }
  \item{anglethreshold}{
    See \link{g.sib.det}
  }
  \item{timethreshold}{
    See \link{g.sib.det}
  }
  \item{time}{
    Vector with time per short epoch
  }
  \item{anglez}{
    Vector with z-angle per short epoch
  }
  \item{ws3}{
    See \link{g.getmeta}
  }
  \item{zeroCrossingCount}{
    Vector with zero crossing counts per epoch as required for Sadeh algortihm
  }
  \item{BrondCount}{
    Vector with Brond counts per epoch to be used by the Sadeh algortihm
  }
}
\value{
  Vector with binary indicator of sustained inactivity bout, 1 is yes, 0 is no.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{g.shell.GGIR}
\alias{g.shell.GGIR}
\title{
Shell function for analysing an accelerometer dataset.
}
\description{
This function is designed to help users operate all steps of the
analysis. It helps to generate and structure milestone data, 
produces user-friendly reports. The function acts as a shell with
calls to \link{g.part1}, \link{g.part2}, \link{g.part3} and
\link{g.part4}. Please see these specific functions for clarification
on optional input arguments.
}
\usage{
g.shell.GGIR(mode=1:5,datadir=c(),outputdir=c(),studyname=c(),
f0=1,f1=0, do.report=c(2), configfile=c(),myfun=c(), ...)
}
\arguments{
  \item{mode}{
    Specify which of the four parts need to be run, e.g. mode = 1 makes
    that \link{g.part1} is run. Default setting, mode = c(1,2), makes
    that both part1 and part2 are ran. Note that if mode = c(1,3) then
    the code will also set do.anglez = TRUE in order to enable sleep
    detection. If you run part 1 and 3 seperatedly then you need to 
    remember to set argument do.anglez to TRUE when running part1.
  }
  \item{datadir}{
    Directory where the accelerometer files are stored or list, e.g. 
    "C:/mydata" of accelerometer filenames and directories, e.g. 
    c("C:/mydata/myfile1.bin", "C:/mydata/myfile2.bin").
  }
  \item{outputdir}{
    Directory where the output needs to be stored. Note that this 
    function will attempt to create folders in this directory and uses
    those folder to keep output
  }
  \item{studyname}{
    If the datadir is a folder then the study will be given the name of the
    data directory. If datadir is a list of filenames then the studyname as specified
    by this input argument will be used as name for the study
  }  
  \item{f0}{
    File index to start with (default = 1). Index refers to the filenames sorted
    in increasing order
  }
  \item{f1}{
    File index to finish with (defaults to number of files available)
  }
  \item{do.report}{
    For which parts to generate a summary spreadsheet: 2 and/or 4. Default is c(2).
    A report will be generated based on the available milestone data. When creating 
    milestone data with multiple machines it is advisable to turn the report
    generation off when generating the milestone data, value = c(),
    and then to merge the milestone data and turn report generation back
    on while setting overwrite to FALSE.  
  }
  \item{configfile}{
    Configuration file previously generated by g.shell.GGIR. See details.
  }
  \item{myfun}{
    External function object to be applied to raw data, see \link{g.getmeta}.
  }
  \item{...}{
    Any of the parameters used GGIR. Given the large number of parameters used in GGIR
    we have grouped them in objects that start with 'params_' these are documented in the 
    function documentation for \link{g.part1}, \link{g.part2}, and \link{g.part3}. You cannot provide 
    these objects as argument to g.shell.GGIR, but you can provide the parameters inside them as input
    to g.shell.GGIR.
  }  
}
\value{
  The function provides no values, it only ensures that other functions are called
  and that their output is stored. Further, a configuration file is stored containing
  all the argument values used to facilitate reproducibility.
}
\details{
  Once you have used g.shell.GGIR and the output directory (outputdir) will be filled
  with milestone data and results. Function g.shell.GGIR stores all the explicitely 
  entered argument values and default values for the argument that are not explicitely
  provided in a csv-file named config.csv stored in the root of the output folder. 
  The config.csv file is accepted as input to g.shell.GGIR with argument `configfile`
  to replace the specification of all the arguments, except `datadir` and `outputdir`.
  
  The practical value of this is that it eases the replication of analysis, because
  instead of having to share you R script, sharing your config.csv file will be 
  sufficient. Further, the config.csv file contribute to the reproducibility
  of your data analysis.
  
  Note: When combining a configuration file with explicitely provided argument
  values, the explicitely provided argument values will overrule
  the argument values in the configuration file. If a parameter is neither provided
  via the configuration file nor as input then GGIR uses its default paramter values which
  can be inspected with command \code{print(load_params())}, and if you are specifically 
  interested in a certain subgroup of parameters, let's say physical activity, then you
  can do \code{print(load_params()$params_phyact)}. These defaults are part of the GGIR
  code and cannot be changed by the user.
}
\examples{
\dontrun{
  mode= c(1,2,3,4)
  datadir= "C:/myfolder/mydata"
  outputdir= "C:/myresults"
  studyname="test"
  f0 = 1 
  f1 = 2
  g.shell.GGIR(#-------------------------------
               # General parameters
               #-------------------------------
               mode=mode, 
               datadir=datadir, 
               outputdir=outputdir, 
               studyname=studyname, 
               f0=f0, 
               f1=f1,
               overwrite = FALSE, 
               do.imp=TRUE,
               idloc=1, 
               print.filename=FALSE,
               storefolderstructure = FALSE,
               #-------------------------------
               # Part 1 parameters:
               #-------------------------------
               windowsizes = c(5,900,3600),
               do.cal=TRUE, 
               do.enmo = TRUE,
               do.anglez=TRUE,
               chunksize=1,
               printsummary=TRUE,
               #-------------------------------
               # Part 2 parameters:
               #-------------------------------
               strategy = 1,
               ndayswindow=7,
               hrs.del.start = 1,
               hrs.del.end = 1, 
               maxdur = 9,
               includedaycrit = 16, 
               L5M5window = c(0,24),
               M5L5res = 10,
               winhr = c(5,10),
               qlevels = c(c(1380/1440),c(1410/1440)),
               qwindow=c(0,24), 
               ilevels = c(seq(0,400,by=50),8000), 
               mvpathreshold =c(100,120),
               #-------------------------------
               # Part 3 parameters:
               #-------------------------------
               timethreshold= c(5,10),
               anglethreshold=5,
               ignorenonwear = TRUE,
               #-------------------------------
               # Part 4 parameters:
               #-------------------------------
               excludefirstlast = FALSE,
               includenightcrit = 16,
               def.noc.sleep = 1,
               loglocation= "D:/sleeplog.csv",
               outliers.only = FALSE,
               criterror = 4,
               relyonsleeplog = FALSE,
               sleeplogidnum = TRUE,
               colid=1, 
               coln1=2, 
               do.visual = TRUE,
               nnights = 9,
               #-------------------------------
               # Part 5 parameters:
               #-------------------------------
               # Key functions: Merging physical activity with sleep analyses
               threshold.lig = c(30,40,50),
               threshold.mod = c(100,120),
               threshold.vig = c(400,500),
               excludefirstlast = FALSE,
               boutcriter = 0.8,
               boutcriter.in = 0.9,
               boutcriter.lig = 0.8,
               boutcriter.mvpa = 0.8,
               boutdur.in = c(10,20,30),
               boutdur.lig = c(1,5,10),
               boutdur.mvpa = c(1,5,10),
               timewindow = c("WW"),
               #-----------------------------------
               # Report generation
               #-------------------------------
               do.report=c(2,4,5))
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\references{
  \itemize{
    \item van Hees VT, Gorzelniak L, Dean Leon EC, Eder M, Pias M, et al. (2013) Separating
    Movement and Gravity Components in an Acceleration Signal and Implications for the
    Assessment of Human Daily Physical Activity. PLoS ONE 8(4): e61691.
    doi:10.1371/journal.pone.0061691
    \item van Hees VT, Fang Z, Langford J, Assah F, Mohammad A, da Silva IC, Trenell MI, 
    White T, Wareham NJ, Brage S. Auto-calibration of accelerometer data for
    free-living physical activity assessment using local gravity and temperature: 
    an evaluation on four continents. J Appl Physiol (1985). 2014 Aug 7
    \item van Hees VT, Sabia S, et al. (2015) A novel, open access method to
    assess sleep duration using a wrist-worn accelerometer, PLoS ONE, November 2015
  }
}
\name{g.part4}
\alias{g.part4}
\title{
  Labels detected sustained inactivity periods by g.part3 as either
  part of the Sleep Period Time window or not
}
\description{
  Combines output from g.part3 and guider information to estimate
  sleep variables. See vignette paragraph "Sleep and full day
  time-use analysis in GGIR" for an elaborate descript of the sleep detection.
}
\usage{
 g.part4(datadir = c(), metadatadir = c(), f0 = f0, f1 = f1, params_sleep = c(), 
    params_metrics = c(),  params_cleaning = c(), params_output = c(),
                   params_general = c(), ...)
}
\arguments{
  \item{datadir}{
    Directory where the accelerometer files are stored or list
    of accelerometer filenames and directories
  }
  \item{metadatadir}{
    Directory that holds a folders 'meta' and inside this a folder
    'basic' which contains the milestone data produced by g.part1.
    The folderstructure is normally created by g.part3. When using
    g.part4 via g.shell.GGIR then g.shell.GGIR will automatically
    recognise what the value of metadatadir is, so the user does not
    need to specify this.
  }
  \item{f0}{
    File index to start with (default = 1). Index refers to the
    filenames sorted in increasing order
  }
  \item{f1}{
    File index to finish with (defaults to number of files available)
  }
  \item{params_sleep}{
    List of parameters used for sleep analysis (GGIR part 3, 4, and 5): 
    see documentation \link{g.part3}.
  }
  \item{params_metrics}{
    List of parameters used for metrics extraction (GGIR part 1): 
    see documentation \link{g.part1}.
  }
  \item{params_cleaning}{
    List, see \link{g.part1}
  }
  \item{params_output}{
    List, see \link{g.part2}
  }
  \item{params_general}{
    List, see \link{g.part1}
  }
  \item{...}{
    To enable compatibility with R scripts written for older GGIR versions,
    the user can also provide parameters listed in the params_ objects as direct argument.
  }
}
\value{
  The function does not produce values but generates an RData file
  in the milestone subfolder ms4.out which incudes a dataframe
  named \code{nightsummary}. This dataframe is used in g.report.part4 to create
  two reports one per night and one per person. See package vignette 
  paragraph "Output part 4" for description of all the variables.
}

\examples{
  \dontrun{
    metadatadir = "C:/myfolder/meta" # assumes that there is a subfolder in
    # metadatadir named 'ms3.out' containing the output from g.part3
    g.part4(metadatadir=metadatadir)
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\references{
  \itemize{
    \item van Hees VT, Sabia S, et al. (2018) AEstimating sleep parameters
    using an accelerometer without sleep diary, Scientific Reports.
    \item van Hees VT, Sabia S, et al. (2015) A novel, open access method
    to assess sleep duration using a wrist-worn accelerometer, PLoS ONE.
  }
}
\name{createConfigFile}
\alias{createConfigFile}
\title{
Creates Config File based on variables in g.shell.GGIR environment
}
\description{
Only used inside \link{g.shell.GGIR}. Not intended for direct use by user.
}
\usage{
createConfigFile(config.parameters = c())	
}
\arguments{
  \item{config.parameters}{
List with all arguments used in g.shell.GGIR.
}
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}\name{load_params}
\alias{load_params}
\title{
  Load default parameters
}
\description{
  Loads default paramter values
  Not intended for direct use by GGIR users.
}
\usage{
  load_params(group = c("sleep", "metrics", "rawdata", "247",
                        "phyact", "cleaning", "output", "general"))
}
\arguments{
  \item{group}{
    Character vector with parameter groups to be loaded.
  }
}
\value{
  Lists of parameter objects 
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\name{is_this_a_dst_night}
\alias{is_this_a_dst_night}
\title{
 Check whether the night starting on a calendar date has DST.
}
\description{
Tests whether the night that follows the input calendar date is a
night with day saving time (DST) and on what hour the time moved.
}
\usage{
  is_this_a_dst_night(calendar_date=c(),tz="Europe/London")
}
\arguments{
  \item{calendar_date}{
    Character in the format dd/mm/yyyy 
  }
  \item{tz}{
    Time zone in "Europe/London" format.
  }
}
\value{
  \item{dst_night_or_not}{If value=0 no DST, if value=1 time moved forward,
  if value=-1 time moved forward}
  \item{dsthour}{Either the double hour or the hour that was skipped, 
  this differs between countries}
}
  
\examples{
  test4dst = is_this_a_dst_night("23/03/2014",tz="Europe/London")
}
\name{g.getbout}
\alias{g.getbout}
\title{
  function to calculate bouts from vector of binary classes
}
\description{
  To detect bouts of behaviour in time series. The function is used by \link{g.analyse}
}
\usage{
g.getbout(x, boutduration, boutcriter = 0.8, closedbout = FALSE, 
bout.metric=6, ws3=5)
}
\arguments{
  \item{x}{vector of zeros and/or ones to be screened for bouts of ones
}
  \item{boutduration}{duration of bout in epochs
}
  \item{boutcriter}{ Minimum percentage of boutduration for which the epoch values
  are expected to meet the threshold criterium
}
  \item{closedbout}{TRUE if you want breaks in bouts to be counted towards
  time spent in bouts (argument only active for bout.metric 1 and 2) 
}
  \item{bout.metric}{If value=1 the code uses the MVPA bout definition as has
  been available since 2014 (see papers by Sabia AJE 2014 and da Silva IJE 2014).
  Here, the algorithm looks for 10 minute windows in which more than XX percent
  of the epochs are above mvpathreshold, and then counts the entire window as mvpa.
  If value=2 the code looks for groups of epochs with a value above
  mvpathreshold that span a time window of at least mvpadur minutes in  which
  more than boutcriter percent of the epochs are above the threshold. The motivation 
  for the defition 1 was: A person who spends 10 minutes in MVPA with a 2 minute
  break in the middle is equally active as a person who spends 8 minutes in MVPA
  without taking a break. Therefore, both should be counted equal and counted as
  10 minute MVPA bout. The motivation for the definition 2 is: not counting breaks
  towards MVPA may simplify interpretation and still counts the two persons in
  the example as each others equal. If value=3, using sliding window across the
  data to test bout criteria per window and do not allow for breaks larger than 1 minute
  and with fraction of time larger than the boutcriter threshold.
  If value=4, same as 3 but also requires the first and last epoch to 
  meet the threshold criteria. If value=5, same as 4, but now looks for breaks larger
  than a minute such that 1 minute breaks are allowe, and the fraction of time that meets
  the threshold should be equal than or greater than the bout.criter threshold.
  If value=6, algorithm improved (2021) to check for first and last epoch.
}
  \item{ws3}{epoch length in seconds, only needed for bout.metric =3, because
  it needs to measure how many epochs equal 1 minute breaks
}

}
\value{
  \item{x}{Vector with binary numbers indicator where bouts where detected}
  \item{boutcount}{Vector with binary numbers indicator where bouts where detected
  and counted towards time spent in bouts, see  argument closedbout for clarification}
  
}
\examples{
y = g.getbout(x=round(runif(1000,0.4,1)),boutduration = 120,boutcriter=0.9,
    closedbout=FALSE,bout.metric=3,ws3=5)
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.part5.addfirstwake}
\alias{g.part5.addfirstwake}
\title{
  Adds first wake if it is missing in part 4 output.
}
\description{
  Not intended for direct use by GGIR users.
  Adds first wake if it is missing in part 4 output
  as part of \link{g.part5}.
}
\usage{
  g.part5.addfirstwake(ts, summarysleep_tmp2, nightsi, sleeplog,
  ID, Nepochsinhour, Nts, SPTE_end, ws3new)
}
\arguments{
  \item{ts}{
  }
  \item{summarysleep_tmp2}{
  }
  \item{nightsi}{
  }
  \item{sleeplog}{
  }
  \item{ID}{
  }
  \item{Nepochsinhour}{
  }
  \item{Nts}{
  }
  \item{SPTE_end}{
  }
  \item{ws3new}{
  }
}
\value{
  Data.frame ts
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.getmeta}
\alias{g.getmeta}
\title{
  Function to extract meta-data (features) from data in accelerometer file
}
\description{
  Reads a accelerometer file in blocks, extracts various features and stores average  feature
  value per short or long epoch. Acceleration and angle metrics are stored at short
  epoch length. The non-wear indication score, the clipping score, temperature
  (if available), light (if available), and Euclidean norm are stored at long epoch
  length. The function has been designed and thoroughly tested with accelerometer files
  from GENEA and GENEActiv. Further, the function should be able to cope with csv-format
  data procuded by GENEActiv and Actigraph
}
\usage{
  g.getmeta(datafile, params_metrics = c(), params_rawdata = c(),
                     params_general = c(), daylimit = FALSE, 
                     offset = c(0, 0, 0), scale = c(1, 1, 1), tempoffset = c(0, 0, 0),
                     meantempcal = c(), selectdaysfile = c(), myfun = c(), ...)
}
\arguments{
  \item{datafile}{
    name of accelerometer file
  }
  \item{params_metrics}{
    See \link{g.part1}
  }
  \item{params_rawdata}{
    See \link{g.part1}
  }
  \item{params_general}{
    See \link{g.part1}
  }
  \item{daylimit}{
    number of days to limit (roughly), if set to FALSE no daylimit
    will be applied
  }
  \item{offset}{
    offset correction value per axis, usage:
    value = scale(value,center = -offset, scale = 1/scale)
  }
  \item{scale}{
    scaling correction value per axis, usage:
    value = scale(value,center = -offset, scale = 1/scale)
  }
  \item{tempoffset}{
    temperature offset correction value per axis, usage:
    value = scale(value,center = -offset, scale = 1/scale)
    + scale(temperature, center = rep(averagetemperate,3), scale = 1/tempoffset)
  }
  \item{meantempcal}{
    mean temperature corresponding to the data as used for
    autocalibration. If autocalibration is not done or if temperature was not
    available then leave blank (default)
  }
  \item{selectdaysfile}{
    see \link{g.part1}
  }
  \item{myfun}{
    External function object to be applied to raw data.
    See details \link{applyExtFunction}.
  }
  \item{...}{
    Any argument used in the previous version of g.getmeta, which will now
    be used to overrule the arguments specified with the parameter objects.
  }
}
\value{
  \item{metalong}{dataframe with long epoch meta-data: EN, non-wear score,
  clipping score, temperature}
  \item{metashort}{dataframe with short epoch meta-data: timestamp and metric}
  \item{tooshort}{indicator of whether file was too short for processing (TRUE or FALSE)}
  \item{corrupt}{indicator of whether file was considered corrupt (TRUE or FALSE)}
}
\examples{
  \dontrun{
    datafile = "C:/myfolder/testfile.bin"
    
    #Extract meta-data:
    M = g.getmeta(datafile)
    
    #Inspect first couple of rows of long epoch length meta data:
    print(M$metalong[1:5,])
    
    #Inspect first couple of rows of short epoch length meta data:
    print(M$metalong[1:5,])
  }
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}
\references{
  \itemize{
    \item van Hees VT, Gorzelniak L, Dean Leon EC, Eder M, Pias M, et al. (2013) Separating
    Movement and Gravity Components in an Acceleration Signal and Implications for the
    Assessment of Human Daily Physical Activity. PLoS ONE 8(4): e61691.
    doi:10.1371/journal.pone.0061691
    \item Aittasalo M, Vaha-Ypya H, Vasankari T, Husu P, Jussila AM, and Sievanen H. Mean
    amplitude deviation calculated from raw acceleration data: a novel method for
    classifying the intensity of adolescents physical activity irrespective of accelerometer
    brand. BMC Sports Science, Medicine and Rehabilitation (2015).
  }
}
\name{g.detecmidnight}
\alias{g.detecmidnight}
\title{
  Detect all midnights in a time series
}
\description{
  Detect all midnights in a time series
}
\usage{
  g.detecmidnight(time,desiredtz, dayborder)
}
\arguments{
  \item{time}{
    Vector of timestamps, either in iso8601 or in POSIX format
  }
  \item{desiredtz}{
    See \link{g.part2}
  }
    \item{dayborder}{
  see \link{g.analyse}
  }
  
}
\value{
Output of the function is list containing the following objects:\cr
\itemize{
\item firstmidnight = timestamp of first midnight
\item firstmidnighti = index of first midnight
\item lastmidnight = timestamp of last midnight
\item lastmidnighti = index of last midnight
\item midnights = timestamps of midnights
\item midnightsi = indeces of midnights
}
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}\name{g.plot5}
\alias{g.plot5}
\title{
Generate user-friendly visual report. The first part of the report summarizes important daily metrics in bar plot format.
The second part of the report shows the raw data and annotations in 24-hr periods.
Angle-z is shown with sleep annotations during the SPT (sleep period time) window.
ENMO is shown with daytime inactivity and PA (physical activity) annotations in the lower
section of each 24-hr plot. The PA annotations are based on a 10 minute bout metric and
80% bout criteria. Moderate PA is a short window of time above threshold.mod that is part
of a 10 minute bout of MVPA. Vigorous PA is a short window of time above threshold.vig that
is part of a bout of MVPA. Light PA is a short window of time above threshold.lig that is
part of a bout of light PA.
}
\description{
  Function called by \link{g.shell.GGIR} to generate report. Not intended
  for direct use by user
}
\usage{
  g.plot5(metadatadir = c(), dofirstpage = FALSE, viewingwindow = 1,
  f0 = c(), f1 = c(), overwrite = FALSE, metric="ENMO",desiredtz = "Europe/London",
  threshold.lig = 30, threshold.mod = 100, threshold.vig = 400)
}
\arguments{
  \item{metadatadir}{
    See \link{g.part2}
  }
  \item{dofirstpage}{
    Boolean to indicate whether a first page with historgrams summarizing the whole
    measurement should be added
  }
  \item{viewingwindow}{
    See \link{g.shell.GGIR}
  }
  \item{f0}{
    See \link{g.shell.GGIR}
  }
  \item{f1}{
    See \link{g.shell.GGIR}
  }
  \item{overwrite}{
    See \link{g.shell.GGIR}
  }
  \item{metric}{
    Which one of the metrics do you want to consider to describe behaviour. The
    metric of interest need to be calculated in M (see \link{g.part1})
  }
  \item{desiredtz}{
    See \link{g.getmeta}
  }
  \item{threshold.lig}{
    See \link{g.part5}
  }
  \item{threshold.mod}{
    See \link{g.part5}
  }
  \item{threshold.vig}{
    See \link{g.part5}
  }
}
\value{
No values, this function only generates a plot
}

\examples{
  \dontrun{
  # generate plots for the first 10 files:
  g.plot5(metadatadir="C:/output_mystudy/meta/basic",dofirstpage=TRUE,
  viewingwindow = 1,f0=1,f1=10,overwrite=FALSE,desiredtz = "Europe/London",
  threshold.lig,threshold.mod,threshold.vig)
  }
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
Matthew R Patterson <mpatterson@shimmersensing.com>
}
\name{g.sib.sum}
\alias{g.sib.sum}
\title{
  sustiained inactivty bouts detection
}
\description{
  Detects sustiained inactivty bouts. Function not intended
  for direct use by package user
}
\usage{
  g.sib.sum(SLE,M,ignorenonwear=TRUE,desiredtz="")
}
\arguments{
  \item{SLE}{
    Output from \link{g.sib.det}
  }
  \item{M}{
    Object produced by \link{g.getmeta}
  }
  \item{ignorenonwear}{
    If TRUE then ignore detected monitor non-wear periods to avoid
  confusion between monitor non-wear time and sustained inactivity
  (default = TRUE)
  }
  \item{desiredtz}{
    See \link{g.part3}
  }
}
\value{
Dataframe with per night and per definition of sustained inactivity bouts
the start and end time of each sustained inactivity bout
}
\author{
Vincent T van Hees <v.vanhees@accelting.com>
}\name{getfolderstructure}
\alias{getfolderstructure}
\title{
Extracts folderstructure based on data directory.
}
\description{
Extracts folderstructure based on data directory. This is used when
accelerometer files are stored in a hierarchical folder structure and
the user likes to have a reference to the exact position in the folder
tree, rather than just the filename.  Function not
intended for direct use by package user.
}
\usage{
getfolderstructure(datadir=c(),referencefnames=c())	
}
\arguments{
  \item{datadir}{
Argument datadir as used in various other functions in GGIR
}
\item{referencefnames}{
vector with filename to filter on
}
}
\value{
List with items:
item{fullfilenames}{vector with all full paths to the folders 
  including the name of the file itself}
item{foldername}{vector with only the names of the folder in which each
file is stroed (so only the most distal folder in the folder tree)}
}
  
\examples{
\dontrun{
folderstructure = getfolderstructure(datadir)
}
}
\name{g.part5.fixmissingnight}
\alias{g.part5.fixmissingnight}
\title{
  Fix missing night in part 4 output
}
\description{
  Not intended for direct use by GGIR users.
  If a night is missing in the part4 output then this function
  tries to fix as part of \link{g.part5}.
}
\usage{
  g.part5.fixmissingnight(summarysleep_tmp2, sleeplog=c(), ID)
}
\arguments{
  \item{summarysleep_tmp2}{
    Object produced by \link{g.part4}
  }
  \item{sleeplog}{
  }
  \item{ID}{
  }
}
\value{
  Corrected summarysleep_tmp2 object.
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{read.myacc.csv}
\alias{read.myacc.csv}
\title{
  Read custom csv files with accelerometer data
}
\description{
  Loads csv files with accelerometer data and standardises 
  the output format (incl. unit of measurement, timestamp format,
  header format, and column locations) to make the data compatible
  with other GGIR functions.
}
\usage{
  read.myacc.csv(rmc.file=c(), rmc.nrow=c(), rmc.skip = c(), rmc.dec=".",
                          rmc.firstrow.acc = 1, rmc.firstrow.header=c(),
                          rmc.header.length = c(),
                          rmc.col.acc = 1:3, rmc.col.temp = c(), 
                          rmc.col.time=c(),
                          rmc.unit.acc = "g", rmc.unit.temp = "C", 
                          rmc.unit.time = "POSIX",
                          rmc.format.time = "\%Y-\%m-\%d \%H:\%M:\%OS",
                          rmc.bitrate = c(), rmc.dynamic_range = c(), 
                          rmc.unsignedbit = TRUE,
                          rmc.origin = "1970-01-01",
                          rmc.desiredtz = "Europe/London",
                          rmc.sf = c(),
                          rmc.headername.sf = c(),
                          rmc.headername.sn = c(),
                          rmc.headername.recordingid = c(),
                          rmc.header.structure = c(),
                          rmc.check4timegaps = FALSE,
                          rmc.col.wear = c(),
                          rmc.doresample=FALSE,
                          interpolationType=1)	
}
\arguments{
  \item{rmc.file}{
    Filename of file to be read.
  }
  \item{rmc.nrow}{
    Number of rows to read, same as nrow argument in read.csv and in fread.
  }
  \item{rmc.skip}{
    Number of rows to skip, same as skip argument in read.csv and in fread.
  }
  \item{rmc.dec}{
    Decimal used for numbers, same as skip argument in read.csv and in fread.
  }
  \item{rmc.firstrow.acc}{
    First row (number) of the acceleration data.
  }
  \item{rmc.firstrow.header}{
    First row (number) of the header. Leave blank if the file does not have
    a header.
  }
  \item{rmc.header.length}{
    If file has header, specify header length (numeric).
  }
  \item{rmc.col.acc}{
    Vector with three column (numbers) in which the acceleration signals
    are stored
  }
  \item{rmc.col.temp}{
    Scalar with column (number) in which the temperature is stored.
    Leave in default setting if no temperature is avaible. The temperature
    will be used by g.calibrate.
  }
  \item{rmc.col.time}{
    Scalar with column (number) in which the timestamps are stored.
    Leave in default setting if timestamps are not stored. 
  }
  \item{rmc.unit.acc}{
    Character with unit of acceleration values: "g", "mg", or "bit"
  }
  \item{rmc.unit.temp}{
    Character with unit of temperature values: (K)elvin, (C)elsius, or (F)ahrenheit
  }
  \item{rmc.unit.time}{
    Character with unit of timestamps: "POSIX",
    "UNIXsec" (seconds since origin, see argument origin), "character", or
    "ActivPAL" (exotic timestamp format only used in the ActivPAL
    activity monitor).
    
  }
  \item{rmc.format.time}{
    Format of timestamp, only used for rmc.unit.time: character and POSIX.
  }
  \item{rmc.bitrate}{
    Numeric: If unit of acceleration is a bit then provide bit rate, e.g. 12 bit.
  }
  \item{rmc.dynamic_range}{
    Numeric, if unit of acceleration is a bit then provide dynamic range deviation
    in g from zero, e.g. +/-6g would mean this argument needs to be 6. If you give this
    argument a character value the code will search the file header for elements with
    a name equal to the character value and use the corresponding numeric value
    next to it as dynamic range.
  }
  \item{rmc.unsignedbit}{
    Boolean, if unsignedbit = TRUE means that bits are only positive numbers.
   if unsignedbit = FALSE then bits are both positive and negative.
  }
  \item{rmc.origin}{
    Origin of time when unit of time is UNIXsec, e.g. 1970-1-1
  }
  \item{rmc.desiredtz}{
    Timezone in which device was configured and expriments took place.
    If experiments took place in a different timezone, then use this
    argument for the timezone in whcih the experiments took place and 
    argument configtz to specify where the device was configured (not implemented yet).
  }
  \item{rmc.sf}{
    Sample rate in Hertz, if this is stored in the file header then the that will be used
    instead.
  }
  \item{rmc.headername.sf}{
    If file has a header: Row name (character) under which the sample
    frequency can be found.
  }
  \item{rmc.headername.sn}{
    If file has a header: Row name (character) under which the
    serial number can be found.
  }
  \item{rmc.headername.recordingid}{
    If file has a header: Row name (character) under which the
    recording ID can be found.
  }
  
  \item{rmc.header.structure}{
    Character used to split the header name from the header
    value, e.g. ":" or " "
  }
  \item{rmc.check4timegaps}{
    Boolean to indicate whether gaps in time should be imputed with zeros.
    Some sensing equipment provides accelerometer with gaps in time. The rest of 
    GGIR is not designed for this, by setting this argument to TRUE the the gaps
    in time will be filled with zeros.
  }
  \item{rmc.col.wear}{
    If external wear detection outcome is stored as part of the data then this can be used by GGIR.
    This argument specifies the column in which the wear detection (Boolean) is stored.
  }
  \item{rmc.doresample}{
    Boolean to indicate whether to resample the data based on the available timestamps and extracted 
   sample rate from the file header
  }
  \item{interpolationType}{
    See \link{g.getmeta}
  }

}
\value{
  List with objects data holding the time series of acceleration, and
  header if it was available in the orignal file.
}
\examples{
  # create test files: No header, with temperature, with time
  N = 30
  sf = 30
  timestamps = as.POSIXlt(Sys.time()+((0:(N-1))/sf),origin="1970-1-1",tz = "Europe/London")
  mydata = data.frame(x=rnorm(N), time=timestamps,y=rnorm(N),z=rnorm(N),temp=rnorm(N)+20)
  testfile = "testcsv1.csv"
  write.csv(mydata, file= testfile, row.names = FALSE)
  loadedData = read.myacc.csv(rmc.file=testfile, rmc.nrow=20, rmc.dec=".",
                      rmc.firstrow.acc = 1, rmc.firstrow.header=c(),
                      rmc.col.acc = c(1,3,4), rmc.col.temp = 5, rmc.col.time=2,
                      rmc.unit.acc = "g", rmc.unit.temp = "C", rmc.origin = "1970-01-01")
  if (file.exists(testfile)) file.remove(testfile)
}
\author{
  Vincent T van Hees <v.vanhees@accelting.com>
}\name{isfilelist}
\alias{isfilelist}
\title{
Checks whether datadir is a directory or a vector with
filenames
}
\description{
Checks whether argument datadir used in various other functions in 
GGIR is the name of a directory that includes data files or whether
it is a vector with the full paths to one or more data files}
\usage{
isfilelist(datadir)	
}
\arguments{
  \item{datadir}{
Argument datadir as used in various other functions in GGIR
}
}
\value{
Boolean whether it is a list of files (TRUE) or not (FALSE)
}
  
\examples{
\dontrun{
isitafilelist = isfilelist(datadir)
}
}
