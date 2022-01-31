<!-- README.md is generated from README.Rmd. Please edit that file -->

dataaimsr <img src="man/figures/logo.png" width = 180 alt="dataaimsr Logo" align="right" />
===========================================================================================

<!-- badges: start -->

[![](https://badges.ropensci.org/428_status.svg)](https://github.com/ropensci/software-review/issues/428)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03282/status.svg)](https://doi.org/10.21105/joss.03282)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R build
status](https://github.com/ropensci/dataaimsr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/dataaimsr/actions)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/dataaimsr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/dataaimsr?branch=master)
![pkgdown](https://github.com/ropensci/dataaimsr/workflows/pkgdown/badge.svg)
[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)
[![packageversion](https://img.shields.io/badge/Package%20version-1.1.0-orange.svg)](commits/master)
[![Ask Us Anything
!](https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg)](https://github.com/ropensci/dataaimsr/issues/new)
![Open Source
Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)
<!-- badges: end -->

**Barneche DR, Coleman G, Fermor D, Klein E, Robinson T, Smith J,
Sheehan JL, Dowley S, Ditton D, Gunn K, Ericson G, Logan M, Rehbein M**
(2021). dataaimsr: An R Client for the Australian Institute of Marine
Science Data Platform API which provides easy access to AIMS Data
Platform. *Journal of Open Source Software*, **6:** 3282. doi:
[10.21105/joss.03282](https://doi.org/10.21105/joss.03282).

Overview
--------

The Australian Institute of Marine Science (AIMS) has a long tradition
in measuring and monitoring a series of environmental parameters along
the tropical coast of Australia. These parameters include long-term
record of sea surface temperature, wind characteristics, atmospheric
temperature, pressure, chlorophyll-a data, among many others. The AIMS
Data Centre team has recently developed the [AIMS Data Platform
API](https://open-aims.github.io/data-platform/) which is a *REST API*
providing JSON-formatted data to users. `dataaimsr` is an **R package**
written to allow users to communicate with the AIMS Data Platform API
using an API key and a few convenience functions to interrogate and
understand the datasets that are available to download. In doing so, it
allows the user to fully explore these datasets in R in whichever
capacity they want (e.g. data visualisation, statistical analyses, etc).
The package itself contains a `plot` method which allows the user to
plot summaries of the different types of dataset made available by the
API. Below we provide a brief context about the existing
[Datasets](#datasets) that can be explored through `dataaimsr`.

Installation
------------

### Requesting an AIMS Data Platform API Key

**AIMS Data Platform** requires an API Key for data requests, [get a key
here](https://open-AIMS.github.io/data-platform/key-request).

The API Key can be passed to the package functions as an additional
`api_key = "XXXX"` argument. **However**, we strongly encourage users to
maintain their API key as a private locally hidden environment variable
(`AIMS_DATAPLATFORM_API_KEY`) in the `.Renviron` file for automatic
loading at the start of an R session. Please read this
[article](https://CRAN.R-project.org/package=httr/vignettes/secrets.html)
which details why keeping your API private is extremely important.

Users can modify their `.Renviron` file by adding the following line:

    AIMS_DATAPLATFORM_API_KEY=XXXXXXXXXXXXX

The `.Renviron` file is usually stored in each users home directory:

<table>
<colgroup>
<col style="width: 35%" />
<col style="width: 64%" />
</colgroup>
<thead>
<tr class="header">
<th>System</th>
<th>.Renviron file locations</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>MS Windows</td>
<td><code>C:\Users\‹username›\.Renviron</code> or <code>C:\Users\‹username›\Documents\.Renviron</code></td>
</tr>
<tr class="even">
<td>Linux / MacOs</td>
<td><code>/home/‹username›/.Renviron</code></td>
</tr>
</tbody>
</table>

### Package

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr class="header">
<th>Type</th>
<th>Source</th>
<th>Command</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Release</td>
<td>CRAN</td>
<td>Not yet available</td>
</tr>
<tr class="even">
<td>Development</td>
<td>GitHub</td>
<td><code>remotes::install_github("ropensci/dataaimsr")</code></td>
</tr>
<tr class="odd">
<td>Development</td>
<td>rOpenSci</td>
<td><code>install.packages("dataaimsr", repos = "https://dev.ropensci.org")</code></td>
</tr>
</tbody>
</table>

Usage
-----

    # assumes that user already has API key saved to
    # .Renviron
    library(dataaimsr)

    # summarised by series
    # for all sites that contain data
    # within a defined date range
    sdf_b <- aims_data("temp_loggers", api_key = NULL,
                       summary = "summary-by-series",
                       filters = list("from_date" = "2018-01-01",
                                      "thru_date" = "2018-12-31"))

    # downloads weather data from site Yongala
    # within a defined date range
    wdf_a <- aims_data("weather", api_key = NULL,
                       filters = list(site = "Yongala",
                                      from_date = "2018-01-01",
                                      thru_date = "2018-01-02"))

More comprehensive examples about how to navigate `dataaimsr` and
interrogate the datasets can be found on our [online
vignettes](https://ropensci.github.io/dataaimsr/articles/).

Datasets
--------

Currently, there are two AIMS long-term monitoring datasets available to
be downloaded through `dataaimsr`:

### Northern Australia Automated Marine Weather And Oceanographic Stations

Automatic weather stations have been deployed by AIMS since 1980. Most
of the stations are along the Great Barrier Reef (GBR) including the
Torres Strait in North-Eastern Australia but there is also a station in
Darwin and one at Ningaloo Reef in Western Australia. Many of the
stations are located on the reef itself either on poles located in the
reef lagoon or on tourist pontoons or other structures. A list of the
weather stations which have been deployed by AIMS and the period of time
for which data may be available can be found on the
[metadata](https://apps.aims.gov.au/metadata/view/0887cb5b-b443-4e08-a169-038208109466)
webpage. **NB:** Records may not be continuous for the time spans given.

### AIMS Sea Water Temperature Observing System (AIMS Temperature Logger Program)

The data provided here are from a number of sea water temperature
monitoring programs conducted in tropical and subtropical coral reefs
environments around Australia. Data are available from approximately 80
GBR sites, 16 Coral Sea sites, 7 sites in North West Western Australia
(WA), 8 Queensland regional ports, 13 sites in the Solitary Islands, 4
sites in Papua New Guinea and 10 sites in the Cocos (Keeling) Islands.
Data are obtained from in-situ data loggers deployed on the reef.
Temperature instruments sample water temperatures every 5-10 minutes
(typically) and are exchanged and downloaded approximately every 12
months. Temperature loggers on the reef-flat are generally placed just
below Lowest Astronomical Tide level. Reef-slope (or where specified as
Upper reef-slope) generally refers to depths 5–9 m while Deep reef-slope
refers to depths of ~20 m. For more information on the dataset and its
usage, please visit the
[metadata](https://apps.aims.gov.au/metadata/view/4a12a8c0-c573-11dc-b99b-00008a07204e)
webpage.

License
-------

`dataaimsr` is provided by the [Australian Institute of Marine
Science](https://www.aims.gov.au) under the MIT License
([MIT](https://opensource.org/licenses/MIT)).

Code of Conduct
---------------

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

AIMS R package logos
--------------------

Our R package logos use a watercolour map of Australia, obtained with
the [ggmap](https://CRAN.R-project.org/package=ggmap) R package, which
downloads original map tiles provided by [Stamen
Design](https://stamen.com/), under [CC BY
3.0](https://creativecommons.org/licenses/by/3.0), with data from
[OpenStreetMap](https://www.openstreetmap.org/), under [CC BY
SA](https://creativecommons.org/licenses/by-sa/3.0).
# dataaimsr 1.1.0

* Added the capacity to download daily aggregated Sea Water Temperature Loggers dataset

# dataaimsr 1.0.3

* Added `aims_*` prefix to all exported functions

* `aims_data`, `aims_filter_values` and `aims_expose_attributes` now take a
string target rather than a string DOI as input argument

* `aims_data` and `aims_filter_values` fail gracefully now if the input
filter parameters are incorrect or there is no internet connection

* `aims_data_doi` is now the non-exported `data_doi`

* Important changes pertaining to `aims_data`:
 - Always returns a data.frame
 - Contains it's own class called `aimsdf` which contains print, summary and plot methods
 - Contains three additional exposed helper functions which allow the user to
 extract metadata/citation/parameter attributes.
 - Example set has been reduced to a minimal amount

* plot method for class `aimsdf` displays either a map or a time series

* improved test coverage

# dataaimsr 1.0.2

* Implemented `summary` datasets for the Temperature Loggers dataset via
`aims_data`

* Implemented `expose_attributes()` to show which filters are accepted
by the different datasets

* restrict `filter_values` to expel info on sites, series and parameters 
only

* Using `parsedate` to standardise date strings and account for time zone

* Created new vignette explaining basic usage of package
## How to contribute
Government employees, public and members of the private sector are encouraged to contribute to the repository by **forking and submitting a pull request**. 

(If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git) and  check out a more detailed guide to [pull requests](https://help.github.com/articles/using-pull-requests/).)

Pull requests will be evaluated by the repository guardians on a schedule and if deemed beneficial will be committed to the master.

All contributors retain the original copyright to their stuff, but by contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users **under the terms of the license under which this project is distributed.**
---
title: 'dataaimsr: An R Client for the Australian Institute of Marine Science Data Platform API which provides easy access to AIMS Data Platform'
tags:
  - R
  - sea surface temperature
  - weather stations
  - long-term environmental monitoring
  - Australia
  - API
authors:
  - name: Diego R. Barneche
    orcid: 0000-0002-4568-2362
    affiliation: "1, 2"
  - name: Greg Coleman
    affiliation: "3"
  - name: Duncan Fermor
    affiliation: "3"
  - name: Eduardo Klein
    affiliation: "3"
  - name: Tobias Robinson
    affiliation: "3"
  - name: Jason Smith
    affiliation: "3"
  - name: Jeffrey L. Sheehan
    affiliation: "3"
  - name: Shannon Dowley
    affiliation: "3"
  - name: Dean Ditton
    affiliation: "3"
  - name: Kevin Gunn
    affiliation: "3"
  - name: Gavin Ericson
    affiliation: "3"
  - name: Murray Logan
    affiliation: "3"
  - name: Mark Rehbein
    affiliation: "3"
affiliations:
 - name: Australian Institute of Marine Science, Crawley, WA 6009, Australia
   index: 1
 - name: The Indian Ocean Marine Research Centre, The University of Western Australia, Crawley, WA 6009, Australia
   index: 2
 - name: Australian Institute of Marine Science, Townsville, Qld 4810, Australia
   index: 3

citation_author: Barneche et al.
date: "2021-06-02"
bibliography: paper.bib
output:
  my_modified_joss:
    fig_caption: yes
csl: apa.csl
journal: JOSS
---

# Summary

`dataaimsr` is an **R package** written to provide open access to decades
of field measurements of atmospheric and oceanographic parameters around the
coast of Australia, conducted by the
[Australian Institute of Marine Science][0] (AIMS). The package communicates
with the recently-developed AIMS Data Platform API via an API key. Here we
describe the available datasets as well as example usage cases.

[0]: https://www.aims.gov.au/

# Statement of Need

The Australian Institute of Marine Science (AIMS) has a long tradition in
measuring and monitoring a series of environmental parameters along the
tropical coast of Australia. These parameters include long-term record of sea
surface temperature, wind characteristics, atmospheric temperature, pressure,
chlorophyll-a data, among many others. The AIMS Data Centre team has recently
developed the [AIMS Data Platform API][1] which is a *REST API* providing
JSON-formatted data to users. `dataaimsr` is an **R package** written to
allow users to communicate with the AIMS Data Platform API using an API key
and a few convenience functions to interrogate and understand the datasets
that are available to download. In doing so, it allows the user to
fully explore these datasets in R in whichever capacity they want (e.g.
data visualisation, statistical analyses, etc). The package itself contains
a `plot` method which allows the user to plot summaries of the different types
of dataset made available by the API.

[1]: https://open-aims.github.io/data-platform/

Currently, there are two AIMS long-term monitoring datasets available to be
downloaded through `dataaimsr`: 1) the Northern Australia Automated Marine
Weather And Oceanographic Stations---a list of the weather stations which have
been deployed by AIMS and the period of time for which data may be available
can be found on the [AIMS metadata][2] webpage; 2) AIMS Sea Water Temperature Observing System (AIMS Temperature Logger Program)---for more information on
the dataset and its usage, please visit the [AIMS metadata][3] webpage.

[2]: https://apps.aims.gov.au/metadata/view/0887cb5b-b443-4e08-a169-038208109466

[3]: https://apps.aims.gov.au/metadata/view/4a12a8c0-c573-11dc-b99b-00008a07204e 

# Technical details and Usage

Before loading the package, a user needs to download and store their personal
[AIMS Data Platform API Key][4]---we strongly encourage users to
maintain their API key as a private, locally hidden environment variable
(`AIMS_DATAPLATFORM_API_KEY`) in the `.Renviron` file for
automatic loading at the start of an R session.

[4]: https://open-aims.github.io/data-platform/key-request

`dataaimsr` imports the packages *httr* [@httrcit], *jsonlite* [@jsonlitecit],
*parsedate* [@parsedatecit], *dplyr* [@dplyrcit], *tidyr* [@tidyrcit],
*rnaturalearth* [@rnaturalearthcit], *sf* [@sfcit], *ggplot2* [@ggplot2cit],
*ggrepel* [@ggrepelcit] and *curl* [@curlcit].

The [Weather Station][2] and [Sea Water Temperature Loggers][3] datasets are
very large (terabytes in size), and as such they are not locally stored.
They are instead downloaded via the API and unique DOI identifiers. The 
datasets are structured by sites, series and parameters. A series is a 
continuing time-series, i.e. a collection of deployments measuring the 
same parameter (e.g. Air Temperature, Air Pressure, Chlorophyll) at the 
same subsite. So, for a given site and parameter, there might exist multiple
subsites and therefore series, in which case they are most likely 
distinguishable by depth.

For the Sea Water Temperature Loggers dataset, series is synonymous 
with the variable called subsite. For the Weather Station dataset, it 
is the combination of subsite and parameter.

## Discover a dataset

The [AIMS Data Platform API][1] points to the full metadata of each
dataset. We are currently working on ways to facilitate the 
visualisation of both datasets and their multiple features directly
through the R package. So please consult our [on-line vignettes][35] to obtain
the most up-to-date instructions on how to navigate the different datasets.
Future versions of this package might even provide more of AIMS monitoring
datasets.

[35]: https://docs.ropensci.org/dataaimsr/articles/

### Data summaries

The first step would be to visualise the dataset. We do this by
mapping all available sites. For example, we download the summary information
for the Sea Water Temperature Loggers dataset using the main function called
`aims_data`. Setting the argument `api_key = NULL` means that `dataaimsr` will
automatically search for the user's API key stored in `.Renviron`.
The `summary` argument should only be used when the
user wants an overview of the available data---this is currently
implemented for the Sea Water Temperature Loggers dataset only. One can
visualise `summary-by-series` or `summary-by-deployment`. The output of
`aims_data` is a `data.frame` of class `aimsdf` with its own plotting
method \autoref{fig:summary}:





![Distribution of all temperature logger series around Australian waters.\label{fig:summary}](summary.png)

For summary data such as `sdata`, plot will always generate a map with the
points around Australia and associated regions, coloured by the number of
calibrated observations.
Observations in a series can be: `uncal_obs`, `cal_obs` and `qc_obs`, which 
respectively stand for uncalibrated, calibrated, and quality-controlled 
observations. Calibrated and quality-controlled are generally the same.
Instruments are routinely calibrated (mostly once a year) in a
temperature-controlled water bath and corrections applied to the data. After
calibration, all data records are quality controlled based on the following
tests: 1) clip to in-water only data, using deployment's metadata, 2)
impossible value check: data outside a fixed temperature range (14˚C – 40˚C)
is flagged as bad data, 3) spike test: individual extreme values are flagged
as probably bad according to the algorithm presented in @morelo2014methods and
4) Excessive gradient test: pairs of data that present a sudden change in the
slope are flagged as probably bad [@toma2016acta]. If any data record fails at
least one of the tests, a QC flag equal to 2 is returned, otherwise, the QC
flag is set to 1. Please refer to our on-line [on-line vignettes][35] to learn
details about the entire structure of an `aimsdf` object.

In the case of the Weather Station dataset, the user can call a
the `aims_filter_values` function which allows one to query what
sites, series and parameters are available for both datasets:


```r
head(aims_filter_values("weather", filter_name = "series"))
```

```
##   series_id                                                                 series
## 1    104918        Myrmidon Reef Weather Station Wind Speed (scalar avg b 10 min) 
## 2    100686                            Saibai Island Weather Station Hail Duration
## 3       266 Orpheus Island Relay Pole 3 Wind Direction (Vector Average 30 Minutes)
## 4      2639 Hardy Reef Weather Station Wind Direction (Vector Standard 10 Minutes)
## 5     10243                           Raine Island Weather Station Air Temperature
## 6       258             Orpheus Island Relay Pole 3 Wind Speed (Scalar avg 10 min)
```

The downside is that one cannot know what time window is available
for each one of those, nor how they are nested (i.e. series /
parameter / site). In a way though the series name generally
gives that information anyway (see code output above). If knowing the 
available observation window is absolutely crucial, then as mentioned 
above the user should refer to the [on-line metadata][3].

## Download slices of datasets

We recommend slicing the datasets because AIMS monitoring datasets are of very 
high temporal resolution and if one tries to download an entire series
it might take a few hours. To slice the datasets properly, the user
needs to apply filters to their query.

### Data filters

Filters are the last important information the user needs to know to 
master the navigation and download of AIMS monitoring datasets. Each 
dataset can filtered by attributes which can be exposed with the function
`aims_expose_attributes`:


```r
aims_expose_attributes("weather")
```

```
## $summary
## [1] NA
## 
## $filters
##  [1] "site"      "subsite"   "series"    "series_id" "parameter" "size"      "min_lat"   "max_lat"   "min_lon"   "max_lon"   "from_date" "thru_date" "version"   "cursor"
```

```r
aims_expose_attributes("temp_loggers")
```

```
## $summary
## [1] "summary-by-series"     "summary-by-deployment"
## 
## $filters
##  [1] "site"      "subsite"   "series"    "series_id" "parameter" "size"      "min_lat"   "max_lat"   "min_lon"   "max_lon"   "from_date" "thru_date" "version"   "cursor"
```

The help file (see `?aims_expose_attributes`) contains the details about what
each filter targets. So, having an understanding of the summaries and what
filters are available provide the user with a great head start.

Downloading the data is achieved using the same `aims_data` function, 
however now the `summary` argument is omitted, and instead 
implement filters. For example, to download all the data collected at the
[Yongala wreck][6] for a specific time window:

[6]: https://en.wikipedia.org/wiki/SS_Yongala


```r
wdata_a <- aims_data("weather", api_key = NULL,
                     filters = list(site = "Yongala",
                                    from_date = "2018-01-01",
                                    thru_date = "2018-01-02"))
```

The returned `aimsdf` object in this case has attributes which give us
summary crucial information:

- `metadata` a doi link containing the metadata record for the data series

- `citation` the citation information for the particular dataset

- `parameters` an output `data.frame`

These can be directly extracted using the convenience functions
`aims_metadata`, `aims_citation` and `aims_parameters`, e.g.:


```r
aims_metadata(wdata_a)
```

```
## [1] "Metadata record https://doi.org/10.25845/5c09bf93f315d"
```

This example data contains multiple parameters available for this site at the
specified time, and the actual measurements are either raw or
quality-controlled. For monitoring data (i.e. when `summary = NA` in a
`aims_data` call), we can either visualise the data as a time series broken
down by parameter, or a map showing the sites with some summary info. If the
parameters are not specified, then `dataaimsr` will plot a maximum of 4
parameters chosen at random for a time series plot. Alternatively the user can
specify which parameters are to be plotted \autoref{fig:wind}.



![Yongala wreck profiles for water pressure and chlorophyll-a between the first and second of January 2018.\label{fig:wind}](wind.png)

The filters `from_date` and `thru_date` can be further refined by including a
time window to download the data:


```r
wdata_b <- aims_data("weather", api_key = NULL,
                     filters = list(series_id = 64,
                                    from_date = "1991-10-18T06:00:00",
                                    thru_date = "1991-10-18T12:00:00"))
range(wdata_b$time)
```

```
## [1] "1991-10-18 06:00:00 UTC" "1991-10-18 12:00:00 UTC"
```

### Methods

Objects of class `aimsdf` have associated `plot`, `print` and `summary`
methods.

### Data citation

Whenever using `dataaimsr`, we ask the user to not only cite this paper, but
also any data used in an eventual publication. Citation data can be extracted
from a dataset using the function `aims_citation`:


```r
aims_citation(wdata_b)
```

```
## [1] "Australian Institute of Marine Science (AIMS). 2009, Australian Institute of Marine Science Automatic Weather Stations, https://doi.org/10.25845/5c09bf93f315d, accessed 1 June 2021.  Time period: 1991-10-18T06:00:00 to 1991-10-18T12:00:00.  Series: Davies Reef Weather Station Air Temperature"
```

## Sister web tool

The Time Series Explorer (https://apps.aims.gov.au/ts-explorer/) is an
interactive web-based application that visualises large time series datasets.
The application utilises the AIMS Data Platform API to dynamically query data
according to user selection and visualise the data as line graphs. Series are
able to be compared visually. For large series, data are aggregated to daily
averages and displayed as minimum, maximum and mean. When the user 'zooms in'
sufficiently, the data will be displayed as non-aggregate values
\autoref{fig:tssa}. This technique is being used to ensure the application
performs well with large time series.

![Interactive discovery and visualisation of data series.\label{fig:tssa}](tssa.png)

The user can then download the displayed data as CSV or obtain a R code
snippet that shows how to obtain the data using the dataaimsr package
\autoref{fig:tssb}. In this way, a user can easily explore and discover
datasets and then quickly and easily have this data in their R environment for
additional analysis.

![Download/Export displayed data via R snippet.\label{fig:tssb}](tssb.png)

# Future directions

The API is still a work in progress. We are working on ways to better
facilitate data visualisation and retrieval, and also we are trying to
standardise the outputs from the different datasets as much as possible. In the
future, we envision that `dataaimsr` will also provide access to other
monitoring datasets collected by AIMS.

# References
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
stopifnot(require(knitr))
options(width = 90)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/",
  out.width = "100%"
)
```

```{r, echo = FALSE}
version <- as.vector(read.dcf("DESCRIPTION")[, "Version"])
version <- gsub("-", ".", version)
```

# dataaimsr <img src="man/figures/logo.png" width = 180 alt="dataaimsr Logo" align="right" />

<!-- badges: start -->
[![](https://badges.ropensci.org/428_status.svg)](https://github.com/ropensci/software-review/issues/428)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03282/status.svg)](https://doi.org/10.21105/joss.03282)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R build status](https://github.com/ropensci/dataaimsr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/dataaimsr/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/dataaimsr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/dataaimsr?branch=master)
![pkgdown](https://github.com/ropensci/dataaimsr/workflows/pkgdown/badge.svg)
[![license](https://img.shields.io/badge/license-MIT + file LICENSE-lightgrey.svg)](https://choosealicense.com/)
[![packageversion](https://img.shields.io/badge/Package%20version-`r version`-orange.svg)](commits/master)
[![Ask Us Anything
\!](https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg)](https://github.com/ropensci/dataaimsr/issues/new)
![Open Source
Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)
<!-- badges: end -->

**Barneche DR, Coleman G, Fermor D, Klein E, Robinson T, Smith J, Sheehan JL, Dowley S, Ditton D, Gunn K, Ericson G, Logan M, Rehbein M** (2021). dataaimsr: An R Client for the Australian Institute of Marine Science Data Platform API which provides easy access to AIMS Data Platform. *Journal of Open Source Software*, **6:** 3282. doi: [10.21105/joss.03282](https://doi.org/10.21105/joss.03282).

## Overview 

The Australian Institute of Marine Science (AIMS) has a long tradition in
measuring and monitoring a series of environmental parameters along the
tropical coast of Australia. These parameters include long-term record of sea
surface temperature, wind characteristics, atmospheric temperature, pressure,
chlorophyll-a data, among many others. The AIMS Data Centre team has recently
developed the [AIMS Data Platform API][1] which is a *REST API* providing
JSON-formatted data to users. `dataaimsr` is an **R package** written to
allow users to communicate with the AIMS Data Platform API using an API key
and a few convenience functions to interrogate and understand the datasets
that are available to download. In doing so, it allows the user to
fully explore these datasets in R in whichever capacity they want (e.g.
data visualisation, statistical analyses, etc). The package itself contains
a `plot` method which allows the user to plot summaries of the different types
of dataset made available by the API. Below we provide a brief context about
the existing [Datasets](#datasets) that can be explored through `dataaimsr`.

[1]: https://open-aims.github.io/data-platform/

## Installation

### Requesting an AIMS Data Platform API Key

**AIMS Data Platform** requires an API Key for data requests, [get a key here](https://open-AIMS.github.io/data-platform/key-request).

The API Key can be passed to the package functions as an additional
`api_key = "XXXX"` argument. **However**, we strongly encourage users to
maintain their API key as a private locally hidden environment variable
(`AIMS_DATAPLATFORM_API_KEY`) in the `.Renviron` file for
automatic loading at the start of an R session. Please read this
[article](https://CRAN.R-project.org/package=httr/vignettes/secrets.html)
which details why keeping your API private is extremely important.

Users can modify their `.Renviron` file by adding the following line:

```
AIMS_DATAPLATFORM_API_KEY=XXXXXXXXXXXXX
```

The `.Renviron` file is usually stored in each users home directory:

System        | .Renviron file locations
--------------|-------------------------
MS Windows    | <code>C:&#92;Users&#92;&#8249;username&#8250;&#92;.Renviron</code>  or <code>C:&#92;Users&#92;&#8249;username&#8250;&#92;Documents&#92;.Renviron</code>
Linux / MacOs | <code>&#47;home&#47;&#8249;username&#8250;&#47;.Renviron</code>

### Package

Type | Source | Command
---|---|---
Release | CRAN | Not yet available
Development | GitHub | `remotes::install_github("ropensci/dataaimsr")`
Development | rOpenSci | `install.packages("dataaimsr", repos = "https://dev.ropensci.org")`

## Usage

```{r, eval = FALSE}
# assumes that user already has API key saved to
# .Renviron
library(dataaimsr)

# summarised by series
# for all sites that contain data
# within a defined date range
sdf_b <- aims_data("temp_loggers", api_key = NULL,
                   summary = "summary-by-series",
                   filters = list("from_date" = "2018-01-01",
                                  "thru_date" = "2018-12-31"))

# downloads weather data from site Yongala
# within a defined date range
wdf_a <- aims_data("weather", api_key = NULL,
                   filters = list(site = "Yongala",
                                  from_date = "2018-01-01",
                                  thru_date = "2018-01-02"))
```

More comprehensive examples about how to navigate `dataaimsr` and interrogate
the datasets can be found on our [online vignettes][4].

[4]: https://ropensci.github.io/dataaimsr/articles/

## Datasets

Currently, there are two AIMS long-term monitoring datasets available to be
downloaded through `dataaimsr`:

### Northern Australia Automated Marine Weather And Oceanographic Stations

Automatic weather stations have been deployed by AIMS since 1980. Most of the
stations are along the Great Barrier Reef (GBR) including the Torres Strait in
North-Eastern Australia but there is also a station in Darwin and one at
Ningaloo Reef in Western Australia. Many of the stations are located on the
reef itself either on poles located in the reef lagoon or on tourist pontoons
or other structures. A list of the weather stations which have been deployed
by AIMS and the period of time for which data may be available can be
found on the [metadata][2] webpage. **NB:** Records may not be continuous for
the time spans given.

[2]: https://apps.aims.gov.au/metadata/view/0887cb5b-b443-4e08-a169-038208109466

### AIMS Sea Water Temperature Observing System (AIMS Temperature Logger Program)

The data provided here are from a number of sea water temperature monitoring
programs conducted in tropical and subtropical coral reefs environments around
Australia. Data are available from approximately 80 GBR sites, 16 Coral Sea
sites, 7 sites in North West Western Australia (WA), 8 Queensland regional
ports, 13 sites in the Solitary Islands, 4 sites in Papua New Guinea and 10
sites in the Cocos (Keeling) Islands. Data are obtained from in-situ data
loggers deployed on the reef. Temperature instruments sample water
temperatures every 5-10 minutes (typically) and are exchanged and downloaded
approximately every 12 months. Temperature loggers on the reef-flat are
generally placed just below Lowest Astronomical Tide level. Reef-slope (or
where specified as Upper reef-slope) generally refers to depths 5--9 m while
Deep reef-slope refers to depths of ~20 m. For more information on the dataset
and its usage, please visit the [metadata][3] webpage.

[3]: https://apps.aims.gov.au/metadata/view/4a12a8c0-c573-11dc-b99b-00008a07204e 

## License

`dataaimsr` is provided by the [Australian Institute of Marine Science](https://www.aims.gov.au) under the MIT License ([MIT](https://opensource.org/licenses/MIT)).

## Code of Conduct

Please note that this package is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/).
By contributing to this project, you agree to abide by its terms.

## AIMS R package logos

Our R package logos use a watercolour map of Australia, obtained with the [ggmap](https://CRAN.R-project.org/package=ggmap) R package, which downloads original map tiles provided by [Stamen Design](https://stamen.com/), under [CC BY 3.0](https://creativecommons.org/licenses/by/3.0), with data from [OpenStreetMap](https://www.openstreetmap.org/), under [CC BY SA](https://creativecommons.org/licenses/by-sa/3.0).
---
title: 'dataaimsr: An R Client for the Australian Institute of Marine Science Data Platform API which provides easy access to AIMS Data Platform'
tags:
  - R
  - sea surface temperature
  - weather stations
  - long-term environmental monitoring
  - Australia
  - API
authors:
  - name: Diego R. Barneche
    orcid: 0000-0002-4568-2362
    affiliation: "1, 2"
  - name: Greg Coleman
    affiliation: "3"
  - name: Duncan Fermor
    affiliation: "3"
  - name: Eduardo Klein
    affiliation: "3"
  - name: Tobias Robinson
    affiliation: "3"
  - name: Jason Smith
    affiliation: "3"
  - name: Jeffrey L. Sheehan
    affiliation: "3"
  - name: Shannon Dowley
    affiliation: "3"
  - name: Dean Ditton
    affiliation: "3"
  - name: Kevin Gunn
    affiliation: "3"
  - name: Gavin Ericson
    affiliation: "3"
  - name: Murray Logan
    affiliation: "3"
  - name: Mark Rehbein
    affiliation: "3"
affiliations:
 - name: Australian Institute of Marine Science, Crawley, WA 6009, Australia
   index: 1
 - name: The Indian Ocean Marine Research Centre, The University of Western Australia, Crawley, WA 6009, Australia
   index: 2
 - name: Australian Institute of Marine Science, Townsville, Qld 4810, Australia
   index: 3

citation_author: Barneche et al.
date: "`r Sys.Date()`"
bibliography: paper.bib
output:
  my_modified_joss:
    fig_caption: yes
csl: apa.csl
journal: JOSS
---

# Summary

`dataaimsr` is an **R package** written to provide open access to decades
of field measurements of atmospheric and oceanographic parameters around the
coast of Australia, conducted by the
[Australian Institute of Marine Science][0] (AIMS). The package communicates
with the recently-developed AIMS Data Platform API via an API key. Here we
describe the available datasets as well as example usage cases.

[0]: https://www.aims.gov.au/

# Statement of Need

The Australian Institute of Marine Science (AIMS) has a long tradition in
measuring and monitoring a series of environmental parameters along the
tropical coast of Australia. These parameters include long-term record of sea
surface temperature, wind characteristics, atmospheric temperature, pressure,
chlorophyll-a data, among many others. The AIMS Data Centre team has recently
developed the [AIMS Data Platform API][1] which is a *REST API* providing
JSON-formatted data to users. `dataaimsr` is an **R package** written to
allow users to communicate with the AIMS Data Platform API using an API key
and a few convenience functions to interrogate and understand the datasets
that are available to download. In doing so, it allows the user to
fully explore these datasets in R in whichever capacity they want (e.g.
data visualisation, statistical analyses, etc). The package itself contains
a `plot` method which allows the user to plot summaries of the different types
of dataset made available by the API.

[1]: https://open-aims.github.io/data-platform/

Currently, there are two AIMS long-term monitoring datasets available to be
downloaded through `dataaimsr`: 1) the Northern Australia Automated Marine
Weather And Oceanographic Stations---a list of the weather stations which have
been deployed by AIMS and the period of time for which data may be available
can be found on the [AIMS metadata][2] webpage; 2) AIMS Sea Water Temperature Observing System (AIMS Temperature Logger Program)---for more information on
the dataset and its usage, please visit the [AIMS metadata][3] webpage.

[2]: https://apps.aims.gov.au/metadata/view/0887cb5b-b443-4e08-a169-038208109466

[3]: https://apps.aims.gov.au/metadata/view/4a12a8c0-c573-11dc-b99b-00008a07204e 

# Technical details and Usage

Before loading the package, a user needs to download and store their personal
[AIMS Data Platform API Key][4]---we strongly encourage users to
maintain their API key as a private, locally hidden environment variable
(`AIMS_DATAPLATFORM_API_KEY`) in the `.Renviron` file for
automatic loading at the start of an R session.

[4]: https://open-aims.github.io/data-platform/key-request

`dataaimsr` imports the packages *httr* [@httrcit], *jsonlite* [@jsonlitecit],
*parsedate* [@parsedatecit], *dplyr* [@dplyrcit], *tidyr* [@tidyrcit],
*rnaturalearth* [@rnaturalearthcit], *sf* [@sfcit], *ggplot2* [@ggplot2cit],
*ggrepel* [@ggrepelcit] and *curl* [@curlcit].

The [Weather Station][2] and [Sea Water Temperature Loggers][3] datasets are
very large (terabytes in size), and as such they are not locally stored.
They are instead downloaded via the API and unique DOI identifiers. The 
datasets are structured by sites, series and parameters. A series is a 
continuing time-series, i.e. a collection of deployments measuring the 
same parameter (e.g. Air Temperature, Air Pressure, Chlorophyll) at the 
same subsite. So, for a given site and parameter, there might exist multiple
subsites and therefore series, in which case they are most likely 
distinguishable by depth.

For the Sea Water Temperature Loggers dataset, series is synonymous 
with the variable called subsite. For the Weather Station dataset, it 
is the combination of subsite and parameter.

## Discover a dataset

The [AIMS Data Platform API][1] points to the full metadata of each
dataset. We are currently working on ways to facilitate the 
visualisation of both datasets and their multiple features directly
through the R package. So please consult our [on-line vignettes][35] to obtain
the most up-to-date instructions on how to navigate the different datasets.
Future versions of this package might even provide more of AIMS monitoring
datasets.

[35]: https://docs.ropensci.org/dataaimsr/articles/

### Data summaries

The first step would be to visualise the dataset. We do this by
mapping all available sites. For example, we download the summary information
for the Sea Water Temperature Loggers dataset using the main function called
`aims_data`. Setting the argument `api_key = NULL` means that `dataaimsr` will
automatically search for the user's API key stored in `.Renviron`.
The `summary` argument should only be used when the
user wants an overview of the available data---this is currently
implemented for the Sea Water Temperature Loggers dataset only. One can
visualise `summary-by-series` or `summary-by-deployment`. The output of
`aims_data` is a `data.frame` of class `aimsdf` with its own plotting
method \autoref{fig:summary}:

```{r, echo = FALSE, message = FALSE, warning = FALSE}
library(dataaimsr)
sdata <- aims_data("temp_loggers", api_key = NULL,
                   summary = "summary-by-series")
```

```{r echo = FALSE, eval = FALSE}
library(dataaimsr)
sdata <- aims_data("temp_loggers", api_key = NULL,
                   summary = "summary-by-series")
out <- plot(sdata, ptype = "map")
ggsave("paper/summary.png", out, width = 7, height = 4.36, units = "in",
       device = "png", dpi = 300)
```

![Distribution of all temperature logger series around Australian waters.\label{fig:summary}](summary.png)

For summary data such as `sdata`, plot will always generate a map with the
points around Australia and associated regions, coloured by the number of
calibrated observations.
Observations in a series can be: `uncal_obs`, `cal_obs` and `qc_obs`, which 
respectively stand for uncalibrated, calibrated, and quality-controlled 
observations. Calibrated and quality-controlled are generally the same.
Instruments are routinely calibrated (mostly once a year) in a
temperature-controlled water bath and corrections applied to the data. After
calibration, all data records are quality controlled based on the following
tests: 1) clip to in-water only data, using deployment's metadata, 2)
impossible value check: data outside a fixed temperature range (14˚C – 40˚C)
is flagged as bad data, 3) spike test: individual extreme values are flagged
as probably bad according to the algorithm presented in @morelo2014methods and
4) Excessive gradient test: pairs of data that present a sudden change in the
slope are flagged as probably bad [@toma2016acta]. If any data record fails at
least one of the tests, a QC flag equal to 2 is returned, otherwise, the QC
flag is set to 1. Please refer to our on-line [on-line vignettes][35] to learn
details about the entire structure of an `aimsdf` object.

In the case of the Weather Station dataset, the user can call a
the `aims_filter_values` function which allows one to query what
sites, series and parameters are available for both datasets:

```{r, message = FALSE, warning = FALSE}
head(aims_filter_values("weather", filter_name = "series"))
```

The downside is that one cannot know what time window is available
for each one of those, nor how they are nested (i.e. series /
parameter / site). In a way though the series name generally
gives that information anyway (see code output above). If knowing the 
available observation window is absolutely crucial, then as mentioned 
above the user should refer to the [on-line metadata][3].

## Download slices of datasets

We recommend slicing the datasets because AIMS monitoring datasets are of very 
high temporal resolution and if one tries to download an entire series
it might take a few hours. To slice the datasets properly, the user
needs to apply filters to their query.

### Data filters

Filters are the last important information the user needs to know to 
master the navigation and download of AIMS monitoring datasets. Each 
dataset can filtered by attributes which can be exposed with the function
`aims_expose_attributes`:

```{r, message = FALSE, warning = FALSE}
aims_expose_attributes("weather")
aims_expose_attributes("temp_loggers")
```

The help file (see `?aims_expose_attributes`) contains the details about what
each filter targets. So, having an understanding of the summaries and what
filters are available provide the user with a great head start.

Downloading the data is achieved using the same `aims_data` function, 
however now the `summary` argument is omitted, and instead 
implement filters. For example, to download all the data collected at the
[Yongala wreck][6] for a specific time window:

[6]: https://en.wikipedia.org/wiki/SS_Yongala

```{r, message = FALSE, warning = FALSE}
wdata_a <- aims_data("weather", api_key = NULL,
                     filters = list(site = "Yongala",
                                    from_date = "2018-01-01",
                                    thru_date = "2018-01-02"))
```

The returned `aimsdf` object in this case has attributes which give us
summary crucial information:

- `metadata` a doi link containing the metadata record for the data series

- `citation` the citation information for the particular dataset

- `parameters` an output `data.frame`

These can be directly extracted using the convenience functions
`aims_metadata`, `aims_citation` and `aims_parameters`, e.g.:

```{r}
aims_metadata(wdata_a)
```

This example data contains multiple parameters available for this site at the
specified time, and the actual measurements are either raw or
quality-controlled. For monitoring data (i.e. when `summary = NA` in a
`aims_data` call), we can either visualise the data as a time series broken
down by parameter, or a map showing the sites with some summary info. If the
parameters are not specified, then `dataaimsr` will plot a maximum of 4
parameters chosen at random for a time series plot. Alternatively the user can
specify which parameters are to be plotted \autoref{fig:wind}.

```{r echo = FALSE, eval = FALSE}
# check parameters with aims_parameters(wdata_a)
out <- plot(wdata_a, ptype = "time_series",
            pars = c("Water Pressure", "Chlorophyll"))
ggsave("paper/wind.png", out, width = 8.5, height = 3.9, units = "in",
       device = "png", dpi = 300)
```

![Yongala wreck profiles for water pressure and chlorophyll-a between the first and second of January 2018.\label{fig:wind}](wind.png)

The filters `from_date` and `thru_date` can be further refined by including a
time window to download the data:

```{r, message = FALSE, warning = FALSE}
wdata_b <- aims_data("weather", api_key = NULL,
                     filters = list(series_id = 64,
                                    from_date = "1991-10-18T06:00:00",
                                    thru_date = "1991-10-18T12:00:00"))
range(wdata_b$time)
```

### Methods

Objects of class `aimsdf` have associated `plot`, `print` and `summary`
methods.

### Data citation

Whenever using `dataaimsr`, we ask the user to not only cite this paper, but
also any data used in an eventual publication. Citation data can be extracted
from a dataset using the function `aims_citation`:

```{r, message = FALSE, warning = FALSE}
aims_citation(wdata_b)
```

## Sister web tool

The Time Series Explorer (https://apps.aims.gov.au/ts-explorer/) is an
interactive web-based application that visualises large time series datasets.
The application utilises the AIMS Data Platform API to dynamically query data
according to user selection and visualise the data as line graphs. Series are
able to be compared visually. For large series, data are aggregated to daily
averages and displayed as minimum, maximum and mean. When the user 'zooms in'
sufficiently, the data will be displayed as non-aggregate values
\autoref{fig:tssa}. This technique is being used to ensure the application
performs well with large time series.

![Interactive discovery and visualisation of data series.\label{fig:tssa}](tssa.png)

The user can then download the displayed data as CSV or obtain a R code
snippet that shows how to obtain the data using the dataaimsr package
\autoref{fig:tssb}. In this way, a user can easily explore and discover
datasets and then quickly and easily have this data in their R environment for
additional analysis.

![Download/Export displayed data via R snippet.\label{fig:tssb}](tssb.png)

# Future directions

The API is still a work in progress. We are working on ways to better
facilitate data visualisation and retrieval, and also we are trying to
standardise the outputs from the different datasets as much as possible. In the
future, we envision that `dataaimsr` will also provide access to other
monitoring datasets collected by AIMS.

# References
---
title: "Navigating dataaimsr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Navigating dataaimsr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



The very first thing to do is read the documentation on our
[README](../index.html) page. Make sure you have the package properly 
installed, and that your personal [AIMS Data Platform API Key][1] has 
been downloaded.

[0]: https://www.aims.gov.au/
[1]: https://open-aims.github.io/data-platform/key-request

As per the installation instructions, we strongly suggest that you hide your
API Key permanently in your `.Renviron` file and set the object `my_api_key` to
`NULL` in the chunk below. You can read more about why that is important
[here](https://CRAN.R-project.org/package=httr/vignettes/secrets.html).


```r
# set my_api_key to NULL after successfully placing it in .Renviron
my_api_key <- NULL
```

We now load `dataaimsr`:


```r
library(dataaimsr)
```

## How this package works

`dataaimsr` contains two sets of monitoring data collected by AIMS---[the
Australian Institute of Marine Science][0]---since the 1980's: the
[Weather Station][2] dataset which contains 
encompasses data for different parameters (e.g. Air Temperature, Air 
Pressure, Chlorophyll, and many others); and the
[Sea Water Temperature Loggers][3] dataset which contains records of 
(you guessed it!) sea water temperature at different sites and water 
depths.

[2]: https://doi.org/10.25845/5c09bf93f315d
[3]: https://doi.org/10.25845/5b4eb0f9bb848

The datasets are very large, and as such they are not locally stored.
They are instead downloaded via the API and unique DOI identifier (just 
hover over the data links above to see the actual DOI codes). The 
datasets are structured by sites, series and parameters. A series is a 
continuing time-series, i.e. a collection of deployments measuring the 
same parameter (e.g. Air Temperature, Air Pressure, Chlorophyll) at the 
same subsite. So, for a given site and parameter, there might exist multiple
subsites and therefore series, in which case they are most likely 
distinguishable by depth.

For the Sea Water Temperature Loggers dataset, series is synonymous 
with the variable called subsite. For the Weather Station dataset, it 
is the combination of subsite and parameter.

This vignette gives an overview of how one would go about discovering
the overall information contained in the datasets. For dataset-specific 
vignettes, see our other [vignette pages][4].

[4]: https://ropensci.github.io/dataaimsr/articles/

## Discover a dataset

The [AIMS Data Platform API][5] points to the full metadata of each
dataset. We are currently working on ways to facilitate the 
visualisation of both datasets and their multiple features directly
through the R package. At the moment though it is only possible
to visualise summary information for the Sea Water Temperature Loggers
dataset. A similar feature for the Weather Station dataset will be 
implemented in the near future (likely early 2021)---so for now, please
refer to the online metadata to discover from where (and when) you can 
download data.

[5]: https://open-aims.github.io/data-platform/

### Data summary

The first step would be to visualise the dataset. Let's do this by
mapping all available sites for the Sea Water Temperature 
Loggers dataset using the main function called `aims_data`:


```r
sdata <- aims_data(target = "temp_loggers", api_key = my_api_key,
                   summary = "summary-by-series")
head(sdata)
#>   site_id                    site subsite_id    subsite series_id     series
#> 1       1 Agincourt Reef Number 3       2687     AG3FL1      2687     AG3FL1
#> 2       1 Agincourt Reef Number 3      14276  AG3SL1old     14276  AG3SL1old
#> 3       3           Cleveland Bay       3007 CLEVAWSSL1      3007 CLEVAWSSL1
#> 4       3           Cleveland Bay       3069 CLEVAWSFL1      3069 CLEVAWSFL1
#> 5       4             Davies Reef       2629     DAVFL1      2629     DAVFL1
#> 6       4             Davies Reef       2630     DAVSL1      2630     DAVSL1
#>           parameter parameter_id time_coverage_start time_coverage_end      lat      lon
#> 1 Water Temperature            1          1996-03-30        2008-12-11 -15.9903 145.8212
#> 2 Water Temperature            1          1996-03-30        2011-07-21 -15.9905 145.8213
#> 3 Water Temperature            1          2004-05-13        2008-05-03 -19.1557 146.8813
#> 4 Water Temperature            1          2005-09-15        2005-12-22 -19.1557 146.8813
#> 5 Water Temperature            1          1997-08-26        2019-06-10 -18.8065 147.6688
#> 6 Water Temperature            1          1996-05-02        2021-03-29 -18.8060 147.6686
#>   depth uncal_obs cal_obs qc_obs
#> 1     0     23130  110480 110480
#> 2     5    114450  216794 216794
#> 3     7     11951   53231  53231
#> 4     1         0    4656   4656
#> 5     1    437544  566585 566585
#> 6     8    463146  589492 589437
```

The `summary` argument here is key. It should be either `"summary-by-series"`
or `"summary-by-deployment"` when the user wants an overview of the available
data. Again, this is currently implemented for the Sea Water Temperature Loggers
dataset only. The output of `aims_data` is a `data.frame` of class `aimsdf`.

Notice that `sdata` contains a lot of information, most of which is
related to site / series / parameter ID. Each row corresponds to a
unique series, and a certain site may contain multiple series; in such
cases, series generally differ from one another by depth. The columns 
`time_coverage_start` and `time_coverage_end` are probably one of the most
valuable pieces of information. They provide the user with the window of data
collection for a particular series, which is probably crucial to decide
whether that particular series is of relevance to the specific question in
hand.

Also note that there are three columns containing the total number of 
observations in a series: `uncal_obs`, `cal_obs` and `qc_obs`, which 
respectively stand for uncalibrated, calibrated, and quality-controlled 
observations. Calibrated and quality-controlled are generally the same.
Instruments are routinely calibrated (mostly once a year) in a
temperature-controlled water bath and corrections applied to the data. After
calibration, all data records are quality controlled based on the following
tests: 1) clip to in-water only data, using deployment's metadata, 2)
impossible value check: data outside a fixed temperature range (14˚C – 40˚C)
is flagged as bad data, 3) spike test: individual extreme values are flagged
as probably bad and 4) Excessive gradient test: pairs of data that present a
sudden change in the slope are flagged as probably bad. If any data record
fails at least one of the tests, a QC flag equal to 2 is returned, otherwise,
the QC flag is set to 1.

`aimsdf` objects can be plotted using the `plot` function. For `summary-by-...`
data such as `sdata`, `plot` will always generate a map with the points around
Australia and associated regions, coloured by the number of calibrated
observations:


```r
plot(sdata, ptype = "map")
```

<img src="vignette-fig-summap-1.png" title="plot of chunk summap" alt="plot of chunk summap" width="100%" />

### Filter values

In the case of the Weather Station dataset, knowing what sites are
out there is a bit tricky. However, currently we have a convenience
function called `aims_filter_values` which allows one to query what
sites, series and parameters are available for both datasets:


```r
head(aims_filter_values("weather", filter_name = "series"))
#>   series_id                                                                 series
#> 1    104918        Myrmidon Reef Weather Station Wind Speed (scalar avg b 10 min) 
#> 2    100686                            Saibai Island Weather Station Hail Duration
#> 3       266 Orpheus Island Relay Pole 3 Wind Direction (Vector Average 30 Minutes)
#> 4      2639 Hardy Reef Weather Station Wind Direction (Vector Standard 10 Minutes)
#> 5     10243                           Raine Island Weather Station Air Temperature
#> 6       258             Orpheus Island Relay Pole 3 Wind Speed (Scalar avg 10 min)
```

The downside is that one cannot know what time window is available
for each one of those, nor how they are nested (i.e. series /
parameter / site). In a way though the series name generally
gives that information anyway (see code output above). If knowing the 
available observation window is absolutely crucial, then as mentioned 
above the user should refer to the [online metadata][5].

## Download slices of datasets

Now that we know how to explore the datasets and what data is out there,
we finish this vignette by showing an example of how one would go about
downloading actual monitoring data.

We say slices of datasets because AIMS monitoring datasets are of very 
high temporal resolution and if one tries to download the entire thing
it might take hours. Generally that is why we download slices of data at a
time, and for that we need filters (see below).

On the other hand, if all we are interested in are aggregated values
(daily means), then we can set `summary = "daily"` in `aims_data` in
combination with a list of `filters` to download more concise datasets.
So far, aggregated data is only available for the Sea Water Temperature
Loggers dataset.

### Data filters

Filters are the last important information the user needs to know to 
master the navigation and download of AIMS monitoring datasets. Each 
dataset can filtered by attributes which can be exposed with the function `aims_expose_attributes`:


```r
aims_expose_attributes("weather")
#> $summary
#> [1] NA
#> 
#> $filters
#>  [1] "site"      "subsite"   "series"    "series_id" "parameter" "size"      "min_lat"  
#>  [8] "max_lat"   "min_lon"   "max_lon"   "from_date" "thru_date" "version"   "cursor"
aims_expose_attributes("temp_loggers")
#> $summary
#> [1] "summary-by-series"     "summary-by-deployment" "daily"                
#> 
#> $filters
#>  [1] "site"      "subsite"   "series"    "series_id" "parameter" "size"      "min_lat"  
#>  [8] "max_lat"   "min_lon"   "max_lon"   "from_date" "thru_date" "version"   "cursor"
```

The help file (see `?aims_expose_attributes`) contains the details about what
each filter targets. So, having an understanding of the summaries and what
filters are available provide the user with a great head start.

Downloading the raw, high-resolution data is achieved using the same `aims_data`
function, however now we ignore the `summary` argument, and instead 
implement filters. For example, let's say we want to download all the
data collected at the [Yongala](https://en.wikipedia.org/wiki/SS_Yongala) for
a specific time window:


```r
wdata_a <- aims_data("weather", api_key = my_api_key,
                     filters = list(site = "Yongala",
                                    from_date = "2018-01-01",
                                    thru_date = "2018-01-02"))
```

The returned `aimsdf` object in this case has attributes which give us
summary crucial information:

- `metadata` a doi link containing the metadata record for the data series

- `citation` the citation information for the particular dataset

- `parameters` an output `data.frame`

These can be directly extracted using the convenience functions
`aims_metadata`, `aims_citation` and `aims_parameters`, e.g.:


```r
aims_metadata(wdata_a)
#> [1] "Metadata record https://doi.org/10.25845/5c09bf93f315d"
```

This example data contains multiple parameters available for this site at the
specified time, and the actual measurements are either raw or
quality-controlled. For monitoring data (both aggregated or non-aggregated,
i.e. when `summary = NA` or `summary = "daily"` in an `aims_data` call), we can
either visualise the data as a time series broken down by parameter, or a map
showing the sites with some summary info. If the parameters are not specified,
then `dataaimsr` will plot a maximum of 4 parameters chosen at random for a
time series plot. Alternatively the user can specify which parameters are to be
plotted.


```r
# check parameters with aims_parameters(wdata_a)
plot(wdata_a, ptype = "time_series",
     pars = c("Water Pressure", "Chlorophyll"))
```

<img src="vignette-fig-monts-1.png" title="plot of chunk monts" alt="plot of chunk monts" width="100%" />

We can also refine even further by including a time window to download the
data:


```r
wdata_b <- aims_data("weather",
                     api_key = my_api_key,
                     filters = list(series_id = 64,
                                    from_date = "1991-10-18T06:00:00",
                                    thru_date = "1991-10-18T12:00:00"))
range(wdata_b$time)
#> [1] "1991-10-18 06:00:00 UTC" "1991-10-18 12:00:00 UTC"
```

Or simply download and plot daily aggregated data:


```r
sdata_c <- aims_data("temp_loggers", api_key = my_api_key, summary = "daily",
                     filters = list(series = "DAVFL1",
                                    from_date = "2018-01-01",
                                    thru_date = "2018-12-31"))
plot(sdata_c, ptype = "time_series", pars = c("Water Temperature"))
```

<img src="vignette-fig-unnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" width="100%" />

## More info

See our other [vignette pages][4] for further dataset-specific 
explorations.
---
title: "Sea Water Temperature Loggers time series dataset"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Sea Water Temperature Loggers time series dataset}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



Please check our [intro vignette][1] first to implement the installation
requirements, and to learn the general approach to navigating the different
datasets. This vignette assumes you have obtained an
[AIMS Data Platform API Key][2].

[1]: https://ropensci.github.io/dataaimsr/articles/navigating.html
[2]: https://open-aims.github.io/data-platform/key-request

As per the installation instructions, we strongly suggest that you hide your
API Key permanently in your `.Renviron` file and set the object `my_api_key` to
`NULL` in the chunk below. You can read more about why that is important
[here](https://CRAN.R-project.org/package=httr/vignettes/secrets.html).


```r
# set my_api_key to NULL after successfully placing it in .Renviron
my_api_key <- NULL
```

Let's start by loading the packages needed for this vignette:


```r
library(purrr)
library(dataaimsr)
library(ggplot2)
```

## Discovering the dataset

The [Sea Water Temperature Loggers][3] dataset is less extensive than the
[AIMS Weather Station][4] dataset because it comprises one single
*"parameter"*---water temperature---that is measured at multiple sites. Not 
all sites have the same temporal coverage; some loggers are still actively 
collecting data, others have been discontinued. So the key distinctive 
variables in this instance are the "site", and the "series". A "series" 
represents a continuing time-series, i.e. a collection of deployments 
measuring the same parameter at the same subsite. Because there is only one
parameter (water temperature), subsite and series are synonymous in the
[Sea Water Temperature Loggers][3] dataset. So a series will comprise a
continuing time-series at a specific site and depth.

Essentially, for the user who has limited knowledge about where the data are,
and of what they are consisted, they would need to do some prior exploration 
to learn more about what can be downloaded. Suppose the goal is to download all
time-series from a particular site. The general procedure would be:

1. Examine documentation and establish query filters
2. Perform data download using `aims_data`
3. Create an exploratory time-series chart

For all datasets, a list of available filters can be retrieved with the 
function `aims_expose_attributes`. Knowing the filters is important because
some time series are quite extensive, with parameters being measured at very
high frequency (e.g. every 5 minutes), so downloading the dataset for an
entire year or more my take quite some time (it's possible though if that is
the true goal of the user). Otherwise, the Sea Water Temperature Loggers dataset
can be downloaded at daily average aggregations, which can reduce the size of
the download by many fold.


```r
aims_expose_attributes("temp_loggers")
#> $summary
#> [1] "summary-by-series"     "summary-by-deployment" "daily"                
#> 
#> $filters
#>  [1] "site"      "subsite"   "series"    "series_id" "parameter" "size"      "min_lat"  
#>  [8] "max_lat"   "min_lon"   "max_lon"   "from_date" "thru_date" "version"   "cursor"
```

In the [Sea Water Temperature Loggers][3] dataset, as demonstrated in our
[intro vignette][1], we have a convenience `summary` argument which facilitates 
learning more about what data is available. We can download the summary 
information for all sites using the main function called `aims_data`:

[3]: https://doi.org/10.25845/5b4eb0f9bb848
[4]: https://doi.org/10.25845/5c09bf93f315d


```r
sdata <- aims_data("temp_loggers", api_key = my_api_key,
                   summary = "summary-by-series")
head(sdata)
#>   site_id                    site subsite_id    subsite series_id     series
#> 1       1 Agincourt Reef Number 3       2687     AG3FL1      2687     AG3FL1
#> 2       1 Agincourt Reef Number 3      14276  AG3SL1old     14276  AG3SL1old
#> 3       3           Cleveland Bay       3007 CLEVAWSSL1      3007 CLEVAWSSL1
#> 4       3           Cleveland Bay       3069 CLEVAWSFL1      3069 CLEVAWSFL1
#> 5       4             Davies Reef       2629     DAVFL1      2629     DAVFL1
#> 6       4             Davies Reef       2630     DAVSL1      2630     DAVSL1
#>           parameter parameter_id time_coverage_start time_coverage_end      lat      lon
#> 1 Water Temperature            1          1996-03-30        2008-12-11 -15.9903 145.8212
#> 2 Water Temperature            1          1996-03-30        2011-07-21 -15.9905 145.8213
#> 3 Water Temperature            1          2004-05-13        2008-05-03 -19.1557 146.8813
#> 4 Water Temperature            1          2005-09-15        2005-12-22 -19.1557 146.8813
#> 5 Water Temperature            1          1997-08-26        2019-06-10 -18.8065 147.6688
#> 6 Water Temperature            1          1996-05-02        2021-03-29 -18.8060 147.6686
#>   depth uncal_obs cal_obs qc_obs
#> 1     0     23130  110480 110480
#> 2     5    114450  216794 216794
#> 3     7     11951   53231  53231
#> 4     1         0    4656   4656
#> 5     1    437544  566585 566585
#> 6     8    463146  589492 589437
```

`summary` should be set to either `summary-by-series` or `summary-by-deployment`
when the user wants an overview of the available data.


```r
ddata <- aims_data("temp_loggers", api_key = my_api_key,
                   summary = "summary-by-deployment")
head(ddata)
#>   deployment_id serial_num site_id             site subsite_id subsite series_id  series
#> 1          2616 SST-905084     860 Black Rocks Reef       2616  BLAFL1      2616  BLAFL1
#> 2          2612 SST-905053     856       Cattle Bay       2612   CBFL2      2612   CBFL2
#> 3          2613 SST-905054     857     Raine Island       2613 RAIDSL1      2613 RAIDSL1
#> 4          2620 SST-905189     863       Kelso Reef       2620  KELSL1      2620  KELSL1
#> 5          2619 SST-905091     863       Kelso Reef       2619  KELFL1      2619  KELFL1
#> 6          2621 SST-905068     865    Hayman Island       2621  HAYFL1      2621  HAYFL1
#>           parameter parameter_id time_coverage_start time_coverage_end      lat      lon
#> 1 Water Temperature            1          1996-07-20        1997-01-19 -16.2445 145.4872
#> 2 Water Temperature            1          1998-11-19        1999-03-02 -18.5719 146.4833
#> 3 Water Temperature            1          1996-11-28        1997-10-09 -11.5898 144.0309
#> 4 Water Temperature            1          1997-08-25        1998-04-18 -18.4221 146.9846
#> 5 Water Temperature            1          1997-08-25        1998-04-18 -18.4448 146.9933
#> 6 Water Temperature            1          1997-05-12        1998-04-21 -20.0413 148.8819
#>   depth uncal_obs cal_obs qc_obs
#> 1   0.1      8728    8728   8728
#> 2   4.0      4896    4896   4896
#> 3  20.0         0   15062  15062
#> 4   7.0     11278   11278  11278
#> 5   2.0     11278   11278  11278
#> 6   2.0     16460   16460  16460
```

Notice that `sdata` contains a lot of information, most of which is
related to site / series / parameter ID. Each row corresponds to a
unique series. The columns `time_coverage_start` and `time_coverage_end` are
probably one of the most valuable pieces of information. They provide the user
with the window of data collection for a particular series, which is probably
crucial to decide whether that particular series is of relevance to the
specific question in hand.

The benefits to choosing a data `series` (or the numeric equivalent,
`series_id`) is that it comes from one location and parameter type (here only
water temperature), making the data easy to plot. If we did not choose a
data series from the [Sea Water Temperature Loggers][4] dataset, we would have
to specify additional arguments to ensure the data is downloaded as expected.

Our values and filters might look like the following:

Variable  | Value                  | Description
----------|------------------------|-------------------------------------------------------
series_id | 2687                   | Found [here][6], Agincourt Reef Number 3
from_date | "2005-01-01"           | We want to start charting on 1/1/2005
thru_date | "2005-01-10"           | We are plotting 10 days of data

[5]: https://open-aims.github.io/data-platform
[6]: https://apps.aims.gov.au/metadata/view/4a12a8c0-c573-11dc-b99b-00008a07204e

## Query and Plot Dataset

After deciding on query parameters, we plug the series id into a `aims_data` function:


```r
agincourt <- aims_data("temp_loggers", api_key = my_api_key,
                       filters = list(series_id = 2687,
                                      from_date = "2005-01-01",
                                      thru_date = "2005-01-10"))
```

We can check that the query filters worked:


```r
range(agincourt$time)
#> [1] "2005-01-01 UTC" "2005-01-10 UTC"
```

We can then visualise where in Australia that data is placed:


```r
plot(agincourt, ptype = "map")
```

<img src="vignette-fig-tlfa-1.png" title="plot of chunk tlfa" alt="plot of chunk tlfa" width="100%" />

We can also visually compare multiple series at once. For instance, let's
compare the air temperature data from Davies Reef and Bramble Cay for the
same period of time:


```r
target_series <- c("Agincourt" = 2687, "Cleveland Bay" = 3007)
aims_data_per_series <- function(series_number, my_api_key, ...) {
  aims_data("temp_loggers", api_key = my_api_key,
            filters = list(series_id = series_number, ...))
}
results <- purrr::map(target_series, aims_data_per_series,
                      my_api_key = my_api_key,
                      from_date = "2005-01-01",
                      thru_date = "2005-01-10")
sst_data <- purrr::map_dfr(results, rbind)
plot(sst_data, ptype = "time_series")
```

<img src="vignette-fig-tlfb-1.png" title="plot of chunk tlfb" alt="plot of chunk tlfb" width="100%" />

One could also download data for a particular time of day throughout
the year, e.g. for Davies Reef at 1 m of depth (`series_id` is 2629):


```r
days <- seq(as.Date("2005-01-01"), as.Date("2005-12-31"), by = "month")
out <- numeric(length = length(days))
for (i in seq_along(days)) {
  hour_in <- paste0(days[i], "T06:00:00")
  hour_out <- paste0(days[i], "T12:00:00")
  df <- aims_data("temp_loggers", api_key = my_api_key,
                  filters = list(series_id = 2629, from_date = hour_in,
                                 thru_date = hour_out))
  out[i] <- mean(df$qc_val)
}

ggplot(data = data.frame(date = days, temps = out)) +
  geom_line(mapping = aes(x = date, y = temps)) +
  labs(x = "Date",
       y = "Water temperature (˚C)",
       title = "Davies Reef @ 1 m (2005)",
       subtitle = "mean 6 A.M. – 12 P.M.") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "bottom")
```

<img src="vignette-fig-tlfc-1.png" title="plot of chunk tlfc" alt="plot of chunk tlfc" width="100%" />

Or simply plot the daily aggregated averages:


```r
df <- aims_data("temp_loggers", api_key = my_api_key, summary = "daily",
                filters = list(series_id = 2629, from_date = "2005-01-01",
                               thru_date = "2005-12-31"))
plot(df, ptype = "time_series", pars = c("Water Temperature"))
```

<img src="vignette-fig-unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" width="100%" />

## Bibliography


```r
purrr::map_chr(results, aims_citation) %>%
  unlist %>%
  unname
#> [1] "Australian Institute of Marine Science (AIMS). 2017, AIMS Sea Temperature Observing System (AIMS Temperature Logger Program), Time period:2005-01-01 to 2005-01-10. https://doi.org/10.25845/5b4eb0f9bb848, accessed 29 Oct 2021."
#> [2] "Australian Institute of Marine Science (AIMS). 2017, AIMS Sea Temperature Observing System (AIMS Temperature Logger Program), Time period:2005-01-01 to 2005-01-10. https://doi.org/10.25845/5b4eb0f9bb848, accessed 29 Oct 2021."
```
---
title: "AIMS Weather Station time series datasets"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{AIMS Weather Station time series datasets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



Please check our [intro vignette][1] first to implement the installation
requirements, and to learn the general approach to navigating the different
datasets. This vignette assumes you have obtained an
[AIMS Data Platform API Key][2].

[1]: https://ropensci.github.io/dataaimsr/articles/navigating.html
[2]: https://open-aims.github.io/data-platform/key-request

As per the installation instructions, we strongly suggest that you hide your
API Key permanently in your `.Renviron` file and set the object `my_api_key` to
`NULL` in the chunk below. You can read more about why that is important
[here](https://CRAN.R-project.org/package=httr/vignettes/secrets.html).


```r
# set my_api_key to NULL after successfully placing it in .Renviron
my_api_key <- NULL
```

Let's start by loading the packages needed for this vignette:


```r
library(purrr)
library(dataaimsr)
library(ggplot2)
```

## Discover a dataset

The [AIMS Weather Station][3] dataset consists of a series of *"parameters"* 
that are measured at multiple sites. Those could be, for instance, minimum 
wind speed, peak wave period, wind direction, water pressure, dissolved 
oxygen, chlorophyll concentration, etc. Not all parameters are measured for 
all sites and at all times. Some loggers are still actively collecting data, 
others have been discontinued.

Essentially, for the user who has limited knowledge about where the data are,
and of what they are consisted, they would need to do some prior exploration 
to learn more about what can downloaded. Suppose the goal is to download all
time-series from a particular site. The general procedure would be:

1. Examine documentation and establish query filters
2. Perform data download using `aims_data`
3. Use the R package `ggplot2` to create an exploratory time-series chart

For all datasets, a list of available filters can be retrieved with the 
function `aims_expose_attributes`. Knowing the filters is important because
some time series are quite extensive, with parameters being measured at very
high frequency (e.g. every 5 minutes), so downloading the dataset for an
entire year or more my take quite some time (it's possible though if that is
the true goal of the user).


```r
aims_expose_attributes("weather")
#> $summary
#> [1] NA
#> 
#> $filters
#>  [1] "site"      "subsite"   "series"    "series_id" "parameter" "size"      "min_lat"  
#>  [8] "max_lat"   "min_lon"   "max_lon"   "from_date" "thru_date" "version"   "cursor"
```

In the [Sea Water Temperature Loggers][4] dataset, as demonstrated in our
[intro vignette][1], we have a convenience `summary` method which facilitates
learning more about what data is available, or downloading daily aggregated
average data. The back-end for these is currently
being developed for the AIMS Weather Station as well. In the meantime, to
explore the [AIMS Weather Station][3] dataset, we use the function
`aims_filter_values`. This function takes a target dataset and a given filter,
and returns all the available information regarding the filter. We recommend
exploring the `series` filter---a series is a continuing time-series, i.e. a
collection of deployments measuring the same parameter (e.g. air temperature,
air pressure, chlorophyll) at the same subsite. So, for a given site and
parameter, there might exist multiple subsites and therefore series, in which
case they are most likely distinguishable by depth.

[3]: https://doi.org/10.25845/5c09bf93f315d
[4]: https://doi.org/10.25845/5b4eb0f9bb848


```r
head(aims_filter_values("weather", filter_name = "series"))
#>   series_id                                                                 series
#> 1    104918        Myrmidon Reef Weather Station Wind Speed (scalar avg b 10 min) 
#> 2    100686                            Saibai Island Weather Station Hail Duration
#> 3       266 Orpheus Island Relay Pole 3 Wind Direction (Vector Average 30 Minutes)
#> 4      2639 Hardy Reef Weather Station Wind Direction (Vector Standard 10 Minutes)
#> 5     10243                           Raine Island Weather Station Air Temperature
#> 6       258             Orpheus Island Relay Pole 3 Wind Speed (Scalar avg 10 min)
```

The benefits to choosing a data `series` is that it comes from one location
and parameter type, making the data easy to plot. If we did not choose a
data series from the [AIMS Weather Station][3] dataset, we would have to 
specify additional arguments to ensure the data is downloaded as expected.

Our values and filters might look like the following:

Variable  | Value                  | Description
----------|------------------------|-------------------------------------------------------
series_id | 64                     | Found [here][6], Davies Reef Air Temperature data series
from_date | "2018-01-01"           | We want to start charting on 1/1/2018
thru_date | "2018-01-10"           | We are plotting 10 days of data

[5]: https://open-aims.github.io/data-platform
[6]: https://apps.aims.gov.au/metadata/view/5fc91100-4ade-11dc-8f56-00008a07204e

## Query and Plot Dataset

After deciding on query parameters, we plug the series id into a `aims_data` function:


```r
davies <- aims_data("weather", api_key = my_api_key,
                    filters = list(series_id = 64,
                                   from_date = "2018-01-01",
                                   thru_date = "2018-01-10"))
```

We can even visually compare multiple series at once. For instance, let's
compare the air temperature data from Davies Reef and Bramble Cay for the
same period of time:


```r
target_series <- c("Davies Reef" = 64, "Bramble Cay" = 87929)
aims_data_per_series <- function(series_number, my_api_key, ...) {
  aims_data("weather", api_key = my_api_key,
            filters = list(series_id = series_number, ...))
}
results <- purrr::map(target_series, aims_data_per_series,
                      my_api_key = my_api_key,
                      from_date = "2018-01-01",
                      thru_date = "2018-01-10")
weather_data <- purrr::map_dfr(results, rbind)
plot(weather_data, ptype = "time_series")
```

<img src="vignette-fig-wfa-1.png" title="plot of chunk wfa" alt="plot of chunk wfa" width="100%" />

One could also download data for a particular time of day throughout
the year, e.g. for Heron Island Relay Pole 8 at 5.4 m of depth (series
10394):


```r
days <- seq(as.Date("2018-01-01"), as.Date("2018-12-31"), by = "month")
out <- numeric(length = length(days))
for (i in seq_along(days)) {
  hour_in <- paste0(days[i], "T06:00:00")
  hour_out <- paste0(days[i], "T12:00:00")
  df <- aims_data("weather",
                  api_key = my_api_key,
                  filters = list(series_id = 10394,
                                 from_date = hour_in,
                                 thru_date = hour_out))
  out[i] <- mean(df$qc_val)
}
ggplot(data = data.frame(date = days, temps = out)) +
  geom_line(mapping = aes(x = date, y = temps)) +
  labs(x = "Date",
       y = "Air temperature (˚C)",
       title = "Heron Island (2018)",
       subtitle = "mean 6 A.M. – 12 P.M.") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "bottom")
```

<img src="vignette-fig-wfb-1.png" title="plot of chunk wfb" alt="plot of chunk wfb" width="100%" />

## Bibliography


```r
purrr::map_chr(results, aims_citation) %>%
  unlist %>%
  unname
#> [1] "Australian Institute of Marine Science (AIMS). 2009, Australian Institute of Marine Science Automatic Weather Stations, https://doi.org/10.25845/5c09bf93f315d, accessed 29 October 2021.  Time period: 2018-01-01 to 2018-01-10.  Series: Davies Reef Weather Station Air Temperature"
#> [2] "Australian Institute of Marine Science (AIMS). 2009, Australian Institute of Marine Science Automatic Weather Stations, https://doi.org/10.25845/5c09bf93f315d, accessed 29 October 2021.  Time period: 2018-01-01 to 2018-01-10.  Series: Bramble Cay Weather Station Air Temperature"
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataaimsr-package.R
\docType{package}
\name{dataaimsr-package}
\alias{dataaimsr-package}
\alias{dataaimsr}
\title{The 'dataaimsr' package.}
\description{
dataaimsr is the Australian Institute of Marine Science (AIMS)
Data Platform R package, and provides the user with easy access to datasets
from the AIMS Data Platform API. Please see ?aims_data for more details.
}
\references{
Australian Institute of Marine Science (AIMS). (2017). AIMS Sea Water
Temperature Observing System (AIMS Temperature Logger Program)
https://doi.org/10.25845/5b4eb0f9bb848

Australian Institute of Marine Science (AIMS). (2017). Northern Australia
Automated Marine Weather and Oceanographic Stations,
https://doi.org/10.25845/5c09bf93f315d
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{process_request}
\alias{process_request}
\title{Format \code{\link{json_results}} output}
\usage{
process_request(dt_req, next_page = FALSE, ...)
}
\arguments{
\item{dt_req}{An URL \code{\link[httr]{GET}} output}

\item{next_page}{Logical. Is this a multi-url request?}

\item{...}{Additional arguments to be passed to internal function
\code{\link{update_format}}}
}
\value{
\code{aims_data} returns a \code{\link[base]{data.frame}} of class
\code{\link{aimsdf}}.

If \code{summary \%in\% c("summary-by-series", "summary-by-deployment")},
the output shows the summary information for the target dataset (i.e.
weather or temperature loggers)
(NB: currently, \code{summary} only works for the temperature logger
database). If \code{summary} is \emph{not} passed as an additional argument, then
the output contains \strong{raw} monitoring data. If \code{summary = "daily"},
then the output contains \strong{mean daily aggregated} monitoring data.
The output also contains five attributes (empty strings if
\code{summary} is passed as an additional argument):
\itemize{
\item{\code{metadata}}{a \href{https://www.doi.org/}{DOI} link
containing the metadata record for the data series.}
\item{\code{citation}}{the citation information for the particular
dataset.}
\item{\code{parameters}}{The measured parameters comprised in the
output.}
\item{\code{type}}{The type of dataset. Either "monitoring" if
\code{summary} is not specified, "monitoring (daily aggregation)" if
\code{summary = "daily"}, or a "summary-by-" otherwise.}
\item{\code{target}}{The input target.}
}
}
\description{
Wrapper function
}
\details{
This function checks for errors in \code{dt_req}
data request and processes result via
\code{\link{json_results}}.
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aims_data.R
\name{aims_data}
\alias{aims_data}
\title{Request data via the AIMS Data Platform API}
\usage{
aims_data(target, filters = NULL, summary = NA, ...)
}
\arguments{
\item{target}{A \code{\link[base]{character}} vector of length 1 specifying
the dataset. Only \code{weather} or \code{temp_loggers} are currently
allowed.}

\item{filters}{A \code{\link[base]{list}} containing a set of
filters for the data query (see Details).}

\item{summary}{Should summary tables (\code{"summary-by-series"} or
\code{"summary-by-deployment"}) or daily aggregated data ("daily") be
returned instead of full data (see Details)?}

\item{...}{Currently unused. Additional arguments to be passed to
non-exported internal functions.}
}
\value{
\code{aims_data} returns a \code{\link[base]{data.frame}} of class
\code{\link{aimsdf}}.

If \code{summary \%in\% c("summary-by-series", "summary-by-deployment")},
the output shows the summary information for the target dataset (i.e.
weather or temperature loggers)
(NB: currently, \code{summary} only works for the temperature logger
database). If \code{summary} is \emph{not} passed as an additional argument, then
the output contains \strong{raw} monitoring data. If \code{summary = "daily"},
then the output contains \strong{mean daily aggregated} monitoring data.
The output also contains five attributes (empty strings if
\code{summary} is passed as an additional argument):
\itemize{
\item{\code{metadata}}{a \href{https://www.doi.org/}{DOI} link
containing the metadata record for the data series.}
\item{\code{citation}}{the citation information for the particular
dataset.}
\item{\code{parameters}}{The measured parameters comprised in the
output.}
\item{\code{type}}{The type of dataset. Either "monitoring" if
\code{summary} is not specified, "monitoring (daily aggregation)" if
\code{summary = "daily"}, or a "summary-by-" otherwise.}
\item{\code{target}}{The input target.}
}
}
\description{
A function that communicates with the
the \href{https://open-aims.github.io/data-platform/}{AIMS Data Platform}
via the AIMS Data Platform API
}
\details{
The AIMS Data Platform R Client provides easy access to
data sets for R applications to the
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}.
The AIMS Data Platform requires an API Key for requests, which can
be obtained at this
\href{https://open-aims.github.io/data-platform/key-request}{link}.
It is preferred that API Keys are not stored in code. We recommend
storing the environment variable \code{AIMS_DATAPLATFORM_API_KEY}
permanently under the user's \code{.Renviron} file in order to load
the API Key automatically.

There are two types of data currently available through the
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}:
\href{https://weather.aims.gov.au/#/overview}{Weather} and
\href{https://tinyurl.com/h93mcojk}{Sea Water Temperature Loggers}.
They are searched internally via unique DOI identifiers.
Only one data type at a time can be passed to the argument \code{target}.

A list of arguments for \code{filters} can be exposed for both
\href{https://weather.aims.gov.au/#/overview}{Weather} and
\href{https://weather.aims.gov.au/#/overview}{Sea Water Temperature Loggers}
using function \code{\link{aims_expose_attributes}}.

Note that at present the user can inspect the range of dates for
the temperature loggers data only (see usage of argument \code{summary} in
the examples below). For that, the argument \code{summary} must be either
the string \code{"summary-by-series"} or \code{"summary-by-deployment"}.
In those cases, time filters will be ignored.

Details about available dates for each dataset and time series can be
accessed via Metadata on
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}.
We raise this caveat here because these time boundaries are very important;
data are collected at very small time intervals, a window of just a few days
can yield very large datasets. The query will return and error
if it reaches the system's memory capacity.

For that same reason, from version 1.1.0 onwards, we are offering the
possibility of downloading a mean daily aggregated version. For that, the
user must set \code{summary = "daily"}. In this particular case, query filter
will be taken into account.
}
\examples{
\dontrun{
library(dataaimsr)
# assumes that user already has API key saved to
# .Renviron

# start downloads:
# 1. downloads weather data from
# site Yongala
# within a defined date range
wdf_a <- aims_data("weather", api_key = NULL,
                   filters = list(site = "Yongala",
                                  from_date = "2018-01-01",
                                  thru_date = "2018-01-02"))

# 2. downloads weather data from all sites
# under series_id 64 from Davies Reef
# within a defined date range
wdf_b <- aims_data("weather", api_key = NULL,
                   filters = list(series_id = 64,
                                  from_date = "1991-10-18",
                                  thru_date = "1991-10-19"))
head(wdf_b)
range(wdf_b$time)

# 3. downloads weather data from all sites
# under series_id 64 from Davies Reef
# within defined date AND time range
wdf_c <- aims_data("weather", api_key = NULL,
                   filters = list(series_id = 64,
                                  from_date = "1991-10-18T06:00:00",
                                  thru_date = "1991-10-18T12:00:00"))
head(wdf_c)
range(wdf_c$time)

# 4. downloads all parameters from all sites
# within a defined date range
wdf_d <- aims_data("weather", api_key = NULL,
                   filters = list(from_date = "2003-01-01",
                                  thru_date = "2003-01-02"))
# note that there are multiple sites and series
# so in this case, because we did not specify a specific
# parameter, series within sites could differ by both
# parameter and depth
head(wdf_d)
unique(wdf_d[, c("site", "series_id", "series")])
unique(wdf_d$parameter)
range(wdf_d$time)

# 5. downloads chlorophyll from all sites
# within a defined date range
wdf_e <- aims_data("weather", api_key = NULL,
                   filters = list(parameter = "Chlorophyll",
                                  from_date = "2018-01-01",
                                  thru_date = "2018-01-02"))
# note again that there are multiple sites and series
# however in this case because we did specify a specific
# parameter, series within sites differ by depth only
head(wdf_e)
unique(wdf_e[, c("site", "series_id", "series", "depth")])
unique(wdf_e$parameter)
range(wdf_e$time)

# 6. downloads temperature data
# summarised by series
sdf_a <- aims_data("temp_loggers", api_key = NULL,
                   summary = "summary-by-series")
head(sdf_a)
dim(sdf_a)

# 7. downloads temperature data
# summarised by series
# for all sites that contain data
# within a defined date range
sdf_b <- aims_data("temp_loggers", api_key = NULL,
                   summary = "summary-by-series",
                   filters = list("from_date" = "2018-01-01",
                                  "thru_date" = "2018-12-31"))
head(sdf_b)
dim(sdf_b) # a subset of sdf_a

# 8. downloads temperature data
# summarised by deployment
sdf_c <- aims_data("temp_loggers", api_key = NULL,
                   summary = "summary-by-deployment")
head(sdf_c)
dim(sdf_c)

# 9. downloads temperature data
# within a defined date range, averaged by day
sdf_d <- aims_data("temp_loggers", api_key = NULL, summary = "daily",
                   filters = list(series = "DAVFL1",
                                  from_date = "2018-01-01",
                                  thru_date = "2018-01-10"))
# note again that there are multiple sites and series
# however in this case because we did specify a specific
# parameter, series within sites differ by depth only
head(sdf_d)
unique(sdf_d[, c("site", "series_id", "series", "depth")])
unique(sdf_d$parameter)
range(sdf_d$time)
}

}
\seealso{
\code{\link{aims_citation}}, \code{\link{aims_metadata}},
\code{\link{aims_parameters}}
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{make_base_end_pt}
\alias{make_base_end_pt}
\title{Expose available query filters}
\usage{
make_base_end_pt(doi, aims_version = NA)
}
\arguments{
\item{doi}{A \href{https://www.doi.org/}{Digital Object Identifier}
for a chosen
\href{https://open-aims.github.io/data-platform/}{AIMS data series}}

\item{aims_version}{A \code{\link[base]{character}} string
defining the version of database. Must be "/v1.0" or "-v2.0".
If none is provided, then "-v2.0" (the most recent) is used.}
}
\description{
Expose available query filters which are allowed to be parsed either
via argument \code{summary} or \code{filters} in \code{\link{aims_data}}
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aims_parameters.R
\name{aims_parameters}
\alias{aims_parameters}
\title{Extracts parameters attribute from object of class aimsdf}
\usage{
aims_parameters(df_)
}
\arguments{
\item{df_}{A data.frame of class \code{\link{aimsdf}} created by
function \code{\link{aims_data}}}
}
\value{
A \code{\link[base]{character}} vector.
}
\description{
Extracts parameters attribute from object of class aimsdf
}
\details{
This function retrieves the parameters attribute from an
\code{\link{aimsdf}} object. If the input \code{\link{aimsdf}} object is
a summary data.frame (see ?\code{\link{aims_data}}), then output will be
an empty string.
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aimsdf-methods.R
\name{plot.aimsdf}
\alias{plot.aimsdf}
\title{plot.aimsdf}
\usage{
\method{plot}{aimsdf}(x, ..., ptype, pars)
}
\arguments{
\item{x}{An object of class \code{\link{aimsdf}} as
returned by \code{\link{aims_data}}.}

\item{...}{Not used.}

\item{ptype}{Type of plot. Can either be "time_series" or "map".}

\item{pars}{Which parameters to plot? Only relevant if ptype is
"time_series"}
}
\value{
An object of class \code{\link[ggplot2]{ggplot}}.
}
\description{
Plotting options for aimsdf objects
}
\details{
Currently plots cannot be customised. Summary datasets can only
be represented by maps.
}
\examples{
\dontrun{
library(dataaimsr)
wdf <- aims_data("weather", api_key = NULL,
                 filters = list(site = "Yongala",
                                from_date = "2018-01-01",
                                thru_date = "2018-01-02"))
plot(wdf, ptype = "map")
plot(wdf, ptype = "time_series")
# summary-by- datasets can only return maps
sdf <- aims_data("temp_loggers", api_key = NULL,
                 summary = "summary-by-deployment")
plot(sdf, ptype = "map")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{handle_error}
\alias{handle_error}
\title{\code{\link[httr]{GET}} error handler}
\usage{
handle_error(dt_req)
}
\arguments{
\item{dt_req}{An URL \code{\link[httr]{GET}} output}
}
\value{
A \code{\link[base]{character}} vector conveying the error message.
}
\description{
Displays error status
}
\details{
This function retrieves the status and content of \code{dt_req}
via the \pkg{httr} package.
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/page_data.R
\name{page_data}
\alias{page_data}
\title{Request data via the AIMS Data Platform API}
\usage{
page_data(
  doi,
  filters = NULL,
  api_key = NULL,
  summary = NA,
  aims_version = NA,
  verbose = FALSE
)
}
\arguments{
\item{doi}{A \href{https://www.doi.org/}{Digital Object Identifier}
for a chosen
\href{https://open-aims.github.io/data-platform/}{AIMS data series}}

\item{filters}{A \code{\link[base]{list}} containing a set of
filters for the data query (see Details).}

\item{api_key}{An AIMS Data Platform
\href{https://open-aims.github.io/data-platform/key-request}{API Key}}

\item{summary}{Should summary tables (\code{"summary-by-series"} or
\code{"summary-by-deployment"}) or daily aggregated data ("daily") be
returned instead of full data (see Details)?}

\item{aims_version}{A \code{\link[base]{character}} string
defining the version of database. Must be "/v1.0" or "-v2.0".
If none is provided, then "-v2.0" (the most recent) is used.}

\item{verbose}{Should links be printed to screen? Used for debugging only}
}
\value{
\code{aims_data} returns a \code{\link[base]{data.frame}} of class
\code{\link{aimsdf}}.

If \code{summary \%in\% c("summary-by-series", "summary-by-deployment")},
the output shows the summary information for the target dataset (i.e.
weather or temperature loggers)
(NB: currently, \code{summary} only works for the temperature logger
database). If \code{summary} is \emph{not} passed as an additional argument, then
the output contains \strong{raw} monitoring data. If \code{summary = "daily"},
then the output contains \strong{mean daily aggregated} monitoring data.
The output also contains five attributes (empty strings if
\code{summary} is passed as an additional argument):
\itemize{
\item{\code{metadata}}{a \href{https://www.doi.org/}{DOI} link
containing the metadata record for the data series.}
\item{\code{citation}}{the citation information for the particular
dataset.}
\item{\code{parameters}}{The measured parameters comprised in the
output.}
\item{\code{type}}{The type of dataset. Either "monitoring" if
\code{summary} is not specified, "monitoring (daily aggregation)" if
\code{summary = "daily"}, or a "summary-by-" otherwise.}
\item{\code{target}}{The input target.}
}
}
\description{
A function that communicates with the
the \href{https://open-aims.github.io/data-platform/}{AIMS Data Platform}
via the AIMS Data Platform API
}
\details{
The AIMS Data Platform R Client provides easy access to
data sets for R applications to the
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}.
The AIMS Data Platform requires an API Key for requests, which can
be obtained at this
\href{https://open-aims.github.io/data-platform/key-request}{link}.
It is preferred that API Keys are not stored in code. We recommend
storing the environment variable \code{AIMS_DATAPLATFORM_API_KEY}
permanently under the user's \code{.Renviron} file in order to load
the API Key automatically.

There are two types of data currently available through the
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}:
\href{https://weather.aims.gov.au/#/overview}{Weather} and
\href{https://tinyurl.com/h93mcojk}{Sea Water Temperature Loggers}.
They are searched internally via unique DOI identifiers.
Only one data type at a time can be passed to the argument \code{target}.

A list of arguments for \code{filters} can be exposed for both
\href{https://weather.aims.gov.au/#/overview}{Weather} and
\href{https://weather.aims.gov.au/#/overview}{Sea Water Temperature Loggers}
using function \code{\link{aims_expose_attributes}}.

Note that at present the user can inspect the range of dates for
the temperature loggers data only (see usage of argument \code{summary} in
the examples below). For that, the argument \code{summary} must be either
the string \code{"summary-by-series"} or \code{"summary-by-deployment"}.
In those cases, time filters will be ignored.

Details about available dates for each dataset and time series can be
accessed via Metadata on
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}.
We raise this caveat here because these time boundaries are very important;
data are collected at very small time intervals, a window of just a few days
can yield very large datasets. The query will return and error
if it reaches the system's memory capacity.

For that same reason, from version 1.1.0 onwards, we are offering the
possibility of downloading a mean daily aggregated version. For that, the
user must set \code{summary = "daily"}. In this particular case, query filter
will be taken into account.
}
\seealso{
\code{\link{aims_expose_attributes}},
\code{\link{aims_filter_values}}, \code{\link{aims_data}}
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{find_api_key}
\alias{find_api_key}
\title{AIMS API Key retriever}
\usage{
find_api_key(api_key)
}
\arguments{
\item{api_key}{An AIMS Data Platform
\href{https://open-aims.github.io/data-platform/key-request}{API Key}}
}
\value{
Either a \code{\link[base]{character}} vector API Key found
in .Renviron or, if missing entirely, an error message.
}
\description{
This function tries to search for an API Key
}
\details{
The AIMS Data Platform R Client provides easy access to
data sets for R applications to the
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}.
The AIMS Data Platform requires an API Key for requests, which can
be obtained at this
\href{https://open-aims.github.io/data-platform/key-request}{link}.
It is preferred that API Keys are not stored in code. We recommend
storing the environment variable \code{AIMS_DATAPLATFORM_API_KEY}
permanently under the user's \code{.Renviron} file in order to load
the API Key automatically.

There are two types of data currently available through the
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}:
\href{https://weather.aims.gov.au/#/overview}{Weather} and
\href{https://tinyurl.com/h93mcojk}{Sea Water Temperature Loggers}.
They are searched internally via unique DOI identifiers.
Only one data type at a time can be passed to the argument \code{target}.

A list of arguments for \code{filters} can be exposed for both
\href{https://weather.aims.gov.au/#/overview}{Weather} and
\href{https://weather.aims.gov.au/#/overview}{Sea Water Temperature Loggers}
using function \code{\link{aims_expose_attributes}}.

Note that at present the user can inspect the range of dates for
the temperature loggers data only (see usage of argument \code{summary} in
the examples below). For that, the argument \code{summary} must be either
the string \code{"summary-by-series"} or \code{"summary-by-deployment"}.
In those cases, time filters will be ignored.

Details about available dates for each dataset and time series can be
accessed via Metadata on
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}.
We raise this caveat here because these time boundaries are very important;
data are collected at very small time intervals, a window of just a few days
can yield very large datasets. The query will return and error
if it reaches the system's memory capacity.

For that same reason, from version 1.1.0 onwards, we are offering the
possibility of downloading a mean daily aggregated version. For that, the
user must set \code{summary = "daily"}. In this particular case, query filter
will be taken into account.
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{data_doi}
\alias{data_doi}
\title{AIMS Dataset DOI retriever}
\usage{
data_doi(target)
}
\arguments{
\item{target}{A \code{\link[base]{character}} vector of length 1 specifying
the dataset. Only \code{weather} or \code{temp_loggers} are currently
allowed.}
}
\value{
A \code{\link[base]{character}} vector
containing the dataset DOI string.
}
\description{
Returns DOI for a given dataset
}
\examples{
\dontrun{
library(dataaimsr)
weather_doi <- data_doi("weather")
ssts_doi <- data_doi("temp_loggers")
}
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aims_filter_values.R
\name{aims_filter_values}
\alias{aims_filter_values}
\title{Retrieve vector of existing filter values}
\usage{
aims_filter_values(target, filter_name)
}
\arguments{
\item{target}{A \code{\link[base]{character}} vector of length 1 specifying
the dataset. Only \code{weather} or \code{temp_loggers} are currently
allowed.}

\item{filter_name}{A \code{\link[base]{character}} string containing the
name of the filter. Must be "site", "subsite", "series", or
"parameter". See details.}
}
\value{
Either a \code{\link[base]{data.frame}} if
\code{filter_name = "series"}, else a \code{\link[base]{character}}
vector.
}
\description{
This is a utility function which allows to user
to query about the existing possibilities of
a given filter name
}
\details{
For a full description of each valid filter_name see
?\code{\link{aims_expose_attributes}}. In the temperature logger dataset,
"subsite" is equivalent to "series"; moreover, note that there is only one
parameter being measured (i.e. water temperature), so the "parameter" filter
contains one single value.
}
\examples{
\dontrun{
library(dataaimsr)
aims_filter_values("weather", filter_name = "site")
aims_filter_values("temp_loggers", filter_name = "subsite")
}

}
\seealso{
\code{\link{aims_data}}, \code{\link{aims_expose_attributes}}
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{capitalise}
\alias{capitalise}
\title{capitalise}
\usage{
capitalise(x)
}
\arguments{
\item{x}{A character}
}
\description{
Internal
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aims_expose_attributes.R
\name{aims_expose_attributes}
\alias{aims_expose_attributes}
\title{Expose available query filters}
\usage{
aims_expose_attributes(target)
}
\arguments{
\item{target}{A \code{\link[base]{character}} vector of length 1 specifying
the dataset. Only \code{weather} or \code{temp_loggers} are currently
allowed.}
}
\value{
A \code{\link[base]{list}} of two \code{\link[base]{character}}
vectors: one detailing summary modes, another detailing filters.
}
\description{
Expose available query filters which are allowed to be parsed either
via argument \code{summary} or \code{filters} in \code{\link{aims_data}}
}
\details{
Use this function to learn which summary modes and
filters are allowed.

We are working on implementing summary visualisation methods for weather
station data. So, for the moment, the options below are only available
for temperature logger data. Three options are available:

\itemize{
\item{summary-by-series}{Expose summary for all available series;
a series is a continuing time-series, i.e. a collection of
deployments measuring the same parameter at the same site.
For temperature loggers, series is synonymous with sub-site.
For weather stations, it is the combination of sub-site and
parameter.}
\item{summary-by-deployment}{Expose summary for all available
deployments.}
\item{daily}{Return mean daily aggregated monitoring data .}
}

We offer a list of valid filter names:

\itemize{
\item{site}{Filter by a particular site.}
\item{subsite}{Filter by a particular subsite.}
\item{series}{Filter by a particular series.}
\item{series_id}{A unique identifier for the series - it should be
unique within a dataset. An alternative to looking up a series by name.}
\item{parameter}{Parameter of interest. Only relevant for
weather station data because temperature logger is always water
temperature.}
\item{min_lat}{Minimum latitude; used to filter by a
lat-lon box.}
\item{max_lat}{Maximum latitude; used to filter by a
lat-lon box.}
\item{min_lon}{Minimum longitude; used to filter by a
lat-lon box.}
\item{max_lon}{Maximum longitude; used to filter by a
lat-lon box.}
\item{from_date}{Filter from time (string of format
YYYY-MM-DD).}
\item{thru_date}{Filter until time (string of format
YYYY-MM-DD).}
}

Some additional options for the actual download, which should be passed as
additional arguments to the function, are:
\itemize{
\item{size}{Set a page size for large queries
(only for the \code{data} and \code{data-no-key} endpoints).}
\item{cursor}{Used for pagination on / data").}
\item{version}{Request the data as recorded at a particular time
(a version history).}
}
}
\examples{
\dontrun{
library(dataaimsr)
aims_expose_attributes("weather")
aims_expose_attributes("temp_loggers")
}

}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/page_data.R
\name{next_page_data}
\alias{next_page_data}
\title{Further data requests via the AIMS Data Platform API}
\usage{
next_page_data(url, api_key = NULL, ...)
}
\arguments{
\item{url}{A data retrieval URL}

\item{api_key}{An AIMS Data Platform
\href{https://open-aims.github.io/data-platform/key-request}{API Key}}

\item{...}{Additional arguments to be passed to internal function
\code{\link{update_format}}}
}
\value{
\code{aims_data} returns a \code{\link[base]{data.frame}} of class
\code{\link{aimsdf}}.

If \code{summary \%in\% c("summary-by-series", "summary-by-deployment")},
the output shows the summary information for the target dataset (i.e.
weather or temperature loggers)
(NB: currently, \code{summary} only works for the temperature logger
database). If \code{summary} is \emph{not} passed as an additional argument, then
the output contains \strong{raw} monitoring data. If \code{summary = "daily"},
then the output contains \strong{mean daily aggregated} monitoring data.
The output also contains five attributes (empty strings if
\code{summary} is passed as an additional argument):
\itemize{
\item{\code{metadata}}{a \href{https://www.doi.org/}{DOI} link
containing the metadata record for the data series.}
\item{\code{citation}}{the citation information for the particular
dataset.}
\item{\code{parameters}}{The measured parameters comprised in the
output.}
\item{\code{type}}{The type of dataset. Either "monitoring" if
\code{summary} is not specified, "monitoring (daily aggregation)" if
\code{summary = "daily"}, or a "summary-by-" otherwise.}
\item{\code{target}}{The input target.}
}
}
\description{
Similar to \code{\link{page_data}}, but for cases #' where there are
multiple URLs for data retrieval
}
\details{
The AIMS Data Platform R Client provides easy access to
data sets for R applications to the
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}.
The AIMS Data Platform requires an API Key for requests, which can
be obtained at this
\href{https://open-aims.github.io/data-platform/key-request}{link}.
It is preferred that API Keys are not stored in code. We recommend
storing the environment variable \code{AIMS_DATAPLATFORM_API_KEY}
permanently under the user's \code{.Renviron} file in order to load
the API Key automatically.

There are two types of data currently available through the
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}:
\href{https://weather.aims.gov.au/#/overview}{Weather} and
\href{https://tinyurl.com/h93mcojk}{Sea Water Temperature Loggers}.
They are searched internally via unique DOI identifiers.
Only one data type at a time can be passed to the argument \code{target}.

A list of arguments for \code{filters} can be exposed for both
\href{https://weather.aims.gov.au/#/overview}{Weather} and
\href{https://weather.aims.gov.au/#/overview}{Sea Water Temperature Loggers}
using function \code{\link{aims_expose_attributes}}.

Note that at present the user can inspect the range of dates for
the temperature loggers data only (see usage of argument \code{summary} in
the examples below). For that, the argument \code{summary} must be either
the string \code{"summary-by-series"} or \code{"summary-by-deployment"}.
In those cases, time filters will be ignored.

Details about available dates for each dataset and time series can be
accessed via Metadata on
\href{https://open-aims.github.io/data-platform/}{AIMS Data Platform API}.
We raise this caveat here because these time boundaries are very important;
data are collected at very small time intervals, a window of just a few days
can yield very large datasets. The query will return and error
if it reaches the system's memory capacity.

For that same reason, from version 1.1.0 onwards, we are offering the
possibility of downloading a mean daily aggregated version. For that, the
user must set \code{summary = "daily"}. In this particular case, query filter
will be taken into account.
}
\seealso{
\code{\link{aims_filter_values}}, \code{\link{page_data}},
\code{\link{aims_data}}
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aimsdf-methods.R
\name{summary.aimsdf}
\alias{summary.aimsdf}
\title{summary.aimsdf}
\usage{
\method{summary}{aimsdf}(object, ...)
}
\arguments{
\item{object}{An object of class \code{\link{aimsdf}} as
returned by \code{\link{aims_data}}.}

\item{...}{Unused.}
}
\value{
A list containing summary info from the input data.frame.
}
\description{
summary.aimsdf
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aimsdf-class.R
\name{is_aimsdf}
\alias{is_aimsdf}
\title{Checks if argument is a \code{aimsdf} object}
\usage{
is_aimsdf(x)
}
\arguments{
\item{x}{An \R object}
}
\description{
Checks if argument is a \code{aimsdf} object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{extract_map_coord}
\alias{extract_map_coord}
\title{extract_map_coord}
\usage{
extract_map_coord(x, ...)
}
\arguments{
\item{x}{An sfc_POINT}

\item{...}{Additional argument "pos" to internal function}
}
\description{
Internal
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{make_pretty_data_label}
\alias{make_pretty_data_label}
\title{make_pretty_data_label}
\usage{
make_pretty_data_label(x)
}
\arguments{
\item{x}{A character}
}
\description{
Internal
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{json_results}
\alias{json_results}
\title{\code{\link[jsonlite]{fromJSON}} data request}
\usage{
json_results(dt_req)
}
\arguments{
\item{dt_req}{An URL \code{\link[httr]{GET}} output}
}
\value{
\code{aims_data} returns a \code{\link[base]{data.frame}} of class
\code{\link{aimsdf}}.

If \code{summary \%in\% c("summary-by-series", "summary-by-deployment")},
the output shows the summary information for the target dataset (i.e.
weather or temperature loggers)
(NB: currently, \code{summary} only works for the temperature logger
database). If \code{summary} is \emph{not} passed as an additional argument, then
the output contains \strong{raw} monitoring data. If \code{summary = "daily"},
then the output contains \strong{mean daily aggregated} monitoring data.
The output also contains five attributes (empty strings if
\code{summary} is passed as an additional argument):
\itemize{
\item{\code{metadata}}{a \href{https://www.doi.org/}{DOI} link
containing the metadata record for the data series.}
\item{\code{citation}}{the citation information for the particular
dataset.}
\item{\code{parameters}}{The measured parameters comprised in the
output.}
\item{\code{type}}{The type of dataset. Either "monitoring" if
\code{summary} is not specified, "monitoring (daily aggregation)" if
\code{summary = "daily"}, or a "summary-by-" otherwise.}
\item{\code{target}}{The input target.}
}
}
\description{
Wrapper function
}
\details{
This function submits a \code{dt_req} data request via
\code{\link[jsonlite]{fromJSON}}.
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aimsdf-methods.R
\name{print.aimsdf}
\alias{print.aimsdf}
\title{print.aimsdf}
\usage{
\method{print}{aimsdf}(x, ...)
}
\arguments{
\item{x}{An object of class \code{\link{aimsdf}} as
returned by \code{\link{aims_data}}.}

\item{...}{Not used.}
}
\value{
A list containing a summary of the model fit as returned a
brmsfit for each model.
}
\description{
print.aimsdf
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aims_citation.R
\name{aims_citation}
\alias{aims_citation}
\title{Extracts citation attribute from object of class aimsdf}
\usage{
aims_citation(df_)
}
\arguments{
\item{df_}{A data.frame of class \code{\link{aimsdf}} created by
function \code{\link{aims_data}}}
}
\value{
A \code{\link[base]{character}} vector.
}
\description{
Extracts citation attribute from object of class aimsdf
}
\details{
This function retrieves the citation attribute from an
\code{\link{aimsdf}} object. If the input \code{\link{aimsdf}} object is
a summary data.frame (see ?\code{\link{aims_data}}), then output will be
an empty string.
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{update_format}
\alias{update_format}
\title{Format \code{\link[jsonlite]{fromJSON}} output list}
\usage{
update_format(results, doi)
}
\arguments{
\item{results}{A \code{\link[jsonlite]{fromJSON}} list
generated by \code{\link{json_results}}.}

\item{doi}{A \href{https://www.doi.org/}{Digital Object Identifier}
for a chosen
\href{https://open-aims.github.io/data-platform/}{AIMS data series}}
}
\value{
\code{aims_data} returns a \code{\link[base]{data.frame}} of class
\code{\link{aimsdf}}.

If \code{summary \%in\% c("summary-by-series", "summary-by-deployment")},
the output shows the summary information for the target dataset (i.e.
weather or temperature loggers)
(NB: currently, \code{summary} only works for the temperature logger
database). If \code{summary} is \emph{not} passed as an additional argument, then
the output contains \strong{raw} monitoring data. If \code{summary = "daily"},
then the output contains \strong{mean daily aggregated} monitoring data.
The output also contains five attributes (empty strings if
\code{summary} is passed as an additional argument):
\itemize{
\item{\code{metadata}}{a \href{https://www.doi.org/}{DOI} link
containing the metadata record for the data series.}
\item{\code{citation}}{the citation information for the particular
dataset.}
\item{\code{parameters}}{The measured parameters comprised in the
output.}
\item{\code{type}}{The type of dataset. Either "monitoring" if
\code{summary} is not specified, "monitoring (daily aggregation)" if
\code{summary = "daily"}, or a "summary-by-" otherwise.}
\item{\code{target}}{The input target.}
}
}
\description{
When \code{\link[jsonlite]{fromJSON}} returns a list, format list names
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aims_metadata.R
\name{aims_metadata}
\alias{aims_metadata}
\title{Extracts metadata attribute from object of class aimsdf}
\usage{
aims_metadata(df_)
}
\arguments{
\item{df_}{A data.frame of class \code{\link{aimsdf}} created by
function \code{\link{aims_data}}}
}
\value{
A \code{\link[base]{character}} vector.
}
\description{
Extracts metadata attribute from object of class aimsdf
}
\details{
This function retrieves the metadata attribute from an
\code{\link{aimsdf}} object. If the input \code{\link{aimsdf}} object is
a summary data.frame (see ?\code{\link{aims_data}}), then output will be
an empty string.
}
\author{
AIMS Datacentre \email{adc@aims.gov.au}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aimsdf-class.R
\docType{class}
\name{aimsdf-class}
\alias{aimsdf-class}
\alias{aimsdf}
\title{Class \code{aimsdf} of data.frame downloaded by the \pkg{dataaimsr} package}
\description{
Datasets downloaded by the
\code{\link[dataaimsr:dataaimsr-package]{dataaimsr}} package inherit
the \code{aimsdf} class, which is data.frame with three attributes.
}
\details{
See \code{methods(class = "aimsdf")} for an overview of available methods.
}
\seealso{
\code{\link{aims_data}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{make_pretty_colour}
\alias{make_pretty_colour}
\title{make_pretty_colour}
\usage{
make_pretty_colour(x, alpha_ = 0.55)
}
\arguments{
\item{x}{A character}

\item{alpha_}{A numeric}
}
\description{
Internal
}
\keyword{internal}
