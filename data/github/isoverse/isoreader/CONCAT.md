---
title: 'Isoreader: An R package to read stable isotope data files for reproducible research'
tags:
  - R
  - stable isotopes
  - earth science
  - ecology
  - isotope ratio mass spectrometry
  - reproducible research
authors:
  - name: Sebastian Kopf^[corresponding author]
    orcid: 0000-0002-2044-0201
    affiliation: 1
  - name: Brett Davidheiser-Kroll
    orcid: 0000-0002-6153-7851
    affiliation: 1
  - name: Ilja Kocken
    orcid: 0000-0003-2196-8718
    affiliation: 2
affiliations:
 - name: Department of Geological Sciences, University of Colorado Boulder, Colorado, USA
   index: 1
 - name: Department of Earth Sciences, Utrecht University, the Netherlands
   index: 2
date: 11 November 2020
bibliography: paper.bib
---

# Summary

The measurement and interpretation of the stable isotope composition of any material or molecule has widespread application in disciplines ranging from the earth sciences to ecology, anthropology, and forensics. The naturally occurring differences in the abundance of the stable isotopes of carbon, nitrogen, oxygen, and many other elements provide valuable insight into environmental conditions and sources, fluxes, and mechanisms of material transfer. Because isotopic variations in nature are very small, the measurement itself requires cutting edge analytical instrumentation using isotope ratio mass spectrometry (IRMS) as well as rigorous data reduction procedures for calibration and quality control. The `isoreader` package implements an easily extendable interface for IRMS data from common instrument vendor file formats and thus enables the reading and processing of stable isotope data directly from the source. This provides a foundational tool for platform-independent, efficient and reproducible data reduction.

# Statement of need

Reproducible data processing is a key prerequisite for efficient data exchange, methodological progress, and productive discourse in scientific research. However, generating a record of every step of a data processing pipeline in a format that is transparent and easy to understand is not an easy task. In the world of stable isotopes, many data processing steps require proprietary software for data access and depend on point-and-click interactions. This makes it challenging to share and discuss one’s approach, review others’ and compare calculations and datasets across laboratories. Moreover, it severely restricts opportunities for iteration, exchange of ideas, and data aggregation.

The `isoreader` package enables efficient and reproducible reading and processing of stable isotope data directly from the data files no matter which operating system (Windows, Mac, Linux). It is already being used for stable isotope data processing in several laboratories and recent publications including @Silverman2019, @Cheng2019, @Ingalls2020, and @Suarez2020. The `isoreader` package was designed to be easily extendable with readers for new file formats, and provides data export functionality to Python using the shared R/Python feather file format. This will enable the development, sharing and vetting of open-source data processing pipelines for stable isotope data across scientific disciplines.

# Acknowledgements

We thank Max Lloyd for his valuable contributions to parsing dual inlet file formats. We thank Seth Newsome and the executive committee of the IsoBank initiative for organizing a workshop that helped improve this software. We also thank all the laboratories that shared test files including the Bender, Bergmann, Bradley, Eiler, Kim, Ono, Pearson, Sessions, Sigman, and Snell Labs (all in the USA), and the Ziegler Lab (in the Netherlands), as well as the Stable Isotope Facilities at the University of California Davis (USA), University of Colorado (USA), University of New Brunswick (Canada), University of New Mexico (USA), University of Ottawa (Canada), University of Utah (USA), University of Washington (USA), University of Wyoming (USA), and the United States Geological Survey (USA). Lastly, we thank the members of the stable isotope community on the ISOGEOCHEM listserv for their contributions and bug reports that improved this software. This project was supported, in part, by grants to SHK from the U.S. National Science Foundation (EAR 1928303) and the University of Colorado Boulder.

# References

<!-- README.md is generated from README.Rmd. Please edit that file -->

# isoreader <a href='https://isoreader.isoverse.org'><img src='man/figures/isoreader_logo_thumb.png' align="right" height="138.5"/></a>

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/isoreader)](https://cran.r-project.org/package=isoreader)
[![Documentation](https://img.shields.io/badge/docs-online-green.svg)](https://isoreader.isoverse.org/)
[![R build
status](https://github.com/isoverse/isoreader/workflows/R-CMD-check/badge.svg)](https://github.com/isoverse/isoreader/actions?workflow=R-CMD-check)
[![Binder](https://img.shields.io/badge/explore%20online-in%20RStudio-blue.svg)](https://mybinder.org/v2/gh/isoverse/isoreader/binder?urlpath=rstudio)
[![Binder](https://img.shields.io/badge/explore%20online-in%20Jupyter-orange.svg)](https://mybinder.org/v2/gh/isoverse/isoreader/binder?urlpath=lab)

## About

This package is intended as a unified one-stop command line interface to
all common IRMS (isotope ratio mass spectrometry) file formats used in
stable isotope geochemistry. It is an extension and highly stream-lined
re-implementation of the proof-of-concept
[isoread](https://github.com/sebkopf/isoread) package and is designed to
fit into a larger framework of IRMS data tools that includes the
web-based graphical user interface package
[isoviewer](https://github.com/isoverse/isoviewer) and the data
processing and visualization pipeline
[isoprocessor](https://github.com/isoverse/isoprocessor).

[isoreader](https://isoreader.isoverse.org/) enables the reading and
processing of stable isotope data directly from the data files and thus
provides a tool for platform-independent (Windows, Mac, Linux),
efficient and reproducible data reduction. Although implemented in R, it
can be used in both RMarkdown as well as Jupyter data processing
notebooks and also provides functionality for easy export to Python
using the shared R/Python feather file format. At present, it can read
most Thermo dual inlet (.did, .caf) and continuous flow (.dxf, .cf) data
files as well as Elementar continuous flow data archives (.iarc) with
additional extensions for other file formats in the works. Due to the
dynamic implementation and design based on the popular
[tidyverse](https://www.tidyverse.org/) style of R programming,
isoreader is easily extendable, takes care of error catching to avoid
pipeline breaks due to problems encountered in source data files
(modeled after [readr](https://readr.tidyverse.org/)) and works great
with [tidyverse](https://www.tidyverse.org/) packages such as
[tidyr](https://tidyr.tidyverse.org/),
[dplyr](https://dplyr.tidyverse.org/) and
[ggplot](https://ggplot2.tidyverse.org/).

## Installation

You can install the latest release of isoreader from
[CRAN](https://cran.r-project.org/package=isoreader):

``` r
install.packages("isoreader")
```

Some isoreader features including Excel and feather export depend on
optional packages that are not required for the core functionality of
isoreader. To use this functionality, please install the following
packages manually if not already installed (isoreader will throw an
informative warning if they are needed but missing):

``` r
# optional extensions
install.packages(c("feather", "openxlsx", "xml2", "BiocManager"))
BiocManager::install("rhdf5")
```

To install the current development version of isoreader directly from
GitHub, please use the devtools package:

``` r
# installs the development tools package if not yet installed
if(!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("isoverse/isoreader")
```

Troubleshooting note: depending on your workspace and operating system,
you may have to re-start your R session or manually install some
dependencies. For example, the `digest` package sometimes causes trouble
- re-install with
`remove.packages("digest"); install.packages("digest")`.

## Show me some code

You can, for example, automatically read the data from **all** supported
scan files in a directory (and all its subdirectories) simply by
providing the path to the folder. The following code demonstrates this
with the example data files bundled with the `isoreader` package. For a
more detailed example including continuous flow and dual inlet file
reads, check out our [**Quick Start
Vignette**](https://isoreader.isoverse.org/articles/quick_start.html).

``` r
library(isoreader)
data_folder <- iso_get_reader_examples_folder()
iso_files <- iso_read_scan(data_folder)
#> Info: preparing to read 4 data files...
#> Info: reading file 'background_scan_example.scn' with '.scn' reader...
#> Info: reading file 'full_scan_example.scn' with '.scn' reader...
#> Info: reading file 'peak_shape_scan_example.scn' with '.scn' reader...
#> Info: reading file 'time_scan_example.scn' with '.scn' reader...
#> Info: finished reading 4 files in 0.82 secs

iso_files
#> Data from 4 scan iso files: 
#> # A tibble: 4 × 5
#>   file_id                     raw_data        file_info method_info file_path   
#>   <chr>                       <glue>          <chr>     <chr>       <chr>       
#> 1 background_scan_example.scn 525 measuremen… 8 entries resistors   background_…
#> 2 full_scan_example.scn       799 measuremen… 8 entries resistors   full_scan_e…
#> 3 peak_shape_scan_example.scn 220 measuremen… 8 entries resistors   peak_shape_…
#> 4 time_scan_example.scn       5532 measureme… 8 entries resistors   time_scan_e…
```

## Supported File Types

Currently supported file types:

| type            | extension | software  | description                         |
|:----------------|:----------|:----------|:------------------------------------|
| continuous flow | .cf       | Isodat    | Continuous Flow file format (older) |
| continuous flow | .cf.rds   | isoreader | R Data Storage                      |
| continuous flow | .dxf      | Isodat    | Continuous Flow file format (newer) |
| continuous flow | .iarc     | ionOS     | Continuous Flow data archive        |
| dual inlet      | .caf      | Isodat    | Dual Inlet file format (older)      |
| dual inlet      | .di.rds   | isoreader | R Data Storage                      |
| dual inlet      | .did      | Isodat    | Dual Inlet file format (newer)      |
| dual inlet      | .txt      | Nu        | Dual Inlet file format              |
| scan            | .scan.rds | isoreader | R Data Storage                      |
| scan            | .scn      | Isodat    | Scan file format                    |

## Documentation

-   for a quick introduction, check out the aforementioned [**Quick
    Start
    Vignette**](https://isoreader.isoverse.org/articles/quick_start.html)
-   for a full reference of all available functions, see the **[Function
    Reference](https://isoreader.isoverse.org/reference/)**
-   for function help within RStudio, simply start typing `?iso_` in the
    console and a list of available function will appear (all functions
    share the `iso_` prefix)
-   for a detailed example of how to work with continuous flow data
    files, see the vignette on **[Continuous
    Flow](https://isoreader.isoverse.org/articles/continuous_flow.html)**
-   for a detailed example of how to work with dual inlet data files,
    see the vignette on **[Dual
    Inlet](https://isoreader.isoverse.org/articles/dual_inlet.html)**
-   for a detailed example of how to work with scan data files, see the
    vignette on
    **[Scans](https://isoreader.isoverse.org/articles/scan.html)**

## Troubleshooting

If you run into a file format that is not currently supported or any
issues with supported formats, please file a request/bug report in the
[issue tracker](https://github.com/isoverse/isoreader/issues). Likewise
if you run into any unexpected behavior or uncaught errors. Most
isoreader functionality is continuously tested on Unix and Windows
systems using [GitHub
Actions](https://github.com/isoverse/isoreader/actions?workflow=R-CMD-check).
This makes it possible to ensure proper functionality and catch issues
quickly, however, sometimes something slips through or is not yet
automatically tested. We try to make sure to fix such errors as soon as
possible but ask for patience due to the small development team. If you
have the skills and are willing to fix problems yourself, that’s great,
please take a look at the development section below.

## Development

If you are interested in helping with development, that’s fantastic!
Please fork the repository and branch off from the [dev
branch](https://github.com/isoverse/isoreader/tree/dev) since it
contains the most up-to-date development version of
[isoreader](https://isoreader.isoverse.org/). Make sure to write
[`testthat` tests](https://r-pkgs.org/tests.html) for your work (stored
in the tests/testthat directory). All tests can be run automatically and
continuously during development to make it easier to spot any code
problems on the go. The easiest way to run them is by running
`make auto_test` in the [isoreader](https://isoreader.isoverse.org/)
directory from command line (it will test everything automatically in a
completely separate R session).

## Open Source

[isoreader](https://isoreader.isoverse.org/) is and will always be fully
open-source (i.e. free as in **freedom** and free as in **free beer**)
and is provided as is. The source code is released under GPL-2.

## isoverse <a href='https://www.isoverse.org'><img src='man/figures/isoverse_logo_thumb.png' align="right" height="138.5"/></a>

This package is part of the isoverse suite of data tools for stable
isotopes. If you like the functionality that isoverse packages provide
to the geochemical community, please help us spread the word and include
an isoverse or individual package logo on one of your posters or slides.
All logos are posted in high resolution in [this
repository](https://github.com/isoverse/logos).
# isoreader 1.3.0

## Major changes

 - functions for previously deprecated `.rda` file format removed
 - R version (>= 4.0.0) and dependency version requirements clarified
 - dependencies simplified (`openxlsx`, `feather`, `xml2` and `rhdf5` are now optional extensions)

## Bug fixes

 - dependency on pandoc removed

# isoreader 1.2.3

## Major Features

#### Reading isotope ratio mass spectrometry (IRMS) data files

This package currently supports reading the following raw data files:

 - continuous flow data files from mass spectrometers run by Isodat and IonOS software (`?iso_read_continuous_flow`)
 - dual inlet data files from mass spectrometers run by Isodat and Nu software (`?iso_read_dual_inlet`)
 - scan data files from mass spectrometers run by Isodat software (`?iso_read_scan`)

#### Aggregating data from files

This package provides the following data aggregation and data processing functionality for all supported data files:

 - aggregating file information including sequence line data (`?iso_get_file_info`)
 - aggregating raw mass spectrometric data (`?iso_get_raw_data`)
 - aggregating mass spectrometric background data (`?iso_get_bgrd_data`)
 - aggregating vendor-specific data tables if provided in raw file formats (`?iso_get_vendor_data_table`)
 - aggregating information on detector resistors if included in raw file formats (`?iso_get_resistors`)
 - aggregating information on internal isotope standards if included in raw file formats (`?iso_get_standards`)
 - processing file information with tidyverse-equivalent `select`, `rename`, `mutate` and `filter` functions (`?iso_filter_files`, `?iso_select_file_info`, `?iso_rename_file_info`, `?iso_mutate_file_info`)
 - adding file information with `?iso_add_file_info`to simplify sequential join operations


#### Exporting information

This package provides the following data export functionality for all supported data files:

 - export to open Excel (.xslx) with `?iso_export_to_excel`
 - export to the Python/R cross-over feather file format with `?iso_export_to_feather`
 Update of an existing CRAN package to address problems during vignette re-building that have arisen on systems without pandoc (e.g Solaris) as pointed out by Prof. Brian Ripley in a message on Feb. 12 2021 to the maintainers of all affected packages. The package has been tested on GitHub without pandoc and all vignette errors have been addressed (mostly stemming from the use of `rmarkdown::paged_table()`).

## Test environments

* Local OS X install, R 4.0.2 (release)
* Mac OS X 10.15.6 (on GitHub), R 4.0.2 (release)
* Ubuntu 16.04 (on GitHub), R 4.0.2 (release)
* Windows Server 2019 (on GitHub), R 4.0.2 (release)
* Win-builder (release and devel)

## R CMD check results

There were no ERRORs, no WARNINGs, and no NOTEs. 

## Downstream dependencies

There are currently no downstream dependencies for this package.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)

version <- as.character(packageVersion("isoreader"))
```

# isoreader <a href='https://isoreader.isoverse.org'><img src='man/figures/isoreader_logo_thumb.png' align="right" height="138.5"/></a>

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/isoreader)](https://cran.r-project.org/package=isoreader)
[![Documentation](https://img.shields.io/badge/docs-online-green.svg)](https://isoreader.isoverse.org/)
[![R build status](https://github.com/isoverse/isoreader/workflows/R-CMD-check/badge.svg)](https://github.com/isoverse/isoreader/actions?workflow=R-CMD-check)
[![Binder](https://img.shields.io/badge/explore%20online-in%20RStudio-blue.svg)](https://mybinder.org/v2/gh/isoverse/isoreader/binder?urlpath=rstudio)
[![Binder](https://img.shields.io/badge/explore%20online-in%20Jupyter-orange.svg)](https://mybinder.org/v2/gh/isoverse/isoreader/binder?urlpath=lab)

## About

This package is intended as a unified one-stop command line interface to all common IRMS (isotope ratio mass spectrometry) file formats used in stable isotope geochemistry. It is an extension and highly stream-lined re-implementation of the proof-of-concept [isoread](https://github.com/sebkopf/isoread) package and is designed to fit into a larger framework of IRMS data tools that includes the web-based graphical user interface package [isoviewer](https://github.com/isoverse/isoviewer) and the data processing and visualization pipeline [isoprocessor](https://github.com/isoverse/isoprocessor).

[isoreader](https://isoreader.isoverse.org/) enables the reading and processing of stable isotope data directly from the data files and thus provides a tool for platform-independent (Windows, Mac, Linux), efficient and reproducible data reduction. Although implemented in R, it can be used in both RMarkdown as well as Jupyter data processing notebooks and also provides functionality for easy export to Python using the shared R/Python feather file format. At present, it can read most Thermo dual inlet (.did, .caf) and continuous flow (.dxf, .cf) data files as well as Elementar continuous flow data archives (.iarc) with additional extensions for other file formats in the works. Due to the dynamic implementation and design based on the popular [tidyverse](https://www.tidyverse.org/) style of R programming, isoreader is easily extendable, takes care of error catching to avoid pipeline breaks due to problems encountered in source data files (modeled after  [readr](https://readr.tidyverse.org/)) and works great with [tidyverse](https://www.tidyverse.org/) packages such as [tidyr](https://tidyr.tidyverse.org/), [dplyr](https://dplyr.tidyverse.org/) and [ggplot](https://ggplot2.tidyverse.org/).

## Installation

You can install the latest release of isoreader from [CRAN](https://cran.r-project.org/package=isoreader):

```{r cran-installation, eval = FALSE}
install.packages("isoreader")
```

Some isoreader features including Excel and feather export depend on optional packages that are not required for the core functionality of isoreader. To use this functionality, please install the following packages manually if not already installed (isoreader will throw an informative warning if they are needed but missing):

```{r, optional-installation, eval = FALSE}
# optional extensions
install.packages(c("feather", "openxlsx", "xml2", "BiocManager"))
BiocManager::install("rhdf5")
```

To install the current development version of isoreader directly from GitHub, please use the devtools package:

```{r gh-installation, eval = FALSE}
# installs the development tools package if not yet installed
if(!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("isoverse/isoreader")
```

Troubleshooting note: depending on your workspace and operating system, you may have to re-start your R session or manually install some dependencies. For example, the `digest` package sometimes causes trouble - re-install with `remove.packages("digest"); install.packages("digest")`.

## Show me some code

You can, for example, automatically read the data from **all** supported scan files in a directory (and all its subdirectories) simply by providing the path to the folder. The following code demonstrates this with the example data files bundled with the `isoreader` package. For a more detailed example including continuous flow and dual inlet file reads, check out our [**Quick Start Vignette**](https://isoreader.isoverse.org/articles/quick_start.html).

```{r, message = -c(1:3), echo = -c(2:3)}
library(isoreader)
iso_turn_reader_caching_off() # make sure reading fresh
setwd(tempdir()) # make sure no wd artifacts
data_folder <- iso_get_reader_examples_folder()
iso_files <- iso_read_scan(data_folder)

iso_files
```

## Supported File Types

Currently supported file types:

```{r, echo=FALSE, warning=FALSE, message=FALSE}
library(isoreader)
iso_get_supported_file_types() %>% 
  dplyr::select(-call) %>%
  knitr::kable()
```

## Documentation

 - for a quick introduction, check out the aforementioned [**Quick Start Vignette**](https://isoreader.isoverse.org/articles/quick_start.html)
 - for a full reference of all available functions, see the **[Function Reference](https://isoreader.isoverse.org/reference/)**
 - for function help within RStudio, simply start typing `?iso_` in the console and a list of available function will appear (all functions share the `iso_` prefix)
 - for a detailed example of how to work with continuous flow data files, see the vignette on **[Continuous Flow](https://isoreader.isoverse.org/articles/continuous_flow.html)**
 - for a detailed example of how to work with dual inlet data files, see the vignette on **[Dual Inlet](https://isoreader.isoverse.org/articles/dual_inlet.html)**
 - for a detailed example of how to work with scan data files, see the vignette on **[Scans](https://isoreader.isoverse.org/articles/scan.html)**

## Troubleshooting

If you run into a file format that is not currently supported or any issues with supported formats, please file a request/bug report in the [issue tracker](https://github.com/isoverse/isoreader/issues). Likewise if you run into any unexpected behavior or uncaught errors. Most isoreader functionality is continuously tested on Unix and Windows systems using [GitHub Actions](https://github.com/isoverse/isoreader/actions?workflow=R-CMD-check). This makes it possible to ensure proper functionality and catch issues quickly, however, sometimes something slips through or is not yet automatically tested. We try to make sure to fix such errors as soon as possible but ask for patience due to the small development team. If you have the skills and are willing to fix problems yourself, that's great, please take a look at the development section below.

## Development

If you are interested in helping with development, that's fantastic! Please fork the repository and branch off from the [dev branch](https://github.com/isoverse/isoreader/tree/dev) since it contains the most up-to-date development version of [isoreader](https://isoreader.isoverse.org/). Make sure to write [```testthat``` tests](https://r-pkgs.org/tests.html) for your work (stored in the tests/testthat directory). All tests can be run automatically and continuously during development to make it easier to spot any code problems on the go. The easiest way to run them is by running ```make auto_test``` in the [isoreader](https://isoreader.isoverse.org/) directory from command line (it will test everything automatically in a completely separate R session).

## Open Source

[isoreader](https://isoreader.isoverse.org/) is and will always be fully open-source (i.e. free as in **freedom** and free as in **free beer**) and is provided as is. The source code is released under GPL-2.

## isoverse <a href='https://www.isoverse.org'><img src='man/figures/isoverse_logo_thumb.png' align="right" height="138.5"/></a>

This package is part of the isoverse suite of data tools for stable isotopes. If you like the functionality that isoverse packages provide to the geochemical community, please help us spread the word and include an isoverse or individual package logo on one of your posters or slides. All logos are posted in high resolution in [this repository](https://github.com/isoverse/logos).
---
title: "Continuous Flow Examples"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
  html_document:
    code_folding: show
    df_print: paged
    number_sections: yes
    toc: yes
    toc_depth: 3
    toc_float: yes
editor_options:
  chunk_output_type: console
vignette: >
  %\VignetteIndexEntry{Continuous Flow Examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction

Isoreader supports several continuous flow IRMS data formats. This vignette shows some of the functionality for continuous flow files. For additional information on operations more generally (caching, combining read files, data export, etc.), please consult the [operations vignette](https://isoreader.isoverse.org/articles/operations.html). For details on downstream data processing and visualization, see the [isoprocessor package](https://isoprocessor.isoverse.org).

```{r, message=FALSE}
# load isoreader package
library(isoreader)
```

# Reading files

Reading continuous flow files is as simple as passing one or multiple file or folder paths to the `iso_read_continuous_flow()` function. If folders are provided, any files that have a recognized continuous flow file extensions within those folders will be processed (e.g. all `.dxf`, `.cf` and `.iarc`). Here we read several files that are bundled with the package as examples (and whose paths can be retrieved using the `iso_get_reader_example()` function). Note that some of the files (.cf, .dxf) are individual analysis files whereas others (.iarc) are collections of several files.

```{r, message = FALSE}
# all available examples
iso_get_reader_examples() %>% knitr::kable()
```

```{r}
# read a few of the continuous flow examples
cf_files <-
  iso_read_continuous_flow(
    iso_get_reader_example("continuous_flow_example.cf"),
    iso_get_reader_example("continuous_flow_example.iarc"),
    iso_get_reader_example("continuous_flow_example.dxf")
  )
```

# File summary

The `cf_files` variable now contains a set of isoreader objects, one for each file. Take a look at what information was retrieved from the files using the `iso_get_data_summary()` function.

```{r}
cf_files %>% iso_get_data_summary() %>% knitr::kable()
```

## Problems

In case there was any trouble with reading any of the files, the following functions provide an overview summary as well as details of all errors and warnings, respectively. The examples here contain no errors but if you run into any unexpected file read problems, please file a bug report in the [isoreader issue tracker](https://github.com/isoverse/isoreader/issues).

```{r}
cf_files %>% iso_get_problems_summary() %>% knitr::kable()
cf_files %>% iso_get_problems() %>% knitr::kable()
```

# File Information

Detailed file information can be aggregated for all isofiles using the `iso_get_file_info()` function which supports the full [select syntax](https://dplyr.tidyverse.org/reference/select.html) of the [dplyr](https://dplyr.tidyverse.org/) package to specify which columns are of interest (by default, all file information is retrieved). Additionally, file information from different file formats can be renamed to the same column name for easy of downstream processing. The following provides a few examples for how this can be used (the names of the interesting info columns may vary between different file formats):

```{r}
# all file information
cf_files %>% iso_get_file_info(select = c(-file_root)) %>% knitr::kable()
# select file information
cf_files %>%
  iso_get_file_info(
    select = c(
       # rename sample id columns from the different file types to a new ID column
      ID = `Identifier 1`, ID = `Name`,
      # select columns without renaming
      Analysis, `Peak Center`, `H3 Factor`,
      # select the time stamp and rename it to `Date & Time`
      `Date & Time` = file_datetime
    ),
    # explicitly allow for file specific rename (for the new ID column)
    file_specific = TRUE
  ) %>% knitr::kable()
```

## Select/Rename

Rather than retrieving specific file info columns using the above example of `iso_get_file_info(select = ...)`, these information can also be modified across an entire collection of isofiles using the `iso_select_file_info()` and `iso_rename_file_info()` functions. For example, the above example could be similarly achieved with the following use of `iso_select_file_info()`:

```{r}
# select + rename specific file info columns
cf_files2 <- cf_files %>%
  iso_select_file_info(
    ID = `Identifier 1`, ID = `Name`, Analysis, `Peak Center`, `H3 Factor`,
    `Date & Time` = file_datetime,
    # recode to the same name in different files
    `Sample Weight` = `Identifier 2`, `Sample Weight` = `EA Sample Weight`,
    file_specific = TRUE
  )

# fetch all file info
cf_files2 %>% iso_get_file_info() %>% knitr::kable()
```

## Filter

Any collection of isofiles can also be filtered based on the available file information using the function `iso_filter_files`. This function can operate on any column available in the file information and supports full [dplyr](https://dplyr.tidyverse.org/reference/filter.html) syntax.

```{r}
# find files that have 'acetanilide' in the new ID field
cf_files2 %>% iso_filter_files(grepl("acetanilide", ID)) %>%
  iso_get_file_info() %>%
  knitr::kable()

# find files that were run since 2015
cf_files2 %>%
  iso_filter_files(`Date & Time` > "2015-01-01") %>%
  iso_get_file_info() %>%
  knitr::kable()
```

## Mutate

The file information in any collection of isofiles can also be mutated using the function `iso_mutate_file_info`. This function can introduce new columns and operate on any existing columns available in the file information (even if it does not exist in all files) and supports full [dplyr](https://dplyr.tidyverse.org/reference/mutate.html) syntax. It can also be used in conjunction with `iso_with_unit` to generate values with implicit units.

```{r}
cf_files3 <-
  cf_files2 %>%
  iso_mutate_file_info(
    # update existing column
    ID = paste("ID:", ID),
    # introduce new column
    `Run since 2015?` = `Date & Time` > "2015-01-01",
    # parse weight as a number and turn into a column with units
    `Sample Weight` = `Sample Weight` %>% parse_number() %>% iso_with_units("mg")
  )

cf_files3 %>%
  iso_get_file_info() %>%
  iso_make_units_explicit() %>%
  knitr::kable()
```

## Add

Additionally, a wide range of new file information can be added in the form of a data frame with any number of columns (usually read from a comma-separated-value/csv file or an Excel/xlsx file) using the function `iso_add_file_info` and specifying which existing file information should be used to merge in the new information. It is similar to [dplyr's left_join](https://dplyr.tidyverse.org/reference/mutate-joins.html) but with additional safety checks and the possibility to join the new information sequentially as illustrated below.

```{r}
# this kind of information data frame is frequently read in from a csv or xlsx file
new_info <-
  dplyr::bind_rows(
    # new information based on new vs. old samples
    dplyr::tribble(
      ~file_id, ~`Run since 2015?`,  ~process,  ~info,
       NA,       TRUE,                "yes",     "new runs",
       NA,       FALSE,               "yes",     "old runs"
    ),
    # new information for a single specific file
    dplyr::tribble(
      ~file_id,        ~process,  ~note,
       "6617_IAEA600",  "no",      "did not inject properly"
    )
  )
new_info %>% knitr::kable()

# adding it to the isofiles
cf_files3 %>%
  iso_add_file_info(new_info, by1 = "Run since 2015?", by2 = "file_id") %>%
  iso_get_file_info(select = !!names(new_info)) %>%
  knitr::kable()
```


## Parse

Most file information is initially read as text to avoid cumbersome specifications during the read process and compatibility issues between different IRMS file formats. However, many file info columns are not easily processed as text. The isoreader package therefore provides several parsing and data extraction functions to facilitate processing the text-based data (some via functionality implemented by the [readr](https://readr.tidyverse.org) package). See code block below for examples. For a complete overview, see the `?extract_data` and `?iso_parse_file_info` documentation.

```{r}
# use parsing and extraction in iso_mutate_file_info
cf_files2 %>%
  iso_mutate_file_info(
    # change type of Peak Center to logical
    `Peak Center` = parse_logical(`Peak Center`),
    # retrieve first word of file_id
    file_id_1st = extract_word(file_id),
    # retrieve second word of ID column
    file_id_2nd = extract_word(file_id, 2),
    # retrieve file extension from the file_id using regular expression
    name = extract_substring(ID, "(\\w+)-?(.*)?", capture_bracket = 1)
  ) %>%
  iso_get_file_info(select = c(matches("file_id"), ID, name, `Peak Center`)) %>%
  knitr::kable()

# use parsing in iso_filter_file_info
cf_files2 %>%
  iso_filter_files(parse_number(`H3 Factor`) > 2) %>%
  iso_get_file_info() %>%
  knitr::kable()

# use iso_parse_file_info for simplified parsing of column data types
cf_files2 %>%
  iso_parse_file_info(
    integer = Analysis,
    number = `H3 Factor`,
    logical = `Peak Center`
  ) %>%
  iso_get_file_info() %>%
  knitr::kable()
```

# Resistors

Additionally, some IRMS data files contain resistor information that are useful for downstream calculations (see e.g. section on signal conversion later in this vignette):

```{r}
cf_files %>% iso_get_resistors() %>% knitr::kable()
```

# Reference values

As well as isotopic reference values for the different gases:

```{r}
# reference delta values without ratio values
cf_files %>% iso_get_standards(file_id:reference) %>% knitr::kable()
# reference values with ratios
cf_files %>% iso_get_standards() %>% knitr::kable()
```


# Raw Data

The raw data read from the IRMS files can be retrieved similarly using the `iso_get_raw_data()` function. Most data aggregation functions also allow for inclusion of file information using the `include_file_info` parameter, which functions identically to the `select` parameter of the `iso_get_file_info` function discussed earlier.

```{r}
# get raw data with default selections (all raw data, no additional file info)
cf_files %>% iso_get_raw_data() %>% head(n=10) %>% knitr::kable()
# get specific raw data and add some file information
cf_files %>%
  iso_get_raw_data(
    # select just time and the m/z 2 and 3 ions
    select = c(time.s, v2.mV, v3.mV),
    # include the Analysis number fron the file info and rename it to 'run'
    include_file_info = c(run = Analysis)
  ) %>%
  # look at first few records only
  head(n=10) %>% knitr::kable()
```

# Data Processing

The isoreader package is intended to make raw stable isotope data easily accessible. However, as with most analytical data, there is significant downstream processing required to turn these raw intensity chromatograms into peak-specific, properly referenced isotopic measurements. This and similar functionality as well as data visualization is part of the [isoprocessor package](https://isoprocessor.isoverse.org) which takes isotopic data through the various corrections in a transparent, efficient and reproducible manner.

That said, most vendor software also performs some of these calculations and it can be useful to be able to compare new data reduction procedures against those implemented in the vendor software. For this purpose, isoreader retrieves vendor computed data tables whenever possible, as illustrated below.

## Vendor Data Table

As with most data retrieval functions, the `iso_get_vendor_data_table()` function also allows specific column selection (by default, all columns are selected) and easy addition of file information via the `include_file_info` parameter (by default, none is included).

```{r}
# entire vendor data table
cf_files %>% iso_get_vendor_data_table() %>% knitr::kable()
# get specific parts and add some file information
cf_files %>%
  iso_get_vendor_data_table(
    # select peak number, ret. time, overall intensity and all H delta columns
    select = c(Nr., Rt, area = `rIntensity All`, matches("^d \\d+H")),
    # include the Analysis number fron the file info and rename it to 'run'
    include_file_info = c(run = Analysis)
  ) %>%
  knitr::kable()

# the data table also provides units if included in the original data file
# which can be made explicit using the function iso_make_units_explicit()
cf_files %>%
  iso_get_vendor_data_table(
    # select peak number, ret. time, overall intensity and all H delta columns
    select = c(Nr., Rt, area = `rIntensity All`, matches("^d \\d+H")),
    # include the Analysis number fron the file info and rename it to 'run'
    include_file_info = c(run = Analysis)
  ) %>%
  # make column units explicit
  iso_make_units_explicit() %>%
  knitr::kable()
```

# For expert users: retrieving all data

For users familiar with the nested data frames from the [tidyverse](https://www.tidyverse.org/) (particularly [tidyr](https://tidyr.tidyverse.org/)'s `nest` and `unnest`), there is an easy way to retrieve all data from the iso file objects in a single nested data frame:

```{r}
all_data <- cf_files %>% iso_get_all_data()
# not printed out because this data frame is very big
```

# Saving collections

Saving entire collections of isofiles for retrieval at a later point is easily done using the `iso_save` function which stores collections or individual isoreader file objects in the efficient R data storage format `.rds` (if not specified, the extension `.cf.rds` will be automatically appended). These saved collections can be conveniently read back using the same `iso_read_continuous_flow` command used for raw data files.

```{r}
# export to R data archive
cf_files %>% iso_save("cf_files_export.cf.rds")

# read back the exported R data archive
cf_files <- iso_read_continuous_flow("cf_files_export.cf.rds")
cf_files %>% iso_get_data_summary() %>% knitr::kable()
```


# Data Export

At the moment, isoreader supports export of all data to Excel and the [Feather file format](https://blog.rstudio.com/2016/03/29/feather/) (a Python/R cross-over format). Note that both export methods have similar syntax and append the appropriate file extension for each type of export file (`.cf.xlsx` and `.cf.feather`, respectively).

```{r, eval = FALSE}
# export to excel
cf_files %>% iso_export_to_excel("cf_files_export")

# data sheets available in the exported data file:
readxl::excel_sheets("cf_files_export.cf.xlsx")
```

```{r, eval=FALSE}
# export to feather
cf_files %>% iso_export_to_feather("cf_files_export")

# exported feather files
list.files(pattern = ".cf.feather")
```
---
title: "Quick Start Guide"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
  html_document:
    code_folding: show
    df_print: paged
    number_sections: yes
    toc: yes
    toc_depth: 3
    toc_float: yes
editor_options:
  chunk_output_type: console
vignette: >
  %\VignetteIndexEntry{Quick Start Guide}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction

Isoreader supports various dual inlet, continuous flow, and scan file formats. This vignette shows how to get started reading these raw IRMS data files and exporting the information to Excel. For more details on isoreader functionality and each file type, please read the **Full Examples** vignettes. For more information on downstream processing with isoverse, check out the [isoprocessor](https://isoprocessor.isoverse.org) package.

```{r, message=FALSE}
# load isoreader package
library(isoreader)
```

```{r, include=FALSE, message=FALSE}
iso_turn_reader_caching_off()
```


# Data files

For demonstration purposes, this vignette simply reads all supported dual inlet, continuous flow, and scan files that are bundled with the isoreader package. For a more detailed introduction to each file type, check out the specific vignettes for each:

 - [Continuous Flow](https://isoreader.isoverse.org/articles/continuous_flow.html)
 - [Dual Inlet](https://isoreader.isoverse.org/articles/dual_inlet.html)
 - [Scans](https://isoreader.isoverse.org/articles/scan.html)


```{r}
# all available examples
iso_get_reader_examples() %>% knitr::kable()
```

# Dual Inlet Files

```{r}
# read all available examples
di_files <- iso_read_dual_inlet(iso_get_reader_examples_folder())
```

```{r}
# save as r data storage (read back in with iso_read_dual_inlet)
iso_save(di_files, filepath = "di_save")
```

```{r, eval = FALSE}
# export to excel
iso_export_to_excel(di_files, filepath = "di_export")
```

# Continuous Flow Files

```{r}
# read all available examples
cf_files <- iso_read_continuous_flow(iso_get_reader_examples_folder())
```

```{r}
# save as r data storage (read back in with iso_read_continuous_flow)
iso_save(cf_files, filepath = "cf_save")
```

```{r, eval = FALSE}
# export to excel
iso_export_to_excel(cf_files, filepath = "cf_export")
```

# Scan Files

```{r}
# read all available examples
scan_files <- iso_read_scan(iso_get_reader_examples_folder())
```

```{r}
# save as r data storage (read back in with iso_read_scan)
iso_save(scan_files, filepath = "scan_save")
```

```{r, eval = FALSE}
# export to excel
iso_export_to_excel(scan_files, filepath = "scan_export")
```

---
title: "Dual Inlet Examples"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
  html_document:
    code_folding: show
    df_print: paged
    number_sections: yes
    toc: yes
    toc_depth: 3
    toc_float: yes
editor_options:
  chunk_output_type: console
vignette: >
  %\VignetteIndexEntry{Dual Inlet Examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction

Isoreader supports several dual inlet IRMS data formats. This vignette shows some of the functionality for dual inlet data files. For additional information on operations more generally (caching, combining read files, data export, etc.), please consult the [operations vignette](https://isoreader.isoverse.org/articles/operations.html). For details on downstream data processing and visualization, see the [isoprocessor package](https://isoprocessor.isoverse.org).


```{r, message=FALSE}
# load isoreader package
library(isoreader)
```


# Reading files

Reading dual inlet files is as simple as passing one or multiple file or folder paths to the `iso_read_dual_inlet()` function. If folders are provided, any files that have a recognized continuous flow file extensions within those folders will be processed (e.g. all `.did` and `.caf`). Here we read several files that are bundled with the package as examples (and whose paths can be retrieved using the `iso_get_reader_example()` function).

```{r, message=FALSE}
# all available examples
iso_get_reader_examples() %>% knitr::kable()
```

```{r}
# read dual inlet examples
di_files <-
  iso_read_dual_inlet(
    iso_get_reader_example("dual_inlet_example.did"),
    iso_get_reader_example("dual_inlet_example.caf"),
    iso_get_reader_example("dual_inlet_nu_example.txt"),
    nu_masses = 49:44
  )
```

# File summary

The `di_files` variable now contains a set of isoreader objects, one for each file. Take a look at what information was retrieved from the files using the `iso_get_data_summary()` function.

```{r}
di_files %>% iso_get_data_summary() %>% knitr::kable()
```

## Problems

In case there was any trouble with reading any of the files, the following functions provide an overview summary as well as details of all errors and warnings, respectively. The examples here contain no errors but if you run into any unexpected file read problems, please file a bug report in the [isoreader issue tracker](https://github.com/isoverse/isoreader/issues).

```{r}
di_files %>% iso_get_problems_summary() %>% knitr::kable()
di_files %>% iso_get_problems() %>% knitr::kable()
```

# File Information

Detailed file information can be aggregated for all isofiles using the `iso_get_file_info()` function which supports the full [select syntax](https://dplyr.tidyverse.org/reference/select.html) of the [dplyr](https://dplyr.tidyverse.org/) package to specify which columns are of interest (by default, all file information is retrieved). Additionally, file information from different file formats can be renamed to the same column name for easy of downstream processing. The following provides a few examples for how this can be used (the names of the interesting info columns may vary between different file formats):

```{r}
# all file information
di_files %>% iso_get_file_info(select = c(-file_root)) %>% knitr::kable()
# select file information
di_files %>%
  iso_get_file_info(
    select = c(
      # rename sample id columns from the different file types to a new ID column
      ID = `Identifier 1`, ID = `Sample Name`,
      # select columns without renaming
      Analysis, Method, `Peak Center`,
      # select the time stamp and rename it to `Date & Time`
      `Date & Time` = file_datetime,
      # rename weight columns from the different file types
      `Sample Weight`, `Sample Weight` = `Weight [mg]`
    ),
    # explicitly allow for file specific rename (for the new ID column)
    file_specific = TRUE
  ) %>% knitr::kable()
```

## Select/Rename

Rather than retrieving specific file info columns using the above example of `iso_get_file_info(select = ...)`, these information can also be modified across an entire collection of isofiles using the `iso_select_file_info()` and `iso_rename_file_info()` functions. For example, the above example could be similarly achieved with the following use of `iso_select_file_info()`:

```{r}
# select + rename specific file info columns
di_files2 <- di_files %>%
  iso_select_file_info(
    ID = `Identifier 1`, ID = `Sample Name`, Analysis, Method,
    `Peak Center`, `Date & Time` = file_datetime,
    `Sample Weight`, `Sample Weight` = `Weight [mg]`,
    file_specific = TRUE
  )

# fetch all file info
di_files2 %>% iso_get_file_info() %>% knitr::kable()
```

## Filter

Any collection of isofiles can also be filtered based on the available file information using the function `iso_filter_files`. This function can operate on any column available in the file information and supports full [dplyr](https://dplyr.tidyverse.org/reference/filter.html) syntax.

```{r}
# find files that have 'CIT' in the new ID field
di_files2 %>% iso_filter_files(grepl("CIT", ID)) %>%
  iso_get_file_info() %>%
  knitr::kable()

# find files that were run in 2017
di_files2 %>%
  iso_filter_files(`Date & Time` > "2017-01-01" & `Date & Time` < "2018-01-01") %>%
  iso_get_file_info() %>%
  knitr::kable()
```

## Mutate

The file information in any collection of isofiles can also be mutated using the function `iso_mutate_file_info`. This function can introduce new columns and operate on/overwrite any existing columns available in the file information (even if it does not exist in all files) and supports full [dplyr](https://dplyr.tidyverse.org/reference/mutate.html) syntax. It can also be used in conjunction with `iso_with_unit` to generate values with implicit units.

```{r}
di_files3 <- di_files2 %>%
  iso_mutate_file_info(
    # update existing column
    ID = paste("ID:", ID),
    # introduce new column
    `Run in 2017?` = `Date & Time` > "2017-01-01" & `Date & Time` < "2018-01-01",
    # parse weight as a number and turn into a column with units
    `Sample Weight` = `Sample Weight` %>% parse_number() %>% iso_with_units("mg")
  )

di_files3 %>%
  iso_get_file_info() %>%
  iso_make_units_explicit() %>%
  knitr::kable()
```

## Add

Additionally, a wide range of new file information can be added in the form of a data frame with any number of columns (usually read from a comma-separated-value/csv file or an Excel/xlsx file) using the function `iso_add_file_info` and specifying which existing file information should be used to merge in the new information. It is similar to [dplyr's left_join](https://dplyr.tidyverse.org/reference/mutate-joins.html) but with additional safety checks and the possibility to join the new information sequentially as illustrated below.

```{r}
# this kind of information data frame is frequently read in from a csv or xlsx file
new_info <-
  dplyr::bind_rows(
    # new information based on new vs. old samples
    dplyr::tribble(
      ~Analysis, ~`Run in 2017?`,  ~process,  ~info,
       NA,       TRUE,              "yes",     "2017 runs",
       NA,       FALSE,             "yes",     "other runs"
    ),
    # new information for a single specific file
    dplyr::tribble(
      ~Analysis, ~process,  ~note,
       "16068",   "no",      "did not inject properly"
    )
  )
new_info %>% knitr::kable()

# adding it to the isofiles
di_files3 %>%
  iso_add_file_info(new_info, by1 = "Run in 2017?", by2 = "Analysis") %>%
  iso_get_file_info(select = !!names(new_info)) %>%
  knitr::kable()
```


## Parse

Most file information is initially read as text to avoid cumbersome specifications during the read process and compatibility issues between different IRMS file formats. However, many file info columns are not easily processed as text. The isoreader package therefore provides several parsing and data extraction functions to facilitate processing the text-based data (some via functionality implemented by the [readr](https://readr.tidyverse.org) package). See code block below for examples. For a complete overview, see the `?extract_data` and `?iso_parse_file_info` documentation.

```{r}
# use parsing and extraction in iso_mutate_file_info
di_files2 %>%
  iso_mutate_file_info(
    # change type of Peak Center to logical
    `Peak Center` = parse_logical(`Peak Center`),
    # retrieve first word of Method column
    Method_1st = extract_word(Method),
    # retrieve second word of Method column
    Method_2nd = extract_word(Method, 2),
    # retrieve file extension from the file_id using regular expression
    extension = extract_substring(file_id, "\\.(\\w+)$", capture_bracket = 1)
  ) %>%
  iso_get_file_info(select = c(extension, `Peak Center`, matches("Method"))) %>%
  knitr::kable()

# use parsing in iso_filter_file_info
di_files2 %>%
  iso_filter_files(parse_integer(Analysis) > 1500) %>%
  iso_get_file_info() %>%
  knitr::kable()

# use iso_parse_file_info for simplified parsing of column data types
di_files2 %>%
  iso_parse_file_info(
    integer = Analysis,
    number = `Sample Weight`,
    logical = `Peak Center`
  ) %>%
  iso_get_file_info() %>%
  knitr::kable()
```

# Resistors

Additionally, some IRMS data files contain resistor information that are useful for downstream calculations (see e.g. section on signal conversion later in this vignette):

```{r}
di_files %>% iso_get_resistors() %>% knitr::kable()
```

# Reference values

As well as isotopic reference values for the different gases:

```{r}
# reference delta values without ratio values
di_files %>% iso_get_standards(file_id:reference) %>% knitr::kable()
# reference values with ratios
di_files %>% iso_get_standards() %>% knitr::kable()
```

# Raw Data

The raw data read from the IRMS files can be retrieved similarly using the `iso_get_raw_data()` function. Most data aggregation functions also allow for inclusion of file information using the `include_file_info` parameter, which functions identically to the `select` parameter of the `iso_get_file_info` function discussed earlier.

```{r}
# get raw data with default selections (all raw data, no additional file info)
di_files %>% iso_get_raw_data() %>% head(n=10) %>% knitr::kable()
# get specific raw data and add some file information
di_files %>%
  iso_get_raw_data(
    # select just time and the two ions
    select = c(type, cycle, v44.mV, v45.mV),
    # include the Analysis number fron the file info and rename it to 'run'
    include_file_info = c(run = Analysis)
  ) %>%
  # look at first few records only
  head(n=10) %>% knitr::kable()
```


# Data Processing

The isoreader package is intended to make raw stable isotope data easily accessible. However, as with most analytical data, there is significant downstream processing required to turn these raw signal intensities into properly referenced isotopic measurement. This and similar functionality as well as data visualization is part of the [isoprocessor package](https://isoprocessor.isoverse.org) which takes isotopic data through the various corrections in a transparent, efficient and reproducible manner.

That said, most vendor software also performs some of these calculations and it can be useful to be able to compare new data reduction procedures against those implemented in the vendor software. For this purpose, isoreader retrieves vendor computed data tables whenever possible, as illustrated below.

## Vendor Data Table

As with most data retrieval functions, the `iso_get_vendor_data_table()` function also allows specific column selection (by default, all columns are selected) and easy addition of file information via the `include_file_info` parameter (by default, none is included).

```{r}
# entire vendor data table
di_files %>% iso_get_vendor_data_table() %>% knitr::kable()
# get specific parts and add some file information
di_files %>%
  iso_get_vendor_data_table(
    # select cycle and all carbon columns
    select = c(cycle, matches("C")),
    # include the Identifier 1 fron the file info and rename it to 'id'
    include_file_info = c(id = `Identifier 1`)
  ) %>% knitr::kable()
```

# For expert users: retrieving all data

For users familiar with the nested data frames from the [tidyverse](https://www.tidyverse.org/) (particularly [tidyr](https://tidyr.tidyverse.org/)'s `nest` and `unnest`), there is an easy way to retrieve all data from the iso file objects in a single nested data frame:

```{r}
all_data <- di_files %>% iso_get_all_data()
# not printed out because this data frame is very big
```

# Saving collections

Saving entire collections of isofiles for retrieval at a later point is easily done using the `iso_save` function which stores collections or individual isoreader file objects in the efficient R data storage format `.rds` (if not specified, the extension `.di.rds` will be automatically appended). These saved collections can be conveniently read back using the same `iso_read_dual_inlet` command used for raw data files.

```{r}
# export to R data archive
di_files %>% iso_save("di_files_export.di.rds")

# read back the exported R data storage
iso_read_dual_inlet("di_files_export.di.rds")
```


# Data Export

At the moment, isoreader supports export of all data to Excel and the [Feather file format](https://blog.rstudio.com/2016/03/29/feather/) (a Python/R cross-over format). Note that both export methods have similar syntax and append the appropriate file extension for each type of export file (`.di.xlsx` and `.di.feather`, respectively).

```{r, eval = FALSE}
# export to excel
di_files %>% iso_export_to_excel("di_files_export")

# data sheets available in the exported data file:
readxl::excel_sheets("di_files_export.di.xlsx")
```

```{r, eval=FALSE}
# export to feather
di_files %>% iso_export_to_feather("di_files_export")

# exported feather files
list.files(pattern = ".di.feather")
```
---
title: "Operations"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
  html_document:
    code_folding: show
    df_print: paged
    number_sections: yes
    toc: yes
    toc_depth: 3
    toc_float: yes
editor_options:
  chunk_output_type: console
vignette: >
  %\VignetteIndexEntry{Operations}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction

Isoreader provides a number of general purpose operations that work on all supported IRMS data formats such as caching of files, parallel processing and catching read errors. This vignette demonstrates some of these general operations.

```{r, message=FALSE}
# load isoreader package
library(isoreader)
```

# Supported file types

```{r}
# list all suported file types
iso_get_supported_file_types() %>%
  dplyr::select(extension, software, description, type) %>%
  knitr::kable()
```

# Messages

By default, isoreader is quite verbose to let the user know what is happening. However, most functions can be silenced by adding the parameter `quiet = TRUE` to the function call. This can also be done globally using `iso_turn_info_messages_off()`

```{r}
# read a file in the default verbose mode
iso_get_reader_example("dual_inlet_example.did") %>%
  iso_read_dual_inlet() %>%
  iso_select_file_info(file_datetime, `Identifier 1`) %>%
  iso_get_file_info() %>%
  knitr::kable()

# read the same file but make the read process quiet
iso_get_reader_example("dual_inlet_example.did") %>%
  iso_read_dual_inlet(quiet = TRUE) %>%
  iso_select_file_info(file_datetime, `Identifier 1`) %>%
  iso_get_file_info() %>%
  knitr::kable()

# read the same file but turn all isoreader messages off
iso_turn_info_messages_off()
iso_get_reader_example("dual_inlet_example.did") %>%
  iso_read_dual_inlet(quiet = TRUE) %>%
  iso_select_file_info(file_datetime, `Identifier 1`) %>%
  iso_get_file_info() %>%
  knitr::kable()

# turn message back on
iso_turn_info_messages_on()
```

# Caching

By default, isoreader caches files as R objects to make access faster in the future. This feature can be turned off if you want to force a fresh read from the source file. Alternatively, you can clear the entire isoreader cache in your working directory to clean up previous file reads.

```{r}
# cleanup reader cache
iso_cleanup_reader_cache()

# read a new file (notice the time elapsed)
cf_file <- iso_get_reader_example("continuous_flow_example.dxf") %>%
  iso_read_continuous_flow()

# re-read the same file much faster (it will be read from cache)
cf_file <- iso_get_reader_example("continuous_flow_example.dxf") %>%
    iso_read_continuous_flow()

# turn reader caching off
iso_turn_reader_caching_off()

# re-read the same file (it will NOT be read from cache)
cf_file <- iso_get_reader_example("continuous_flow_example.dxf") %>%
  iso_read_continuous_flow()

# turn reader caching back on
iso_turn_reader_caching_on()
```

# Parallel processing

Isoreader supports parallel processing of data files based on the number of processors available in a computer simply by setting the `parallel = TRUE` flag in any file read operation. This makes it possible to read large quantities of data files much more quickly on a multi-core system (i.e. most modern laptops).

However, whether parallel processing yields significant improvements in read speeds depends on the number of available processors, file types and operating system. In theory, parallel processing always reduces computation time but in practice this is offset by various factors including the size of the data that needs to be sent back and forth between the processors, file system read/write speed, and the spin-up time for new processes. Generally speaking, parallel processing can provide significant improvements in speed with larger number of files (~10+) and more complex read operations (e.g. continuous flow > dual inlet > scan file). Reading from cache is so efficient that there are rarely gains from parallel processing and it is usually faster NOT to read in parallel once a set of files is already cached.

```{r}
# read 3 files in parallel (note that this is usually not a large enough file number to be worth it)
di_files <-
  iso_read_dual_inlet(
    iso_get_reader_example("dual_inlet_example.did"),
    iso_get_reader_example("dual_inlet_example.caf"),
    iso_get_reader_example("dual_inlet_nu_example.txt"),
    nu_masses = 49:44,
    parallel = TRUE
  )
```


# Combining / subsetting isofiles

All isoreader objects are lists that can be combined or subset to work with only specific files or create a larger collection.

```{r}
# all 3 di_files read above
di_files

# only one of the files (by index)
di_files[[2]]

# only one of the files (by file_id)
di_files$dual_inlet_example.did

# a subset of the files (by index)
di_files[c(1,3)]

# a subset of the files (by file_id)
di_files[c("dual_inlet_example.did", "dual_inlet_example.caf")]

# same result using iso_filter_files (more flexible + verbose output)
di_files %>% iso_filter_files(
  file_id %in% c("dual_inlet_example.did", "dual_inlet_example.caf")
)

# recombining subset files
c(
  di_files[3],
  di_files[1]
)
```

# Dealing with file read problems

Isoreader is designed to catch problems during file reading without crashing the read pipeline. It keeps track of all problems encountered along the way to make it easy to see what went wrong and remove erroneous files. Most times, files that were only partly saved because of an interrupted instrument analysis will have errors. If you encounter a file that should have intact data in it but has an error in isoreader, please file a bug report and submit your file at https://github.com/isoverse/isoreader/issues

```{r}
# read two files, one of which is erroneous
iso_files <-
  iso_read_continuous_flow(
    iso_get_reader_example("continuous_flow_example.dxf"),
    system.file("errdata", "cf_without_data.dxf", package = "isoreader")
  )

# retrieve problem summary
iso_files %>% iso_get_problems_summary() %>% knitr::kable()

# retrieve problem details
iso_files %>% iso_get_problems() %>% knitr::kable()

# filter out erroneous files
iso_files <- iso_files %>% iso_filter_files_with_problems()
```

# Re-reading files

If a file has changed (e.g. is edited through the vendor software) and the changes should be loaded in isoreader, it is easy to re-read and update just those files within a file collection by using the `iso_reread_changed_files()` function. If some of the files are no longer accessible at their original location, it will throw a warning. If the location for all files has changed, it can be easily adjusted by modifying the `file_root` file info parameter using `iso_set_file_root()`.

Similar functions can be used to re-read outdated files from an older isoreader version (`iso_reread_outdated_files()`), attempt to re-read problematic files that had read errors/warnings (`iso_reread_problem_files()`), or simply re-read all files in a collection (`iso_reread_all_files()`).

```{r}
# re-read the 3 dual inlet files from their original location if any have changed
di_files %>%
  iso_reread_changed_files()

# update the file_root for the files before re-read (in this case to a location
# that does not hold these files and hence will lead to a warning)
di_files %>%
  iso_set_file_root(root = ".") %>%
  iso_reread_all_files()
```

# Units

Isoreader provides a built in data type with units (`iso_with_units`) that can be used to easily keep track of units inside data frame. These units can be made explicit (=included in the column header), stripped altogether, or turned back to be implicit.

```{r}
# strip all units
cf_file %>%
  iso_get_vendor_data_table(select = c(`Ampl 28`, `rIntensity 28`, `d 15N/14N`)) %>%
  iso_strip_units() %>% head(3)

# make units explicit
cf_file %>%
  iso_get_vendor_data_table(select = c(`Ampl 28`, `rIntensity 28`, `d 15N/14N`)) %>%
  iso_make_units_explicit() %>% head(3)

# introduce new unit columns e.g. in the file info
cf_file %>%
  iso_mutate_file_info(weight = iso_with_units(0.42, "mg")) %>%
  iso_get_vendor_data_table(select = c(`Ampl 28`, `rIntensity 28`, `d 15N/14N`),
                            include_file_info = weight) %>%
  iso_make_units_explicit() %>% head(3)

# or turn a column e.g. with custom format units in the header into implicit units
cf_file %>%
  iso_mutate_file_info(weight.mg = 0.42) %>%
  iso_get_vendor_data_table(select = c(`Ampl 28`, `rIntensity 28`, `d 15N/14N`),
                            include_file_info = weight.mg) %>%
  iso_make_units_implicit(prefix = ".", suffix = "") %>% head(3)
```

# Formatting

Formatting data into text is easily achieved with the built in R function `sprintf` but this package also provides a convenience function that knows how to incorporate units information from `iso_with_units` values. Use `iso_format` to format and concatenate any single values or entire columns inside a data frame.

```{r}
# concatenation example with single values
iso_format(
   pi = 3.14159,
   x = iso_with_units(42, "mg"),
   ID = "ABC",
   signif = 4,
   sep = " | "
)

# example inside a data frame
cf_file %>%
  iso_get_vendor_data_table(select = c(`Nr.`, `Ampl 28`, `d 15N/14N`)) %>%
  dplyr::select(-file_id) %>%
  head(3) %>%
  # introduce new label columns using iso_format
  dplyr::mutate(
    # default concatenation of values
    label_default = iso_format(
      `Nr.`, `Ampl 28`, `d 15N/14N`,
      sep = ", "
    ),
    # concatenate with custom names for each value
    label_named = iso_format(
      `#` = `Nr.`, A = `Ampl 28`, d15 = `d 15N/14N`,
      sep = ", "
    ),
    # concatenate just the values and increase significant digits
    label_value = iso_format(
      `Nr.`, `Ampl 28`, `d 15N/14N`,
      sep = ", ", format_names = NULL, signif = 6
    )
  )
```
---
title: "Development features of isoreader"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
  html_document:
    code_folding: show
    df_print: paged
    number_sections: yes
    toc: yes
    toc_depth: 3
    toc_float: yes
editor_options:
  chunk_output_type: console
vignette: >
  %\VignetteIndexEntry{Development features of isoreader}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(isoreader)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This vignette introduces some of the development features of the isoreader package and is aimed primarily at code contributors interested in expanding its functionality or helping with bug fixes.

# Adding new file format readers

Testing out new file format readers is easiest by registering a new reader function for a specific file extension using `iso_register_dual_inlet_file_reader` and `iso_register_continuous_flow_file_reader`, respectively. Both require an extension (e.g. `".ext"`), name of the new reader function (`"new_reader"`), and optionally a description. Both functions automatically return a data frame with a list of all registered reader. Overwriting of existing readers with a different function requires an explicit `overwrite = TRUE` flag. All reader functions must accept an isoreader data structure object (`ds`) as the first argument, a list of reader specific options as the second argument (`options`), and should return the structure with data filled in for downstream isoreader operations to work smoothly. The following minimal example illustrates how to do this with the `new_reader` function simply printing out the layout of the provided data structure skeleton `ds`.

```{r}
new_reader <- function(ds, options = list()) {
  isoreader:::log_message("this is the new reader!")
  str(ds)
  return(ds)
}

# register new reader
readers <- iso_register_dual_inlet_file_reader(".new.did", "new_reader")
knitr::kable(readers)

# copy an example file from the package with the new extension
iso_get_reader_example("dual_inlet_example.did") %>% file.copy(to = "example.new.did")

# read the file
iso_read_dual_inlet("example.new.did", read_cache = FALSE)
file.remove("example.new.did")
```

Note that for parallel processing to work during the read process (`parallel = TRUE`), isoreader needs to know where to find the new reader function. It will figure this out automatically as long as the function name is unique but if this fails (or to be on the safe side), please specify e.g. `env = "R_GlobalEnv"` or `env = "newpackage"` during the reader registration. Also note that isoreader will not automatically know where to find all functions called from within the new reader function if they are not part of base R and it is recommended to make all outside calls explicit (e.g. `dplyr::filter(...)`) to preempt this potential problem. For info messages and warnings to work with the progress bar and in parallel reads, make sure to use `isoreader:::log_message(...)` and `isoreader:::log_warning(...)` instead of base R's `message(...)` and `warning(...)`.

If you have designed and tested a new reader, please consider contributing it to the `isoreader` github repository via pull request.

# Processing hooks

Isoreader defines two processing hooks at the beginning and end of reading an individual file. This is useful for integration into pipelines that require additional output (such as GUIs) but is also sometimes useful for debugging purposes. The expressions are evaluated in the context of the `isoreader:::read_iso_file` function and have access to all parameters passed to this function, such as e.g. `file_n` and `path`. Same as for new readers: for info messages and warnings to work with the progress bar and in parallel reads, make sure to use `isoreader:::log_message(...)` and `isoreader:::log_warning(...)` instead of base R's `message(...)` and `warning(...)`. The main difference between the two is that `log_message()` will honor the `quiet = TRUE` flag passed to the main `iso_read...()` call whereas `log_warning()` will always show its message no matter the `quiet` setting.

```{r}
isoreader:::set_read_file_event_expr({
  isoreader:::log_message(sprintf("starting file #%.d, named '%s'", file_n, basename(path)))
})
isoreader:::set_finish_file_event_expr({
  isoreader:::log_message(sprintf("finished file #%.d", file_n))
})

c(
  iso_get_reader_example("dual_inlet_example.did"),
  iso_get_reader_example("dual_inlet_example.caf")
) %>% iso_read_dual_inlet(read_cache = FALSE)

isoreader:::initialize_options() # reset all isoreader options
```


# Debugging isoreader

The best way to start debugging an isoreader call is to switch the package into debug mode. This is done using the internal `iso_turn_debug_on()` function. This enables debug messages, turns caching off by default so files are always read anew, and makes the package keep more information in the isofile objects. It continues to catch errors inside file readers (keeping track of them in the [problems](operations.html#dealing-with-file-read-problems)) unless you set `iso_turn_debug_on(catch_errors = FALSE)`, in which case no errors are caught and stop the processing so you get the full traceback and debugging options of your IDE.

## Debugging binary file reads (Isodat)

Errors during the binary file reads usually indicate the approximate position in the file where the error was encountered. The easiest way to get started on figuring out what the file looks like at that position is to use a binary file editor and jump to the position. For a sense of the interpreted structure around that position, one can use the internal function `map_binary_structure` which tries to apply all frequently occurring binary patterns recognized by isoreader. The binary representation of the source file is only available if in debug mode but if debug mode is ON, it can be accessed as follows:

```{r}
# turn on debug mode
isoreader:::iso_turn_debug_on()
# read example file
ex <- iso_get_reader_example("dual_inlet_example.did") %>%  
  iso_read_dual_inlet(quiet = TRUE)
# access binary
bin <- ex$binary
# use structure mapping
bin %>%
  isoreader:::move_to_pos(1340) %>%
  isoreader:::map_binary_structure(length = 200)
```

This structure representation shows recognized control elements in `<...>` and data elements in `{...}` which are converted to text or numeric representation if the interpretation is unambiguous, or plain hexadecimal characters if the nature of the data cannot be determined with certainty. Because this function tries all possible control elements and data interpretations, it is quite slow and may take a while if run for large stretches of binary code (i.e. if the `length` parameter is very long).

For an overview of all the control elements that are currently consider, use the internal `get_ctrl_blocks_config_df()` function.

```{r}
isoreader:::get_ctrl_blocks_config_df()
```

Additional information can be gleaned from the so-called control blocks, which are larger structural elements of Isodat binary files and are kept in a data frame within the binary object (again only available in debug mode).

```{r}
bin$C_blocks
```

Same as for specific byte positions, one can use the control blocks to navigate the file and `map_binary_structure`.

```{r}
bin %>%
  isoreader:::move_to_C_block("CMethod") %>%
  isoreader:::map_binary_structure(length = 200)
```
---
title: "Scan Examples"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
  html_document:
    code_folding: show
    df_print: paged
    number_sections: yes
    toc: yes
    toc_depth: 3
    toc_float: yes
editor_options:
  chunk_output_type: console
vignette: >
  %\VignetteIndexEntry{Scan Examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r setup, include = FALSE}
# global knitting options for code rendering
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")
```

# Introduction

Isoreader supports several dual inlet IRMS data formats. This vignette shows some of the functionality for scan data files. For additional information on operations more generally (caching, combining read files, data export, etc.), please consult the [operations vignette](https://isoreader.isoverse.org/articles/operations.html). For details on downstream data processing and visualization, see the [isoprocessor package](https://isoprocessor.isoverse.org).

Note: this vignette is still a work in progress.


```{r, message=FALSE}
# load isoreader package
library(isoreader)
```


# Reading files

Reading scan files is as simple as passing one or multiple file or folder paths to the `iso_read_scan()` function. If folders are provided, any files that have a recognized scan file extensions within those folders will be processed (e.g. all `.scn`). Here we read several files that are bundled with the package as examples (and whose paths can be retrieved using the `iso_get_reader_example()` function).

```{r, message=FALSE}
# all available examples
iso_get_reader_examples() %>% knitr::kable()
```

```{r}
# read scan examples
scan_files <-
  iso_read_scan(
    iso_get_reader_example("peak_shape_scan_example.scn"),
    iso_get_reader_example("background_scan_example.scn"),
    iso_get_reader_example("full_scan_example.scn"),
    iso_get_reader_example("time_scan_example.scn")
  )
```

# File summary

The `scan_files` variable now contains a set of isoreader objects, one for each file. Take a look at what information was retrieved from the files using the `iso_get_data_summary()` function.

```{r}
scan_files %>% iso_get_data_summary() %>% knitr::kable()
```

## Problems

In case there was any trouble with reading any of the files, the following functions provide an overview summary as well as details of all errors and warnings, respectively. The examples here contain no errors but if you run into any unexpected file read problems, please file a bug report in the [isoreader issue tracker](https://github.com/isoverse/isoreader/issues).

```{r}
scan_files %>% iso_get_problems_summary() %>% knitr::kable()
scan_files %>% iso_get_problems() %>% knitr::kable()
```

# File Information

Detailed file information can be aggregated for all isofiles using the `iso_get_file_info()` function which supports the full [select syntax](https://dplyr.tidyverse.org/reference/select.html) of the [dplyr](https://dplyr.tidyverse.org/) package to specify which columns are of interest (by default, all file information is retrieved).

```{r}
# all file information
scan_files %>% iso_get_file_info(select = c(-file_root)) %>% knitr::kable()
```

## Select/Rename

File information can also be modified across an entire collection of isofiles using the `iso_select_file_info()` and `iso_rename_file_info()` functions:

```{r}
# select + rename specific file info columns
scan_files2 <- scan_files %>%
  iso_select_file_info(-file_root) %>%
  iso_rename_file_info(`Date & Time` = file_datetime)

# fetch all file info
scan_files2 %>% iso_get_file_info() %>% knitr::kable()
```


## Filter

Any collection of isofiles can also be filtered based on the available file information using the function `iso_filter_files`. This function can operate on any column available in the file information and supports full [dplyr](https://dplyr.tidyverse.org/reference/filter.html) syntax.

```{r}
# find files that have 'CIT' in the new ID field
scan_files2 %>%
  iso_filter_files(type == "High Voltage") %>%
  iso_get_file_info() %>%
  knitr::kable()
```

## Mutate

The file information in any collection of isofiles can also be mutated using the function `iso_mutate_file_info`. This function can introduce new columns and operate on any existing columns available in the file information (even if it does not exist in all files) and supports full [dplyr](https://dplyr.tidyverse.org/reference/mutate.html) syntax.

```{r}
scan_files3 <- scan_files2 %>%
  iso_mutate_file_info(
    # introduce new column
    `Run in 2019?` = `Date & Time` > "2019-01-01" & `Date & Time` < "2020-01-01"
  )

scan_files3 %>%
  iso_get_file_info() %>%
  knitr::kable()
```

# Resistors

Additionally, some IRMS data files contain resistor information that are useful for downstream calculations (see e.g. section on signal conversion later in this vignette):

```{r}
scan_files %>% iso_get_resistors() %>% knitr::kable()
```

# Raw Data

The raw data read from the scan files can be retrieved similarly using the `iso_get_raw_data()` function. Most data aggregation functions also allow for inclusion of file information using the `include_file_info` parameter, which functions identically to the `select` parameter of the `iso_get_file_info` function discussed earlier.

```{r}
# get raw data with default selections (all raw data, no additional file info)
scan_files %>% iso_get_raw_data() %>% head(n=10) %>% knitr::kable()
# get specific raw data and add some file information
scan_files %>%
  iso_get_raw_data(
    # select just time and the two ions
    select = c(x, x_units, v44.mV, v45.mV),
    # include the scan type and rename the column
    include_file_info = c(`Scan Type` = type)
  ) %>%
  # look at first few records only
  head(n=10) %>% knitr::kable()
```

# For expert users: retrieving all data

For users familiar with the nested data frames from the [tidyverse](https://www.tidyverse.org/) (particularly [tidyr](https://tidyr.tidyverse.org/)'s `nest` and `unnest`), there is an easy way to retrieve all data from the iso file objects in a single nested data frame:

```{r}
all_data <- scan_files %>% iso_get_all_data()
# not printed out because this data frame is very big
```


# Saving collections

Saving entire collections of isofiles for retrieval at a later point is easily done using the `iso_save` function which stores collections or individual isoreader file objects in the efficient R data storage format `.rds` (if not specified, the extension `.scan.rds` will be automatically appended). These saved collections can be conveniently read back using the same `iso_read_scan` command used for raw data files.

```{r}
# export to R data archive
scan_files %>% iso_save("scan_files_export.scan.rds")

# read back the exported R data storage
iso_read_scan("scan_files_export.scan.rds")
```

# Data Export

At the moment, isoreader supports export of all data to Excel and the [Feather file format](https://blog.rstudio.com/2016/03/29/feather/) (a Python/R cross-over format). Note that both export methods have similar syntax and append the appropriate file extension for each type of export file (`.scan.xlsx` and `.scan.feather`, respectively).

```{r, eval=FALSE}
# export to excel
scan_files %>% iso_export_to_excel("scan_files_export")

# data sheets available in the exported data file:
readxl::excel_sheets("scan_files_export.scan.xlsx")
```

```{r, eval=FALSE}
# export to feather
scan_files %>% iso_export_to_feather("scan_files_export")

# exported feather files
list.files(pattern = ".scan.feather")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_info_operations.R
\name{iso_parse_file_info}
\alias{iso_parse_file_info}
\title{Parse file info}
\usage{
iso_parse_file_info(
  iso_files,
  number = c(),
  double = c(),
  integer = c(),
  logical = c(),
  datetime = c(),
  text = c(),
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{number}{dplyr-style \link[dplyr]{select} condition to choose columns that should be converted to a number using \link[readr:parse_atomic]{parse_number}. Use \code{c(...)} to select multiple columns.}

\item{double}{dplyr-style \link[dplyr]{select} condition to choose columns that should be converted to a double using \link[readr:parse_atomic]{parse_double}. Use \code{c(...)} to select multiple columns.}

\item{integer}{dplyr-style \link[dplyr]{select} condition to choose columns that should be converted to an integer using \link[readr:parse_atomic]{parse_integer}. Use \code{c(...)} to select multiple columns.}

\item{logical}{dplyr-style \link[dplyr]{select} condition to choose columns that should be converted to a boolean (TRUE/FALSE) using \link[readr:parse_atomic]{parse_logical}. Use \code{c(...)} to select multiple columns.}

\item{datetime}{dplyr-style \link[dplyr]{select} condition to choose columns that should be converted to a date-time using \link[readr:parse_atomic]{parse_datetime}. Use \code{c(...)} to select multiple columns.}

\item{text}{dplyr-style \link[dplyr]{select} condition to choose columns that should be converted to text using \link[base]{as.character}. Use \code{c(...)} to select multiple columns.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Convenience function to batch parse file info (\code{\link{iso_get_file_info}}) columns in isofile objects for the most common parsing calls. Uses the \code{parse_} functions exported from \link{readr} and described in \link{extract_data}. Note that for less common parsing calls or calls that require additional parameters to the parsing function, it is better to parse columns one-by-one using \code{\link{iso_mutate_file_info}} instead.
}
\seealso{
Other file_info operations: 
\code{\link{iso_add_file_info.iso_file_list}()},
\code{\link{iso_filter_files}()},
\code{\link{iso_mutate_file_info}()},
\code{\link{iso_rename_file_info}()},
\code{\link{iso_select_file_info}()},
\code{\link{iso_set_file_root}()}
}
\concept{file_info operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_info_operations.R
\name{iso_select_file_info}
\alias{iso_select_file_info}
\title{Select file info columns}
\usage{
iso_select_file_info(
  iso_files,
  ...,
  file_specific = FALSE,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{...}{dplyr-style \link[dplyr]{select} conditions applied based on each file's file_info (see \code{\link{iso_get_file_info}}). Note that the \code{file_id} column will always be kept, no matter the selection criteria, and cannot be renamed to protect from unexpected behavior.}

\item{file_specific}{whether to run the select criteria (\code{...}) specifically within each individual file rather than on all files jointly. This is a lot slower but makes it possible to  select different columns in different iso_files depending on what exists in each file and is mostly of use when working with data from multiple instruments.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Select which file info columns (\code{\link{iso_get_file_info}}) to keep within isofile objects. Works just like dplyr's \link[dplyr]{select} and can rename columns on-the-fly. You can also use \link[dplyr]{select} directly but it will not provide summary information on the operation. To rename columns without removing all other information, use \link{iso_rename_file_info} instead. Set \code{file_specific = TRUE} to select different columns in different iso_files depending on what exists in each file. This is very useful when working with data from multiple instruments that may have the same information (e.g. sample name) stored in different columns.
}
\seealso{
Other file_info operations: 
\code{\link{iso_add_file_info.iso_file_list}()},
\code{\link{iso_filter_files}()},
\code{\link{iso_mutate_file_info}()},
\code{\link{iso_parse_file_info}()},
\code{\link{iso_rename_file_info}()},
\code{\link{iso_set_file_root}()}
}
\concept{file_info operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{iso_plot_dual_inlet_data}
\alias{iso_plot_dual_inlet_data}
\title{moved to isoprocessor}
\usage{
iso_plot_dual_inlet_data(...)
}
\arguments{
\item{...}{deprecated}
}
\description{
moved to isoprocessor
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/problems.R
\name{iso_filter_files_with_problems}
\alias{iso_filter_files_with_problems}
\title{Filter out problematic files}
\usage{
iso_filter_files_with_problems(
  iso_files,
  remove_files_with_errors = TRUE,
  remove_files_with_warnings = FALSE,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{remove_files_with_errors}{whether to remove files with errors (default is TRUE)}

\item{remove_files_with_warnings}{whether to remove files with warnings (default is FALSE)}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Use this function to filter out files that have encountered problems, either errors, warnings or both and returns the remaining iso_files. For additional functions available to check for and deal with problems, see the \link{iso_problem_functions}.
}
\seealso{
Other problem functions: 
\code{\link{iso_get_problems_summary}()},
\code{\link{iso_get_problems}()},
\code{\link{iso_has_problems}()},
\code{\link{iso_problem_functions}}
}
\concept{problem functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{vec_cast.iso_double_with_units}
\alias{vec_cast.iso_double_with_units}
\title{vec_cast for iso_double_with_units}
\usage{
\method{vec_cast}{iso_double_with_units}(x, to, ...)
}
\arguments{
\item{x}{Vectors to cast.}

\item{to}{Type to cast to. If \code{NULL}, \code{x} will be returned as is.}

\item{...}{For \code{vec_cast_common()}, vectors to cast. For
\code{vec_cast()}, \code{vec_cast_default()}, and \code{vec_restore()}, these
dots are only for future extensions and should be empty.}
}
\description{
vec_cast for iso_double_with_units
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_binary_files.R
\name{map_binary_structure}
\alias{map_binary_structure}
\title{Map isodat file binary structure.}
\usage{
map_binary_structure(
  bfile,
  length = 100,
  start = bfile$pos,
  ctrl_blocks = get_ctrl_blocks_config()
)
}
\arguments{
\item{bfile}{the binary file, stored in each iso_file under \code{$binary} if (and only if) the file was read with \link{iso_turn_debug_on} activated before.}

\item{length}{how many bytes to map}

\item{start}{at which byte position to start mapping (index 1 based)}

\item{ctrl_blocks}{named list of block patterns with size, regexp and [optional] replace function}
}
\description{
Map out binary structure for easy visualization (used mostly for error messages and debugging). See the development vignette for details and example application.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_data}
\alias{iso_get_data}
\title{DEPRECATED}
\usage{
iso_get_data(...)
}
\arguments{
\item{...}{forwarded to \link{iso_get_all_data}}
}
\description{
Please use \link{iso_get_all_data} instead.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{iso_read_dual_inlet}
\alias{iso_read_dual_inlet}
\title{Load dual inlet data}
\usage{
iso_read_dual_inlet(
  ...,
  root = ".",
  read_raw_data = default(read_raw_data),
  read_file_info = default(read_file_info),
  read_method_info = default(read_method_info),
  read_vendor_data_table = default(read_vendor_data_table),
  nu_masses = c(),
  discard_duplicates = TRUE,
  parallel = FALSE,
  parallel_plan = future::multisession,
  parallel_cores = future::availableCores(),
  cache = default(cache),
  read_cache = default(cache),
  reread_outdated_cache = FALSE,
  quiet = default(quiet),
  cache_files_with_errors = TRUE
)
}
\arguments{
\item{...}{one or multiple file/folder paths. All files must have a supported file extension. All folders are expanded and searched for files with supported file extensions (which are then included in the read).}

\item{root}{root directory for the isofiles. Can be relative to the current working directory (e.g. \code{"data"}) or an absolute path on the file system (e.g. \code{"/Users/..."} or \code{"C:/Data/.."}). The default is the current working directory (\code{"."}). Can be supplied as a vector of same length as the provided paths if the paths have different roots.}

\item{read_raw_data}{whether to read the raw mass/ion data from the file}

\item{read_file_info}{whether to read auxiliary file information (file id, sequence information, etc.)}

\item{read_method_info}{whether to read methods information (standards, processing info)}

\item{read_vendor_data_table}{whether to read the vendor computed data table}

\item{nu_masses}{list of masses (e.g. \code{c("46","45","44")}) to map the collector channels (interpreted in order, i.e. the first channel will be linked to the first mass, the second channel to the second mass, etc.). This parameter is only used for reading Nu data files.}

\item{discard_duplicates}{whether to automatically discard files with duplicate file IDs (i.e. duplicate file names). If \code{TRUE} (the default), only the first files are kept and any files with the same file ID are discarded. If \code{FALSE}, all duplicate files are kept but their file IDs are appended with suffix \code{#1}, \code{#2}, etc.}

\item{parallel}{whether to process in parallel based on the number of available CPU cores. This may yield performance increases for files that are slow to parse such as continuous flow isodat files but usually provides little benefit for efficient data formats such as reading from R Data Archives.}

\item{parallel_plan}{which parallel processing strategy to use, see \link[future]{plan}, typically \code{future::multisession} for compatibility with RStudio interactive mode. If supported by the operating system and running in detached mode (not interactively in RStudio) can also use \code{future::multicore}.}

\item{parallel_cores}{how many processor cores to use for parallel processing. By default the maximum available number of cores (\link[future]{availableCores}), which will allow maximal processing speed but may slow other programs running on your machine. Choose a smaller number if you want some processing resources to remain available for other processes. Will issue a warning if too many cores are requested and reset to the maximum available.}

\item{cache}{whether to cache iso_files. Note that R Data Storage files (.rds, see \link{iso_save}) are never cached since they are already essentially in cached form.}

\item{read_cache}{whether to reload from cache if a cached version exists. Note that it will only read from cache if the raw data file has not been modified since. Files that have been modified on disc (e.g. edited in the vendor software) will always be read anew. To automatically reread cached files that were cached by an outdated version of the isoreader package, set the \code{reread_outdated_cache} flag.}

\item{reread_outdated_cache}{whether to re-read outdated cache files whenever they are encountered.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}

\item{cache_files_with_errors}{deprecated. Please use \link{iso_reread_problem_files} instead to selectively re-read all files in a collection of iso files that had been previously read with errors or warnings.}
}
\description{
Load dual inlet data
}
\seealso{
Other isoread functions for different types of IRMS data: 
\code{\link{iso_read_continuous_flow}()},
\code{\link{iso_read_scan}()}
}
\concept{isoread functions for different types of IRMS data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.R
\docType{package}
\name{isoreader-package}
\alias{isoreader}
\alias{isoreader-package}
\title{isoreader: Read Stable Isotope Data Files}
\description{
Interface to the raw data file formats commonly encountered in scientific disciplines that make use of stable isotopes.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/isoverse/isoreader}
  \item Report bugs at \url{https://github.com/isoverse/isoreader/issues}
}

}
\author{
\strong{Maintainer}: Sebastian Kopf \email{sebastian.kopf@colorado.edu} (\href{https://orcid.org/0000-0002-2044-0201}{ORCID})

Authors:
\itemize{
  \item Brett Davidheiser-Kroll \email{brett.davidheiserkroll@colorado.edu} (\href{https://orcid.org/0000-0002-6153-7851}{ORCID})
  \item Ilja Kocken \email{i.j.kocken@uu.nl} (\href{https://orcid.org/0000-0003-2196-8718}{ORCID})
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{isoread}
\alias{isoread}
\title{Read isotope data file}
\usage{
isoread(...)
}
\arguments{
\item{...}{original isoread parameters}
}
\description{
This function from the original isoread package is deprecated, please use \link{iso_read_dual_inlet}, \link{iso_read_continuous_flow} and \link{iso_read_scan} instead.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/problems.R
\name{iso_omit_files_with_problems}
\alias{iso_omit_files_with_problems}
\title{Renamed to iso_filter_files_with_problems}
\usage{
iso_omit_files_with_problems(...)
}
\arguments{
\item{...}{deprecated}
}
\description{
This function has been renamed to \link{iso_filter_files_with_problems} for naming consistency.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_info_operations.R
\name{iso_mutate_file_info}
\alias{iso_mutate_file_info}
\title{Mutate file info}
\usage{
iso_mutate_file_info(iso_files, ..., quiet = default(quiet))
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{...}{dplyr-style \link[dplyr]{mutate} conditions applied to the combined file info (see \code{\link{iso_get_file_info}})}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Mutate the file info (\code{\link{iso_get_file_info}}) within isofile objects by changing existing columns or introducing new ones. Works just like dplyr's \link[dplyr]{mutate}. You can also use \link[dplyr]{mutate} directly but it will not provide summary information on the operation. Note that this will create missing columns that exist in some but not all of the passed in isofile objects in all isofile objects (filling them with NAs) the same way that \code{\link{iso_get_file_info}} does.
}
\seealso{
Other file_info operations: 
\code{\link{iso_add_file_info.iso_file_list}()},
\code{\link{iso_filter_files}()},
\code{\link{iso_parse_file_info}()},
\code{\link{iso_rename_file_info}()},
\code{\link{iso_select_file_info}()},
\code{\link{iso_set_file_root}()}
}
\concept{file_info operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_data_summary}
\alias{iso_get_data_summary}
\title{Get data summary}
\usage{
iso_get_data_summary(iso_files, quiet = default(quiet))
}
\arguments{
\item{iso_files}{single iso file or collection of iso_file objects}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\value{
a \code{\link[tibble]{tibble}} that summarizes the data in the \code{iso_files}
}
\description{
Summarize the data information from one or multiple iso files.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleanup.R
\name{extract_substring}
\alias{extract_substring}
\title{Extract a substring from text}
\usage{
extract_substring(
  string,
  pattern,
  capture_n = 1,
  capture_bracket = 0,
  missing = NA_character_
)
}
\arguments{
\item{string}{string to extract}

\item{pattern}{regular expression pattern to search for}

\item{capture_n}{within each string, which match of the \code{pattern} should be extracted? e.g. if the pattern searches for words, should the first, second or third word be captured?}

\item{capture_bracket}{for the captured match, which capture group should be extracted? i.e. which parentheses-enclosed segment of the \code{pattern}?
by default captures the whole pattern (\code{capture_bracket = 0}).}

\item{missing}{what to replace missing values with? Note that values can be missing because there are not enough captured matches or because the actual capture_bracket is empty.}
}
\value{
character vector of same length as \code{string} with the extracted substrings
}
\description{
This is a convenience function to capture substrings from textual data.
Uses \code{\link[stringr:str_match]{str_match_all}} internally but instead of returning everything, always returns only one single part of the match, depending on parameters \code{capture_n} and \code{capture_group}.
}
\seealso{
Other data extraction functions: 
\code{\link{extract_data}},
\code{\link{extract_word}()}
}
\concept{data extraction functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/settings.R
\name{iso_set_default_read_parameters}
\alias{iso_set_default_read_parameters}
\title{Set default read options}
\usage{
iso_set_default_read_parameters(
  data = NULL,
  read_raw_data,
  read_file_info,
  read_method_info,
  read_vendor_data_table,
  quiet = default(quiet)
)
}
\arguments{
\item{data}{a data frame - returned invisibly as is if provided (e.g. in the middle of a pipeline)}

\item{read_raw_data}{if provided, set as the default for `read_raw_data` parameters}

\item{read_file_info}{if provided, set as the default for `read_file_info` parameters}

\item{read_method_info}{if provided, set as the default for `read_method_info` parameters}

\item{read_vendor_data_table}{if provided, set as the default for `read_vendor_data_table` parameters}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Set default read options
}
\seealso{
Other settings functions: 
\code{\link{iso_caching}},
\code{\link{iso_get_default_reader_parameters}()},
\code{\link{iso_info_messages}}
}
\concept{settings functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleanup.R
\name{extract_word}
\alias{extract_word}
\title{Extract words from text}
\usage{
extract_word(
  string,
  capture_n = 1,
  include_numbers = TRUE,
  include_underscore = FALSE,
  include_dash = FALSE,
  include_space = FALSE,
  include_colon = FALSE,
  missing = NA_character_
)
}
\arguments{
\item{string}{string to extract}

\item{capture_n}{which word to extract? 1st, 2nd, 3rd?}

\item{include_numbers}{whether to include numbers (0-9) as part of the word (if FALSE, numbers will work as a word separator)}

\item{include_underscore}{whether to include the underscore character (_) as part of a word (if FALSE, it will work as a word separator)}

\item{include_dash}{whether to include the dash character (-) as part of a word (if FALSE, it will work as a word separator)}

\item{include_space}{whether to include the space character ( ) as part of a word (if FALSE, it will work as a word separator)}

\item{include_colon}{whether to include the colon character (.) as part of a word (if FALSE, it will work as a word separator)}

\item{missing}{what to replace missing values with? Note that values can be missing because there are not enough captured matches or because the actual capture_bracket is empty.}
}
\description{
This extracts words from text, by default looks for continuous sequences of numbers and/or letters.
Can adjust whether characters such as "_", "-", " ", and "." should be counted as part of a word or separate them and whether numbers should be included.
}
\examples{
x_text <- extract_word(c("sample number16.2", "sample number7b"),
                       capture_n = 2, include_colon = TRUE)
# "number16.2" "number7b"
x_num <- parse_number(x_text)
# 16.2 7.0
}
\seealso{
Other data extraction functions: 
\code{\link{extract_data}},
\code{\link{extract_substring}()}
}
\concept{data extraction functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{iso_format}
\alias{iso_format}
\title{Format values}
\usage{
iso_format(
  ...,
  signif = 3,
  format_names = "\%s: ",
  format_units = "\%s",
  replace_permil = TRUE,
  sep = "\\n"
)
}
\arguments{
\item{...}{variable names with data. Must have the same dimensions if multiple are supplied. Can be named to rename variable name output. Will include units in output for all \link{iso_with_units}.}

\item{signif}{number of significant digits for numbered data}

\item{format_names}{how to format the variable names, set to \code{NULL} to remove names}

\item{format_units}{how to format the units from \code{\link{iso_double_with_units}} variables, set to \code{NULL} to omit units}

\item{replace_permil}{whether to replace the term 'permil' with the permil symbol (\\u2030)}

\item{sep}{separator between variables if multiple are provided in \code{...}}
}
\description{
Convenience function to easily format and concatenate text and numeric values. Can be used with any test and number data. Automatically detects \code{\link{iso_with_units}} values and incorporates the units into the formatting.
}
\examples{
x <- iso_with_units(1:5, "V")
y <- iso_with_units(1:5, "permil")
iso_format(x, y)
iso_format(amplitude = x, d13C = y)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{iso_read_continuous_flow}
\alias{iso_read_continuous_flow}
\title{Load continuous flow data}
\usage{
iso_read_continuous_flow(
  ...,
  root = ".",
  read_raw_data = default(read_raw_data),
  read_file_info = default(read_file_info),
  read_method_info = default(read_method_info),
  read_vendor_data_table = default(read_vendor_data_table),
  discard_duplicates = TRUE,
  parallel = FALSE,
  parallel_plan = future::multisession,
  parallel_cores = future::availableCores(),
  cache = default(cache),
  read_cache = default(cache),
  reread_outdated_cache = FALSE,
  quiet = default(quiet),
  cache_files_with_errors = TRUE
)
}
\arguments{
\item{...}{one or multiple file/folder paths. All files must have a supported file extension. All folders are expanded and searched for files with supported file extensions (which are then included in the read).}

\item{root}{root directory for the isofiles. Can be relative to the current working directory (e.g. \code{"data"}) or an absolute path on the file system (e.g. \code{"/Users/..."} or \code{"C:/Data/.."}). The default is the current working directory (\code{"."}). Can be supplied as a vector of same length as the provided paths if the paths have different roots.}

\item{read_raw_data}{whether to read the raw mass/ion data from the file}

\item{read_file_info}{whether to read auxiliary file information (file id, sequence information, etc.)}

\item{read_method_info}{whether to read methods information (standards, processing info)}

\item{read_vendor_data_table}{whether to read the vendor computed data table}

\item{discard_duplicates}{whether to automatically discard files with duplicate file IDs (i.e. duplicate file names). If \code{TRUE} (the default), only the first files are kept and any files with the same file ID are discarded. If \code{FALSE}, all duplicate files are kept but their file IDs are appended with suffix \code{#1}, \code{#2}, etc.}

\item{parallel}{whether to process in parallel based on the number of available CPU cores. This may yield performance increases for files that are slow to parse such as continuous flow isodat files but usually provides little benefit for efficient data formats such as reading from R Data Archives.}

\item{parallel_plan}{which parallel processing strategy to use, see \link[future]{plan}, typically \code{future::multisession} for compatibility with RStudio interactive mode. If supported by the operating system and running in detached mode (not interactively in RStudio) can also use \code{future::multicore}.}

\item{parallel_cores}{how many processor cores to use for parallel processing. By default the maximum available number of cores (\link[future]{availableCores}), which will allow maximal processing speed but may slow other programs running on your machine. Choose a smaller number if you want some processing resources to remain available for other processes. Will issue a warning if too many cores are requested and reset to the maximum available.}

\item{cache}{whether to cache iso_files. Note that R Data Storage files (.rds, see \link{iso_save}) are never cached since they are already essentially in cached form.}

\item{read_cache}{whether to reload from cache if a cached version exists. Note that it will only read from cache if the raw data file has not been modified since. Files that have been modified on disc (e.g. edited in the vendor software) will always be read anew. To automatically reread cached files that were cached by an outdated version of the isoreader package, set the \code{reread_outdated_cache} flag.}

\item{reread_outdated_cache}{whether to re-read outdated cache files whenever they are encountered.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}

\item{cache_files_with_errors}{deprecated. Please use \link{iso_reread_problem_files} instead to selectively re-read all files in a collection of iso files that had been previously read with errors or warnings.}
}
\description{
Load continuous flow data
}
\seealso{
Other isoread functions for different types of IRMS data: 
\code{\link{iso_read_dual_inlet}()},
\code{\link{iso_read_scan}()}
}
\concept{isoread functions for different types of IRMS data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_raw_data}
\alias{iso_get_raw_data}
\title{Aggregate raw data}
\usage{
iso_get_raw_data(
  iso_files,
  select = everything(),
  gather = FALSE,
  include_file_info = NULL,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{select}{which data columns to select - use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax. By default, all columns are selected.}

\item{gather}{whether to gather raw data into long format (e.g. for ease of use in plotting). Not that the \code{select} parameter applies to the data columns BEFORE gathering.}

\item{include_file_info}{which file information to include (see \code{\link{iso_get_file_info}}). Use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Aggregate the raw ion data from the provided iso_files. Can aggregate either in a wide table (for easy overview) or a gathered long table (for plotting and further data processing). The raw data is only available if the iso_files were read with parameter \code{read_raw_data=TRUE}.
}
\seealso{
Other data retrieval functions: 
\code{\link{iso_get_all_data}()},
\code{\link{iso_get_bgrd_data}()},
\code{\link{iso_get_file_info}()},
\code{\link{iso_get_resistors}()},
\code{\link{iso_get_standards}()},
\code{\link{iso_get_vendor_data_table}()}

Other data retrieval functions: 
\code{\link{iso_get_all_data}()},
\code{\link{iso_get_bgrd_data}()},
\code{\link{iso_get_file_info}()},
\code{\link{iso_get_resistors}()},
\code{\link{iso_get_standards}()},
\code{\link{iso_get_vendor_data_table}()}
}
\concept{data retrieval functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_standards_info}
\alias{iso_get_standards_info}
\title{DEPRECATED}
\usage{
iso_get_standards_info(...)
}
\arguments{
\item{...}{forwarded to \link{iso_get_standards}}
}
\description{
Please use \link{iso_get_standards} instead.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{vec_ptype2.iso_double_with_units}
\alias{vec_ptype2.iso_double_with_units}
\title{vec_ptype2 for iso_double_with_units}
\usage{
\method{vec_ptype2}{iso_double_with_units}(x, y, ...)
}
\arguments{
\item{x}{Vector types.}

\item{y}{Vector types.}

\item{...}{These dots are for future extensions and must be empty.}
}
\description{
vec_ptype2 for iso_double_with_units
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{iso_read_files}
\alias{iso_read_files}
\title{Core function to read isotope data files}
\usage{
iso_read_files(
  paths,
  root,
  supported_extensions,
  data_structure,
  read_options = c(),
  reader_options = list(),
  discard_duplicates = TRUE,
  cache_files_with_errors = TRUE,
  parallel = FALSE,
  parallel_plan = future::multisession,
  parallel_cores = future::availableCores(),
  cache = default(cache),
  read_cache = default(cache),
  reread_outdated_cache = FALSE,
  quiet = default(quiet)
)
}
\arguments{
\item{paths}{one or multiple file/folder paths. All files must have a supported file extension. All folders are expanded and searched for files with supported file extensions (which are then included in the read). Paths can be absolute paths or relative to the provided file \code{root} (which is the current working directory by default). For absolute paths, a common root directory will be guessed using \link{iso_find_absolute_path_roots}. The root portion of paths will never be displayed in info messages.}

\item{root}{root directory for the isofiles. Can be relative to the current working directory (e.g. \code{"data"}) or an absolute path on the file system (e.g. \code{"/Users/..."} or \code{"C:/Data/.."}). The default is the current working directory (\code{"."}). Can be supplied as a vector of same length as the provided paths if the paths have different roots.}

\item{supported_extensions}{data frame with supported extensions and corresponding reader functions (columns 'extension', 'func', 'cacheable')}

\item{data_structure}{the basic data structure for the type of iso_file}

\item{read_options}{vector of read options to be stored in the data structure (e.g. \code{c(read_vendor_data_table = FALSE)}). The \code{read_} prefix is optional.}

\item{reader_options}{list of parameters to be passed on to the reader}

\item{discard_duplicates}{whether to automatically discard files with duplicate file IDs (i.e. duplicate file names). If \code{TRUE} (the default), only the first files are kept and any files with the same file ID are discarded. If \code{FALSE}, all duplicate files are kept but their file IDs are appended with suffix \code{#1}, \code{#2}, etc.}

\item{cache_files_with_errors}{deprecated. Please use \link{iso_reread_problem_files} instead to selectively re-read all files in a collection of iso files that had been previously read with errors or warnings.}

\item{parallel}{whether to process in parallel based on the number of available CPU cores. This may yield performance increases for files that are slow to parse such as continuous flow isodat files but usually provides little benefit for efficient data formats such as reading from R Data Archives.}

\item{parallel_plan}{which parallel processing strategy to use, see \link[future]{plan}, typically \code{future::multisession} for compatibility with RStudio interactive mode. If supported by the operating system and running in detached mode (not interactively in RStudio) can also use \code{future::multicore}.}

\item{parallel_cores}{how many processor cores to use for parallel processing. By default the maximum available number of cores (\link[future]{availableCores}), which will allow maximal processing speed but may slow other programs running on your machine. Choose a smaller number if you want some processing resources to remain available for other processes. Will issue a warning if too many cores are requested and reset to the maximum available.}

\item{cache}{whether to cache iso_files. Note that R Data Storage files (.rds, see \link{iso_save}) are never cached since they are already essentially in cached form.}

\item{read_cache}{whether to reload from cache if a cached version exists. Note that it will only read from cache if the raw data file has not been modified since. Files that have been modified on disc (e.g. edited in the vendor software) will always be read anew. To automatically reread cached files that were cached by an outdated version of the isoreader package, set the \code{reread_outdated_cache} flag.}

\item{reread_outdated_cache}{whether to re-read outdated cache files whenever they are encountered.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\value{
single iso_file object (if single file) or list of iso_files (iso_file_list)
}
\description{
This function takes care of extracting basic information about iso_files, dealing with problems and making sure only valid fire formats are processed.
This function is not typically called directly but indirectly by calling \link{iso_read_dual_inlet}, \link{iso_read_continuous_flow} and \link{iso_read_scan}.
It is made available outside the package because it can be very useful for testing new file readers.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isodata_structures.R
\name{print.iso_file_list}
\alias{print.iso_file_list}
\alias{print.iso_file}
\alias{print.dual_inlet}
\alias{print.continuous_flow}
\alias{print.scan}
\title{Isofile printing}
\usage{
\method{print}{iso_file_list}(x, ...)

\method{print}{iso_file}(x, ..., show_problems = TRUE)

\method{print}{dual_inlet}(x, ..., show_problems = TRUE)

\method{print}{continuous_flow}(x, ..., show_problems = TRUE)

\method{print}{scan}(x, ..., show_problems = TRUE)
}
\arguments{
\item{x}{Object to show.}

\item{...}{additional parameters passed to print.default}

\item{show_problems}{whether to show encountered problems}
}
\description{
Print summary of individual iso_files (dual inlet or continuous flow) or collection of iso_files.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{iso_make_units_implicit}
\alias{iso_make_units_implicit}
\title{Make units implicit}
\usage{
iso_make_units_implicit(df, prefix = " [", suffix = "]")
}
\arguments{
\item{df}{the data frame in which to make the units implicit/explicit}

\item{prefix}{the prefix for the units}

\item{suffix}{the suffix for the units}
}
\description{
This function is intended for data frames /tibbles only and tries to figure out which numeric columns have units in the column names and makes those units implicit using \code{\link{iso_double_with_units}}. The reverse function is \code{\link{iso_make_units_explicit}}.
}
\examples{
# generate implicit units
df <- tibble(peak = 1:5, `height [V]` = 1:5)
iso_make_units_implicit(df)

# convert back and forth
iso_make_units_implicit(df) \%>\% iso_make_units_explicit()

# implicit units from custom prefix & suffix
df <- tibble(peak = 1:5, height.V = 1:5)
iso_make_units_implicit(df, prefix = ".", suffix = "")
}
\seealso{
Other functions for values with units: 
\code{\link{iso_get_units}()},
\code{\link{iso_is_double_with_units}()},
\code{\link{iso_make_units_explicit}()},
\code{\link{iso_strip_units}()},
\code{\link{iso_with_units}()}
}
\concept{functions for values with units}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_info_operations.R
\name{iso_set_file_root}
\alias{iso_set_file_root}
\title{Set iso file directory root}
\usage{
iso_set_file_root(
  iso_files,
  root = ".",
  remove_embedded_root = NULL,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{root}{new root directory for the isofiles. Can be relative to the current working directory (e.g. \code{"data"}) or an absolute path on the file system (e.g. \code{"/Users/..."} or \code{"C:/Data/.."}). Can be supplied as a vector of same length as the \code{iso_files} if the files have different roots. Use \code{root = "."} to set the root to the current working directory (the default).}

\item{remove_embedded_root}{set this parameter to a root path that is embedded in the isofiles' \code{file_path}. Will warn about any paths that cannot be simplified by removing the specified \code{remove_embedded_root}.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Sets the root directory for a set of iso_files (property \code{file_root} in the file information), which is particularly useful for re-reading files (\link{reread_iso_files}) after they have changed location. Can optionally remove the previous root (\code{remove_embedded_root}) if it is still embedded in the isofiles' \code{file_path} instead of \code{file_root}. Will warn about any paths that cannot be simplified by removing the embedded root.
}
\seealso{
Other file_info operations: 
\code{\link{iso_add_file_info.iso_file_list}()},
\code{\link{iso_filter_files}()},
\code{\link{iso_mutate_file_info}()},
\code{\link{iso_parse_file_info}()},
\code{\link{iso_rename_file_info}()},
\code{\link{iso_select_file_info}()}
}
\concept{file_info operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleanup.R
\name{extract_data}
\alias{extract_data}
\title{Overview of text data extraction functions}
\description{
The following functions are intended to make it easy to extract relevant information from textual data.
These functions are primarily intended for use in \code{\link{iso_mutate_file_info}} and inside the filtering conditions passed to \code{\link{iso_filter_files}}. However, they can of course also be used stand-alone and in regular \code{\link[dplyr]{mutate}} or \code{\link[dplyr]{filter}} calls on the data frames returned by the data retrieval functions (\code{\link{iso_get_raw_data}}, \code{\link{iso_get_file_info}}, \code{\link{iso_get_vendor_data_table}}, etc.). Not that all the \code{parse_} functions are used in \code{\link{iso_parse_file_info}} for easy type conversions.
}
\details{
For simultaneous extraction of pure text data into multiple columns, please see the \code{\link[tidyr]{extract}} function from the \link{tidyr} package.

\itemize{
\item \code{\link{extract_substring}} is a generic convenience function to extract parts of textual data (based on regular expression matches).
Can be used in combination with the parsing functions to turn extracted substrings into numerical or logical data.

\item \code{\link{extract_word}} is a more specific convenience function to extract the 1st/2nd/3rd word from textual data.

\item \code{\link[readr:parse_atomic]{parse_number}} is a convenience function to extract a number even if it is surrounded by text (re-exported from the \link{readr} package).

\item \code{\link[readr:parse_atomic]{parse_double}} parses text that holds double (decimal) numerical values without any extraneous text around -
use \code{\link[readr:parse_atomic]{parse_number}} instead if this is not the case (re-exported from the \link{readr} package)

\item \code{\link[readr:parse_atomic]{parse_integer}} parses text that holds integer (whole number) numerical values without any extraneous text around -
use \code{\link[readr:parse_atomic]{parse_number}} instead if this is not the case (re-exported from the \link{readr} package)

\item \code{\link[readr:parse_atomic]{parse_logical}} parses text that holds logical (boolean, i.e. TRUE/FALSE) values (re-exported from the \link{readr} package)

\item \code{\link[readr:parse_atomic]{parse_datetime}} parses text that holds date and time information (re-exported from the \link{readr} package)

}
}
\seealso{
Other data extraction functions: 
\code{\link{extract_substring}()},
\code{\link{extract_word}()}
}
\concept{data extraction functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{iso_cleanup_reader_cache}
\alias{iso_cleanup_reader_cache}
\title{Cleanup cached files}
\usage{
iso_cleanup_reader_cache(all = FALSE)
}
\arguments{
\item{all}{deprecated}
}
\description{
Removes all cached files.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_bgrd_data}
\alias{iso_get_bgrd_data}
\title{Aggregate background data}
\usage{
iso_get_bgrd_data(
  iso_files,
  select = everything(),
  gather = FALSE,
  include_file_info = NULL,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{select}{which data columns to select - use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax. By default, all columns are selected.}

\item{gather}{whether to gather raw data into long format (e.g. for ease of use in plotting). Not that the \code{select} parameter applies to the data columns BEFORE gathering.}

\item{include_file_info}{which file information to include (see \code{\link{iso_get_file_info}}). Use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Aggregate the background data from the provided iso_files. Can aggregate either in a wide table (for easy overview) or a gathered long table (for plotting and further data processing). The background data is only available if the iso_files were read with parameter \code{read_raw_data=TRUE}.
}
\seealso{
Other data retrieval functions: 
\code{\link{iso_get_all_data}()},
\code{\link{iso_get_file_info}()},
\code{\link{iso_get_raw_data}()},
\code{\link{iso_get_resistors}()},
\code{\link{iso_get_standards}()},
\code{\link{iso_get_vendor_data_table}()}
}
\concept{data retrieval functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{iso_root_paths}
\alias{iso_root_paths}
\title{Root paths}
\usage{
iso_root_paths(path, root = ".", check_existence = TRUE)
}
\arguments{
\item{path}{vector of file/folder paths, mixed relative and absolute paths are allowed.}

\item{root}{root directory for the isofiles. Can be relative to the current working directory (e.g. \code{"data"}) or an absolute path on the file system (e.g. \code{"/Users/..."} or \code{"C:/Data/.."}). The default is the current working directory (\code{"."}). Can be supplied as a vector of same length as the provided paths if the paths have different roots.}

\item{check_existence}{whether to check for the existence of the paths}
}
\value{
a data frame with the root directories and paths relative to the root - order of input paths is preserved
}
\description{
Function to root both relative and absolute paths to a root directory (or directories) commonly relative to current working directory. Determines the best way to shorten relative paths and put absolute paths in a relative context (if possible) using \link{iso_shorten_relative_paths} and \link{iso_find_absolute_path_roots}, respectively.
}
\seealso{
Other file system functions: 
\code{\link{iso_expand_paths}()},
\code{\link{iso_find_absolute_path_roots}()},
\code{\link{iso_shorten_relative_paths}()}
}
\concept{file system functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{iso_reread_files}
\alias{iso_reread_files}
\alias{iso_reread_all_files}
\alias{iso_reread_changed_files}
\alias{iso_reread_outdated_files}
\alias{iso_reread_problem_files}
\alias{iso_reread_storage}
\alias{iso_reread_archive}
\title{Re-read iso_files}
\usage{
iso_reread_files(iso_files, ...)

iso_reread_all_files(
  iso_files,
  ...,
  stop_if_missing = FALSE,
  quiet = default(quiet)
)

iso_reread_changed_files(
  iso_files,
  ...,
  stop_if_missing = FALSE,
  quiet = default(quiet)
)

iso_reread_outdated_files(
  iso_files,
  ...,
  stop_if_missing = FALSE,
  quiet = default(quiet)
)

iso_reread_problem_files(
  iso_files,
  ...,
  stop_if_missing = FALSE,
  reread_files_with_errors = TRUE,
  reread_files_with_warnings = FALSE,
  quiet = default(quiet)
)

iso_reread_storage(...)

iso_reread_archive(...)
}
\arguments{
\item{iso_files}{collection of iso_files}

\item{...}{additional read parameters that should be used for re-reading the iso_files, see \code{\link{iso_read_dual_inlet}}, \code{\link{iso_read_continuous_flow}} and \code{\link{iso_read_scan}} for details (except \code{read_cache} which is always set to \code{FALSE} to force re-reads).}

\item{stop_if_missing}{whether to stop re-reading if any of the original data files are missing (if FALSE, will warn about the missing files adding a warning to them, but also re-read those that do exist)}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}

\item{reread_files_with_errors}{whether to re-read files that had read in with errors the last time (default TRUE)}

\item{reread_files_with_warnings}{whether to re-read files that had read in with warnings the last time (default TRUE)}
}
\description{
Sometimes it is useful to reload isotope files from their original data files (e.g. after modifying raw data files in vendor software, or after upgrading to a newer version of the isoreader package that provides new functionality). The functions described below are intended to make this very easy. However, re-reading files from disc is only possible if file paths still point to the original raw data files. If they have moved, please use \code{\link{iso_set_file_root}} first to change the root directory of your \code{iso_files}.
}
\details{
To re-read files that have been modified on disc, please use \code{iso_reread_changed_files()}. To re-read files because of an isoreader version upgrade, please use \code{iso_reread_outdated_files()}. To try re-reading files that previously had warnings and/or errors, please use \code{iso_reread_problem_files()}.

\code{iso_reread_all_files} re-reads all files in the collection.

\code{iso_reread_changed_files} re-reads all files that have been modified (e.g. in the vendor software) since they were last read by isoreader.

\code{iso_reread_outdated_files} re-reads all files that were read with an outdated version of isoreader.

\code{iso_reread_problem_files} re-reads all files that have had errors the last time they were read by isoreader (set \code{reread_files_with_warnings = TRUE} to also re-read those that have warnings).

\code{iso_reread_storage} is deprecated.

\code{iso_reread_archive} is deprecated.
}
\examples{
# example for re-reading a saved isofile collection
iso_turn_reader_caching_off()
saved_files_path <- "saved_isofile.scan.rds"

# create saved collection
iso_get_reader_examples_folder() \%>\% 
 iso_read_scan() \%>\%
 iso_save(saved_files_path)
 
# load collection
iso_read_scan(saved_files_path) \%>\%
 # reread outdated files (alternatively "_all_" or "_changed_")
 iso_reread_outdated_files() \%>\%
 # re-save collection to its original location
 iso_save(saved_files_path)

# cleanup
unlink(saved_files_path)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_info_operations.R
\name{iso_rename_file_info}
\alias{iso_rename_file_info}
\title{Rename file info columns}
\usage{
iso_rename_file_info(
  iso_files,
  ...,
  file_specific = FALSE,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{...}{dplyr-style \link[dplyr]{rename} conditions applied based on each file's file_info (see \code{\link{iso_get_file_info}})}

\item{file_specific}{whether to run the select criteria (\code{...}) specifically within each individual file rather than on all files jointly. This is a lot slower but makes it possible to  select different columns in different iso_files depending on what exists in each file and is mostly of use when working with data from multiple instruments.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Rename file info columns (\code{\link{iso_get_file_info}}) within isofile objects. Works just like dplyr's \link[dplyr]{rename}. You can also use \link[dplyr]{rename} directly but it will not provide summary information on the operation. To select specific columns to keep (discarding all others), use \link{iso_select_file_info} instead. Set \code{file_specific = TRUE} to rename different columns in different iso_files depending on what exists in each file. This is very useful when working with data from multiple instruments that may have the same information (e.g. sample name) stored in different columns.
}
\seealso{
Other file_info operations: 
\code{\link{iso_add_file_info.iso_file_list}()},
\code{\link{iso_filter_files}()},
\code{\link{iso_mutate_file_info}()},
\code{\link{iso_parse_file_info}()},
\code{\link{iso_select_file_info}()},
\code{\link{iso_set_file_root}()}
}
\concept{file_info operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{vec_arith.iso_double_with_units}
\alias{vec_arith.iso_double_with_units}
\title{vec_arith for iso_double_with_units}
\usage{
\method{vec_arith}{iso_double_with_units}(op, x, y, ...)
}
\arguments{
\item{op}{An arithmetic operator as a string}

\item{x}{A pair of vectors. For \code{!}, unary \code{+} and unary \code{-}, \code{y} will be
a sentinel object of class \code{MISSING}, as created by \code{MISSING()}.}

\item{y}{A pair of vectors. For \code{!}, unary \code{+} and unary \code{-}, \code{y} will be
a sentinel object of class \code{MISSING}, as created by \code{MISSING()}.}

\item{...}{These dots are for future extensions and must be empty.}
}
\description{
vec_arith for iso_double_with_units
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{iso_strip_units}
\alias{iso_strip_units}
\title{Strip units from variables}
\usage{
iso_strip_units(x)
}
\arguments{
\item{x}{variable to strip units from (vector or data frame)}
}
\description{
This function converts numbers with units back into unitless numbers both for single variables and data frames / tibbles. For single variables, this is equivalent to the \code{as.numeric} function.
}
\seealso{
Other functions for values with units: 
\code{\link{iso_get_units}()},
\code{\link{iso_is_double_with_units}()},
\code{\link{iso_make_units_explicit}()},
\code{\link{iso_make_units_implicit}()},
\code{\link{iso_with_units}()}
}
\concept{functions for values with units}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/settings.R
\name{iso_get_default_reader_parameters}
\alias{iso_get_default_reader_parameters}
\title{Get the current default parameters}
\usage{
iso_get_default_reader_parameters()
}
\description{
Retrieve a table with all default function parameters for this package. 
To set read parameters, see \code{\link{iso_set_default_read_parameters}}. 
To set messaging and caching parameters see \code{\link{iso_info_messages}} and see \code{\link{iso_caching}}.
}
\seealso{
Other settings functions: 
\code{\link{iso_caching}},
\code{\link{iso_info_messages}},
\code{\link{iso_set_default_read_parameters}()}
}
\concept{settings functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{iso_expand_paths}
\alias{iso_expand_paths}
\title{Expand file paths}
\usage{
iso_expand_paths(path, extensions = c(), root = ".")
}
\arguments{
\item{path}{vector of file/folder paths, mixed relative and absolute paths are allowed.}

\item{extensions}{which extensions to look for? (with or without leading .) - this is typically one or more of the extensions listed by \code{\link{iso_get_supported_file_types}}}

\item{root}{root directory for the isofiles. Can be relative to the current working directory (e.g. \code{"data"}) or an absolute path on the file system (e.g. \code{"/Users/..."} or \code{"C:/Data/.."}). The default is the current working directory (\code{"."}). Can be supplied as a vector of same length as the provided paths if the paths have different roots.}
}
\value{
data frame with columns \code{root} (\code{root} as provided) and \code{path} of all the found files.
}
\description{
Helper function to expand the provided paths to find data files in folders and subfolders that match any of the specified extensions. Filepaths will be kept as is, only folders will be expanded. Note that this function is rarely called directly. It is used automatically by \code{\link{iso_read_dual_inlet}} and \code{\link{iso_read_continuous_flow}} to identify files of interest based on the file paths provided.
}
\seealso{
Other file system functions: 
\code{\link{iso_find_absolute_path_roots}()},
\code{\link{iso_root_paths}()},
\code{\link{iso_shorten_relative_paths}()}
}
\concept{file system functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export.R
\name{iso_export_to_feather}
\alias{iso_export_to_feather}
\title{Export to feather}
\usage{
iso_export_to_feather(
  iso_files,
  filepath_prefix,
  include_file_info = everything(),
  include_raw_data = everything(),
  include_standards = !!enexpr(include_method_info),
  include_resistors = !!enquo(include_method_info),
  include_vendor_data_table = everything(),
  include_problems = everything(),
  with_explicit_units = FALSE,
  include_method_info = everything(),
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{filepath_prefix}{what to use as the prefix for the feather file names (e.g. name of the data collection or current date)}

\item{include_file_info}{which file information to include (see \code{\link{iso_get_file_info}}). Use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns.}

\item{include_raw_data}{which columns from the raw data to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes all columns by default. Set to NULL to include no raw data.}

\item{include_standards}{which columns from the standards info to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. By default, everything is included (both standards and ratios). To omit the ratios, change to \code{select = file_id:reference}. Set to NULL to include no standards info.}

\item{include_resistors}{which columns from the resistors info to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes all columns by default. Set to NULL to include no resistors info.}

\item{include_vendor_data_table}{which columns from the vendor data table to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes all columns by default. Set parameter \code{with_explicit_units = TRUE} to make column units explicit (keep in mind that this will require specific \code{include_vendor_data_table} column selections to reflect the column names including the units). Set to NULL to include no vendor data table.}

\item{include_problems}{which columns from problems to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes none of the read problems by default. Set to \code{include_problems = everything()} to include all columns.}

\item{with_explicit_units}{whether to include units in the column headers of the returned data frame instead of the column data types (see \code{\link{iso_double_with_units}}). Note that any \code{select} conditions have to refer to the column names including the full units.}

\item{include_method_info}{deprecated in favor of the more specific include_standards and include_resistors}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\value{
returns the iso_files object invisibly for use in pipelines
}
\description{
This function exports the passed in iso_files to the Python and R shared feather file format. The different kinds of data (raw data, file info, methods info, etc.) are exported to separate feather files that are saved with the provided \code{filepath_prefix} as prefix. All are only exported if the corresponding \code{include_} parameter is set to \code{TRUE} and only for data types for which this type of data is available and was read (see \code{\link{iso_read_dual_inlet}}, \code{\link{iso_read_continuous_flow}} for details on read parameters). Note that in rare instances where vectorized data columns exist in the file information (e.g. measurement_info), they are concatenated with ', ' in feather output. Note that the feather package required for this export is not installed automatically as part of isoreader. Please install it manually if missing using \code{install.packages("feather")}.
}
\seealso{
Other export functions: 
\code{\link{iso_export_to_excel}()},
\code{\link{iso_save}()}
}
\concept{export functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/problems.R
\name{iso_get_problems_summary}
\alias{iso_get_problems_summary}
\title{Retrieve a summary of the problems}
\usage{
iso_get_problems_summary(
  iso_files,
  problem_files_only = TRUE,
  include_file_info = NULL
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{problem_files_only}{whether to list only problem files or all files}

\item{include_file_info}{which file information to include (see \code{\link{iso_get_file_info}}). Use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns.}
}
\value{
data frame with file_id and number of encountered errors and warnings
}
\description{
Returns a data frame listing how many errors and warnings were encountered for each file. For details on each error/warning, see \link[readr]{problems} and the \link{iso_problem_functions}.
}
\seealso{
Other problem functions: 
\code{\link{iso_filter_files_with_problems}()},
\code{\link{iso_get_problems}()},
\code{\link{iso_has_problems}()},
\code{\link{iso_problem_functions}}
}
\concept{problem functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/problems.R
\name{iso_has_problems}
\alias{iso_has_problems}
\title{Check for parsing problems}
\usage{
iso_has_problems(iso_files)
}
\arguments{
\item{iso_files}{collection of iso_file objects}
}
\value{
boolean
}
\description{
Check for parsing problems
}
\seealso{
Other problem functions: 
\code{\link{iso_filter_files_with_problems}()},
\code{\link{iso_get_problems_summary}()},
\code{\link{iso_get_problems}()},
\code{\link{iso_problem_functions}}
}
\concept{problem functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{iso_shorten_relative_paths}
\alias{iso_shorten_relative_paths}
\title{Shorten relative paths}
\usage{
iso_shorten_relative_paths(path, root = ".")
}
\arguments{
\item{path}{vector of file/folder paths, mixed relative and absolute paths are allowed.}

\item{root}{root directory for the isofiles. Can be relative to the current working directory (e.g. \code{"data"}) or an absolute path on the file system (e.g. \code{"/Users/..."} or \code{"C:/Data/.."}). The default is the current working directory (\code{"."}). Can be supplied as a vector of same length as the provided paths if the paths have different roots.}
}
\value{
a data frame with the root directories and paths relative to the root - order of input paths is preserved
}
\description{
Convenience function to shorten relative paths based on overlap with the provided root(s). Also simplifies current directory repeats (e.g. "././." becomes ".") for better legibility. Does not check whether the original or resulting paths point to valid files or folders. Relative paths that do not start with the supplied \code{root} default back to the current working directory (\code{.}). Absolute paths are allowed but are returned as is without attempts at shortening. See \code{iso_find_absolute_path_roots} for rooting absolute paths.
}
\examples{
iso_shorten_relative_paths(file.path("A", "B", "C"), "A") # root = "A", path = B/C
iso_shorten_relative_paths(file.path("A", "B", "C"), file.path("A", "B")) # root = "A/B", path = "C"
iso_shorten_relative_paths(file.path("A", "C", "D"), file.path("A", "B")) # root = "A", path = "C/D"
iso_shorten_relative_paths(file.path("A", "B", "C"), "B") # root = ".", path stays "A/B/C"
}
\seealso{
Other file system functions: 
\code{\link{iso_expand_paths}()},
\code{\link{iso_find_absolute_path_roots}()},
\code{\link{iso_root_paths}()}
}
\concept{file system functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_info_operations.R
\name{iso_add_file_info.iso_file_list}
\alias{iso_add_file_info.iso_file_list}
\alias{iso_add_file_info.data.frame}
\alias{iso_add_file_info}
\title{Add additional file information}
\usage{
\method{iso_add_file_info}{iso_file_list}(iso_files, new_file_info, ..., quiet = default(quiet))

\method{iso_add_file_info}{data.frame}(df, new_file_info, ..., quiet = default(quiet))

iso_add_file_info(...)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{new_file_info}{data frame with new file information to add to the isofiles}

\item{...}{each parameter specifies a set of \code{join_by} column(s) to add the \code{new_file_info} to the existing file information. The provided parameters are applied sequentially. At least one must be specified.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}

\item{df}{a data frame of iso files data retrieved by any of the data retrieval functions (e.g. \code{\link{iso_get_file_info}}, \code{\link{iso_get_raw_data}, etc.}}
}
\value{
the original iso files or data frame with the new file info added in.
}
\description{
This function makes it easy to add additional file info (\code{\link{iso_get_file_info}}) to isofile objects and data frames by a single \code{\link[dplyr:mutate-joins]{left_join}} or multiple sequential \code{\link[dplyr:mutate-joins]{left_join}} operations. The function provides a detailed summary of the information that was added unless \code{quiet = TRUE}. Note that one-to-many joins are not permitted (and will fail with an informative error) since this would lead to likely unintended data duplication in the isofiles. However, one-to-one and many-to-one joins are fully supported and should cover all needed use cases for this function. Also note that for each join, only the \code{new_file_info} rows that have defined non-NA, non-empty ("") values in all \code{join_by} columns will be considered for the join and that only \code{new_file_info} columns that do NOT already exist in ANY file information will be added. For changing the values of existing file information, please use \code{\link{iso_mutate_file_info}} instead.
}
\details{
Single \code{\link[dplyr:mutate-joins]{left_join}}: this is the most common use of this function and basically a simple left join operation (with some additional safety checks). Specify a single \code{join_by} in the \code{...}, such as e.g. \code{c("file_id")} to add additional file information joining by the \code{file_id} column.

Multiple sequential \code{\link[dplyr:mutate-joins]{left_join}}: this use case is for applying a set of increasingly more specific \code{join_by} rules. For example, \code{... = c("Identifier 1", "Identifier 2"), c("file_id")} would serve to first add one set of new file information for all isofiles based on their \code{Identifier 1} and \code{Identifier 2} columns and then overwrite the new information with more specific details for a subset of isofiles based on their \code{file_id} column, all based on a single overview \code{new_file_info} data frame. Basically, each set of \code{join_by} conditions specified in \code{...} must describe a valid \code{\link[dplyr:mutate-joins]{left_join}} \code{join_by} parameter to merge the \code{new_file_info} with the existing file info. Each set of \code{new_file_info} data can overwrite the previous \code{join_by} matches such that the last set of \code{join_by} column(s) provided in \code{...} will overwrite all previous matches for which it applies, even if they have already been a match for a previous column.
}
\seealso{
Other file_info operations: 
\code{\link{iso_filter_files}()},
\code{\link{iso_mutate_file_info}()},
\code{\link{iso_parse_file_info}()},
\code{\link{iso_rename_file_info}()},
\code{\link{iso_select_file_info}()},
\code{\link{iso_set_file_root}()}
}
\concept{file_info operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/settings.R
\name{set_temp}
\alias{set_temp}
\title{Set temporary option}
\usage{
set_temp(name, value)
}
\arguments{
\item{name}{name of the temporary option}

\item{value}{value of the temporary option}
}
\description{
Set a temporary option for parallel processing in isoprocessor.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{iso_get_units}
\alias{iso_get_units}
\title{Retrieve number units}
\usage{
iso_get_units(x)
}
\arguments{
\item{x}{variable to get the units for (vector or data frame)}
}
\description{
This function returns the units of a numerical value generated by \code{\link{iso_double_with_units}}. It returns \code{NA}) for unitless variables. Returns a column-named vector of units if \code{x} is a data frame / tibble. Returns the direct units of \code{x} in all other cases.
}
\seealso{
Other functions for values with units: 
\code{\link{iso_is_double_with_units}()},
\code{\link{iso_make_units_explicit}()},
\code{\link{iso_make_units_implicit}()},
\code{\link{iso_strip_units}()},
\code{\link{iso_with_units}()}
}
\concept{functions for values with units}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/problems.R
\name{iso_problem_functions}
\alias{iso_problem_functions}
\title{Problem Functions Overview}
\description{
The following functions to check for and deal with problems are available.
}
\details{
\itemize{
\item \code{iso_get_problems} is a re-export of \code{\link[readr]{problems}}

\item \code{\link{iso_get_problems_summary}}

\item \code{\link{iso_has_problems}}

\item \code{\link[readr:problems]{stop_for_problems}}

\item \code{\link{iso_filter_files_with_problems}}
}
}
\seealso{
Other problem functions: 
\code{\link{iso_filter_files_with_problems}()},
\code{\link{iso_get_problems_summary}()},
\code{\link{iso_get_problems}()},
\code{\link{iso_has_problems}()}
}
\concept{problem functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/settings.R
\name{iso_info_messages}
\alias{iso_info_messages}
\alias{iso_turn_info_messages_on}
\alias{iso_turn_info_messages_off}
\alias{iso_turn_datetime_warnings_on}
\alias{iso_turn_datetime_warnings_off}
\title{Control information messages}
\usage{
iso_turn_info_messages_on(data = NULL)

iso_turn_info_messages_off(data = NULL)

iso_turn_datetime_warnings_on(data = NULL)

iso_turn_datetime_warnings_off(data = NULL)
}
\arguments{
\item{data}{a data frame - returned invisibly as is if provided (e.g. in the middle of a pipeline)}
}
\description{
These functions control the global settings for information messages.
}
\details{
\code{iso_turn_info_messages_on()} and \code{iso_turn_info_messages_off()} turn information messages on/off in all subsequent function calls by changing the global settings for the \code{quiet} parameter of most isoreader functions. These functions can be called stand alone or within a pipeline to turn messages on/off at a certain point during the pipeline.

\code{iso_turn_datetime_warnings_on()} and \code{iso_turn_datetime_warnings_off()} turn datetime warnings that occur on some platforms (mostly linux distributions) on/off for all subsequent isoreader functions. These warnings inform the user that file creation dates are not available from the operating system.
}
\seealso{
Other settings functions: 
\code{\link{iso_caching}},
\code{\link{iso_get_default_reader_parameters}()},
\code{\link{iso_set_default_read_parameters}()}
}
\concept{settings functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/problems.R
\name{iso_get_problems}
\alias{iso_get_problems}
\title{Retrieve parsing problems}
\usage{
iso_get_problems(iso_files, select = everything())
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{select}{which data columns to select - use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax. By default, all columns are selected.}
}
\description{
This function retrieves parsing problems encountered during the reading of a set of iso files.
}
\seealso{
Other problem functions: 
\code{\link{iso_filter_files_with_problems}()},
\code{\link{iso_get_problems_summary}()},
\code{\link{iso_has_problems}()},
\code{\link{iso_problem_functions}}
}
\concept{problem functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{iso_is_double_with_units}
\alias{iso_is_double_with_units}
\title{Check if a value has units}
\usage{
iso_is_double_with_units(x)
}
\arguments{
\item{x}{vector to check for whether it is a double with units}
}
\description{
Check if a variable is a double with units. That is if it has been generated by \code{\link{iso_double_with_units}}.
}
\seealso{
Other functions for values with units: 
\code{\link{iso_get_units}()},
\code{\link{iso_make_units_explicit}()},
\code{\link{iso_make_units_implicit}()},
\code{\link{iso_strip_units}()},
\code{\link{iso_with_units}()}
}
\concept{functions for values with units}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isosave.R
\name{iso_save}
\alias{iso_save}
\title{Save data to R Data Storage (.rds)}
\usage{
iso_save(iso_files, filepath, quiet = default(quiet))
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{filepath}{the path (folder and filename) to the export file. The correct file extension is automatically added if not already in the filename, i.e. filename can be provided with or without extension.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\value{
returns the iso_files object invisibly for use in pipelines
}
\description{
This function saves the passed in iso_files to an R Data Storage (.rds) file, which is an efficient compressed data storage format. Data exported this way can be easily read back into isoreader using the standard \code{\link{iso_read_continuous_flow}} and \code{\link{iso_read_dual_inlet}} functions.
}
\seealso{
Other export functions: 
\code{\link{iso_export_to_excel}()},
\code{\link{iso_export_to_feather}()}
}
\concept{export functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/settings.R
\name{iso_debug_mode}
\alias{iso_debug_mode}
\alias{iso_turn_debug_on}
\alias{iso_turn_debug_off}
\alias{set_read_file_event_expr}
\alias{set_finish_file_event_expr}
\title{Debugging functions}
\usage{
iso_turn_debug_on(data = NULL, catch_errors = TRUE, cache = FALSE)

iso_turn_debug_off(data = NULL)

set_read_file_event_expr(event_expr = NULL)

set_finish_file_event_expr(event_expr = NULL)
}
\arguments{
\item{data}{a data frame - returned invisibly as is if provided (e.g. in the middle of a pipeline)}

\item{catch_errors}{whether to still catch errors in debug mode or whether to throw them}

\item{cache}{whether to cache or read anything from cache}

\item{event_expr}{an expression to evaluate in the context of reading individual iso files (evaluated in the local environment at the beginning of a file read)}
}
\description{
For troubleshooting. Not exported.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_standards}
\alias{iso_get_standards}
\title{Aggregate standards from methods info}
\usage{
iso_get_standards(
  iso_files,
  select = everything(),
  include_file_info = NULL,
  with_ratios = NULL,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{select}{which data columns to select - use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax. By default, everything is included (both standards and ratios). To omit the ratios, change to \code{select = file_id:reference}.}

\item{include_file_info}{which file information to include (see \code{\link{iso_get_file_info}}). Use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns.}

\item{with_ratios}{deprecated, please use the \code{select} parameter to explicitly include or exclude ratio columns}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Aggregates the isotopic standard information recovered from the provided iso_files. Can aggregate just the standards' delta values or combine the delta values with the recovered ratios (if any). Use parameter \code{select} to exclude/include the ratios. All standards info is only available if the iso_files were read with parameter \code{read_method_info=TRUE}.
}
\seealso{
Other data retrieval functions: 
\code{\link{iso_get_all_data}()},
\code{\link{iso_get_bgrd_data}()},
\code{\link{iso_get_file_info}()},
\code{\link{iso_get_raw_data}()},
\code{\link{iso_get_resistors}()},
\code{\link{iso_get_vendor_data_table}()}
}
\concept{data retrieval functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{iso_register_dual_inlet_file_reader}
\alias{iso_register_dual_inlet_file_reader}
\alias{iso_register_continuous_flow_file_reader}
\alias{iso_register_scan_file_reader}
\title{Register file readers}
\usage{
iso_register_dual_inlet_file_reader(
  extension,
  func,
  description = NA_character_,
  software = NA_character_,
  cacheable = TRUE,
  post_read_check = TRUE,
  overwrite = FALSE,
  env = find_func(func)
)

iso_register_continuous_flow_file_reader(
  extension,
  func,
  description = NA_character_,
  software = NA_character_,
  cacheable = TRUE,
  post_read_check = TRUE,
  overwrite = FALSE,
  env = find_func(func)
)

iso_register_scan_file_reader(
  extension,
  func,
  description = NA_character_,
  software = NA_character_,
  cacheable = TRUE,
  post_read_check = TRUE,
  overwrite = FALSE,
  env = find_func(func)
)
}
\arguments{
\item{extension}{the file extension (e.g. \code{.dxf}) of the data file. Must be unique otherwise different files can not automatically be matched with the appropriate file reader based on their extension.}

\item{func}{the name of the function that should be used a filter reader. All file reader functions must accept a data structure argument as the first argument and return the same data structure with added data.}

\item{description}{what is this file type about?}

\item{software}{what is the software program that creates this file type?}

\item{cacheable}{whether this file type is cacheable. If \code{TRUE} (the default), user requests to cache the file will be honored. If \code{FALSE}, this file type will never be cached no matter what the user requests.}

\item{post_read_check}{whether isoreader should conduct a data integrity check after reading the file. Should always be \code{TRUE} unless there is independent data integrity checking already taking place inside the reader.}

\item{overwrite}{whether to overwrite an existing file reader for the same extension}

\item{env}{the environment where to find the function, by default this will be determined automatically and will throw an error if there is any ambiguity (e.g. the same function name in multiple packages) in which case it should be set manually}
}
\description{
Register file extensions and reader functions for different data files. Isoreader automatically registers all built-in file readers so this function is usually only needed when registering additional readers provided for testing purposes from outside of the isoreader package. Note that file extensions are case-insensitive, i.e. a reader for \code{.ext} will also recognize \code{.Ext} and \code{.EXT}
}
\details{
\code{iso_register_dual_inlet_file_reader}: use this function to register file readers for dual inlet files.

\code{iso_register_continuous_flow_file_reader}: use this function to register file readers for continuous flow files.

\code{iso_register_scan_file_reader}: use this function to register file readers for scan files.
}
\seealso{
Other file_types: 
\code{\link{iso_get_supported_file_types}()}

Other file_types: 
\code{\link{iso_get_supported_file_types}()}

Other file_types: 
\code{\link{iso_get_supported_file_types}()}
}
\concept{file_types}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{iso_get_reader_example}
\alias{iso_get_reader_example}
\alias{iso_get_reader_examples}
\alias{iso_get_reader_examples_folder}
\title{Example files}
\usage{
iso_get_reader_example(filename)

iso_get_reader_examples()

iso_get_reader_examples_folder()
}
\arguments{
\item{filename}{the name of the example file for which to retrieve the system path}
}
\description{
The isoreader package comes with a few example files to make it easy to illustrate the functionality.
}
\details{
\code{iso_get_reader_example}: retrieve the path to an isoreader example file

\code{iso_get_reader_examples}: list of all available isoreader example files

\code{iso_get_reader_examples_folder}: path to the location of the reader examples
}
\examples{
iso_get_reader_examples()
iso_get_reader_examples_folder()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_vendor_data_table}
\alias{iso_get_vendor_data_table}
\title{Aggregate vendor computed table data}
\usage{
iso_get_vendor_data_table(
  iso_files,
  with_units = FALSE,
  select = everything(),
  include_file_info = NULL,
  with_explicit_units = with_units,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{with_units}{this parameter has been DEPRECATED with the introduction of unit-data types (see \code{\link{iso_double_with_units}}) and will be removed in future versions of isoreader. Please use \code{with_explicit_units} instead if you really want columns to have units explicitly in the column name. Alternatively, consider working with the new implicit unit system and convert vendor data tables as needed with \code{\link{iso_make_units_explicit}} and \code{\link{iso_make_units_implicit}}.}

\item{select}{which data columns to select - use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax. By default, all columns are selected.}

\item{include_file_info}{which file information to include (see \code{\link{iso_get_file_info}}). Use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns.}

\item{with_explicit_units}{whether to include units in the column headers of the returned data frame instead of the column data types (see \code{\link{iso_double_with_units}}). Note that any \code{select} conditions have to refer to the column names including the full units.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Aggregate data from the vendor-computed data table. This information is only available if the iso_files were read with parameter \code{read_vendor_data_table=TRUE}.
}
\seealso{
Other data retrieval functions: 
\code{\link{iso_get_all_data}()},
\code{\link{iso_get_bgrd_data}()},
\code{\link{iso_get_file_info}()},
\code{\link{iso_get_raw_data}()},
\code{\link{iso_get_resistors}()},
\code{\link{iso_get_standards}()}
}
\concept{data retrieval functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/settings.R
\name{iso_caching}
\alias{iso_caching}
\alias{iso_turn_reader_caching_on}
\alias{iso_turn_reader_caching_off}
\title{Turn caching on/off}
\usage{
iso_turn_reader_caching_on(data = NULL)

iso_turn_reader_caching_off(data = NULL)
}
\arguments{
\item{data}{a data frame - returned invisibly as is if provided (e.g. in the middle of a pipeline)}
}
\description{
These functions turn caching of data files (and reading from cache) on/off in all subsequent isoread calls by changing the global settings for the \code{cache} parameter. Can be called stand alone or within a pipeline.
}
\seealso{
Other settings functions: 
\code{\link{iso_get_default_reader_parameters}()},
\code{\link{iso_info_messages}},
\code{\link{iso_set_default_read_parameters}()}
}
\concept{settings functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unit_conversion.R
\name{iso_convert_time}
\alias{iso_convert_time}
\title{moved to isoprocessor}
\usage{
iso_convert_time(...)
}
\arguments{
\item{...}{deprecated}
}
\description{
moved to isoprocessor
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isodata_structures.R
\name{iso_is_file}
\alias{iso_is_file}
\alias{iso_is_file_list}
\alias{iso_is_object}
\alias{iso_is_dual_inlet}
\alias{iso_is_continuous_flow}
\alias{iso_is_scan}
\alias{iso_as_file_list}
\title{Isoreader data structure functions}
\usage{
iso_is_file(x)

iso_is_file_list(x)

iso_is_object(x)

iso_is_dual_inlet(x)

iso_is_continuous_flow(x)

iso_is_scan(x)

iso_as_file_list(..., discard_duplicates = TRUE)
}
\arguments{
\item{x}{an object to test whether it has the specific class}

\item{...}{iso_file and iso_file_list objects to concatenate}

\item{discard_duplicates}{whether to automatically discard files with duplicate file IDs (i.e. duplicate file names). If \code{TRUE} (the default), only the first files are kept and any files with the same file ID are discarded. If \code{FALSE}, all duplicate files are kept but their file IDs are appended with suffix \code{#1}, \code{#2}, etc.}
}
\description{
\code{iso_is_file} tests if the object is an iso_file

\code{iso_is_file_list} tests if the object is an iso_file list (collection of iso_files)

\code{iso_is_object} test if the object is an iso-object (iso_file or iso_file list)

\code{iso_is_dual_inlet} tests if an iso_file or iso_file list consists exclusively of dual inlet file objects

\code{iso_is_continuous_flow} tests if an iso_file or iso_file list consists exclusively of continuous flow file objects

\code{iso_is_scan} tests if an iso_file or iso_file list consists exclusively of scan file objects

\code{iso_as_file_list} concatenates iso_file and iso_file list object(s) into one combined iso_file list (equivalent to calling \code{c(...)}), flattens all passed lists into one list structure, all individual objects and objects within iso_file lists have to be the same type of iso_file, issues warnings if there are duplicate file ids and summarizes all problems in the iso_file list. If duplicates are allowed (\code{discard_duplicates = FALSE}), their file IDs will append a #1, #2, #3, etc. to preserve unique file IDs (important for many data aggregation operations).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export.R
\name{iso_export_to_excel}
\alias{iso_export_to_excel}
\title{Export data to Excel}
\usage{
iso_export_to_excel(
  iso_files,
  filepath,
  include_file_info = everything(),
  include_raw_data = everything(),
  include_standards = !!enexpr(include_method_info),
  include_resistors = !!enquo(include_method_info),
  include_vendor_data_table = everything(),
  include_problems = everything(),
  with_explicit_units = FALSE,
  include_method_info = everything(),
  with_ratios = NULL,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{filepath}{the path (folder and filename) to the export file. The correct file extension is automatically added if not already in the filename, i.e. filename can be provided with or without extension.}

\item{include_file_info}{which file information to include (see \code{\link{iso_get_file_info}}). Use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns.}

\item{include_raw_data}{which columns from the raw data to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes all columns by default. Set to NULL to include no raw data.}

\item{include_standards}{which columns from the standards info to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. By default, everything is included (both standards and ratios). To omit the ratios, change to \code{select = file_id:reference}. Set to NULL to include no standards info.}

\item{include_resistors}{which columns from the resistors info to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes all columns by default. Set to NULL to include no resistors info.}

\item{include_vendor_data_table}{which columns from the vendor data table to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes all columns by default. Set parameter \code{with_explicit_units = TRUE} to make column units explicit (keep in mind that this will require specific \code{include_vendor_data_table} column selections to reflect the column names including the units). Set to NULL to include no vendor data table.}

\item{include_problems}{which columns from problems to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes none of the read problems by default. Set to \code{include_problems = everything()} to include all columns.}

\item{with_explicit_units}{whether to include units in the column headers of the returned data frame instead of the column data types (see \code{\link{iso_double_with_units}}). Note that any \code{select} conditions have to refer to the column names including the full units.}

\item{include_method_info}{deprecated in favor of the more specific include_standards and include_resistors}

\item{with_ratios}{deprecated, please use the \code{select} parameter to explicitly include or exclude ratio columns}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\value{
returns the iso_files object invisibly for use in pipelines
}
\description{
This function exports the passed in iso_files to Excel. The different kinds of data (raw data, file info, methods info, etc.) are exported to separate tabs within the excel file. Use the various \code{include_...} parameters to specify what information to include. Note that in rare instances where vectorized data columns exist in the file information (e.g. measurement_info), they are concatenated with ', ' in the excel export. Note that the openxlsx package required for this export is not installed automatically as part of isoreader. Please install it manually if missing using \code{install.packages("openxlsx")}.
}
\seealso{
Other export functions: 
\code{\link{iso_export_to_feather}()},
\code{\link{iso_save}()}
}
\concept{export functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleanup.R, R/package.R, R/problems.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{parse_number}
\alias{parse_logical}
\alias{parse_integer}
\alias{parse_double}
\alias{parse_datetime}
\alias{!!}
\alias{!!!}
\alias{\%>\%}
\alias{everything}
\alias{starts_with}
\alias{ends_with}
\alias{matches}
\alias{filter}
\alias{tibble}
\alias{problems}
\alias{stop_for_problems}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{dplyr}{\code{\link[dplyr]{filter}}}

  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}

  \item{readr}{\code{\link[readr]{parse_datetime}}, \code{\link[readr:parse_atomic]{parse_double}}, \code{\link[readr:parse_atomic]{parse_integer}}, \code{\link[readr:parse_atomic]{parse_logical}}, \code{\link[readr]{parse_number}}, \code{\link[readr]{problems}}, \code{\link[readr:problems]{stop_for_problems}}}

  \item{rlang}{\code{\link[rlang:nse-force]{!!}}, \code{\link[rlang:nse-force]{!!!}}}

  \item{tibble}{\code{\link[tibble]{tibble}}}

  \item{tidyselect}{\code{\link[tidyselect:starts_with]{ends_with}}, \code{\link[tidyselect]{everything}}, \code{\link[tidyselect:starts_with]{matches}}, \code{\link[tidyselect]{starts_with}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{iso_get_supported_file_types}
\alias{iso_get_supported_file_types}
\title{Supported file types}
\usage{
iso_get_supported_file_types()
}
\description{
Get an overview of all the file types currently supported by the isoreader package. To register additional file readers, use the \code{\link{iso_register_dual_inlet_file_reader}} and \code{\link{iso_register_continuous_flow_file_reader}} functions.
}
\seealso{
Other file_types: 
\code{\link{iso_register_dual_inlet_file_reader}()}
}
\concept{file_types}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_all_data}
\alias{iso_get_all_data}
\title{Aggregate all isofiles data}
\usage{
iso_get_all_data(
  iso_files,
  include_file_info = everything(),
  include_raw_data = everything(),
  include_standards = everything(),
  include_resistors = everything(),
  include_vendor_data_table = everything(),
  include_problems = NULL,
  gather = FALSE,
  with_explicit_units = with_units,
  with_units = FALSE,
  with_ratios = NULL,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{include_file_info}{which file information to include (see \code{\link{iso_get_file_info}}). Use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns.}

\item{include_raw_data}{which columns from the raw data to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes all columns by default. Set to NULL to include no raw data.}

\item{include_standards}{which columns from the standards info to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. By default, everything is included (both standards and ratios). To omit the ratios, change to \code{select = file_id:reference}. Set to NULL to include no standards info.}

\item{include_resistors}{which columns from the resistors info to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes all columns by default. Set to NULL to include no resistors info.}

\item{include_vendor_data_table}{which columns from the vendor data table to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes all columns by default. Set parameter \code{with_explicit_units = TRUE} to make column units explicit (keep in mind that this will require specific \code{include_vendor_data_table} column selections to reflect the column names including the units). Set to NULL to include no vendor data table.}

\item{include_problems}{which columns from problems to include. Use \code{c(...)} to select multiple columns, supports all \link[dplyr]{select} syntax including renaming columns. Includes none of the read problems by default. Set to \code{include_problems = everything()} to include all columns.}

\item{gather}{whether to gather raw data into long format (e.g. for ease of use in plotting). Not that the \code{select} parameter applies to the data columns BEFORE gathering.}

\item{with_explicit_units}{whether to include units in the column headers of the returned data frame instead of the column data types (see \code{\link{iso_double_with_units}}). Note that any \code{select} conditions have to refer to the column names including the full units.}

\item{with_units}{this parameter has been DEPRECATED with the introduction of unit-data types (see \code{\link{iso_double_with_units}}) and will be removed in future versions of isoreader. Please use \code{with_explicit_units} instead if you really want columns to have units explicitly in the column name. Alternatively, consider working with the new implicit unit system and convert vendor data tables as needed with \code{\link{iso_make_units_explicit}} and \code{\link{iso_make_units_implicit}}.}

\item{with_ratios}{deprecated, please use the \code{select} parameter to explicitly include or exclude ratio columns}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\value{
data_frame with file_ids, file_types and nested data frames for each data type (file_info, raw_data, vendor_data_table, etc.)
}
\description{
This function aggregates all isofiles data and returns it in a large data frame with nested columns for each type of information (file_info, raw_data, etc.). For targeted retrieval of specific data \code{\link{iso_get_raw_data}}, \code{\link{iso_get_file_info}}, \code{\link{iso_get_vendor_data_table}}, etc. are much faster and easier to work with. This function is primarily useful for downstream processing pipelines that want to carry all information along. To \code{\link[tidyr:nest]{unnest}} any of the specific data types (e.g. \code{raw_data}), make sure to filter first for the files that have this data type available (e.g. \code{filter(has_raw_data)}). Exclude specific types of information by setting its \code{include...} parameter to \code{NULL} (Note: for historical reasons, setting it to \code{FALSE} will also include the information).
}
\seealso{
Other data retrieval functions: 
\code{\link{iso_get_bgrd_data}()},
\code{\link{iso_get_file_info}()},
\code{\link{iso_get_raw_data}()},
\code{\link{iso_get_resistors}()},
\code{\link{iso_get_standards}()},
\code{\link{iso_get_vendor_data_table}()}
}
\concept{data retrieval functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{iso_with_units}
\alias{iso_with_units}
\alias{iso_double_with_units}
\title{Generate values with units}
\usage{
iso_with_units(x, units = "undefined units")

iso_double_with_units(x = double(), units = "undefined units")
}
\arguments{
\item{x}{the values (single value or vector)}

\item{units}{the units for the value, by default "undefined units" but this parameter should always be supplied when working with real data that has units}
}
\description{
These functions generate values with units that work well within data frames and tibbles and implement safety checks on operations that combine values with different units. To retrieve the value without units, use \code{\link{iso_strip_units}} (works for single variables and data frames/tibbles). To retrieve the unit use \code{\link{iso_get_units}}. Note that to correctly combine data frames / tibbles that have values with units in them, use \link[vctrs:vec_bind]{vec_rbind} instead of \link{rbind} or \link[dplyr:bind]{bind_rows}. \link[vctrs:vec_bind]{vec_rbind} will combine columns that have values with units if they have the same unit and otherwise convert back to plain values without units with a warning. The other functions will either fail or reduce the unit values to plain values with a cryptic warning message about not preserving attributes.
}
\details{
\code{iso_with_units} is the primary function to generate values with units. At present, only numeric values are supported so this function is just a shorter alias for the number-specific \code{iso_double_with_units}. It is not clear yet whether any non-numeric values with units make sense to be supported at a later point or whether integer and decimal numbers should be treated differently when they have units.
}
\seealso{
Other functions for values with units: 
\code{\link{iso_get_units}()},
\code{\link{iso_is_double_with_units}()},
\code{\link{iso_make_units_explicit}()},
\code{\link{iso_make_units_implicit}()},
\code{\link{iso_strip_units}()}

Other functions for values with units: 
\code{\link{iso_get_units}()},
\code{\link{iso_is_double_with_units}()},
\code{\link{iso_make_units_explicit}()},
\code{\link{iso_make_units_implicit}()},
\code{\link{iso_strip_units}()}
}
\concept{functions for values with units}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_resistors_info}
\alias{iso_get_resistors_info}
\title{DEPRECATED}
\usage{
iso_get_resistors_info(...)
}
\arguments{
\item{...}{forwarded to \link{iso_get_resistors}}
}
\description{
Please use \link{iso_get_resistors} instead.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{iso_read_scan}
\alias{iso_read_scan}
\title{Load scan data}
\usage{
iso_read_scan(
  ...,
  root = ".",
  read_raw_data = default(read_raw_data),
  read_file_info = default(read_file_info),
  read_method_info = default(read_method_info),
  discard_duplicates = TRUE,
  parallel = FALSE,
  parallel_plan = future::multisession,
  parallel_cores = future::availableCores(),
  cache = default(cache),
  read_cache = default(cache),
  reread_outdated_cache = FALSE,
  quiet = default(quiet),
  cache_files_with_errors = TRUE
)
}
\arguments{
\item{...}{one or multiple file/folder paths. All files must have a supported file extension. All folders are expanded and searched for files with supported file extensions (which are then included in the read).}

\item{root}{root directory for the isofiles. Can be relative to the current working directory (e.g. \code{"data"}) or an absolute path on the file system (e.g. \code{"/Users/..."} or \code{"C:/Data/.."}). The default is the current working directory (\code{"."}). Can be supplied as a vector of same length as the provided paths if the paths have different roots.}

\item{read_raw_data}{whether to read the raw mass/ion data from the file}

\item{read_file_info}{whether to read auxiliary file information (file id, sequence information, etc.)}

\item{read_method_info}{whether to read methods information (standards, processing info)}

\item{discard_duplicates}{whether to automatically discard files with duplicate file IDs (i.e. duplicate file names). If \code{TRUE} (the default), only the first files are kept and any files with the same file ID are discarded. If \code{FALSE}, all duplicate files are kept but their file IDs are appended with suffix \code{#1}, \code{#2}, etc.}

\item{parallel}{whether to process in parallel based on the number of available CPU cores. This may yield performance increases for files that are slow to parse such as continuous flow isodat files but usually provides little benefit for efficient data formats such as reading from R Data Archives.}

\item{parallel_plan}{which parallel processing strategy to use, see \link[future]{plan}, typically \code{future::multisession} for compatibility with RStudio interactive mode. If supported by the operating system and running in detached mode (not interactively in RStudio) can also use \code{future::multicore}.}

\item{parallel_cores}{how many processor cores to use for parallel processing. By default the maximum available number of cores (\link[future]{availableCores}), which will allow maximal processing speed but may slow other programs running on your machine. Choose a smaller number if you want some processing resources to remain available for other processes. Will issue a warning if too many cores are requested and reset to the maximum available.}

\item{cache}{whether to cache iso_files. Note that R Data Storage files (.rds, see \link{iso_save}) are never cached since they are already essentially in cached form.}

\item{read_cache}{whether to reload from cache if a cached version exists. Note that it will only read from cache if the raw data file has not been modified since. Files that have been modified on disc (e.g. edited in the vendor software) will always be read anew. To automatically reread cached files that were cached by an outdated version of the isoreader package, set the \code{reread_outdated_cache} flag.}

\item{reread_outdated_cache}{whether to re-read outdated cache files whenever they are encountered.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}

\item{cache_files_with_errors}{deprecated. Please use \link{iso_reread_problem_files} instead to selectively re-read all files in a collection of iso files that had been previously read with errors or warnings.}
}
\description{
Load scan data
}
\seealso{
Other isoread functions for different types of IRMS data: 
\code{\link{iso_read_continuous_flow}()},
\code{\link{iso_read_dual_inlet}()}
}
\concept{isoread functions for different types of IRMS data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_ratios.R
\name{iso_calculate_ratios}
\alias{iso_calculate_ratios}
\title{moved to isoprocessor}
\usage{
iso_calculate_ratios(...)
}
\arguments{
\item{...}{deprecated}
}
\description{
moved to isoprocessor
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{iso_plot_continuous_flow_data}
\alias{iso_plot_continuous_flow_data}
\title{moved to isoprocessor}
\usage{
iso_plot_continuous_flow_data(...)
}
\arguments{
\item{...}{deprecated}
}
\description{
moved to isoprocessor
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{reread_iso_files}
\alias{reread_iso_files}
\title{Re-read iso_files}
\usage{
reread_iso_files(
  iso_files,
  ...,
  stop_if_missing = FALSE,
  reread_only_changed_files = FALSE,
  reread_only_outdated_files = FALSE,
  reread_files_without_problems = TRUE,
  reread_files_with_errors = TRUE,
  reread_files_with_warnings = TRUE,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_files}

\item{...}{additional read parameters that should be used for re-reading the iso_files, see \code{\link{iso_read_dual_inlet}}, \code{\link{iso_read_continuous_flow}} and \code{\link{iso_read_scan}} for details (except \code{read_cache} which is always set to \code{FALSE} to force re-reads).}

\item{stop_if_missing}{whether to stop re-reading if any of the original data files are missing (if FALSE, will warn about the missing files adding a warning to them, but also re-read those that do exist)}

\item{reread_only_changed_files}{whether to re-read only files that have since be changed on disc (i.e. have no valid cache file), default FALSE i.e. re-read ALL files}

\item{reread_only_outdated_files}{whether to re-read only files that were read by an outdated version of isoreader (default FALSE, i.e. re-read ALL files)}

\item{reread_files_without_problems}{whether to re-read files that had read in without problems the last time (default TRUE)}

\item{reread_files_with_errors}{whether to re-read files that had read in with errors the last time (default TRUE)}

\item{reread_files_with_warnings}{whether to re-read files that had read in with warnings the last time (default TRUE)}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Actual multi-purpose file-reread function (not exported) that powers \link{iso_reread_files}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_info_operations.R
\name{iso_filter_files}
\alias{iso_filter_files}
\title{Filter iso_files}
\usage{
iso_filter_files(iso_files, ..., quiet = default(quiet))
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{...}{dplyr-style \link[dplyr]{filter} conditions applied based on each file's file_info (see \code{\link{iso_get_file_info}})}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Filter for specific isofiles using file info columns (\code{\link{iso_get_file_info}}). Works just like dplyr's \link[dplyr]{filter} except that it provides the user with some information on what has been filtered. Returns \code{NULL} if none of the isofiles' file info matches the filter criteria. You can also use \link[dplyr]{filter} directly to filter collections of \code{iso_file} objects.
}
\seealso{
Other file_info operations: 
\code{\link{iso_add_file_info.iso_file_list}()},
\code{\link{iso_mutate_file_info}()},
\code{\link{iso_parse_file_info}()},
\code{\link{iso_rename_file_info}()},
\code{\link{iso_select_file_info}()},
\code{\link{iso_set_file_root}()}
}
\concept{file_info operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_resistors}
\alias{iso_get_resistors}
\title{Aggregate resistors from methods info}
\usage{
iso_get_resistors(
  iso_files,
  select = everything(),
  include_file_info = NULL,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{select}{which data columns to select - use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax. By default, all columns are selected.}

\item{include_file_info}{which file information to include (see \code{\link{iso_get_file_info}}). Use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Aggregates the resistor information recovered from the provided iso_files. This information is only available if the iso_files were read with parameter \code{read_method_info=TRUE} and only linked to specific masses if the iso_files were additionally read with parameter \code{read_raw_data=TRUE}.
}
\seealso{
Other data retrieval functions: 
\code{\link{iso_get_all_data}()},
\code{\link{iso_get_bgrd_data}()},
\code{\link{iso_get_file_info}()},
\code{\link{iso_get_raw_data}()},
\code{\link{iso_get_standards}()},
\code{\link{iso_get_vendor_data_table}()}
}
\concept{data retrieval functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_binary_files.R
\name{print.binary_structure_map}
\alias{print.binary_structure_map}
\title{Print binary structure map}
\usage{
\method{print}{binary_structure_map}(
  x,
  ...,
  data_as_raw = FALSE,
  line_break_blocks = c("cblock", "stx", "etx"),
  pos_info = TRUE
)
}
\arguments{
\item{x}{object to show.}

\item{...}{additional parameters passed to print.default}

\item{data_as_raw}{whether to show data as raw}

\item{line_break_blocks}{at which blocks to introduce a line break}

\item{pos_info}{whether to include position information}
}
\description{
Print binary structure map
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/units.R
\name{iso_make_units_explicit}
\alias{iso_make_units_explicit}
\title{Make units explicit}
\usage{
iso_make_units_explicit(df, prefix = " [", suffix = "]")
}
\arguments{
\item{df}{the data frame in which to make the units explicit}

\item{prefix}{the prefix for the units}

\item{suffix}{the suffix for the units}
}
\description{
This function is intended for data frames / tibbles only and makes the units of columns that have numbers with units explicit in the column name. It also strips the units attribute from those columns using \code{\link{iso_strip_units}}. The reverse function is \code{\link{iso_make_units_implicit}}.
}
\examples{
# a data frame with implicit units
df <- tibble(peak = 1:5, height = iso_double_with_units(1:5, "V"))
df

# show with explicit units
iso_make_units_explicit(df)

# show with explicit units (custom prefix & suffix)
iso_make_units_explicit(df, prefix = ".", suffix = "")
}
\seealso{
Other functions for values with units: 
\code{\link{iso_get_units}()},
\code{\link{iso_is_double_with_units}()},
\code{\link{iso_make_units_implicit}()},
\code{\link{iso_strip_units}()},
\code{\link{iso_with_units}()}
}
\concept{functions for values with units}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{iso_find_absolute_path_roots}
\alias{iso_find_absolute_path_roots}
\title{Find roots for absolute paths}
\usage{
iso_find_absolute_path_roots(path, root = ".", check_existence = TRUE)
}
\arguments{
\item{path}{vector of file/folder paths, mixed relative and absolute paths are allowed.}

\item{root}{root directory for the isofiles. Can be relative to the current working directory (e.g. \code{"data"}) or an absolute path on the file system (e.g. \code{"/Users/..."} or \code{"C:/Data/.."}). The default is the current working directory (\code{"."}). Can be supplied as a vector of same length as the provided paths if the paths have different roots.}

\item{check_existence}{whether to check for the existence of the paths}
}
\value{
a data frame with the root directories and paths relative to the root - order of input paths is preserved
}
\description{
Helper function to find the roots of absolute paths. Tries to put absolute paths into the context of the relative root. For those that this is not possible (because they are not in fact a sub-path of the relative roots), identifies the greatest common denominator for absolute paths as their root. Does not change relative paths but does check whether they do exist if \code{check_existence = TRUE} (the default). To modify relative paths, use \link{iso_shorten_relative_paths} prior to calling this function.
}
\seealso{
Other file system functions: 
\code{\link{iso_expand_paths}()},
\code{\link{iso_root_paths}()},
\code{\link{iso_shorten_relative_paths}()}
}
\concept{file system functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_data.R
\name{iso_get_file_info}
\alias{iso_get_file_info}
\title{Aggregate file info}
\usage{
iso_get_file_info(
  iso_files,
  select = everything(),
  file_specific = FALSE,
  simplify = TRUE,
  quiet = default(quiet)
)
}
\arguments{
\item{iso_files}{collection of iso_file objects}

\item{select}{which columns to select - use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns. File id is always included and cannot be renamed.}

\item{file_specific}{whether to run the select criteria (\code{...}) specifically within each individual file rather than on all files jointly. This is a lot slower but makes it possible to  select different columns in different iso_files depending on what exists in each file and is mostly of use when working with data from multiple instruments.}

\item{simplify}{if set to TRUE (the default), nested value columns in the file info will be unnested as long as they are compatible across file types. Note that file info entries with multiple values still remain nested multi-value (=list) columns even with \code{simplify=TRUE}. These can be unnested using \link[tidyr:nest]{unnest}.}

\item{quiet}{whether to display (quiet=FALSE) or silence (quiet = TRUE) information messages. Set parameter to overwrite global defaults for this function or set global defaults with calls to \link[=iso_info_messages]{iso_turn_info_messages_on} and \link[=iso_info_messages]{iso_turn_info_messages_off}}
}
\description{
Combine file information from multiple iso_files. By default all information is included but specific columns can be targeted using the \code{select} parameter to select and/or rename columns. File information beyond \code{file_id}, \code{file_root}, \code{file_path}, \code{file_datetime} and \code{file_size} (in bytes) is only available if the \code{iso_files} were read with parameter \code{read_file_info=TRUE}.
}
\note{
this function used to allow selecting/renaming different file_info_columns in different files to the same column. This was a significant speed impediment and only covered very rare use cases. It is still available in the related function \code{\link{iso_select_file_info}} with a special flag but is no longer the default and not encouraged for use in the frequently called \code{iso_get_file_info}.
}
\seealso{
Other data retrieval functions: 
\code{\link{iso_get_all_data}()},
\code{\link{iso_get_bgrd_data}()},
\code{\link{iso_get_raw_data}()},
\code{\link{iso_get_resistors}()},
\code{\link{iso_get_standards}()},
\code{\link{iso_get_vendor_data_table}()}
}
\concept{data retrieval functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unit_conversion.R
\name{iso_convert_signals}
\alias{iso_convert_signals}
\title{moved to isoprocessor}
\usage{
iso_convert_signals(...)
}
\arguments{
\item{...}{deprecated}
}
\description{
moved to isoprocessor
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{iso_plot_raw_data}
\alias{iso_plot_raw_data}
\title{moved to isoprocessor}
\usage{
iso_plot_raw_data(...)
}
\arguments{
\item{...}{deprecated}
}
\description{
moved to isoprocessor
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoread.R
\name{read_iso_file}
\alias{read_iso_file}
\title{Read individual iso file}
\usage{
read_iso_file(
  ds,
  root,
  path,
  file_n,
  files_n,
  read_from_cache,
  read_from_old_cache,
  reread_outdated_cache,
  write_to_cache,
  cachepath,
  old_cachepath,
  post_read_check,
  ext,
  reader_fun,
  reader_options,
  reader_fun_env
)
}
\arguments{
\item{ds}{the basic data structure for the type of iso_file}

\item{root}{root directory for the isofiles. Can be relative to the current working directory (e.g. \code{"data"}) or an absolute path on the file system (e.g. \code{"/Users/..."} or \code{"C:/Data/.."}). The default is the current working directory (\code{"."}). Can be supplied as a vector of same length as the provided paths if the paths have different roots.}

\item{path}{file path}

\item{file_n}{number of processed file for info messages}

\item{files_n}{total number of files for info messages}

\item{read_from_cache}{whether to read from cache}

\item{read_from_old_cache}{whether to read from old cache files (to be deprecated in isoreader 2.0)}

\item{reread_outdated_cache}{whether to reread outdated cache files}

\item{write_to_cache}{whether to write to cache}

\item{cachepath}{path for the cache file}

\item{old_cachepath}{path for the old cache files}

\item{post_read_check}{whether to run data integrity checks after a file read}

\item{ext}{file extension}

\item{reader_fun}{file reader function}

\item{reader_options}{list of parameters to be passed on to the reader}

\item{reader_fun_env}{where to find the reader function}
}
\description{
Low level read function for an individual iso file. Usually not called directly but available for methods development.
}
