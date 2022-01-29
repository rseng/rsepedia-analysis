
# historydata

[![Build Status](https://travis-ci.org/ropensci/historydata.svg)](https://travis-ci.org/ropensci/historydata)

## Overview

These sample data sets are intended for historians learning R. They
include population, institutional, religious, military, and
prosopographical data suitable for mapping, quantitative analysis, and
network analysis.

## Installation

To install the package from CRAN:

    install.packages("historydata")

To install the development version, you will first have to install
[devtools][] and then install this package from GitHub.

To install:

    devtools::install_github("ropensci/historydata")

## Use

To list all the datasets in the package with their documentation:

    library(historydata)
    help(package = historydata)

To load a dataset:

    data(catholic_dioceses)

## Contributing

If you have a dataset that you think would be good for this package,
feel free to contribute it. You can send the dataset via e-mail if you
can provide a citation and guarantee that the data is available under an
open license.

But it is much preferred that you contribute the dataset via a pull
request. To add a dataset:

1.  Add the raw data and an R script to load and save it as an R data
    object to `data-raw/`. The `.rda` file should be saved to `data/`.
    See `data-raw/sarna.R` as a model. Please try keeping the data
    [tidy][].
2.  Add the documentation using the [Roxygen][] format to a file in the
    `R/`. Use `R/sarna.R` as a model. Be sure to include a citation.

## License

This project is released under the MIT License:
<http://lmullen.mit-license.org/>

  [devtools]: https://github.com/hadley/devtools
  [tidy]: https://doi.org/10.18637/jss.v059.i10
  [Roxygen]: http://roxygen.org/
# historydata 0.2 (in development)

-   Added `dijon_prices` and `dijon_prices_wide`
-   Added `presbyterians`
-   Fixed data formatting of `paulist_missions`
-   Added `methodists`

# historydata 0.1

-   Added `sarna`: Population estimates for American Jews.
-   Added us\_state\_populations\`: Populations of US states and territories.
-   Added `catholic_dioceses`: Locations and dates founded for Roman Catholic dioceses in the United States, Mexico, and Canada.
-   Added `naval_promotions`: Career data for officers of the line in the early U.S. Navy.
-   Added `us_national_population`: Populations of the United States of America.
-   Added `early_colleges`: Colleges founded in the United States before 1848.
-   Added `paulist_missions`: Missions held by the Paulists Fathers,1851-1893.
-   Added `judges_people` and `judges_appointments`: Federal judges, 1789--present.
-   Added `tudors`: relationships among the Tudor dynasty.
This is a resubmission of the historydata package version 0.1. The package description has been updated as suggested by CRAN maintainers.

This package has been checked on Mac, Linux (Ubuntu 14.04), and Windows  without warnings or notes. It is a new release of a data package which  will be infrequently updated.  
