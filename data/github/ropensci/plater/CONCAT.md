---
title: 'plater: Read, Tidy, and Display Data from Microtiter Plates'
tags:
 - data import
 - R
authors:
- name: Sean M Hughes
  orcid: 0000-0002-9409-9405
  affiliation: University of Washington
date: 27 September 2016
bibliography: paper.bib
---

# Summary

plater is an R [@R] package that makes it easy to work with data from experiments performed in microtiter plates.

Many scientific instruments (such as plate readers and qPCR machines) produce data in tabular form that mimics a microtiter plate: each cell corresponds to a well as physically laid out on the plate. For experiments like this, it's often easiest to keep records of what was what (control vs. treatment, concentration, sample type, etc.) in a similar plate layout form. 

plater defines a simple, plate-shaped file format for data storage, so it's easy to remember the experimental design, and provides functions to seamlessly convert between that format and a tidy [@tidy] data frame that's optimal for analysis. When the instrument produces data that's already tidy, plater helps combine that data with plate-shaped experimental metadata. Once the data is tidy, it's sometimes useful to look back at it in plate shape, so plater makes that easy, too. 

# References
# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
<!-- README.md is generated from README.Rmd. Please edit that file -->
plater
======

[![Travis-CI Build Status](https://travis-ci.org/ropensci/plater.svg?branch=master)](https://travis-ci.org/ropensci/plater) [![CRAN version](http://www.r-pkg.org/badges/version/plater)](https://cran.r-project.org/package=plater) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/plater)](http://cran.rstudio.com/web/packages/plater/index.html) [![DOI](https://zenodo.org/badge/32951641.svg)](https://zenodo.org/badge/latestdoi/32951641) [![](https://badges.ropensci.org/60_status.svg)](https://github.com/ropensci/onboarding/issues/60)

plater makes it easy to work with data from experiments performed in plates. It is aimed at scientists and analysts who deal with microtiter plate-based instruments.

Installation
------------

plater is available through CRAN. Just run:

``` r
install.packages("plater") 
```

Getting your data in
--------------------

Many scientific instruments (such as plate readers and qPCR machines) produce data in tabular form that mimics a microtiter plate: each cell corresponds to a well as physically laid out on the plate. For experiments like this, it's often easiest to keep records of what was what (control vs. treatment, concentration, sample type, etc.) in a similar plate layout form.

But data in those dimensions aren't ideal for analysis. That's where `read_plate()` and `add_plate()` come in.

-   `read_plate()` takes data in plate layout form and converts it to a data frame, with one well per row, identified by well name.
-   `add_plate()` does the same thing, but merges the new columns into an existing data frame you provide.

In other words, these functions seamlessly convert plate-shaped data (easy to think about) into tidy data (easy to analyze).

To make it even easier, if you have multiple plates in an experiment, use `read_plates()` to read them all in and combine them into a single data frame.

Seeing your data
----------------

Sometimes it's useful to map your data back onto a plate (are the weird outliers all from the same corner of the plate?). For that, there's `view_plate()`, which takes a data frame with one well per row, and lays it out like it's on a plate.

Vignette
--------

For a detailed example of how to use `plater`, check out [the vignette.](https://cran.r-project.org/web/packages/plater/vignettes/plater-basics.html)

Contributing to `plater`
------------------------

`plater` is developed under a [Contributor Code of Conduct](CONDUCT.md). To contribute to its development, you must agree to abide by its terms. Pull requests for changes are accepted with gratitude. Please include tests as appropriate with any pull requests.

Requests for new features and reports of bugs or security vulnerabilities can be made [here](https://github.com/ropensci/plater/issues) or emailed to the address listed [here](https://github.com/ropensci/plater/blob/master/DESCRIPTION).

[![ropensci footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# plater 1.0.3
* Change `add_plate()` to return a tibble rather than trying to preserve initial class
* Remove use of deprecated `select_` function

# plater 1.0.2
* Changes to tests to comply with new CRAN policy on `data.frame(..., stringsAsFactors = FALSE)`
* Add support for 6- and 1536-well plates
* Change behavior of add_plate so that when the plate layout contains more wells than the input data frame, those wells are appended to the end of the data frame instead of erroring. 

# plater 1.0.1
* Eliminate warnings from readLines on files without EOF (Mac issue)
* Fix issue with numeric formatting in mixed numeric/character layouts
* Fix issue with grouped tibbles and view_plate

# plater 1.0.0 (5 Oct 2016)
* Changes in response to rOpenSci reviewers
* Reorder arguments of `add_plate()` for better pipelining
* add `check_plater_format()` to help with preparing files
* rename all lowercase

# plateR 0.2.1
* Reorganize parameters for consistency
* Add defaults for parameters
* Add `read_plates()`

# plateR 0.2
* Introduce new data format with multiple plate layouts per .csv file (replacing multiple files at once)

# plateR 0.1
* Add support for reading multiple files at once# Update, version 1.0.3, 4 Jan 2021

This is a minor update: 

* Fix bug where class of some objects was mishandled
* Replace internal use of deprecated function

## Test environments
* ubuntu 18.04 on travis-ci:  devel   2021-01-02 r79767 
                              release 4.0.2 (2020-06-22)
* win-builder:                devel   2021-01-02 r79767
                              release 4.0.3 (2020-10-10)
* rhub                        
    * ubuntu 16.04            3.6.1
    * fedora                  2020-10-24 r79367
    * windows                 2020-12-14 r79633

## R CMD check results
There were no ERRORs or WARNINGs. 

There was one NOTE:

   * checking CRAN incoming feasibility ... NOTE_to_CRAN_maintainers
   Maintainer: 'Sean Hughes <smhughes@uw.edu>'

## Downstream dependencies
There are currently no downstream dependencies for this package.