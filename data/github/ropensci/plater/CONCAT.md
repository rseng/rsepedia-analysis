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
There are currently no downstream dependencies for this package.---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# plater

[![Travis-CI Build Status](https://travis-ci.org/ropensci/plater.svg?branch=master)](https://travis-ci.org/ropensci/plater)
[![CRAN version](http://www.r-pkg.org/badges/version/plater)](https://cran.r-project.org/package=plater)
[![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/plater)](http://cran.rstudio.com/web/packages/plater/index.html) [![DOI](https://zenodo.org/badge/32951641.svg)](https://zenodo.org/badge/latestdoi/32951641)
[![](https://badges.ropensci.org/60_status.svg)](https://github.com/ropensci/onboarding/issues/60)

plater makes it easy to work with data from experiments performed in plates. It is aimed at scientists and analysts who deal with microtiter plate-based instruments.  

## Installation

plater is available through CRAN. Just run: 

```{r, eval = FALSE}
install.packages("plater") 
```

## Getting your data in

Many scientific instruments (such as plate readers and qPCR machines) produce data in tabular form that mimics a microtiter plate: each cell corresponds to a well as physically laid out on the plate. For experiments like this, it's often easiest to keep records of what was what (control vs. treatment, concentration, sample type, etc.) in a similar plate layout form. 

But data in those dimensions aren't ideal for analysis. That's where `read_plate()` and `add_plate()` come in. 

* `read_plate()` takes data in plate layout form and converts it to a data frame, with one well per row, identified by well name.
* `add_plate()` does the same thing, but merges the new columns into an existing data frame you provide. 

In other words, these functions seamlessly convert plate-shaped data (easy to think about) into tidy data (easy to analyze).  

To make it even easier, if you have multiple plates in an experiment, use `read_plates()` to read them all in and combine them into a single data frame. 

## Seeing your data

Sometimes it's useful to map your data back onto a plate (are the weird outliers all from the same corner of the plate?). For that, there's `view_plate()`, which takes a data frame with one well per row, and lays it out like it's on a plate. 

## Vignette

For a detailed example of how to use `plater`, check out [the vignette.](https://cran.r-project.org/web/packages/plater/vignettes/plater-basics.html)

## Contributing to `plater`

`plater` is developed under a [Contributor Code of Conduct](CONDUCT.md). To contribute to its development, you must agree to abide by its terms. Pull requests for changes are accepted with gratitude. Please include tests as appropriate with any pull requests.  

Requests for new features and reports of bugs or security vulnerabilities can be made [here](https://github.com/ropensci/plater/issues) or emailed to the address listed [here](https://github.com/ropensci/plater/blob/master/DESCRIPTION).

[![ropensci footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)---
title: "Getting started with `plater`"
author: "Sean Hughes"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting started with `plater`}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```{r, echo = FALSE}
# for github flavored markdown, use 
# output:
#   md_document:
#     variant: markdown_github
# then switch back to output: rmarkdown::html_vignette. Just open plater-basics.md and resave to update time stamp

library(plater)

# print results of code using #>
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

## How `plater` helps you

`plater` makes it easy to work with data from experiments performed in plates.

Many scientific instruments (such as plate readers and qPCR machines) produce data in tabular form that mimics a microtiter plate: each cell corresponds to a well as physically laid out on the plate. For experiments like this, it's often easiest to keep records of what was what (control vs. treatment, concentration, sample type, etc.) in a similar plate layout form. 
But while plate-shaped data is easy to think about, it's not easy to analyze. The point of `plater` is to seamlessly convert plate-shaped data (easy to think about) into tidy data (easy to analyze). It does this by defining a simple, systematic format for storing information in plate layouts. Then it painlessly rearranges data that intuitive format into a tidy data frame.

There are just two steps: 

1. Put the data in a file in `plater` format 
2. Read in the data `plater` functions

## The example

Imagine you've invented two new antibiotics. To show how well they work, you filled up a 96-well plate with dilutions of the antibiotics and mixed in four different types of bacteria. Then, you measured how many of the bacteria got killed. So for each well in the plate you know: 

* The drug (A or B)
* The concentration of drug (100 uM to 0.01 nM and no drug)
* The bacterial species (E. coli, S. enterocolitis, C. trachomatis, and N. gonorrhoeae)
* The amount of killing in the well

The first three items are variables you chose in setting up the experiment. The fourth item is what you measured.

## Step 1: Put the data in `plater` format

The first step is to create a file for the experiment. `plater` format is designed to store all the information about an experiment in one file. It's simply a .csv file representing a single plate, containing one or more plate layouts. Each layout maps to a variable, so for the example experiment, there are four layouts in the file: Drug, Concentration, Bacteria, and Killing. 

A `plater` format file for the example experiment came with the package. Load `plater` (i.e. run `library(plater)`) and then run `system.file("extdata", package = "plater")`. Open the folder listed there and then open `example-1.csv` in a spreadsheet editor. 

An abbreviated version of that file is shown below: 

![plater format example](plater-format-image.png)

The format is pretty simple: 

* .csv file
* Top left cell of each layout is the name
* The rest of the top row of each layout is the column numbers (1:12 for a 96-well plate)
* The rest of the left column is the row names (A:H for a 96-well plate)
* One line between layouts (This row should appear as blank in a spreadsheet editor, but as a row of commas when viewed as plain text.)

You can use `plater` format with any standard plate size (6 to 1536 wells). Not every well has to be filled. If a well is blank in every layout in a file, it's omitted. If it's blank in some but not others, it'll get `NA` where it's blank.

While creating a file in `plater` format, it can be helpful to check whether you're doing it right. For that purpose, you can pass the path of the file to `check_plater_format()`, which will check that the format is correct and diagnose any problems.  

## Step 2: Read in the data

Now that your file is set up, you're ready to read in the data. 

We will analyze this experiment two different ways to illustrate two common data analysis scenarios: 

1. Assuming the instrument gives back the killing data shaped like a plate, we'll create one file with all four variables and read it in with `read_plate()`.
2. Assuming the instrument gives back tidy data (one-well-per-row), we'll create two files--one with the data and one with the three variables--and then combine the files with `add_plate()`.  

### Step 2: Read a single `plater` format file with `read_plate()`

Here is how it works. (Note that below we use `system.file()` here to get the file path of the example file, but for your own files you would specify the file path without using `system.file()`).

```{r}
file_path <- system.file("extdata", "example-1.csv", package = "plater")
   
data <- read_plate(
      file = file_path,             # full path to the .csv file
      well_ids_column = "Wells"     # name to give column of well IDs (optional)
)
str(data)

head(data)
```

So what happened? `read_plate()` read in the `plater` format file you created and turned each layout into a column, using the name of the layout specified in the file. So you have four columns: Drug, Concentration, Bacteria, and Killing. It additionally creates a column named "Wells" with the well identifiers for each well. Now, each well is represented by a single row, with the values indicated in the file for each column. 

### Step 2 (again): Combine a one-well-per-row file and a `plater` format file with `add_plate()`

In the previous example, we assumed that the killing data was provided by the instrument in plate-shaped form, so it could just be pasted into the `plater` format file. Sometimes, though, you'll get data back formatted with one well per row. 

`add_plate()` is set up to help in this situation. You provide a tidy data frame including well IDs and then you provide a `plater` format file with the other information and `add_plate()` knits them together well-by-well. Here's an example using the other two files installed along with `plater`. 

```{r}
file2A <- system.file("extdata", "example-2-part-A.csv", package = "plater")
data2 <- read.csv(file2A)

str(data2)

head(data2)

meta <- system.file("extdata", "example-2-part-B.csv", package = "plater")
data2 <- add_plate(
      data = data2,               # data frame to add to 
      file = meta,                # full path to the .csv file
      well_ids_column = "Wells"   # name of column of well IDs in data frame
)

str(data2)

head(data2)
```

`add_plate` then makes it easy to store data in a mix of formats, in some cases tidy and in some cases plate-shaped, which is the reality of many experiments. 

## Multiple plates

Say you were happy with the tests of you antibiotics, so you decided to do a second experiment, testing some other common pathogenic bacteria. Now you have data from two separate plates. Rather than handling them separately, you can combine them all into a common data frame with the `read_plates()` function. 

Just like before, you create one `plater` file per plate, with all the information describing the experiment. In this case, you'll have two files, one from each experiment. Then, just read them in with `read_plates()`. You can specify names for each plate, which will become a column in the output identifying which plate the well was on. By default it'll use the file names. 

```{r}
# same file as above
file1 <- system.file("extdata", "example-1.csv", package = "plater")

# new file
file2 <- system.file("extdata", "more-bacteria.csv", package = "plater")

data <- read_plates(
   files = c(file1, file2),
   plate_names = c("Experiment 1", "Experiment 2"),
   well_ids_column = "Wells") # optional

str(data)

head(data)
```

## Viewing plate-shaped data

Sometimes it's useful to look back at the data in plate shape. Was there something weird about that one column? Was there contamination all in one corner of the plate? 

For this, use `view_plate()` which takes a tidy data frame and displays columns from it as plate layouts. 

```{r}
view_plate(
  data = data2, 
  well_ids_column = "Wells", 
  columns_to_display = c("Concentration", "Killing")
)
```% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_plater_format.R
\name{check_plater_format}
\alias{check_plater_format}
\title{Check whether a file is in plater format.}
\usage{
check_plater_format(file)
}
\arguments{
\item{file}{The path of the file to check}
}
\value{
Displays a number of messages as it checks the file. Will stop with
a descriptive error message if the file is not formatted correctly.
}
\description{
Runs the provided file through a number of diagnostics to determine whether
it is a valid plater format file and displays information about any 
deficiencies found.
}
\examples{
file_path <- system.file("extdata", "example-1.csv", package = "plater")

data <- check_plater_format(file_path)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plater.R
\docType{package}
\name{plater}
\alias{plater}
\title{Tools to Make it Easy to Work with Microtiter Plate-Shaped Data}
\description{
plater defines a simple, plate-shaped file format for data storage, so it's 
easy to remember the experimental design. The package provides functions to 
seamlessly convert between that format and a tidy data frame that's optimal 
for analysis. \code{\link[plater]{check_plater_format}} is provided to help 
you manage plate-shaped files. \cr\cr You can work with purely plate-shaped 
data (\code{\link[plater]{read_plate}} and 
\code{\link[plater]{read_plates}}), as well as with a combination of 
plate-shaped data and tidy data (\code{\link[plater]{add_plate}}). It further
allows easy plate-shaped visualization of tidy data 
(\code{\link[plater]{view_plate}}).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/view_plate.R
\name{view_plate}
\alias{view_plate}
\title{Displays the data in the form of a microtiter plate.}
\usage{
view_plate(data, well_ids_column, columns_to_display, plate_size = 96)
}
\arguments{
\item{data}{A data frame containing the data}

\item{well_ids_column}{The name of the column in \code{data} containing the 
well IDs.}

\item{columns_to_display}{A vector of the names of one or more columns you'd
like to display.}

\item{plate_size}{The number of wells in the plate. Must be 6, 12, 24, 48, 96
384, or 1536. Default 96.}
}
\value{
A depiction of the data in \code{columns_to_display} as 
though laid out on a microtiter plate with \code{plate_size} wells.
}
\description{
Displays the data in the form of a microtiter plate.
}
\examples{
# Generate some tidy data
data <- data.frame(Wells = paste0(LETTERS[1:3], 0, rep(1:4, each = 3)), 
Species = rep(c("Alien", "Human", "Cat"), 4), 
OxygenProduction = round(rnorm(12), 3))
head(data)

# See which wells had cells from which species and the amount of oxygen 
# produced for each well
view_plate(data, "Wells", c("Species", "OxygenProduction"), 12)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_plates.R
\name{read_plates}
\alias{read_plates}
\title{Read multiple plater-formatted files and combine result into one data frame.}
\usage{
read_plates(files, plate_names = NULL, well_ids_column = "Wells")
}
\arguments{
\item{files}{A character vector with the paths of one or more plater-formatted
.csv files.}

\item{plate_names}{A character vector the same length as \code{files} with the
names to give the individual plates in the resulting data frame. Defaults to
the file names (stripped of path and .csv).}

\item{well_ids_column}{The name to give the column that will contain the well
identifiers. Default "Wells".}
}
\value{
Returns a data frame like that returned by \code{read_plate}, 
containing the data from all of the plates. The plates will be identified 
with a column called "Plate" containing the names given in 
\code{plate_names}.
}
\description{
A wrapper around \code{read_plate} that handles multiple plates and combines
them all into a single data frame.
}
\examples{
# Combine multiple files into one tidy data frame
file1 <- system.file("extdata", "example-1.csv", package = "plater")
file2 <- system.file("extdata", "more-bacteria.csv", package = "plater")

# Data are stored in plate-shaped form
data <- read_plates(
   files = c(file1, file2),
   plate_names = c("Experiment 1", "Experiment 2"),
   well_ids_column = "Wells")

# Data from both plates are tidy and in the same data frame
head(data)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_plate.R
\name{read_plate}
\alias{read_plate}
\title{Read a plater-formatted file and turn it into a tidy data frame.}
\usage{
read_plate(file, well_ids_column = "Wells")
}
\arguments{
\item{file}{The path of a .csv file formatted as described below.}

\item{well_ids_column}{The name to give the column that will contain the well
identifiers. Default "Wells".}
}
\value{
Returns a data frame with each well as a row. One column will be 
named with \code{well_ids_column} and contain the well names (A01, A02..). 
There will be as many additional columns as layouts in \code{file}. Empty 
wells are omitted.
}
\description{
Converts data from \code{plater} format to a data frame with one well 
per row identified by well name.
}
\section{\code{plater} format}{

The .csv file should be formatted as a microtiter plate. The top-left most 
cell contains the name to use for the column representing that plate. For 
example, for a 96-well plate, the subsequent wells in the top row should be 
labeled 1-12. The subsequent cells in the first column should be labeled A-H. 
That is:

\tabular{ccccc}{
ColName      \tab \strong{1} \tab \strong{2} \tab \strong{3} \tab \strong{...}\cr
\strong{A}   \tab A01        \tab A02        \tab A03        \tab ... \cr
\strong{B}   \tab B01        \tab B02        \tab B03        \tab ... \cr
\strong{...} \tab ...        \tab ...        \tab ...        \tab ... \cr
}

In this example, the cells within the plate contain the well IDs ("A01", 
"A02"), but they may contain arbitrary characters: numbers, letters, or 
punctuation. Any cell may also be blank.

Note that Microsoft Excel will sometimes include cells that appear to be 
blank in the .csv files it produces, so the files may have spurious columns
or rows outside of the plate, causing errors. To solve this problem, copy and
paste just the cells within the plate to a fresh worksheet and save it.
}

\section{Multiple columns}{
 
Multiple columns of information about a plate can be included in a single 
file. After the first plate, leave one row blank, and then add another plate
formatted as described above. (The "blank" row should appear as blank in a 
spreadsheet editor, but as a row of commas when viewed as plain text.) As 
many plates as necessary can be included in a single file (e.g. data 
measured, subject, treatment, replicate, etc.).
}

\examples{
file_path <- system.file("extdata", "example-1.csv", package = "plater")

# Data are stored in plate-shaped form
data <- read_plate(
   file = file_path,
   well_ids_column = "Wells")

# Now data are tidy
head(data)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_plate.R
\name{add_plate}
\alias{add_plate}
\title{Read a plater-formatted file and combine it with an existing data frame.}
\usage{
add_plate(data, file, well_ids_column)
}
\arguments{
\item{data}{The data frame to merge the file into. Must contain a column with
well names.}

\item{file}{The path of a .csv file formatted as described in 
\code{\link{read_plate}}.}

\item{well_ids_column}{The name of the column in \code{data} containing the 
well IDs.}
}
\value{
Returns data as a tibble with as many new columns as plates in  
\code{file}. Empty wells are indicated with NA.
}
\description{
Converts data from \code{plater} format to a data frame with one well 
per row and merges it into an existing data frame by well name.
}
\details{
If data contains more wells than in \code{file}, NA will be added to the 
merged column(s) for those wells. If the file contains more wells than 
\code{data}, those wells will be added to the bottom of the result with NA
for the columns in \code{data}.
}
\examples{
# Part of the data is tidy
file <- system.file("extdata", "example-2-part-A.csv", package = "plater")
data <- read.csv(file)

# Part of the data is plate-shaped
plate_shaped <- system.file("extdata", "example-2-part-B.csv", package = "plater")

# Combine the two
data <- add_plate(
   data = data, 
   file = plate_shaped,
   well_ids_column = "Wells")

# Now data are tidy
head(data)
}
