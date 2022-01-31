
<!-- README.md is generated from README.Rmd. Please edit that file -->

# baRcodeR

<!-- badges: start -->

[![Status](https://www.r-pkg.org/badges/version/baRcodeR)](https://CRAN.R-project.org/package=baRcodeR)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/baRcodeR)](https://CRAN.R-project.org/package=baRcodeR)
[![Travis build
status](https://travis-ci.org/ropensci/baRcodeR.svg?branch=master)](https://travis-ci.org/ropensci/baRcodeR)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/baRcodeR/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/baRcodeR?branch=master)
[![Project Status: Active â€“ The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN
checks](https://cranchecks.info/badges/worst/baRcodeR)](https://cran.r-project.org/web/checks/check_results_baRcodeR.html)
<!-- badges: end -->

baRcodeR generates labels for more repeatable workflows with biological
samples

## Installation

You can install the released version of baRcodeR from
[CRAN](https://CRAN.R-project.org) with:

    install.packages("baRcodeR")

And the development version from [GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("ropensci/baRcodeR", build_vignettes = T)
    # for windows users to build vignettes
    # install_github("ropensci/baRcodeR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

> NOTE: Restarting RStudio is necessary for the addin for baRcodeR to
> appear.

## Quick Start

Text identifiers can be created in a sequential or hierarchical pattern.

``` r
library(baRcodeR)
```

    ## Loading required package: qrcode

    ## Registered S3 method overwritten by 'R.oo':
    ##   method        from       
    ##   throw.default R.methodsS3

``` r
example_labels <- uniqID_maker(user = FALSE, string = "Example", level = 1:80)
head(example_labels)
```

    ##        label ind_string ind_number
    ## 1 Example001    Example        001
    ## 2 Example002    Example        002
    ## 3 Example003    Example        003
    ## 4 Example004    Example        004
    ## 5 Example005    Example        005
    ## 6 Example006    Example        006

Then the text identifiers can be printed out with a laser printer on
sticker sheets.

``` r
pdf_file_name <- tempfile()
create_PDF(Labels = example_labels, name = pdf_file_name)
```

![](man/figures/example.png)<!-- -->

Th particular layout above defaults to ULINE 1.75" \* 0.5" labels but
other layouts can be specified through parameters in the
`custom_create_PDF` function.

## Introduction

`baRcodeR` is a R package for generating unique identifier strings and
printable 2D (QR) barcodes, with the aim of improving repeatability of
labelling, tracking and curating data from biological samples.
Specifically, users can:

  - generate simple ID codes (Ex001, Ex002, Ex003 …),
  - generate hierarchical (i.e. nested) ID codes (A01-B01, A01-B02,
    A02-B01, A02-B02, A03-B01 …),
  - generate printable PDF files of paired ID codes and QR barcodes with
    default spacing for ULINE 1.75" \* 0.5" WEATHER RESISTANT LABEL for
    laser printer; item \# S-19297 (uline.ca)
  - customize the PDF layout for any type of printable format (e.g,
    vinyl stickers, waterproof paper)
  - generate reproducible code for archival purposes (e.g. in
    publications or online repositories)
  - create CSV files to link unique IDs and sampling hierarchy with
    downstream data collection workflows. For example, the PyTrackDat
    pipeline can be used to set up a web-based data collection platform:
    <https://github.com/pytrackdat/pytrackdat>

Creating unique, scannable barcodes generally involves two steps:

1.  Generate unique ID codes with `uniqID_maker()` or
    `uniqID_hier_maker()`
2.  Create a PDF file containing unique ID codes coupled with 2D barcode
    using `create_PDF()`

If you already have ID codes saved in a CSV file, the csv can be read
into a `data.frame()` in R. The `label` column, if it exists will be
used as input to generate barcodes. Otherwise, the first column in the
data frame will be used.

> NOTE: When printing from pdf, ensure that ‘anti-aliasing’ or
> ‘smoothing’ options are turned OFF, and that you are not using ‘fit
> to page’ or similar options that will re-scale the output.

![Flowchart of major functions](man/figures/Flowchart.png)

### Cheat Sheet

A 2-page, quick-reference guide is available via
[Figshare](https://dx.doi.org/10.6084/m9.figshare.7043309)

## Usage with RStudio addin

Please load the vignette “Use Addin”.

``` r
library(baRcodeR)
vignette("use-addin")
```

## Usage from the console

Please load the vignette “Using-baRcodeR” for console use.

``` r
vignette("Using-baRcodeR")
```

# Contribution

Please note that the ‘baRcodeR’ project is released with a [Contributor
Code of
Conduct](https://github.com/ropensci/baRcodeR/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

Please document issues with a description, a minimal reproducible
example, and the `sessionInfo()`.

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 17763)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252   
    ## [3] LC_MONETARY=English_Canada.1252 LC_NUMERIC=C                   
    ## [5] LC_TIME=English_Canada.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] baRcodeR_0.1.4 qrcode_0.1.1  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.3        png_0.1-7         digest_0.6.23     R.methodsS3_1.7.1
    ##  [5] magrittr_1.5      evaluate_0.14     rlang_0.4.1       stringi_1.4.3    
    ##  [9] rstudioapi_0.10   R.oo_1.22.0       R.utils_2.9.0     rmarkdown_1.17   
    ## [13] tools_3.6.1       stringr_1.4.0     xfun_0.11         yaml_2.2.0       
    ## [17] compiler_3.6.1    htmltools_0.4.0   knitr_1.26

# See also:

  - [zintr](https://github.com/carlganz/zintr) is an R interface to the C
    zint library. Use zintr if you want to create single barcode images.
    zintr does not include functions for (i) automating the creation of
    biologically-relevant, unique ID codes or (ii) customizable layouts
    for printing multiple barcodes.

  - [zint](http://zint.org.uk/) is a C library that generates a variety
    of different barcodes. Just like zintr, zint produces single barcode
    images.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# baRcodeR 0.1.6

At the request of a COVID-2019 research group, we have added an option to allow non-encoded text to appear with linear & 2D barcodes. These have been added to custom_create_PDF():

* alt_text -- adds human-readable text that is NOT encoded in the digital barcode
* denote -- characters used to denote non-encoded text

These are considered advanced features that should be used cautiously, and therefore they are not made available through the Addins GUI. 

# baRcodeR 0.1.5

Bugs:

- fixed errors in tests due to r-devel switching to using stringsAsFactor = FALSE as default
- fixed broken link in readme


# baRcodeR 0.1.4

Bugs and Improvements:

- tests have been added for user prompts and the RStudio addin with shinytest
- vignettes describing how to use the add-in and the command line are now organized separately from the restructured README
- `custom_create_PDF()` page generation should now be faster (helpful when making sheets for hundreds of labels)
- the command prompts have been restructured for a more menu-like selection for yes/no questions.
- other minor changes in addition to the major ones outlined above as suggested by reviews at rOpenSci documented [here](https://github.com/ropensci/software-review/issues/338)
- the baRcodeR package is now a part of [rOpenSci](https://ropensci.org/) and the documentation is online (here)[https://docs.ropensci.org/baRcodeR/]

# baRcodeR 0.1.3

Bugs and Improvements:

- major bug fix for linear barcodes that occasionally created unscannable barcodes.
- added documentation on how to create alternative formatting of labels (e.g. spaces, line breaks)
- added padding for labels which were single character or blank
- 2-page cheatsheet now available as addin 

# baRcodeR 0.1.2

New Feature:

- In response to a user request, there is now an option to print linear (code 128 set B) barcodes. 

# baRcodeR 0.1.1

Bugs and Improvements: 

- create_PDF() function will replace all underscores in text with dashes. Underscores are not specified in the encoding dictionary of `qcrode` and will throw errors.
- x_space and y_space parameters are now limited between 0 and 1 for easier use. These parameters are used to position text on the printed labels.
- Font size is no longer limited and is now measured as points. Font size is automatically reduced if text code is too long for the printed labels.

New Features:

- label_width and label_height parameters specify the width and height of the label to enable alleys (i.e. gaps) between physical labels.


-----------------


# baRcodeR 0.1.0

This is the first official release of the package.# Contributor Code of Conduct

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
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
### baRcodeR 0.1.6 - new parameters

At the request of a COVID-2019 research group, we have added an option to allow non-encoded text to appear with linear & 2D barcodes:

* alt_text -- adds human-readable text that is NOT encoded in the digital barcode
* denote -- characters used to denote non-encoded text

## baRcodeR 0.1.5 R CMD check results

This is a minor update to make sure tests pass on r-devel and fix a broken link. 


## Test environments

* local Windows 10 install, R 3.6.1, r-devel
* win-builder (devel and release)
* Ubuntu 16.04 (on travis-ci), R-oldrel, R-release, R-devel
* R-hub 
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran


## Amendment to previous rejected submission

We have removed the license file from the description file and the package file itself as suggested by Uwe Ligges.


  
## baRcodeR 0.1.4 R CMD check results

There were no ERRORS or WARNINGS. 1 NOTE regarding change of package maintainer from Robert Colautti to Yihan Wu. 

This is our third submission. From our last submission, we made

1. changes to clarify documentation 

2. modify command line prompts with menu-like choices

3. added a vignette tutorial to show how to use the package with the RStudio addin and removed similar content from the package README

4. corrected error in current CRAN version of package relating to R-devel changes

### April 26, 2019

This is our second submission. From first submission on September 10, 2018, we made three changes:

1. We have changed the example for create_PDF() to write to tempdir() as suggested by Uwe Liggs.

2. We renamed functions containing label_ to uniqID_ to clarify that the functions generate unique identifiers, in contrast to create_PDF, which generates labels for printing. We found this was already a source of confusions among early testers.

3. We fixed a minor bug in the preview of one of the Addins tabs and we clarified explanations in the vignette and documentation.

### Oct 4, 2018 

additional changes advised by Swetlana Herbrandt

1. Removed 'Open-Source R Tools for' in the title
2. Replaced T and F with TRUE and FALSE throughout
3. Repaced \dontrun{} with if(interactive()){} or \donttest{}

### Oct 12, 2018

A few changes from Oct 4 were not included in the last upload.
Additionally, fixed a minor text errors:
1. 'collection' occurred twice in DESCRIPTION
2. 'unique' was mis-spelled in documentation of uniqID_maker

### Nov 29, 2018

Bug fixes and minor features added, as outlined in NEWS.md

### Jan 10, 2019

- Added the ability to create linear (1D) barcodes using the type='linear' parameter in create_pdf()
- fixed a few minor spelling/grammar changes

### April 26, 2019

Bug fix for linear barcodes and minor changes outlined in NEWS.md (ver. 0.1.3)

---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

# baRcodeR

<!-- badges: start -->
[![Status](https://www.r-pkg.org/badges/version/baRcodeR)](https://CRAN.R-project.org/package=baRcodeR)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/baRcodeR)](https://CRAN.R-project.org/package=baRcodeR)
[![Travis build status](https://travis-ci.org/ropensci/baRcodeR.svg?branch=master)](https://travis-ci.org/ropensci/baRcodeR)
[![Codecov test coverage](https://codecov.io/gh/ropensci/baRcodeR/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/baRcodeR?branch=master)
[![Project Status: Active â€“ The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN checks](https://cranchecks.info/badges/worst/baRcodeR)](https://cran.r-project.org/web/checks/check_results_baRcodeR.html)
[![rOpenSci peer-review](https://badges.ropensci.org/338_status.svg)](https://github.com/ropensci/software-review/issues/338)
<!-- badges: end -->

baRcodeR generates labels for more repeatable workflows with biological samples

## Installation

You can install the released version of baRcodeR from [CRAN](https://CRAN.R-project.org) with:

```
install.packages("baRcodeR")
```

And the development version from [GitHub](https://github.com/) with:

```
# install.packages("devtools")
devtools::install_github("ropensci/baRcodeR", build_vignettes = T)
# for windows users to build vignettes
# install_github("ropensci/baRcodeR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```

> NOTE: Restarting RStudio is necessary for the addin for baRcodeR to appear.

## Quick Start

Text identifiers can be created in a sequential or hierarchical pattern.

```{r}
library(baRcodeR)

example_labels <- uniqID_maker(user = FALSE, string = "Example", level = 1:80)
head(example_labels)
```

Then the text identifiers can be printed out with a laser printer on sticker sheets.


```{r include = F, eval = F}
# generate the image shown in vignettes as png rather than pdf so it shows up in markdown
# pdf_file_name <- "man/figures/example_file"
# create_PDF(Labels = example_labels, name = pdf_file_name)
# pdf_file <- magick::image_read_pdf("man/figures/example_file.pdf")
# magick::image_write(pdf_file, path = "man/figures/example.png", format = "png")
```


```{r eval = F}
pdf_file_name <- tempfile()
create_PDF(Labels = example_labels, name = pdf_file_name)
```

```{r echo = F}
knitr::include_graphics("man/figures/example.png")
```


Th particular layout above defaults to  ULINE 1.75" * 0.5"  labels but other layouts can be specified through parameters in the `custom_create_PDF` function. 


## Introduction

`baRcodeR` is a R package for generating unique identifier strings and printable 2D (QR) barcodes, with the aim of improving repeatability of labelling, tracking and curating data from biological samples. Specifically, users can:

* generate simple ID codes (Ex001, Ex002, Ex003 ...),
* generate hierarchical (i.e. nested) ID codes (A01-B01, A01-B02, A02-B01, A02-B02, A03-B01 ...),
* generate printable PDF files of paired ID codes and QR barcodes with default spacing for ULINE 1.75" * 0.5" WEATHER RESISTANT LABEL for laser printer; item # S-19297 (uline.ca)
* customize the PDF layout for any type of printable format (e.g, vinyl stickers, waterproof paper)
* generate reproducible code for archival purposes (e.g. in publications or online repositories)
* create CSV files to link unique IDs and sampling hierarchy with downstream data collection workflows. For example, the PyTrackDat pipeline can be used to set up a web-based data collection platform: https://github.com/pytrackdat/pytrackdat

Creating unique, scannable barcodes generally involves two steps:

  1. Generate unique ID codes with `uniqID_maker()` or `uniqID_hier_maker()`
  2. Create a PDF file containing unique ID codes coupled with 2D barcode using `create_PDF()`

If you already have ID codes saved in a CSV file, the csv can be read into a `data.frame()` in R. The `label` column, if it exists will be used as input to generate barcodes. Otherwise, the first column in the data frame will be used. 

> NOTE: When printing from pdf, ensure that 'anti-aliasing' or 'smoothing' options are turned OFF, and that you are not using 'fit to page' or similar options that will re-scale the output.

![Flowchart of major functions](man/figures/Flowchart.png)

### Cheat Sheet

A 2-page, quick-reference guide is available via [Figshare](https://dx.doi.org/10.6084/m9.figshare.7043309)

## Usage with RStudio addin

Please load the vignette "Use Addin".

```{r eval = F}
library(baRcodeR)
vignette("use-addin")
```


## Usage from the console

Please load the vignette "Using-baRcodeR" for console use.

```{r eval = F} 
vignette("Using-baRcodeR")
```
# Contribution

Please note that the 'baRcodeR' project is released with a
  [Contributor Code of Conduct](https://github.com/ropensci/baRcodeR/blob/master/CODE_OF_CONDUCT.md).
  By contributing to this project, you agree to abide by its terms.
  
Please document issues with a description, a minimal reproducible example, and the `sessionInfo()`. 

```{r}
sessionInfo()
```

# See also:

- [zintr](https://github.com/carlganz/zintr)is an R interface to the C zint library. Use zintr if you want to create single barcode images. zintr does not include functions for (i) automating the creation of biologically-relevant, unique ID codes or (ii) customizable layouts for printing multiple barcodes.

- [zint](http://zint.org.uk/) is a C library that generates a variety of different barcodes. Just like zintr, zint produces single barcode images.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

---
title: "Using baRcodeR with command line prompts"
author: "Yihan Wu and Robert I. Colautti"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Quick Start to Use baRcodeR}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

BaRcodeR is an open-source program to facilitate repeatable label generation and management for labelling, tracking and curating data from biological samples.

Flowchart of major functions
<img src="https://raw.githubusercontent.com/yihanwu/baRcodeR/master/man/figures/Flowchart.png" alt="drawing" width="650"/>

For a quick start, see the [introduction](https://docs.ropensci.org/baRcodeR).


# Cheat Sheet

A 2-page, quick-reference guide is available via [Figshare](https://dx.doi.org/10.6084/m9.figshare.7043309)

# Overview

Creating unique, scannable 2-D barcodes involves two steps:

  1. Generate unique ID codes with `uniqID_maker()` or `uniqID_hier_maker()`
  2. Create a PDF file containing unique ID codes coupled with 2D barcode using `create_PDF()`

If you already have ID codes saved in a CSV file, skip step #1 and import directly as a `data.frame()` object for step #2: see [Create barcodes](#create-barcodes).

For using the baRcodeR RStudio addin, refer to the package vignette *Use baRcodeR addin* through `vignette("use-addin", package="baRcodeR")`. 

The functions in `baRcodeR` will accept parameters in two ways: 1) On the command line as part of the function, and 2) interactively through user input when `user=TRUE` is specified as part of the function. This vignette will cover the two different ways to obtain the same text identifiers and PDF sheets needed for printing.


# Create simple ID codes {#uniqID-maker}

Simple ID codes can be generated using the `uniqID_maker()` function. One-level ID codes consist of two parts, a leading string shared by all samples and a number unique to each sample. 

For example, if we want to generate basic ID codes for 5 individuals:

```
example005
example006
example007
example008
example009
example010
```

First, load the `baRcodeR` library.

```{r eval=T, include=T}

library(baRcodeR)

```

## By user prompt {#uniqID-maker-user-prompt}

Run the `uniqID_make(r))` function in interactive mode.

```{r eval=F, include=T}
IDcodes <- uniqID_maker(user = TRUE)
```

> NOTE: When typing into the console for user prompts, strings should not be quoted. 

User prompts will appear in the console. Inputs received from the user are passed as parameters to the `uniqID_maker()` to create the object `IDcodes`, which is a `data.frame()` object containing a vector of unique IDs in the first column, and additional columns for strings and individual numbers provided by the user. 


> `Please enter string for level: `


The string entered here defines the leading part of the ID code that will be the same for each individual code. This can be the name of a population, a species or a location. In this example, the string entered by the user into the console as denoted by the `>` symbol is "example".


> `Please enter string for level: ` example


The second user prompt is:

> `Please enter the starting number for level (integer): `

The prompt asks for the starting number to generate unique IDs. These codes do not have to start from 1. 

> `Please enter the starting number for level (integer):` 5

The third user prompt is:

> `Enter the ending number for level: `

The prompt asks for the ending number. Unique IDs are generated sequentially from the starting number to the ending number. Note that a higher number can be used as the starting number to generate a reverse order. It is possible to generate IDs that are not sequentially numbered by passing a vector through `uniqID_maker()` (e.g. `seq()`) (see [By arguments](#IDcode-maker-arguments)) 

> `Enter the ending number for level: ` 10

After the starting and ending numbers are set, the function generates a series of numbers. The user is then asked for the number of digits in the unique ID. 

> `Number of digits to print for level:`

This number must be >= the maximum number of digits in the unique ID code and will add leading zeros as needed. This is particularly useful for generating a smaller number ID codes that are expected to be part of a much larger sample set.

> `Number of digits to print for level: ` 3


## By arguments {#IDcode-maker-arguments}

It is also possible to create unique ID codes directly, without user prompts, by defining parameters in `uniqID_maker`. The interactive example above can be reproduced a single line of code:

```{r eval=T, include=T}
IDcodes <- uniqID_maker(string = "example", level = 5:10, digits = 3)
IDcodes
```

More complicated, non-sequential ID codes can be generated via the `levels` parameter:

```{r eval=T, include=T}
number_sequence <- seq(1, 10, 2)
IDcodes <- uniqID_maker(string = "example", level = number_sequence, digits = 3)
```

The output will then be:

```{r eval=T, include=F}
IDcodes
```

ID codes can be saved to a text file for use with other programs (e.g. spreadsheet):

```{r eval=F, include=T}
write.csv(IDcodes, "IDcodes.csv")
```


# Create hierarchical ID codes {#uniqID_hier_maker}

`uniqID_hier_maker` is used to make unique ID codes that follow a hierarchical (i.e. nested) structure, for example we might sample B individuals from A populations at each of C time points. Similar to `uniqID_maker`, this function can be run interactively or by directly defining parameters. In contrast to `uniqID_maker`, `uniqID_hier_maker` is used to generate nested pairs of strings and unique ID codes. Below is an example of a list of hierarchical identifier codes with three levels (a, b, c) and varying numbers of individuals for each level (a=3, b=2, c=2). 

```
a1-b1-c1
a1-b1-c2
a1-b2-c1
a1-b2-c2
a2-b1-c1
a2-b1-c2
a2-b2-c1
a2-b2-c2
a3-b1-c1
a3-b1-c2
a3-b2-c1
a3-b2-c2
```



## By user prompts {#uniqID-hier-maker-user-prompt}

To create hierarchical ID codes in interactive mode, start with the argument `user=T` in the function.

```{r eval=F, include=T}
IDcodes <- uniqID_hier_maker(user = TRUE)
```

The first prompt that appears in the console is:

> `What is the number of levels in hierarchy: `

In this example, we have levels a, b, and c; so three levels in total. 

> `What is the # of levels in hierarchy: ` 3

The second prompt asks if a string should be appended to the end of the IDs. 

> `String at end of label:`

> `1 : Yes`

> `2 : No` 

There are only two possible inputs, 1 or 2. Typing any other number will result in an "invalid input" warning.

In this example, there is no ending string.

> `String at end of label:  `

>` 1 : Yes ` 

>` 2 : No ` 

> 2

A series of prompts will repeat for each level, allowing the user to set the number of digits to be printed, the leading string, the starting number, and the ending number. These are similar to `uniqID_maker`(see [uniqID_maker](#uniqID-maker-user-prompt) for step by step instructions). The number of digits to print applies to all levels.

## By argument

Instead of interactive mode, it is possible to define arguments directly into `uniqID_hier_maker`. Unlike `uniqID_maker`, a `list` object is required to specify the parameters in each level of the hierarchy. First, define a vector for each level with three parameters: the leading string, the start value, and the end value. Second, combine vectors into a `list` object in the order of their hierarchy. For the example below, 'c' is nested within 'b' within 'a'. 

```{r}
level_one <- c("a", 1, 3)
level_two <- c("b", 1, 2)
level_three <- c("c", 10, 12)
hier_list <- list(level_one, level_two, level_three)
```

You can specify a custom suffix string for all ID codes using the `end` argument, and the number of digits to print for all levels through the `digits` argument. It is not possible at this time to vary the number of digits printed for each level, however this can be done using interactive mode (i.e. user=T).

The list can then be passed into `uniqID_hier_maker` to generate the unique ID codes.

```{r, eval=T, include=T}
IDcodes <- uniqID_hier_maker(hierarchy = hier_list, digits = 1)
```

The output will be:

```{r eval=T, include=T}
IDcodes
```

The data frame contains ID codes in the first column, and a separate column for each level of the hierarchy, with the user-defined string as the header. This can be saved to a CSV:

```{r eval=F, include=T}
write.csv(IDcodes, "IDcodes.csv")
```

This file is useful for archiving ID codes and as a starting point for data entry. For example, it can be opened in a spreadsheet program to add data measurement columns. It is also the input for creating printable, QR-coded labels with `create_PDF`.


# Create barcodes {#create-barcodes}

2D barcodes (i.e. QR codes) paired with ID code strings are created from an input vector of text labels. Users can manually create their own ID codes as a vector, as the first column of an existing `data.frame()` object, or as a `data.frame()` from `<-uniqID_maker` or `<-uniqID_hier_maker`. 

The function `create_PDF()` produces a pdf file containing barcodes that fit to the dimensions of ULINE 1.75" * 0.5" WEATHER RESISTANT LABEL for laser printer; item # S-19297 (uline.ca). If needed, the page setup can be modified using advanced options in [custom_create_PDF](#custom-create-pdf). 

The first step is to read in the vector of ID codes, in this example from a CSV file: 

```{r eval=F, include=T}
# Reading in from a csv file
IDcodes<-read.csv("IDcodes.csv")
```

## By user prompt {#create-barcodes-user-prompt}

In the following example, the IDcodes data.frame object is used to create a PDF file called "example.pdf", with a font size of 3.5, and an error correction level "Q" meaning the barcode can tolerate up to 25% damage. The parameter `user=T` will prompt the user to guide creation of the pdf file containing scannable barcodes.

If `IDcodes` is a vector, the vector will be directly used to generate barcodes. If `IDcodes` is a data frame, the function will use the column called `label` or else the first column in the data frame. 

```{r eval=F, include=T}
create_PDF(user=TRUE, Labels=IDcodes)
```
A user prompt is printed into the console. For example:

> `Please enter name for PDF output file: `

Any combination of letters and numbers can be used for the name of the pdf file. Here, the file name is set to "example."

> `Please enter name for PDF output file: ` example

The next user prompt is to set the size of the text printed on each barcode. 

> `Please enter a font size: `

This font size is the point size and will apply to the size of the text on the PDF. For example, entering 12 will create 12 point font on the PDF.

> `Please enter a font size: ` 12

The last basic parameter to set is the error correction level. There are four possible levels: L, M, Q, and H. 

Level "L" - up to 7% damage -- ideal for very small labels (large pixels)
Level "M" - up to 15% damage
Level "Q" - up to 25% damage
Level "H" - up to 30% damage -- good for bigger labels (small pixels)

The user prompt for error correction level is similar to previous prompts.

> `Select an error correction level. `

> `1 : L (up to 7% damage)`

> `2 : M (up to 15% damage)`

> `3 : Q (up to 25% damage)`

> `4 : H (up to 30% damage) `

This example uses an error correction level "Q" so `3` should be entered.

> `Select an error correction level. `

> `1 : L (up to 7% damage)`

> `2 : M (up to 15% damage)`

> `3 : Q (up to 25% damage)`

> `4 : H (up to 30% damage) `

> 3

The last user prompt asks whether the user wants to modify the advanced parameters. 

> `Edit advanced parameters?`

> `1 : Yes`

> `2 : No`


In this example, no advanced parameters are modified (input `2`). Using advanced parameters are covered in [advanced options](#custom-create-pdf)

> `Edit advanced parameters?`

> `1 : Yes`

> `2 : No`

> 2

## By arguments

The same example above can be reproduced directly with the following parameters:

```{r eval=F, include=T}
create_PDF(Labels = IDcodes, name = "example", ErrCorr = "Q", Fsz = 2.5)
```


# Advanced Options for pdf output {#custom-create-pdf}

There are advanced options for the pdf output which can be accessed interactively or by specifying additional arguments in `create_PDF`. The user prompts are similar to the ones shown above but allow customization of the output document for other printing formats. Documentation of the advanced options can be found using through the man page `?custom_create_PDF`.

Arguments can be passed from `create_PDF` to `custom_create_PDF` as `create_PDF` is just a wrapper for `custom_create_PDF`. 

```{r eval=F, include=T}

## This will create a pdf sheet where the labels are printed in columns then rows. It will skip 3 rows from the top and 1 column from the left. 
create_PDF(Labels = Labels, name = "example_advanced", ErrCorr = "Q", Fsz = 2.5, Across = F, ERows = 3, ECol = 1)
```


# Formatting of text labels

Trying to format text labels in the baRcodeR GUI is not recommended. 

## New Lines

It is possible to force formatting of labels by inserting `\n` as line breaks.

## Tabs

Using `\t` will create error due to the underlying `qrcode` library used in this package. One solution is to use the space character in place of tabs.

```{r}
# original label
X <- "text\ttext"
cat(X)
cat(gsub("\\t", "\x20\x20\x20\x20", X))
```
## When reading from file

R automatically escapes the escape characters when reading from file. For example, `text\ntext` will be read in as `text\\ntext` from `read.csv`. 

Use a global substitution to replace the double slashes.

```{r}
X <- "text\\ntext"
cat(X)
cat(gsub("\\\\n","\n",X))
```

---
title: "Use baRcodeR addin"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Use baRcodeR addin}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r}
knitr::opts_chunk$set(echo = F)
```

BaRcodeR is an open-source program to facilitate repeatable label generation and management for labelling, tracking and curating data from biological samples.

Flowchart of major functions
<img src="https://raw.githubusercontent.com/ropensci/baRcodeR/master/man/figures/Flowchart.png" alt="drawing" width="650"/>

For a quick start, see the [introduction](https://docs.ropensci.org/baRcodeR).

# Cheat Sheet

A 2-page, quick-reference guide is available via [Figshare](https://dx.doi.org/10.6084/m9.figshare.7043309)

If RStudio is not available, see the [introduction](https://docs.ropensci.org/baRcodeR) and `vignette("Using-baRcodeR)"` for command line use. 

## Using the RStudio addin 

The main baRcodeR functions for unique identifiers and QR code generation can be performed interactive via the RStudio addin found on the toolbar. 

### Find the addin

Make sure to restart RStudio after installing. Then the addin should appear in the toolbar. 

Click on the add-in, and a popup window will appear. 

```{r, fig.cap = "Screenshot of RStudio addins bar"}
knitr::include_graphics("add-in-screenshot.png")
```

Note the 3 tabs along the bottom, corresponding to the three main baRcodeR commands: `uniqID_maker`, `uniqID_hier_maker` and `create_PDF`.


```{r, fig.cap = "Screenshot of the simple ID Code tab"}
knitr::include_graphics("tab-1-screenshot.png")
```

### Generate simple ID codes

The first tab generates basic ID codes with user input as seen below:

```{r, fig.cap = "Active simple ID code tab"}
knitr::include_graphics("tab-1-screenshot-2.png")
```


As you fill in the fields, a preview of the ID codes will appear on the right-hand side along with reproducible code, which can be copied for archival purposes. Clicking 'Create Label.csv' will create a CSV file called 'Label_YYYY-MM-DD.csv', which contains a data frame with the full unique ID strings as the first column, the user-defined prefix string in the second column, and the unique ID number in the third column. This file is useful for archiving ID codes and as a starting point for data entry. For example, it can be opened in a spreadsheet program to add data measurement columns. It is also the input for creating printable, QR-coded labels with `create_PDF`.

```{r, fig.cap = "Screenshot of the hierarchical ID code tab"}
knitr::include_graphics("tab-2-screenshot.png")
```


### Generate Hierarchical ID codes

You can switch from the simple ID code generation tab to the hierarchical ID code generation or QR code creation tabs at the bottom.

Hierarchical ID codes have a nested structure (e.g. X subsamples from Y individuals at Z time points), the information for each level is saved under the "Hierarchy" section. The "Add level" button is used to add more levels to the hierarchy, and the "Remove level" button will remove the most recently added level. The data frame output will contain ID codes in the first column, and a separate column for each level of the hierarchy, with the user-defined string as the header; as shown under 'Preview'. As with the simple ID code tab, the output of Hierarchical ID codes is a CSV file "Labels_YYYY-MM-DD.csv", saved in the working directory. This file is useful for archiving ID codes and as a starting point for data entry. For example, it can be opened in a spreadsheet program to add data measurement columns. It is also the input for creating printable, QR-coded labels with `create_PDF`.

```{r, fig.cap = "Sceenshot of PDF creation tab"}
knitr::include_graphics("tab-3-screenshot.png")
```

### Create the PDF for sticker printing

The Barcode Creation tab contains all the advanced options for page layout. The default options fit a specific format: ULINE 1.75" * 0.5" WEATHER RESISTANT LABEL for laser printer; item # S-19297 (uline.ca). A text file containing ID codes is imported by clicking the "Browse" button and selecting the CSV text file in the file browser. The file is be previewed by clicking "Import File". 

After importing a CSV file, the preview shows part of the expected output PDF file based on font size and other layout options. The first column is highlighted by default and defines the column to use for the labels. Clicking on a different column will set it as the ID code column, as shown in the preview.

```{r, fig.cap = "Screenshot of Column Selection"}
knitr::include_graphics("tab-3-screenshot-2.png")
```


Clicking "Make PDF" will generate a printable PDF of all barcodes provided. This can take several minutes for >100 barcodes, depending on computer speed. The text "Done" will appear upon completion of the PDF file.

> NOTE: When printing from pdf, ensure that 'anti-aliasing' or 'smoothing' options are turned OFF, and that you are not using 'fit to page' or similar options that will re-scale the output.

```{r}
sessionInfo()
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniqID_maker_addin.R
\name{cheatsheet}
\alias{cheatsheet}
\title{baRcodeR Cheatsheet}
\usage{
cheatsheet()
}
\value{
Opens webpage of PDF
}
\description{
This addin links to a downloadable PDF version of the baRcodeR cheatsheet.
}
\examples{
if(interactive()){
baRcodeR::cheatsheet()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createPDF.R
\name{create_PDF}
\alias{create_PDF}
\title{Make barcodes and print labels}
\usage{
create_PDF(
  user = FALSE,
  Labels = NULL,
  name = "LabelsOut",
  type = "matrix",
  ErrCorr = "H",
  Fsz = 12,
  ...
)
}
\arguments{
\item{user}{logical. Run function using interactive mode (prompts user for
parameter values) Default is \code{FALSE}}

\item{Labels}{vector or data frame object containing label names (i.e. unique
ID codes) with either UTF-8 or ASCII encoding.}

\item{name}{character. Name of the PDF output file. Default is
\code{"LabelsOut"}. A file named \code{name.pdf} will be saved to the
working directory by default. Use \code{"dirname/name"} to produce a file
called \code{name.pdf} in the \code{dirname} directory.}

\item{type}{character. Choice of \code{"linear"} code 128 or \code{"matrix"}
QR code labels. Default is \code{"matrix"}.}

\item{ErrCorr}{error correction value for matrix labels only. Level of damage
from low to high: \code{"L"}, \code{"M"}, \code{"Q"}, \code{"H"}. Default
is \code{"H"}. See details for explanation of values.}

\item{Fsz}{numerical. Sets font size in points. Longer ID codes may be shrunk
to fit if truncation is not used for matrix labels. Default font size is
\code{5}. ID codes are also shrunk automatically to fit on the label if
actual size is bigger than label dimensions.}

\item{...}{advanced arguments to modify the PDF layout. See
\code{\link{custom_create_PDF}} for arguments. The advanced options can be
 accessed interactively with \code{user = TRUE} and then entering TRUE when prompted to
  modify advanced options.}
}
\value{
a PDF file containing QR-coded labels, saved to the default directory.
}
\description{
Input vector or data.frame of ID codes to produce a PDF of QR codes which can
be printed. This is a wrapper function for \code{\link{custom_create_PDF}}.
See details of \code{\link{custom_create_PDF}} on how to format text labels
if needed.
}
\details{
The default PDF setup is for ULINE 1.75" * 0.5" WEATHER RESISTANT LABEL for laser
printer; item # S-19297 (uline.ca). The page format can be modified using
the \code{...} (advanced arguments) for other label types.
}
\examples{
## data frame
example_vector <- as.data.frame(c("ao1", "a02", "a03"))

\dontrun{
## run with default options
## pdf file will be "example.pdf" saved into a temp directory

temp_file <- tempfile()

create_PDF(Labels = example_vector, name = temp_file)

## view example output from temp folder
system2("open", paste0(temp_file, ".pdf"))
}

## run interactively. Overrides default pdf options
if(interactive()){
    create_PDF(user = TRUE, Labels = example_vector)
}

\dontrun{
## run using a data frame, automatically choosing the "label" column
example_df <- data.frame("level1" = c("a1", "a2"), "label" = c("a1-b1",
"a1-b2"), "level2" = c("b1", "b1"))
create_PDF(user = FALSE, Labels = example_df, name = file.path(tempdir(), "example_2"))
}

\dontrun{
## run using an unnamed data frame
example_df <- data.frame(c("a1", "a2"), c("a1-b1", "a1-b2"), c("b1", "b1"))
## specify column from data frame
create_PDF(user = FALSE, Labels = example_df[,2], name = file.path(tempdir(), "example_3"))
}
\dontrun{
## create linear (code128) label rather than matrix (2D/QR) labels
example_df <- data.frame(c("a1", "a2"), c("a1-b1", "a1-b2"), c("b1", "b1"))
## specify column from data frame
create_PDF(user = FALSE, Labels = example_df, name = file.path(tempdir(),
"example_4", type = "linear"))
}
}
\seealso{
\code{\link{custom_create_PDF}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hidden_createPDF.R
\name{custom_create_PDF}
\alias{custom_create_PDF}
\alias{qrcode_make}
\alias{code_128_make}
\title{Make barcodes and print labels}
\usage{
custom_create_PDF(
  user = FALSE,
  Labels = NULL,
  name = "LabelsOut",
  type = "matrix",
  ErrCorr = "H",
  Fsz = 12,
  Across = TRUE,
  ERows = 0,
  ECols = 0,
  trunc = TRUE,
  numrow = 20,
  numcol = 4,
  page_width = 8.5,
  page_height = 11,
  width_margin = 0.25,
  height_margin = 0.5,
  label_width = NA,
  label_height = NA,
  x_space = 0,
  y_space = 0.5,
  alt_text = NULL,
  replace_label = FALSE,
  denote = c("\\n(", ")")
)

qrcode_make(Labels, ErrCorr)

code_128_make(Labels)
}
\arguments{
\item{user}{logical. Run function using interactive mode (prompts user for
parameter values) Default is \code{FALSE}}

\item{Labels}{vector or data frame object containing label names (i.e. unique
ID codes) with either UTF-8 or ASCII encoding.}

\item{name}{character. Name of the PDF output file. Default is
\code{"LabelsOut"}. A file named \code{name.pdf} will be saved to the
working directory by default. Use \code{"dirname/name"} to produce a file
called \code{name.pdf} in the \code{dirname} directory.}

\item{type}{character. Choice of \code{"linear"} code 128 or \code{"matrix"}
QR code labels. Default is \code{"matrix"}.}

\item{ErrCorr}{error correction value for matrix labels only. Level of damage
from low to high: \code{"L"}, \code{"M"}, \code{"Q"}, \code{"H"}. Default
is \code{"H"}. See details for explanation of values.}

\item{Fsz}{numerical. Sets font size in points. Longer ID codes may be shrunk
to fit if truncation is not used for matrix labels. Default font size is
\code{5}. ID codes are also shrunk automatically to fit on the label if
actual size is bigger than label dimensions.}

\item{Across}{logical. When \code{TRUE}, print labels across rows, left to
right. When \code{FALSE}, print labels down columns, top to bottom. Default
is \code{TRUE}.}

\item{ERows}{number of rows to skip. Default is \code{0}. Example: setting
ERows to 6 will begin printing at row 7. ERows and ECols are useful for
printing on partially-used label sheets.}

\item{ECols}{number of columns to skip. Default is \code{0}. Example: setting
ECols to 2 will put the first label at column 3. ERows and ECols are useful
for printing on partially-used label sheets.}

\item{trunc}{logical. Text is broken into multiple lines for longer ID codes,
to prevent printing off of the label area. Default is \code{TRUE}. If
\code{trunc = FALSE}, and text is larger than the physical label, the text will
be shrunk down automatically.}

\item{numrow}{numerical. Number of rows per page. Default is \code{20}.}

\item{numcol}{numerical. Number of columns per page. Default is \code{4}.}

\item{page_width}{numerical. Width of page (in inches). Default is set to
\code{8.5}.}

\item{page_height}{numerical. Height of page (in inches). Default is set to
\code{11}.}

\item{width_margin}{numerical. The width margin of the page (in inches).
Default is \code{0.25}.}

\item{height_margin}{numerical. The height margin of the page (in inches).
Default is \code{0.5}.}

\item{label_width}{numerical. The width of label (in inches). Will be
calculated as \code{(page_width - 2 * width_margin)/numcol} if
\code{label_width} is set as \code{NULL}.}

\item{label_height}{numerical. The height of the label (in inches). Will be
calculated as \code{(page_height - 2 * height_margin)/numrow} if
\code{label_height} is set as \code{NULL}.}

\item{x_space}{numerical. A value between \code{0} and \code{1}. This sets
the distance between the QR code and text of each label. Only applies when
\code{type = "matrix"}. Default is \code{0}.}

\item{y_space}{numerical. The height position of the text on the physical
label as a proportion of the label height. Only applies when \code{type =
"matrix"}. A value between \code{0} and \code{1}. Default is \code{0.5}.}

\item{alt_text}{vector containing alternative names that are printed along with 
Labels BUT ARE NOT ENCODED in the barcode image. Use with caution!}

\item{replace_label}{logical. Replace label text with \code{alt_text}.
Generated barcode will contain more information than text label. Use with
caution!}

\item{denote}{character (prefix) or vector of length 2 (prefix, suffix). 
Denotes alt_text that is not encoded in the barcode image. 
Default is brackets before and after ().}
}
\value{
a PDF file containing QR-coded labels, saved to the default
  directory.
}
\description{
Input a vector or data frame of ID codes to produce a PDF of barcode labels
that can then be printed. The PDF setup is for the ULINE 1.75" * 0.5" WEATHER
RESISTANT LABEL for laser printer; item # S-19297 (uline.ca). See details for
how to format text labels properly.
}
\details{
\code{qrcode_make} is the helper function for generating a QR code matrix.
\code{code_128_make} is the helper function for generating a linear barcode
according to code 128 set B. \code{custom_create_PDF} is the main function
which sets page layout, and creates the PDF file.

Correction levels for QR codes refer to the level of damage a label can
tolerate before the label become unreadable by a scanner (L = Low (7\%), M =
Medium (15\%), Q = Quantile (25\%), H = High (30\%)). So a label with L
correction can lose up to at most 7% of the code before it is unreadable
while a H label can lose up to 30% of the code. This also means that L codes
can be printed at smaller sizes compared to H codes.

The escape characters \code{\\n} and \code{\\s} (and the hex equivalents
\code{\\x0A} and \code{\\x20} can be used to format text labels. Tab character
\code{\\t} (\code{\\x09}) does not work for QR codes and should be replaced by
a number of space characters. See the package vignette for examples.

If \code{ECol} or \code{ERow} is greater than \code{numcol} and \code{numrow}, 
the labels will be printed starting on the second page.
}
\examples{

## this is the same examples used with create_PDF
## data frame
example_vector <- as.data.frame(c("ao1", "a02", "a03"))

\dontrun{
## run with default options
## pdf file will be "example.pdf" saved into a temp directory
temp_file <- tempfile()

custom_create_PDF(Labels = example_vector, name = temp_file)

## view example output from temp folder
system2("open", paste0(temp_file, ".pdf"))
}

## run interactively. Overrides default pdf options
if(interactive()){
    custom_create_PDF(user = TRUE, Labels = example_vector)
}

\dontrun{
## run using a data frame, automatically choosing the "label" column
example_df <- data.frame("level1" = c("a1", "a2"), "label" = c("a1-b1",
 "a1-b2"), "level2" = c("b1", "b1"))

custom_create_PDF(user = FALSE, Labels = example_df, name = file.path(tempdir(), 
 "example_2"))
 }
\dontrun{
## run using an unnamed data frame
example_df <- data.frame(c("a1", "a2"), c("a1-b1", "a1-b2"), c("b1", "b1"))
## specify column from data frame
custom_create_PDF(user = FALSE, Labels = example_df[,2], name = file.path(tempdir(), "example_3"))
}
\dontrun{
## create linear (code128) label rather than matrix (2D/QR) labels
example_df <- data.frame(c("a1", "a2"), c("a1-b1", "a1-b2"), c("b1", "b1"))
## specify column from data frame
custom_create_PDF(user = FALSE, Labels = example_df, name = file.path(tempdir(),
"example_4", type = "linear"))
}
\dontrun{
## Include text for the user that is NOT encoded into the barcode image
## Excluded text is denoted with brackets by default
example_df <- data.frame(ID = floor(runif(3) * 10000), name = c("A", "B", "C"),
 dob = c("1/1/2020", "12/6/2001", "2/8/1986")
 
## linear (1d) barcodes with custom denote parameter
custom_create_PDF(Labels = example_df$ID, alt_text = paste(example_df$name,
 example_df$dob), type = "linear", denote=".")
}  
}
\seealso{
\code{\link{create_PDF}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniqID_maker_addin.R
\name{make_labels}
\alias{make_labels}
\title{baRcodeR GUI}
\usage{
make_labels()
}
\value{
Opens RStudio addin gadget window for making labels and barcodes in a GUI
}
\description{
This addin will allow you to interactive create ID codes and generate PDF files
of QR codes.
}
\examples{
if(interactive()){
library(baRcodeR)
make_labels()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniqIDMaker.R
\name{uniqID_maker}
\alias{uniqID_maker}
\title{Generate a list of ID codes}
\usage{
uniqID_maker(
  user = FALSE,
  string = NULL,
  level,
  digits = 3,
  ending_string = NULL
)
}
\arguments{
\item{user}{logical. Run function using interactive mode (prompts user for 
parameter values). Default is \code{FALSE}}

\item{string}{character. Text string for label. Default \code{null}.}

\item{level}{integer vector. Defines the numerical values to be appended
to the character string. Can be any sequence of numbers (see examples).}

\item{digits}{numerical. Default is \code{2}. Number of digits to be printed, 
adding leading 0s as needed. This will apply to all levels when \code{user=FALSE}. 
When the numeric value of the label has a greater number of digits than 
\code{digits}, \code{digits} is automatically increased for the entire level. 
Default is \code{3}.}

\item{ending_string}{a character string or vector of strings to attach to the label.
If a vector is used, all combinations of that vector with a unique label will be produced.}
}
\value{
data.frame with text labels in the first column, along with string
and numeric values in two additional columns.
}
\description{
Create ID codes consisting of a text string and unique numbers (string001, string002, ...). 
Can be run in interactive mode, prompting user for input. The data.frame 
output can be saved as CSV for (i) the \code{\link{create_PDF}} function 
to generate printable QR-coded labels; and (ii) to downstream data 
collection software (spreadsheets, relational databases, etc.)
}
\details{
When the function is called with \code{user = TRUE}, a sequence of 
numbers is generated between the starting and ending number provided by the 
user. When \code{user = FALSE}, a vector of custom numbers can be provided. 
See example below.
}
\examples{


## sequential string of numbers in label
Labels <- uniqID_maker(string = "string", level = c(1:5), digits = 2)
Labels

## can also use nonsequential strings in input for levels
level <- c(1:5, 8:10, 999:1000)
Labels <- uniqID_maker(string = "string", level = level, digits = 4)
Labels


## Using the ending_string to produce labels with unique endings
## this is different from hierarchical labels with two levels as there 
## is no numbering, just the text string


Labels <- uniqID_maker(string = "string", level = c(1:5), digits = 2, ending_string = "A")
Labels

Labels <- uniqID_maker(string = "string", level = c(1:5), 
                       digits = 2, ending_string = c("A", "B"))
Labels


if(interactive()){
## function using user prompt does not use any of the other parameters
Labels <- uniqID_maker(user = TRUE)
Labels
}
}
\seealso{
\code{\link{uniqID_hier_maker}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniqIDHierarchy.R
\name{uniqID_hier_maker}
\alias{uniqID_hier_maker}
\title{Make hierarchical ID codes}
\usage{
uniqID_hier_maker(user = FALSE, hierarchy, end = NULL, digits = 2)
}
\arguments{
\item{user}{logical. Run function using interactive mode (prompts user for 
parameter values). Default is \code{FALSE}}

\item{hierarchy}{list. A list with each element consisting of three members
a vector of three elements (string, beginning value, end value). See examples.
Used only when \code{user=FALSE})}

\item{end}{character. A string to be appended to the end of each label.}

\item{digits}{numerical. Default is \code{2}. Number of digits to be printed, 
adding leading 0s as needed. This will apply to all levels when \code{user=FALSE}. 
When the max number of digits in the ID code is greater than number of digits 
defined in \code{digits}, then \code{digits} is automatically increased 
to avoid errors.}
}
\value{
data.frame of text labels in the first column, with additional columns 
for each level in the hierarchy list, as defined by the user.
}
\description{
Generate hierarchical ID codes for barcode labels. 
Hierarchical codes have a nested structure: e.g. Y subsamples from 
each of X individuals. Use \code{\link{uniqID_maker}} 
for sequential single-level labels. Can be run in interactive mode, 
prompting user for input. The data.frame can be saved as CSV for 
(i) the \code{\link{create_PDF}} function to generate printable 
QR-coded labels; and (ii) to downstream data collection using spreadsheet, 
relational database, etc.
}
\examples{
if(interactive()){
## for interactive mode
uniqID_hier_maker(user = TRUE)
}

## how to make hierarchy list

## create vectors for each level in the order string_prefix, beginning_value,
## end_value and combine in list

a <- c("a", 3, 6)
b <- c("b", 1, 3)
c <- list(a, b)
Labels <- uniqID_hier_maker(hierarchy = c)
Labels

## add string at end of each label
Labels <- uniqID_hier_maker(hierarchy = c, end = "end")
Labels

}
\seealso{
\code{\link{uniqID_maker}}
}
