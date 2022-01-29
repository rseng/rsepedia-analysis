
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

