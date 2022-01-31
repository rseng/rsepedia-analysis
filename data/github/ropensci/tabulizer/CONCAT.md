
# Extract Tables from PDFs

[![CRAN](https://www.r-pkg.org/badges/version/tabulizer)](https://cran.r-project.org/package=tabulizer)
[![Downloads](https://cranlogs.r-pkg.org/badges/tabulizer)](https://cran.r-project.org/package=tabulizer)
[![Build
Status](https://travis-ci.org/ropensci/tabulizer.png?branch=master)](https://travis-ci.org/ropensci/tabulizer)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ropensci/tabulizer?branch=master&svg=true)](https://ci.appveyor.com/project/tpaskhalis/tabulizer)
[![codecov.io](https://codecov.io/github/ropensci/tabulizer/coverage.svg?branch=master)](https://codecov.io/github/ropensci/tabulizer?branch=master)
[![](https://badges.ropensci.org/42_status.svg)](https://github.com/ropensci/onboarding/issues/42)

**tabulizer** provides R bindings to the [Tabula java
library](https://github.com/tabulapdf/tabula-java/), which can be used
to computationaly extract tables from PDF documents.

Note: tabulizer is released under the MIT license, as is Tabula itself.

## Installation

tabulizer depends on [rJava](https://cran.r-project.org/package=rJava),
which implies a system requirement for Java. This can be frustrating,
especially on Windows. The preferred Windows workflow is to use
[Chocolatey](https://chocolatey.org/) to obtain, configure, and update
Java. You need do this before installing rJava or attempting to use
tabulizer. More on [this](#installing-java-on-windows-with-chocolatey)
and [troubleshooting](#troubleshooting) below.

To install the latest CRAN version:

``` r
install.packages("tabulizer")
```

To install the latest development version:

``` r
if (!require("remotes")) {
    install.packages("remotes")
}
# on 64-bit Windows
remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")
# elsewhere
remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"))
```

## Code Examples

The main function, `extract_tables()` provides an R clone of the Tabula
command line application:

``` r
library("tabulizer")
f <- system.file("examples", "data.pdf", package = "tabulizer")
out1 <- extract_tables(f)
str(out1)
## List of 4
##  $ : chr [1:32, 1:10] "mpg" "21.0" "21.0" "22.8" ...
##  $ : chr [1:7, 1:5] "Sepal.Length " "5.1 " "4.9 " "4.7 " ...
##  $ : chr [1:7, 1:6] "" "145 " "146 " "147 " ...
##  $ : chr [1:15, 1] "supp" "VC" "VC" "VC" ...
```

By default, it returns the most table-like R structure available: a
matrix. It can also write the tables to disk or attempt to coerce them
to data.frames using the `output` argument. It is also possible to
select tables from only specified pages using the `pages`
argument.

``` r
out2 <- extract_tables(f, pages = 1, guess = FALSE, output = "data.frame")
str(out2)
## List of 1
##  $ :'data.frame':       33 obs. of  13 variables:
##   ..$ X   : chr [1:33] "Mazda RX4 " "Mazda RX4 Wag " "Datsun 710 " "Hornet 4 Drive " ...
##   ..$ mpg : num [1:33] 21 21 22.8 21.4 18.7 18.1 14.3 24.4 22.8 19.2 ...
##   ..$ cyl : num [1:33] 6 6 4 6 8 6 8 4 4 6 ...
##   ..$ X.1 : int [1:33] NA NA NA NA NA NA NA NA NA NA ...
##   ..$ disp: num [1:33] 160 160 108 258 360 ...
##   ..$ hp  : num [1:33] 110 110 93 110 175 105 245 62 95 123 ...
##   ..$ drat: num [1:33] 3.9 3.9 3.85 3.08 3.15 2.76 3.21 3.69 3.92 3.92 ...
##   ..$ wt  : num [1:33] 2.62 2.88 2.32 3.21 3.44 ...
##   ..$ qsec: num [1:33] 16.5 17 18.6 19.4 17 ...
##   ..$ vs  : num [1:33] 0 0 1 1 0 1 0 1 1 1 ...
##   ..$ am  : num [1:33] 1 1 1 0 0 0 0 0 0 0 ...
##   ..$ gear: num [1:33] 4 4 4 3 3 3 3 4 4 4 ...
##   ..$ carb: int [1:33] 4 4 1 1 2 1 4 2 2 4 ...
```

It is also possible to manually specify smaller areas within pages to
look for tables using the `area` and `columns` arguments to
`extract_tables()`. This facilitates extraction from smaller portions of
a page, such as when a table is embeded in a larger section of text or
graphics.

Another function, `extract_areas()` implements this through an
interactive style in which each page of the PDF is loaded as an R
graphic and the user can use their mouse to specify upper-left and
lower-right bounds of an area. Those areas are then extracted
auto-magically (and the return value is the same as for
`extract_tables()`). Here’s a shot of it in action:

![extract\_areas()](https://i.imgur.com/USTyQl7.gif)

`locate_areas()` handles the area identification process without
performing the extraction, which may be useful as a debugger.

`extract_text()` simply returns text, possibly separately for each
(specified) page:

``` r
out3 <- extract_text(f, page = 3)
cat(out3, sep = "\n")
## len supp dose
## 4.2 VC 0.5
## 11.5 VC 0.5
## 7.3 VC 0.5
## 5.8 VC 0.5
## 6.4 VC 0.5
## 10.0 VC 0.5
## 11.2 VC 0.5
## 11.2 VC 0.5
## 5.2 VC 0.5
## 7.0 VC 0.5
## 16.5 VC 1.0
## 16.5 VC 1.0
## 15.2 VC 1.0
## 17.3 VC 1.0
## 22.5 VC 1.0
## 3
```

Note that for large PDF files, it is possible to run up against Java
memory constraints, leading to a `java.lang.OutOfMemoryError: Java heap
space` error message. Memory can be increased using
`options(java.parameters = "-Xmx16000m")` set to some reasonable amount
of memory.

Some other utility functions are also provided (and made possible by the
Java [Apache PDFBox library](https://pdfbox.apache.org/)):

  - `extract_text()` converts the text of an entire file or specified
    pages into an R character vector.
  - `split_pdf()` and `merge_pdfs()` split and merge PDF documents,
    respectively.
  - `extract_metadata()` extracts PDF metadata as a list.
  - `get_n_pages()` determines the number of pages in a document.
  - `get_page_dims()` determines the width and height of each page in pt
    (the unit used by `area` and `columns` arguments).
  - `make_thumbnails()` converts specified pages of a PDF file to image
    files.

### Installing Java on Windows with Chocolatey

In command prompt, install Chocolately if you don’t already have
    it:

    @powershell -NoProfile -ExecutionPolicy Bypass -Command "iex ((new-object net.webclient).DownloadString('https://chocolatey.org/install.ps1'))" && SET PATH=%PATH%;%ALLUSERSPROFILE%\chocolatey\bin

Then, install java using Chocolately’s `choco install` command:

    choco install jdk7 -y

You may also need to then set the `JAVA_HOME` environment variable to
the path to your Java installation (e.g., `C:\Program
Files\Java\jdk1.8.0_92`). This can be done:

1.  within R using `Sys.setenv(JAVA_HOME = "C:/Program
    Files/Java/jdk1.8.0_92")` (note slashes), or
2.  from command prompt using the `setx` command: `setx JAVA_HOME
    C:\Program Files\Java\jdk1.8.0_92`, or
3.  from PowerShell, using the .NET framework:
    `[Environment]::SetEnvironmentVariable("JAVA_HOME", "C:\Program
    Files\Java\jdk1.8.0_92", "User")`, or
4.  from the Start Menu, via `Control Panel » System » Advanced »
    Environment Variables` ([instructions
    here](http://superuser.com/a/284351/221772)).

You should now be able to safely open R, and use rJava and tabulizer.
Note, however, that some users report that rather than setting this
variable, they instead need to delete it (e.g., with
`Sys.setenv(JAVA_HOME = "")`), so if the above instructions fail, that
is the next step in troubleshooting.

### Troubleshooting

Some notes for troubleshooting common installation problems:

  - On Mac OS, you may need to install [a particular version of
    Java](https://support.apple.com/kb/DL1572?locale=en_US) prior to
    attempting to install tabulizer.
  - On a Unix-like, you need to ensure that R has been installed with
    Java support. This can often be fixed by running `R CMD javareconf`
    on the command line (possibly with `sudo`, etc. depending on your
    system setup).
  - On Windows, make sure you have permission to write to and install
    packages to your R directory before trying to install the package.
    This can be changed from “Properties” on the right-click context
    menu. Alternatively, you can ensure write permission by choosing
    “Run as administrator” when launching R (again, from the
    right-click context menu).

## Meta

  - Please [report any issues or
    bugs](https://github.com/ropensci/tabulizer/issues).
  - License: MIT
  - Get citation information for `tabulizer` in R doing
    `citation(package =
'tabulizer')`

[![rofooter](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# CHANGES TO tabulizer 0.2.2

* `extract_tables()` gets `outdir` argument for writing out CSV, TSV and JSON
files.
* Fixes in vignette.

# CHANGES TO tabulizer 0.2.1

* `make_thumbnails()` and `split_pdf()` now use `tempdir()` as the default
output directory.
* `extract_` functions get `copy` argument for copying original local files to 
R session's temporary directory.
* PDF files used do not get copied to temporary directory by default.
* General clean-up for CRAN submission.

# CHANGES TO tabulizer 0.2.0

* Upgrade to PDFBox 2/Tabula 1.0.1 ([#48](https://github.com/ropensci/tabulizer/issues/48))
* `method` argument is changed to `output` in `extract_tables()`.
* New `method` argument reflects method of extraction as in Tabula command-line Java utility.
* `extract_text()` accepts `area` as argument.

# CHANGES TO tabulizer 0.1.24

* Switch to using new document loading algorithm and localize all remote URLs. ([#40](https://github.com/ropensci/tabulizer/issues/40))
* Removed kind of annoying message about overwriting a temporary file. ([#36](https://github.com/ropensci/tabulizer/issues/36))

# CHANGES TO tabulizer 0.1.23

* Transferred source code repository to ropensci: https://github.com/ropensci/tabulizer

# CHANGES TO tabulizer 0.1.22

* Transferred source code repository to ropenscilabs: https://github.com/ropenscilabs/tabulizer

# CHANGES TO tabulizer 0.1.21

* Exposed a new option `widget` in `locate_areas()` to control which widget is used in locating areas. 

# CHANGES TO tabulizer 0.1.20

* Fixed bug in internal function `try_area_full()` introduced by changes in8.

# CHANGES TO tabulizer 0.1.19

* Further expand the `locate_areas()` interface to use a Shiny gadget when working within RStudio, or otherwise rely on the full functionality interface (based on graphics device events) or reduced functionality interface (relying on `locator()`). (#8)

# CHANGES TO tabulizer 0.1.18

* Completely rewrite the `locate_areas()` interface to rely on graphics device event handling where possible. This may behave differently across platforms or in RStudio. (#8)

# CHANGES TO tabulizer 0.1.17

* Fixed a bug in `extract_tables()` such that when no tables are found, an empty list is returned (for `method` values with list response structures). (h/t Lincoln Mullen)

# CHANGES TO tabulizer 0.1.16

* `split_pdfs()` and `make_thumbnails()` gain an `outdir` argument to specify where to save the output. The file numbering of output files is also now zero-padded.
* An warning in `merge_pdfs()` has been fixed.
* `stop_logging()` is called when the package is attached to the search path.
* `get_page_dims()` earns a `doc` argument and argument order in `get_n_pages()` is reversed.

# CHANGES TO tabulizer 0.1.15

* Added better support for specifying character encoding. (#10)

# CHANGES TO tabulizer 0.1.14

* Added support for password-protected PDF files. (#11)

# CHANGES TO tabulizer 0.1.13

* Expand file paths where needed. (h/t David Gohel)

# CHANGES TO tabulizer 0.1.12

* Improve handling of URLs when using `extract_areas()` by downloading PDF to temporary directory.

# CHANGES TO tabulizer 0.1.11

* Exposed new functions `split_pdf()` and `merge_pdfs()` to split and merge PDFs, respectively. (#9)
* Exposed a `get_n_pages()` to determine the page length of a PDF document.
* Moved tabula .jar file to separate package, tabulizerjars. (#2)

# CHANGES TO tabulizer 0.1.10

* Added a new function `extract_metadata()` to extract PDF metadata as a list.
* Added a new function `extract_text()` to convert PDF contents to an R character vector.
* Changed the internal `localize_file()` function to use PDFBox to natively read from a URL.
* Removed illogical default `file` argument value in `extract_tables()`.

# CHANGES TO tabulizer 0.1.9

* Expanded test suite to cover `areas` and `columns` arguments and utilities. (#3)
* Fixed the same bug in `make_columns()` as was corrected for `make_areas()`. (#5)

# CHANGES TO tabulizer 0.1.8

* Fixed a bug in `make_areas()` internal when `area` was specified as a length 1 list for a multi-page document. (#5, h/t Tony Hirst)

# CHANGES TO tabulizer 0.1.7

* Added a function, `extract_areas()`, to interactively identify and extract page areas. Another new function, `locates_areas()` implements the locator functionality without performing any extraction.
* Added a function, `make_thumbnails()`, to convert pages into individual image files.
* Added a function, `get_page_dims()`, to extract page dimensions.

# CHANGES TO tabulizer 0.1.6

* Fixed a bug in the repeating of the `area` argument when `length(area) == 1 & length(pages) > 1`. (#5, #6)

# CHANGES TO tabulizer 0.1.5

* Fixed a bug in the handling of the `area` argument. (#5, #6)

# CHANGES TO tabulizer 0.1.4

* Added vignette. (#4)
* Added tests of guess parameter. (#3)
* Added a `spreadsheet` argument, a la Tabula itself.
* Fixed bugs in parsing of `area` and `columns` arguments.

# CHANGES TO tabulizer 0.1.3

* Added multiple table writing options beyond the default list of matrices. (#1)

# CHANGES TO tabulizer 0.1.1

* Initial release.
Contributions to **tabulizer** are welcome from anyone and are best sent as pull requests on [the GitHub repository](https://github.com/leeper/tabulizer/). This page provides some instructions to potential contributors about how to add to the package.

 1. Contributions can be submitted as [a pull request](https://help.github.com/articles/creating-a-pull-request/) on GitHub by forking or cloning the [repo](https://github.com/leeper/tabulizer/), making changes and submitting the pull request.
 
 2. Pull requests should involve only one commit per substantive change. This means if you change multiple files (e.g., code and documentation), these changes should be committed together. If you don't know how to do this (e.g., you are making changes in the GitHub web interface) just submit anyway and the maintainer will clean things up.
 
 3. All contributions must be submitted consistent with the package license ([MIT](https://opensource.org/licenses/MIT)).
 
 4. Non-trivial contributions need to be noted in the `Authors@R` field in the [DESCRIPTION](https://github.com/leeper/tabulizer/blob/master/DESCRIPTION). Just follow the format of the existing entries to add your name (and, optionally, email address). Substantial contributions should also be noted in [`inst/CITATION`](https://github.com/leeper/tabulizer/blob/master/inst/CITATION).
 
 5. The project uses royxgen code and documentation markup, so changes should be made to roxygen comments in the source code `.R` files. If changes are made, roxygen needs to be run. The easiest way to do this is a command line call to: `Rscript -e devtools::document()`. Please resolve any roxygen errors before submitting a pull request.
 
 6. Please run `R CMD BUILD tabulizer` and `R CMD CHECK tabulizer_VERSION.tar.gz` before submitting the pull request to check for any errors.
 
Some specific types of changes that you might make are:

 1. Bug fixes. Great!
 
 2. Documentation-only changes (e.g., to Rd files, README, vignettes). This is great! All contributions are welcome.
 
 3. New functionality. This is fine, but should be discussed on [the GitHub issues page](https://github.com/leeper/tabulizer/issues) before submitting a pull request.
 
 3. Changes requiring a new package dependency should also be discussed on [the GitHub issues page](https://github.com/leeper/tabulizer/issues) before submitting a pull request.
 
 4. Message translations. These are very appreciated! The format is a pain, but if you're doing this I'm assuming you're already familiar with it.

Any questions you have can be opened as GitHub issues or directed to thosjleeper (at) gmail.com.
\pagenumbering{gobble}

To cite R in publications use:

  R Core Team (2018). R: A language and environment for statistical computing. R
  Foundation for Statistical Computing, Vienna, Austria. URL

  https://www.R-project.org/.

A BibTeX entry for LaTeX users is

```
  @Manual{,
    title = {R: A Language and Environment for Statistical Computing},
    author = {{R Core Team}},
    organization = {R Foundation for Statistical Computing},
    address = {Vienna, Austria},
    year = {2018},
    url = {https://www.R-project.org/},
  }
```

We have invested a lot of time and effort in creating R, please cite it when
using it for data analysis. See also ‘citation("pkgname")’ for citing R
packages.

\newpage

To cite R in publications use:

  R Core Team (2018). R: A language and environment for statistical computing. R
  Foundation for Statistical Computing, Vienna, Austria. URL

  https://www.R-project.org/.

A BibTeX entry for LaTeX users is

```
  @Manual{,
    title = {R: A Language and Environment for Statistical Computing},
    author = {{R Core Team}},
    organization = {R Foundation for Statistical Computing},
    address = {Vienna, Austria},
    year = {2018},
    url = {https://www.R-project.org/},
  }
```

We have invested a lot of time and effort in creating R, please cite it when
using it for data analysis. See also ‘citation("pkgname")’ for citing R
packages.
Please specify whether your issue is about:

 - [ ] a possible bug
 - [ ] a question about package functionality
 - [ ] a suggested code or documentation change, improvement to the code, or feature request

If you are reporting (1) a bug or (2) a question about code, please supply:

 - ensure that you can install and successfully load [**rJava**](https://cran.r-project.org/package=rJava)
 - [a fully reproducible example](http://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example) using a publicly available dataset (or provide your data)
 - if an error is occurring, include the output of `traceback()` run immediately after the error occurs
 - the output of `sessionInfo()`

Put your code here:

```R
## rJava loads successfully
# install.packages("rJava")
library("rJava")

## load package
library("tabulizer")

## code goes here


## session info for your system
sessionInfo()
```

Please ensure the following before submitting a PR:

 - [ ] if suggesting code changes or improvements, [open an issue](https://github.com/leeper/responserates/issues/new) first
 - [ ] for all but trivial changes (e.g., typo fixes), add your name to [DESCRIPTION](https://github.com/leeper/responserates/blob/master/DESCRIPTION)
 - [ ] for all but trivial changes (e.g., typo fixes), documentation your change in [NEWS.md](https://github.com/leeper/responserates/blob/master/NEWS.md) with a parenthetical reference to the issue number being addressed
 - [ ] if changing documentation, edit files in `/R` not `/man` and run `devtools::document()` to update documentation
 - [ ] add code or new test files to [`/tests`](https://github.com/leeper/responserates/tree/master/tests/testthat) for any new functionality or bug fix
 - [ ] make sure `R CMD check` runs without error before submitting the PR

---
output: github_document
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##"
)
```

# Extract Tables from PDFs

[![CRAN](https://www.r-pkg.org/badges/version/tabulizer)](https://cran.r-project.org/package=tabulizer)
[![Downloads](https://cranlogs.r-pkg.org/badges/tabulizer)](https://cran.r-project.org/package=tabulizer)
[![Build Status](https://travis-ci.org/ropensci/tabulizer.png?branch=master)](https://travis-ci.org/ropensci/tabulizer)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/tabulizer?branch=master&svg=true)](https://ci.appveyor.com/project/tpaskhalis/tabulizer)
[![codecov.io](https://codecov.io/github/ropensci/tabulizer/coverage.svg?branch=master)](https://codecov.io/github/ropensci/tabulizer?branch=master)
[![](https://badges.ropensci.org/42_status.svg)](https://github.com/ropensci/onboarding/issues/42)

**tabulizer** provides R bindings to the [Tabula java library](https://github.com/tabulapdf/tabula-java/), which can be used to computationaly extract tables from PDF documents.

Note: tabulizer is released under the MIT license, as is Tabula itself.

## Installation

tabulizer depends on [rJava](https://cran.r-project.org/package=rJava), which implies a system requirement for Java. This can be frustrating, especially on Windows. The preferred Windows workflow is to use [Chocolatey](https://chocolatey.org/) to obtain, configure, and update Java. You need do this before installing rJava or attempting to use tabulizer. More on [this](#installing-java-on-windows-with-chocolatey) and [troubleshooting](#troubleshooting) below.

To install the latest CRAN version:
```{r eval = FALSE}
install.packages("tabulizer")
```

To install the latest development version:
```{r eval = FALSE}
if (!require("remotes")) {
    install.packages("remotes")
}
# on 64-bit Windows
remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")
# elsewhere
remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"))
```

## Code Examples

The main function, `extract_tables()` provides an R clone of the Tabula command line application:

```{r eval = FALSE}
library("tabulizer")
f <- system.file("examples", "data.pdf", package = "tabulizer")
out1 <- extract_tables(f)
str(out1)
## List of 4
##  $ : chr [1:32, 1:10] "mpg" "21.0" "21.0" "22.8" ...
##  $ : chr [1:7, 1:5] "Sepal.Length " "5.1 " "4.9 " "4.7 " ...
##  $ : chr [1:7, 1:6] "" "145 " "146 " "147 " ...
##  $ : chr [1:15, 1] "supp" "VC" "VC" "VC" ...
```

By default, it returns the most table-like R structure available: a matrix. It can also write the tables to disk or attempt to coerce them to data.frames using the `output` argument. It is also possible to select tables from only specified pages using the `pages` argument.

```{r eval = FALSE}
out2 <- extract_tables(f, pages = 1, guess = FALSE, output = "data.frame")
str(out2)
## List of 1
##  $ :'data.frame':       33 obs. of  13 variables:
##   ..$ X   : chr [1:33] "Mazda RX4 " "Mazda RX4 Wag " "Datsun 710 " "Hornet 4 Drive " ...
##   ..$ mpg : num [1:33] 21 21 22.8 21.4 18.7 18.1 14.3 24.4 22.8 19.2 ...
##   ..$ cyl : num [1:33] 6 6 4 6 8 6 8 4 4 6 ...
##   ..$ X.1 : int [1:33] NA NA NA NA NA NA NA NA NA NA ...
##   ..$ disp: num [1:33] 160 160 108 258 360 ...
##   ..$ hp  : num [1:33] 110 110 93 110 175 105 245 62 95 123 ...
##   ..$ drat: num [1:33] 3.9 3.9 3.85 3.08 3.15 2.76 3.21 3.69 3.92 3.92 ...
##   ..$ wt  : num [1:33] 2.62 2.88 2.32 3.21 3.44 ...
##   ..$ qsec: num [1:33] 16.5 17 18.6 19.4 17 ...
##   ..$ vs  : num [1:33] 0 0 1 1 0 1 0 1 1 1 ...
##   ..$ am  : num [1:33] 1 1 1 0 0 0 0 0 0 0 ...
##   ..$ gear: num [1:33] 4 4 4 3 3 3 3 4 4 4 ...
##   ..$ carb: int [1:33] 4 4 1 1 2 1 4 2 2 4 ...
```

It is also possible to manually specify smaller areas within pages to look for tables using the `area` and `columns` arguments to `extract_tables()`. This facilitates extraction from smaller portions of a page, such as when a table is embeded in a larger section of text or graphics.

Another function, `extract_areas()` implements this through an interactive style in which each page of the PDF is loaded as an R graphic and the user can use their mouse to specify upper-left and lower-right bounds of an area. Those areas are then extracted auto-magically (and the return value is the same as for `extract_tables()`). Here's a shot of it in action:

![extract_areas()](https://i.imgur.com/USTyQl7.gif)

`locate_areas()` handles the area identification process without performing the extraction, which may be useful as a debugger.

`extract_text()` simply returns text, possibly separately for each (specified) page:

```{r eval = FALSE}
out3 <- extract_text(f, page = 3)
cat(out3, sep = "\n")
## len supp dose
## 4.2 VC 0.5
## 11.5 VC 0.5
## 7.3 VC 0.5
## 5.8 VC 0.5
## 6.4 VC 0.5
## 10.0 VC 0.5
## 11.2 VC 0.5
## 11.2 VC 0.5
## 5.2 VC 0.5
## 7.0 VC 0.5
## 16.5 VC 1.0
## 16.5 VC 1.0
## 15.2 VC 1.0
## 17.3 VC 1.0
## 22.5 VC 1.0
## 3
```

Note that for large PDF files, it is possible to run up against Java memory constraints, leading to a `java.lang.OutOfMemoryError: Java heap space` error message. Memory can be increased using `options(java.parameters = "-Xmx16000m")` set to some reasonable amount of memory.

Some other utility functions are also provided (and made possible by the Java [Apache PDFBox library](https://pdfbox.apache.org/)):

 - `extract_text()` converts the text of an entire file or specified pages into an R character vector.
 - `split_pdf()` and `merge_pdfs()` split and merge PDF documents, respectively.
 - `extract_metadata()` extracts PDF metadata as a list.
 - `get_n_pages()` determines the number of pages in a document.
 - `get_page_dims()` determines the width and height of each page in pt (the unit used by `area` and `columns` arguments).
 - `make_thumbnails()` converts specified pages of a PDF file to image files.

### Installing Java on Windows with Chocolatey

In command prompt, install Chocolately if you don't already have it:

```
@powershell -NoProfile -ExecutionPolicy Bypass -Command "iex ((new-object net.webclient).DownloadString('https://chocolatey.org/install.ps1'))" && SET PATH=%PATH%;%ALLUSERSPROFILE%\chocolatey\bin
```

Then, install java using Chocolately's `choco install` command:

```
choco install jdk7 -y
```

You may also need to then set the `JAVA_HOME` environment variable to the path to your Java installation (e.g., `C:\Program Files\Java\jdk1.8.0_92`). This can be done:

 1. within R using `Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jdk1.8.0_92")` (note slashes), or
 2. from command prompt using the `setx` command: `setx JAVA_HOME C:\Program Files\Java\jdk1.8.0_92`, or
 3. from PowerShell, using the .NET framework: `[Environment]::SetEnvironmentVariable("JAVA_HOME", "C:\Program Files\Java\jdk1.8.0_92", "User")`, or
 4. from the Start Menu, via `Control Panel » System » Advanced » Environment Variables` ([instructions here](http://superuser.com/a/284351/221772)).

You should now be able to safely open R, and use rJava and tabulizer. Note, however, that some users report that rather than setting this variable, they instead need to delete it (e.g., with `Sys.setenv(JAVA_HOME = "")`), so if the above instructions fail, that is the next step in troubleshooting.

### Troubleshooting

Some notes for troubleshooting common installation problems:

 - On Mac OS, you may need to install [a particular version of Java](https://support.apple.com/kb/DL1572?locale=en_US) prior to attempting to install tabulizer.
 - On a Unix-like, you need to ensure that R has been installed with Java support. This can often be fixed by running `R CMD javareconf` on the command line (possibly with `sudo`, etc. depending on your system setup).
 - On Windows, make sure you have permission to write to and install packages to your R directory before trying to install the package. This can be changed from "Properties" on the right-click context menu. Alternatively, you can ensure write permission by choosing "Run as administrator" when launching R (again, from the right-click context menu).

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/tabulizer/issues).
* License: MIT
* Get citation information for `tabulizer` in R doing `citation(package = 'tabulizer')`

[![rofooter](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "Introduction to tabulizer"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to tabulizer}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

**tabulizer** provides R bindings to the [Tabula java library](https://github.com/tabulapdf/tabula-java/), which can be used to computationally extract tables from PDF documents. The main function `extract_tables()` mimics the command-line behavior of the Tabula, by extracting all tables from a PDF file and, by default, returns those tables as a list of character matrices in R.

```{r}
library("tabulizer")
f <- system.file("examples", "data.pdf", package = "tabulizer")

# extract table from first page of example PDF
tab <- extract_tables(f, pages = 1)
head(tab[[1]])
```

The `pages` argument allows you to select which pages to attempt to extract tables from. By default, Tabula (and thus tabulizer) checks every page for tables using a detection algorithm and returns all of them. `pages` can be an integer vector of any length; pages are indexed from 1.

It is possible to specify a remote file, which will be copied to R's temporary directory before processing:

```{r}
f2 <- "https://github.com/leeper/tabulizer/raw/master/inst/examples/data.pdf"
extract_tables(f2, pages = 2)
```

## Changing the Method of Extraction

The default method used by `extract_tables()` mimics the behaviour of Tabula.
For each page the algorithm decides whether it contains one consistent table
and then extracts it by using spreadsheet-tailored algorithm
`method = "lattice"`. The correct recognition of a table depends on whether the
page contains a table grid. If it doesn't and the table is a matrix of cells
with values without borders, it might not be able to recognise it. This also
happens when multiple tables with different number of columns are present on the
same page. In those cases another, more general, algorithm `method = "stream"`
is used, which relies on the distances between text characters on the page.

```{r}
# Extract tables by deciding for each page individually
extract_tables(f2, method = "decide")
```

It is possible to specify the preferred algorithm which might be a better option
for more difficult cases.

```{r}
# Extract tables by using "lattice" method
extract_tables(f2, pages = 2, method = "lattice")
```

```{r}
# Extract tables by using "stream" method
extract_tables(f2, pages = 2, method = "stream")
```

## Modifying the Return Value ##

By default, `extract_tables()` returns a list of character matrices. This is because many tables might be malformed or irregular and thus not be easily coerced to an R data.frame. This can easily be changed by specifying the `output` argument:

```{r}
# attempt to coerce tables to data.frames
extract_tables(f, pages = 2, output = "data.frame")
```

Tabula itself implements three "writer" methods that write extracted tables to disk as CSV, TSV, or JSON files. These can be specified by `output = "csv"`, `output = "tsv"`, and `output = "json"`, respectively. For CSV and TSV, one file is written to disk for each table and R session's temporary directory `tempdir()` is used by default (alternatively, the directory can be specified through `output` argument). For JSON, one file is written containing information about all tables. For these methods, `extract_tables()` returns a path to the directory containing the output files.

```{r}
# extract tables to CSVs
extract_tables(f, output = "csv")
```

If none of the standard methods works well, you can specify `output = "asis"` to return an rJava "jobjRef" object, which is a pointer to a Java ArrayList of Tabula Table objects. Working with that object might be quite awkward as it requires knowledge of Java and Tabula's internals, but might be useful to advanced users for debugging purposes.

## Extracting Areas ##

By default, tabulizer uses Tabula's table detection algorithm to automatically identify tables within each page of a PDF. This automatic detection can be toggled off by setting `guess = FALSE` and specifying an "area" within each PDF page to extract the table from. Here is a comparison of the default settings, versus extracting from two alternative areas within a page:

```{r}
str(extract_tables(f, pages = 2, guess = TRUE, output = "data.frame"))
str(extract_tables(f, pages = 2, area = list(c(126, 149, 212, 462)), guess = FALSE, output = "data.frame"))
str(extract_tables(f, pages = 2, area = list(c(126, 284, 174, 417)), guess = FALSE, output = "data.frame"))
```

The `area` argument should be a list either of length 1 (to use the same area for each specified page) or equal to the number of pages specified. This also means that you can extract multiple areas from one page, but specifying the page twice and indicating the two areas separately:

```{r}
a2 <- list(c(126, 149, 212, 462), c(126, 284, 174, 417))
str(extract_tables(f, pages = c(2,2), area = a2, guess = FALSE, output = "data.frame"))
```

## Interactive Table Extraction ##

In addition to the programmatic extraction offered by `extract_tables()`, it is also possible to work interactively with PDFs via the `extract_areas()` function. This function triggers a process by which each (specified) page of a PDF is converted to a PNG image file and then loaded as an R graphic. From there, you can use your mouse to specify upper-left and lower-right bounds of an area on each page. Pages are cycled through automatically and, after selecting areas for each page, those areas are extracted auto-magically (and the return value is the same as for `extract_tables()`). Here's a shot of it in action:

[![extract_areas()](http://i.imgur.com/USTyQl7.gif)](http://i.imgur.com/USTyQl7.gif)

`locate_areas()` handles the area identification process without performing the extraction, which may be useful as a debugger, or simply to define areas to be used in a programmatic extraction.

## Miscellaneous Functionality ##

Tabula is built on top of the [Java PDFBox library](https://pdfbox.apache.org/)), which provides low-level functionality for working with PDFs. A few of these tools are exposed through tabulizer, as they might be useful for debugging or generally for working with PDFs. These functions include:


 - `extract_text()` converts the text of an entire file or specified pages into an R character vector.
 - `split_pdf()` and `merge_pdfs()` split and merge PDF documents, respectively.
 - `extract_metadata()` extracts PDF metadata as a list.
 - `get_n_pages()` determines the number of pages in a document.
 - `get_page_dims()` determines the width and height of each page in pt (the unit used by `area` and `columns` arguments).
 - `make_thumbnails()` converts specified pages of a PDF file to image files.

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_metadata.R
\name{extract_metadata}
\alias{extract_metadata}
\title{extract_metadata}
\usage{
extract_metadata(file, password = NULL, copy = FALSE)
}
\arguments{
\item{file}{A character string specifying the path or URL to a PDF file.}

\item{password}{Optionally, a character string containing a user password to access a secured PDF.}

\item{copy}{Specifies whether the original local file(s) should be copied to
\code{tempdir()} before processing. \code{FALSE} by default. The argument is
ignored if \code{file} is URL.}
}
\value{
A list.
}
\description{
Extract metadata from a file
}
\details{
This function extracts metadata from a PDF
}
\examples{
\dontrun{
# simple demo file
f <- system.file("examples", "data.pdf", package = "tabulizer")

extract_metadata(f)
}
}
\seealso{
\code{\link{extract_tables}}, \code{\link{extract_areas}}, \code{\link{extract_text}}, \code{\link{split_pdf}}
}
\author{
Thomas J. Leeper <thosjleeper@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_page_dims.R
\name{get_page_dims}
\alias{get_page_dims}
\alias{get_n_pages}
\title{Page length and dimensions}
\usage{
get_page_dims(file, doc, pages = NULL, password = NULL, copy = FALSE)

get_n_pages(file, doc, password = NULL, copy = FALSE)
}
\arguments{
\item{file}{A character string specifying the path or URL to a PDF file.}

\item{doc}{Optionally,, in lieu of \code{file}, an rJava reference to a PDDocument Java object.}

\item{pages}{An optional integer vector specifying pages to extract from.}

\item{password}{Optionally, a character string containing a user password to access a secured PDF.}

\item{copy}{Specifies whether the original local file(s) should be copied to
\code{tempdir()} before processing. \code{FALSE} by default. The argument is
ignored if \code{file} is URL.}
}
\value{
For \code{get_n_pages}, an integer. For \code{get_page_dims}, a list of two-element numeric vectors specifying the width and height of each page, respectively.
}
\description{
Get Page Length and Dimensions
}
\details{
\code{get_n_pages} returns the page length of a PDF document. \code{get_page_dims} extracts the dimensions of specified pages in a PDF document. This can be useful for figuring out how to specify the \code{area} argument in \code{\link{extract_tables}}
}
\examples{
\dontrun{
# simple demo file
f <- system.file("examples", "data.pdf", package = "tabulizer")

get_n_pages(file = f)
get_page_dims(f)
}
}
\references{
\href{http://tabula.technology/}{Tabula}
}
\seealso{
\code{\link{extract_tables}}, \code{\link{extract_text}}, \code{\link{make_thumbnails}}
}
\author{
Thomas J. Leeper <thosjleeper@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/locate_area.R
\name{locate_areas}
\alias{locate_areas}
\alias{extract_areas}
\title{extract_areas}
\usage{
locate_areas(file, pages = NULL, resolution = 60L, widget = c("shiny",
  "native", "reduced"), copy = FALSE)

extract_areas(file, pages = NULL, guess = FALSE, copy = FALSE, ...)
}
\arguments{
\item{file}{A character string specifying the path to a PDF file. This can also be a URL, in which case the file will be downloaded to the R temporary directory using \code{download.file}.}

\item{pages}{An optional integer vector specifying pages to extract from. To extract multiple tables from a given page, repeat the page number (e.g., \code{c(1,2,2,3)}).}

\item{resolution}{An integer specifying the resolution of the PNG images conversions. A low resolution is used by default to speed image loading.}

\item{widget}{A one-element character vector specifying the type of \dQuote{widget} to use for locating the areas. The default (\dQuote{shiny}) is a shiny widget. The alternatives are a widget based on the native R graphics device (\dQuote{native}, where available), or a very reduced functionality model (\dQuote{reduced}).}

\item{copy}{Specifies whether the original local file(s) should be copied to
\code{tempdir()} before processing. \code{FALSE} by default. The argument is
ignored if \code{file} is URL.}

\item{guess}{See \code{\link{extract_tables}} (note the different default value).}

\item{\dots}{Other arguments passed to \code{\link{extract_tables}}.}
}
\value{
For \code{extract_areas}, see \code{\link{extract_tables}}. For \code{locate_areas}, a list of four-element numeric vectors (top,left,bottom,right), one per page of the file.
}
\description{
Interactively identify areas and extract
}
\details{
\code{extract_areas} is an interactive mode for \code{\link{extract_tables}} allowing the user to specify areas of each PDF page in a file that they would like extracted. When used, each page is rendered to a PNG file and displayed in an R graphics window sequentially, pausing on each page to call \code{\link[graphics]{locator}} so the user can click and highlight an area to extract.

The exact behaviour is a somewhat platform-dependent, and depends on the value of \code{widget} (and further, whether you are working in RStudio or the R console). In RStudio (where \code{widget = "shiny"}), a Shiny gadget is provided which allows the user to click and drag to select areas on each page of a file, clicking \dQuote{Done} on each page to advance through them. It is not possible to return to previous pages. In the R console, a Shiny app will be launched in a web browser.

For other values of \code{widget}, functionality is provided through the graphics device. If graphics events are supported, then it is possibly to interactively highlight a page region, make changes to that region, and navigate through the pages of the document while retaining the area highlighted on each page. If graphics events are not supported, then some of this functionality is not available (see below).

In \emph{full functionality mode} (\code{widget = "native"}), areas are input in a native graphics device. For each page, the first mouse click on a page initializes a highlighting rectangle; the second click confirms it. If unsatisfied with the selection, the process can be repeated. The window also responds to keystrokes. \kbd{PgDn}, \kbd{Right}, and \kbd{Down} advance to the next page image, while \kbd{PgUp}, \kbd{Left}, and \kbd{Up} return to the previous page image. \kbd{Home} returns to the first page image and \kbd{End} advances to the final page image. \kbd{Q} quits the interactive mode and proceeds with extraction. When navigating between pages, any selected areas will be displayed and can be edited. \kbd{Delete} removes a highlighted area from a page (and then displays it again). (This mode may not work correctly from within RStudio.)

In \emph{reduced functionality mode} (where \code{widget = "reduced"} or on platforms where graphics events are unavailable), the interface requires users to indicate the upper-left and lower-right (or upper-right and lower-left) corners of an area on each page, this area will be briefly confirmed with a highlighted rectangle and the next page will be displayed. Dynamic page navigation and area editing are not possible.

In any of these modes, after the areas are selected, \code{extract_areas} passes these user-defined areas to \code{\link{extract_tables}}. \code{locate_areas} implements the interactive component only, without actually extracting; this might be useful for interactive work that needs some modification before executing \code{extract_tables} computationally.
}
\examples{
\dontrun{
# simple demo file
f <- system.file("examples", "data.pdf", package = "tabulizer")

# locate areas only, using Shiny app
locate_areas(f)

# locate areas only, using native graphics device
locate_areas(f, widget = "shiny")

# locate areas and extract
extract_areas(f)
}
}
\seealso{
\code{\link{extract_tables}}, \code{\link{make_thumbnails}}, , \code{\link{get_page_dims}}
}
\author{
Thomas J. Leeper <thosjleeper@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.R
\docType{package}
\name{tabulizer-package}
\alias{tabulizer-package}
\alias{tabulizer}
\title{tabulizer}
\description{
Bindings for \dQuote{Tabula} PDF Table Extractor Library
}
\details{
Tabula is a Java library designed to computationally extract tables from PDF documents. tabulizer provides a thin R package with bindings to the library. It presently offers two principal functions: \code{\link{extract_tables}}, which mimics the command line functionality of Tabula, and \code{\link{extract_areas}} which provides an interactive interface to the former.
}
\references{
\href{http://tabula.technology/}{tabula}
}
\seealso{
\code{\link{extract_tables}}, \code{\link{extract_areas}}
}
\author{
Thomas J. Leeper <thosjleeper@gmail.com>
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logging.R
\name{stop_logging}
\alias{stop_logging}
\title{rJava logging}
\usage{
stop_logging()
}
\value{
\code{NULL}, invisibly.
}
\description{
Toggle verbose rJava logging
}
\details{
This function turns off the somewhat verbose rJava logging, most of which is uninformative. It is called automatically when tabulizer is attached via \code{library()}, \code{require}, etc. To keep logging on, load the package namespace using \code{requireNamespace("tabulizer")} and reference functions in using fully qualified references (e.g., \code{tabulizer::extract_tables()}.
}
\note{
This resets a global Java setting and may affect logging of other rJava operations, requiring a restart of R.
}
\examples{
\dontrun{
stop_logging()
}
}
\author{
Thomas J. Leeper <thosjleeper@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_thumbnails.R
\name{make_thumbnails}
\alias{make_thumbnails}
\title{make_thumbnails}
\usage{
make_thumbnails(file, outdir = NULL, pages = NULL, format = c("png",
  "jpeg", "bmp", "gif"), resolution = 72, password = NULL, copy = FALSE)
}
\arguments{
\item{file}{A character string specifying the path or URL to a PDF file.}

\item{outdir}{An optional character string specifying a directory into which
to split the resulting files. If \code{NULL}, the \code{outdir} is
\code{tempdir()}. If \code{file} is a URL, both file and thumbnails are
stored in the R session's temporary directory.}

\item{pages}{An optional integer vector specifying pages to extract from.}

\item{format}{A character string specifying an image file format.}

\item{resolution}{A numeric value specifying the image resolution in DPI.}

\item{password}{Optionally, a character string containing a user password
to access a secured PDF.}

\item{copy}{Specifies whether the original local file(s) should be copied to
\code{tempdir()} before processing. \code{FALSE} by default. The argument is
ignored if \code{file} is URL.}
}
\value{
A character vector of file paths.
}
\description{
Convert Pages to Image Thumbnails
}
\details{
This function save each (specified) page of a document as an image
with 720 dpi resolution. Images are saved in the same directory as the
original file, with file names specified by the original file name,
a page number, and the corresponding file format extension.
}
\note{
This may generate Java \dQuote{INFO} messages in the console,
which can be safely ignored.
}
\examples{
\dontrun{
# simple demo file
f <- system.file("examples", "data.pdf", package = "tabulizer")

make_thumbnails(f)
}
}
\references{
\href{http://tabula.technology/}{Tabula}
}
\seealso{
\code{\link{extract_tables}}, \code{\link{extract_text}},
\code{\link{make_thumbnails}}
}
\author{
Thomas J. Leeper <thosjleeper@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_merge.R
\name{split_pdf}
\alias{split_pdf}
\alias{merge_pdfs}
\title{Split and merge PDFs}
\usage{
split_pdf(file, outdir = NULL, password = NULL, copy = FALSE)

merge_pdfs(file, outfile, copy = FALSE)
}
\arguments{
\item{file}{For \code{merge_pdfs}, a character vector specifying the path to
one or more \emph{local} PDF files. For \code{split_pdf}, a character string
specifying the path or URL to a PDF file.}

\item{outdir}{For \code{split_pdf}, an optional character string specifying a
directory into which to split the resulting files. If \code{NULL}, the
\code{outdir} is \code{tempdir()}. If \code{file} is a URL, both the original
file and separate pages are stored in the R session's temporary directory.}

\item{password}{Optionally, a character string containing a user password to
access a secured PDF. Currently, encrypted PDFs cannot be merged with
\code{merge_pdfs}.}

\item{copy}{Specifies whether the original local file(s) should be copied to
\code{tempdir()} before processing. \code{FALSE} by default. The argument is
ignored if \code{file} is URL.}

\item{outfile}{For \code{merge_pdfs}, a character string specifying the path
to the PDF file to create from the merged documents.}
}
\value{
For \code{split_pdfs}, a character vector specifying the output file
names, which are patterned after the value of \code{file}. For
\code{merge_pdfs}, the value of \code{outfile}.
}
\description{
Split PDF into separate pages or merge multiple PDFs into one.
}
\details{
\code{\link{split_pdf}} splits the file listed in \code{file} into
separate one-page doucments. \code{\link{merge_pdfs}} creates a single PDF
document from multiple separate PDF files.
}
\examples{
\dontrun{
# simple demo file
f <- system.file("examples", "data.pdf", package = "tabulizer")
get_n_pages(file = f)

# split PDF by page
sf <- split_pdf(f)

# merge pdf
mf <- file.path(tempdir(), "merged.pdf")
merge_pdfs(sf, mf)
get_n_pages(mf)
}
}
\seealso{
\code{\link{extract_areas}}, \code{\link{get_page_dims}},
\code{\link{make_thumbnails}}
}
\author{
Thomas J. Leeper <thosjleeper@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_text.R
\name{extract_text}
\alias{extract_text}
\title{extract_text}
\usage{
extract_text(file, pages = NULL, area = NULL, password = NULL,
  encoding = NULL, copy = FALSE)
}
\arguments{
\item{file}{A character string specifying the path or URL to a PDF file.}

\item{pages}{An optional integer vector specifying pages to extract from.}

\item{area}{An optional list, of length equal to the number of pages specified, where each entry contains a four-element numeric vector of coordinates (top,left,bottom,right) containing the table for the corresponding page. As a convenience, a list of length 1 can be used to extract the same area from all (specified) pages.}

\item{password}{Optionally, a character string containing a user password to access a secured PDF.}

\item{encoding}{Optionally, a character string specifying an encoding for the text, to be passed to the assignment method of \code{\link[base]{Encoding}}.}

\item{copy}{Specifies whether the original local file(s) should be copied to
\code{tempdir()} before processing. \code{FALSE} by default. The argument is
ignored if \code{file} is URL.}
}
\value{
If \code{pages = NULL} (the default), a length 1 character vector, otherwise a vector of length \code{length(pages)}.
}
\description{
Extract text from a file
}
\details{
This function converts the contents of a PDF file into a single unstructured character string.
}
\examples{
\dontrun{
# simple demo file
f <- system.file("examples", "text.pdf", package = "tabulizer")

# extract all text
extract_text(f)

# extract all text from page 1 only
extract_text(f, pages = 1)

# extract text from selected area only
extract_text(f, area = list(c(209.4, 140.5, 304.2, 500.8)))

}
}
\seealso{
\code{\link{extract_tables}}, \code{\link{extract_areas}}, \code{\link{split_pdf}}
}
\author{
Thomas J. Leeper <thosjleeper@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_tables.R
\name{extract_tables}
\alias{extract_tables}
\title{extract_tables}
\usage{
extract_tables(file, pages = NULL, area = NULL, columns = NULL,
  guess = TRUE, method = c("decide", "lattice", "stream"),
  output = c("matrix", "data.frame", "character", "asis", "csv", "tsv",
  "json"), outdir = NULL, password = NULL, encoding = NULL,
  copy = FALSE, ...)
}
\arguments{
\item{file}{A character string specifying the path or URL to a PDF file.}

\item{pages}{An optional integer vector specifying pages to extract from.}

\item{area}{An optional list, of length equal to the number of pages specified, where each entry contains a four-element numeric vector of coordinates (top,left,bottom,right) containing the table for the corresponding page. As a convenience, a list of length 1 can be used to extract the same area from all (specified) pages. Only specify \code{area} xor \code{columns}.}

\item{columns}{An optional list, of length equal to the number of pages specified, where each entry contains a numeric vector of horizontal (x) coordinates separating columns of data for the corresponding page. As a convenience, a list of length 1 can be used to specify the same columns for all (specified) pages. Only specify \code{area} xor \code{columns}.}

\item{guess}{A logical indicating whether to guess the locations of tables on each page. If \code{FALSE}, \code{area} or \code{columns} must be specified; if \code{TRUE}, columns is ignored.}

\item{method}{A string identifying the prefered method of table extraction.
\itemize{
  \item \code{method = "decide"} (default) automatically decide (for each page) whether spreadsheet-like formatting is present and "lattice" is appropriate
  \item \code{method = "lattice"} use Tabula's spreadsheet extraction algorithm
  \item \code{method = "stream"} use Tabula's basic extraction algorithm
}}

\item{output}{A function to coerce the Java response object (a Java ArrayList of Tabula Tables) to some output format. The default method, \dQuote{matrices}, returns a list of character matrices. See Details for other options.}

\item{outdir}{Output directory for files if \code{output} is set to
\code{"csv"}, \code{"tsv"} or \code{"json"}, ignored otherwise. If equals
\code{NULL} (default), uses R sessions temporary directory \code{tempdir()}.}

\item{password}{Optionally, a character string containing a user password to access a secured PDF.}

\item{encoding}{Optionally, a character string specifying an encoding for the text, to be passed to the assignment method of \code{\link[base]{Encoding}}.}

\item{copy}{Specifies whether the original local file(s) should be copied to
\code{tempdir()} before processing. \code{FALSE} by default. The argument is
ignored if \code{file} is URL.}

\item{\dots}{These are additional arguments passed to the internal functions dispatched by \code{method}.}
}
\value{
By default, a list of character matrices. This can be changed by specifying an alternative value of \code{method} (see Details).
}
\description{
Extract tables from a file
}
\details{
This function mimics the behavior of the Tabula command line utility. It returns a list of R character matrices containing tables extracted from a file by default. This response behavior can be changed by using the following options.
\itemize{
  \item \code{output = "character"} returns a list of single-element character vectors, where each vector is a tab-delimited, line-separate string of concatenated table cells.
  \item \code{output = "data.frame"} attempts to coerce the structure returned by \code{method = "character"} into a list of data.frames and returns character strings where this fails.
  \item \code{output = "csv"} writes the tables to comma-separated (CSV) files using Tabula's CSVWriter method in the same directory as the original PDF. \code{method = "tsv"} does the same but with tab-separated (TSV) files using Tabula's TSVWriter and \code{method = "json"} does the same using Tabula's JSONWriter method. Any of these three methods return the path to the directory containing the extract table files. 
  \item \code{output = "asis"} returns the Java object reference, which can be useful for debugging or for writing a custom parser.
}
\code{\link{extract_areas}} implements this functionality in an interactive mode allowing the user to specify extraction areas for each page.
}
\examples{
\dontrun{
# simple demo file
f <- system.file("examples", "data.pdf", package = "tabulizer")

# extract all tables
extract_tables(f)

# extract tables from only second page
extract_tables(f, pages = 2)

# extract areas from a page
## full table
extract_tables(f, pages = 2, area = list(c(126, 149, 212, 462)))
## part of the table
extract_tables(f, pages = 2, area = list(c(126, 284, 174, 417)))

# return data.frames
extract_tables(f, pages = 2, output = "data.frame")
}
}
\references{
\href{http://tabula.technology/}{Tabula}
}
\seealso{
\code{\link{extract_areas}}, \code{\link{get_page_dims}}, \code{\link{make_thumbnails}}, \code{\link{split_pdf}}
}
\author{
Thomas J. Leeper <thosjleeper@gmail.com>, Tom Paskhalis <tpaskhalis@gmail.com>
}
