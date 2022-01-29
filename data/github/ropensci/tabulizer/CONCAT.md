
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

