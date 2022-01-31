
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epubr <img src="man/figures/logo.png" style="margin-left:10px;margin-bottom:5px;" width="120" align="right">

**Author:** [Matthew Leonawicz](https://github.com/leonawicz)
<a href="https://orcid.org/0000-0001-9452-2771" target="orcid.widget">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
<br/> **License:** [MIT](https://opensource.org/licenses/MIT)<br/>

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ropensci/epubr?branch=master&svg=true)](https://ci.appveyor.com/project/leonawicz/epubr)
[![codecov](https://codecov.io/gh/ropensci/epubr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/epubr)
[![](https://badges.ropensci.org/222_status.svg)](https://github.com/ropensci/software-review/issues/222)

[![CRAN
status](http://www.r-pkg.org/badges/version/epubr)](https://cran.r-project.org/package=epubr)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/epubr)](https://cran.r-project.org/package=epubr)
[![Github
Stars](https://img.shields.io/github/stars/ropensci/epubr.svg?style=social&label=Github)](https://github.com/ropensci/epubr)

## Read EPUB files in R

Read EPUB text and metadata.

The `epubr` package provides functions supporting the reading and
parsing of internal e-book content from EPUB files. E-book metadata and
text content are parsed separately and joined together in a tidy, nested
tibble data frame.

E-book formatting is not completely standardized across all literature.
It can be challenging to curate parsed e-book content across an
arbitrary collection of e-books perfectly and in completely general
form, to yield a singular, consistently formatted output. Many EPUB
files do not even contain all the same pieces of information in their
respective metadata.

EPUB file parsing functionality in this package is intended for
relatively general application to arbitrary EPUB e-books. However,
poorly formatted e-books or e-books with highly uncommon formatting may
not work with this package. There may even be cases where an EPUB file
has DRM or some other property that makes it impossible to read with
`epubr`.

Text is read ‘as is’ for the most part. The only nominal changes are
minor substitutions, for example curly quotes changed to straight
quotes. Substantive changes are expected to be performed subsequently by
the user as part of their text analysis. Additional text cleaning can be
performed at the user’s discretion, such as with functions from packages
like `tm` or `qdap`.

## Installation

Install `epubr` from CRAN with:

``` r
install.packages("epubr")
```

Install the development version from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/epubr")
```

## Example

Bram Stoker’s Dracula novel sourced from Project Gutenberg is a good
example of an EPUB file with unfortunate formatting. The first thing
that stands out is the naming convention using `item` followed by some
ordered digits does not differentiate sections like the book preamble
from the chapters. The numbering also starts in a weird place. But it is
actually worse than this. Notice that sections are not broken into
chapters; they can begin and end in the middle of chapters\!

These annoyances aside, the metadata and contents can still be read into
a convenient table. Text mining analyses can still be performed on the
overall book, if not so easily on individual chapters. See the [package
vignette](https://docs.ropensci.org/epubr/articles/epubr.html) for
examples on how to further improve the structure of an e-book with
formatting like this.

``` r
file <- system.file("dracula.epub", package = "epubr")
(x <- epub(file))
#> # A tibble: 1 x 9
#>   rights         identifier           creator   title  language subject                                                                       date                source                 data       
#>   <chr>          <chr>                <chr>     <chr>  <chr>    <chr>                                                                         <chr>               <chr>                  <list>     
#> 1 Public domain~ http://www.gutenber~ Bram Sto~ Dracu~ en       Horror tales|Epistolary fiction|Gothic fiction (Literary genre)|Vampires -- ~ 1995-10-01|2017-10~ http://www.gutenberg.~ <tibble[,4~

x$data[[1]]
#> # A tibble: 15 x 4
#>    section          text                                                                                                                                                                 nword nchar
#>    <chr>            <chr>                                                                                                                                                                <int> <int>
#>  1 item6            "The Project Gutenberg EBook of Dracula, by Bram StokerThis eBook is for the use of anyone anywhere at no cost and withalmost no restrictions whatsoever.  You may ~ 11446 60972
#>  2 item7            "But I am not in heart to describe beauty, for when I had seen the view I explored further; doors, doors, doors everywhere, and all locked and bolted. In no place ~ 13879 71798
#>  3 item8            "\" 'Lucy, you are an honest-hearted girl, I know. I should not be here speaking to you as I am now if I did not believe you clean grit, right through to the very ~ 12474 65522
#>  4 item9            "CHAPTER VIIIMINA MURRAY'S JOURNAL\nSame day, 11 o'clock p. m.-Oh, but I am tired! If it were not that I had made my diary a duty I should not open it to-night. We~ 12177 62724
#>  5 item10           "CHAPTER X\nLetter, Dr. Seward to Hon. Arthur Holmwood.\n\"6 September.\n\"My dear Art,-\n\"My news to-day is not so good. Lucy this morning had gone back a bit. T~ 12806 66678
#>  6 item11           "Once again we went through that ghastly operation. I have not the heart to go through with the details. Lucy had got a terrible shock and it told on her more than~ 12103 62949
#>  7 item12           "CHAPTER XIVMINA HARKER'S JOURNAL\n23 September.-Jonathan is better after a bad night. I am so glad that he has plenty of work to do, for that keeps his mind off t~ 12214 62234
#>  8 item13           "CHAPTER XVIDR. SEWARD'S DIARY-continued\nIT was just a quarter before twelve o'clock when we got into the churchyard over the low wall. The night was dark with oc~ 13990 72903
#>  9 item14           "\"Thus when we find the habitation of this man-that-was, we can confine him to his coffin and destroy him, if we obey what we know. But he is clever. I have asked~ 13356 69779
#> 10 item15           "\"I see,\" I said. \"You want big things that you can make your teeth meet in? How would you like to breakfast on elephant?\"\n\"What ridiculous nonsense you are ~ 12866 66921
#> 11 item16           "CHAPTER XXIIIDR. SEWARD'S DIARY\n3 October.-The time seemed terrible long whilst we were waiting for the coming of Godalming and Quincey Morris. The Professor tri~ 11928 61550
#> 12 item17           "CHAPTER XXVDR. SEWARD'S DIARY\n11 October, Evening.-Jonathan Harker has asked me to note this, as he says he is hardly equal to the task, and he wants an exact re~ 13119 68564
#> 13 item18           " \nLater.-Dr. Van Helsing has returned. He has got the carriage and horses; we are to have some dinner, and to start in an hour. The landlady is putting us up a h~  8435 43464
#> 14 item19           "End of the Project Gutenberg EBook of Dracula, by Bram Stoker*** END OF THIS PROJECT GUTENBERG EBOOK DRACULA ******** This file should be named 345-h.htm or 345-h~  2665 18541
#> 15 coverpage-wrapp~ ""                                                                                                                                                                       0     0
```

## Related packages

[tesseract](https://github.com/ropensci/tesseract) by @jeroen for more
direct control of the OCR process.

[pdftools](https://github.com/ropensci/pdftools) for extracting metadata
and text from PDF files (therefore more specific to PDF, and without a
Java dependency)

[tabulizer](https://github.com/ropensci/tabulizer) by @leeper and
@tpaskhalis, Bindings for Tabula PDF Table Extractor Library, to extract
tables, therefore not text, from PDF files.

[rtika](https://github.com/ropensci/rtika) by @goodmansasha for more
general text parsing.

[gutenbergr](https://github.com/ropensci/gutenbergr) by @dgrtwo for
searching and downloading public domain texts from Project Gutenberg.

-----

Please note that the `epubr` project is released with a [Contributor
Code of
Conduct](https://github.com/ropensci/epubr/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# epubr 0.6.3

* Allow multiple values for a metadata field, pipe-separated.
* Update documentation and tests.

# epubr 0.6.2

* Documentation updates.

# epubr 0.6.1

* Fix warnings caused by updates to other packages like `tidyr`.

# epubr 0.6.0

* Added `count_words` helper function.
* Improved word count accuracy in `epub`. Now also splitting words on new line characters rather than only on spaces. Now also ignoring vector elements in the split result that are most likely to not be words, such as stranded pieces of punctuation.
* Added `epub_recombine` for breaking apart and recombining text sections into new data frame rows using alternative breaks based on a regular expression pattern.
* Added `epub_sift` function for filtering out small text sections based on low word or character count. This function can also be used directly inside calls to `epub_recombine` through an argument list.
* Added `epub_reorder` for reordering a specified (by index) subset of text section data frame rows according to a text parsing function (several template functions are available for convenient use to address common cases).
* Refactored code to remove `purrr` dependency.
* Added unit tests.
* Updated function documentation, readme and vignette.

# epubr 0.5.0

* Added `epub_cat` function for pretty printing to console as a helpful way to quickly inspect the parsed text in a more easily readable format than looking at the quoted strings in the table entries. `epub_cat` can take an EPUB filename string (may be a vector) as its first argument or a data frame already returned by `epub`.
* Like `epub_cat`, `epub_head` accepts EPUB character filenames or now also a data frame already returned by `epub` based on those files. Because of this change, the first argument has been renamed from `file` to `x`.
* Added `encoding` argument to `epub` function, defaulting to UTF-8.
    * This helps significantly with reading EPUB archive files properly, e.g., providing ability to parse and substitute all the curly single and double quotes, apostrophes, various forms of hyphens and ellipses.
    * Previously, these were not substituted (e.g., replacing curly quotes with straight quotes), but attempting to do so would have failed anyway because they were not initially read correctly due to the lack of encoding specification.
    * Now non-standard characters are more likely to be read correctly, and those mentioned above are substituted with standard versions. If necessary, the encoding can be changed from UTF-8 via the new argument.
    * It appears that the EPUB format *requires* UTF encoding. Currently the only permissible option other than UTF-8 is UTF-16. This keeps things very simple and straightforward. Users should not encounter EPUB files in other encodings.
* Added unit tests and updated documentation.

# epubr 0.4.1

* Improved handling of errors and better messages.
* More robust handling of `title` field when missing, redundant or requiring remapping/renaming. All outputs of `epub` now include a `title` as well as `data` field, even if the e-book does not have a metadata field named `title`.
* Minor improvements to e-book section handling.
* Added `epub_head` function for previewing the opening text of each e-book section.
* Removed R version from Depends field of DESCRIPTION. Package Imports that necessitated a higher R version were previously removed.
* Minor fixes.
* Updated documentation, vignette, unit tests.

# epubr 0.4.0

* Enhanced function documentation details.
* Added `epub_meta` for strictly parsing EPUB metadata without reading the full file contents.
* When working with a vector of EPUB files, functions now clean up each unzipped archive temp directory with `unlink` immediately after use, rather than after all files are read into memory or by overwriting files in a single temp directory.
* Added initial introduction vignette content.
* Minor function refactors.
* Minor bug fixes.
* Added unit tests.

# epubr 0.3.0

* Refactor functions.
* Further reduce package dependencies.
* Update unit tests and documentation.

# epubr 0.2.0

* Refactor functions.
* Move contextual and e-book collection-specific functionality to other packages.
* Make any other remaining edge-case related options hidden arguments so that general usage of `epubr` functions is not too inflexible.
* Reduce package dependencies.
* Add basic unit tests.
* Add example public domain EPUB book for examples and testing.
* Update readme and documentation.

# epubr 0.1.0

* Added initial package scaffolding and function.
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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
## Test environments

* local Windows 10 install, R 4.0.5
* Windows 10 (AppVeyor), R 4.0.5
* Ubuntu 16.04 (Travis CI), R-devel, R-release, R-oldrel
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Update release

* Minor fixes.
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 3.6.1 (2019-07-05) |
|os       |Windows 10 x64               |
|system   |x86_64, mingw32              |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |English_United States.1252   |
|ctype    |English_United States.1252   |
|tz       |America/Denver               |
|date     |2019-11-27                   |

# Dependencies

|package    |old      |new      |<U+0394>  |
|:----------|:--------|:--------|:--|
|epubr      |0.6.0    |0.6.1    |*  |
|assertthat |0.2.1    |0.2.1    |   |
|backports  |1.1.5    |1.1.5    |   |
|BH         |1.69.0-1 |1.69.0-1 |   |
|cli        |1.1.0    |1.1.0    |   |
|crayon     |1.3.4    |1.3.4    |   |
|digest     |0.6.23   |0.6.23   |   |
|dplyr      |0.8.3    |0.8.3    |   |
|ellipsis   |0.3.0    |0.3.0    |   |
|fansi      |0.4.0    |0.4.0    |   |
|glue       |1.3.1    |1.3.1    |   |
|lifecycle  |0.1.0    |0.1.0    |   |
|magrittr   |1.5      |1.5      |   |
|pillar     |1.4.2    |1.4.2    |   |
|pkgconfig  |2.0.3    |2.0.3    |   |
|plogr      |0.2.0    |0.2.0    |   |
|purrr      |0.3.3    |0.3.3    |   |
|R6         |2.4.1    |2.4.1    |   |
|Rcpp       |1.0.3    |1.0.3    |   |
|rlang      |0.4.2    |0.4.2    |   |
|stringi    |1.4.3    |1.4.3    |   |
|tibble     |2.1.3    |2.1.3    |   |
|tidyr      |1.0.0    |1.0.0    |   |
|tidyselect |0.2.5    |0.2.5    |   |
|utf8       |1.1.4    |1.1.4    |   |
|vctrs      |0.2.0    |0.2.0    |   |
|xml2       |1.2.2    |1.2.2    |   |
|xslt       |1.3      |1.3      |   |
|zeallot    |0.1.0    |0.1.0    |   |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", fig.path = "man/figures/README-",
  message = FALSE, warning = FALSE, error = FALSE
)
library(epubr)
```

# epubr <img src="man/figures/logo.png" style="margin-left:10px;margin-bottom:5px;" width="120" align="right">
**Author:** [Matthew Leonawicz](https://github.com/leonawicz) <a href="https://orcid.org/0000-0001-9452-2771" target="orcid.widget">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
<br/>
**License:** [MIT](https://opensource.org/licenses/MIT)<br/>

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/epubr?branch=master&svg=true)](https://ci.appveyor.com/project/leonawicz/epubr)
[![codecov](https://codecov.io/gh/ropensci/epubr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/epubr)
[![](https://badges.ropensci.org/222_status.svg)](https://github.com/ropensci/software-review/issues/222)

[![CRAN status](http://www.r-pkg.org/badges/version/epubr)](https://cran.r-project.org/package=epubr)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/epubr)](https://cran.r-project.org/package=epubr)
[![Github Stars](https://img.shields.io/github/stars/ropensci/epubr.svg?style=social&label=Github)](https://github.com/ropensci/epubr)

## Read EPUB files in R

Read EPUB text and metadata. 

The `epubr` package provides functions supporting the reading and parsing of internal e-book content from EPUB files. E-book metadata and text content are parsed separately and joined together in a tidy, nested tibble data frame. 

E-book formatting is not completely standardized across all literature. It can be challenging to curate parsed e-book content across an arbitrary collection of e-books perfectly and in completely general form, to yield a singular, consistently formatted output. Many EPUB files do not even contain all the same pieces of information in their respective metadata.

EPUB file parsing functionality in this package is intended for relatively general application to arbitrary EPUB e-books. However, poorly formatted e-books or e-books with highly uncommon formatting may not work with this package.
There may even be cases where an EPUB file has DRM or some other property that makes it impossible to read with `epubr`.

Text is read 'as is' for the most part. The only nominal changes are minor substitutions, for example curly quotes changed to straight quotes. Substantive changes are expected to be performed subsequently by the user as part of their text analysis. Additional text cleaning can be performed at the user's discretion, such as with functions from packages like `tm` or `qdap`.

## Installation

Install `epubr` from CRAN with:

``` r
install.packages("epubr")
```

Install the development version from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/epubr")
```

## Example

Bram Stoker's Dracula novel sourced from Project Gutenberg is a good example of an EPUB file with unfortunate formatting.
The first thing that stands out is the naming convention using `item` followed by some ordered digits does not differentiate sections like the book preamble from the chapters.
The numbering also starts in a weird place. But it is actually worse than this. Notice that sections are not broken into chapters; they can begin and end in the middle of chapters!

These annoyances aside, the metadata and contents can still be read into a convenient table. Text mining analyses can still be performed on the overall book, if not so easily on individual chapters. See the [package vignette](https://docs.ropensci.org/epubr/articles/epubr.html) for examples on how to further improve the structure of an e-book with formatting like this.

```{r ex}
file <- system.file("dracula.epub", package = "epubr")
(x <- epub(file))

x$data[[1]]
```

## Related packages

[tesseract](https://github.com/ropensci/tesseract) by @jeroen for more direct control of the OCR process.

[pdftools](https://github.com/ropensci/pdftools) for extracting metadata and text from PDF files (therefore more specific to PDF, and without a Java dependency)

[tabulizer](https://github.com/ropensci/tabulizer) by @leeper and @tpaskhalis, Bindings for Tabula PDF Table Extractor Library, to extract tables, therefore not text, from PDF files.

[rtika](https://github.com/ropensci/rtika) by @goodmansasha for more general text parsing.

[gutenbergr](https://github.com/ropensci/gutenbergr) by @dgrtwo for searching and downloading public domain texts from Project Gutenberg.

---

Please note that the `epubr` project is released with a [Contributor Code of Conduct](https://github.com/ropensci/epubr/blob/master/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Introduction to epubr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to epubr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE, error = FALSE, tidy = TRUE
)
library(dplyr)
```

The `epubr` package provides functions supporting the reading and parsing of internal e-book content from EPUB files. E-book metadata and text content are parsed separately and joined together in a tidy, nested tibble data frame. 

E-book formatting is not completely standardized across all literature. It can be challenging to curate parsed e-book content across an arbitrary collection of e-books perfectly and in completely general form, to yield a singular, consistently formatted output. Many EPUB files do not even contain all the same pieces of information in their respective metadata.

EPUB file parsing functionality in this package is intended for relatively general application to arbitrary EPUB e-books. However, poorly formatted e-books or e-books with highly uncommon formatting may not work with this package.
There may even be cases where an EPUB file has DRM or some other property that makes it impossible to read with `epubr`.

Text is read 'as is' for the most part. The only nominal changes are minor substitutions, for example curly quotes changed to straight quotes. Substantive changes are expected to be performed subsequently by the user as part of their text analysis. Additional text cleaning can be performed at the user's discretion, such as with functions from packages like `tm` or `qdap`.

## Read EPUB files

Bram Stoker's Dracula novel sourced from Project Gutenberg is a good example of an EPUB file with unfortunate formatting.
The first thing that stands out is the naming convention using `item` followed by some ordered digits does not differentiate sections like the book preamble from the chapters.
The numbering also starts in a weird place. But it is actually worse than this. Notice that sections are not broken into chapters; they can begin and end in the middle of chapters!

These annoyances aside, the metadata and contents can still be read into a convenient table. Text mining analyses can still be performed on the overall book, if not so easily on individual chapters. See the section below on restructuring for examples of `epubr` functions that help get around these issues.

Here a single file is read with `epub()`. The output of the returned primary data frame and the book text data frame that is nested within its `data` column are shown.

```{r read}
library(epubr)
file <- system.file("dracula.epub", package = "epubr")
(x <- epub(file))

x$data[[1]]
```

The `file` argument may be a vector of EPUB files. There is one row for each book.

## EPUB metadata

The above examples jump right in, but it can be helpful to inspect file metadata before reading a large number of books into memory. Formatting may differ across books. It can be helpful to know what fields to expect, the degree of consistency, and what content you may want to drop during the file reading process. `epub_meta()` strictly parses file metadata and does not read the e-book text.

```{r metadata}
epub_meta(file)
```

This provides the big picture, though it will not reveal the internal breakdown of book section naming conventions that were seen in the first `epub()` example.

`file` can also be a vector for `epub_meta()`. Whenever `file` is a vector, the fields (columns) returned are the union of all fields detected across all EPUB files. Any books (rows) that do not have a field found in another book return `NA` for that row and column.

## Additional arguments

There are three optional arguments that can be provided to `epub()` to:

*    select fields, or columns of the primary data frame.
*    filter sections, or rows of the nested data frame.
*    attempt to detect which rows or sections in the nested data frame identify book chapters.

Unless you have a collection of well-formatted and similarly formatted EPUB files, these arguments may not be helpful and can be ignored, especially chapter detection.

### Select fields

Selecting fields is straightforward. All fields found are returned unless a vector of fields is provided.

```{r fields}
epub(file, fields = c("title", "creator", "file"))
```

Note that `file` was not a field identified in the metadata. This is a special case. Including `file` will include the `basename` of the input file. This is helpful when you want to retain file names and `source` is included in the metadata but may represent something else. Some fields like `data` and `title` are always returned and do not need to be specified in `fields`.

Also, if your e-book does not have a metadata field named `title`, you can pass an additional argument to `...` to map a different, known metadata field to `title`. E.g., `title = "BookTitle"`. The resulting table always has a `title` field, but in this case `title` would be populated with information from the `BookTitle` metadata field. If the default `title` field or any other field name passed to the additional `title` argument does not exist in the file metadata, the output `title` column falls back on filling in with the same unique file names obtained when requesting the `file` field.

### Drop sections

Filtering out unwanted sections, or rows of the nested data frame, uses a regular expression pattern. Matched rows are dropped. This is where knowing the naming conventions used in the e-books in `file`, or at least knowing they are satisfactorily consistent and predictable for a collection, helps with removing extraneous clutter.

One section that can be discarded is the cover. For many books it can be helpful to use a pattern like `"^(C|c)ov"` to drop any sections whose IDs begin with `Cov`, `cov`, and may be that abbreviation or the full word. For this book, `cov` suffices. The nested data frame has one less row than before.

```{r drop_sections}
epub(file, drop_sections = "cov")$data[[1]]
```

### Guess chapters

This e-book unfortunately does not have great formatting. For the sake of example, pretend that chapters are known to be sections beginning with `item` and followed by *two* digits, using the pattern `^item\\d\\d`. This does two things. It adds a new metadata column to the primary data frame called `nchap` giving the estimated number of chapters in the book. In the nested data frame containing the parsed e-book text, the `section` column is conditionally mutated to reflect a new, consistent chapter naming convention for the identified chapters and a logical `is_chapter` column is added.

```{r chapters}
x <- epub(file, drop_sections = "cov", chapter_pattern = "^item\\d\\d")
x

x$data[[1]]
```

This renaming of sections is helpful if the e-book content is split into meaningful sections, but the sections are not named in a similarly meaningful way. It is not helpful in cases where the sections are poorly or arbitrarily defined as in the above example where they begin with `item*` regardless of content. See the examples below on restructuring parsed content for an approach to dealing with this more problematic case.

Another problematic formatting that you may encounter is the random order of book sections as defined in the metadata, even if the sections themselves are well-defined. For example, there may be clear breaks between chapters, but the chapters may occur out of order such as `ch05, ch11, ch02, ...` and so on. Merely identifying which sections are chapters does not help in this case. They are still be labeled incorrectly since they are already out of order. Properly reordering sections is also addressed below.

Also note that not all books have chapters or something like them. Make sure an optional argument like `chapter_pattern` makes sense to use with a given e-book in the first place.

Ultimately, everything depends on the quality of the EPUB file. Some publishers are better than others. Formatting standards may also change over time.

## Restructure parsed content

When reading EPUB files it is ideal to be able to identify meaningful sections to retain via a regular expression pattern, as well as to drop extraneous sections in a similar manner. Using pattern matching as shown above is a convenient way to filter rows of the nested text content data frame.

For e-books with poor metadata formatting this is not always possible, or may be possible only after some other pre-processing. `epubr` provides other functions to assist in restructuring the text table. The Dracula EPUB file included in `epubr` is a good example to continue with here.

### Split and recombine into new sections

This book is internally broken into sections at arbitrary break points, hence why several sections begin in the middle of chapters, as seen above. Other chapters begin in the middle of sections. Use `epub_recombine()` along with a regular expression that can match the true section breaks. This function collapses the full text and then rebuilds the text table using new sections with proper break points. In the process it also recalculates the numbers of words and characters and relabels the sections with chapter notation.

Fortunately, a reliable pattern exists, which consists of `CHAPTER` in capital letters followed by a space and some Roman numerals. Recombine the text into a new object.

```{r recombine1}
pat <- "CHAPTER [IVX]+"
x2 <- epub_recombine(x, pat)
x2

x2$data[[1]]
```

But this is not quite as expected. There should be 27 chapters, not 54. What was not initially apparent was that the same pattern matching each chapter name also appears in the first section where every chapter is listed in the table of contents. The new section breaks were successful in keeping chapter text in single, unique sections, but there are now twice as many as needed. Unintentionally, the first 27 "chapters" represent the table of contents being split on each chapter ID. These should be removed.

An easy way to do this is with `epub_sift()`, which sifts, or filters out, small word- or character-count sections from the nested data frame. It's a simple sieve and you can control the size of the holes with `n`. You can choose `type = "word"` (default) or `type = "character"`. This is somewhat of a blunt instrument, but is useful in a circumstance like this one where it is clear it will work as desired.

```{r recombine2}
library(dplyr)
x2 <- epub_recombine(x, pat) %>% epub_sift(n = 200)
x2

x2$data[[1]]
```

This removes the unwanted rows, but one problem remains. Note that sifting the table sections in this case results in a need to re-apply `epub_recombine()` because the sections we removed had nevertheless offset the chapter indexing. Another call to `epub_recombine()` can be chained, but it may be more convenient to use the `sift` argument to `epub_recombine()`, which is applied recursively.

```{r recombine3}
#epub_recombine(x, pat) %>% epub_sift(n = 200) %>% epub_recombine(pat)
x2 <- epub_recombine(x, pat, sift = list(n = 200))
x2

x2$data[[1]]
```

### Reorder sections based on pattern in text

Some poorly formatted e-books have their internal sections occur in an arbitrary order. This can be frustrating to work with when doing text analysis on each section and where order matters. Just like recombining into new sections based on a pattern, sections that are out of order can be reordered based on a pattern. This requires a bit more work. In this case the user must provide a function that will map something in the matched pattern to an integer representing the desired row index.

Continue with the Dracula example, but with one difference. Even though the sections were originally broken at arbitrary points, they were in chronological order. To demonstrate the utility of `epub_reorder()`, first randomize the rows so that chronological order can be recovered.

```{r reorder1}
set.seed(1)
x2$data[[1]] <- sample_frac(x2$data[[1]]) # randomize rows for example
x2$data[[1]]
```

It is clear above that sections are now out of order. It is common enough to load poorly formatted EPUB files and yield this type of result. If all you care about is the text in its entirely, this does not matter, but if your analysis involves trends over the course of a book, this is problematic.

For this book, you need a function that will convert an annoying Roman numeral to an integer. You already have the pattern for finding the relevant information in each text section. You only need to tweak it for proper substitution. Here is an example:

```{r reorder2}
f <- function(x, pattern) as.numeric(as.roman(gsub(pattern, "\\1", x)))
```

This function is passed to `epub_reorder()`. It takes and returns scalars. It must take two arguments: the first is a text string. The second is the regular expression. It must return a single number representing the index of that row. For example, if the pattern matches `CHAPTER IV`, the function should return a `4`.

`epub_reorder()` takes care of the rest. It applies your function to every row in the the nested data frame and then reorders the rows based on the full set of indices. Note that it also repeats this for every row (book) in the primary data frame, i.e., for every nested table. This means that the same function will be applied to every book. Therefore, you should only use this in bulk on a collection of e-books if you know the pattern does not change and the function will work correctly in each case.

The pattern has changed slightly. Parentheses are used to retain the important part of the matched pattern, the Roman numeral. The function `f` here substitutes the entire string (because now it begins with `^` and ends with `.*`) with only the part stored in parentheses (In `f`, this is the `\\1` substitution). `epub_reorder()` applies this to all rows in the nested data frame:

```{r reorder3}
x2 <- epub_reorder(x2, f, "^CHAPTER ([IVX]+).*")
x2$data[[1]]
```

It is important that this is done on a nested data frame that has already been cleaned to the point of not containing extraneous rows that cannot be matched by the desired pattern. If they cannot be matched, then it is unknown where those rows should be placed relative to the others.

If sections are both out of order and use arbitrary break points, it would be necessary to reorder them before you split and recombine. If you split and recombine first, this would yield new sections that contain text from different parts of the e-book. However, the two are not likely to occur together; in fact it may be impossible for an EPUB file to be structured this way. In developing `epubr`, no such examples have been encountered. In any event, reordering out of order sections essentially requires a human-identifiable pattern near the beginning of each section text string, so it does not make sense to perform this operation unless the sections have meaningful break points. 

## Unzip EPUB file

Separate from using `epub_meta()` and `epub()`, you can call `epub_unzip()` directly if all you want to do is extract the files from the `.epub` file archive. By default the archive files are extracted to `tempdir()` so you may want to change this with the `exdir` argument.

```{r unzip}
bookdir <- file.path(tempdir(), "dracula")
epub_unzip(file, exdir = bookdir)
list.files(bookdir, recursive = TRUE)
```

## Other functions

### Word count

The helper function `count_words()` provides word counts for strings, but allows you to control the regular expression patterns used for both splitting the string and conditionally counting the resulting character elements. This is the same function used internally by `epub()` and `epub_recombine()`. It is exported so that it can be used directly.

By default, `count_words()` splits on spaces and new line characters. It counts as a word any element containing at least one alphanumeric character or the ampersand. It ignores everything else as noise, such as extra spaces, empty strings and isolated bits of punctuation.

```{r count_words}
x <- " This   sentence will be counted to have:\n\n10 (ten) words."
count_words(x)
```

### Inspection

Helper functions for inspecting the text in the R console include `epub_head()` and `epub_cat()`. 

`epub_head()` provides an overview of the text by section for each book in the primary data frame. The nested data frames are unnested and row bound to one another and returned as a single data frame. The text is shortened to only the first few characters (defaults to `n = 50`).

`epub_cat()` can be used to `cat` the text of an e-book to the console for quick inspection in a more readable form. It can take several arguments that help slice out a section of the text and customize how it is printed.

Both functions can take an EPUB filename or a data frame of an already loaded EPUB file as their first argument.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sift.R
\name{epub_sift}
\alias{epub_sift}
\title{Sift EPUB sections}
\usage{
epub_sift(data, n, type = c("word", "char"))
}
\arguments{
\item{data}{a data frame created by \code{epub}.}

\item{n}{integer, minimum number of words or characters to retain a section.}

\item{type}{character, \code{"word"} or \code{"character"}.}
}
\value{
a data frame
}
\description{
Sift out EPUB sections that have suspiciously low word or character count.
}
\details{
This function is like a sieve that lets small section rows fall through.
Choose the minimum number of words or characters to accept as a meaningful section in the e-book worth retaining in the nested data frame, e.g., book chapters.
Data frame rows pertaining to smaller sections are dropped.

This function is helpful for isolating meaningful content by removing extraneous e-book sections that may be difficult to remove by other methods when working with poorly formatted e-books.
The EPUB file included in \code{epubr} is a good example of this. It does not contain meaningful section identifiers in its metadata.
This creates a need to restructure the text table after reading it with \code{epub} by subsequently calling \code{epub_recombine}.
However, some unavoidable ambiguity in this leads to many small sections appearing from the table of contents.
These can then be dropped with \code{epub_sift}. See a more comprehensive in the \code{\link{epub_recombine}} documentation.
A simpler example is shown below.
}
\examples{
\donttest{
file <- system.file("dracula.epub", package = "epubr")
x <- epub(file) # parse entire e-book
x$data[[1]]

x <- epub_sift(x, n = 3000) # drops last two sections
x$data[[1]]
}
}
\seealso{
\code{\link{epub_recombine}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{epub_head}
\alias{epub_head}
\title{Preview the first n characters}
\usage{
epub_head(x, n = 50)
}
\arguments{
\item{x}{a data frame returned by \code{\link{epub}} or a character string giving the EPUB filename(s).}

\item{n}{integer, first n characters to retain from each e-book section.}
}
\value{
a data frame.
}
\description{
Preview the first n characters of each EPUB e-book section.
}
\details{
This function returns a simplified data frame of only the unnested \code{section} and \code{text} columns of a data frame returned by \code{\link{epub}}, with the text included only up to the first \code{n} characters.
This is useful for previewing the opening text of each e-book section to inspect for possible useful regular expression patterns to use for text-based section identification.
For example, an e-book may not have meaningful section IDs that distinguish one type of book section from another, such as chapters from non-chapter sections,
but the text itself may contain this information at or near the start of a section.
}
\examples{
\donttest{
file <- system.file("dracula.epub", package = "epubr")
epub_head(file)
}
}
\seealso{
\code{\link{epub_cat}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reorder.R
\name{epub_reorder}
\alias{epub_reorder}
\title{Reorder sections}
\usage{
epub_reorder(data, .f, pattern)
}
\arguments{
\item{data}{a data frame created by \code{epub}.}

\item{.f}{a scalar function to determine a single row index based on a matched regular expression. It must take two strings, the text and the pattern, and return a single number. See examples.}

\item{pattern}{regular expression passed to \code{.f}.}
}
\value{
a data frame
}
\description{
Reorder text sections in an e-book based on a user-provided function.
}
\details{
Many e-books have chronologically ordered sections based on quality metadata.
This results in properly book sections in the nested data frame.
However, some poorly formatted e-books have their internal sections occur in an arbitrary order.
This can be frustrating to work with when doing text analysis on each section and where order matters.

This function addresses this case by reordering the text sections in the nested data frame based on a user-provided function that re-indexes the data frame rows based on their content.
In general, the approach is to find something in the content of each section that describes the section order.
For example, \code{epub_recombine} can use a regular expression to identify chapters.
Taking this a step further, \code{epub_reorder} can use a function that works with the same information to reorder the rows.

It is enough in the former case to identify where in the text the pattern occurs. There is no need to extract numeric ordering from it.
The latter takes more effort. In the example EPUB file included in \code{epubr}, chapters can be identified using a pattern of the word CHAPTER in capital letters followed by a space and then some Roman numerals.
The user must provide a function that would parse the Roman numerals in this pattern so that the rows of the data frame can be reordered properly.
}
\examples{
\donttest{
file <- system.file("dracula.epub", package = "epubr")
x <- epub(file) # parse entire e-book
x <- epub_recombine(x, "CHAPTER [IVX]+", sift = list(n = 1000)) # clean up

library(dplyr)
set.seed(1)
x$data[[1]] <- sample_frac(x$data[[1]]) # randomize rows for example
x$data[[1]]

f <- function(x, pattern) as.numeric(as.roman(gsub(pattern, "\\\\1", x)))
x <- epub_reorder(x, f, "^CHAPTER ([IVX]+).*")
x$data[[1]]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recombine.R
\name{epub_recombine}
\alias{epub_recombine}
\title{Recombine text sections}
\usage{
epub_recombine(data, pattern, sift = NULL)
}
\arguments{
\item{data}{a data frame created by \code{epub}.}

\item{pattern}{character, a regular expression.}

\item{sift}{\code{NULL} or a named list of parameters passed to \code{\link{epub_sift}}. See details.}
}
\value{
a data frame
}
\description{
Split and recombine EPUB text sections based on regular expression pattern matching.
}
\details{
This function takes a regular expression and uses it to determine new break points for the full e-book text.
This is particularly useful when sections pulled from EPUB metadata have arbitrary breaks and the text contains meaningful breaks at random locations in various sections.
\code{epub_recombine} collapses the text and then creates a new nested data frame containing new chapter/section labels, word counts and character counts,
associated with the text based on the new break points.

Usefulness depends on the quality of the e-book. While this function exists to improve the data structure of e-book content parsed from e-books with poor metadata formatting,
it still requires original formatting that will at least allow such an operation to be successful, specifically a consistent, non-ambiguous regular expression pattern.
See examples below using the built in e-book dataset.

When used in conjunction with \code{epub_sift} via the \code{sift} argument, recombining and resifting is done recursively.
This is because it is possible that sifting can create a need to rerun the recombine step in order to regenerate proper chapter indexing for the section column.
However, recombining a second time does not lead to a need to resift, so recursion ends after one round regardless.

This is a convenient way to avoid the syntax:

\code{epub_recombine([args]) \%>\% epub_sift([args]) \%>\% epub_recombine([args])}.
}
\examples{
\donttest{
file <- system.file("dracula.epub", package = "epubr")
x <- epub(file) # parse entire e-book
x$data[[1]] # note arbitrary section breaks (not between chapters)

pat <- "CHAPTER [IVX]+" # but a reliable pattern exists for new breaks
epub_recombine(x, pat) # not as expected; pattern also in table of contents

epub_recombine(x, pat, sift = list(n = 1000)) # sift low word-count sections
}
}
\seealso{
\code{\link{epub_sift}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epub.R
\name{epub}
\alias{epub}
\alias{epub_meta}
\alias{epub_unzip}
\title{Extract and read EPUB e-books}
\usage{
epub(
  file,
  fields = NULL,
  drop_sections = NULL,
  chapter_pattern = NULL,
  encoding = "UTF-8",
  ...
)

epub_meta(file)

epub_unzip(file, exdir = tempdir())
}
\arguments{
\item{file}{character, input EPUB filename. May be a vector for \code{epub} and \code{epub_meta}. Always a single file for \code{epub_unzip}.}

\item{fields}{character, vector of metadata fields (data frame columns) to parse from metadata, if they exist. See details.}

\item{drop_sections}{character, a regular expression pattern string to identify text sections (rows of nested text data frame) to drop.}

\item{chapter_pattern}{character, a regular expression pattern string to attempt distinguishing nested data frame rows of chapter text entries from other types of entries.}

\item{encoding}{character, defaults to \code{"UTF-8"}.}

\item{...}{additional arguments. With the exception of passing \code{title} (see details), currently developmental/unsupported.}

\item{exdir}{for \code{epub_unzip}, extraction directory to place archive contents (files). It will be created if necessary.}
}
\value{
\code{epub} returns a data frame. \code{epub_unzip} returns nothing but extracts files from an EPUB file archive.
}
\description{
Read EPUB format e-books into a data frame using \code{epub} or extract EPUB archive files for direct use with \code{epub_unzip}.
}
\details{
The primary function here is \code{epub}. It parses EPUB file metadata and textual content into a data frame.
The output data frame has one row for each file in \code{file}.
It has metadata in all columns except the \code{data} column, which is a column of nested data frames containing e-book text by book section (e.g., chapters).
Both the primary and nested data frames are tibbles and safe to print to the console "as is".

Be careful if \code{file} is a long vector of many EPUB files.
This could take a long time to process as well as could potentially use up all of your system RAM if you have far too many large books in one call to \code{epub}.

On a case by case basis, you can always select columns and filter rows of a resulting data frame for a single e-book subsequent to visual inspection.
However, the optional arguments \code{fields}, \code{drop_sections} and \code{chapter_pattern} allow you to do some of this as part of the EPUB file reading process.
You can ignore these arguments and do all your own post-processing of the resulting data frame, but if using these arguments,
they are most likely to be useful for bulk e-book processing where \code{file} is a vector of like-formatted files.

\subsection{Main columns}{
The \code{fields} argument can be used to limit the columns returned in the primary data frame.
E.g., \code{fields = c("title", "creator", "date", "identifier", "publisher", "file")}. Some fields will be returned even if not in \code{fields}, such as \code{data} and \code{title}.
\cr\cr
Ideally, you should already know what metadata fields are in the EPUB file. This is not possible for large collections with possibly different formatting.
Note that when \code{"file"} is included in \code{fields}, the output will include a column of the original file names, in case this is different from the content of a \code{source} field that may be present in the metadata.
So this field is always available even if not part of the file metadata.
\cr\cr
Additionally, if there is no \code{title} field in the metadata, the output data frame will include a \code{title} column filled in with the same file names,
unless you pass the additional optional title argument, e.g. \code{title = "TitleFieldID"} so that another field can me mapped to \code{title}.
If supplying a \code{title} argument that also does not match an existing field in the e-book, the output \code{title} column will again default to file names.
File names are the fallback option because unlike e-book metadata fields, file names always exist and should also always be unique when performing vectorized reads over multiple books,
ensuring that \code{title} can be a column in the output data frame that uniquely identifies different e-books even if the books did not have a \code{title} field in their metadata.
\cr\cr
Columns of the nested data frames in \code{data} are fixed. Select from these in subsequent data frame manipulations.
}
\subsection{Nested rows}{
The \code{chapter_pattern} argument may be helpful for bulk processing of similarly formatted EPUB files. This should be ignored for poorly formatted EPUB files or where there is inconsistent naming across an e-book collection.
Like with \code{fields}, you should explore file metadata in advance or this argument will not be useful. If provided, a column \code{nchap} is added to the output data frame giving the guessed number of chapters.
In the \code{data} column, the \code{section} column of the nested data frames will also be updated to reflect guessed chapters with new, consistent chapter IDs, always beginning with \code{ch} and ending with digits.
\cr\cr
The \code{drop_sections} argument also uses regular expression pattern matching like \code{chapter_pattern} and operates on the same \code{section} column. It simply filters out any matched rows.
This is useful for dropping rows that may pertain to book cover, copyright and acknowledgements pages, and other similar, clearly non-chapter text, e-book sections.
An example that might work for many books could be \code{drop_sections = "^co(v|p)|^ack"}
\cr\cr
Rows of the primary data frame are fixed. Filter or otherwise manipulate these in subsequent data frame manipulations. There is one row per file so filtering does not make sense to do as part of the initial file reading.
}
\subsection{EPUB metadata}{
Use \code{epub_meta} to return a data frame of only the metadata for each file in \code{file}. This skips the reading of each file's text contents, strictly parsing the metadata.
It returns a data frame with one row for each file and \code{n} columns where \code{n} is equal to the union of all fields identified across all files in \code{file}.
Fields available for at least one e-book in \code{file} will return \code{NA} in that column for any row pertaining to an e-book that does not have that field in its metadata.
If the metadata contains multiple entries for a field, such as multiple subjects or publication dates, they are collapsed using the pipe character.
}
\subsection{Unzipping EPUB files}{
If using \code{epub_unzip} directly on individual EPUB files, this gives you control over where to extract archive files to and what to do with them subsequently.
\code{epub} and \code{epub_meta} use \code{epub_unzip} internally to extract EPUB archive files to the R session temp directory (with \code{tempdir()}).
You do not need to use \code{epub_unzip} directly prior to using these other functions. It is only needed if you want the internal files for some other purpose in or out of R.
}
}
\examples{
# Use a local example EPUB file included in the package
file <- system.file("dracula.epub", package = "epubr")
bookdir <- file.path(tempdir(), "dracula")
epub_unzip(file, exdir = bookdir) # unzip to directly inspect archive files
list.files(bookdir, recursive = TRUE)

\donttest{
epub_meta(file) # parse EPUB file metadata only

x <- epub(file) # parse entire e-book
x
x$data[[1]]

epub(file, fields = c("title", "creator"), drop_sections = "^cov")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{epub_cat}
\alias{epub_cat}
\title{Pretty printing of EPUB text}
\usage{
epub_cat(
  x,
  max_paragraphs = 10,
  skip = 0,
  paragraph_spacing = 1,
  paragraph_indent = 2,
  section_sep = "====",
  book_sep = "====\\n===="
)
}
\arguments{
\item{x}{a data frame returned by \code{\link{epub}} or a character string giving the EPUB filename(s).}

\item{max_paragraphs}{integer, maximum number of paragraphs (non-empty lines) to \code{cat} to console.}

\item{skip}{integer, number of paragraphs to skip.}

\item{paragraph_spacing}{integer, number of empty lines between paragraphs.}

\item{paragraph_indent}{integer, number of spaces to indent paragraphs.}

\item{section_sep}{character, a string to indicate section breaks.}

\item{book_sep}{character, separator shown between books when \code{x} has multiple rows (books).}
}
\value{
nothing is returned but a more readable format of the text content for books in \code{x} is printed to the console.
}
\description{
Print EPUB text to the console in a more readable format.
}
\details{
This function prints text from EPUB files to the console using \code{cat}.
This is useful for quickly obtaining an overview of the book text parsed by \code{\link{epub}} that is easier to read that looking at strings in the table.
\code{max_paragraphs} is set low by default to prevent accidentally printing entire books to the console.
To print everything in \code{x}, set \code{max_paragraphs = NULL}.
}
\examples{
\donttest{
file <- system.file("dracula.epub", package = "epubr")
d <- epub(file)
epub_cat(d, max_paragraphs = 2, skip = 147)
}
}
\seealso{
\code{\link{epub_head}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{count_words}
\alias{count_words}
\title{Word count}
\usage{
count_words(x, word_pattern = "[A-Za-z0-9&]", break_pattern = " |\\n")
}
\arguments{
\item{x}{character, a string containing words to be counted. May be a vector.}

\item{word_pattern}{character, regular expression to match words. Elements not matched are not counted.}

\item{break_pattern}{character, regular expression to split a string between words.}
}
\value{
an integer
}
\description{
Count the number of words in a string.
}
\details{
This function estimates the number of words in strings. Words are first separated using \code{break_pattern}.
Then the resulting character vector elements are counted, including only those that are matched by \code{word_pattern}.
The approach taken is meant to be simple and flexible.

\code{epub} uses this function internally to estimate the number of words for each e-book section alongside the use of \code{nchar} for counting individual characters.
It can be used directly on character strings and is convenient for applying with different regular expression pattern arguments as needed.

These two arguments are provided for control, but the defaults are likely good enough.
By default, strings are split only on spaces and new line characters.
The "words" that are counted in the resulting vector are those that contain any alphanumeric characters or the ampersand.
This means for example that hyphenated words, acronyms and numbers displayed with digits, are all counted as words.
The presence of any other characters does not negate that a word has been found.
}
\examples{
x <- " This   sentence will be counted to have:\n\n10 (ten) words."
count_words(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epubr.R
\docType{package}
\name{epubr}
\alias{epubr}
\title{epubr: Read EPUB file metadata and text.}
\description{
\code{epubr} provides functions supporting the reading and parsing of internal e-book content from EPUB files.
E-book metadata and textual contents are parsed separately.
}
\details{
E-book formatting is not completely standardized across all literature.
It can be challenging to curate parsed e-book content across an arbitrary collection of e-books perfectly and in completely general form, to yield a singular, consistently formatted output.
Many EPUB files do not even contain all the same pieces of information in their respective metadata.

EPUB file parsing functionality in this package is intended for relatively general application to arbitrary EPUB e-books.
However, poorly formatted e-books or e-books with highly uncommon formatting may not work with this package.
There may even be cases where an EPUB file has DRM or some other property that makes it impossible to read with \code{epubr}.

Text is read as is for the most part. The only nominal changes are minor substitutions, for example curly quotes changed to straight quotes.
Substantive changes are expected to be performed subsequently by the user as part of their text analysis.
Additional text cleaning can be performed at the user's discretion, such as with functions from packages like \code{tm} or \code{qdap}.
}
