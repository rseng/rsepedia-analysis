---
title: 'jstor: Import and Analyse Data from Scientific Texts'
authors:
 - affiliation: 1
   name: Thomas Klebel
   orcid: 0000-0002-7331-4751
date: 2018-08-07
output:
  html_document:
    df_print: paged
bibliography: paper.bib
tags:
 - JSTOR
 - DfR
 - Data for Research
 - scientometrics
 - bibliometrics
 - text mining
 - text analysis
 - citation analysis
affiliations:
 - index: 1
   name: Department of Sociology, University of Graz
---

# Summary
The interest in text as data has seen a sharp increase in the 
past few years, mostly due to the advent of methods for automated text analysis.
At the same time, researches within the field of scientometrics have analysed
citations and other aspects of the scholarly literature with great sophistication.
The archival content of [JSTOR](http://www.jstor.org) offers a rich and diverse
set of primary sources like research articles or book chapters for both 
approaches. 
[Data for Research (DfR)](http://www.jstor.org/dfr/) by JSTOR gives all 
researchers, regardless of whether they have access to JSTOR or not, the
opportunity to analyse metadata,
n-grams and, upon special request, full-text materials about all available
articles and books from JSTOR. The package `jstor` [@jstor] helps in
analysing these datasets by enabling researchers to easily import the metadata
to R [@r_core], a task, for which no other integrated solution exists to date.

The metadata from DfR
can either be analysed on their own or be used in conjunction with n-grams
or full-text data. Commonly, metadata from DfR include information
on the articles' authors, their title, journal, date of publishing, and quite
frequently all footnotes and references. All this information can be of interest
for specific research questions. For the analysis of n-grams or full-texts,
the metadata imported with `jstor`
allow the researchers to
filter articles based on specific journals, the dates of publication, the
authors, keywords in titles and other aspects.

`jstor` provides functions for three main tasks within the research process:

- Importing different parts of metadata, either from XML-files or directly from
the .zip-archive provided by DfR.
- Importing n-gram and full-text files.
- Performing common tasks of cleaning metadata like unifying the journal id or
cleaning page numbers.


Full documentation of `jstor`, including a comprehensive 
case study about analysing 
n-grams from DfR, is available at 
https://ropensci.github.io/jstor/. The package can be obtained from 
CRAN (https://CRAN.R-project.org/package=jstor)
or from GitHub (https://github.com/ropensci/jstor). 
Archived versions of all releases are available at Zenodo 
(https://doi.org/10.5281/zenodo.1169861). 





# Acknowledgements
I am indebted to Jason Becker and Elin Waring for their helpful requests during 
[package review](https://github.com/ropensci/onboarding/issues/189). 
Benjamin Klebel, Antonia Schirgi and Matthias Duller provided helpful comments
and requests for clarifications when writing the case study.

Work on `jstor` benefited from financial support for the project "Academic
Super-Elites in Sociology and Economics" by the Austrian Science Fund (FWF), 
project number "P 29211 Einzelprojekte". `jstor` is currently being used by
the project team to analyse academic elites in sociology and economics.



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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/

<!-- README.md is generated from README.Rmd. Please edit that file -->

# jstor: Import and Analyse Data from Scientific Articles

**Author:** [Thomas Klebel](https://thomasklebel.eu) <br> **License:**
[GPL v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)

[![R-CMD-check](https://github.com/ropensci/jstor/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/ropensci/jstor/actions/workflows/check-standard.yaml)
[![AppVeyorBuild
status](https://ci.appveyor.com/api/projects/status/sry2gtwam7qyfw6l?svg=true)](https://ci.appveyor.com/project/tklebel/jstor)
[![Coverage
status](https://codecov.io/gh/ropensci/jstor/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/jstor?branch=master)
[![lifecycle](https://img.shields.io/badge/lifecycle-superseded-blue.svg)](https://www.tidyverse.org/lifecycle/#superseded)
[![CRAN
status](http://www.r-pkg.org/badges/version/jstor)](https://cran.r-project.org/package=jstor)
[![CRAN\_Download\_Badge](http://cranlogs.r-pkg.org/badges/grand-total/jstor)](https://CRAN.R-project.org/package=jstor)
[![rOpenSci
badge](https://badges.ropensci.org/189_status.svg)](https://github.com/ropensci/onboarding/issues/189)
[![JOSS
badge](http://joss.theoj.org/papers/ba29665c4bff35c37c0ef68cfe356e44/status.svg)](http://joss.theoj.org/papers/ba29665c4bff35c37c0ef68cfe356e44)
[![Zenodo
DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1169861.svg)](https://doi.org/10.5281/zenodo.1169861)

The tool [Data for Research (DfR)](http://www.jstor.org/dfr/) by JSTOR
is a valuable source for citation analysis and text mining. `jstor`
provides functions and suggests workflows for importing datasets from
DfR. It was developed to deal with very large datasets which require an
agreement, but can be used with smaller ones as well.

**Note**: As of 2021, JSTOR has moved changed the way they provide data to a new 
platform called [Constellate](https://constellate.org/). The package `jstor` has
not been adapted to this change, and might therefore only be used for legacy 
data that was optained from the old DfR platform.

The most important set of functions is a group of `jst_get_*` functions:

  - `jst_get_article`
  - `jst_get_authors`
  - `jst_get_references`
  - `jst_get_footnotes`
  - `jst_get_book`
  - `jst_get_chapters`
  - `jst_get_full_text`
  - `jst_get_ngram`

All functions which are concerned with meta data (therefore excluding
`jst_get_full_text` and `jst_get_ngram`) operate along the same lines:

1.  The file is read with `xml2::read_xml()`.
2.  Content of the file is extracted via XPATH or CSS-expressions.
3.  The resulting data is returned in a `tibble`.

## Installation

To install the package use:

``` r
install.packages("jstor")
```

You can install the development version from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/jstor")
```

## Usage

In order to use `jstor`, you first need to load it:

``` r
library(jstor)
library(magrittr)
```

The basic usage is simple: supply one of the `jst_get_*`-functions with
a path and it will return a tibble with the extracted
information.

``` r
jst_get_article(jst_example("article_with_references.xml")) %>% knitr::kable()
```

| file\_name                | journal\_doi | journal\_jcode   | journal\_pub\_id | journal\_title                                     | article\_doi    | article\_pub\_id | article\_jcode | article\_type    | article\_title                     | volume | issue | language | pub\_day | pub\_month | pub\_year | first\_page | last\_page | page\_range |
| :------------------------ | :----------- | :--------------- | :--------------- | :------------------------------------------------- | :-------------- | :--------------- | :------------- | :--------------- | :--------------------------------- | :----- | :---- | :------- | :------- | :--------- | --------: | :---------- | :--------- | :---------- |
| article\_with\_references | NA           | tranamermicrsoci | NA               | Transactions of the American Microscopical Society | 10.2307/3221896 | NA               | NA             | research-article | On the Protozoa Parasitic in Frogs | 41     | 2     | eng      | 1        | 4          |      1922 | 59          | 76         | 59-76       |

``` r

jst_get_authors(jst_example("article_with_references.xml")) %>% knitr::kable()
```

| file\_name                | prefix | given\_name | surname | string\_name | suffix | author\_number |
| :------------------------ | :----- | :---------- | :------ | :----------- | :----- | -------------: |
| article\_with\_references | NA     | R.          | Kudo    | NA           | NA     |              1 |

Further explanations, especially on how to use jstor’s functions for
importing many files, can be found in the vignettes.

## Getting started

In order to use `jstor`, you need some data from DfR. From the [main
page](http://www.jstor.org/dfr/) you can create a dataset by searching
for terms and restricting the search regarding time, subject and content
type. After you created an account, you can download your selection.
Alternatively, you can download [sample
datasets](http://www.jstor.org/dfr/about/sample-datasets) with documents
from before 1923 for the US, and before 1870 for all other countries.

## Supported Elements

In their [technical
specifications](http://www.jstor.org/dfr/about/technical-specifications),
DfR lists fields which should be reliably present in all articles and
books.

The following table gives an overview, which elements are supported by
`jstor`.

### Articles

| `xml`-field                      | reliably present | supported in `jstor` |
| :------------------------------- | :--------------- | :------------------- |
| journal-id (type=“jstor”)        | x                | x                    |
| journal-id (type=“publisher-id”) | x                | x                    |
| journal-id (type=“doi”)          |                  | x                    |
| issn                             | x                |                      |
| journal-title                    | x                | x                    |
| publisher-name                   | x                |                      |
| article-id (type=“doi”)          | x                | x                    |
| article-id (type=“jstor”)        | x                | x                    |
| article-id (type=“publisher-id”) |                  | x                    |
| article-type                     |                  | x                    |
| volume                           |                  | x                    |
| issue                            |                  | x                    |
| article-categories               | x                |                      |
| article-title                    | x                | x                    |
| contrib-group                    | x                | x                    |
| pub-date                         | x                | x                    |
| fpage                            | x                | x                    |
| lpage                            |                  | x                    |
| page-range                       |                  | x                    |
| product                          | x                |                      |
| self-uri                         | x                |                      |
| kwd-group                        | x                |                      |
| custom-meta-group                | x                | x                    |
| fn-group (footnotes)             |                  | x                    |
| ref-list (references)            |                  | x                    |

### Books

| `xml`-field            | reliably present | supported in `jstor` |
| :--------------------- | :--------------- | :------------------- |
| book-id (type=“jstor”) | x                | x                    |
| discipline             | x                | x                    |
| call-number            | x                |                      |
| lcsh                   | x                |                      |
| book-title             | x                | x                    |
| book-subtitle          |                  | x                    |
| contrib-group          | x                | x                    |
| pub-date               | x                | x                    |
| isbn                   | x                | x                    |
| publisher-name         | x                | x                    |
| publisher-loc          | x                | x                    |
| permissions            | x                |                      |
| self-uri               | x                |                      |
| counts                 | x                | x                    |
| custom-meta-group      | x                | x                    |

### Book Chapters

| `xml`-field            | reliably present | supported in `jstor` |
| :--------------------- | :--------------- | :------------------- |
| book-id (type=“jstor”) | x                | x                    |
| part\_id               | x                | x                    |
| part\_label            | x                | x                    |
| part-title             | x                | x                    |
| part-subtitle          |                  | x                    |
| contrib-group          | x                | x                    |
| fpage                  | x                | x                    |
| abstract               | x                | x                    |

## Code of conduct

Please note that this project is released with a [Contributor Code of
Conduct](CONDUCT.md). By participating in this project you agree to
abide by its terms.

## Citation

To cite `jstor`, please refer to `citation(package =
    "jstor")`:

    Klebel (2018). jstor: Import and Analyse Data from Scientific Texts. Journal of 
    Open Source Software, 3(28), 883, https://doi.org/10.21105/joss.00883

## Acknowledgements

Work on `jstor` benefited from financial support for the project
“Academic Super-Elites in Sociology and Economics” by the Austrian
Science Fund (FWF), project number “P 29211 Einzelprojekte”.

Some internal functions regarding file paths and example files were
adapted from the package
`readr`.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# jstor 0.3.10
This is a small release to fix compatibility with `readr v2.0.0`. In addition,
various minor improvements have been made across to package to get it back to
CRAN.

# jstor 0.3.9
This is a small release to fix compatibility with `dplyr v1.0.0`.

# jstor 0.3.8
* The package homepage has switched from https://ropensci.github.io/jstor to
https://docs.ropensci.org/jstor. This simplifies building the site and aligns it
with the rest of rOpenSci's packages. #80
* Compatibility fix for tibble 3.0.0.

# jstor 0.3.7
This is a small release to fix compatibility with `tidyr v1.0.0`. Furthermore,
the formerly defunct functions following the old naming conventions (like
`find_article()`, `find_references()`, etc.) have been removed.

# jstor 0.3.6
This is another small release to fix compatibility with `readr v1.3.0` and
`tibble v2.0.0`. There are no other changes.


# jstor 0.3.5
This is a small release, mainly to fix compatibility with version `1.2.0` of 
`readr`. There is one breaking change however:

## Breaking changes
* Output column names for `jst_get_refernces` have been renamed to avoid 
ambiguity when matching with output from `jst_get_article`. All columns now have
a `ref_*` prefix.


# jstor 0.3.4
* Added option to parse references, if the information is available. #32
* Output from error messages is "re-imported" properly as well.
* Example files have been renamed for clarity. `sample_` has been replaced with
`article_` or removed altogether.
* Reworked internals for catching misspecifications in `jst_define_import` due 
to upcoming release of `rlang v0.3.0`.

# jstor 0.3.3

## Removed functionality
* The option to download the most recent data on journals directly from JSTOR
(as in `jst_get_journal_overview(most_recent = T)`)
had to be removed due changes on their server. I will try to find a solution
with JSTOR support so we can add the functionality again.

## New features
* `jst_define_import` now prints the specification in a pretty and informative
way.
* `jst_define_import` now checks the definition more extensively: 
`jst_define_import(article = jst_get_book)` or similar mis-specifications 
will raise an error.
* Import the crayon package for more colourful error messages in some places,
which are easier to read.

## Bug fixes
* Removed an outdated function from the vignette on batch importing files.

## Other changes
* The old group of `find_*` functions is now defunct (they raise an error).
* Updated cached version of journal data.

# jstor 0.3.2
This is a hotfix to resolve an issue with writing to other directories than
temporary folders during tests, which should not have happend in the first 
place.


# jstor 0.3.1

* `jstor` is now part of rOpenSci.
* removed arguments `new_col` for `jst_unify_journal_id` and
`jst_add_total_pages`, since both built on the dev version of rlang. Once this
version is on CRAN, they will be re-introduced.
* fixed a few issues in README and man files regarding changes introduced with 
0.3.0.

# jstor 0.3.0

## Breaking changes
`jst_import` and `jst_import_zip` now use futures as a backend for parallel 
processing. This makes internals more compact and reduces dependencies. 
Furthermore this reduces the number of arguments, since the argument `cores` 
has been removed. By default, the functions run sequentially. If you want them
to execute in parallel, use futures:
```
library(future)
plan(multiprocess)

jst_import_zip("zip-archive.zip",
               import_spec = jst_define_import(article = jst_get_article),
               out_file = "outfile")
```
If you want to terminate the proceses, at least on *nix-systems you need to kill
them manually (once again).

* All functions have been renamed to use a unified naming scheme: `jst_*`.
The former group of `find_*` functions is now called `jst_get_*`, as in
`jst_get_article()`. The previous functions have been deprecated and will be 
removed before submission to CRAN.
* The unique identifier for matching across files has been renamed to 
`file_name`, and the corresponding helper to get this file name from
`get_basename` to `jst_get_file_name`.

## Importing data directly from zip-files
There is a new set of functions which lets you directly import files from 
.zip-archives: `jst_import_zip()` and `jst_define_import()`.

In the following example, we have a zip-archive from DfR and want to import
metadata on books and articles. For all articles we want to apply
`jst_get_article()` and `jst_get_authors()`, for books only `jst_get_book()`, 
and we want to read unigrams (ngram1).

First we specify what we want, and then we apply it to our zip-archive:
```r
# specify definition
import_spec <- jst_define_import(article = c(jst_get_article, jst_get_authors),
                                 book = jst_get_book,
                                 ngram1 = jst_get_ngram)

# apply definition to archive
jst_import_zip("zip_archive.zip",
               import_spec = import_spec,
               out_file = "out_path")
```

If the archive contains also research reports, pamphlets or other ngrams, they
will not be imported. We could however change our specification, if we wanted
to import all kinds of ngrams (given that we originally requested them from
DfR):

```r
# import multiple forms of ngrams
import_spec <- jst_define_import(article = c(jst_get_article, jst_get_authors),
                                 book = jst_get_book,
                                 ngram1 = jst_get_ngram,
                                 ngram2 = jst_get_ngram,
                                 ngram3 = jst_get_ngram)
```

Note however that for larger archives, importing all ngrams takes a very long
time. It is thus advisable to only import ngrams for articles which you
want to analyse, i.e. most likely a subset of the initial request. The new 
function `jst_subset_ngrams()` helps you with this (see also the section on
importing bigrams in the 
[case study](https://docs.ropensci.org/jstor/articles/analysing-n-grams.html#importing-bigrams).

Before importing all files from a zip-archive, you can get a quick overview with
`jst_preview_zip()`.

## New vignette
The new `vignette("known-quirks")` lists common problems with data from
JSTOR/DfR. Contributions with further cases are welcome!


## New functions
* New function `jst_get_journal_overview()` supplies a tibble with contextual
information about the journals in JSTOR.
* New function `jst_combine_outputs()` applies `jst_re_import()` to a whole 
directory and lets you combine all related files in one go. It uses the file 
structure that `jst_import()` and `jst_import_zip()` provide as a heuristic: a
filename with a dash and one or multiple digits at its end (`filename-1.csv`). 
All files
with identical names (disregarding dash and digits) are combined into one file.
* New function `jst_re_import()` lets you re_import a `.csv` file that
`jstor_import()` or `jst_import_zip()` had exported. It tries to guess the type
of
content based on the column names or, if column names are not available, from
the number of columns, raising a warning if guessing fails and reverting to a
generic import.
* A new function `jst_subset_ngrams()` lets you create a subset of ngram files
within a zip-file which you can import with `jst_get_ngram()`.
* A new set of convenience functions for taking a few cleaning steps:
`jst_clean_page()` tries to turn a character vector with pages into a numeric
one, `jst_unify_journal_id()` merges different specifications of journals into
one, `jst_add_total_pages()` adds a total count of pages per article, and
`jst_augment()` calls all three functions to clean the data set in one go.


## Minor changes
* Improved documentation regarding endnotes (thanks @elinw)
* jst_import and jst_import_zip have a new argument: `n_batches` which lets you 
specify the number of batches directly


# jstor 0.2.6

* added lengthy case study at https://tklebel.github.io/jstor/articles/analysing-n-grams.html
* added a pkgdown site at https://tklebel.github.io/jstor/
* changed implementation of parallel execution in `jstor_import` from 
`parallel::mclapply` to `foreach::foreach` with `snow` as a backend for
`%dopar%`. 
* added support for progress bars #34
* `jstor_import` now writes column names by default #29
* new helper `get_basename` helps to get the basename of a file without its
extension
* `find_article` does not coerce days and months to integer any more, since there
might be information stored as text.
* Added a `NEWS.md` file to track changes to the package.
## Test environments
* local MacOS 11.6.1, R 4.0.4
* ubuntu 20.04 (GitHub actions), devel and release
* macOS-latest (GitHub actions), release
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a fix addressing changes in readr (v2.0.0). Changes in readr where the
reason for the check errors that led to the removal of `jstor` from CRAN.
As the tests and checks noted above indicate, all related issues have been
resolved.

* The notes regarding the links are spurious (jstor.org blocking access), as in
previous submissions.

* Changed all occurences of `T` or `F` to `TRUE` and `FALSE`, as requested.

* Removed file LICENSE and mention of it from DESCRIPTION, as requested.

* Added `return` statements to the documentation of all requested files.


# extracting footnotes works

     [1] "[Footnotes]"                                                           
     [2] "1\nJudge Edwards's (1998)\nD.C. Circuit (1997)\nRevesz's reply (1999)."
     [3] "2\nRaz (1994),\nSchauer (1991)\nSoper (1996, 1998)."                   
     [4] "3\nCooter (1983)\nPosner (1993)."                                      
     [5] "4\nHacker (1973)"                                                      
     [6] "5\nSavage (1954)."                                                     
     [7] "6\nKornhauser (1998)."                                                 
     [8] "9\nRasmusen (1996)."                                                   
     [9] "10\nRabin (1995)"                                                      
    [10] "12\nHiggins and Rubin (1980)\nCohen (1991,\n1992)."                    

---

    [1] "[Footnotes]"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    [2] "9\nR. E. Bellman, Dynamic Programming (Princeton, N.J.: Princeton University Press,\n1957)Bellman\n                  Dynamic Programming\n                  1957\n               "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
    [3] "10\nC.\nR. Carr and C. W. Howe, Quantitative Decision Procedures in Management and Economics\n(New York: McGraw-Hill 1964), Ch. 7Carr\n                  Quantitative Decision Procedures in Management and Economics\n                  1964\n               \nW. J. Fabrycky and P. E. Torgersen, Operations\nEconomy: Industrial Applications of Operations Research (Englewood Cliffs, N. J.: Prentice-\nHall, 1966), Ch. 16Fabrycky\n                  Operations Economy: Industrial Applications of Operations Research\n                  1966\n               \nG. Hadley, Nonlinear and Dynamic Programming (Reading: Addison-\nWesley, 1964), Chs. 10 and 11Hadley\n                  Nonlinear and Dynamic Programming\n                  1964\n               \nR. A. Howard, \"Dynamic Programming,\" Management\nScience, XII, No. 5 (Jan., 1966), p. 31710.2307/2627818\n                  317\n               \nG. L. Nemhauser, Dynamic Programming (New York:\nJohn Wiley and Sons, 1966)Nemhauser\n                  Dynamic Programming\n                  1966\n               \nM. Sasieni, A. Yaspan, L. Friedman, Operations Research\n(New York: John Wiley and Sons, 1959), Ch. 10Sasieni\n                  Operations Research\n                  1959\n               "

# full-text is read correctly

    Code
      text[["full_text"]]
    Output
      [1] "<plain_text> <page sequence=\"1\">Lorem ipsum dolor sit amet, \r\nconsectetur adipisici elit, sed eiusmod tempor incidunt ut labore et dolore \r\nmagna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris\r\nnisi ut aliquid ex ea commodi consequat. Quis aute iure reprehenderit in \r\nvoluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint\r\nobcaecat cupiditat non proident, sunt in culpa qui officia deserunt mollit anim\r\nid est laborum.\r\n</page>  </plain_text>\r\n"

# full-text is read correctly

    Code
      text[["full_text"]]
    Output
      [1] "<plain_text> <page sequence=\"1\">Lorem ipsum dolor sit amet, \nconsectetur adipisici elit, sed eiusmod tempor incidunt ut labore et dolore \nmagna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris\nnisi ut aliquid ex ea commodi consequat. Quis aute iure reprehenderit in \nvoluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint\nobcaecat cupiditat non proident, sunt in culpa qui officia deserunt mollit anim\nid est laborum.\n</page>  </plain_text>\n"

---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# jstor: Import and Analyse Data from Scientific Articles

**Author:** [Thomas Klebel](https://thomasklebel.eu) <br>
**License:** [GPL v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)

[![R-CMD-check](https://github.com/ropensci/jstor/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/ropensci/jstor/actions/workflows/check-standard.yaml)
[![AppVeyorBuild status](https://ci.appveyor.com/api/projects/status/sry2gtwam7qyfw6l?svg=true)](https://ci.appveyor.com/project/tklebel/jstor)
[![Coverage status](https://codecov.io/gh/ropensci/jstor/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/jstor?branch=master)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN status](http://www.r-pkg.org/badges/version/jstor)](https://cran.r-project.org/package=jstor)
[![CRAN\_Download\_Badge](http://cranlogs.r-pkg.org/badges/grand-total/jstor)](https://CRAN.R-project.org/package=jstor)
[![rOpenSci badge](https://badges.ropensci.org/189_status.svg)](https://github.com/ropensci/onboarding/issues/189)
[![JOSS badge](http://joss.theoj.org/papers/ba29665c4bff35c37c0ef68cfe356e44/status.svg)](http://joss.theoj.org/papers/ba29665c4bff35c37c0ef68cfe356e44)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1169861.svg)](https://doi.org/10.5281/zenodo.1169861)
  
The tool [Data for Research (DfR)](http://www.jstor.org/dfr/) by JSTOR is a
valuable source for citation analysis and text mining. `jstor`
provides functions and suggests workflows for importing
datasets from DfR. It was developed to deal with very large datasets which
require an agreement, but can be used with smaller ones as well.

**Note**: As of 2021, JSTOR has moved changed the way they provide data to a new 
platform called [Constellate](https://constellate.org/). The package `jstor` has
not been adapted to this change, and might therefore only be used for legacy 
data that was optained from the old DfR platform.

The most important set of functions is a group of `jst_get_*` functions:

- `jst_get_article`
- `jst_get_authors`
- `jst_get_references`
- `jst_get_footnotes`
- `jst_get_book`
- `jst_get_chapters`
- `jst_get_full_text`
- `jst_get_ngram`


All functions which are concerned with meta data (therefore excluding
`jst_get_full_text` and `jst_get_ngram`) operate along the same lines:

1. The file is read with `xml2::read_xml()`.
2. Content of the file is extracted via XPATH or CSS-expressions.
3. The resulting data is returned in a `tibble`.


## Installation

To install the package use: 

```{r, eval=FALSE}
install.packages("jstor")
```


You can install the development version from GitHub with:

```{r gh-installation, eval = FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/jstor")
```

## Usage
In order to use `jstor`, you first need to load it:
```{r}
library(jstor)
library(magrittr)
```

The basic usage is simple: supply one of the `jst_get_*`-functions with a path
and it will return a tibble with the extracted information.
```{r, results='asis'}
jst_get_article(jst_example("article_with_references.xml")) %>% knitr::kable()

jst_get_authors(jst_example("article_with_references.xml")) %>% knitr::kable()
```

Further explanations, especially on how to use jstor's functions for importing
many files, can be found in the vignettes.

## Getting started
In order to use `jstor`, you need some data from DfR. From the
[main page](http://www.jstor.org/dfr/) you can create a dataset by searching for
terms and restricting the search regarding time, subject and content type. After
you created an account, you can download your selection. Alternatively,
you can download 
[sample datasets](http://www.jstor.org/dfr/about/sample-datasets) with documents
from before 1923 for the US, and before 1870 for all other countries. 

## Supported Elements
In their [technical specifications](http://www.jstor.org/dfr/about/technical-specifications), 
DfR lists fields which should be reliably present in all articles and books.

The following table gives an overview, which elements are supported by `jstor`.

### Articles
|`xml`-field                       |reliably present |supported in `jstor`|
|:---------------------------------|:----------------|:-------------------|
|journal-id (type="jstor")         |x                |x                   |
|journal-id (type="publisher-id")  |x                |x                   |
|journal-id (type="doi")           |                 |x                   |
|issn                              |x                |                    |
|journal-title                     |x                |x                   |
|publisher-name                    |x                |                    |
|article-id (type="doi")           |x                |x                   |
|article-id (type="jstor")         |x                |x                   |
|article-id (type="publisher-id")  |                 |x                   |
|article-type                      |                 |x                   |
|volume                            |                 |x                   |
|issue                             |                 |x                   |
|article-categories                |x                |                    |
|article-title                     |x                |x                   |
|contrib-group                     |x                |x                   |
|pub-date                          |x                |x                   |
|fpage                             |x                |x                   |
|lpage                             |                 |x                   |
|page-range                        |                 |x                   |
|product                           |x                |                    |
|self-uri                          |x                |                    |
|kwd-group                         |x                |                    |
|custom-meta-group                 |x                |x                   |
|fn-group (footnotes)              |                 |x                   |
|ref-list (references)             |                 |x                   |



### Books
|`xml`-field                       |reliably present |supported in `jstor`|
|:---------------------------------|:----------------|:-------------------|
|book-id (type="jstor")            |x                |x                   |
|discipline                        |x                |x                   |
|call-number                       |x                |                    |
|lcsh                              |x                |                    |
|book-title                        |x                |x                   |
|book-subtitle                     |                 |x                   |
|contrib-group                     |x                |x                   |
|pub-date                          |x                |x                   |
|isbn                              |x                |x                   |
|publisher-name                    |x                |x                   |
|publisher-loc                     |x                |x                   |
|permissions                       |x                |                    |
|self-uri                          |x                |                    |
|counts                            |x                |x                   |
|custom-meta-group                 |x                |x                   |


### Book Chapters
|`xml`-field                       |reliably present |supported in `jstor`|
|:---------------------------------|:----------------|:-------------------|
|book-id (type="jstor")            |x                |x                   |
|part_id                           |x                |x                   |
|part_label                        |x                |x                   |
|part-title                        |x                |x                   |
|part-subtitle                     |                 |x                   |
|contrib-group                     |x                |x                   |
|fpage                             |x                |x                   |
|abstract                          |x                |x                   |


## Code of conduct
Please note that this project is released with a 
[Contributor Code of Conduct](CONDUCT.md).
By participating in this project you agree to abide by its terms.

## Citation
To cite `jstor`, please refer to `citation(package = "jstor")`:

```
Klebel (2018). jstor: Import and Analyse Data from Scientific Texts. Journal of 
Open Source Software, 3(28), 883, https://doi.org/10.21105/joss.00883
```

## Acknowledgements
Work on `jstor` benefited from financial support for the project "Academic
Super-Elites in Sociology and Economics" by the Austrian Science Fund (FWF), 
project number "P 29211 Einzelprojekte".

Some internal functions regarding file paths and example files were adapted from
the package `readr`.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

---
title: "Automating File Import"
author: "Thomas Klebel"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Automating File Import}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r}
library(jstor)
library(purrr)
library(stringr)
```


# Intro
The `find_*` functions from `jstor` all work on a single file. Data from DfR
however contains many single files, from up to 25,000 when using the self-service
functions, up to several hundreds of thousands of files when requesting a large
dataset via an agreement.

Currently `jstor` offers two implementations to import many files: 
`jst_import_zip()` and `jst_import()`. The first one lets you import data directly
from zip archives, the second works for file paths, so you need to unzip
the archive first. I will first introduce `jst_import_zip()` and discuss the 
approach of with `jst_import()` afterwards.

# Importing data directly from zip-archives
Unpacking and working with many files directly is unpractical for at least three
reasons: 

1. If you unzip the archive, the single files will occupy a lot more disk space
than the single archive.
2. Before you can import the files, you need to list them  via `list.files` or
`system("find...")` on UNIX. Depending on the size of your data, this can take
some time.
3. There might be different types of data in your sample (journal articles, book
chapters, etc.). You need to manage matching the paths to the appropriate 
functions, which is extra work.

Importing directly from the zip archive simplifies all those tasks with a
single function: `jst_import_zip()`. For the following demonstration, we will
use a small sample archive that comes with the package.

As a first step, we should take a look at the archive and its content. This
is made easy with `jst_preview_zip()`:

```{r}
jst_preview_zip(jst_example("pseudo_dfr.zip")) %>% knitr::kable()
```

We can see that we have a simple archive with three metadata files and one 
ngram file. Before we can use `jst_import_zip()`, we first need to think about,
what we actually want to import: all of the data, or just parts? What kind of
data do we want to extract from articles, books and pamphlets? We can specify
this via `jst_define_import()`:

```{r}
import_spec <- jst_define_import(
  article = c(jst_get_article, jst_get_authors),
  book = jst_get_book,
  ngram1 = jst_get_ngram
)
```

In this case, we want to import data on articles (standard metadata plus
information on the authors), general data on books and unigrams (ngram1). This
specification can then be used with `jst_import_zip()`:

```{r}
# set up a temporary folder for output files
tmp <- tempdir()

# extract the content and write output to disk
jst_import_zip(jst_example("pseudo_dfr.zip"),
               import_spec = import_spec,
               out_file = "my_test",
               out_path = tmp)

```

We can take a look at the files within `tmp` with `list.files()`:
```{r}
list.files(tmp, pattern = "csv")
```

As you can see, `jst_import_zip()` automatically creates file names based on
the string you supplied to `out_file` to delineate the different types of
output.

If we want to re-import the data for further analysis, we can either use 
functions like `readr::read_csv()`, or a small helper function from the package
which determines and sets the column types correctly:

```{r}
jst_re_import(
  file.path(tmp, "my_test_journal_article_jst_get_article-1.csv")
) %>% 
  knitr::kable()
```

A side note on ngrams: For larger archives, importing all ngrams can take a very
long time. It is thus advisable to only import ngrams for articles which you
want to analyse, i.e. most likely a subset of the initial request. The 
function `jst_subset_ngrams()` helps you with this (see also the section on
importing bigrams in the 
[case study](https://docs.ropensci.org/jstor/articles/analysing-n-grams.html#importing-bigrams)).

## Parallel processing
Since the above process might take a while for larger archives (files have to
be unpacked, read and parsed), there might be a benefit of executing the 
function in parallel. `jst_import_zip()` and `jst_import()` use `furrr` at their
core, therefore adding parallelism is very easy. Just add the following lines
at the beginning of your script, and the import will use all available cores:

```{r, eval=FALSE}
library(future)

plan(multisession)
```

You can find out more about futures by reading the package vignette:

```{r, eval=FALSE}
vignette("future-1-overview", package = "future")
```



# Working with file paths
The above approach of importing directly from zip archives is very convenient,
but in some cases you might want to have more control over how data is imported.
For example, if you run into problems because the output from any of the 
functions provided with the package looks corrupted, you could
want to look at the original files. Alongside this, you could unzip the archive
and work with the files directly, which I will demonstrate in the following
sections.


## Unzip containers
For simple purposes it might be sensible
to unzip to a temporary directory (with `temp()` and `unzip()`) but for my 
research I simply extracted files to an external SSD, since I a) lacked disk space,
b) needed to read them fast, and c) wanted to be able to look at specific files for
debugging.

## List files
There are many ways to generate a list of all files: 
`list.files()` or using `system()` in conjunction with `find` on unix-like systems
are common options. 

For demonstration purposes I use files contained in `jstor` which can be accessed
via `system.file`:
```{r, echo=TRUE}
example_dir <- system.file("extdata", package = "jstor")
```

## `list.files`

```{r}
file_names_listed <- list.files(path = example_dir, full.names = TRUE,
                                pattern = "*.xml")
file_names_listed
```

### `system` and `find`
```{r, eval=FALSE}
file_names <- system(paste0("cd ", example_dir, "; find . -name '*.xml' -type f"), intern = TRUE)
```

```{r, eval=FALSE}
library(stringr)

file_names_system <- file_names %>%
  str_replace("^\\.\\/", "") %>%
  str_c(example_dir, "/", .)

file_names_system
#> [1] "/Library/Frameworks/R.framework/Versions/3.4/Resources/library/jstor/extdata/sample_with_footnotes.xml" 
#> [2] "/Library/Frameworks/R.framework/Versions/3.4/Resources/library/jstor/extdata/sample_book.xml"
#> [3] "/Library/Frameworks/R.framework/Versions/3.4/Resources/library/jstor/extdata/sample_with_references.xml"
```

In this case the two approaches give the same result.
The main difference seems to be though, that `list.files` sorts the output, 
whereas `find` does not. For a large amount of files (200,000) this makes
`list.files` slower, for smaller datasets the difference shouldn't make an
impact.

## Batch import
Once the file list is generated, we can apply any of the `jst_get_*`-functions 
to the list. A good and simple way for small to moderate amounts of files is to
use `purrr::map_df()`:

```{r}
# only work with journal articles
article_paths <- file_names_listed %>% 
  keep(str_detect, "with")

article_paths %>% 
  map_df(jst_get_article) %>% 
  knitr::kable()

```

This works well if 1) there are no errors and 2) if there is only a moderate
size of files. For larger numbers of files, `jst_import()` can streamline the 
process for you. This function works very similar to `jst_import_zip()`, the 
main difference being that it needs file paths as input and can only handle
one type of output.

```{r}
jst_import(article_paths, out_file = "my_second_test", .f = jst_get_article, 
           out_path = tmp)
```

The result is again written to disk, as can be seen below:

```{r}
list.files(tmp, pattern = "csv")
```

---
title: "Known Quirks of JSTOR/DfR Data"
author: "Thomas Klebel"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{Known Quirks of JSTOR Data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Collecting all the quirks
Data from JSTOR/DfR is unlike most other data you encounter when doing text 
analysis. First and foremost, the data about articles and books come from a
wide variety of journals and publishers. The level of detail and certain formats
vary because of this. `jstor` tries to deal with this situation with two
strategies:

- try to recognise the format and read data accordingly
- if this is not possible, read data as "raw" as possible, i.e. without any 
conversions

An example for the first case are references. Four different ways how
references can be specified are known at this time, and all are imported in 
specific ways to deal this variation. There might however be other formats,
which should lead to an informative error when trying to import them via 
`jst_get_references()`.

An example for the latter case are page numbers. Most of the time, the entries
for page numbers are simply `42`, or `61`. This is as expected, and could be
parsed as integers. Sometimes, there are characters present, like `M77`. This
would pose no problem either, we could simply extract all digits via regex and
parse as character. Unfortunately, sometimes the page is specified like this:
`v75i2p84`. Extracting all digits would result in `75284`, which is wrong by
a long shot. Since there might be other ways of specifying pages, `jstor` does
not attempt to parse the pages to integers when importing. However, it
offers a set of convenience functions which deal with a few common cases
(see `jst_augment()` and below).

There are many other problems or peculiarities like this. This vignette tries to
list as many as possible, and offer solutions for dealing with them.
Unfortunately I have neither the time nor the interest to wade through all the
data which you could get from DfR in order to find all possible quirks. The
following list is thus inevitably incomplete. If you encounter new 
quirks/peculiarities, it would be greatly appreciated if you sent me an email,
or [opened an issue at GitHub](https://github.com/ropensci/jstor/issues).
I will then include your findings in future version of this vignette, so this 
vignette can
be a starting point for everybody who conducts text analysis with data from
JSTOR/DfR.


# Data augmentation
After importing data via `jst_get_article()`, there are at least two tasks you
might typically want to undertake:

- Merge different identifiers for journals into one, so you can filter journals.
- Convert pages from character into integers and calculate the total number
of pages per article.

There are four functions which help you to streamline this process:

- `jst_clean_page()` attempts to turn a character vector with pages into an integer
vector.
- `jst_add_total_pages()` adds a column with the total number of pages per
article.
- `jst_unify_journal_id()` merges different identifiers for journals into one.
- `jst_augment()` wraps the above functions for convenience.

# Known quirks

In the following sections, known issues with data from DfR are described in
greater detail.

## Page numbers
Page numbers are a mess. Besides the issues mentioned above, page numbers might
sometimes be specified as "pp. 1234-83" as in
[this article from the American Journal of 
Sociology](https://www.jstor.org/stable/10.1086/659100).
Of course, this results in `first_page = 1234` and `last_page = 83`, and the
computed number of total pages from `jst_get_total_pages()` will be
negative. There is currently no general solution for this issue.


### Calculating total pages
As outlined above, page numbers come in very different forms. Besides this
problem, there is actually another issue. Imagine you would like to quantify 
the lengths of articles. Obviously you will need information on the first
and the last page of the articles. Furthermore, the pages need to be
parsed properly: you will run into troubles if you calculate page numbers like
`75284 - 42 + 1`, in case the number was parsed badly. `jst_clean_page()` tries
to do this properly, based on a few known possibilities:

- "2" -> 2
- "A2" -> 2
- "v75i2p84" -> 84


Parsing correctly is unfortunately not enough. Things like "Errata" might come
to haunt you. For example there might be an article with `first_page = 42` and 
`last_page = 362`, which
would leave you puzzled as to if this can be true^[Although it sounds absurd,
this can actually be true. There are some articles which are 200 pages long.
Obviously, they are not standard research articles. You will need to decide
if they fall into your sample or not.]. There could be a simple explanation:
the article might start on page 42, and end on page 65, and there is furthermore
an erratum on page 362. Technically, `last_page = 362` is true then, but it
will cause problems for calculating the total number of pages. Quite often,
there is information in another column which could resolve this: `page_range`, 
which in this case would look like `42 - 65, 362`. 

A small helper to deal with those situations is `jst_get_total_pages()`. It
works for page ranges, but also for first and last pages:

```{r, message=FALSE}
library(jstor)
library(dplyr)

input <- tibble::tribble(
  ~first_page,   ~last_page,    ~page_range,
  NA_real_,      NA_real_,      NA_character_, 
  1,             10,            "1 - 10",
  1,             10,            NA_character_,
  1,             NA_real_,      NA_character_, 
  1,             NA_real_,      "1-10",
  NA_real_,       NA_real_,      "1, 5-10",
  NA_real_,       NA_real_,      "1-4, 5-10",
  NA_real_,       NA_real_,      "1-4, C5-C10"
)

input %>% 
  mutate(n_pages = jst_get_total_pages(first_page, last_page, page_range))
```
This is actually identical to using `jst_add_total_pages()`:

```{r}
input %>% 
  jst_add_total_pages()
```


## Journal identifiers
Identifiers for the journal usually appear in three columns:

- `journal_doi`
- `journal_jcode`
- `journal_pub_id`

```{r, results='asis'}
sample_article <- jst_get_article(jst_example("article_with_references.xml")) 

knitr::kable(sample_article)
```

From my samples, it seems that the information in `journal_pub_id` is 
often missing, as is journal_doi. The most important identifier is thus
`journal_jcode`. In cases where both `journal_jcode` and `journal_pub_id` are
present, at least in my samples, the format of `journal_jcode` was different.
For consistency, `jst_unify_journal_id()` thus takes content of
`journal_pub_id` if it is present, and that of `journal_jcode` otherwise.

With this algorithm, it should be possible to reliably match them to 
general information about the respective journals, which are available from
`jst_get_journal_overview()`:

```{r}
sample_article %>% 
  jst_unify_journal_id() %>% 
  left_join(jst_get_journal_overview()) %>% 
  tidyr::gather(variable, value) %>% 
  knitr::kable()
```



## Duplicated ngrams
|Source                            |time span          |Part              |
|:---------------------------------|-----------------:|-----------------:|
|American Journal of Sociology     |Unknown           |Book Reviews      |

For the AJS, ngrams for book reviews are calculated per issue. There are 
numerous reviews per issue, and each of them has an identical file of ngrams,
containing ngrams for all book reviews of this issue. 

A possible strategy for dealing with this is either not to use those ngrams,
since they are calculated on all reviews in the issue, irrespective of 
whether actually all reviews of the given issue are in the sample or not.
Alternatively, one could group by issues, and only take one set of ngrams
per issue.

## Language codes
Information on langues is not consistent. For the sample article, `language` is
`eng`. 

```{r}
sample_article %>% 
  pull(language)
```

In other cases it might be `en`. It is thus advisable to take a quick look at
different variants via `distinct(meta_data, language)` or 
`count(meta_data, language)`. 


## Incorrect/odd references
When analysing data about references and footnotes, you will encounter many
inconsistencies and errors. Most of them are not due to errors from DfR, but 
stem simply from the fact, that humans make mistakes when creating manuscripts,
and not all errors with references are caught before printing.

### Problems with non-english characters
A common problem are names with non-english characters like german umlauts
(Ferdinand Tönnies) or nordic names (Gøsta Esping-Andersen). These will appear
in many different variations: Tonnies, Tönnies, Gosta, Gösta, etc.

### OCR-Issues
For older articles, you might encounter issues that stem from digitising text 
with OCR-software. A common problem is distinguishing `I` from `l`, like in
the phrase "In love". Depending on which names appear in your data, this might
lead to inconsistencies.

### Errors by article authors
There are many examples where authors make mistakes and your summary statistics
end up being skewed. 
[This article](https://www.jstor.org/stable/25074331?seq=11&refreqid=excelsior%3A6f9520d8aecff945ab2033fa66d3438e#page_scan_tab_contents) 
about "Ethics Education in the Workplace" cites the same items multiple times,
which is possibly an artifact. The advantage of using JSTOR/DfR data is, that 
you can inspect all sources and
check, if a specific pattern you see is an artifact or genuine.
---
title: Introduction to jstor
author: "Thomas Klebel"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
      toc: true
vignette: >
  %\VignetteIndexEntry{Introduction to jstor}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The tool [Data for Research (DfR)](https://www.jstor.org/dfr/) by JSTOR is a
valuable source for citation analysis and text mining. `jstor`
provides functions and suggests workflows for importing
datasets from DfR.

When using DfR, requests for datasets can
be made for small excerpts (max. 25,000 records) or large ones, which require an
agreement between the researcher and JSTOR. `jstor` was developed to deal with
very large datasets which require an agreement, but can be used with smaller
ones as well.

The most important set of functions is a group of `jst_get_*` functions:

- `jst_get_article` (for journal documents and pamphlets)
- `jst_get_authors` (for all sources)
- `jst_get_references` (for journal documents)
- `jst_get_footnotes` (for journal documents)
- `jst_get_book` (for books and research reports)
- `jst_get_chapters` (for books and possibly research reports)

I will demonstrate their usage using the
[sample dataset](https://www.jstor.org/dfr/about/sample-datasets)
which is provided by JSTOR on their website.

# General Concept
All functions from the `jst_get_*` family which are concerned with meta 
data operate along the same lines:

1. The file is read with `xml2::read_xml()`.
2. Content of the file is extracted via XPATH or CSS-expressions.
3. The resulting data is returned in a tidy `tibble`.

The functions are similar in that all operate on single files
(article, book, research report or pamphlet). Depending on the content of the
file, the output of the functions might have one or multiple rows.
`jst_get_article` always returns a `tibble` with one row: the core meta data
(like title, id, or first page of the article) are single items,
and only one article is processed at a time. 
Running `jst_get_authors` for the same article might give you a tibble with one
or multiple rows, depending on the number of authors the article has. The same
is true for `jst_get_references` and `jst_get_footnotes`. If a file has no data on
references (they might still exist, but JSTOR might not have parsed them), the
output is only one row, with missing references. If there is data on references,
each entry gets its own row. Note however, that the number of rows does not
equal the number of references. References usually start with a title like 
"References", which is obviously not a reference to another article.
Be sure to think carefully about your assumptions and to check the content of 
your data before you make inferences.

Books work a bit differently. Searching for data on
https://www.jstor.org/dfr/results lets you filter for books, which are actually
book chapters. If you receive data from DfR on a book chapter, you always get 
one xml-file with the whole book, including data on all chapters. Ngram or 
full-text data for the same entry however is processed only from single
chapters^[See the
[technical specifications](https://www.jstor.org/dfr/about/technical-specifications)
for more detail.]. Thus, the output of `jst_get_book` for a single file is
similar to the one from `jst_get_article`: it is one row with general data about
the book. `jst_get_chapters` gives you data on all chapters, and the resulting 
tibble therefore might have multiple rows.


The following sections showcase the different functions separately.

# Application
Apart from `jstor` we only need to load `dplyr` for matching records and `knitr`
for printing nice tables.

```{r, message=FALSE, warning=FALSE}
library(jstor)
library(dplyr)
library(knitr)
```



## jst_get_article
The basic usage of the `jst_get_*` functions is very simple. They take only one 
argument, the path to the file to import:
```{r}
meta_data <- jst_get_article(file_path = jst_example("article_with_references.xml"))
```

The resulting object is a `tibble` with one row and 17 columns. The columns
correspond to most of the elements documented here: https://www.jstor.org/dfr/about/technical-specifications.

The columns are:

- file_name *(chr)*: The file name of the original .xml-file. Can be used 
 for joining with other parts (authors, references, footnotes, full-texts).
- journal_doi *(chr)*: A registered identifier for the journal.
- journal_jcode *(chr)*: A identifier for the journal like "amerjsoci" for
 the "American Journal of Sociology".
- journal_pub_id *(chr)*: Similar to journal_jcode. Most of the time either
 one is present.
- article_doi *(chr)*: A registered unique identifier for the article.
- article_jcode *(chr)*: A unique identifier for the article (not a DOI).
- article_pub_id *(chr)*: Infrequent, either part of the DOI or the 
 article_jcode.
- article_type *(chr)*: The type of article (research-article, book-review,
 etc.).
- article_title *(chr)*: The title of the article.
- volume *(chr)*: The volume the article was published in.
- issue *(chr)*: The issue the article was published in.
- language *(chr)*: The language of the article.
- pub_day *(chr)*: Publication day, if specified.
- pub_month *(chr)*: Publication month, if specified.
- pub_year *(int)*: Year of publication.
- first_page *(int)*: Page number for the first page of the article.
- last_page *(int)*: Page number for the last page of the article.

Since the output from all functions are tibbles, the result is nicely formatted:

```{r, results='asis'}
meta_data %>% kable()
```

## jst_get_authors
Extracting the authors works in similar fashion:

```{r, results='asis'}
authors <- jst_get_authors(jst_example("article_with_references.xml"))
kable(authors)
```

Here we have the following columns:

- *file_name*: The same as above, used for matching articles.
- *prefix*: A prefix to the name.
- *given_name*: The given name of the author (i.e. `Albert` or `A.`).
- *surname*: The surname of the author (i.e. `Einstein`).
- *string_name*: Sometimes instead of given_name and surname, only a full string is
supplied, i.e.: `Albert Einstein`, or `Einstein, Albert`.
- *suffix*: A suffix to the name, as in `Albert Einstein, II.`.
- *author_number*: An integer representing the order of how the authors appeared in the data.

The number of rows matches the number of authors -- each author get its' own row.

## jst_get_references
```{r}
references <- jst_get_references(jst_example("article_with_references.xml"))

# # we need to remove line breaks for knitr::kable() to work properly for printing
references <- references %>%
  mutate(ref_unparsed = stringr::str_remove_all(ref_unparsed, "\\\n"))
```

We have two columns:

- *file_name*: Identifier, can be used for matching.
- *ref_title*: The title of the references sections.
- *ref_authors*: A string of authors. Several authors are separated with `;`.
- *ref_editors*: A string of editors, if present.
- *ref_collab*: A field that may contain information on the authors, if authors
            are not available.
- *ref_item_title*: The title of the cited entry.
- *ref_year*: A year, often the article's publication year, but not always. 
- *ref_source*: The source of the cited entry. For books often the title of the book,
            for articles the publisher of the journal.
- *ref_volume*: The volume of the journal article.
- *ref_first_page*: The first page of the article/chapter.
- *ref_last_page*: The last page of the article/chapter.
- *ref_publisher*: For books the publisher, for articles often missing.
- *ref_publication_type*: Known types: book, journal, web, other.
- *ref_unparsed*: The full references entry in unparsed form.


Here I display 5 random entries:

```{r, echo=FALSE}
set.seed(1234)
```

```{r}
references %>% 
  sample_n(5) %>% 
  kable()
```

This example shows several things: `file_name` is identical among rows, since
it identifies the article and all references came from one article. The the
sample file doesn't follow a typical convention (it was published in 1922), 
therefore there are several different headings (`ref_title`). Usually, this is
only "Bibliography" or "References".

Since the references were not parsed by JSTOR, we only get an unparsed version.
In general, the content of references (`unparsed_refs`) is in quite a raw state, 
quite often the result of digitising scans via OCR. For example, the last entry 
reads like this: 
`MACHADO, A.1911 Zytologische Untersuchungen fiber Trypanosoma rotatorium ...`.
There is an error here: `fiber` should be `über`. The language of the source
is German, but the OCR-software assumed English. Therefore, it didn't
recognize the *Umlaut*. Similar errors are common for text read via OCR.

For other files, we can set `parse_refs = TRUE`, so references will be imported 
in their parsed form, whenever they are available. 

```{r}
jst_get_references(
  jst_example("parsed_references.xml"),
  parse_refs = TRUE
) %>% 
  kable()
```


Note, that there might be other content present like endnotes, in case the
article used endnotes rather than footnotes.


## jst_get_footnotes
```{r, results='asis'}
jst_get_footnotes(jst_example("article_with_references.xml")) %>% 
  kable()
```

Very commonly, articles either have footnotes or references. The sample file
used here does not have footnotes, therefore a simple `tibble` with missing
footnotes is returned.

I will use another file to demonstrate footnotes.

```{r, results='asis'}
footnotes <- jst_get_footnotes(jst_example("article_with_footnotes.xml"))

footnotes %>% 
  mutate(footnotes = stringr::str_remove_all(footnotes, "\\\n")) %>% 
  kable()

```

In general, you might need to combine `jst_get_footnotes()` with
`jst_get_references()` to get all available information on citation data.

## jst_get_full_text
The function to extract full texts can't be demonstrated with proper data, since
the full texts are only supplied upon special request with DfR. The function
guesses the encoding of the specified file via `readr::guess_encoding()`, reads
the whole file and returns a `tibble` with `file_name`, `full_text` and
`encoding`.

I created a file that looks similar to files supplied by DfR with sample text:

```{r, results='asis'}
full_text <- jst_get_full_text(jst_example("full_text.txt"))
full_text %>% 
  mutate(full_text = stringr::str_remove_all(full_text, "\\\n")) %>% 
  kable()
```


## Combining results
Different parts of meta-data can be combined by using `dplyr::left_join()`.

### Matching with authors

```{r, results='asis'}
meta_data %>% 
  left_join(authors) %>%
  select(file_name, article_title, pub_year, given_name, surname) %>% 
  kable()
```

### Matching with references

```{r, results='asis'}
meta_data %>% 
  left_join(references) %>% 
  select(file_name, article_title, volume, pub_year, ref_unparsed) %>%
  head(5) %>% 
  kable()
```


# Books
Quite recently DfR added book chapters to their stack. To import metadata about
the books and chapters, jstor supplies `jst_get_book` and `jst_get_chapters`.

`jst_get_book` is very similar to `jst_get_article`. We obtain general information
about the complete book:

```{r, results='asis'}
jst_get_book(jst_example("book.xml")) %>% knitr::kable()
```

A single book might contain many chapters. `jst_get_chapters` extracts all of them.
Due to this, the function is a bit slower than most of jstor's other functions.

```{r}
chapters <- jst_get_chapters(jst_example("book.xml"))

str(chapters)
```

Without the abstracts (they are rather long) the first 10 chapters look like
this:

```{r, results='asis'}
chapters %>% 
  select(-abstract) %>% 
  head(10) %>% 
  kable()
```



Since extracting all authors for all chapters needs considerably
more time, by default authors are not extracted. You can import them like so:

```{r}
author_chap <- jst_get_chapters(jst_example("book.xml"), authors = TRUE) 
```

The authors are supplied in a list column:
```{r}
class(author_chap$authors)
```

You can expand this list with `tidyr::unnest`:

```{r}
author_chap %>% 
  tidyr::unnest(authors) %>% 
  select(part_id, given_name, surname) %>% 
  head(10) %>% 
  kable()
```

You can learn more about the concept of list-columns in Hadley Wickham's book
[R for Data Science](https://r4ds.had.co.nz/many-models.html).
---
title: Analysing n-grams with `jstor` for R
author: "Thomas Klebel"
date: "2019-12-22"
output:
  rmarkdown::html_vignette:
      toc: true
---




The service [DfR](http://www.jstor.org/dfr/) by JSTOR offers several ways for
text analysis of scientific articles. 
In this vignette I will demonstrate how to analyse n-grams which DfR delivers.

Let's suppose, we are interested in the topic of "inequality" within the 
discipline of sociology. Social inequality can be considered a prime subject of
sociological inquiry. In order to gain some context on the subject, we might be
interested to analyse frequently occurring terms. DfR offers different grades of
tokenization: n-grams for 1-3 words, e.g. unigrams, bigrams and trigrams. In
case you are unfamiliar with the analysis of tokenized text, you could read
the first few paragraphs of chapters 1 and 4 in https://www.tidytextmining.com
as an introduction.

Our analysis starts at the main page of [DfR](http://www.jstor.org/dfr/).
We create a dataset^[An introduction on how to create datasets can be found
on the page of DfR: http://www.jstor.org/dfr/about/creating-datasets]
by searching for
"inequality" and selecting "sociology" as our subject. To trim down the number
of articles, we only select articles from 1997 to 2017. After
logging in/creating an account, we select unigrams and bigrams. The resulting 
`.zip`-file can be downloaded from the email which DfR sends.


Up-front, we need to load some packages. `jstor` is  available from
CRAN, so it can be installed via `install.packages("jstor")`.



```r
library(jstor)
library(tidyverse)
library(future)

# set a lighter theme for plots
theme_set(theme_bw())

# import files in parallel
plan(multisession)
```


With the latest version of `jstor` we can now directly import files from a zip
archive. We only need to specify, where the zip archive is located, which
parts we want to extract, and where the resulting files should be saved. 

Before importing data, we can take a quick look at the contents of each
archive with `jst_preview_zip()`:


```r
jst_preview_zip("part1.zip")

## # A tibble: 4 x 3
##   type     meta_type           n
##   <chr>    <chr>           <int>
## 1 metadata journal_article 19782
## 2 ngram1   ngram1          19759
## 3 ngram2   ngram2          19739
## 4 ngram3   ngram3          19727
```


```r
jst_preview_zip("part2.zip")

## # A tibble: 4 x 3
##   type     meta_type           n
##   <chr>    <chr>           <int>
## 1 metadata journal_article  4127
## 2 ngram1   ngram1           4151
## 3 ngram2   ngram2           4172
## 4 ngram3   ngram3           4184
```

Although we could import all ngrams at this point, this would be extremely 
inefficient. It is thus best first to decide, which articles we want to analyze,
and import the corresponding ngram-files afterwards.

The following code assumes that you follow a workflow
organised around projects within RStudio (refer to
http://r4ds.had.co.nz/workflow-projects.html for further information). 


```r
import_spec <- jst_define_import(article = jst_get_article)

jst_import_zip("part1.zip", import_spec = import_spec, out_file = "part1")
jst_import_zip("part2.zip", import_spec = import_spec, out_file = "part2")

#> Processing files for journal_article with functions jst_get_article
#> Processing chunk 1/1
#>  Progress: ───────────────────────────────────────────────────────────── 100%
```

Since `jst_import_zip` writes the results to disk, we need to read the metadata
from the newly created file. This is made easy by `jst_re_import` which ensures,
that the data are read with the right column types.


```r
imported_metadata <- c("part1_journal_article_jst_get_article-1.csv",
                       "part2_journal_article_jst_get_article-1.csv") %>% 
  map_df(jst_re_import)
imported_metadata
```

```
## # A tibble: 23,909 x 19
##    file_name journal_doi journal_jcode journal_pub_id journal_title
##    <chr>     <chr>       <chr>         <chr>          <chr>        
##  1 journal-… <NA>        <NA>          amerjsoci      American Jou…
##  2 journal-… <NA>        <NA>          amerjsoci      American Jou…
##  3 journal-… <NA>        <NA>          amerjsoci      American Jou…
##  4 journal-… <NA>        <NA>          amerjsoci      American Jou…
##  5 journal-… <NA>        <NA>          amerjsoci      American Jou…
##  6 journal-… <NA>        <NA>          amerjsoci      American Jou…
##  7 journal-… <NA>        <NA>          amerjsoci      American Jou…
##  8 journal-… <NA>        <NA>          amerjsoci      American Jou…
##  9 journal-… <NA>        <NA>          amerjsoci      American Jou…
## 10 journal-… <NA>        <NA>          amerjsoci      American Jou…
## # … with 23,899 more rows, and 14 more variables: article_doi <chr>,
## #   article_pub_id <chr>, article_jcode <chr>, article_type <chr>,
## #   article_title <chr>, volume <chr>, issue <chr>, language <chr>,
## #   pub_day <chr>, pub_month <chr>, pub_year <int>, first_page <chr>,
## #   last_page <chr>, page_range <chr>
```

# Cleaning the data
Data from DfR is inherently messy. To fix a few common issues, we can use 
`jst_augment()`:


```r
imported_metadata <- jst_augment(imported_metadata, quietly = TRUE)
```

For more information on common quirks with data from DfR and how to deal with
them, take a look at the `vignette("known-quirks")`.


# Exploration
Before diving into the analysis of n-grams, we might wish to take an exploratory
look at our metadata. The first thing to look at are the types of articles.

```r
ggplot(imported_metadata, aes(article_type)) +
  geom_bar() +
  coord_flip()
```

![plot of chunk article-types](figure/article-types-1.png)

We can see, that the majority of articles are proper "research-articles", which
together with book-reviews and miscellaneous articles amount to ~99% of all
articles.


```r
imported_metadata %>% 
  count(article_type, sort = T) %>% 
  mutate(perc = scales::percent(n/sum(n)))
```

```
## # A tibble: 14 x 3
##    article_type         n perc 
##    <chr>            <int> <chr>
##  1 research-article 16289 68.1%
##  2 book-review       4552 19.0%
##  3 misc              2850 11.9%
##  4 other               89 0.4% 
##  5 in-brief            38 0.2% 
##  6 discussion          31 0.1% 
##  7 review-article      25 0.1% 
##  8 announcement        12 0.1% 
##  9 index                9 0.0% 
## 10 editorial            7 0.0% 
## 11 introduction         3 0.0% 
## 12 news                 2 0.0% 
## 13 bibliography         1 0.0% 
## 14 letter               1 0.0%
```

We must be cautious, however, when using this variable to distinguish articles
into categories. In this instance, we have "research-articles" which are
actually book-reviews:


```r
imported_metadata %>% 
  filter(article_type == "research-article" & str_detect(article_title, "Book")) %>% 
  select(file_name, article_title, pub_year)
```

```
## # A tibble: 190 x 3
##    file_name                      article_title pub_year
##    <chr>                          <chr>            <int>
##  1 journal-article-10.1086_210272 Book Reviews      1999
##  2 journal-article-10.1086_210273 Book Reviews      1999
##  3 journal-article-10.1086_210274 Book Reviews      1999
##  4 journal-article-10.1086_210275 Book Reviews      1999
##  5 journal-article-10.1086_210276 Book Reviews      1999
##  6 journal-article-10.1086_210278 Book Reviews      1999
##  7 journal-article-10.1086_210279 Book Reviews      1999
##  8 journal-article-10.1086_210280 Book Reviews      1999
##  9 journal-article-10.1086_210281 Book Reviews      1999
## 10 journal-article-10.1086_210283 Book Reviews      1999
## # … with 180 more rows
```

For the current demonstration, we want to restrict the type of articles to
research articles, therefore we need to take steps to remove book reviews and other
miscellaneous articles: First, filter by `article_type`, then remove articles
where the title starts with "Book Review".


```r
research_articles <- imported_metadata %>% 
  filter(article_type == "research-article") %>% 
  filter(!str_detect(article_title, "^Book Review"))
```


## The moving wall - filtering articles by time
Since JSTOR has a [moving wall](https://support.jstor.org/hc/en-us/articles/115004879547-JSTOR-s-Moving-Wall-Archive-vs-Current-Definitions),
we should take a look at the number of articles per year in our dataset.

```r
research_articles %>% 
  ggplot(aes(pub_year)) +
  geom_bar() 
```

![plot of chunk publications-per-year](figure/publications-per-year-1.png)

From this graph we can see an increase in research articles until 2010, after
which the number of articles
first tapers off, and then drops off sharply. For this reason we should
exclude articles at least
from 2015 onward, since the sample might get quite biased toward specific
journals.


```r
without_wall <- research_articles %>% 
  filter(pub_year < 2015)
```

## Flagship journals - filtering articles by journal

Since the amount of articles is still rather large for this demonstration, we
could select only a few journals. Here, we will look at articles from two
leading journals within the discipline, "American Journal
of Sociology" and "American Sociological Review". 


Since we cleaned the identifiers for journals with `jst_augment` earlier,
we can select our two flagship-journals very easily.


```r
flagship_journals <- without_wall %>%
  filter(journal_id %in% c("amerjsoci", "amersocirevi"))
```


# Importing bigrams
> Disclaimer: Much of the following analysis was inspired by the book
"Text Mining with R" by Julia Silge and David Robinson:
https://www.tidytextmining.com

For this demonstration we will look at bigrams to find the most common pairs of
words. Until now, we were only dealing with the metadata, therefore we need a
way to link our reduced dataset to the bigram files from DfR. The directory
structure for deliveries from DfR looks something like this:

```r
receipt-id-123456-part-001
  -- metadata
    -- journal_article_foo.xml
       .
       .
       .
  -- ngram2
    -- journal_article_foo.txt
       .
       .
       .
receipt-id-123456-part-002
  -- metadata
    -- journal_article_bar.xml
       .
       .
       .
  -- ngram2
    -- journal_article_bar.txt
       .
       .
       .
```

From this structure we can see, that the file name can serve as an identifier to
match articles and n-grams, since it is similar between metadata and
n-grams. 


To make importing a subset of ngrams more convenient, we can use
`jst_subset_ngrams`. This function returns a list of "zip-locations", which
`jst_get_ngram` can read.


```r
ngram_selection <- jst_subset_ngrams(c("part1.zip", "part2.zip"), "ngram2",
                                    flagship_journals)
```

```
## Error in utils::unzip(zip_archive, list = TRUE): zip Datei 'part2.zip' kann nicht geöffnet werden
```

```r
head(ngram_selection, 2)
```

```
## Error in head(ngram_selection, 2): Objekt 'ngram_selection' nicht gefunden
```




```r
ngram_selection <- jst_subset_ngrams(c("part1.zip", "part2.zip"), "ngram2",
                                    flagship_journals)

imported_bigrams <- ngram_selection %>% 
  furrr::future_map_dfr(jst_get_ngram)
```




From the 872 articles in our two flagship journals we
now have 6,729,813 bigrams. The 
bigrams are calculated by
JSTOR for each article independently. In order to reduce the sample to the most 
common bigrams, we have two choices: either to include only terms which occur
*within each* article a given amount of times, or to include terms which occur 
*within all* articles a given amount of times. By only including terms which 
occur more than 5 times in each article, we can drastically reduce the number of
terms. However, we might miss some important ones: there might be terms which
do not occur repeatedly within articles, but are present in all of them.

For demonstration purposes we are a bit restrictive and include only those
terms, which occur at least three times per article. 


```r
top_bigrams <- imported_bigrams %>%
  filter(n >= 3)
```



## Cleaning up bigrams
When constructing n-grams, DfR uses a stop-word list, which is quite limited 
^[for more information see the 
[technical specifications](http://www.jstor.org/dfr/about/technical-specifications)
on their page]. If we would like to restrict the terms a bit further, we could
use stopwords from `tidytext`:



```r
library(tidytext)
bigrams_separated <- top_bigrams %>%
  separate(bigrams, c("word1", "word2"), sep = " ")

bigrams_filtered <- bigrams_separated %>%
  filter(!word1 %in% stop_words$word) %>%
  filter(!word2 %in% stop_words$word)
```

After removing the stopwords we need to consider the fact, that our bigrams
were created for each article on its own. In order to analyse them together,
we need to count the terms for all articles in combination.


```r
bigram_counts <- bigrams_filtered %>%
  group_by(word1, word2) %>%
  summarise(n = sum(n)) %>%
  arrange(desc(n))

bigram_counts
```

```
## # A tibble: 106,706 x 3
## # Groups:   word1 [18,152]
##    word1        word2            n
##    <chr>        <chr>        <int>
##  1 american     sociological  9593
##  2 sociological review        9198
##  3 university   press         4603
##  4 labor        market        3555
##  5 american     journal       3273
##  6 9            7             3270
##  7 7            6             3260
##  8 10           9             3230
##  9 amsmath      amsxtra       3192
## 10 begin        document      3192
## # … with 106,696 more rows
```

From the first few terms we can see, that there are still many terms which are
not very interesting for our analysis. The terms "american" and "sociological"
are simply part of the title of a journal we selected (American Sociological 
Review). To clean the terms up, we can employ different approaches. One is to
simply filter the terms we wish to exclude:


```r
bigram_counts_clean <- bigram_counts %>%
  unite(bigram, word1, word2, sep = " ") %>%
  filter(!bigram %in% c("american sociological", "sociological review",
                        "university press", "american journal",
                        "journal sociology")) %>%
  separate(bigram, c("word1", "word2"))
```

We will look at another approach after plotting our bigrams.


# Visualize relationships
When analyzing bigrams, we might want to look at the relationships between
common terms. For this we can leverage the power of
[igraph](http://igraph.org/r/) and
[ggraph](https://cran.r-project.org/web/packages/ggraph/index.html).


```r
library(igraph)
library(ggraph)
```


First, we only keep the most common terms and then convert our `data.frame` to
an `igraph`-object. ^[If you are unfamiliar with graph theory, just take a look
at Wikipedia: [Graph Theory](https://en.wikipedia.org/wiki/Graph_theory).]



```r
bigram_graph <- bigram_counts_clean %>%
  filter(n > 500) %>%
  graph_from_data_frame()

bigram_graph
```

```
## IGRAPH 755e99b DN-- 170 161 -- 
## + attr: name (v/c), n (e/n)
## + edges from 755e99b (vertex names):
##  [1] labor                 ->market   9                     ->7       
##  [3] 7                     ->6        10                    ->9       
##  [5] amsmath               ->amsxtra  begin                 ->document
##  [7] declaremathsizes      ->10       declaretextfontcommand->textcyr 
##  [9] documentclass         ->aastex   encodingdefault       ->ot2     
## [11] newcommand            ->cyr      ot1                   ->fontenc 
## [13] ot2                   ->ot1      pagestyle             ->empty   
## [15] portland              ->xspace  
## + ... omitted several edges
```

For plotting, we will use a simple plotting function, adapted from
https://www.tidytextmining.com/ngrams.html#visualizing-a-network-of-bigrams-with-ggraph.


```r
plot_bigrams <- function(igraph_df, seed = 2016) {
  set.seed(seed)
  
  a <- grid::arrow(type = "closed", length = unit(.15, "inches"))

  ggraph(igraph_df, layout = "fr") +
    geom_edge_link(aes(edge_alpha = n), show.legend = FALSE,
                   arrow = a, end_cap = circle(.07, 'inches')) +
    geom_node_point(color = "lightblue", size = 4) +
    geom_node_text(aes(label = name), repel = T) +
    theme_graph()
}
```


```r
plot_bigrams(bigram_graph)
```

![plot of chunk first-bigram-graph](figure/first-bigram-graph-1.png)
Very obvious is a group of nodes which are not relevant to the topic of
inequality. They come from LaTeX documents and somehow made their way into the
original dataset. However, since they are more common than most of the other
terms, they are quite easy to remove. We can look at the nodes/vertices of our
graph with `V(bigram_graph)`.


```r
V(bigram_graph)
```

```
## + 170/170 vertices, named, from 755e99b:
##   [1] labor                  9                      7                     
##   [4] 10                     amsmath                begin                 
##   [7] declaremathsizes       declaretextfontcommand documentclass         
##  [10] encodingdefault        newcommand             ot1                   
##  [13] ot2                    pagestyle              portland              
##  [16] renewcommand           rmdefault              sfdefault             
##  [19] textcyr                usepackage             6                     
##  [22] aastex                 amsbsy                 amsfonts              
##  [25] amssymb                amsxtra                bm                    
##  [28] cyr                    document               empty                 
## + ... omitted several vertices
```

The first node, "labor", is relevant to us, but all other nodes from 2 to at 
least 40 are clearly irrelevant. We can remove them by simple subtraction:


```r
bigram_graph_clean <- bigram_graph - 2:40
bigram_graph_clean
```

```
## IGRAPH 3272619 DN-- 131 105 -- 
## + attr: name (v/c), n (e/n)
## + edges from 3272619 (vertex names):
##  [1] labor     ->market     labor     ->force      0         ->0         
##  [4] table     ->2          1         ->1          income    ->inequality
##  [7] black     ->white      table     ->1          table     ->3         
## [10] social    ->capital    model     ->1          model     ->2         
## [13] african   ->american   0         ->1          human     ->capital   
## [16] african   ->americans  table     ->4          1         ->2         
## [19] model     ->3          racial    ->ethnic     individual->level     
## [22] civil     ->rights     cross     ->national  
## + ... omitted several edges
```

Another apparent group is a combination of "table" or "figure" with digits. This
evidently comes from tables or figures in the papers and might suggest, that the
articles in our sample quite frequently employ quantitative methods, where
figures and tables are very common. For the analysis at hand however, we might
remove them, along with a few other irrelevant terms.



```r
bigram_graph_clean <- bigram_graph_clean - c("table", "model",
                                             as.character(0:5),
                                             "xd", "rh", "landscape", "00",
                                             "figure", "review", "79",
                                             "http", "www", "000", "01")
```


After cleaning up a bit, we can take a fresh look at our bigrams.

```r
plot_bigrams(bigram_graph_clean, 234)
```

![plot of chunk second-bigram-graph](figure/second-bigram-graph-1.png)

The figure is still far from perfect ("eco" -> "nomic" should clearly be one
term), but we can begin to analyse our network.

The most frequent bigrams are now "labor market", "labor force", and "income
inequality", which are not very surprising given that most individuals in
capitalist societies need to supply their work in exchange for income. For this
reason, the labor market and its stratification is a prime subject of the 
sociological inquiry into inequality. A few further key dimensions of
sociological analysis are apparent from the
graph: gender, race/ethnicity, occupational and socioeconomic status. That we
find many terms to be associated with the term "social" seems quite likely
given the discipline's subject.

At least two surprising results should be pointed out. First, it is not evident
how the terms "ethnic" and "racial" are connected. They do not form a typical
term like "social capital", "middle class" or similar, nor could they be 
considered a dichotomy like "black" and "white" which are often included in
tables from regressions. From a theoretical point of view, they have slightly
different meanings but are frequently being used as synonyms. 
Second, there is a group of nodes around the term
"university": university -> chicago, university -> california,
harvard -> university, etc. At least two explanations seem plausible: either,
many books are being cited which are in some way associated with those 
universities ("The University of Chicago Press" is the largest university press
in the United States), or many researchers who publish in the two
flagship-journals we selected are affiliated with those four universities: 
Harvard, Chicago, Cambridge and California. At least partly the prominence of
university -> chicago -> press might be due to the fact, that it is the
publisher of the American Journal of Sociology, and therefore included in each
article by this journal.

Altogether, these findings will not be surprising to well-educated sociologists.
Almost all bigrams in the graph are common concepts or terms within the
discipline. Most importantly, very often those concepts and terms comprise two
single words, as in "social science", "social structure", "social networks".
Two approaches might be useful to examine, *how* these concepts are being used:

1. Analyzing trigrams. It could be the case, that many of the above concepts 
would show up in combination with more interesting terms, if we were to analyze
combinations of three words. 
2. Another approach would be to first analyze data for bigrams as above, 
determining the core concepts of a field. In a second step, one could search for
those core concepts in the raw data files and extract adjacent terms. This would
not be possible with the default deliveries from DfR however, since full-text
content is only available through a dedicated agreement.

# Comparison over time
Besides looking at the overall relationship of bigrams, we could be interested
in the development over time of specific terms. Here, we want to look at how
often "labor market" and "income inequality" appear from year to year.

For this, we need to join our bigrams with the metadata.



```r
time_bigrams <- top_bigrams %>% 
  left_join(flagship_journals, by = "file_name") %>% 
  select(bigrams, n, pub_year)

head(time_bigrams)
```

```
##             bigrams  n pub_year
## 1    private sector 92     1998
## 2 market transition 68     1998
## 3               3 t 38     1998
## 4              1 00 37     1998
## 5 journal sociology 34     1998
## 6      state sector 34     1998
```

Again, we need to sum up the counts, but this time grouped by year:

```r
time_bigrams <- time_bigrams %>%
  group_by(bigrams, pub_year) %>%
  summarise(n = sum(n)) %>%
  arrange(desc(n))

time_bigrams
```

```
## # A tibble: 248,725 x 3
## # Groups:   bigrams [157,266]
##    bigrams               pub_year     n
##    <chr>                    <int> <int>
##  1 0 0                       2004  1071
##  2 et al                     2014   916
##  3 women s                   2006   885
##  4 american sociological     2014   860
##  5 sociological review       2014   814
##  6 u s                       2014   793
##  7 et al                     2011   792
##  8 et al                     2013   748
##  9 et al                     2010   691
## 10 et al                     2012   687
## # … with 248,715 more rows
```

We now only keep the two terms of interest and plot them in a simple chart.


```r
# filter the terms of interest
time_comparison <- time_bigrams %>% 
  filter(bigrams == "labor market" | bigrams == "income inequality")

ggplot(time_comparison, aes(pub_year, n, colour = bigrams)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = scales::pretty_breaks(7))
```

![plot of chunk bigrams-over-time](figure/bigrams-over-time-1.png)

In this instance, the plot does not reveal trends over time-the frequency of 
the terms is fluctuating a lot but staying on a similar level. Single spikes
of term frequency for specific years (for example income inequality in 2011)
could stem from special issues being explicitly concerned with income
inequality, although a quick glance at the corresponding issues invalidates this
hypothesis.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zip.R
\name{jst_preview_zip}
\alias{jst_preview_zip}
\title{Preview content of zip files}
\usage{
jst_preview_zip(zip_archive)
}
\arguments{
\item{zip_archive}{A path to a .zip-file from DfR}
}
\value{
The function returns a tibble with three columns:
\itemize{
\item \emph{type}: metadata or some form of ngram
\item \emph{meta_type}: which type of metadata (book_chapter, journal article, ...)
\item \emph{n}: a count for each category
}
}
\description{
This function gives you a quick preview about what a .zip-file from DfR
contains.
}
\examples{
jst_preview_zip(jst_example("pseudo_dfr.zip"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/augmentations.R
\name{jst_clean_page}
\alias{jst_clean_page}
\title{Clean a character vector of pages}
\usage{
jst_clean_page(page)
}
\arguments{
\item{page}{A character vector for pages.}
}
\value{
An integer vector, cleaned and converted from the input vector.
}
\description{
This function tries to convert character vectors into integers. This function
should not be called on page ranges.
}
\examples{
jst_clean_page("2")

# anything that is not a digit gets removed
jst_clean_page("A2-")

# a weird format from the American Journal of Sociology is convered correctly
jst_clean_page("AJSv104p126")
# this is done by searching for "p", and if it is found, extracting the
# content after "p".
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_import_fun.R
\name{jst_import}
\alias{jst_import}
\alias{jst_import_zip}
\title{Wrapper for file import}
\usage{
jst_import(
  in_paths,
  out_file,
  out_path = NULL,
  .f,
  col_names = TRUE,
  n_batches = NULL,
  files_per_batch = NULL,
  show_progress = TRUE
)

jst_import_zip(
  zip_archive,
  import_spec,
  out_file,
  out_path = NULL,
  col_names = TRUE,
  n_batches = NULL,
  files_per_batch = NULL,
  show_progress = TRUE,
  rows = NULL
)
}
\arguments{
\item{in_paths}{A character vector to the \code{xml}-files which should be
imported}

\item{out_file}{Name of files to export to. Each batch gets appended by an
increasing number.}

\item{out_path}{Path to export files to (combined with filename).}

\item{.f}{Function to use for import. Can be one of \code{jst_get_article},
\code{jst_get_authors}, \code{jst_get_references}, \code{jst_get_footnotes}, \code{jst_get_book}
or \code{jst_get_chapter}.}

\item{col_names}{Should column names be written to file? Defaults to \code{TRUE}.}

\item{n_batches}{Number of batches, defaults to 1.}

\item{files_per_batch}{Number of files for each batch. Can be used instead of
n_batches, but not in conjunction.}

\item{show_progress}{Displays a progress bar for each batch, if the session
is interactive.}

\item{zip_archive}{A path to a .zip-archive from DfR}

\item{import_spec}{A specification from \link{jst_define_import}
for which parts of a .zip-archive should be imported via which functions.}

\item{rows}{Mainly used for testing, to decrease the number of files which
are imported (i.e. 1:100).}
}
\value{
Writes \code{.csv}-files to disk.
}
\description{
This function applies an import function to a list of \code{xml}-files
or a .zip-archive in case of \code{jst_import_zip} and saves
the output in batches of \code{.csv}-files to disk.
}
\details{
Along the way, we wrap three functions, which make the process of converting
many files easier:
\itemize{
\item \code{\link[purrr:safely]{purrr::safely()}}
\item \code{\link[furrr:future_map]{furrr::future_map()}}
\item \code{\link[readr:write_delim]{readr::write_csv()}}
}

When using one of the \verb{find_*} functions, there should usually be no errors.
To avoid the whole computation to fail in the unlikely event that an error
occurs, we use \code{safely()} which let's us
continue the process, and catch the error along the way.

If you have many files to import, you might benefit from executing the
function in parallel. We use futures for this to give you maximum
flexibility. By default the code is executed sequentially. If you want to
run it in parallel, simply call \code{\link[future:plan]{future::plan()}} with
\code{\link[future:multiprocess]{future::multiprocess()}} as an argument before
running \code{jst_import} or \code{jst_import_zip}.

After importing all files, they are written to disk with
\code{\link[readr:write_delim]{readr::write_csv()}}.

Since you might run out of memory when importing a large quantity of files,
you can split up the files to import  into batches. Each batch is being
treated separately, therefore for each batch multiple processes from
\code{\link[future:multiprocess]{future::multiprocess()}} are spawned, if you added this plan.
For this reason, it is not recommended to have very small batches,
as there is an overhead for starting and ending the processes. On the other
hand, the batches should not be too large, to not exceed memory limitations.
A value of 10000 to 20000 for \code{files_per_batch} should work fine on most
machines. If the session is interactive and \code{show_progress} is \code{TRUE}, a
progress bar is displayed for each batch.
}
\examples{
\dontrun{
# read from file list --------
# find all files
meta_files <- list.files(pattern = "xml", full.names = TRUE)

# import them via `jst_get_article`
jst_import(meta_files, out_file = "imported_metadata", .f = jst_get_article,
           files_per_batch = 25000)
           
# do the same, but in parallel
library(future)
plan(multiprocess)
jst_import(meta_files, out_file = "imported_metadata", .f = jst_get_article,
           files_per_batch = 25000)

# read from zip archive ------ 
# define imports
imports <- jst_define_import(article = c(jst_get_article, jst_get_authors))

# convert the files to .csv
jst_import_zip("my_archive.zip", out_file = "my_out_file", 
                 import_spec = imports)
} 
}
\seealso{
\code{\link[=jst_combine_outputs]{jst_combine_outputs()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/re-import.R
\name{jst_combine_outputs}
\alias{jst_combine_outputs}
\title{Combine outputs from converted files}
\usage{
jst_combine_outputs(
  path,
  write_to_file = TRUE,
  out_path = NULL,
  overwrite = FALSE,
  clean_up = FALSE,
  warn = TRUE
)
}
\arguments{
\item{path}{A path to a directory, containing .csv-files from
\code{\link[=jst_import]{jst_import()}} or \code{\link[=jst_import_zip]{jst_import_zip()}}, or a vector of files which are to be
imported.}

\item{write_to_file}{Should combined data be written to a file?}

\item{out_path}{A directory where to write the combined files. If no
directory is supplied and \code{write_to_file} is \code{TRUE}, the combined files are
written to \code{path}.}

\item{overwrite}{Should files be overwritten?}

\item{clean_up}{Do you want to remove the original batch files? Use with
caution.}

\item{warn}{Should warnings be raised, if the file type cannot be determined?}
}
\value{
Either writes to disk, or returns a list with all combined files.
}
\description{
\code{jst_combine_outputs()} helps you to manage the multitude of files you might
receive after running \code{\link[=jst_import]{jst_import()}} or \code{\link[=jst_import_zip]{jst_import_zip()}} with more than
one batch.
}
\details{
Splitting the output of \code{\link[=jst_import]{jst_import()}} or \code{\link[=jst_import_zip]{jst_import_zip()}} might be done
for multiple reasons, but in the end you possibly want to combine all outputs
into one file/data.frame. This function makes a few assumptions in order to
combine files:
\itemize{
\item Files with similar names (except for trailing dashes with numbers) belong
together and will be combined into one file.
\item The names of the combined files can be determined from the original files.
If you want to combine \verb{foo-1.csv} and \verb{foo-2.csv}, the combined file will
be \code{combined_foo.csv}.
\item The directory only contains files which were imported via
\code{\link[=jst_import]{jst_import()}} or \code{\link[=jst_import_zip]{jst_import_zip()}}. If the directory contains other
\code{.csv} files, you should supply a character vector with paths to only those
files, which you want to import.
}
}
\examples{
# set up a temporary directory
tmp <- tempdir()

# find multiple files
file_list <- rep(jst_example("article_with_references.xml"), 2)

# convert and write to file
jst_import(file_list, "article", out_path = tmp, .f = jst_get_article,
             n_batches = 2, show_progress = FALSE)
             
# combine outputs
jst_combine_outputs(tmp)
list.files(tmp, "csv")

\dontrun{
# Trying to combine the files again raises an error.
jst_combine_outputs(tmp)
}

# this doesn't
jst_combine_outputs(tmp, overwrite = TRUE)

# we can remove the original files too
jst_combine_outputs(tmp, overwrite = TRUE, clean_up = TRUE)
list.files(tmp, "csv")

}
\seealso{
\code{\link[=jst_re_import]{jst_re_import()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/journal_info.R
\name{jst_get_journal_overview}
\alias{jst_get_journal_overview}
\title{Get table with information on journals}
\usage{
jst_get_journal_overview(most_recent = FALSE, quiet = FALSE)
}
\arguments{
\item{most_recent}{Should the most recent version be downloaded from DfR?
(Currently disabled due to changes on the JSTOR-servers).}

\item{quiet}{Should status messages about the download be printed?}
}
\value{
A \code{tibble} with various information about journals.
}
\description{
Download most recent or display cached version of data on journals.
}
\details{
When analysing your sample of articles from DfR, it might be helpful to have
some context about the journals in your sample. This function provides a
\code{tibble} with various information like the full name of the journal, the
short version of the name (sometimes referred to as \code{JCODE}), dates on where
the first
and last (available) issues were published, etc.

The data on journals might change. Therefore this function provides two
sources of data: a cached version which gets updated with every release, and
the ability to pull the most recent version directly from DfR (this had to
be temporarily disabled.)

The cached version was updated on 2020-04-03.
}
\examples{
# use the function without arguments to get a tibble from disk
jst_get_journal_overview()

\dontrun{
# download the most recent version from DfR
jst_get_journal_overview(most_recent = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/books.R
\name{jst_get_chapters}
\alias{jst_get_chapters}
\title{Extract information on book chapters}
\usage{
jst_get_chapters(file_path, authors = FALSE)
}
\arguments{
\item{file_path}{The path to a \code{.xml}-file for a book or research report.}

\item{authors}{Extracting the authors is an expensive operation which makes
the function ~3 times slower, depending on the number of chapters and
the number of authors. Defaults to \code{FALSE}. Use \code{authors = TRUE} to
import the authors too.}
}
\value{
A \code{tibble} containing the extracted meta-data with the following
columns:
\itemize{
\item book_id \emph{(chr)}: The book id of type "jstor", which is not a registered
DOI.
\item file_name \emph{(chr)}: The filename of the original .xml-file. Can be used
for joining with other data for the same file.
\item part_id \emph{(chr)}: The id of the part.
\item part_label \emph{(chr)}: A label for the part, if specified.
\item part_title \emph{(chr)}: The title of the part.
\item part_subtitle \emph{(chr)}: The subtitle of the part, if specified.
\item authors \emph{(list)}: A list-column with information on the authors. Can be
unnested with \code{\link[tidyr:nest]{tidyr::unnest()}}. See the examples and \code{\link[=jst_get_authors]{jst_get_authors()}}.
\item abstract \emph{(chr)}: The abstract to the part.
\item part_first_page \emph{(chr)}: The page where the part begins.
}
}
\description{
\code{jst_get_chapters()} extracts meta-data from JSTOR-XML files for book
chapters.
}
\details{
Currently, \code{jst_get_chapters()} is quite a lot slower than most of the other
functions. It is roughly 10 times slower than \code{jst_get_book}, depending on
the number of chapters to extract.
}
\examples{
# extract parts without authors
jst_get_chapters(jst_example("book.xml"))

# import authors too
parts <- jst_get_chapters(jst_example("book.xml"), authors = TRUE)
parts

tidyr::unnest(parts)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngram.R
\name{jst_subset_ngrams}
\alias{jst_subset_ngrams}
\title{Define a subset of ngrams}
\usage{
jst_subset_ngrams(zip_archives, ngram_type, selection, by = file_name)
}
\arguments{
\item{zip_archives}{A character vector of one or multiple zip-files.}

\item{ngram_type}{One of \code{"ngram1"}, \code{"ngram2"} or \code{"ngram3"}}

\item{selection}{A data.frame with the articles/books which are to be
selected.}

\item{by}{A column name for matching.}
}
\value{
A list of zip-locations which can be read via \code{\link[=jst_get_ngram]{jst_get_ngram()}}.
}
\description{
This function helps in defining a subset of ngram files which should be
imported, since importing all ngrams at once can be very expensive (in
terms of cpu and memory).
}
\examples{
# create sample output
tmp <- tempdir()
jst_import_zip(jst_example("pseudo_dfr.zip"),
               import_spec = jst_define_import(book = jst_get_book),
               out_file = "test", out_path = tmp)

# re-import as our selection for which we would like to import ngrams
selection <- jst_re_import(file.path(tmp, 
                                     "test_book_chapter_jst_get_book-1.csv"))

# get location of file
zip_loc <- jst_subset_ngrams(jst_example("pseudo_dfr.zip"), "ngram1",
                             selection) 

# import ngram
jst_get_ngram(zip_loc[[1]])
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngram.R
\name{jst_get_ngram}
\alias{jst_get_ngram}
\title{Read ngram data}
\usage{
jst_get_ngram(file)
}
\arguments{
\item{file}{A path to a file or a zip location from \code{\link[=jst_subset_ngrams]{jst_subset_ngrams()}}.}
}
\value{
A \code{\link[tibble:tibble]{tibble::tibble()}} with two columns:
\itemize{
\item \emph{ngram}: the ngram term (unigram, bigram, trigram)
\item \emph{n}: an integer for the number of times the term occurred in the original
file
}
}
\description{
Read in data on ngrams via \code{\link[readr:read_delim]{readr::read_tsv()}}.
}
\details{
This function is mainly useful when it is used in together with
\link{jst_import_zip}, where you can use it to specify reading in ngram data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/footnotes.R
\name{jst_get_footnotes}
\alias{jst_get_footnotes}
\title{Extract all footnotes}
\usage{
jst_get_footnotes(file_path)
}
\arguments{
\item{file_path}{The path to the \code{.xml}-file from which footnotes should be
extracted.}
}
\value{
A \code{tibble} containing the content from \code{fn-group} (usually the
footnotes). If there were no footnotes, \code{NA_character} is returned for the
column \code{footnotes}.
}
\description{
This function extracts the content of \code{fn-group} from journal-articles.
}
\details{
The \code{fn-group} usually contains footnotes corresponding to the article.
However, since footnotes are currently not fully supported by DfR,
there is no comprehensive documentation on the different variants. \code{jstor}
therefore extracts the content of \code{fn-group} exactly as it appears in the
data. Because of this, there might be other content present than footnotes.

In order to get all available information on citation data, you might need to
combine \code{jst_get_footnotes()} with \code{jst_get_references()}.
}
\examples{
jst_get_footnotes(jst_example("article_with_footnotes.xml"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{jst_get_file_name}
\alias{jst_get_file_name}
\title{Extract the basename of a path}
\usage{
jst_get_file_name(file_path)
}
\arguments{
\item{file_path}{A path to a file}
}
\value{
A character vector, containing the basename of the file without an
extension.
}
\description{
This helper simply extracts the basename of a path and removes the extension,
e.g. \code{foo/bar.txt} is shortened to \code{bar}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/augmentations.R
\name{jst_unify_journal_id}
\alias{jst_unify_journal_id}
\title{Unify journal IDs}
\usage{
jst_unify_journal_id(meta_data, remove_cols = TRUE)
}
\arguments{
\item{meta_data}{Data which was processed via \code{\link[=jst_get_article]{jst_get_article()}}.}

\item{remove_cols}{Should the original columns be removed after unifying?}
}
\value{
A modified \code{tibble}.

A modified tibble.
}
\description{
This function is a simple wrapper to unify journal ids.
}
\details{
Date on journal ids can be found in three columns:
\code{journal_pub_id}, \code{journal_jcode} and \code{journal_doi}. From my experience,
most of the time the relevant dat ais present in \code{journal_pub_id} or
\code{journal_jcode}, with \code{journal_jcode} being to most common identifier.
This function takes the value from \code{journal_pub_id}, and if it is missing,
that from \code{journal_jcode}. \code{journal_doi} is currently disregarded.
}
\examples{
article <- jst_get_article(jst_example("article_with_references.xml"))

jst_unify_journal_id(article)


# per default, original columns with data on the journal are removed
library(dplyr)

jst_unify_journal_id(article) \%>\% 
  select(contains("journal")) \%>\% 
  names()
  
# you can keep them by setting `remove_cols` to `FALSE`
jst_unify_journal_id(article, remove_cols = FALSE) \%>\%  
  select(contains("journal")) \%>\%
  names()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example.R
\name{jst_example}
\alias{jst_example}
\title{Get path to jstor example}
\usage{
jst_example(path = NULL)
}
\arguments{
\item{path}{Name of the example file. If \code{NULL}, the example files will be
listed.}
}
\value{
Either a character vector with the names of example files (if
\code{jst_example()} is called without supplying an argument), or a character
vector indicating the path to the example file.
}
\description{
jstor includes several sample files for demonstration purposes. This helper
makes them easy to access.
}
\details{
The code for this function was adapted from the package \code{readr}.
}
\examples{
jst_example()
jst_example("article_with_references.xml") 
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/re-import.R
\name{jst_re_import}
\alias{jst_re_import}
\title{Re-import files}
\usage{
jst_re_import(file, warn = TRUE)
}
\arguments{
\item{file}{A path to a .csv file.}

\item{warn}{Should warnings be emitted, if the type of file cannot be
determined?}
}
\value{
A \code{tibble}, with the columns determined based on heuristics applied
to the input file.
}
\description{
\code{jst_re_import()} lets you re-import a file which was exported via
\code{\link[=jst_import]{jst_import()}} or \code{\link[=jst_import_zip]{jst_import_zip()}}.
}
\details{
When attempting to re-import, a heuristic is applied. If the file has column
names which match the names from any of the \verb{find_*} functions, the file
is read with the corresponding specifications. If no column names are
recognized, files are recognized based on the number of columns. Since both
references and footnotes have only two columns, the first line is inspected
for either \code{"Referenc...|Bilbio...|Endnote..."} or \code{"Footnote..."}.
In case there is still no match, the file is read with
\code{\link[readr:read_delim]{readr::read_csv()}} with \code{guess_max = 5000} and a warning is raised.
}
\seealso{
\code{\link[=jst_combine_outputs]{jst_combine_outputs()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/full_text.R
\name{jst_get_full_text}
\alias{jst_get_full_text}
\title{Import full-text}
\usage{
jst_get_full_text(filename)
}
\arguments{
\item{filename}{The path to the file.}
}
\value{
A \code{tibble}, containing the file-path as id, the full content of
the file, and the encoding which was used to read it.
}
\description{
This function imports the full_text contents of a JSTOR-article.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/article.R
\name{jst_get_article}
\alias{jst_get_article}
\title{Extract meta information for articles}
\usage{
jst_get_article(file_path)
}
\arguments{
\item{file_path}{A \code{.xml}-file for a journal-article.}
}
\value{
A \code{tibble} containing the extracted meta-data with the following
columns:
\itemize{
\item file_name \emph{(chr)}: The file_name of the original .xml-file. Can be used
for joining with other parts (authors, references, footnotes, full-texts).
\item journal_doi \emph{(chr)}: A registered identifier for the journal.
\item journal_jcode \emph{(chr)}: A identifier for the journal like "amerjsoci" for
the "American Journal of Sociology".
\item journal_pub_id \emph{(chr)}: Similar to journal_jcode. Most of the time either
one is present.
\item journal_title \emph{(chr)}: The title of the journal.
\item article_doi \emph{(chr)}: A registered unique identifier for the article.
\item article_jcode \emph{(chr)}: A unique identifier for the article (not a DOI).
\item article_pub_id \emph{(chr)}: Infrequent, either part of the DOI or the
article_jcode.
\item article_type \emph{(chr)}: The type of article (research-article, book-review,
etc.).
\item article_title \emph{(chr)}: The title of the article.
\item volume \emph{(chr)}: The volume the article was published in.
\item issue \emph{(chr)}: The issue the article was published in.
\item language \emph{(chr)}: The language of the article.
\item pub_day \emph{(chr)}: Publication day, if specified.
\item pub_month \emph{(chr)}: Publication month, if specified.
\item pub_year \emph{(int)}: Year of publication.
\item first_page \emph{(int)}: Page number for the first page of the article.
\item last_page \emph{(int)}: Page number for the last page of the article.
\item page_range \emph{(chr)}: The range of pages for the article.
}

A note about publication dates: always the first entry is being extracted,
which should correspond to the oldest date, in case there is more than one
date.
}
\description{
\code{jst_get_article()} extracts meta-data from JSTOR-XML files for journal
articles.
}
\examples{
jst_get_article(jst_example("article_with_references.xml"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/augmentations.R
\name{jst_get_total_pages}
\alias{jst_get_total_pages}
\title{Calculate total pages}
\usage{
jst_get_total_pages(first_page, last_page, page_range, quietly = FALSE)
}
\arguments{
\item{first_page}{The first page of an article (numeric).}

\item{last_page}{The last page of an article (numeric).}

\item{page_range}{The page range of an article (character).}

\item{quietly}{Sometimes page ranges contain roman numerals like \code{xiv}. These
are not recognized, return \code{NA} and raise a warning. If set to \code{TRUE}, this
warning not raised.}
}
\value{
A vector with the calculated total pages.
}
\description{
This function is a simple helper to calculate the total number of pages of
an article.
}
\details{
This function deals with four cases:
\itemize{
\item if all three arguments are missing, NA is returned.
\item if page_range is supplied, the number of pages is calculated from it.
\item if only the first page is supplied, NA is returned.
\item if first and last page are supplied, the number of pages is calculated as
\code{last_page - first_page + 1}.
}

The algorithm to parse page ranges works as follows: A typical page range is
\verb{1-10, 200} where the article starts at page 1, ends at page 10, and has an
erratum at page 200. For this case, the range is calculated as
\code{range + single_page}, as in\code{(10 - 1 + 1) + 1 = 11}. Sometimes multiple
ranges are given: \verb{1-10, 11-20}. For those cases all ranges are summed:
\code{(10 - 1 + 1) + (20 - 11 + 1) = 20}. Another specification for multiple
ranges is \code{1-10+11-20}, which is treated similarly.
}
\examples{
# calculate pages from first and last page
first_pages <- sample(30:50, 10)
last_pages <- first_pages + sample(5:20, 10)
page_ranges <- rep(NA_character_, 10)

jst_get_total_pages(first_pages, last_pages, page_ranges)

# get pages from page range
jst_get_total_pages(NA_real_, NA_real_, "51 - 70")
jst_get_total_pages(NA_real_, NA_real_, "51 - 70, 350")
jst_get_total_pages(NA_real_, NA_real_, "350, 51 - 70")
jst_get_total_pages(NA_real_, NA_real_, "51 - 70, 80-100")
jst_get_total_pages(NA_real_, NA_real_, "51-70+350")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/books.R
\name{jst_get_book}
\alias{jst_get_book}
\title{Extract meta information for books}
\usage{
jst_get_book(file_path)
}
\arguments{
\item{file_path}{A \code{.xml}-file for a book or research report.}
}
\value{
A \code{tibble} containing the extracted meta-data with the following
columns:
\itemize{
\item file_name \emph{(chr)}: The filename of the original .xml-file. Can be used
for joining with other data for the same file.
\item discipline \emph{(chr)}: The discipline from the discipline names used on JSTOR.
\item book_id \emph{(chr)}: The book id of type "jstor", which is not a registered
DOI.
\item book_title \emph{(chr)}: The title of the book.
\item book_subtitle \emph{(chr)}: The subtitle of the book.
\item pub_day \emph{(int)}: Publication day, if specified.
\item pub_month \emph{(int)}: Publication month, if specified.
\item pub_year \emph{(int)}: Year of publication.
\item isbn \emph{(chr)}: One or more entries for the book's ISBN. If two or more,
separated by \code{"; "}.
\item publisher_name \emph{(chr)}: The name of the publisher.
\item publisher_loc \emph{(chr)}: The location of the publisher.
\item n_pages \emph{(int)}: The number of pages.
\item language \emph{(chr)}: The language of the book.
}

A note about publication dates: always the first entry is being extracted,
which should correspond to the oldest date, in case there is more than one
date.
}
\description{
\code{jst_get_book()} extracts meta-data from JSTOR-XML files for book chapters.
}
\examples{
jst_get_book(jst_example("book.xml"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/augmentations.R
\name{jst_add_total_pages}
\alias{jst_add_total_pages}
\title{Add total count of pages}
\usage{
jst_add_total_pages(meta_data, quietly = FALSE)
}
\arguments{
\item{meta_data}{Data which was processed via \code{\link[=jst_get_article]{jst_get_article()}}.}

\item{quietly}{Should warnings from converting page ranges be suppressed?}
}
\value{
A \code{tibble}, as provided with in \code{meta_data}, with an additional
column on total number of pages.
}
\description{
This function adds a column with the total count of pages. It calls
\code{\link[=jst_get_total_pages]{jst_get_total_pages()}} which does the main work.
}
\seealso{
\code{\link[=jst_get_total_pages]{jst_get_total_pages()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jstor_package.R
\docType{package}
\name{jstor}
\alias{jstor}
\title{jstor: Read Data from JSTOR/DfR}
\description{
The tool \href{https://www.jstor.org/dfr/}{Data for Research (DfR)} by JSTOR is a
valuable source for citation analysis and text mining. \code{jstor}
provides functions and suggests workflows for importing
datasets from DfR.
}
\details{
Please refer to the vignettes for information on how to use the package:

\code{browseVignettes("jstor")}

If you encounter any issues or have ideas for new features, please file an
issue at \url{https://github.com/ropensci/jstor/issues}.
}
\author{
Thomas Klebel
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/references.R
\name{jst_get_references}
\alias{jst_get_references}
\title{Extract all references}
\usage{
jst_get_references(file_path, parse_refs = FALSE)
}
\arguments{
\item{file_path}{The path to the \code{.xml}-file from which references should be
extracted.}

\item{parse_refs}{Should references be parsed, if available?}
}
\value{
A \code{tibble} with the following columns:
\itemize{
\item \code{file_name}: the identifier for the article the references come from.
\item \code{ref_title}: the title of the references sections.
\item \code{ref_authors}: a string of authors. Several authors are separated with \verb{;}.
\item \code{ref_editors}: a string of editors, if available.
\item \code{ref_collab}: a field that may contain information on the authors, if authors
are not available.
\item \code{ref_item_title}: the title of the cited entry. For books this is often
empty, with the title being in \code{ref_source}.
\item \code{ref_year}: a year, often the article's publication year, but not always.
\item \code{ref_source}: the source of the cited entry. For books often the title of the
book, for articles the publisher of the journal.
\item \code{ref_volume}: the volume of the journal article.
\item \code{ref_first_page}: the first page of the article/chapter.
\item \code{ref_last_page}: the last page of the article/chapter.
\item \code{ref_publisher}: For books the publisher, for articles often missing.
\item \code{ref_publication_type}: Known types: \code{book}, \code{journal}, \code{web}, \code{other}.
\item \code{ref_unparsed}: The full references entry in unparsed form.
}
}
\description{
This function extracts the content of \code{ref-list} from the \code{xml}-file.
}
\details{
This content may contain references or endnotes, depending on how the article
used citations. Since references are currently not fully supported by DfR,
there is no comprehensive documentation on the different variants. \code{jstor}
therefore extracts the content of \code{ref-list} exactly as it appears in the
data. Because of this, there might be other content present than references.

In order to get all available information on citation data, you might need to
combine \code{jst_get_references()} with \code{jst_get_footnotes()}.

For newer \code{xml}-files, there would be the option to extract single elements
like authors, title or date of the source, but this is not yet implemented.

In general, the implementation is not as fast as \code{jst_get_article()} -
articles with many references slow the process down.
}
\examples{
jst_get_references(jst_example("article_with_references.xml"))

# import parsed references
jst_get_references(
  jst_example("parsed_references.xml"),
  parse_refs = TRUE
) 
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/augmentations.R
\name{jst_augment}
\alias{jst_augment}
\title{Clean data from DfR}
\usage{
jst_augment(meta_data, quietly = FALSE)
}
\arguments{
\item{meta_data}{Data which was processed via \code{\link[=jst_get_article]{jst_get_article()}}.}

\item{quietly}{Should warnings from converting page ranges be suppressed?}
}
\value{
A cleaned tibble.
}
\description{
This function takes data from \code{\link[=jst_get_article]{jst_get_article()}} and
applies helper functions for cleaning the data.
}
\details{
Data from DfR is inherently messy. For many examples see
\code{vignette("known-quirks", package = "jstor")}. \code{jst_augment()} is a
convenience function that tries to deal with a few common tasks to
clean the data.

For journal articles, it calls \code{\link[=jst_clean_page]{jst_clean_page()}} to convert first and last
page, \code{\link[=jst_unify_journal_id]{jst_unify_journal_id()}} and \code{\link[=jst_add_total_pages]{jst_add_total_pages()}}.
}
\seealso{
\code{\link[=jst_clean_page]{jst_clean_page()}} \code{\link[=jst_unify_journal_id]{jst_unify_journal_id()}} \code{\link[=jst_add_total_pages]{jst_add_total_pages()}}
\code{\link[=jst_get_total_pages]{jst_get_total_pages()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_spec.R
\name{jst_define_import}
\alias{jst_define_import}
\title{Define an import specification}
\usage{
jst_define_import(...)
}
\arguments{
\item{...}{Named arguments with bare function names.}
}
\value{
A specification of imports which is necessary for
\code{\link[=jst_import_zip]{jst_import_zip()}}.
}
\description{
Define which parts of a zip file should be converted via which functions.
}
\details{
The function accepts the following names: article, book, report, pamphlet,
ngram1, ngram2, ngram3.
The corresponding files from a .zip-archive will be imported via the supplied
functions.
}
\examples{
# articles will be imported via `jst_get_article()` and `jst_get_authors()`
jst_define_import(article = c(jst_get_article, jst_get_authors))

# define a specification for importing article metadata and unigrams (ngram1)
jst_define_import(article = jst_get_article,
                  ngram1 = jst_get_ngram)
                  
                  
# import all four types with one function each
jst_define_import(article = jst_get_article,
                  book = jst_get_book,
                  report = jst_get_book,
                  pamphlet = jst_get_article)
                  
# import all four types with multiple functions
jst_define_import(article = c(jst_get_article, jst_get_authors, jst_get_references),
                  book = c(jst_get_book, jst_get_chapters),
                  report = jst_get_book,
                  pamphlet = jst_get_article)

# if you want to import chapters with authors, you can use an anonymous
# function

chapters_w_authors <- function(x) jst_get_chapters(x, authors = TRUE)
jst_define_import(book = chapters_w_authors)


\dontrun{
# define imports
imports <- jst_define_import(article = c(jst_get_article, jst_get_authors))

# convert the files to .csv
jst_import_zip("my_archive.zip", out_file = "my_out_file", 
                 import_spec = imports)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/authors.R
\name{jst_get_authors}
\alias{jst_get_authors}
\title{Extract author information}
\usage{
jst_get_authors(file_path)
}
\arguments{
\item{file_path}{A \code{.xml}-file from JSTOR containing meta-data.}
}
\value{
A \code{tibble} containing the extracted authors. All empty fields are
\code{NA_character}.
}
\description{
\code{jst_get_authors()} extracts information about authors from JSTOR-XML files.
}
\details{
The function returns a \code{tibble} with the following six columns:
\itemize{
\item \emph{prefix}: in case there was a prefix to the name, like \code{"Dr."}.
\item \emph{given_name}: The author's given name, like \code{"Albert"}.
\item \emph{surname}: The author's surname like \code{"Einstein"}.
\item \emph{string_name}: In some cases data the name is not available in separate
fields, but just as a complete string: \code{"Albert Einstein"}.
\item \emph{suffix}: a suffix to the name, like \code{"Jr."}.
\item \emph{author_number}: The authors are enumerated in the order they appear in the
data.
}
}
\examples{
jst_get_authors(jst_example("article_with_references.xml"))
}
