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

