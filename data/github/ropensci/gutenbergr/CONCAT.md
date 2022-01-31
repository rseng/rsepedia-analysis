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
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->




gutenbergr: R package to search and download public domain texts from Project Gutenberg
----------------

**Authors:** [David Robinson](http://varianceexplained.org/)<br/>
**License:** [GPL-2](https://opensource.org/licenses/GPL-2.0)

<!-- badges: start -->
[![Build Status](https://travis-ci.org/ropensci/gutenbergr.svg?branch=master)](https://travis-ci.org/ropensci/gutenbergr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/gutenbergr)]( https://CRAN.R-project.org/package=gutenbergr)
[![Build status](https://ci.appveyor.com/api/projects/status/lqb7hngtj5epsmd1?svg=true)](https://ci.appveyor.com/project/ropensci/gutenbergr-dujv9)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/gutenbergr/master.svg)](https://codecov.io/github/ropensci/gutenbergr?branch=master)
[![rOpenSci peer-review](https://badges.ropensci.org/41_status.svg)](https://github.com/ropensci/software-review/issues/41)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

Download and process public domain works from the [Project Gutenberg](https://www.gutenberg.org/) collection. Includes

* A function `gutenberg_download()` that downloads one or more works from Project Gutenberg by ID: e.g., `gutenberg_download(84)` downloads the text of Frankenstein.
* Metadata for all Project Gutenberg works as R datasets, so that they can be searched and filtered:
  * `gutenberg_metadata` contains information about each work, pairing Gutenberg ID with title, author, language, etc
  * `gutenberg_authors` contains information about each author, such as aliases and birth/death year
  * `gutenberg_subjects` contains pairings of works with Library of Congress subjects and topics

### Installation

Install the package with:


```r
install.packages("gutenbergr")
```

Or install the development version using [devtools](https://github.com/r-lib/devtools) with:


```r
devtools::install_github("ropensci/gutenbergr")
```

### Examples

The `gutenberg_works()` function retrieves, by default, a table of metadata for all unique English-language Project Gutenberg works that have text associated with them. (The `gutenberg_metadata` dataset has all Gutenberg works, unfiltered).



Suppose we wanted to download Emily Bronte's "Wuthering Heights." We could find the book's ID by filtering:


```r
library(dplyr)
library(gutenbergr)

gutenberg_works() %>%
  filter(title == "Wuthering Heights")
#> # A tibble: 1 x 8
#>   gutenberg_id title             author        gutenberg_author_id language
#>          <int> <chr>             <chr>                       <int> <chr>   
#> 1          768 Wuthering Heights Brontë, Emily                 405 en      
#>   gutenberg_bookshelf                                 rights                    has_text
#>   <chr>                                               <chr>                     <lgl>   
#> 1 Gothic Fiction/Movie Books/Best Books Ever Listings Public domain in the USA. TRUE

# or just:
gutenberg_works(title == "Wuthering Heights")
#> # A tibble: 1 x 8
#>   gutenberg_id title             author        gutenberg_author_id language
#>          <int> <chr>             <chr>                       <int> <chr>   
#> 1          768 Wuthering Heights Brontë, Emily                 405 en      
#>   gutenberg_bookshelf                                 rights                    has_text
#>   <chr>                                               <chr>                     <lgl>   
#> 1 Gothic Fiction/Movie Books/Best Books Ever Listings Public domain in the USA. TRUE
```

Since we see that it has `gutenberg_id` 768, we can download it with the `gutenberg_download()` function:


```r
wuthering_heights <- gutenberg_download(768)
#>                                          768 
#> "http://aleph.gutenberg.org/7/6/768/768.zip"
wuthering_heights
#> # A tibble: 12,085 x 2
#>    gutenberg_id text                                                                     
#>           <int> <chr>                                                                    
#>  1          768 "WUTHERING HEIGHTS"                                                      
#>  2          768 ""                                                                       
#>  3          768 ""                                                                       
#>  4          768 "CHAPTER I"                                                              
#>  5          768 ""                                                                       
#>  6          768 ""                                                                       
#>  7          768 "1801.--I have just returned from a visit to my landlord--the solitary"  
#>  8          768 "neighbour that I shall be troubled with.  This is certainly a beautiful"
#>  9          768 "country!  In all England, I do not believe that I could have fixed on a"
#> 10          768 "situation so completely removed from the stir of society.  A perfect"   
#> # … with 12,075 more rows
```

`gutenberg_download` can download multiple books when given multiple IDs. It also takes a `meta_fields` argument that will add variables from the metadata.


```r
# 1260 is the ID of Jane Eyre
books <- gutenberg_download(c(768, 1260), meta_fields = "title")
#>                                              768                                             1260 
#>     "http://aleph.gutenberg.org/7/6/768/768.zip" "http://aleph.gutenberg.org/1/2/6/1260/1260.zip"
books
#> # A tibble: 32,744 x 3
#>    gutenberg_id text                                                                     
#>           <int> <chr>                                                                    
#>  1          768 "WUTHERING HEIGHTS"                                                      
#>  2          768 ""                                                                       
#>  3          768 ""                                                                       
#>  4          768 "CHAPTER I"                                                              
#>  5          768 ""                                                                       
#>  6          768 ""                                                                       
#>  7          768 "1801.--I have just returned from a visit to my landlord--the solitary"  
#>  8          768 "neighbour that I shall be troubled with.  This is certainly a beautiful"
#>  9          768 "country!  In all England, I do not believe that I could have fixed on a"
#> 10          768 "situation so completely removed from the stir of society.  A perfect"   
#>    title            
#>    <chr>            
#>  1 Wuthering Heights
#>  2 Wuthering Heights
#>  3 Wuthering Heights
#>  4 Wuthering Heights
#>  5 Wuthering Heights
#>  6 Wuthering Heights
#>  7 Wuthering Heights
#>  8 Wuthering Heights
#>  9 Wuthering Heights
#> 10 Wuthering Heights
#> # … with 32,734 more rows

books %>%
  count(title)
#> # A tibble: 2 x 2
#>   title                           n
#>   <chr>                       <int>
#> 1 Jane Eyre: An Autobiography 20659
#> 2 Wuthering Heights           12085
```

It can also take the output of `gutenberg_works` directly. For example, we could get the text of all Aristotle's works, each annotated with both `gutenberg_id` and `title`, using:


```r
aristotle_books <- gutenberg_works(author == "Aristotle") %>%
  gutenberg_download(meta_fields = "title")
#>                                                 1974 
#>     "http://aleph.gutenberg.org/1/9/7/1974/1974.zip" 
#>                                                 2412 
#>     "http://aleph.gutenberg.org/2/4/1/2412/2412.zip" 
#>                                                 6762 
#>     "http://aleph.gutenberg.org/6/7/6/6762/6762.zip" 
#>                                                 6763 
#>     "http://aleph.gutenberg.org/6/7/6/6763/6763.zip" 
#>                                                 8438 
#>     "http://aleph.gutenberg.org/8/4/3/8438/8438.zip" 
#>                                                12699 
#> "http://aleph.gutenberg.org/1/2/6/9/12699/12699.zip" 
#>                                                26095 
#> "http://aleph.gutenberg.org/2/6/0/9/26095/26095.zip"

aristotle_books
#> # A tibble: 39,950 x 3
#>    gutenberg_id text                                                                    
#>           <int> <chr>                                                                   
#>  1         1974 "THE POETICS OF ARISTOTLE"                                              
#>  2         1974 ""                                                                      
#>  3         1974 "By Aristotle"                                                          
#>  4         1974 ""                                                                      
#>  5         1974 "A Translation By S. H. Butcher"                                        
#>  6         1974 ""                                                                      
#>  7         1974 ""                                                                      
#>  8         1974 "[Transcriber's Annotations and Conventions: the translator left"       
#>  9         1974 "intact some Greek words to illustrate a specific point of the original"
#> 10         1974 "discourse. In this transcription, in order to retain the accuracy of"  
#>    title                   
#>    <chr>                   
#>  1 The Poetics of Aristotle
#>  2 The Poetics of Aristotle
#>  3 The Poetics of Aristotle
#>  4 The Poetics of Aristotle
#>  5 The Poetics of Aristotle
#>  6 The Poetics of Aristotle
#>  7 The Poetics of Aristotle
#>  8 The Poetics of Aristotle
#>  9 The Poetics of Aristotle
#> 10 The Poetics of Aristotle
#> # … with 39,940 more rows
```

### FAQ

#### What do I do with the text once I have it?

* The [Natural Language Processing CRAN View](https://CRAN.R-project.org/view=NaturalLanguageProcessing) suggests many R packages related to text mining, especially around the [tm package](https://cran.r-project.org/package=tm).
* The [tidytext](https://github.com/juliasilge/tidytext) package is useful for tokenization and analysis, especially since gutenbergr downloads books as a data frame already.
* You could match the `wikipedia` column in `gutenberg_author` to Wikipedia content with the [WikipediR](https://cran.r-project.org/package=WikipediR) package or to pageview statistics with the [wikipediatrend](https://cran.r-project.org/package=wikipediatrend) package.
* If you're considering an analysis based on author name, you may find the [humaniformat](https://cran.r-project.org/package=humaniformat) (for extraction of first names) and [gender](https://cran.r-project.org/package=gender) (prediction of gender from first names) packages useful. (Note that humaniformat has a `format_reverse` function for reversing "Last, First" names).

#### How were the metadata R files generated?

See the [data-raw](https://github.com/ropensci/gutenbergr/tree/master/data-raw) directory for the scripts that generate these datasets. As of now, these were generated from [the Project Gutenberg catalog](https://www.gutenberg.org/wiki/Gutenberg:Feeds#The_Complete_Project_Gutenberg_Catalog) on **05 May 2016**.

#### Do you respect the rules regarding robot access to Project Gutenberg?

Yes! The package respects [these rules](https://www.gutenberg.org/wiki/Gutenberg:Information_About_Robot_Access_to_our_Pages) and complies to the best of our ability. Namely:

* Project Gutenberg allows wget to harvest Project Gutenberg using [this list of links](http://www.gutenberg.org/robot/harvest?filetypes[]=html). The gutenbergr package visits that page once to find the recommended mirror for the user's location.
* We retrieve the book text directly from that mirror using links in the same format. For example, Frankenstein (book 84) is retrieved from `http://www.gutenberg.lib.md.us/8/84/84.zip`.
* We retrieve the .zip file rather than txt to minimize bandwidth on the mirror.

Still, this package is *not* the right way to download the entire Project Gutenberg corpus (or all from a particular language). For that, follow [their recommendation](https://www.gutenberg.org/wiki/Gutenberg:Information_About_Robot_Access_to_our_Pages) to use wget or set up a mirror. This package is recommended for downloading a single work, or works for a particular author or topic.

### Code of Conduct

This project is released with a [Contributor Code of Conduct](https://github.com/ropensci/gutenbergr/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
# gutenbergr 0.2.0

* Changed to comply with CRAN policies for API packages. Tests that do connect to Project Gutenberg are skipped on CRAN, and are supplemented with tests that mock the connection.
* This adds a files argument to gutenberg_download that is generally used only for testing.
* Made changes to work with dplyr 1.0.0, removing filter_ and distinct_.
* Fixed links to https

# gutenbergr 0.1.5

* Make compatible with tidyr v1.0.0
* data_frame is deprecated, use tibble (thanks @evanodell for #21)
* ROpenSci updates to README (thanks @maelle for #23)

# gutenbergr 0.1.4

* Added curl to SUGGESTS, since if it's not installed `readr::read_lines` could fail

# gutenbergr 0.1.3

* The Project Gutenberg mirror in Maryland Public Libraries (http://www.gutenberg.lib.md.us) has been broken for months. When it is provided from robot/harvest, replaces with `http://aleph.gutenberg.org`.
* Changed test of .zip capability not to run on CRAN
* Removed rvest dependency

# gutenbergr 0.1.2

* Made compatible with change to `distinct` in dplyr 0.5 (which is about to be submitted to CRAN)
* Removed xml2 dependency

# gutenbergr 0.1.1

* Transferred repo ownership to [ropenscilabs](https://github.com/ropenscilabs)
* The license was changed from MIT to GPL-2. This is based on the realization that the catalog data is licensed under the GPL, and the package includes a processed version of the catalog data. (See [here](https://www.gutenberg.org/wiki/Gutenberg:Feeds#The_Complete_Project_Gutenberg_Catalog)).
* Updated datasets to 5/5/2016 and added a "date_updated" attribute to tell when they were last updated
* Added `all_languages` and `only_languages` arguments to `gutenberg_works`, allowing fine-grained control of languages. (For example, "either English or French" or "both English and French")
* Changed get_gutenberg_mirror to use xml2 directly, in order to handle AppVeyor
* Removed use of data() in `gutenberg_works`, since it slows down `gutenberg_works` about 2X
* Various documentation, vignette, and README adjustments in response to ROpenSci feedback.
* Added AppVeyor for Windows continuous integration
* Added code coverage information through codecov.io and covr, along with tests to improve coverage

# gutenbergr 0.1

* First version of package, including
  * `gutenberg_download` function, for downloading one or more works from Project Gutenberg using Gutenberg IDs
  * Datasets of Project Gutenberg metadata: `gutenberg_metadata`, `gutenberg_subjects`, `gutenberg_authors`
  * `gutenberg_works` function to retrieve filtered version of `gutenberg_metadata`
  * Introductory vignette including basic examples of downloading books
  * Unit tests for `gutenberg_download` and `gutenberg_works`
* Added a `NEWS.md` file to track changes to the package.
# gutenbergr 0.2.0

This is a re-submission after gutenbergr 0.1.5 was archived, which complies with the CRAN policies. My sincere apologies for not fixing the issue sooner, and I hope it can be returned to CRAN.

## Changes

Major changes:

* Fixed to comply with CRAN policies for API packages. Tests that do connect to Project Gutenberg are skipped on CRAN, and are supplemented with tests that mock the connection.

Minor changes:

* Made changes to work with dplyr 1.0.0, removing filter_ and distinct_.
* Fixed links to https

## Test environments

* local OS X install, R 4.0.2
* Ubuntu 16.04.6 LTS (on travis-ci)
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 2 notes

The Notes occur on win-builder and are about examples that take slightly more than 5 seconds to run, as well as UTF strings in the data that are generally part of author names.
This set of scripts produces the three .rda files in the [data](../data) directory. It is only moderately reproducible, as it needs to be run in a particular order and in a few cases even hardcodes paths. Expect it to be updated and improved soon.

If you're interested in parsing and processing the Gutenberg metadata yourself, the only file you really need is [metadata.json.gz](metadata.json.gz) (written by [gitenberg_meta.py](gitenberg_meta.py)), which contains one line with a JSON dictionary for every Project Gutenberg work. The JSON dictionary was produced from the Gutenberg RDF using the [gitberg Python package](https://github.com/gitenberg-dev/gitberg).

You can find out the last time the metadata was updated in the package by running:

    attr(gutenberg_metadata, "date_updated")

There are no guarantees about how often the metadata will be updated in the package. If you're interested in works that have been recently added or had their metadata edited on Project Gutenberg, you may want to run the scripts yourself.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->


```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-",
  message = FALSE,
  warning = FALSE
)
```

gutenbergr: R package to search and download public domain texts from Project Gutenberg
----------------

**Authors:** [David Robinson](http://varianceexplained.org/)<br/>
**License:** [GPL-2](https://opensource.org/licenses/GPL-2.0)

<!-- badges: start -->
[![Build Status](https://travis-ci.org/ropensci/gutenbergr.svg?branch=master)](https://travis-ci.org/ropensci/gutenbergr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/gutenbergr)]( https://CRAN.R-project.org/package=gutenbergr)
[![Build status](https://ci.appveyor.com/api/projects/status/lqb7hngtj5epsmd1?svg=true)](https://ci.appveyor.com/project/ropensci/gutenbergr-dujv9)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/gutenbergr/master.svg)](https://codecov.io/github/ropensci/gutenbergr?branch=master)
[![rOpenSci peer-review](https://badges.ropensci.org/41_status.svg)](https://github.com/ropensci/software-review/issues/41)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

Download and process public domain works from the [Project Gutenberg](https://www.gutenberg.org/) collection. Includes

* A function `gutenberg_download()` that downloads one or more works from Project Gutenberg by ID: e.g., `gutenberg_download(84)` downloads the text of Frankenstein.
* Metadata for all Project Gutenberg works as R datasets, so that they can be searched and filtered:
  * `gutenberg_metadata` contains information about each work, pairing Gutenberg ID with title, author, language, etc
  * `gutenberg_authors` contains information about each author, such as aliases and birth/death year
  * `gutenberg_subjects` contains pairings of works with Library of Congress subjects and topics

### Installation

Install the package with:

```{r eval = FALSE}
install.packages("gutenbergr")
```

Or install the development version using [devtools](https://github.com/r-lib/devtools) with:

```{r eval = FALSE}
devtools::install_github("ropensci/gutenbergr")
```

### Examples

The `gutenberg_works()` function retrieves, by default, a table of metadata for all unique English-language Project Gutenberg works that have text associated with them. (The `gutenberg_metadata` dataset has all Gutenberg works, unfiltered).

```{r echo = FALSE}
options(dplyr.width = 140)
options(width = 100)
```

Suppose we wanted to download Emily Bronte's "Wuthering Heights." We could find the book's ID by filtering:

```{r}
library(dplyr)
library(gutenbergr)

gutenberg_works() %>%
  filter(title == "Wuthering Heights")

# or just:
gutenberg_works(title == "Wuthering Heights")
```

Since we see that it has `gutenberg_id` 768, we can download it with the `gutenberg_download()` function:

```{r}
wuthering_heights <- gutenberg_download(768)
wuthering_heights
```

`gutenberg_download` can download multiple books when given multiple IDs. It also takes a `meta_fields` argument that will add variables from the metadata.

```{r}
# 1260 is the ID of Jane Eyre
books <- gutenberg_download(c(768, 1260), meta_fields = "title")
books

books %>%
  count(title)
```

It can also take the output of `gutenberg_works` directly. For example, we could get the text of all Aristotle's works, each annotated with both `gutenberg_id` and `title`, using:

```{r}
aristotle_books <- gutenberg_works(author == "Aristotle") %>%
  gutenberg_download(meta_fields = "title")

aristotle_books
```

### FAQ

#### What do I do with the text once I have it?

* The [Natural Language Processing CRAN View](https://CRAN.R-project.org/view=NaturalLanguageProcessing) suggests many R packages related to text mining, especially around the [tm package](https://cran.r-project.org/package=tm).
* The [tidytext](https://github.com/juliasilge/tidytext) package is useful for tokenization and analysis, especially since gutenbergr downloads books as a data frame already.
* You could match the `wikipedia` column in `gutenberg_author` to Wikipedia content with the [WikipediR](https://cran.r-project.org/package=WikipediR) package or to pageview statistics with the [wikipediatrend](https://cran.r-project.org/package=wikipediatrend) package.
* If you're considering an analysis based on author name, you may find the [humaniformat](https://cran.r-project.org/package=humaniformat) (for extraction of first names) and [gender](https://cran.r-project.org/package=gender) (prediction of gender from first names) packages useful. (Note that humaniformat has a `format_reverse` function for reversing "Last, First" names).

#### How were the metadata R files generated?

See the [data-raw](https://github.com/ropensci/gutenbergr/tree/master/data-raw) directory for the scripts that generate these datasets. As of now, these were generated from [the Project Gutenberg catalog](https://www.gutenberg.org/wiki/Gutenberg:Feeds#The_Complete_Project_Gutenberg_Catalog) on **`r format(attr(gutenberg_metadata, "date_updated"), '%d %B %Y')`**.

#### Do you respect the rules regarding robot access to Project Gutenberg?

Yes! The package respects [these rules](https://www.gutenberg.org/wiki/Gutenberg:Information_About_Robot_Access_to_our_Pages) and complies to the best of our ability. Namely:

* Project Gutenberg allows wget to harvest Project Gutenberg using [this list of links](http://www.gutenberg.org/robot/harvest?filetypes[]=html). The gutenbergr package visits that page once to find the recommended mirror for the user's location.
* We retrieve the book text directly from that mirror using links in the same format. For example, Frankenstein (book 84) is retrieved from `http://www.gutenberg.lib.md.us/8/84/84.zip`.
* We retrieve the .zip file rather than txt to minimize bandwidth on the mirror.

Still, this package is *not* the right way to download the entire Project Gutenberg corpus (or all from a particular language). For that, follow [their recommendation](https://www.gutenberg.org/wiki/Gutenberg:Information_About_Robot_Access_to_our_Pages) to use wget or set up a mirror. This package is recommended for downloading a single work, or works for a particular author or topic.

### Code of Conduct

This project is released with a [Contributor Code of Conduct](https://github.com/ropensci/gutenbergr/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
---
title: "gutenbergr: Search and download public domain texts from Project Gutenberg"
author: "David Robinson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{gutenbergr: Search and download public domain texts from Project Gutenberg}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo = FALSE}
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
```

The gutenbergr package helps you download and process public domain works from the [Project Gutenberg](http://www.gutenberg.org/) collection. This includes both tools for downloading books (and stripping header/footer information), and a complete dataset of Project Gutenberg metadata that can be used to find words of interest. Includes:

* A function `gutenberg_download()` that downloads one or more works from Project Gutenberg by ID: e.g., `gutenberg_download(84)` downloads the text of Frankenstein.
* Metadata for all Project Gutenberg works as R datasets, so that they can be searched and filtered:
  * `gutenberg_metadata` contains information about each work, pairing Gutenberg ID with title, author, language, etc
  * `gutenberg_authors` contains information about each author, such as aliases and birth/death year
  * `gutenberg_subjects` contains pairings of works with Library of Congress subjects and topics
  
### Project Gutenberg Metadata

This package contains metadata for all Project Gutenberg works as R datasets, so that you can search and filter for particular works before downloading.

The dataset `gutenberg_metadata` contains information about each work, pairing Gutenberg ID with title, author, language, etc:

```{r}
library(gutenbergr)
gutenberg_metadata
```

For example, you could find the Gutenberg ID of Wuthering Heights by doing:

```{r}
library(dplyr)

gutenberg_metadata %>%
  filter(title == "Wuthering Heights")
```

In many analyses, you may want to filter just for English works, avoid duplicates, and include only books that have text that can be downloaded. The `gutenberg_works()` function does this pre-filtering:

```{r}
gutenberg_works()
```

It also allows you to perform filtering as an argument:

```{r}
gutenberg_works(author == "Austen, Jane")

# or with a regular expression

library(stringr)
gutenberg_works(str_detect(author, "Austen"))
```

The meta-data currently in the package was last updated on **`r format(attr(gutenberg_metadata, "date_updated"), '%d %B %Y')`**.

### Downloading books by ID

The function `gutenberg_download()` downloads one or more works from Project Gutenberg based on their ID. For example, we earlier saw that "Wuthering Heights" has ID 768 (see [the URL here](https://www.gutenberg.org/ebooks/768)), so `gutenberg_download(768)` downloads this text.

```{r}
f768 <- system.file("extdata", "768.zip", package = "gutenbergr")
wuthering_heights <- gutenberg_download(768,
                                        files = f768,
                                        mirror = "http://aleph.gutenberg.org")
```


```{r eval = FALSE}
wuthering_heights <- gutenberg_download(768)
```

```{r}
wuthering_heights
```

Notice it is returned as a tbl_df (a type of data frame) including two variables: `gutenberg_id` (useful if multiple books are returned), and a character vector of the text, one row per line. Notice that the header and footer added by Project Gutenberg (visible [here](http://www.gutenberg.org/files/768/768.txt)) have been stripped away.

Provide a vector of IDs to download multiple books. For example, to download Jane Eyre (book [1260](https://www.gutenberg.org/ebooks/1260)) along with Wuthering Heights, do:

```{r}
f1260 <- system.file("extdata", "1260.zip", package = "gutenbergr")
books <- gutenberg_download(c(768, 1260),
                            meta_fields = "title",
                            files = c(f768, f1260),
                            mirror = "http://aleph.gutenberg.org")
```


```{r, eval = FALSE}
books <- gutenberg_download(c(768, 1260), meta_fields = "title")
```

```{r}
books
```

Notice that the `meta_fields` argument allows us to add one or more additional fields from the `gutenberg_metadata` to the downloaded text, such as title or author.

```{r}
books %>%
  count(title)
```

### Other meta-datasets

You may want to select books based on information other than their title or author, such as their genre or topic. `gutenberg_subjects` contains pairings of works with Library of Congress subjects and topics. "lcc" means [Library of Congress Classification](https://www.loc.gov/catdir/cpso/lcco/), while "lcsh" means [Library of Congress subject headings](https://id.loc.gov/authorities/subjects.html):

```{r}
gutenberg_subjects
```

This is useful for extracting texts from a particular topic or genre, such as detective stories, or a particular character, such as Sherlock Holmes. The `gutenberg_id` column can then be used to download these texts or to link with other metadata.

```{r}
gutenberg_subjects %>%
  filter(subject == "Detective and mystery stories")

gutenberg_subjects %>%
  filter(grepl("Holmes, Sherlock", subject))
```

`gutenberg_authors` contains information about each author, such as aliases and birth/death year:

```{r}
gutenberg_authors
```

### Analysis

What's next after retrieving a book's text? Well, having the book as a data frame is especially useful for working with the [tidytext](https://github.com/juliasilge/tidytext) package for text analysis.

```{r}
library(tidytext)

words <- books %>%
  unnest_tokens(word, text)

words

word_counts <- words %>%
  anti_join(stop_words, by = "word") %>%
  count(title, word, sort = TRUE)

word_counts
```

You may also find these resources useful:

* The [Natural Language Processing CRAN View](https://CRAN.R-project.org/view=NaturalLanguageProcessing) suggests many R packages related to text mining, especially around the [tm package](https://cran.r-project.org/package=tm)
* You could match the `wikipedia` column in `gutenberg_author` to Wikipedia content with the [WikipediR](https://cran.r-project.org/package=WikipediR) package or to pageview statistics with the [wikipediatrend](https://cran.r-project.org/package=wikipediatrend) package
* If you're considering an analysis based on author name, you may find the [humaniformat](https://cran.r-project.org/package=humaniformat) (for extraction of first names) and [gender](https://cran.r-project.org/package=gender) (prediction of gender from first names) packages useful. (Note that humaniformat has a `format_reverse` function for reversing "Last, First" names).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gutenberg_download.R
\name{gutenberg_strip}
\alias{gutenberg_strip}
\title{Strip header and footer content from a Project Gutenberg book}
\usage{
gutenberg_strip(text)
}
\arguments{
\item{text}{A character vector with lines of a book}
}
\description{
Strip header and footer content from a Project Gutenberg book. This
is based on some formatting guesses so it may not be perfect. It
will also not strip tables of contents, prologues, or other text
that appears at the start of a book.
}
\examples{

library(dplyr)
book <- gutenberg_works(title == "Pride and Prejudice") \%>\%
  gutenberg_download(strip = FALSE)

head(book$text, 10)
tail(book$text, 10)

text_stripped <- gutenberg_strip(book$text)

head(text_stripped, 10)
tail(text_stripped, 10)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gutenberg_works.R
\name{gutenberg_works}
\alias{gutenberg_works}
\title{Get a filtered table of Gutenberg work metadata}
\usage{
gutenberg_works(
  ...,
  languages = "en",
  only_text = TRUE,
  rights = c("Public domain in the USA.", "None"),
  distinct = TRUE,
  all_languages = FALSE,
  only_languages = TRUE
)
}
\arguments{
\item{...}{Additional filters, given as expressions using the variables
in the \link{gutenberg_metadata} dataset (e.g. \code{author == "Austen, Jane"})}

\item{languages}{Vector of languages to include}

\item{only_text}{Whether the works must have Gutenberg text attached. Works
without text (e.g. audiobooks) cannot be downloaded with
\code{\link{gutenberg_download}}}

\item{rights}{Values to allow in the \code{rights} field. By default allows
public domain in the US or "None", while excluding works under copyright.
NULL allows any value of Rights}

\item{distinct}{Whether to return only one distinct combination of each
title and gutenberg_author_id. If multiple occur (that fulfill the other
conditions), it uses the one with the lowest ID}

\item{all_languages}{Whether, if multiple languages are given, all of them
need to be present in a work. For example, if \code{c("en", "fr")} are given,
whether only \code{en/fr} as opposed to English or French works should be
returned}

\item{only_languages}{Whether to exclude works that have other languages
besides the ones provided. For example, whether to include \code{en/fr}
when English works are requested}
}
\value{
A tbl_df (see the tibble or dplyr packages) with one row for
each work, in the same format as \link{gutenberg_metadata}.
}
\description{
Get a table of Gutenberg work metadata that has been filtered by some common
(settable) defaults, along with the option to add additional filters.
This function is for convenience when working with common conditions
when pulling a set of books to analyze.
For more detailed filtering of the entire Project Gutenberg
metadata, use the \link{gutenberg_metadata} and related datasets.
}
\details{
By default, returns

\itemize{
  \item{English-language works}
  \item{That are in text format in Gutenberg (as opposed to audio)}
  \item{Whose text is not under copyright}
  \item{At most one distinct field for each title/author pair}
}
}
\examples{

library(dplyr)

gutenberg_works()

# filter conditions
gutenberg_works(author == "Shakespeare, William")

# language specifications

gutenberg_works(languages = "es") \%>\%
  count(language, sort = TRUE)

gutenberg_works(languages = c("en", "es")) \%>\%
  count(language, sort = TRUE)

gutenberg_works(languages = c("en", "es"), all_languages = TRUE) \%>\%
  count(language, sort = TRUE)

gutenberg_works(languages = c("en", "es"), only_languages = FALSE) \%>\%
  count(language, sort = TRUE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gutenberg_download.R
\name{gutenberg_get_mirror}
\alias{gutenberg_get_mirror}
\title{Get the recommended mirror for Gutenberg files}
\usage{
gutenberg_get_mirror(verbose = TRUE)
}
\arguments{
\item{verbose}{Whether to show messages about the Project Gutenberg
mirror that was chosen}
}
\description{
Get the recommended mirror for Gutenberg files by accessing
the wget harvest path, which is
\url{http://www.gutenberg.org/robot/harvest?filetypes[]=txt}.
Also sets the global \code{gutenberg_mirror} options.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gutenberg_subjects}
\alias{gutenberg_subjects}
\title{Gutenberg metadata about the subject of each work}
\format{
A tbl_df (see tibble or dplyr) with one row for each pairing
of work and subject, with columns:
\describe{
  \item{gutenberg_id}{ID describing a work that can be joined with
  \link{gutenberg_metadata}}
  \item{subject_type}{Either "lcc" (Library of Congress Classification) or
  "lcsh" (Library of Congress Subject Headings)}
  \item{subject}{Subject}
}
}
\usage{
gutenberg_subjects
}
\description{
Gutenberg metadata about the subject of each work, particularly
Library of Congress Classifications (lcc) and Library of Congress
Subject Headings (lcsh).
}
\details{
Find more information about Library of Congress Categories
here: \url{https://www.loc.gov/catdir/cpso/lcco/}, and about
Library of Congress Subject Headings here:
\url{https://id.loc.gov/authorities/subjects.html}.

To find the date on which this metadata was last updated,
run \code{attr(gutenberg_subjects, "date_updated")}.
}
\examples{

library(dplyr)
library(stringr)

gutenberg_subjects \%>\%
  filter(subject_type == "lcsh") \%>\%
  count(subject, sort = TRUE)

sherlock_holmes_subjects <- gutenberg_subjects \%>\%
  filter(str_detect(subject, "Holmes, Sherlock"))

sherlock_holmes_subjects

sherlock_holmes_metadata <- gutenberg_works() \%>\%
  filter(author == "Doyle, Arthur Conan") \%>\%
  semi_join(sherlock_holmes_subjects, by = "gutenberg_id")

sherlock_holmes_metadata

\dontrun{
holmes_books <- gutenberg_download(sherlock_holmes_metadata$gutenberg_id)

holmes_books
}

# date last updated
attr(gutenberg_subjects, "date_updated")

}
\seealso{
\link{gutenberg_metadata}, \link{gutenberg_authors}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gutenberg_download.R
\name{gutenberg_download}
\alias{gutenberg_download}
\title{Download one or more works using a Project Gutenberg ID}
\usage{
gutenberg_download(
  gutenberg_id,
  mirror = NULL,
  strip = TRUE,
  meta_fields = NULL,
  verbose = TRUE,
  files = NULL,
  ...
)
}
\arguments{
\item{gutenberg_id}{A vector of Project Gutenberg ID, or a data frame
containing a \code{gutenberg_id} column, such as from the results of
a \code{gutenberg_works()} call}

\item{mirror}{Optionally a mirror URL to retrieve the books from. By
default uses the mirror from \code{\link{gutenberg_get_mirror}}}

\item{strip}{Whether to strip suspected headers and footers using the
\code{\link{gutenberg_strip}} function}

\item{meta_fields}{Additional fields, such as \code{title} and \code{author},
to add from \link{gutenberg_metadata} describing each book. This is useful
when returning multiple}

\item{verbose}{Whether to show messages about the Project Gutenberg
mirror that was chosen}

\item{files}{A vector of .zip file paths. If given, this reads from the
files rather than from the site. This is mostly used for testing when
the Project Gutenberg website may not be available.}

\item{...}{Extra arguments passed to \code{\link{gutenberg_strip}}, currently
unused}
}
\value{
A two column tbl_df (a type of data frame; see tibble or
dplyr packages) with one row for each line of the text or texts,
with columns
\describe{
  \item{gutenberg_id}{Integer column with the Project Gutenberg ID of
  each text}
  \item{text}{A character vector}
}
}
\description{
Download one or more works by their Project Gutenberg IDs into
a data frame with one row per line per work. This can be used to download
a single work of interest or multiple at a time. You can look up the
Gutenberg IDs of a work using the \code{gutenberg_works()} function or
the \code{gutenberg_metadata} dataset.
}
\details{
Note that if \code{strip = TRUE}, this tries to remove the
Gutenberg header and footer using the \code{\link{gutenberg_strip}}
function. This is not an exact process since headers and footers differ
between books. Before doing an in-depth analysis you may want to check
the start and end of each downloaded book.
}
\examples{

\dontrun{
library(dplyr)

# download The Count of Monte Cristo
gutenberg_download(1184)

# download two books: Wuthering Heights and Jane Eyre
books <- gutenberg_download(c(768, 1260), meta_fields = "title")
books
books \%>\% count(title)

# download all books from Jane Austen
austen <- gutenberg_works(author == "Austen, Jane") \%>\%
  gutenberg_download(meta_fields = "title")

austen
austen \%>\%
 count(title)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gutenberg_authors}
\alias{gutenberg_authors}
\title{Metadata about Project Gutenberg authors}
\format{
A tbl_df (see tibble or dplyr) with one row for each
author, with the columns
\describe{
  \item{gutenberg_author_id}{Unique identifier for the author that can
  be used to join with the \link{gutenberg_metadata} dataset}
  \item{author}{The \code{agent_name} field from the original metadata}
  \item{alias}{Alias}
  \item{birthdate}{Year of birth}
  \item{deathdate}{Year of death}
  \item{wikipedia}{Link to Wikipedia article on the author. If there
  are multiple, they are "/"-delimited}
  \item{aliases}{Character vector of aliases. If there
  are multiple, they are "/"-delimited}
}
}
\usage{
gutenberg_authors
}
\description{
Data frame with metadata about each author of a Project
Gutenberg work. Although the Project Gutenberg raw data
also includes metadata on contributors, editors, illustrators,
etc., this dataset contains only people who have been the
single author of at least one work.
}
\details{
To find the date on which this metadata was last updated,
run \code{attr(gutenberg_authors, "date_updated")}.
}
\examples{

# date last updated
attr(gutenberg_authors, "date_updated")

}
\seealso{
\link{gutenberg_metadata}, \link{gutenberg_subjects}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gutenberg_metadata}
\alias{gutenberg_metadata}
\title{Gutenberg metadata about each work}
\format{
A tbl_df (see tibble or dplyr) with one row for each work in Project Gutenberg
and the following columns:
\describe{
  \item{gutenberg_id}{Numeric ID, used to retrieve works from
  Project Gutenberg}
  \item{title}{Title}
  \item{author}{Author, if a single one given. Given as last name
  first (e.g. "Doyle, Arthur Conan")}
  \item{author_id}{Project Gutenberg author ID}
  \item{language}{Language ISO 639 code, separated by / if multiple. Two
  letter code if one exists, otherwise three letter. See
  \url{https://en.wikipedia.org/wiki/List_of_ISO_639-2_codes}}
  \item{gutenberg_bookshelf}{Which collection or collections this
  is found in, separated by / if multiple}
  \item{rights}{Generally one of three options: "Public domain in the USA."
  (the most common by far), "Copyrighted. Read the copyright notice inside this book
  for details.", or "None"}
  \item{has_text}{Whether there is a file containing digits followed by
  \code{.txt} in Project Gutenberg for this record (as opposed to, for
  example, audiobooks). If not, cannot be retrieved with
  \code{\link{gutenberg_download}}}
}
}
\usage{
gutenberg_metadata
}
\description{
Selected fields of metadata about each of the Project Gutenberg
works. These were collected using the gitenberg Python package,
particularly the \code{pg_rdf_to_json} function.
}
\details{
To find the date on which this metadata was last updated,
run \code{attr(gutenberg_metadata, "date_updated")}.
}
\examples{

library(dplyr)
library(stringr)

gutenberg_metadata

gutenberg_metadata \%>\%
  count(author, sort = TRUE)

# look for Shakespeare, excluding collections (containing "Works") and translations
shakespeare_metadata <- gutenberg_metadata \%>\%
  filter(author == "Shakespeare, William",
         language == "en",
         !str_detect(title, "Works"),
         has_text,
         !str_detect(rights, "Copyright")) \%>\%
         distinct(title)

\dontrun{
shakespeare_works <- gutenberg_download(shakespeare_metadata$gutenberg_id)
}

# note that the gutenberg_works() function filters for English
# non-copyrighted works and does de-duplication by default:

shakespeare_metadata2 <- gutenberg_works(author == "Shakespeare, William",
                                         !str_detect(title, "Works"))

# date last updated
attr(gutenberg_metadata, "date_updated")

}
\seealso{
\link{gutenberg_works}, \link{gutenberg_authors},
\link{gutenberg_subjects}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{read_zip_url}
\alias{read_zip_url}
\title{Read a file from a .zip URL}
\usage{
read_zip_url(url)
}
\arguments{
\item{url}{URL to a .zip file}
}
\description{
Download, read, and delete a .zip file
}
