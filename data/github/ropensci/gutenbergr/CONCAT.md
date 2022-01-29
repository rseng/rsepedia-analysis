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
