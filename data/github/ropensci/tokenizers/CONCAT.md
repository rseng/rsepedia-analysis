---
title: "Fast, Consistent Tokenization of Natural Language Text"
tags:
- text mining
- tokenization
- natural language processing
authors:
- name: Lincoln A. Mullen
  orcid: 0000-0001-5103-6917
  affiliation: 1
- name: Kenneth Benoit
  orcid: 0000-0002-0797-564X
  affiliation: 2
- name: Os Keyes
  orcid: 0000-0001-5196-609X
  affiliation: 3
- name: Dmitry Selivanov
  affiliation: 4
- name: Jeffrey Arnold
  orcid: 0000-0001-9953-3904
  affiliation: 5
affiliations: 
- name: "Department of History and Art History, George Mason University"
  index: 1
- name: "Department of Methodology, London School of Economics and Political Science"
  index: 2
- name: "Department of Human Centered Design and Engineering, University of Washington"
  index: 3
- name: "Open Data Science"
  index: 4
- name: "Department of Political Science, University of Washington"
  index: 5
date: 12 March 2018
bibliography: paper.bib
...

Computational text analysis usually proceeds according to a series of 
well-defined steps. After importing texts, the usual next step is to turn the 
human-readable text into machine-readable tokens. Tokens are defined as 
segments of a text identified as meaningful units for the purpose of analyzing 
the text. They may consist of individual words or of larger or smaller 
segments, such as word sequences, word subsequences, paragraphs, sentences, or 
lines [@Manningetal2008, 22]. Tokenization is the process of splitting the text 
into these smaller pieces, and it often involves preprocessing the text to 
remove punctuation and transform all tokens into lowercase [@welbers_text_2017, 
250-251]. Decisions made during tokenization have a significant effect on 
subsequent analysis [@denny_text_forthcoming; @guthrie_closer_2006]. Especially 
for large corpora, tokenization can be computationally expensive, and 
tokenization is highly language dependent.  Efficiency and correctness are 
therefore paramount concerns for tokenization.

The [tokenizers](http://lincolnmullen.com/software/tokenizers/) package for 
R provides fast, consistent tokenization for natural language text 
[@tokenizers; @rbase]. (The package is available on 
[GitHub](https://github.com/ropensci/tokenizers) and archived on 
[Zenodo](https://doi.org/10.5281/zenodo.1205017).) Each of the tokenizers 
expects a consistent input and returns a consistent output, so that the 
tokenizers can be used interchangeably with one another or relied on in other 
packages. To ensure the correctness of output, the package depends on the 
stringi package, which implements Unicode support for R [@gagolewski_2018]. 
To ensure the speed of tokenization, key components such as the _n_-gram and 
skip _n_-gram tokenizers are written using the Rcpp package 
[@eddelbuettel_2013; @eddelbuettel_2017]. The tokenizers package is part of 
the [rOpenSci project](https://ropensci.org/).

The most important tokenizers in the current version of the package can be 
grouped as follows:

- tokenizers for characters and shingled characters
- tokenizers for words and word stems, as well as for Penn Treebank tokens 
- tokenizers _n_-grams and skip _n_-grams
- tokenizers for tweets, which preserve formatting of usernames and hashtags

In addition the package provides functions for splitting longer documents into 
sentences and paragraphs, or for splitting a long text into smaller chunks each 
with the same number of words. This allows users to treat parts of very long 
texts as documents in their own right. The package also provides functions for 
counting words, characters, and sentences.

The tokenizers in this package can be used on their own, or they can be wrapped 
by higher-level R packages. For instance, the tokenizers package is a 
dependency for the tidytext [@silge_2016], text2vec [@selivanov_2018], 
and textreuse [@mullen_2016] packages. More broadly, the output of the 
tokenization functions follows the guidelines set by the text-interchange 
format  defined at an rOpenSci Text Workshop in 2017 [@tif_2017]. Other 
packages which buy into the text-interchange format can thus use the 
tokenizers package interchangeably.

The tokenizers package has research applications in any discipline which 
uses computational text analysis. The package was originally created for  
historical research into the use of the Bible in American newspapers 
[@mullen_americas] and into the borrowing of legal codes of civil procedure in 
the nineteenth-century United States [@funkmullen_spine_2018, 
@funkmullen_servile_2016]. The tokenizers package underlies the tidytext 
package [@silge_text_2017], and via that package tokenizers has been used 
in disciplines such as political science [@sanger_2015_], social science 
[@warin_mapping], communication studies [@xu_using_2018], English 
[@ballier_rbased_2017], and the digital humanities more generally.

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

# tokenizers

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/tokenizers)](https://cran.r-project.org/package=tokenizers)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00655/status.svg)](https://doi.org/10.21105/joss.00655)
[![rOpenSci peer
review](https://badges.ropensci.org/33_status.svg)](https://github.com/ropensci/onboarding/issues/33)
[![CRAN\_Downloads](http://cranlogs.r-pkg.org/badges/grand-total/tokenizers)](https://cran.r-project.org/package=tokenizers)
[![Travis-CI Build
Status](https://travis-ci.org/ropensci/tokenizers.svg?branch=master)](https://travis-ci.org/ropensci/tokenizers)
[![Appveyor Build
status](https://ci.appveyor.com/api/projects/status/qx3vh3ukjgo99iu4/branch/master?svg=true)](https://ci.appveyor.com/project/lmullen/tokenizers-dkf3v/branch/master)
[![Coverage
Status](https://img.shields.io/codecov/c/github/ropensci/tokenizers/master.svg)](https://codecov.io/github/ropensci/tokenizers?branch=master)

## Overview

This R package offers functions with a consistent interface to convert
natural language text into tokens. It includes tokenizers for shingled
n-grams, skip n-grams, words, word stems, sentences, paragraphs,
characters, shingled characters, lines, tweets, Penn Treebank, and
regular expressions, as well as functions for counting characters,
words, and sentences, and a function for splitting longer texts into
separate documents, each with the same number of words. The package is
built on the [stringi](http://www.gagolewski.com/software/stringi/) and
[Rcpp](http://www.rcpp.org/) packages for fast yet correct tokenization
in UTF-8.

See the “[Introduction to the tokenizers
Package](http://lincolnmullen.com/software/tokenizers/articles/introduction-to-tokenizers.html)”
vignette for an overview of all the functions in this package.

This package complies with the standards for input and output
recommended by the Text Interchange Formats. The TIF initiative was
created at an rOpenSci meeting in 2017, and its recommendations are
available as part of the [tif package](https://github.com/ropensci/tif).
See the “[The Text Interchange Formats and the tokenizers
Package](http://lincolnmullen.com/software/tokenizers/articles/tif-and-tokenizers.html)”
vignette for an explanation of how this package fits into that
ecosystem.

## Suggested citation

If you use this package for your research, we would appreciate a
citation.

``` r
citation("tokenizers")
#> 
#> To cite the tokenizers package in publications, please cite the
#> paper in the Journal of Open Source Software:
#> 
#>   Lincoln A. Mullen et al., "Fast, Consistent Tokenization of
#>   Natural Language Text," Journal of Open Source Software 3, no.
#>   23 (2018): 655, https://doi.org/10.21105/joss.00655.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Fast, Consistent Tokenization of Natural Language Text},
#>     author = {Lincoln A. Mullen and Kenneth Benoit and Os Keyes and Dmitry Selivanov and Jeffrey Arnold},
#>     journal = {Journal of Open Source Software},
#>     year = {2018},
#>     volume = {3},
#>     issue = {23},
#>     pages = {655},
#>     url = {https://doi.org/10.21105/joss.00655},
#>     doi = {10.21105/joss.00655},
#>   }
```

## Installation

You can install this package from CRAN:

``` r
install.packages("tokenizers")
```

To get the development version from GitHub, use
[devtools](https://github.com/hadley/devtools).

``` r
# install.packages("devtools")
devtools::install_github("ropensci/tokenizers")
```

## Examples

The tokenizers in this package have a consistent interface. They all
take either a character vector of any length, or a list where each
element is a character vector of length one, or a data.frame that
adheres to the [tif corpus format](https://github.com/ropensci/tif). The
idea is that each element (or row) comprises a text. Then each function
returns a list with the same length as the input vector, where each
element in the list contains the tokens generated by the function. If
the input character vector or list is named, then the names are
preserved, so that the names can serve as identifiers. For a
tif-formatted data.frame, the `doc_id` field is used as the element
names in the returned token list.

``` r
library(magrittr)
library(tokenizers)

james <- paste0(
  "The question thus becomes a verbal one\n",
  "again; and our knowledge of all these early stages of thought and feeling\n",
  "is in any case so conjectural and imperfect that farther discussion would\n",
  "not be worth while.\n",
  "\n",
  "Religion, therefore, as I now ask you arbitrarily to take it, shall mean\n",
  "for us _the feelings, acts, and experiences of individual men in their\n",
  "solitude, so far as they apprehend themselves to stand in relation to\n",
  "whatever they may consider the divine_. Since the relation may be either\n",
  "moral, physical, or ritual, it is evident that out of religion in the\n",
  "sense in which we take it, theologies, philosophies, and ecclesiastical\n",
  "organizations may secondarily grow.\n"
)
names(james) <- "varieties"

tokenize_characters(james)[[1]] %>% head(50)
#>  [1] "t" "h" "e" "q" "u" "e" "s" "t" "i" "o" "n" "t" "h" "u" "s" "b" "e"
#> [18] "c" "o" "m" "e" "s" "a" "v" "e" "r" "b" "a" "l" "o" "n" "e" "a" "g"
#> [35] "a" "i" "n" "a" "n" "d" "o" "u" "r" "k" "n" "o" "w" "l" "e" "d"
tokenize_character_shingles(james)[[1]] %>% head(20)
#>  [1] "the" "heq" "equ" "que" "ues" "est" "sti" "tio" "ion" "ont" "nth"
#> [12] "thu" "hus" "usb" "sbe" "bec" "eco" "com" "ome" "mes"
tokenize_words(james)[[1]] %>% head(10)
#>  [1] "the"      "question" "thus"     "becomes"  "a"        "verbal"  
#>  [7] "one"      "again"    "and"      "our"
tokenize_word_stems(james)[[1]] %>% head(10)
#>  [1] "the"      "question" "thus"     "becom"    "a"        "verbal"  
#>  [7] "one"      "again"    "and"      "our"
tokenize_sentences(james) 
#> $varieties
#> [1] "The question thus becomes a verbal one again; and our knowledge of all these early stages of thought and feeling is in any case so conjectural and imperfect that farther discussion would not be worth while."                                               
#> [2] "Religion, therefore, as I now ask you arbitrarily to take it, shall mean for us _the feelings, acts, and experiences of individual men in their solitude, so far as they apprehend themselves to stand in relation to whatever they may consider the divine_."
#> [3] "Since the relation may be either moral, physical, or ritual, it is evident that out of religion in the sense in which we take it, theologies, philosophies, and ecclesiastical organizations may secondarily grow."
tokenize_paragraphs(james)
#> $varieties
#> [1] "The question thus becomes a verbal one again; and our knowledge of all these early stages of thought and feeling is in any case so conjectural and imperfect that farther discussion would not be worth while."                                                                                                                                                                                                                                                                   
#> [2] "Religion, therefore, as I now ask you arbitrarily to take it, shall mean for us _the feelings, acts, and experiences of individual men in their solitude, so far as they apprehend themselves to stand in relation to whatever they may consider the divine_. Since the relation may be either moral, physical, or ritual, it is evident that out of religion in the sense in which we take it, theologies, philosophies, and ecclesiastical organizations may secondarily grow. "
tokenize_ngrams(james, n = 5, n_min = 2)[[1]] %>% head(10)
#>  [1] "the question"                   "the question thus"             
#>  [3] "the question thus becomes"      "the question thus becomes a"   
#>  [5] "question thus"                  "question thus becomes"         
#>  [7] "question thus becomes a"        "question thus becomes a verbal"
#>  [9] "thus becomes"                   "thus becomes a"
tokenize_skip_ngrams(james, n = 5, k = 2)[[1]] %>% head(10)
#>  [1] "the"                  "the question"         "the thus"            
#>  [4] "the becomes"          "the question thus"    "the question becomes"
#>  [7] "the question a"       "the thus becomes"     "the thus a"          
#> [10] "the thus verbal"
tokenize_ptb(james)[[1]] %>% head(10)
#>  [1] "The"      "question" "thus"     "becomes"  "a"        "verbal"  
#>  [7] "one"      "again"    ";"        "and"
tokenize_lines(james)[[1]] %>% head(5)
#> [1] "The question thus becomes a verbal one"                                   
#> [2] "again; and our knowledge of all these early stages of thought and feeling"
#> [3] "is in any case so conjectural and imperfect that farther discussion would"
#> [4] "not be worth while."                                                      
#> [5] "Religion, therefore, as I now ask you arbitrarily to take it, shall mean"
tokenize_tweets("Hey @handle, #rstats is awesome!")[[1]]
#> [1] "hey"     "@handle" "#rstats" "is"      "awesome"
```

The package also contains functions to count words, characters, and
sentences, and these functions follow the same consistent interface.

``` r
count_words(james)
#> varieties 
#>       112
count_characters(james)
#> varieties 
#>       673
count_sentences(james)
#> varieties 
#>        13
```

The `chunk_text()` function splits a document into smaller chunks, each
with the same number of words.

## Contributing

Contributions to the package are more than welcome. One way that you can
help is by using this package in your R package for natural language
processing. If you want to contribute a tokenization function to this
package, it should follow the same conventions as the rest of the
functions whenever it makes sense to do so.

Please note that this project is released with a [Contributor Code of
Conduct](CONDUCT.md). By participating in this project you agree to
abide by its terms.

-----

[![rOpenSCi
logo](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# tokenizers 0.2.1

- Add citation information to JOSS paper.

# tokenizers 0.2.0

## Features

- Add the `tokenize_ptb()` function for Penn Treebank tokenizations (@jrnold) (#12).
- Add a function `chunk_text()` to split long documents into pieces (#30).
- New functions to count words, characters, and sentences without tokenization (#36).
- New function `tokenize_tweets()` preserves usernames, hashtags, and URLS (@kbenoit) (#44).
- The `stopwords()` function has been removed in favor of using the **stopwords** package (#46).
- The package now complies with the basic recommendations of the **Text Interchange Format**. All tokenization functions are now methods. This enables them to take corpus inputs as either TIF-compliant named character vectors, named lists, or data frames. All outputs are still named lists of tokens, but these can be easily coerced to data frames of tokens using the `tif` package. (#49)
- Add a new vignette "The Text Interchange Formats and the tokenizers Package" (#49).

## Bug fixes and performance improvements

- `tokenize_skip_ngrams` has been improved to generate unigrams and bigrams, according to the skip definition (#24).
- C++98 has replaced the C++11 code used for n-gram generation, widening the range of compilers `tokenizers` supports (@ironholds) (#26).
- `tokenize_skip_ngrams` now supports stopwords (#31).
- If tokenisers fail to generate tokens for a particular entry, they return `NA` consistently (#33).
- Keyboard interrupt checks have been added to Rcpp-backed functions to enable users to terminate them before completion (#37).
- `tokenize_words()` gains arguments to preserve or strip punctuation and numbers (#48).
- `tokenize_skip_ngrams()` and `tokenize_ngrams()` to return properly marked UTF8 strings on Windows (@patperry) (#58).
- `tokenize_tweets()` now removes stopwords prior to stripping punctuation, making its behavior more consistent with `tokenize_words()` (#76).

# tokenizers 0.1.4

- Add the `tokenize_character_shingles()` tokenizer.
- Improvements to documentation.

# tokenizers 0.1.3

- Add vignette.
- Improvements to n-gram tokenizers.

# tokenizers 0.1.2

- Add stopwords for several languages.
- New stopword options to `tokenize_words()` and `tokenize_word_stems()`.

# tokenizers 0.1.1

- Fix failing test in non-UTF-8 locales.

# tokenizers 0.1.0

- Initial release with tokenizers for characters, words, word stems, sentences
  paragraphs, n-grams, skip n-grams, lines, and regular expressions.
This is a minor update to fix bugs and add peer review and citation information.

## Test environments

* Local OS X install: R-Release
* Ubuntu (on Travis-CI): R-release, R-devel, R-oldrel
* Local Ubuntu 16.04 install: R-release
* Win-builder: R-devel

## R CMD check results

* One NOTE pertains to non-ASCII strings in test files only, which are
  necessary to ensure the package's functionality on Windows.
---
output: github_document
pagetitle: "tokenizers: Fast, Consistent Tokenization of Natural Language Text"
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# tokenizers

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tokenizers)](https://cran.r-project.org/package=tokenizers)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00655/status.svg)](https://doi.org/10.21105/joss.00655)
[![rOpenSci peer review](https://badges.ropensci.org/33_status.svg)](https://github.com/ropensci/onboarding/issues/33)
[![CRAN_Downloads](http://cranlogs.r-pkg.org/badges/grand-total/tokenizers)](https://cran.r-project.org/package=tokenizers)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/tokenizers.svg?branch=master)](https://travis-ci.org/ropensci/tokenizers)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/qx3vh3ukjgo99iu4/branch/master?svg=true)](https://ci.appveyor.com/project/lmullen/tokenizers-dkf3v/branch/master)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/tokenizers/master.svg)](https://codecov.io/github/ropensci/tokenizers?branch=master)

## Overview

This R package offers functions with a consistent interface to convert natural language text into tokens. It includes tokenizers for shingled n-grams, skip n-grams, words, word stems, sentences, paragraphs, characters, shingled characters, lines, tweets, Penn Treebank, and regular expressions, as well as functions for counting characters, words, and sentences, and a function for splitting longer texts into separate documents, each with the same number of words. The package is built on the [stringi](http://www.gagolewski.com/software/stringi/) and [Rcpp](http://www.rcpp.org/) packages for fast yet correct tokenization in UTF-8. 

See the "[Introduction to the tokenizers Package](http://lincolnmullen.com/software/tokenizers/articles/introduction-to-tokenizers.html)" vignette for an overview of all the functions in this package.

This package complies with the standards for input and output recommended by the Text Interchange Formats. The TIF initiative was created at an rOpenSci meeting in 2017, and its recommendations are available as part of the [tif package](https://github.com/ropensci/tif). See the "[The Text Interchange Formats and the tokenizers Package](http://lincolnmullen.com/software/tokenizers/articles/tif-and-tokenizers.html)" vignette for an explanation of how this package fits into that ecosystem.

## Suggested citation

If you use this package for your research, we would appreciate a citation.

```{r}
citation("tokenizers")
```

## Installation

You can install this package from CRAN:

```{r eval=FALSE}
install.packages("tokenizers")
```

To get the development version from GitHub, use  [devtools](https://github.com/hadley/devtools).

```{r eval=FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/tokenizers")
```

## Examples

The tokenizers in this package have a consistent interface. They all take either a character vector of any length, or a list where each element is a character vector of length one, or a data.frame that adheres to the [tif corpus format](https://github.com/ropensci/tif). The idea is that each element (or row) comprises a text. Then each function returns a list with the same length as the input vector, where each element in the list contains the tokens generated by the function.  If the input character vector or list is named, then the names are preserved, so that the names can serve as identifiers.  For a tif-formatted data.frame, the `doc_id` field is used as the element names in the returned token list.

```{r}
library(magrittr)
library(tokenizers)

james <- paste0(
  "The question thus becomes a verbal one\n",
  "again; and our knowledge of all these early stages of thought and feeling\n",
  "is in any case so conjectural and imperfect that farther discussion would\n",
  "not be worth while.\n",
  "\n",
  "Religion, therefore, as I now ask you arbitrarily to take it, shall mean\n",
  "for us _the feelings, acts, and experiences of individual men in their\n",
  "solitude, so far as they apprehend themselves to stand in relation to\n",
  "whatever they may consider the divine_. Since the relation may be either\n",
  "moral, physical, or ritual, it is evident that out of religion in the\n",
  "sense in which we take it, theologies, philosophies, and ecclesiastical\n",
  "organizations may secondarily grow.\n"
)
names(james) <- "varieties"

tokenize_characters(james)[[1]] %>% head(50)
tokenize_character_shingles(james)[[1]] %>% head(20)
tokenize_words(james)[[1]] %>% head(10)
tokenize_word_stems(james)[[1]] %>% head(10)
tokenize_sentences(james) 
tokenize_paragraphs(james)
tokenize_ngrams(james, n = 5, n_min = 2)[[1]] %>% head(10)
tokenize_skip_ngrams(james, n = 5, k = 2)[[1]] %>% head(10)
tokenize_ptb(james)[[1]] %>% head(10)
tokenize_lines(james)[[1]] %>% head(5)
tokenize_tweets("Hey @handle, #rstats is awesome!")[[1]]
```

The package also contains functions to count words, characters, and sentences, and these functions follow the same consistent interface.

```{r}
count_words(james)
count_characters(james)
count_sentences(james)
```

The `chunk_text()` function splits a document into smaller chunks, each with the same number of words.

## Contributing

Contributions to the package are more than welcome. One way that you can help is by using this package in your R package for natural language processing. If you want to contribute a tokenization function to this package, it should follow the same conventions as the rest of the functions whenever it makes sense to do so. 

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

------------------------------------------------------------------------

[![rOpenSCi logo](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "Introduction to the tokenizers Package"
author: "Lincoln Mullen"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to the tokenizers Package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Package overview

In natural language processing, tokenization is the process of breaking human-readable text into machine readable components. The most obvious way to tokenize a text is to split the text into words. But there are many other ways to tokenize a text, the most useful of which are provided by this package.

The tokenizers in this package have a consistent interface. They all take either a character vector of any length, or a list where each element is a character vector of length one. The idea is that each element comprises a text. Then each function returns a list with the same length as the input vector, where each element in the list contains the tokens generated by the function. If the input character vector or list is named, then the names are preserved, so that the names can serve as identifiers.

Using the following sample text, the rest of this vignette demonstrates the different kinds of tokenizers in this package.

```{r}
library(tokenizers)
options(max.print = 25)

james <- paste0(
  "The question thus becomes a verbal one\n",
  "again; and our knowledge of all these early stages of thought and feeling\n",
  "is in any case so conjectural and imperfect that farther discussion would\n",
  "not be worth while.\n",
  "\n",
  "Religion, therefore, as I now ask you arbitrarily to take it, shall mean\n",
  "for us _the feelings, acts, and experiences of individual men in their\n",
  "solitude, so far as they apprehend themselves to stand in relation to\n",
  "whatever they may consider the divine_. Since the relation may be either\n",
  "moral, physical, or ritual, it is evident that out of religion in the\n",
  "sense in which we take it, theologies, philosophies, and ecclesiastical\n",
  "organizations may secondarily grow.\n"
)
```

## Character and character-shingle tokenizers

The character tokenizer splits texts into individual characters. 

```{r}
tokenize_characters(james)[[1]] 
```

You can also tokenize into character-based shingles.

```{r}
tokenize_character_shingles(james, n = 3, n_min = 3, 
                            strip_non_alphanum = FALSE)[[1]][1:20]
```

## Word and word-stem tokenizers

The word tokenizer splits texts into words. 

```{r}
tokenize_words(james)
```

Word stemming is provided by the [SnowballC](https://cran.r-project.org/package=SnowballC) package.

```{r}
tokenize_word_stems(james)
```

You can also provide a vector of stopwords which will be omitted. The [stopwords package](https://github.com/quanteda/stopwords), which contains stopwords for many languages from several sources, is recommended. This argument also works with the n-gram and skip n-gram tokenizers.

```{r}
library(stopwords)
tokenize_words(james, stopwords = stopwords::stopwords("en"))
```

An alternative word stemmer often used in NLP that preserves punctuation and separates common English contractions is the Penn Treebank tokenizer.

```{r}
tokenize_ptb(james)
```

## N-gram and skip n-gram tokenizers

An n-gram is a contiguous sequence of words containing at least `n_min` words and at most `n` words. This function will generate all such combinations of n-grams, omitting stopwords if desired.

```{r}
tokenize_ngrams(james, n = 5, n_min = 2,
                stopwords = stopwords::stopwords("en"))
```

A skip n-gram is like an n-gram in that it takes the `n` and `n_min` parameters. But rather than returning contiguous sequences of words, it will also return sequences of n-grams skipping words with gaps between `0` and the value of `k`. This function generates all such sequences, again omitting stopwords if desired. Note that the number of tokens returned can be very large.

```{r}
tokenize_skip_ngrams(james, n = 5, n_min = 2, k = 2,
                     stopwords = stopwords::stopwords("en"))
```

## Tweet tokenizer

Tokenizing tweets requires special attention, since usernames (`@whoever`) and hashtags (`#hashtag`) use special characters that might otherwise be stripped away.

```{r}
tokenize_tweets("Welcome, @user, to the tokenizers package. #rstats #forever")
```

## Sentence and paragraph tokenizers

Sometimes it is desirable to split texts into sentences or paragraphs prior to tokenizing into other forms.

```{r, collapse=FALSE}
tokenize_sentences(james) 
tokenize_paragraphs(james)
```

## Text chunking

When one has a very long document, sometimes it is desirable to split the document into smaller chunks, each with the same length. This function chunks a document and gives it each of the chunks an ID to show their order. These chunks can then be further tokenized.

```{r}
chunks <- chunk_text(mobydick, chunk_size = 100, doc_id = "mobydick")
length(chunks)
chunks[5:6]
tokenize_words(chunks[5:6])
```

## Counting words, characters, sentences

The package also offers functions for counting words, characters, and sentences in a format which works nicely with the rest of the functions.

```{r}
count_words(mobydick)
count_characters(mobydick)
count_sentences(mobydick)
```

---
title: "The Text Interchange Formats and the tokenizers Package"
author: "Lincoln Mullen"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{The Text Interchange Formats and the tokenizers Package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The [Text Interchange Formats](https://github.com/ropensci/tif) are a set of standards defined at an [rOpenSci](https://ropensci.org/) sponsored [meeting in London](http://textworkshop17.ropensci.org/) in 2017. The formats allow R text analysis packages to target defined inputs and outputs for corpora, tokens, and document-term matrices. By adhering to these recommendations, R packages can buy into an interoperable ecosystem.

The TIF recommendations are still a draft, but the tokenizers package implements its recommendation to accept both of the corpora formats and to output one of its recommended tokens formats. 

Consider these two recommended forms of a corpus. One (`corpus_c`) is a named character vector; the other (`corpus_d`) is a data frame. They both include a document ID and the full text for each item. The data frame format obviously allows for the use of other metadata fields besides the document ID, whereas the other format does not. Using the coercion functions in the tif package, one could switch back and forth between these formats. Tokenizers also supports a corpus formatted as a named list where each element is a character vector of length one (`corpus_l`), though this is not a part of the draft TIF standards.

```{r}
# Named list
(corpus_l <- list(man_comes_around = "There's a man goin' 'round takin' names",
                  wont_back_down = "Well I won't back down, no I won't back down",
                  bird_on_a_wire = "Like a bird on a wire"))

# Named character vector
(corpus_c <- unlist(corpus_l))

# Data frame
(corpus_d <- data.frame(doc_id = names(corpus_c), text = unname(corpus_c),
                        stringsAsFactors = FALSE))
```

All of the tokenizers in this package can accept any of those formats and will return an identical output for each.

```{r}
library(tokenizers)

tokens_l <- tokenize_ngrams(corpus_l, n = 2)
tokens_c <- tokenize_ngrams(corpus_c, n = 2)
tokens_d <- tokenize_ngrams(corpus_c, n = 2)

# Are all these identical?
all(identical(tokens_l, tokens_c),
    identical(tokens_c, tokens_d),
    identical(tokens_l, tokens_d))
```

The output of all of the tokenizers is a named list, where each element of the list corresponds to a document in the corpus. The names of the list are the document IDs, and the elements are character vectors containing the tokens.

```{r}
tokens_l
```

This format can be coerced to a data frame of document IDs and tokens, one row per token, using the coercion functions in the tif package. That tokens data frame would look like this.

```{r, echo=FALSE}
sample_tokens_df <- structure(list(doc_id = c("man_comes_around", "man_comes_around", 
"man_comes_around", "man_comes_around", "man_comes_around", "man_comes_around", 
"wont_back_down", "wont_back_down", "wont_back_down", "wont_back_down", 
"wont_back_down", "wont_back_down", "wont_back_down", "wont_back_down", 
"wont_back_down", "bird_on_a_wire", "bird_on_a_wire", "bird_on_a_wire", 
"bird_on_a_wire", "bird_on_a_wire"), token = c("there's a", "a man", 
"man goin", "goin round", "round takin", "takin names", "well i", 
"i won't", "won't back", "back down", "down no", "no i", "i won't", 
"won't back", "back down", "like a", "a bird", "bird on", "on a", 
"a wire")), .Names = c("doc_id", "token"), row.names = c(NA, 
-20L), class = "data.frame")
head(sample_tokens_df, 10)
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stem-tokenizers.R
\name{tokenize_word_stems}
\alias{tokenize_word_stems}
\title{Word stem tokenizer}
\usage{
tokenize_word_stems(
  x,
  language = "english",
  stopwords = NULL,
  simplify = FALSE
)
}
\arguments{
\item{x}{A character vector or a list of character vectors to be tokenized.
If \code{x} is a character vector, it can be of any length, and each
element will be tokenized separately. If \code{x} is a list of character
vectors, where each element of the list should have a length of 1.}

\item{language}{The language to use for word stemming. This must be one of
the languages available in the SnowballC package. A list is provided by
\code{\link[SnowballC]{getStemLanguages}}.}

\item{stopwords}{A character vector of stop words to be excluded}

\item{simplify}{\code{FALSE} by default so that a consistent value is
returned regardless of length of input. If \code{TRUE}, then an input with
a single element will return a character vector of tokens instead of a
list.}
}
\value{
A list of character vectors containing the tokens, with one element
  in the list for each element that was passed as input. If \code{simplify =
  TRUE} and only a single element was passed as input, then the output is a
  character vector of tokens.
}
\description{
This function turns its input into a character vector of word stems. This is
just a wrapper around the \code{\link[SnowballC]{wordStem}} function from the
SnowballC package which does the heavy lifting, but this function provides a
consistent interface with the rest of the tokenizers in this package. The
input can be a character vector of any length, or a list of character vectors
where each character vector in the list has a length of 1.
}
\details{
This function will strip all white space and punctuation and make
  all word stems lowercase.
}
\examples{
song <-  paste0("How many roads must a man walk down\n",
                "Before you call him a man?\n",
                "How many seas must a white dove sail\n",
                "Before she sleeps in the sand?\n",
                "\n",
                "How many times must the cannonballs fly\n",
                "Before they're forever banned?\n",
                "The answer, my friend, is blowin' in the wind.\n",
                "The answer is blowin' in the wind.\n")

tokenize_word_stems(song)
}
\seealso{
\code{\link[SnowballC]{wordStem}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ptb-tokenizer.R
\name{tokenize_ptb}
\alias{tokenize_ptb}
\title{Penn Treebank Tokenizer}
\usage{
tokenize_ptb(x, lowercase = FALSE, simplify = FALSE)
}
\arguments{
\item{x}{A character vector or a list of character vectors to be tokenized
into n-grams. If \code{x} is a character vector, it can be of any length,
and each element will be tokenized separately. If \code{x} is a list of
character vectors, each element of the list should have a length of 1.}

\item{lowercase}{Should the tokens be made lower case?}

\item{simplify}{\code{FALSE} by default so that a consistent value is
returned regardless of length of input. If \code{TRUE}, then an input with
a single element will return a character vector of tokens instead of a
list.}
}
\value{
A list of character vectors containing the tokens, with one element
  in the list for each element that was passed as input. If \code{simplify =
  TRUE} and only a single element was passed as input, then the output is a
  character vector of tokens.
}
\description{
This function implements the Penn Treebank word tokenizer.
}
\details{
This tokenizer uses regular expressions to tokenize text similar to
  the tokenization used in the Penn Treebank. It assumes that text has
  already been split into sentences. The tokenizer does the following:

  \itemize{ \item{splits common English contractions, e.g. \verb{don't} is
  tokenized into \verb{do n't} and \verb{they'll} is tokenized into ->
  \verb{they 'll},} \item{handles punctuation characters as separate tokens,}
  \item{splits commas and single quotes off from words, when they are
  followed by whitespace,} \item{splits off periods that occur at the end of
  the sentence.} }

This function is a port of the Python NLTK version of the Penn
  Treebank Tokenizer.
}
\examples{
song <- list(paste0("How many roads must a man walk down\n",
                    "Before you call him a man?"),
             paste0("How many seas must a white dove sail\n",
                    "Before she sleeps in the sand?\n"),
             paste0("How many times must the cannonballs fly\n",
                    "Before they're forever banned?\n"),
             "The answer, my friend, is blowin' in the wind.",
             "The answer is blowin' in the wind.")
tokenize_ptb(song)
tokenize_ptb(c("Good muffins cost $3.88\nin New York. Please buy me\ntwo of them.",
  "They'll save and invest more.",
  "Hi, I can't say hello."))
}
\references{
\href{http://www.nltk.org/_modules/nltk/tokenize/treebank.html#TreebankWordTokenizer}{NLTK
TreebankWordTokenizer}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wordcount.R
\name{count_words}
\alias{count_words}
\alias{count_characters}
\alias{count_sentences}
\title{Count words, sentences, characters}
\usage{
count_words(x)

count_characters(x)

count_sentences(x)
}
\arguments{
\item{x}{A character vector or a list of character vectors. If \code{x} is a
character vector, it can be of any length, and each element will be
tokenized separately. If \code{x} is a list of character vectors, each
element of the list should have a length of 1.}
}
\value{
An integer vector containing the counted elements. If the input
  vector or list has names, they will be preserved.
}
\description{
Count words, sentences, and characters in input texts. These functions use
the \code{stringi} package, so they handle the counting of Unicode strings
(e.g., characters with diacritical marks) in a way that makes sense to people
counting characters.
}
\examples{
count_words(mobydick)
count_sentences(mobydick)
count_characters(mobydick)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic-tokenizers.R, R/tokenize_tweets.R
\name{basic-tokenizers}
\alias{basic-tokenizers}
\alias{tokenize_characters}
\alias{tokenize_words}
\alias{tokenize_sentences}
\alias{tokenize_lines}
\alias{tokenize_paragraphs}
\alias{tokenize_regex}
\alias{tokenize_tweets}
\title{Basic tokenizers}
\usage{
tokenize_characters(
  x,
  lowercase = TRUE,
  strip_non_alphanum = TRUE,
  simplify = FALSE
)

tokenize_words(
  x,
  lowercase = TRUE,
  stopwords = NULL,
  strip_punct = TRUE,
  strip_numeric = FALSE,
  simplify = FALSE
)

tokenize_sentences(x, lowercase = FALSE, strip_punct = FALSE, simplify = FALSE)

tokenize_lines(x, simplify = FALSE)

tokenize_paragraphs(x, paragraph_break = "\\n\\n", simplify = FALSE)

tokenize_regex(x, pattern = "\\\\s+", simplify = FALSE)

tokenize_tweets(
  x,
  lowercase = TRUE,
  stopwords = NULL,
  strip_punct = TRUE,
  strip_url = FALSE,
  simplify = FALSE
)
}
\arguments{
\item{x}{A character vector or a list of character vectors to be tokenized.
If \code{x} is a character vector, it can be of any length, and each element
will be tokenized separately. If \code{x} is a list of character vectors,
where each element of the list should have a length of 1.}

\item{lowercase}{Should the tokens be made lower case? The default value
varies by tokenizer; it is only \code{TRUE} by default for the tokenizers
that you are likely to use last.}

\item{strip_non_alphanum}{Should punctuation and white space be stripped?}

\item{simplify}{\code{FALSE} by default so that a consistent value is
returned regardless of length of input. If \code{TRUE}, then an input with
a single element will return a character vector of tokens instead of a
list.}

\item{stopwords}{A character vector of stop words to be excluded.}

\item{strip_punct}{Should punctuation be stripped?}

\item{strip_numeric}{Should numbers be stripped?}

\item{paragraph_break}{A string identifying the boundary between two
paragraphs.}

\item{pattern}{A regular expression that defines the split.}

\item{strip_url}{Should URLs (starting with \code{http(s)}) be preserved intact, or
removed entirely?}
}
\value{
A list of character vectors containing the tokens, with one element
  in the list for each element that was passed as input. If \code{simplify =
  TRUE} and only a single element was passed as input, then the output is a
  character vector of tokens.
}
\description{
These functions perform basic tokenization into words, sentences, paragraphs,
lines, and characters. The functions can be piped into one another to create
at most two levels of tokenization. For instance, one might split a text into
paragraphs and then word tokens, or into sentences and then word tokens.
}
\examples{
song <-  paste0("How many roads must a man walk down\n",
                "Before you call him a man?\n",
                "How many seas must a white dove sail\n",
                "Before she sleeps in the sand?\n",
                "\n",
                "How many times must the cannonballs fly\n",
                "Before they're forever banned?\n",
                "The answer, my friend, is blowin' in the wind.\n",
                "The answer is blowin' in the wind.\n")

tokenize_words(song)
tokenize_words(song, strip_punct = FALSE)
tokenize_sentences(song)
tokenize_paragraphs(song)
tokenize_lines(song)
tokenize_characters(song)
tokenize_tweets("@rOpenSci and #rstats see: https://cran.r-project.org",
                strip_punct = TRUE)
tokenize_tweets("@rOpenSci and #rstats see: https://cran.r-project.org",
                strip_punct = FALSE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngram-tokenizers.R
\name{ngram-tokenizers}
\alias{ngram-tokenizers}
\alias{tokenize_ngrams}
\alias{tokenize_skip_ngrams}
\title{N-gram tokenizers}
\usage{
tokenize_ngrams(
  x,
  lowercase = TRUE,
  n = 3L,
  n_min = n,
  stopwords = character(),
  ngram_delim = " ",
  simplify = FALSE
)

tokenize_skip_ngrams(
  x,
  lowercase = TRUE,
  n_min = 1,
  n = 3,
  k = 1,
  stopwords = character(),
  simplify = FALSE
)
}
\arguments{
\item{x}{A character vector or a list of character vectors to be tokenized
into n-grams. If \code{x} is a character vector, it can be of any length,
and each element will be tokenized separately. If \code{x} is a list of
character vectors, each element of the list should have a length of 1.}

\item{lowercase}{Should the tokens be made lower case?}

\item{n}{The number of words in the n-gram. This must be an integer greater
than or equal to 1.}

\item{n_min}{The minimum number of words in the n-gram. This must be an
integer greater than or equal to 1, and less than or equal to \code{n}.}

\item{stopwords}{A character vector of stop words to be excluded from the
n-grams.}

\item{ngram_delim}{The separator between words in an n-gram.}

\item{simplify}{\code{FALSE} by default so that a consistent value is
returned regardless of length of input. If \code{TRUE}, then an input with
a single element will return a character vector of tokens instead of a
list.}

\item{k}{For the skip n-gram tokenizer, the maximum skip distance between
words. The function will compute all skip n-grams between \code{0} and
\code{k}.}
}
\value{
A list of character vectors containing the tokens, with one element
  in the list for each element that was passed as input. If \code{simplify =
  TRUE} and only a single element was passed as input, then the output is a
  character vector of tokens.
}
\description{
These functions tokenize their inputs into different kinds of n-grams. The
input can be a character vector of any length, or a list of character vectors
where each character vector in the list has a length of 1. See details for an
explanation of what each function does.
}
\details{
\describe{ \item{\code{tokenize_ngrams}:}{ Basic shingled n-grams. A
contiguous subsequence of \code{n} words. This will compute shingled n-grams
for every value of between \code{n_min} (which must be at least 1) and
\code{n}. } \item{\code{tokenize_skip_ngrams}:}{Skip n-grams. A subsequence
of \code{n} words which are at most a gap of \code{k} words between them. The
skip n-grams will be calculated for all values from \code{0} to \code{k}. } }

These functions will strip all punctuation and normalize all whitespace to a
single space character.
}
\examples{
song <-  paste0("How many roads must a man walk down\n",
                "Before you call him a man?\n",
                "How many seas must a white dove sail\n",
                "Before she sleeps in the sand?\n",
                "\n",
                "How many times must the cannonballs fly\n",
                "Before they're forever banned?\n",
                "The answer, my friend, is blowin' in the wind.\n",
                "The answer is blowin' in the wind.\n")

tokenize_ngrams(song, n = 4)
tokenize_ngrams(song, n = 4, n_min = 1)
tokenize_skip_ngrams(song, n = 4, k = 2)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chunk-text.R
\name{chunk_text}
\alias{chunk_text}
\title{Chunk text into smaller segments}
\usage{
chunk_text(x, chunk_size = 100, doc_id = names(x), ...)
}
\arguments{
\item{x}{A character vector or a list of character vectors to be tokenized
into n-grams. If \code{x} is a character vector, it can be of any length,
and each element will be chunked separately. If \code{x} is a list of
character vectors, each element of the list should have a length of 1.}

\item{chunk_size}{The number of words in each chunk.}

\item{doc_id}{The document IDs as a character vector. This will be taken from
the names of the \code{x} vector if available. \code{NULL} is acceptable.}

\item{...}{Arguments passed on to \code{\link{tokenize_words}}.}
}
\description{
Given a text or vector/list of texts, break the texts into smaller segments
each with the same number of words. This allows you to treat a very long
document, such as a novel, as a set of smaller documents.
}
\details{
Chunking the text passes it through \code{\link{tokenize_words}},
  which will strip punctuation and lowercase the text unless you provide
  arguments to pass along to that function.
}
\examples{
\dontrun{
chunked <- chunk_text(mobydick, chunk_size = 100)
length(chunked)
chunked[1:3]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tokenizers-package.r
\docType{package}
\name{tokenizers}
\alias{tokenizers}
\title{Tokenizers}
\description{
A collection of functions with a consistent interface to convert natural
language text into tokens.
}
\details{
The tokenizers in this package have a consistent interface. They all take
either a character vector of any length, or a list where each element is a
character vector of length one. The idea is that each element comprises a
text. Then each function returns a list with the same length as the input
vector, where each element in the list are the tokens generated by the
function. If the input character vector or list is named, then the names are
preserved.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/character-shingles-tokenizers.R
\name{tokenize_character_shingles}
\alias{tokenize_character_shingles}
\title{Character shingle tokenizers}
\usage{
tokenize_character_shingles(
  x,
  n = 3L,
  n_min = n,
  lowercase = TRUE,
  strip_non_alphanum = TRUE,
  simplify = FALSE
)
}
\arguments{
\item{x}{A character vector or a list of character vectors to be tokenized
into character shingles. If \code{x} is a character vector, it can be of
any length, and each element will be tokenized separately. If \code{x} is a
list of character vectors, each element of the list should have a length of
1.}

\item{n}{The number of characters in each shingle. This must be an integer
greater than or equal to 1.}

\item{n_min}{This must be an integer greater than or equal to 1, and less
than or equal to \code{n}.}

\item{lowercase}{Should the characters be made lower case?}

\item{strip_non_alphanum}{Should punctuation and white space be stripped?}

\item{simplify}{\code{FALSE} by default so that a consistent value is
returned regardless of length of input. If \code{TRUE}, then an input with
a single element will return a character vector of tokens instead of a
list.}
}
\value{
A list of character vectors containing the tokens, with one element
  in the list for each element that was passed as input. If \code{simplify =
  TRUE} and only a single element was passed as input, then the output is a
  character vector of tokens.
}
\description{
The character shingle tokenizer functions like an n-gram tokenizer, except
the units that are shingled are characters instead of words. Options to the
function let you determine whether non-alphanumeric characters like
punctuation should be retained or discarded.
}
\examples{
x <- c("Now is the hour of our discontent")
tokenize_character_shingles(x)
tokenize_character_shingles(x, n = 5)
tokenize_character_shingles(x, n = 5, strip_non_alphanum = FALSE)
tokenize_character_shingles(x, n = 5, n_min = 3, strip_non_alphanum = FALSE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-docs.R
\docType{data}
\name{mobydick}
\alias{mobydick}
\title{The text of Moby Dick}
\format{
A named character vector with length 1.
}
\source{
\url{http://www.gutenberg.org/}
}
\usage{
mobydick
}
\description{
The text of Moby Dick, by Herman Melville, taken from Project Gutenberg.
}
\keyword{datasets}
