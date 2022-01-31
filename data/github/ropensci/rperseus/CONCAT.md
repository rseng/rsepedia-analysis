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
rperseus
--------

------------------------------------------------------------------------

[![Build Status](https://travis-ci.org/ropensci/rperseus.svg?branch=master)](https://travis-ci.org/ropensci/rperseus) [![codecov](https://codecov.io/gh/ropensci/rperseus/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rperseus) [![](https://badges.ropensci.org/145_status.svg)](https://github.com/ropensci/onboarding/issues/145)

![](http://www.infobiblio.es/wp-content/uploads/2015/06/perseus-logo.png)

Author: David Ranzolin

License: MIT

Goal
----

The goal of `rperseus` is to furnish classicists, textual critics, and R enthusiasts with texts from the Classical World. While the English translations of most texts are available through `gutenbergr`, `rperseus`returns these works in their original language--Greek, Latin, and Hebrew.

Description
-----------

`rperseus` provides access to classical texts within the [Perseus Digital Library's](http://www.perseus.tufts.edu/hopper/) CapiTainS environment. A wealth of Greek, Latin, and Hebrew texts are available, from Homer to Cicero to Boetheius. The Perseus Digital Library includes English translations in some cases. The base API url is `http://cts.perseids.org/api/cts`.

Installation
------------

`rperseus` is not on CRAN, but can be installed via:

``` r
devtools::install_github("ropensci/rperseus")
```

Usage
-----

[See the vignette to get started.](https://daranzolin.github.io/rperseus//articles/rperseus-vignette.html)

To obtain a particular text, you must first know its full Uniform Resource Name (URN). URNs can be perused in the `perseus_catalog`, a data frame lazily loaded into the package. For example, say I want a copy of Virgil's *Aeneid*:

``` r
library(dplyr)
library(purrr)
library(rperseus)

aeneid_latin <- perseus_catalog %>% 
  filter(group_name == "Virgil",
         label == "Aeneid",
         language == "lat") %>% 
  pull(urn) %>% 
  get_perseus_text()
```

You can also request an English translation for some texts:

``` r
aeneid_english <- perseus_catalog %>% 
  filter(group_name == "Virgil",
         label == "Aeneid",
         language == "eng") %>% 
  pull(urn) %>% 
  get_perseus_text()
```

Refer to the language variable in `perseus_catalog` for translation availability.

Excerpts
--------

You can also specify excerpts:

``` r
qoheleth <- get_perseus_text(urn = "urn:cts:ancJewLit:hebBible.ecclesiastes.leningrad-pntd", excerpt = "1.1-1.3")
qoheleth$text
#> [1] "דִּבְרֵי֙ קֹהֶ֣לֶת בֶּן־ דָּוִ֔ד מֶ֖לֶךְ בִּירוּשָׁלִָֽם : הֲבֵ֤ל הֲבָלִים֙ אָמַ֣ר קֹהֶ֔לֶת הֲבֵ֥ל הֲבָלִ֖ים הַכֹּ֥ל הָֽבֶל : מַה־ יִּתְר֖וֹן לָֽאָדָ֑ם בְּכָל־ עֲמָל֔וֹ שֶֽׁיַּעֲמֹ֖ל תַּ֥חַת הַשָּֽׁמֶשׁ :"
```

Parsing Excerpts
----------------

You can parse any Greek excerpt, returning a data frame with each word's part of speech, gender, case, mood, voice, tense, person, number, and degree.

``` r
parse_excerpt("urn:cts:greekLit:tlg0031.tlg002.perseus-grc2", "5.1-5.2") %>% 
  head(7) %>% 
  knitr::kable()
```

| word    | form     | verse | part\_of\_speech | person | number   | tense  | mood       | voice  | gender   | case       | degree |
|:--------|:---------|:------|:-----------------|:-------|:---------|:-------|:-----------|:-------|:---------|:-----------|:-------|
| καί     | Καὶ      | 5.1   | conjunction      | NA     | NA       | NA     | NA         | NA     | NA       | NA         | NA     |
| ἔρχομαι | ἦλθον    | 5.1   | verb             | third  | plural   | aorist | indicative | active | NA       | NA         | NA     |
| εἰς     | εἰς      | 5.1   | preposition      | NA     | NA       | NA     | NA         | NA     | NA       | NA         | NA     |
| ὁ       | τὸ       | 5.1   | article          | NA     | singular | NA     | NA         | NA     | neuter   | accusative | NA     |
| πέραν   | πέραν    | 5.1   | adverb           | NA     | NA       | NA     | NA         | NA     | NA       | NA         | NA     |
| ὁ       | τῆς      | 5.1   | article          | NA     | singular | NA     | NA         | NA     | feminine | genative   | NA     |
| θάλασσα | θαλάσσης | 5.1   | noun             | NA     | singular | NA     | NA         | NA     | feminine | genative   | NA     |

tidyverse and tidytext
----------------------

`rperseus` plays well with the `tidyverse` and `tidytext`. Here I obtain all of Plato's works that have English translations available:

``` r
library(purrr)
plato <- perseus_catalog %>% 
  filter(group_name == "Plato",
         language == "eng") %>% 
  pull(urn) %>% 
  map_df(get_perseus_text)
```

And here's how to retrieve the Greek text from Sophocles' underrated *Philoctetes* before unleashing the `tidytext` toolkit:

``` r
library(tidytext)

philoctetes <- perseus_catalog %>% 
  filter(group_name == "Sophocles",
         label == "Philoctetes",
         language == "grc") %>% 
  pull(urn) %>%
  get_perseus_text()

philoctetes %>% 
  unnest_tokens(word, text) %>% 
  count(word, sort = TRUE) %>% 
  anti_join(greek_stop_words)
#> Joining, by = "word"
#> # A tibble: 3,514 x 2
#>           word     n
#>          <chr> <int>
#>  1 νεοπτόλεμος   164
#>  2  φιλοκτήτης   141
#>  3           ὦ   119
#>  4          μʼ    74
#>  5    ὀδυσσεύς    56
#>  6      τέκνον    47
#>  7          τʼ    43
#>  8       χορός    41
#>  9          γʼ    40
#> 10         νῦν    39
#> # ... with 3,504 more rows
```

Rendering Parallels
-------------------

You can render small parallels with `perseus_parallel`:

``` r
tibble(label = c("Colossians", "1 Thessalonians", "Romans"),
              excerpt = c("1.4", "1.3", "8.35-8.39")) %>%
    left_join(perseus_catalog) %>%
    filter(language == "grc") %>%
    select(urn, excerpt) %>%
    pmap_df(get_perseus_text) %>%
    perseus_parallel(words_per_row = 4)
#> Joining, by = "label"
```

![](README-unnamed-chunk-9-1.png)

Meta
----

-   [Report bugs or issues here.](https://github.com/daranzolin/rperseus/issues)
-   If you'd like to contribute to the development of `rperseus`, first get acquainted with the Perseus Digital Library, fork the repo, and send a pull request.
-   This project is released with a [Contributor Code of Conduct.](https://github.com/daranzolin/rperseus/blob/master/CONDUCT.md) By participating in this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# rperseus 0.1.2

* Transferred repo ownership to rOpenSci
* Includes `perseus_parallel`, a function to render parallels with `ggplot2`
* Includes `parse_excerpt`, a function to parse any Greek excerpt
* Includes `greek_stop_words`, a dictionary of Greek stop words
* Updated README, vignette

# rperseus 0.1

First version of package, including
  * `get_perseus_text` function, for getting texts from the Perseus Digital Library
  * `perseus_catalog` data frame to peruse text URNs
  * Introductory vignette including basic examples of retrieving texts
  * Added a `NEWS.md` file to track changes to the package
  * `pkgdown` site
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

## rperseus

***

[![Build Status](https://travis-ci.org/ropensci/rperseus.svg?branch=master)](https://travis-ci.org/ropensci/rperseus)
[![codecov](https://codecov.io/gh/ropensci/rperseus/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rperseus)
[![](https://badges.ropensci.org/145_status.svg)](https://github.com/ropensci/onboarding/issues/145)


![](http://www.infobiblio.es/wp-content/uploads/2015/06/perseus-logo.png)

Author: David Ranzolin

License: MIT

## Goal

The goal of `rperseus` is to furnish classicists, textual critics, and R enthusiasts with texts from the Classical World. While the English translations of most texts are available through `gutenbergr`, `rperseus`returns these works in their original language--Greek, Latin, and Hebrew.

## Description

`rperseus` provides access to classical texts within the [Perseus Digital Library's](http://www.perseus.tufts.edu/hopper/) CapiTainS environment. A wealth of Greek, Latin, and Hebrew texts are available, from Homer to Cicero to Boetheius. The Perseus Digital Library includes English translations in some cases. The base API url is `http://cts.perseids.org/api/cts`. 

## Installation

`rperseus` is not on CRAN, but can be installed via:

```{r eval = FALSE}
devtools::install_github("ropensci/rperseus")
```

## Usage

[See the vignette to get started.](https://daranzolin.github.io/rperseus//articles/rperseus-vignette.html)

To obtain a particular text, you must first know its full Uniform Resource Name (URN). URNs can be perused in the `perseus_catalog`, a data frame lazily loaded into the package. For example, say I want a copy of Virgil's *Aeneid*:

```{r warning = FALSE, message=FALSE}
library(dplyr)
library(purrr)
library(rperseus)

aeneid_latin <- perseus_catalog %>% 
  filter(group_name == "Virgil",
         label == "Aeneid",
         language == "lat") %>% 
  pull(urn) %>% 
  get_perseus_text()
```

You can also request an English translation for some texts:

```{r eval = FALSE}
aeneid_english <- perseus_catalog %>% 
  filter(group_name == "Virgil",
         label == "Aeneid",
         language == "eng") %>% 
  pull(urn) %>% 
  get_perseus_text()
```

Refer to the language variable in `perseus_catalog` for translation availability.

## Excerpts

You can also specify excerpts:

```{r}
qoheleth <- get_perseus_text(urn = "urn:cts:ancJewLit:hebBible.ecclesiastes.leningrad-pntd", excerpt = "1.1-1.3")
qoheleth$text
```

## Parsing Excerpts

You can parse any Greek excerpt, returning a data frame with each word's part of speech, gender, case, mood, voice, tense, person, number, and degree.

```{r}
parse_excerpt("urn:cts:greekLit:tlg0031.tlg002.perseus-grc2", "5.1-5.2") %>% 
  head(7) %>% 
  knitr::kable()
```

## tidyverse and tidytext 

`rperseus` plays well with the `tidyverse` and `tidytext`. Here I obtain all of Plato's works that have English translations available:

```{r eval = FALSE, warning = FALSE}
library(purrr)
plato <- perseus_catalog %>% 
  filter(group_name == "Plato",
         language == "eng") %>% 
  pull(urn) %>% 
  map_df(get_perseus_text)
```

And here's how to retrieve the Greek text from Sophocles' underrated *Philoctetes* before unleashing the `tidytext` toolkit:

```{r warning = FALSE}
library(tidytext)

philoctetes <- perseus_catalog %>% 
  filter(group_name == "Sophocles",
         label == "Philoctetes",
         language == "grc") %>% 
  pull(urn) %>%
  get_perseus_text()

philoctetes %>% 
  unnest_tokens(word, text) %>% 
  count(word, sort = TRUE) %>% 
  anti_join(greek_stop_words)
```


## Rendering Parallels

You can render small parallels with `perseus_parallel`:

```{r,fig.width=8, fig.height=6}
tibble(label = c("Colossians", "1 Thessalonians", "Romans"),
              excerpt = c("1.4", "1.3", "8.35-8.39")) %>%
    left_join(perseus_catalog) %>%
    filter(language == "grc") %>%
    select(urn, excerpt) %>%
    pmap_df(get_perseus_text) %>%
    perseus_parallel(words_per_row = 4)
```

## Meta

* [Report bugs or issues here.](https://github.com/daranzolin/rperseus/issues)
* If you'd like to contribute to the development of `rperseus`, first get acquainted with the Perseus Digital Library, fork the repo, and send a pull request.
* This project is released with a [Contributor Code of Conduct.](https://github.com/daranzolin/rperseus/blob/master/CONDUCT.md) By participating in this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)


---
title: "rperseus Vignette"
author: "David Ranzolin"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

### Introduction

`rperseus` provides tools to get and analyze classical texts. Version 0.1.2 includes:

* `get_perseus_text`, a function to obtain a text from the Perseus Digital Library
* `perseus_parallel`, a function to render a text parallel in `ggplot2`
* `parse_excerpt`, a function to parse any Greek excerpt
* `perseus_catalog`, a data frame of available texts
* `greek_stop_words`, a data frame of Greek pronouns, articles, prepositions, and particles

```{r}
library(rperseus)
head(perseus_catalog)
```

A snapshot of available authors:

```{r warning=FALSE, message=FALSE}
library(dplyr)
count(perseus_catalog, group_name, sort = TRUE)
```

### Getting a Text

Once you've identified the relevant URN, paste it into a call to `get_perseus_text`. Here I've called for the Greek text of Plato's *Crito*:

```{r}
crito <- get_perseus_text(urn = "urn:cts:greekLit:tlg0059.tlg003.perseus-grc2")
crito$text[1]
```

### Getting Multiple Texts with the tidyverse

You can collect all of Plato's available English translations with the `tidyverse:`

```{r eval = FALSE}
plato <- perseus_catalog %>% 
  filter(group_name == "Plato",
         language == "eng") %>% 
  pull(urn) %>% 
  map_df(get_perseus_text)
```

### Rendering Parallels

You can render small parallels with `perseus_parallel`:

```{r,fig.width=8, fig.height=6}
tibble::tibble(label = c("Colossians", "1 Thessalonians", "Romans"),
              excerpt = c("1.4", "1.3", "8.35-8.39")) %>%
    dplyr::left_join(perseus_catalog) %>%
    dplyr::filter(language == "grc") %>%
    dplyr::select(urn, excerpt) %>%
    as.list() %>%
    purrr::pmap_df(get_perseus_text) %>%
    perseus_parallel(words_per_row = 4)
```

### Parsing Excerpts

You can parse any Greek excerpt with `parse_excerpt`. A data frame is returned including part of speech, person, number, tense, mood, voice, gender, case, and degree.

```{r}
parse_excerpt("urn:cts:greekLit:tlg0031.tlg002.perseus-grc2", "5.1-5.3")
```




% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_excerpt.R
\name{parse_excerpt}
\alias{parse_excerpt}
\title{Parse a Greek excerpt}
\usage{
parse_excerpt(urn, excerpt)
}
\arguments{
\item{urn}{a valid urn from the perseus_catalog object.}

\item{excerpt}{a valid excerpt, e.g. 5.1-5.5}
}
\value{
a data frame
}
\description{
This function parses a Greek excerpt from the Perseus Digital Library. Parsing includes
part of speech, gender, case, mood, voice, tense, person, number, and degree.
}
\examples{
parse_excerpt("urn:cts:greekLit:tlg0031.tlg002.perseus-grc2", "5.1-5.4")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_perseus_text.R
\name{get_perseus_text}
\alias{get_perseus_text}
\title{Get a primary text by URN.}
\usage{
get_perseus_text(urn, excerpt = NULL)
}
\arguments{
\item{urn}{Valid uniform resource number (URN) obtained from \code{\link{perseus_catalog}}.}

\item{excerpt}{An index to excerpt the text. For example, the first four "verses" of a text might be 1.1-1.4. If NULL, the entire work is returned.}
}
\value{
A seven column \code{tbl_df} with one row for each "section" (splits vary from text--could be line, chapter, etc.).
Columns:
\describe{
  \item{text}{character vector of text}
  \item{urn}{Uniform Resource Number}
  \item{group_name}{Could refer to author (e.g. "Aristotle") or corpus (e.g. "New Testament")}
  \item{label}{Text label, e.g. "Phaedrus"}
  \item{description}{Text description}
  \item{language}{Text language, e.g. "grc" = Greek, "lat" = Latin, "eng" = English}
  \item{section}{Text section or excerpt if specified}
}
}
\description{
Get a primary text by URN.
}
\examples{
get_perseus_text("urn:cts:greekLit:tlg0013.tlg028.perseus-grc2")
get_perseus_text("urn:cts:latinLit:stoa0215b.stoa003.opp-lat1")
get_perseus_text(urn = "urn:cts:greekLit:tlg0031.tlg009.perseus-grc2", excerpt = "5.1-5.5")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{greek_stop_words}
\alias{greek_stop_words}
\title{A dictionary of Greek stop words}
\format{A data frame with 223 rows and 1 variable:
\describe{
  \item{word}{Greek stop word}
}}
\source{
\url{Compiled manually by filtering prepositions, pronouns, conjunctions, particles, and articles.}
}
\usage{
greek_stop_words
}
\description{
A dictionary of Greek stop words
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rperseus.R
\docType{package}
\name{rperseus}
\alias{rperseus}
\alias{rperseus-package}
\title{rperseus: A package to obtain texts from the Perseus Digital Library}
\description{
The package contains a catalog and a function. The catalog is \code{\link{perseus_catalog}} and can be
perused to locate a particular text and its corresponding URN. Then pass the URN to \code{\link{get_perseus_text}}
to obtain the text.
}
\author{
David Ranzolin \email{daranzolin@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perseus_parallel.R
\name{perseus_parallel}
\alias{perseus_parallel}
\title{Render a text parallel with ggplot2}
\usage{
perseus_parallel(perseus_df, words_per_row = 6)
}
\arguments{
\item{perseus_df}{a data frame obtained from \code{get_perseus_text}. Can contain multiple texts.}

\item{words_per_row}{adjusts the words displayed per "row".}
}
\value{
a ggplot object
}
\description{
Render a text parallel with ggplot2
}
\examples{
\dontrun{
tibble::tibble(label = c("Colossians", rep("1 Thessalonians", 2), "Romans"),
               excerpt = c("1.4", "1.3", "5.8", "8.35-8.39")) \%>\%
 dplyr::left_join(perseus_catalog) \%>\%
 dplyr::filter(language == "grc") \%>\%
 dplyr::select(urn, excerpt) \%>\%
 as.list() \%>\%
 purrr::pmap_df(get_perseus_text) \%>\%
 perseus_parallel()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{perseus_catalog}
\alias{perseus_catalog}
\title{Metadata for texts available via the Perseus Digital Library.}
\format{A data frame with 2291 rows and 5 variables:
\describe{
  \item{urn}{Uniform Resource Number}
  \item{group_name}{Could refer to author (e.g. "Aristotle") or corpus (e.g. "New Testament")}
  \item{label}{Text label, e.g. "Phaedrus"}
  \item{description}{Text description}
  \item{language}{Text language, e.g. "grc" = Greek, "lat" = Latin, "eng" = English, "hpt" = Hebrew pointed text, "hct" = Hebrew consonantal text, "ger" = German, "oth" = other}
}}
\source{
\url{http://cts.perseids.org/api/cts/?request=GetCapabilities}
}
\usage{
perseus_catalog
}
\description{
A dataset containing the texts available from the Perseus Digital Library.
}
\keyword{datasets}
