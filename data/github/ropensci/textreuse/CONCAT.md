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

textreuse
=========

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/textreuse)](https://cran.r-project.org/package=textreuse)
[![CRAN\_Downloads](http://cranlogs.r-pkg.org/badges/grand-total/textreuse)](https://cran.r-project.org/package=textreuse)
[![Build
Status](https://travis-ci.org/ropensci/textreuse.svg?branch=master)](https://travis-ci.org/ropensci/textreuse)
[![Build
status](https://ci.appveyor.com/api/projects/status/9qwf0473xi8cyuoh/branch/master?svg=true)](https://ci.appveyor.com/project/lmullen/textreuse-6xljc/branch/master)
[![Coverage
Status](https://img.shields.io/codecov/c/github/ropensci/textreuse/master.svg)](https://codecov.io/github/ropensci/textreuse?branch=master)
[![rOpenSci
badge](https://badges.ropensci.org/20_status.svg)](https://github.com/ropensci/onboarding/issues/20)

Overview
--------

This [R](https://www.r-project.org/) package provides a set of functions
for measuring similarity among documents and detecting passages which
have been reused. It implements shingled n-gram, skip n-gram, and other
tokenizers; similarity/dissimilarity functions; pairwise comparisons;
minhash and locality sensitive hashing algorithms; and a version of the
Smith-Waterman local alignment algorithm suitable for natural language.
It is broadly useful for, for example, detecting duplicate documents in
a corpus prior to text analysis, or for identifying borrowed passages
between texts. The classes provides by this package follow the model of
other natural language processing packages for R, especially the
[NLP](https://cran.r-project.org/package=NLP) and
[tm](https://cran.r-project.org/package=tm) packages. (However, this
package has no dependency on Java, which should make it easier to
install.)

### Citation

If you use this package for scholarly research, I would appreciate a
citation.

``` r
citation("textreuse")
#> 
#> To cite package 'textreuse' in publications use:
#> 
#>   Lincoln Mullen (2020). textreuse: Detect Text Reuse and Document
#>   Similarity. https://docs.ropensci.org/textreuse,
#>   https://github.com/ropensci/textreuse.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {textreuse: Detect Text Reuse and Document Similarity},
#>     author = {Lincoln Mullen},
#>     year = {2020},
#>     note = {https://docs.ropensci.org/textreuse, https://github.com/ropensci/textreuse},
#>   }
```

Installation
------------

To install this package from CRAN:

``` r
install.packages("textreuse")
```

To install the development version from GitHub, use
[devtools](https://github.com/hadley/devtools).

``` r
# install.packages("devtools")
devtools::install_github("ropensci/textreuse", build_vignettes = TRUE)
```

Examples
--------

There are three main approaches that one may take when using this
package: pairwise comparisons, minhashing/locality sensitive hashing,
and extracting matching passages through text alignment.

See the [introductory
vignette](https://cran.r-project.org/package=textreuse/vignettes/textreuse-introduction.html)
for a description of the classes provided by this package.

``` r
vignette("textreuse-introduction", package = "textreuse")
```

### Pairwise comparisons

In this example we will load a tiny corpus of three documents. These
documents are drawn from Kellen Funk’s
[research](http://kellenfunk.org/field-code/) into the propagation of
legal codes of civil procedure in the nineteenth-century United States.

``` r
library(textreuse)
dir <- system.file("extdata/legal", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, meta = list(title = "Civil procedure"),
                          tokenizer = tokenize_ngrams, n = 7)
```

We have loaded the three documents into a corpus, which involves
tokenizing the text and hashing the tokens. We can inspect the corpus as
a whole or the individual documents that make it up.

``` r
corpus
#> TextReuseCorpus
#> Number of documents: 3 
#> hash_func : hash_string 
#> title : Civil procedure 
#> tokenizer : tokenize_ngrams
names(corpus)
#> [1] "ca1851-match"   "ca1851-nomatch" "ny1850-match"
corpus[["ca1851-match"]]
#> TextReuseTextDocument
#> file : /Users/lmullen/R/library/textreuse/extdata/legal/ca1851-match.txt 
#> hash_func : hash_string 
#> id : ca1851-match 
#> minhash_func : 
#> tokenizer : tokenize_ngrams 
#> content : § 4. Every action shall be prosecuted in the name of the real party
#> in interest, except as otherwise provided in this Act.
#> 
#> § 5. In the case of an assignment of a thing in action, the action by
#> the as
```

Now we can compare each of the documents to one another. The
`pairwise_compare()` function applies a comparison function (in this
case, `jaccard_similarity()`) to every pair of documents. The result is
a matrix of scores. As we would expect, some documents are similar and
others are not.

``` r
comparisons <- pairwise_compare(corpus, jaccard_similarity)
comparisons
#>                ca1851-match ca1851-nomatch ny1850-match
#> ca1851-match             NA              0    0.3842549
#> ca1851-nomatch           NA             NA    0.0000000
#> ny1850-match             NA             NA           NA
```

We can convert that matrix to a data frame of pairs and scores if we
prefer.

``` r
pairwise_candidates(comparisons)
#> # A tibble: 3 x 3
#>   a              b              score
#> * <chr>          <chr>          <dbl>
#> 1 ca1851-match   ca1851-nomatch 0    
#> 2 ca1851-match   ny1850-match   0.384
#> 3 ca1851-nomatch ny1850-match   0
```

See the [pairwise
vignette](https://cran.r-project.org/package=textreuse/vignettes/textreuse-pairwise.html)
for a fuller description.

``` r
vignette("textreuse-pairwise", package = "textreuse")
```

### Minhashing and locality sensitive hashing

Pairwise comparisons can be very time-consuming because they grow
geometrically with the size of the corpus. (A corpus with 10 documents
would require at least 45 comparisons; a corpus with 100 documents would
require 4,950 comparisons; a corpus with 1,000 documents would require
499,500 comparisons.) That’s why this package implements the minhash and
locality sensitive hashing algorithms, which can detect candidate pairs
much faster than pairwise comparisons in corpora of any significant
size.

For this example we will load a small corpus of ten documents published
by the American Tract Society. We will also create a minhash function,
which represents an entire document (regardless of length) by a fixed
number of integer hashes. When we create the corpus, the documents will
each have a minhash signature.

``` r
dir <- system.file("extdata/ats", package = "textreuse")
minhash <- minhash_generator(200, seed = 235)
ats <- TextReuseCorpus(dir = dir,
                       tokenizer = tokenize_ngrams, n = 5,
                       minhash_func = minhash)
```

Now we can calculate potential matches, extract the candidates, and
apply a comparison function to just those candidates.

``` r
buckets <- lsh(ats, bands = 50, progress = FALSE)
candidates <- lsh_candidates(buckets)
scores <- lsh_compare(candidates, ats, jaccard_similarity, progress = FALSE)
scores
#> # A tibble: 2 x 3
#>   a                     b                      score
#>   <chr>                 <chr>                  <dbl>
#> 1 practicalthought00nev thoughtsonpopery00nevi 0.463
#> 2 remember00palm        remembermeorholy00palm 0.701
```

For details, see the [minhash
vignette](https://cran.r-project.org/package=textreuse/vignettes/textreuse-minhash.html).

``` r
vignette("textreuse-minhash", package = "textreuse")
```

### Text alignment

We can also extract the optimal alignment between to documents with a
version of the
[Smith-Waterman](https://en.wikipedia.org/wiki/Smith-Waterman_algorithm)
algorithm, used for protein sequence alignment, adapted for natural
language. The longest matching substring according to scoring values
will be extracted, and variations in the alignment will be marked.

``` r
a <- "'How do I know', she asked, 'if this is a good match?'"
b <- "'This is a match', he replied."
align_local(a, b)
#> TextReuse alignment
#> Alignment score: 7 
#> Document A:
#> this is a good match
#> 
#> Document B:
#> This is a #### match
```

For details, see the [text alignment
vignette](https://cran.r-project.org/package=textreuse/vignettes/textreuse-alignment.html).

``` r
vignette("textreuse-alignment", package = "textreuse")
```

### Parallel processing

Loading the corpus and creating tokens benefit from using multiple
cores, if available. (This works only on non-Windows machines.) To use
multiple cores, set `options("mc.cores" = 4L)`, where the number is how
many cores you wish to use.

### Contributing and acknowledgments

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/textreuse/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.

Thanks to [Noam Ross](http://www.noamross.net/) for his thorough [peer
review](https://github.com/ropensci/onboarding/issues/20) of this
package for [rOpenSci](https://ropensci.org/).

------------------------------------------------------------------------

[![rOpenSCi
logo](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# textreuse 0.1.5

- Updates due to dplyr 1.0.0 release

# textreuse 0.1.4

- Preventative maintenance release to avoid failing tests when new version of
  BH is released.

# textreuse 0.1.3

- Preventative maintenance release to avoid failing tests when new versions of 
  the dplyr and testthat packages are released.

# textreuse 0.1.2

- Fix memory error in `shingle_ngrams()`
- Fix tests for retokenizing on Windows
- More informative error message if using `lsh()` on corpora without minhashes

# textreuse 0.1.1

- Fix progress bars in vignettes

# textreuse 0.1.0

- Initial release
This is a maintenance release to handle changes in the forthcoming
dplyr 1.0.0.

This resubmission handles a problem with a URL to CONDUCT.md.

## Test environments

* local OS X 10.15.4 install: R-release
* Ubuntu 18.04 (on Travis-CI): R-devel, R-release, R-oldrel
* Win-builder: R-devel, R-release

## R CMD check results

There were two NOTEs. 

There was one NOTE about this being a updated version of the package. The NOTE mentions three possible misspelled words, all of which are correctly spelled proper nouns or terms of art.

Win-builder reported one NOTE about large components, but these are just some sample text files for testing and demonstration purposes.
---
output:
  md_document:
    variant: markdown_github
pagetitle: Detect Text Reuse and Document Similarity
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
suppressPackageStartupMessages(library(dplyr))
```

# textreuse

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/textreuse)](https://cran.r-project.org/package=textreuse)
[![CRAN_Downloads](http://cranlogs.r-pkg.org/badges/grand-total/textreuse)](https://cran.r-project.org/package=textreuse)
[![Build Status](https://travis-ci.org/ropensci/textreuse.svg?branch=master)](https://travis-ci.org/ropensci/textreuse) 
[![Build status](https://ci.appveyor.com/api/projects/status/9qwf0473xi8cyuoh/branch/master?svg=true)](https://ci.appveyor.com/project/lmullen/textreuse-6xljc/branch/master)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/textreuse/master.svg)](https://codecov.io/github/ropensci/textreuse?branch=master)
[![rOpenSci badge](https://badges.ropensci.org/20_status.svg)](https://github.com/ropensci/onboarding/issues/20)

## Overview

This [R](https://www.r-project.org/) package provides a set of functions for measuring similarity among documents and detecting passages which have been reused. It implements shingled n-gram, skip n-gram, and other tokenizers; similarity/dissimilarity functions; pairwise comparisons; minhash and locality sensitive hashing algorithms; and a version of the Smith-Waterman local alignment algorithm suitable for natural language. It is broadly useful for, for example, detecting duplicate documents in a corpus prior to text analysis, or for identifying borrowed passages between texts. The classes provides by this package follow the model of other natural language processing packages for R, especially the [NLP](https://cran.r-project.org/package=NLP) and [tm](https://cran.r-project.org/package=tm) packages. (However, this package has no dependency on Java, which should make it easier to install.)

### Citation

If you use this package for scholarly research, I would appreciate a citation.

```{r}
citation("textreuse")
```

## Installation

To install this package from CRAN:

```{r eval=FALSE}
install.packages("textreuse")
```

To install the development version from GitHub, use [devtools](https://github.com/hadley/devtools).  

```{r eval=FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/textreuse", build_vignettes = TRUE)
```

## Examples

There are three main approaches that one may take when using this package: pairwise comparisons, minhashing/locality sensitive hashing, and extracting matching passages through text alignment.

See the [introductory vignette](https://cran.r-project.org/package=textreuse/vignettes/textreuse-introduction.html) for a description of the classes provided by this package.

```{r eval = FALSE}
vignette("textreuse-introduction", package = "textreuse")
```

### Pairwise comparisons

In this example we will load a tiny corpus of three documents. These documents are drawn from Kellen Funk's [research](http://kellenfunk.org/field-code/) into the propagation of legal codes of civil procedure in the nineteenth-century United States.

```{r}
library(textreuse)
dir <- system.file("extdata/legal", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, meta = list(title = "Civil procedure"),
                          tokenizer = tokenize_ngrams, n = 7)
```

We have loaded the three documents into a corpus, which involves tokenizing the text and hashing the tokens. We can inspect the corpus as a whole or the individual documents that make it up.

```{r}
corpus
names(corpus)
corpus[["ca1851-match"]]
```

Now we can compare each of the documents to one another. The `pairwise_compare()` function applies a comparison function (in this case, `jaccard_similarity()`) to every pair of documents. The result is a matrix of scores. As we would expect, some documents are similar and others are not.

```{r}
comparisons <- pairwise_compare(corpus, jaccard_similarity)
comparisons
```

We can convert that matrix to a data frame of pairs and scores if we prefer.

```{r}
pairwise_candidates(comparisons)
```

See the [pairwise vignette](https://cran.r-project.org/package=textreuse/vignettes/textreuse-pairwise.html) for a fuller description.

```{r eval=FALSE}
vignette("textreuse-pairwise", package = "textreuse")
```

### Minhashing and locality sensitive hashing

Pairwise comparisons can be very time-consuming because they grow geometrically with the size of the corpus. (A corpus with 10 documents would require at least 45 comparisons; a corpus with 100 documents would require 4,950 comparisons; a corpus with 1,000 documents would require 499,500 comparisons.) That's why this package implements the minhash and locality sensitive hashing algorithms, which can detect candidate pairs much faster than pairwise comparisons in corpora of any significant size. 

For this example we will load a small corpus of ten documents published by the American Tract Society. We will also create a minhash function, which represents an entire document (regardless of length) by a fixed number of integer hashes. When we create the corpus, the documents will each have a minhash signature.

```{r}
dir <- system.file("extdata/ats", package = "textreuse")
minhash <- minhash_generator(200, seed = 235)
ats <- TextReuseCorpus(dir = dir,
                       tokenizer = tokenize_ngrams, n = 5,
                       minhash_func = minhash)
```

Now we can calculate potential matches, extract the candidates, and apply a comparison function to just those candidates.

```{r}
buckets <- lsh(ats, bands = 50, progress = FALSE)
candidates <- lsh_candidates(buckets)
scores <- lsh_compare(candidates, ats, jaccard_similarity, progress = FALSE)
scores
```

For details, see the [minhash vignette](https://cran.r-project.org/package=textreuse/vignettes/textreuse-minhash.html).

```{r eval=FALSE}
vignette("textreuse-minhash", package = "textreuse")
```

### Text alignment

We can also extract the optimal alignment between to documents with a version of the  [Smith-Waterman](https://en.wikipedia.org/wiki/Smith-Waterman_algorithm) algorithm, used for protein sequence alignment, adapted for natural language. The longest matching substring according to scoring values will be extracted, and variations in the alignment will be marked.

```{r}
a <- "'How do I know', she asked, 'if this is a good match?'"
b <- "'This is a match', he replied."
align_local(a, b)
```

For details, see the [text alignment vignette](https://cran.r-project.org/package=textreuse/vignettes/textreuse-alignment.html).

```{r eval=FALSE}
vignette("textreuse-alignment", package = "textreuse")
```

### Parallel processing

Loading the corpus and creating tokens benefit from using multiple cores, if available. (This works only on non-Windows machines.) To use multiple cores, set `options("mc.cores" = 4L)`, where the number is how many cores you wish to use.

### Contributing and acknowledgments

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/textreuse/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

Thanks to [Noam Ross](http://www.noamross.net/) for his thorough [peer review](https://github.com/ropensci/onboarding/issues/20) of this package for [rOpenSci](https://ropensci.org/).

------------------------------------------------------------------------

[![rOpenSCi logo](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "Text Alignment"
author: "Lincoln Mullen"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Text alignment}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Local alignment is the process of finding taking two documents and finding the best subset of each document that aligns with one another. A commonly used local alignment algorithm for genetics is the [Smith-Waterman algorithm](https://en.wikipedia.org/wiki/Smith-Waterman_algorithm). This package offers a version of the Smith-Waterman algorithm intended to be used for natural language processing.

Consider these two documents. The first is part of Shakespeare's *Measure for Measure*. The second is a made-up piece of literary criticism quoting the play, but our imaginary literary critic has bungled the quotation. This is a common class of problems (not bungling literary critics but) documents which contain pieces, often heavily modified, from other documents.

```{r}
shakespeare <- paste(
  "Haste still pays haste, and leisure answers leisure;",
  "Like doth quit like, and MEASURE still FOR MEASURE.",
  "Then, Angelo, thy fault's thus manifested;",
  "Which, though thou wouldst deny, denies thee vantage.",
  "We do condemn thee to the very block",
  "Where Claudio stoop'd to death, and with like haste.",
  "Away with him!")
critic <- paste(
  "The play comes to its culmination where Duke Vincentio, quoting from",
  "the words of the Sermon on the Mount, says,",
  "'Haste still goes very quickly , and leisure answers leisure;",
  "Like doth cancel like, and measure still for measure.'",
  "These titular words sum up the meaning of the play.")
```

We can uses the local alignment function to extract the part of the text that was borrowed. Notice that the resulting object shows us the changes that have been made.

```{r}
library(textreuse)
align_local(shakespeare, critic)
```

See the documentation for the function to see how to tune the match: `?align_local`. This function works with character vectors or with documents of class `TextReuseTextDocument`.

---
title: "Minhash and locality-sensitive hashing"
author: "Lincoln Mullen"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Minhash and locality-sensitive hashing}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo=FALSE, message=FALSE}
library("dplyr")
```

Performing pairwise comparisons in a corpus is time-consuming because the number of comparisons grows geometrically with the size of the corpus. Most of those comparisons, furthermore, are unnecessary because they do not result in matches. The combination of minhash and locality-sensitive hashing (LSH) seeks to solve these problems. They make it possible to compute possible matches only once for each document, so that the cost of computation grows linearly rather than exponentially. This vignette explains how to use the minhash and locality-sensitive hashing functions in this package. For an explanation of why they work, see Jure Leskovec, Anand Rajaraman, and Jeff Ullman, *[Mining of Massive Datasets](http://www.mmds.org/#book)* (Cambridge University Press, 2011), ch. 3. (This [blog post](http://matthewcasperson.blogspot.com/2013/11/minhash-for-dummies.html) is a more succinct explanation.)

We begin by creating a minhash function. A minhash function converts tokenized text into a set of hash integers, then selects the minimum value. This is the equivalent of randomly selecting a token. The function then does the same thing repeatedly with different hashing functions, in effect selecting `n` random shingles. The additional hashing functions come from a bitwise XOR with random integers. That is why the `minhash_generator()` accepts a seed, so that we can re-create the same minhash function again. In other words, a minhash function converts a set of tokens of any length into `n` randomly selected and hashed tokens.

```{r}
library(textreuse)
minhash <- minhash_generator(n = 240, seed = 3552)
head(minhash(c("turn tokens into", "tokens into hashes", "into hashes fast")))
```

Now when we load our corpus, we will tokenize our texts as usual, but we will use our generated `minhash()` function to compute the hashes. We specify that we want to create a minhash signature by passing our minhash function to the `minhash_func =` parameter.

```{r}
dir <- system.file("extdata/ats", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, tokenizer = tokenize_ngrams, n = 5,
                          minhash_func = minhash, keep_tokens = TRUE,
                          progress = FALSE)
```

We can verify that we have minhashes in our corpus:

```{r}
head(minhashes(corpus[[1]]))
length(minhashes(corpus[[1]]))
```


Now all our documents are represented by `n = 240` randomly selected and hashed shingles. Comparing those shingles should be the equivalent of finding the Jaccard similarity of the two documents. However, we still have the problem of pairwise comparison.

The locality-sensitive hashing algorithm, provided in this package by the `lsh()` function, solves this problem. LSH breaks the minhashes into a series of bands comprised of rows. For example, 200 minhashes might broken into 50 bands of 4 rows each. Each band is hashed to a bucket. If two documents have the exact same minhashes in a band, they will be hashed to the same bucket, and so will be considered candidate pairs. Each pair of documents has as many chances to be considered a candidate as their are bands, and the fewer rows there are in each band, the more likely it is that each document will match another.

How likely is it, then, that we will detect a match? The probability of a match depends on the Jaccard similarity of a pair of documents. The more similar two documents are, the more likely they are to be considered candidates, which is what we want. The probability of a match is an S-curve (see Leskovec, Rajaraman, and Ullman), so there is a threshold Jaccard similarity above which documents are likely to be a match. We can calculate the likely threshold based on the number of minhashes and bands that we are using.

```{r}
lsh_threshold(h = 200, b = 50)
lsh_threshold(h = 240, b = 80)
```

Using 240 minhashes and 80 bands, we will likely detect documents with an actual Jaccard similarity of above 0.232. We can also estimate the probability that a pair of documents with a Jaccard similarity `s` will be marked as potential matches.

```{r}
lsh_probability(h = 240, b = 80, s = 0.25)
lsh_probability(h = 240, b =  80, s = 0.75)
```

These numbers seem reasonable for our purposes, so we will set the number of minhashes at 240 and the number of bands at 80.

Now we can use the `lsh()` function to calculate the locality-sensitive hashes for our documents. 

```{r}
buckets <- lsh(corpus, bands = 80, progress = FALSE)
buckets
```

Note that using the LSH method only requires us to calculate the signatures (or buckets) for each document one time. This implies that we can take several data frames of LSH signatures and bind their rows together (e.g., with `dplyr::bind_rows()`). This permits us to compute the signatures for only part of a corpus at a time, or to continue to add to the corpus. Note, however, that you **must** use the same minhash function, generating the same number of minhashes and using the same seed and you **must** use the same number of bands in order to get valid results.

We can extract the potential matches from the cache using `lsh_query()` or  `lsh_candidates()`. The first function returns matches for only one document, specified by its ID; the second functions returns all potential pairs of matches.

```{r}
baxter_matches <- lsh_query(buckets, "calltounconv00baxt")
baxter_matches
candidates <- lsh_candidates(buckets)
candidates
```

Notice that LSH has identified the same three pairs of documents as potential matches that we found with pairwise comparisons, but did so much faster. But we do not have similarity scores; we only know that these documents are likely to have Jaccard similarity scores above the `r round(lsh_threshold(h = 240, b = 80), 3)` threshold.

Now we can use `lsh_compare()` to apply a similarity function to the candidate pairs of documents. Note that we only have to do 3 comparisons for all the candidates, instead of 28 pairs when comparing all 8 documents in the corpus pairwise.

```{r}
lsh_compare(candidates, corpus, jaccard_similarity, progress = FALSE)
```

Note that these results are identical to what we calculated in the pairwise vignette, but required much less computation.
---
title: "Pairwise comparisons for document similarity"
author: "Lincoln Mullen"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Pairwise comparisons for document similarity}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The most straightforward way to compare documents within a corpus is to compare each document to every other document. 

First we will load the corpus and tokenize it with shingled n-grams.

```{r}
library(textreuse)
dir <- system.file("extdata/ats", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, tokenizer = tokenize_ngrams, n = 5,
                          progress = FALSE)
```

We can use any of the comparison functions to compare two documents in the corpus. (Note that these functions, when applied to documents, compare their hashed tokens and not the tokens directly.)

```{r}
jaccard_similarity(corpus[["remember00palm"]], 
                   corpus[["remembermeorholy00palm"]])
```

The `pairwise_compare()` function applies a comparison function to each pair of documents in a corpus. The result is a matrix with the scores for each comparison.

```{r eval=FALSE}
comparisons <- pairwise_compare(corpus, jaccard_similarity, progress = FALSE)
comparisons[1:4, 1:4]
```

```{r, echo=FALSE}
comparisons <- pairwise_compare(corpus, jaccard_similarity, progress = FALSE)
round(comparisons[1:3, 1:3], digits = 3)
```

If you prefer, you can convert the matrix of all comparisons to a data frame of pairs and scores. Here we create the data frame and keep only the pairs with scores above a significant value.

```{r}
candidates <- pairwise_candidates(comparisons)
candidates[candidates$score > 0.1, ]
```

The pairwise comparison method is inadequate for a corpus of any size, however. For a corpus of size $n$, the number of comparisons (assuming the comparisons are commutative) is $\frac{n^2 - n}{2}$. A corpus of 100 documents would require 4,950 comparisons; a corpus of 1,000 documents would require 499,500 comparisons. A better approach for corpora of any appreciable size is to use the minhash/LSH algorithms described in another vignette: 

```{r eval=FALSE}
vignette("minhash", package = "textreuse")
```
---
title: "Introduction to the textreuse package"
author: "Lincoln Mullen"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to the textreuse packages}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The textreuse package provides classes and functions to detect document similarity and text reuse in text corpora. This introductory vignette provides details on the `TextReuseTextDocument` and `TextReuseCorpus` classes, as well as functions for tokenizing, hashing, and measuring similarity. See the pairwise, minhash/LSH, or alignment vignettes for details on solving text similarity problems.

```{r eval=FALSE}
vignette("textreuse-pairwise", package = "textreuse")
vignette("textreuse-minhash", package = "textreuse")
vignette("textreuse-alignment", package = "textreuse")
```

For these vignette we will use a small corpus of eight documents published by the [American Tract Society](https://en.wikipedia.org/wiki/American_Tract_Society) and available from the Internet Archive. The [full corpus](http://lincolnmullen.com/blog/corpus-of-american-tract-society-publications/) is also available to be downloaded if you wish to test the package.

## TextReuse classes

### TextReuseTextDocument

The most basic class provided by this package is the `TextReuseTextDocument` class. This class contains the text of a document and its metadata. When the document is loaded, the text is also tokenized. (See the section on tokenizers below.) Those tokens are then hashed using a hash function. By default the hashes are retained and the tokens are discarded, since using only hashes results in a significant memory savings. 

Here we load a file into a `TextReuseTextDocument` and tokenize it into shingled n-grams, adding an option to retain the tokens.

```{r}
library(textreuse)
file <- system.file("extdata/ats/remember00palm.txt", 
                    package = "textreuse")
doc <- TextReuseTextDocument(file = file, meta = list("publisher" = "ATS"),
                             tokenizer = tokenize_ngrams, n = 5,
                             keep_tokens = TRUE)
doc
```

We can see details of the document with accessor functions. These are derived from the S3 virtual class `TextDocument ` in the [NLP](https://cran.r-project.org/package=NLP) package. Notice that an ID has been assigned to the document based on the filename (without the extension). The name of the tokenizer and hash functions are also saved in the metadata.

```{r}
meta(doc)
meta(doc, "id")
meta(doc, "date") <- 1865
head(tokens(doc))
head(hashes(doc))
wordcount(doc)
```

The `tokens()` and `hashes()` function return the tokens and hashes associated with the document. The `meta()` function returns a named list of all the metadata fields. If that function is called with a specific ID, as in `meta(doc, "myfield")`, then the value for only that field is returned. You can also assign to the metadata as a whole or a specific field, as in the example above.

In addition the `content()` function provides the unprocessed text of the document.

The assumption is that is that you want to tokenize and hash the tokens from the start. If, however, you wish to do any of those steps yourself, you can load a document with `tokenizer = NULL`, then use `tokenize()` or `rehash()` to recompute the tokens and hashes.

Note that a `TextReuseTextDocument` can actually contain two kinds of hashes. The `hashes()` accessor gives you integer representations of each of the tokens in the document: if there are 100,000 tokens in the document, there will be 100,000 hashes. The `minhashes()` accessor gives you a signature that represents the document as a whole but not the specific tokens within it. See the minhash vignette for details: `vignette("textreuse-minhash")`.

### TextReuseCorpus

The class `TextReuseCorpus` provides a list of `TextReuseTextDocuments`. It derives from the S3 virtual class `Corpus` in the [tm](https://cran.r-project.org/package=tm) package. It can be created from a directory of files (or by providing a vector of paths to files).

```{r}
dir <- system.file("extdata/ats", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, tokenizer = tokenize_ngrams, n = 5,
                          progress = FALSE)
corpus
```

The names of the items in a `TextReuseCorpus` are the IDs of the documents. You can use these IDs to subset the corpus or to retrieve specific documents.

```{r}
names(corpus)
corpus[["remember00palm"]]
corpus[c("calltounconv00baxt", "lifeofrevrichard00baxt")]
```

Accessor functions such as `meta()`, `tokens()`, `hashes()`, and `wordcount()` have methods that work on corpora.

```{r}
wordcount(corpus)
```

Note that when creating a corpus, very short or empty documents will be skipped with a warning. A document must have enough words to create at least two n-grams. For example, if five-grams are desired, then the document must have at least six words.

## Tokenizers

One of the steps that is performed when loading a `TextReuseTextDocument`, either individual or in a corpus, is tokenization. Tokenization breaks up a text into pieces, often overlapping. These pieces are the features which are compared when measuring document similarity.

The textreuse package provides a number of tokenizers.

```{r}
text <- "How many roads must a man walk down\nBefore you'll call him a man?"

tokenize_words(text)
tokenize_sentences(text)
tokenize_ngrams(text, n = 3)
tokenize_skip_ngrams(text, n = 3, k = 2)
```

You can write your own tokenizers or use tokenizers from other packages. They should accept a character vector as their first argument.

As an example, we will write a tokenizer function using the \link[stringr]{stringr} package which splits a text on new lines, perhaps useful for poetry. Notice that the function takes a single string and returns a character vector with one element for each line. (A more robust tokenizer might strip blank lines and punctuation, include an option for lowercasing the text, and check for the validity of arguments.)

```{r}
poem <- "Roses are red\nViolets are blue\nI like using R\nAnd you should too"
cat(poem)

tokenize_lines <- function(string) {
  stringr::str_split(string, "\n+")[[1]]
}

tokenize_lines(poem)
```

## Hash functions

This package provides one function to hash tokens to integers, `hash_string()`. 
```{r}
hash_string(tokenize_words(text))
```

You can write your own hash functions, or use those provided by the [digest](https://cran.r-project.org/package=digest) package.

## Comparison functions

This package provides a number of comparison functions for measuring similarity. These functions take either a set (in which each token is counted one time) or a bag (in which each token is counted as many times as it appears) and compares it to another set or bag.

```{r}
a <- tokenize_words(paste("How does it feel, how does it feel?",
                          "To be without a home",
                          "Like a complete unknown, like a rolling stone"))
b <- tokenize_words(paste("How does it feel, how does it feel?",
                          "To be on your own, with no direction home",
                          "A complete unknown, like a rolling stone"))

jaccard_similarity(a, b)
jaccard_dissimilarity(a, b)
jaccard_bag_similarity(a, b)
ratio_of_matches(a, b)
```

See the documentation for `?similarity-functions` for details on what is measured with these functions.

You can write your own similarity functions, which should accept two sets or bags, `a` and `b`, should work on both character and numeric vectors, since they are used with either tokens or hashes of tokens, and should return a single numeric score for the comparison. You will need to implement a method for the `TextReuseTextDocument` class.

## Parallelization

This package will use multiple cores for a few functions is an option is set. This only benefits the corpus loading and tokenizing functions, which are often the slowest parts of an analysis. This is implemented with the [parallel package](https://cran.r-project.org/view=HighPerformanceComputing), and does not work on Windows machines. (Regardless of the options set, this package will never attempt to parallelize computations on Windows.)

To use the parallel option, you must specify the number of CPU cores that you wish to use:

```{R eval = FALSE}
options("mc.cores" = 4L)
```

If that option is set, this package will use multiple cores when possible.

You can figure out how many cores your computer has with `parallel::detectCores()`. See `help(package = "parallel")` for more details.

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsh_probability.R
\name{lsh_probability}
\alias{lsh_probability}
\alias{lsh_threshold}
\title{Probability that a candidate pair will be detected with LSH}
\usage{
lsh_probability(h, b, s)

lsh_threshold(h, b)
}
\arguments{
\item{h}{The number of minhash signatures.}

\item{b}{The number of LSH bands.}

\item{s}{The Jaccard similarity.}
}
\description{
Functions to help choose the correct parameters for the \code{\link{lsh}} and
\code{\link{minhash_generator}} functions. Use \code{lsh_threshold} to
determine the minimum Jaccard similarity for two documents for them to likely
be considered a match. Use \code{lsh_probability} to determine the
probability that a pair of documents with a known Jaccard similarity will be
detected.
}
\details{
Locality sensitive hashing returns a list of possible matches for
similar documents. How likely is it that a pair of documents will be detected
as a possible match? If \code{h} is the number of minhash signatures,
\code{b} is the number of bands in the LSH function (implying then that the
number of rows \code{r = h / b}), and \code{s} is the actual Jaccard
similarity of the two documents, then the probability \code{p} that the two
documents will be marked as a candidate pair is given by this equation.

\deqn{p = 1 - (1 - s^{r})^{b}}

According to \href{http://infolab.stanford.edu/~ullman/mmds/book.pdf}{MMDS},
that equation approximates an S-curve. This implies that there is a threshold
(\code{t}) for \code{s} approximated by this equation.

\deqn{t = \frac{1}{b}^{\frac{1}{r}}}
}
\examples{
# Threshold for default values
lsh_threshold(h = 200, b = 40)

# Probability for varying values of s
lsh_probability(h = 200, b = 40, s = .25)
lsh_probability(h = 200, b = 40, s = .50)
lsh_probability(h = 200, b = 40, s = .75)
}
\references{
Jure Leskovec, Anand Rajaraman, and Jeff Ullman,
 \href{http://www.mmds.org/#book}{\emph{Mining of Massive Datasets}}
 (Cambridge University Press, 2011), ch. 3.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsh_compare.R
\name{lsh_compare}
\alias{lsh_compare}
\title{Compare candidates identified by LSH}
\usage{
lsh_compare(candidates, corpus, f, progress = interactive())
}
\arguments{
\item{candidates}{A data frame returned by \code{\link{lsh_candidates}}.}

\item{corpus}{The same \code{\link{TextReuseCorpus}} corpus which was used to generate the candidates.}

\item{f}{A comparison function such as \code{\link{jaccard_similarity}}.}

\item{progress}{Display a progress bar while comparing documents.}
}
\value{
A data frame with values calculated for \code{score}.
}
\description{
The \code{\link{lsh_candidates}} only identifies potential matches, but
cannot estimate the actual similarity of the documents. This function takes a
data frame returned by \code{\link{lsh_candidates}} and applies a comparison
function to each of the documents in a corpus, thereby calculating the
document similarity score. Note that since your corpus will have minhash
signatures rather than hashes for the tokens itself, you will probably wish
to use \code{\link{tokenize}} to calculate new hashes. This can be done for
just the potentially similar documents. See the package vignettes for
details.
}
\examples{
dir <- system.file("extdata/legal", package = "textreuse")
minhash <- minhash_generator(200, seed = 234)
corpus <- TextReuseCorpus(dir = dir,
                          tokenizer = tokenize_ngrams, n = 5,
                          minhash_func = minhash)
buckets <- lsh(corpus, bands = 50)
candidates <- lsh_candidates(buckets)
lsh_compare(candidates, corpus, jaccard_similarity)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TextReuseTextDocument.R
\name{TextReuseTextDocument}
\alias{TextReuseTextDocument}
\alias{is.TextReuseTextDocument}
\alias{has_content}
\alias{has_tokens}
\alias{has_hashes}
\alias{has_minhashes}
\title{TextReuseTextDocument}
\usage{
TextReuseTextDocument(
  text,
  file = NULL,
  meta = list(),
  tokenizer = tokenize_ngrams,
  ...,
  hash_func = hash_string,
  minhash_func = NULL,
  keep_tokens = FALSE,
  keep_text = TRUE,
  skip_short = TRUE
)

is.TextReuseTextDocument(x)

has_content(x)

has_tokens(x)

has_hashes(x)

has_minhashes(x)
}
\arguments{
\item{text}{A character vector containing the text of the document. This
argument can be skipped if supplying \code{file}.}

\item{file}{The path to a text file, if \code{text} is not provided.}

\item{meta}{A list with named elements for the metadata associated with this
document. If a document is created using the \code{text} parameter, then
you must provide an \code{id} field, e.g., \code{meta = list(id =
"my_id")}. If the document is created using \code{file}, then the ID will
be created from the file name.}

\item{tokenizer}{A function to split the text into tokens. See
\code{\link{tokenizers}}. If value is \code{NULL}, then tokenizing and
hashing will be skipped.}

\item{...}{Arguments passed on to the \code{tokenizer}.}

\item{hash_func}{A function to hash the tokens. See
\code{\link{hash_string}}.}

\item{minhash_func}{A function to create minhash signatures of the document.
See \code{\link{minhash_generator}}.}

\item{keep_tokens}{Should the tokens be saved in the document that is
returned or discarded?}

\item{keep_text}{Should the text be saved in the document that is returned or
discarded?}

\item{skip_short}{Should short documents be skipped? (See details.)}

\item{x}{An R object to check.}
}
\value{
An object of class \code{TextReuseTextDocument}. This object inherits
  from the virtual S3 class \code{\link[NLP]{TextDocument}} in the NLP
  package. It contains the following elements: \describe{ \item{content}{The
  text of the document.} \item{tokens}{The tokens created from the text.}
  \item{hashes}{Hashes created from the tokens.} \item{minhashes}{The minhash
  signature of the document.} \item{metadata}{The document metadata,
  including the filename (if any) in \code{file}.} }
}
\description{
This is the constructor function for \code{TextReuseTextDocument} objects.
This class is used for comparing documents.
}
\details{
This constructor function follows a three-step process. It reads in
  the text, either from a file or from memory. It then tokenizes that text.
  Then it hashes the tokens. Most of the comparison functions in this package
  rely only on the hashes to make the comparison. By passing \code{FALSE} to
  \code{keep_tokens} and \code{keep_text}, you can avoid saving those
  objects, which can result in significant memory savings for large corpora.

  If \code{skip_short = TRUE}, this function will return \code{NULL} for very
  short or empty documents. A very short document is one where there are two
  few words to create at least two n-grams. For example, if five-grams are
  desired, then a document must be at least six words long. If no value of
  \code{n} is provided, then the function assumes a value of \code{n = 3}. A
  warning will be printed with the document ID of a skipped document.
}
\examples{
file <- system.file("extdata/legal/ny1850-match.txt", package = "textreuse")
doc  <- TextReuseTextDocument(file = file, meta = list(id = "ny1850"))
print(doc)
meta(doc)
head(tokens(doc))
head(hashes(doc))
\dontrun{
content(doc)
}
}
\seealso{
\link[=TextReuseTextDocument-accessors]{Accessors for TextReuse
  objects}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minhash.R
\name{minhash_generator}
\alias{minhash_generator}
\title{Generate a minhash function}
\usage{
minhash_generator(n = 200, seed = NULL)
}
\arguments{
\item{n}{The number of minhashes that the returned function should generate.}

\item{seed}{An option parameter to set the seed used in generating the random
numbers to ensure that the same minhash function is used on repeated
applications.}
}
\value{
A function which will take a character vector and return \code{n}
  minhashes.
}
\description{
A minhash value is calculated by hashing the strings in a character vector to
integers and then selecting the minimum value. Repeated minhash values are
generated by using different hash functions: these different hash functions
are created by using performing a bitwise \code{XOR} operation
(\code{\link{bitwXor}}) with a vector of random integers. Since it is vital
that the same random integers be used for each document, this function
generates another function which will always use the same integers. The
returned function is intended to be passed to the \code{hash_func} parameter
of \code{\link{TextReuseTextDocument}}.
}
\examples{
set.seed(253)
minhash <- minhash_generator(10)

# Example with a TextReuseTextDocument
file <- system.file("extdata/legal/ny1850-match.txt", package = "textreuse")
doc <- TextReuseTextDocument(file = file, hash_func = minhash,
                             keep_tokens = TRUE)
hashes(doc)

# Example with a character vector
is.character(tokens(doc))
minhash(tokens(doc))
}
\references{
Jure Leskovec, Anand Rajaraman, and Jeff Ullman,
  \href{http://www.mmds.org/#book}{\emph{Mining of Massive Datasets}}
  (Cambridge University Press, 2011), ch. 3. See also Matthew Casperson,
  "\href{http://matthewcasperson.blogspot.com/2013/11/minhash-for-dummies.html}{Minhash
   for Dummies}" (November 14, 2013).
}
\seealso{
\code{\link{lsh}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wordcount.R
\name{wordcount}
\alias{wordcount}
\title{Count words}
\usage{
wordcount(x)
}
\arguments{
\item{x}{The object containing a text.}
}
\value{
An integer vector for the word count.
}
\description{
This function counts words in a text, for example, a character vector, a
\code{\link{TextReuseTextDocument}}, some other object that inherits from
\code{\link[NLP]{TextDocument}}, or a all the documents in a
\code{\link{TextReuseCorpus}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rehash.R
\name{rehash}
\alias{rehash}
\title{Recompute the hashes for a document or corpus}
\usage{
rehash(x, func, type = c("hashes", "minhashes"))
}
\arguments{
\item{x}{A \code{\link{TextReuseTextDocument}} or
\code{\link{TextReuseCorpus}}.}

\item{func}{A function to either hash the tokens or to generate the minhash
signature. See \code{\link{hash_string}}, \code{\link{minhash_generator}}.}

\item{type}{Recompute the \code{hashes} or \code{minhashes}?}
}
\value{
The modified \code{\link{TextReuseTextDocument}} or
  \code{\link{TextReuseCorpus}}.
}
\description{
Given a \code{\link{TextReuseTextDocument}} or a
\code{\link{TextReuseCorpus}}, this function recomputes either the hashes or
the minhashes with the function specified. This implies that you have
retained the tokens with the \code{keep_tokens = TRUE} parameter.
}
\examples{
dir <- system.file("extdata/legal", package = "textreuse")
minhash1 <- minhash_generator(seed = 1)
corpus <- TextReuseCorpus(dir = dir, minhash_func = minhash1, keep_tokens = TRUE)
head(minhashes(corpus[[1]]))
minhash2 <- minhash_generator(seed = 2)
corpus <- rehash(corpus, minhash2, type = "minhashes")
head(minhashes(corpus[[2]]))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/similarity.R
\name{similarity-functions}
\alias{similarity-functions}
\alias{jaccard_similarity}
\alias{jaccard_dissimilarity}
\alias{jaccard_bag_similarity}
\alias{ratio_of_matches}
\title{Measure similarity/dissimilarity in documents}
\usage{
jaccard_similarity(a, b)

jaccard_dissimilarity(a, b)

jaccard_bag_similarity(a, b)

ratio_of_matches(a, b)
}
\arguments{
\item{a}{The first set (or bag) to be compared. The origin bag for
directional comparisons.}

\item{b}{The second set (or bag) to be compared. The destination bag for
directional comparisons.}
}
\description{
A set of functions which take two sets or bag of words and measure their
similarity or dissimilarity.
}
\details{
The functions \code{jaccard_similarity} and
  \code{jaccard_dissimilarity} provide the Jaccard measures of similarity or
  dissimilarity for two sets. The coefficients will be numbers between
  \code{0} and \code{1}. For the similarity coefficient, the higher the
  number the more similar the two sets are. When applied to two documents of
  class \code{\link{TextReuseTextDocument}}, the hashes in those documents
  are compared. But this function can be passed objects of any class accepted
  by the set functions in base R. So it is possible, for instance, to pass
  this function two character vectors comprised of word, line, sentence, or
  paragraph tokens, or those character vectors hashed as integers.

  The Jaccard similarity coeffecient is defined as follows:

  \deqn{J(A, B) = \frac{ | A \cap B | }{ | A \cup B | }}{ length(intersect(a,
  b)) / length(union(a, b))}

  The Jaccard dissimilarity is simply

  \deqn{1 - J(A, B)}

  The function \code{jaccard_bag_similarity} treats \code{a} and \code{b} as
  bags rather than sets, so that the result is a fraction where the numerator
  is the sum of each matching element counted the minimum number of times it
  appears in each bag, and the denominator is the sum of the lengths of both
  bags. The maximum value for the Jaccard bag similarity is \code{0.5}.

  The function \code{ratio_of_matches} finds the ratio between the number of
  items in \code{b} that are also in \code{a} and the total number of items
  in \code{b}. Note that this similarity measure is directional: it measures
  how much \code{b} borrows from \code{a}, but says nothing about how much of
  \code{a} borrows from \code{b}.
}
\examples{
jaccard_similarity(1:6, 3:10)
jaccard_dissimilarity(1:6, 3:10)

a <- c("a", "a", "a", "b")
b <- c("a", "a", "b", "b", "c")
jaccard_similarity(a, b)
jaccard_bag_similarity(a, b)
ratio_of_matches(a, b)
ratio_of_matches(b, a)

ny         <- system.file("extdata/legal/ny1850-match.txt", package = "textreuse")
ca_match   <- system.file("extdata/legal/ca1851-match.txt", package = "textreuse")
ca_nomatch <- system.file("extdata/legal/ca1851-nomatch.txt", package = "textreuse")

ny         <- TextReuseTextDocument(file = ny,
                                    meta = list(id = "ny"))
ca_match   <- TextReuseTextDocument(file = ca_match,
                                    meta = list(id = "ca_match"))
ca_nomatch <- TextReuseTextDocument(file = ca_nomatch,
                                    meta = list(id = "ca_nomatch"))

# These two should have higher similarity scores
jaccard_similarity(ny, ca_match)
ratio_of_matches(ny, ca_match)

# These two should have lower similarity scores
jaccard_similarity(ny, ca_nomatch)
ratio_of_matches(ny, ca_nomatch)

}
\references{
Jure Leskovec, Anand Rajaraman, and Jeff Ullman,
  \href{http://www.mmds.org/#book}{\emph{Mining of Massive Datasets}}
  (Cambridge University Press, 2011).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TextReuseCorpus.R
\name{TextReuseCorpus}
\alias{TextReuseCorpus}
\alias{is.TextReuseCorpus}
\alias{skipped}
\title{TextReuseCorpus}
\usage{
TextReuseCorpus(
  paths,
  dir = NULL,
  text = NULL,
  meta = list(),
  progress = interactive(),
  tokenizer = tokenize_ngrams,
  ...,
  hash_func = hash_string,
  minhash_func = NULL,
  keep_tokens = FALSE,
  keep_text = TRUE,
  skip_short = TRUE
)

is.TextReuseCorpus(x)

skipped(x)
}
\arguments{
\item{paths}{A character vector of paths to files to be opened.}

\item{dir}{The path to a directory of text files.}

\item{text}{A character vector (possibly named) of documents.}

\item{meta}{A list with named elements for the metadata associated with this
corpus.}

\item{progress}{Display a progress bar while loading files.}

\item{tokenizer}{A function to split the text into tokens. See
\code{\link{tokenizers}}. If value is \code{NULL}, then tokenizing and
hashing will be skipped.}

\item{...}{Arguments passed on to the \code{tokenizer}.}

\item{hash_func}{A function to hash the tokens. See
\code{\link{hash_string}}.}

\item{minhash_func}{A function to create minhash signatures of the document.
See \code{\link{minhash_generator}}.}

\item{keep_tokens}{Should the tokens be saved in the documents that are
returned or discarded?}

\item{keep_text}{Should the text be saved in the documents that are returned
or discarded?}

\item{skip_short}{Should short documents be skipped? (See details.)}

\item{x}{An R object to check.}
}
\description{
This is the constructor function for a \code{TextReuseCorpus}, modeled on the
virtual S3 class \code{\link[tm]{Corpus}} from the \code{tm} package. The
object is a \code{TextReuseCorpus}, which is basically a list containing
objects of class \code{\link{TextReuseTextDocument}}. Arguments are passed
along to that constructor function. To create the corpus, you can pass either
a character vector of paths to text files using the \code{paths =} parameter,
a directory containing text files (with any extension) using the \code{dir =}
parameter, or a character vector of documents using the \code{text = }
parameter, where each element in the characer vector is a document. If the
character vector passed to \code{text = } has names, then those names will be
used as the document IDs. Otherwise, IDs will be assigned to the documents.
Only one of the \code{paths}, \code{dir}, or \code{text} parameters should be
specified.
}
\details{
If \code{skip_short = TRUE}, this function will skip very short or
  empty documents. A very short document is one where there are two few words
  to create at least two n-grams. For example, if five-grams are desired,
  then a document must be at least six words long. If no value of \code{n} is
  provided, then the function assumes a value of \code{n = 3}. A warning will
  be printed with the document ID of each skipped document. Use
  \code{skipped()} to get the IDs of skipped documents.

  This function will use multiple cores on non-Windows machines if the
  \code{"mc.cores"} option is set. For example, to use four cores:
  \code{options("mc.cores" = 4L)}.
}
\examples{
dir <- system.file("extdata/legal", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, meta = list("description" = "Field Codes"))
# Subset by position or file name
corpus[[1]]
names(corpus)
corpus[["ca1851-match"]]

}
\seealso{
\link[=TextReuseTextDocument-accessors]{Accessors for TextReuse
  objects}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{hash_string}
\alias{hash_string}
\title{Hash a string to an integer}
\usage{
hash_string(x)
}
\arguments{
\item{x}{A character vector to be hashed.}
}
\value{
A vector of integer hashes.
}
\description{
Hash a string to an integer
}
\examples{
s <- c("How", "many", "roads", "must", "a", "man", "walk", "down")
hash_string(s)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsh_candidates.R
\name{lsh_candidates}
\alias{lsh_candidates}
\title{Candidate pairs from LSH comparisons}
\usage{
lsh_candidates(buckets)
}
\arguments{
\item{buckets}{A data frame returned from \code{\link{lsh}}.}
}
\value{
A data frame of candidate pairs.
}
\description{
Given a data frame of LSH buckets returned from \code{\link{lsh}}, this
function returns the potential candidates.
}
\examples{
dir <- system.file("extdata/legal", package = "textreuse")
minhash <- minhash_generator(200, seed = 234)
corpus <- TextReuseCorpus(dir = dir,
                          tokenizer = tokenize_ngrams, n = 5,
                          minhash_func = minhash)
buckets <- lsh(corpus, bands = 50)
lsh_candidates(buckets)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filenames.R
\name{filenames}
\alias{filenames}
\title{Filenames from paths}
\usage{
filenames(paths, extension = FALSE)
}
\arguments{
\item{paths}{A character vector of paths.}

\item{extension}{Should the file extension be preserved?}
}
\description{
This function takes a character vector of paths and returns just the file
name, by default without the extension. A \code{\link{TextReuseCorpus}} uses
the paths to the files in the corpus as the names of the list. This function
is intended to turn those paths into more manageable identifiers.
}
\examples{
paths <- c("corpus/one.txt", "corpus/two.md", "corpus/three.text")
filenames(paths)
filenames(paths, extension = TRUE)
}
\seealso{
\code{\link{basename}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tokenizers.R
\name{tokenizers}
\alias{tokenizers}
\alias{tokenize_words}
\alias{tokenize_sentences}
\alias{tokenize_ngrams}
\alias{tokenize_skip_ngrams}
\title{Split texts into tokens}
\usage{
tokenize_words(string, lowercase = TRUE)

tokenize_sentences(string, lowercase = TRUE)

tokenize_ngrams(string, lowercase = TRUE, n = 3)

tokenize_skip_ngrams(string, lowercase = TRUE, n = 3, k = 1)
}
\arguments{
\item{string}{A character vector of length 1 to be tokenized.}

\item{lowercase}{Should the tokens be made lower case?}

\item{n}{For n-gram tokenizers, the number of words in each n-gram.}

\item{k}{For the skip n-gram tokenizer, the maximum skip distance between
words. The function will compute all skip n-grams between \code{0} and
\code{k}.}
}
\value{
A character vector containing the tokens.
}
\description{
These functions each turn a text into tokens. The \code{tokenize_ngrams}
functions returns shingled n-grams.
}
\details{
These functions will strip all punctuation.
}
\examples{
dylan <- "How many roads must a man walk down? The answer is blowin' in the wind."
tokenize_words(dylan)
tokenize_sentences(dylan)
tokenize_ngrams(dylan, n = 2)
tokenize_skip_ngrams(dylan, n = 3, k = 2)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/textreuse-package.r
\docType{package}
\name{textreuse-package}
\alias{textreuse}
\alias{textreuse-package}
\title{textreuse: Detect Text Reuse and Document Similarity}
\description{
Tools for measuring similarity among documents and detecting
    passages which have been reused. Implements shingled n-gram, skip n-gram,
    and other tokenizers; similarity/dissimilarity functions; pairwise
    comparisons; minhash and locality sensitive hashing algorithms; and a
    version of the Smith-Waterman local alignment algorithm suitable for
    natural language.
}
\details{
The best place to begin with this package in the introductory vignette.

\code{vignette("textreuse-introduction", package = "textreuse")}

After reading that vignette, the "pairwise" and "minhash" vignettes introduce
specific paths for working with the package.

\code{vignette("textreuse-pairwise", package = "textreuse")}

\code{vignette("textreuse-minhash", package = "textreuse")}

\code{vignette("textreuse-alignment", package = "textreuse")}

Another good place to begin with the package is the documentation for loading
documents (\code{\link{TextReuseTextDocument}} and
\code{\link{TextReuseCorpus}}), for \link{tokenizers},
\link[=similarity-functions]{similarity functions}, and
\link[=lsh]{locality-sensitive hashing}.
}
\references{
The sample data provided in the \code{extdata/legal} directory is
  taken from a
  \href{http://lincolnmullen.com/blog/corpus-of-american-tract-society-publications/}{corpus
   of American Tract Society publications} from the nineteen-century,
  gathered from the \href{https://archive.org/}{Internet Archive}.

  The sample data provided in the \code{extdata/legal} directory, are taken
  from the following nineteenth-century codes of civil procedure from
  California and New York.

  \emph{Final Report of the Commissioners on Practice and Pleadings}, in 2
  \emph{Documents of the Assembly of New York}, 73rd Sess., No. 16, (1850):
  243-250, sections 597-613.
  \href{http://books.google.com/books?id=9HEbAQAAIAAJ&pg=PA243#v=onepage&q&f=false}{Google
   Books}.

  \emph{An Act To Regulate Proceedings in Civil Cases}, 1851 \emph{California
  Laws} 51, 51-53 sections 4-17; 101, sections 313-316.
  \href{http://books.google.com/books?id=4PHEAAAAIAAJ&pg=PA51#v=onepage&q&f=false}{Google
   Books}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/textreuse}
  \item \url{https://github.com/ropensci/textreuse}
  \item Report bugs at \url{https://github.com/ropensci/textreuse/issues}
}

}
\author{
\strong{Maintainer}: Lincoln Mullen \email{lincoln@lincolnmullen.com} (\href{https://orcid.org/0000-0001-5103-6917}{ORCID})

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsh.R
\name{lsh}
\alias{lsh}
\title{Locality sensitive hashing for minhash}
\usage{
lsh(x, bands, progress = interactive())
}
\arguments{
\item{x}{A \code{\link{TextReuseCorpus}} or
\code{\link{TextReuseTextDocument}}.}

\item{bands}{The number of bands to use for locality sensitive hashing. The
number of hashes in the documents in the corpus must be evenly divisible by
the number of bands. See \code{\link{lsh_threshold}} and
\code{\link{lsh_probability}} for guidance in selecting the number of bands
and hashes.}

\item{progress}{Display a progress bar while comparing documents.}
}
\value{
A data frame (with the additional class \code{lsh_buckets}),
 containing a column with the document IDs and a column with their LSH
 signatures, or buckets.
}
\description{
Locality sensitive hashing (LSH) discovers potential matches among a corpus of
documents quickly, so that only likely pairs can be compared.
}
\details{
Locality sensitive hashing is a technique for detecting document
 similarity that does not require pairwise comparisons. When comparing pairs
 of documents, the number of pairs grows rapidly, so that only the smallest
 corpora can be compared pairwise in a reasonable amount of computation time.
 Locality sensitive hashing, on the other hand, takes a document which has
 been tokenized and hashed using a minhash algorithm. (See
 \code{\link{minhash_generator}}.) Each set of minhash signatures is then
 broken into bands comprised of a certain number of rows. (For example, 200
 minhash signatures might be broken down into 20 bands each containing 10
 rows.) Each band is then hashed to a bucket. Documents with identical rows
 in a band will be hashed to the same bucket. The likelihood that a document
 will be marked as a potential duplicate is proportional to the number of
 bands and inversely proportional to the number of rows in each band.

 This function returns a data frame with the additional class
 \code{lsh_buckets}. The LSH technique only requires that the signatures for
 each document be calculated once. So it is possible, as long as one uses the
 same minhash function and the same number of bands, to combine the outputs
 from this function at different times. The output can thus be treated as a
 kind of cache of LSH signatures.

 To extract pairs of documents from the output of this function, see
 \code{\link{lsh_candidates}}.
}
\examples{
dir <- system.file("extdata/legal", package = "textreuse")
minhash <- minhash_generator(200, seed = 235)
corpus <- TextReuseCorpus(dir = dir,
                          tokenizer = tokenize_ngrams, n = 5,
                          minhash_func = minhash)
buckets <- lsh(corpus, bands = 50)
buckets
}
\references{
Jure Leskovec, Anand Rajaraman, and Jeff Ullman,
 \href{http://www.mmds.org/#book}{\emph{Mining of Massive Datasets}}
 (Cambridge University Press, 2011), ch. 3. See also Matthew Casperson,
 "\href{http://matthewcasperson.blogspot.com/2013/11/minhash-for-dummies.html}{Minhash
  for Dummies}" (November 14, 2013).
}
\seealso{
\code{\link{minhash_generator}}, \code{\link{lsh_candidates}},
 \code{\link{lsh_query}}, \code{\link{lsh_probability}},
 \code{\link{lsh_threshold}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conversion-functions.R
\name{as.matrix.textreuse_candidates}
\alias{as.matrix.textreuse_candidates}
\title{Convert candidates data frames to other formats}
\usage{
\method{as.matrix}{textreuse_candidates}(x, ...)
}
\arguments{
\item{x}{An object of class \code{\link[=lsh_compare]{textreuse_candidates}}.}

\item{...}{Additional arguments.}
}
\value{
A similarity matrix with row and column names containing document IDs.
}
\description{
These S3 methods convert a \code{textreuse_candidates} object to a matrix.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align_local.R
\name{align_local}
\alias{align_local}
\title{Local alignment of natural language texts}
\usage{
align_local(
  a,
  b,
  match = 2L,
  mismatch = -1L,
  gap = -1L,
  edit_mark = "#",
  progress = interactive()
)
}
\arguments{
\item{a}{A character vector of length one, or a
\code{\link{TextReuseTextDocument}}.}

\item{b}{A character vector of length one, or a
\code{\link{TextReuseTextDocument}}.}

\item{match}{The score to assign a matching word. Should be a positive
integer.}

\item{mismatch}{The score to assign a mismatching word. Should be a negative
integer or zero.}

\item{gap}{The penalty for opening a gap in the sequence. Should be a
negative integer or zero.}

\item{edit_mark}{A single character used for displaying for displaying
insertions/deletions in the documents.}

\item{progress}{Display a progress bar and messages while computing the
alignment.}
}
\value{
A list with the class \code{textreuse_alignment}. This list contains
  several elements: \itemize{ \item \code{a_edit} and \code{b_edit}:
  Character vectors of the sequences with edits marked. \item \code{score}:
  The score of the optimal alignment. }
}
\description{
This function takes two texts, either as strings or as
\code{TextReuseTextDocument} objects, and finds the optimal local alignment
of those texts. A local alignment finds the best matching subset of the two
documents. This function adapts the
\href{https://en.wikipedia.org/wiki/Smith-Waterman_algorithm}{Smith-Waterman
algorithm}, used for genetic sequencing, for use with natural language. It
compare the texts word by word (the comparison is case-insensitive) and
scores them according to a set of parameters. These parameters define the
score for a \code{match}, and the penalties for a \code{mismatch} and for
opening a \code{gap} (i.e., the first mismatch in a potential sequence). The
function then reports the optimal local alignment. Only the subset of the
documents that is a match is included. Insertions or deletions in the text
are reported with the \code{edit_mark} character.
}
\details{
The compute time of this function is proportional to the product of the
lengths of the two documents. Thus, longer documents will take considerably
more time to compute. This function has been tested with pairs of documents
containing about 25 thousand words each.

If the function reports that there were multiple optimal alignments, then it
is likely that there is no strong match in the document.

The score reported for the local alignment is dependent on both the size of
the documents and on the strength of the match, as well as on the parameters
for match, mismatch, and gap penalties, so the scores are not directly
comparable.
}
\examples{
align_local("The answer is blowin' in the wind.",
            "As the Bob Dylan song says, the answer is blowing in the wind.")

# Example of matching documents from a corpus
dir <- system.file("extdata/legal", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, progress = FALSE)
alignment <- align_local(corpus[["ca1851-match"]], corpus[["ny1850-match"]])
str(alignment)

}
\references{
For a useful description of the algorithm, see
  \href{http://etherealbits.com/2013/04/string-alignment-dynamic-programming-dna/}{this
  post}. For the application of the Smith-Waterman algorithm to natural
  language, see David A. Smith, Ryan Cordell, and Elizabeth Maddock Dillon,
  "Infectious Texts: Modeling Text Reuse in Nineteenth-Century Newspapers."
  IEEE International Conference on Big Data, 2013,
  \url{http://hdl.handle.net/2047/d20004858}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsh_query.R
\name{lsh_query}
\alias{lsh_query}
\title{Query a LSH cache for matches to a single document}
\usage{
lsh_query(buckets, id)
}
\arguments{
\item{buckets}{An \code{lsh_buckets} object created by \code{\link{lsh}}.}

\item{id}{The document ID to find matches for.}
}
\value{
An \code{lsh_candidates} data frame with matches to the document specified.
}
\description{
This function retrieves the matches for a single document from an \code{lsh_buckets} object created by \code{\link{lsh}}. See \code{\link{lsh_candidates}} to retrieve all pairs of matches.
}
\examples{
dir <- system.file("extdata/legal", package = "textreuse")
minhash <- minhash_generator(200, seed = 235)
corpus <- TextReuseCorpus(dir = dir,
                          tokenizer = tokenize_ngrams, n = 5,
                          minhash_func = minhash)
buckets <- lsh(corpus, bands = 50)
lsh_query(buckets, "ny1850-match")

}
\seealso{
\code{\link{lsh}}, \code{\link{lsh_candidates}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TextReuseTextDocument.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{meta}
\alias{meta<-}
\alias{content}
\alias{content<-}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{NLP}{\code{\link[NLP]{content}}, \code{\link[NLP]{content<-}}, \code{\link[NLP]{meta}}, \code{\link[NLP]{meta<-}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tokenize.R
\name{tokenize}
\alias{tokenize}
\title{Recompute the tokens for a document or corpus}
\usage{
tokenize(
  x,
  tokenizer,
  ...,
  hash_func = hash_string,
  minhash_func = NULL,
  keep_tokens = FALSE,
  keep_text = TRUE
)
}
\arguments{
\item{x}{A \code{\link{TextReuseTextDocument}} or
\code{\link{TextReuseCorpus}}.}

\item{tokenizer}{A function to split the text into tokens. See
\code{\link{tokenizers}}.}

\item{...}{Arguments passed on to the \code{tokenizer}.}

\item{hash_func}{A function to hash the tokens. See
\code{\link{hash_string}}.}

\item{minhash_func}{A function to create minhash signatures. See
\code{\link{minhash_generator}}.}

\item{keep_tokens}{Should the tokens be saved in the document that is
returned or discarded?}

\item{keep_text}{Should the text be saved in the document that is returned or
discarded?}
}
\value{
The modified \code{\link{TextReuseTextDocument}} or
  \code{\link{TextReuseCorpus}}.
}
\description{
Given a \code{\link{TextReuseTextDocument}} or a
\code{\link{TextReuseCorpus}}, this function recomputes the tokens and hashes
with the functions specified. Optionally, it can also recompute the minhash signatures.
}
\examples{
dir <- system.file("extdata/legal", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir, tokenizer = NULL)
corpus <- tokenize(corpus, tokenize_ngrams)
head(tokens(corpus[[1]]))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairwise_candidates.R
\name{pairwise_candidates}
\alias{pairwise_candidates}
\title{Candidate pairs from pairwise comparisons}
\usage{
pairwise_candidates(m, directional = FALSE)
}
\arguments{
\item{m}{A matrix from \code{\link{pairwise_compare}}.}

\item{directional}{Should be set to the same value as in
\code{\link{pairwise_compare}}.}
}
\value{
A data frame containing all the non-\code{NA} values from \code{m}.
  Columns \code{a} and \code{b} are the IDs from the original corpus as
  passed to the comparison function. Column \code{score} is the score
  returned by the comparison function.
}
\description{
Converts a comparison matrix generated by \code{\link{pairwise_compare}} into a
data frame of candidates for matches.
}
\examples{
dir <- system.file("extdata/legal", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir)

m1 <- pairwise_compare(corpus, ratio_of_matches, directional = TRUE)
pairwise_candidates(m1, directional = TRUE)

m2 <- pairwise_compare(corpus, jaccard_similarity)
pairwise_candidates(m2)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsh_subset.R
\name{lsh_subset}
\alias{lsh_subset}
\title{List of all candidates in a corpus}
\usage{
lsh_subset(candidates)
}
\arguments{
\item{candidates}{A data frame of candidate pairs from
\code{\link{lsh_candidates}}.}
}
\value{
A character vector of document IDs from the candidate pairs, to be
  used to subset the \code{\link{TextReuseCorpus}}.
}
\description{
List of all candidates in a corpus
}
\examples{
dir <- system.file("extdata/legal", package = "textreuse")
minhash <- minhash_generator(200, seed = 234)
corpus <- TextReuseCorpus(dir = dir,
                          tokenizer = tokenize_ngrams, n = 5,
                          minhash_func = minhash)
buckets <- lsh(corpus, bands = 50)
candidates <- lsh_candidates(buckets)
lsh_subset(candidates)
corpus[lsh_subset(candidates)]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairwise_compare.R
\name{pairwise_compare}
\alias{pairwise_compare}
\title{Pairwise comparisons among documents in a corpus}
\usage{
pairwise_compare(corpus, f, ..., directional = FALSE, progress = interactive())
}
\arguments{
\item{corpus}{A \code{\link{TextReuseCorpus}}.}

\item{f}{The function to apply to \code{x} and \code{y}.}

\item{...}{Additional arguments passed to \code{f}.}

\item{directional}{Some comparison functions are commutative, so that
\code{f(a, b) == f(b, a)} (e.g., \code{\link{jaccard_similarity}}). Other
functions are directional, so that \code{f(a, b)} measures \code{a}'s
borrowing from \code{b}, which may not be the same as \code{f(b, a)} (e.g.,
\code{\link{ratio_of_matches}}). If \code{directional} is \code{FALSE},
then only the minimum number of comparisons will be made, i.e., the upper
triangle of the matrix. If \code{directional} is \code{TRUE}, then both
directional comparisons will be measured. In no case, however, will
documents be compared to themselves, i.e., the diagonal of the matrix.}

\item{progress}{Display a progress bar while comparing documents.}
}
\value{
A square matrix with dimensions equal to the length of the corpus,
  and row and column names set by the names of the documents in the corpus. A
  value of \code{NA} in the matrix indicates that a comparison was not made.
  In cases of directional comparisons, then the comparison reported is
  \code{f(row, column)}.
}
\description{
Given a \code{\link{TextReuseCorpus}} containing documents of class
\code{\link{TextReuseTextDocument}}, this function applies a comparison
function to every pairing of documents, and returns a matrix with the
comparison scores.
}
\examples{
dir <- system.file("extdata/legal", package = "textreuse")
corpus <- TextReuseCorpus(dir = dir)
names(corpus) <- filenames(names(corpus))

# A non-directional comparison
pairwise_compare(corpus, jaccard_similarity)

# A directional comparison
pairwise_compare(corpus, ratio_of_matches, directional = TRUE)
}
\seealso{
See these document comparison functions,
  \code{\link{jaccard_similarity}}, \code{\link{ratio_of_matches}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TextReuseTextDocument.R
\name{TextReuseTextDocument-accessors}
\alias{TextReuseTextDocument-accessors}
\alias{tokens}
\alias{tokens<-}
\alias{hashes}
\alias{hashes<-}
\alias{minhashes}
\alias{minhashes<-}
\title{Accessors for TextReuse objects}
\usage{
tokens(x)

tokens(x) <- value

hashes(x)

hashes(x) <- value

minhashes(x)

minhashes(x) <- value
}
\arguments{
\item{x}{The object to access.}

\item{value}{The value to assign.}
}
\value{
Either a vector or a named list of vectors.
}
\description{
Accessor functions to read and write components of
\code{\link{TextReuseTextDocument}} and \code{\link{TextReuseCorpus}}
objects.
}
