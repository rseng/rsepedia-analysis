# rOpenSci: The *hunspell* package

> High-Performance Stemmer, Tokenizer, and Spell Checker for R

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/ropensci/hunspell.svg?branch=master)](https://travis-ci.org/ropensci/hunspell)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/hunspell?branch=master&svg=true)](https://ci.appveyor.com/project/jeroen/hunspell)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/hunspell)](https://cran.r-project.org/package=hunspell)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/hunspell)](https://cran.r-project.org/package=hunspell)
[![Github Stars](https://img.shields.io/github/stars/ropensci/hunspell.svg?style=social&label=Github)](https://github.com/ropensci/hunspell)

Low level spell checker and morphological analyzer based on the 
famous hunspell library <https://hunspell.github.io>. The package can analyze
or check individual words as well as tokenize text, latex, html or xml documents.
For a more user-friendly interface use the 'spelling' package which builds on
this package with utilities to automate checking of files, documentation and 
vignettes in all common formats.

## Installation

This package includes a bundled version of libhunspell and no longer
depends on external system libraries:

```r
install.packages("hunspell")
```


## Documentation

About the R package:
 - Blog post: [Hunspell: Spell Checker and Text Parser for R](https://www.opencpu.org/posts/hunspell-release/)
 - Blog post: [Stemming and Spell Checking in R](https://www.opencpu.org/posts/hunspell-1-2/)

## Hello World

```r
# Check individual words
words <- c("beer", "wiskey", "wine")
correct <- hunspell_check(words)
print(correct)

# Find suggestions for incorrect words
hunspell_suggest(words[!correct])

# Extract incorrect from a piece of text
bad <- hunspell("spell checkers are not neccessairy for langauge ninja's")
print(bad[[1]])
hunspell_suggest(bad[[1]])

# Stemming
words <- c("love", "loving", "lovingly", "loved", "lover", "lovely", "love")
hunspell_stem(words)
hunspell_analyze(words)
```

The [spelling](https://docs.ropensci.org/spelling) package uses this package to spell R package documentation:

```r
# Spell check a package
library(spelling)
spell_check_package("~/mypackage")
```

[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "The hunspell package: High-Performance Stemmer, Tokenizer, and Spell Checker for R"
date: "`r Sys.Date()`"
output:
  html_document:
    fig_caption: false
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    toc_depth: 3
vignette: >
  %\VignetteIndexEntry{The hunspell package: High-Performance Stemmer, Tokenizer, and Spell Checker for R}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE, message = FALSE}
knitr::opts_chunk$set(comment = "")
NOT_CRAN = isTRUE(nchar(Sys.getenv('NOT_CRAN')) || (Sys.getenv('USER') == 'jeroen'))
```

Hunspell is the spell checker library used by LibreOffice, OpenOffice, Mozilla Firefox, Google Chrome, Mac OS-X, InDesign, Opera, RStudio and many others. It provides a system for tokenizing, stemming and spelling in almost any language or alphabet. The R package exposes both the high-level spell-checker as well as low-level stemmers and tokenizers which analyze or extract individual words from various formats (text, html, xml, latex).

Hunspell uses a special dictionary format that defines which characters, words and conjugations are valid in a given language. The examples below use the (default) `"en_US"` dictionary. However each function can be used in another language by setting a custom dictionary in the `dict` parameter. See the section on [dictionaries](#hunspell_dictionaries) below.

## Spell Checking

Spell checking text consists of the following steps: 
 
 1. Parse a document by extracting (**tokenizing**) words that we want to check
 2. Analyze each word by breaking it down in it's root (**stemming**) and conjugation affix
 3. Lookup in a **dictionary** if the word+affix combination if valid for your language
 4. (optional) For incorrect words, **suggest** corrections by finding similar (correct) words in the dictionary
 
We can do each these steps manually or have Hunspell do them automatically.

### Check Individual Words

The `hunspell_check` and `hunspell_suggest` functions can test individual words for correctness, and suggest similar (correct) words that look similar to the given (incorrect) word.

```{r}
library(hunspell)

# Check individual words
words <- c("beer", "wiskey", "wine")
correct <- hunspell_check(words)
print(correct)

# Find suggestions for incorrect words
hunspell_suggest(words[!correct])
```

### Check Documents

In practice we often want to spell check an entire document at once by searching for incorrect words. This is done using the `hunspell` function:


```{r}
bad <- hunspell("spell checkers are not neccessairy for langauge ninjas")
print(bad[[1]])
hunspell_suggest(bad[[1]])
```

Besides plain text, `hunspell` supports various document formats, such as html or latex:

```{r}
download.file("https://arxiv.org/e-print/1406.4806v1", "1406.4806v1.tar.gz",  mode = "wb")
untar("1406.4806v1.tar.gz", "content.tex")
text <- readLines("content.tex", warn = FALSE)
bad_words <- hunspell(text, format = "latex")
sort(unique(unlist(bad_words)))
```

### Check PDF files

Use the text-extraction from the `pdftools` package to spell check text from PDF files!

```{r, eval = require('pdftools')}
text <- pdftools::pdf_text('https://www.gnu.org/licenses/quick-guide-gplv3.pdf')
bad_words <- hunspell(text)
sort(unique(unlist(bad_words)))
```

### Check Manual Pages

The `spelling` package builds on hunspell and has a wrapper to spell-check manual pages from R packages. Results might contain a lot of false positives for technical jargon, but you might also catch a typo or two. Point it to the root of your source package:

```{r, eval=FALSE}
spelling::spell_check_package("~/workspace/V8")
```
```
  WORD          FOUND IN
ECMA          V8.Rd:16, description:2,4
ECMAScript    description:2
emscripten    description:5
htmlwidgets   JS.Rd:16
JSON          V8.Rd:33,38,39,57,58,59,120
jsonlite      V8.Rd:42
Ooms          V8.Rd:41,120
Xie           JS.Rd:26
Yihui         JS.Rd:26
```

## Morphological Analysis

In order to lookup a word in a dictionary, hunspell needs to break it down in a stem (**stemming**) and conjugation affix. The `hunspell` function does this automatically but we can also do it manually.

### Stemming Words

The `hunspell_stem` looks up words from the dictionary which match the root of the given word. Note that the function returns a list because some words can have multiple matches.

```{r}
# Stemming
words <- c("love", "loving", "lovingly", "loved", "lover", "lovely")
hunspell_stem(words)
```

### Analyzing Words

The `hunspell_analyze` function is similar, but it returns both the stem and the affix syntax of the word:

```{r}
hunspell_analyze(words)
```

## Tokenizing

To support spell checking on documents, Hunspell includes parsers for various document formats, including  *text*, *html*, *xml*, *man* or *latex*. The Hunspell package also exposes these tokenizers directly so they can be used for other application than spell checking.


```{r}
text <- readLines("content.tex", warn = FALSE)
allwords <- hunspell_parse(text, format = "latex")

# Third line (title) only
print(allwords[[3]])
```

### Summarizing Text

In text analysis we often want to summarize text via it's stems. For example we can count words for display in a wordcloud:

```{r}
allwords <- hunspell_parse(janeaustenr::prideprejudice)
stems <- unlist(hunspell_stem(unlist(allwords)))
words <- sort(table(stems), decreasing = TRUE)
print(head(words, 30))
```

Most of these are stop words. Let's filter these out:

```{r}
df <- as.data.frame(words)
df$stems <- as.character(df$stems)
stops <- df$stems %in% stopwords::stopwords(source="stopwords-iso")
wcdata <- head(df[!stops,], 150)
print(wcdata, max = 40)
```

```{r, eval = NOT_CRAN}
library(wordcloud2)
names(wcdata) <- c("word", "freq")
wcdata$freq <- (wcdata$freq)^(2/3)
wordcloud2(wcdata)
```



## Hunspell Dictionaries

Hunspell is based on MySpell and is backward-compatible with MySpell and aspell dictionaries. Chances are your dictionaries in your language are already available on your system! 

A Hunspell dictionary consists of two files:

 - The `[lang].aff` file specifies the affix syntax for the language
 - The `[lang].dic` file contains a wordlist formatted using syntax from the aff file.
 
Typically both files are located in the same directory and share the same filename, for example `en_GB.aff` and `en_GB.dic`. The `list_dictionaries()` function lists available dictionaries in the current directory and standard system paths where dictionaries are usually installed.

```{r}
list_dictionaries()
```

The `dictionary` function is then used to load any of these available dictionaries:

```{r}
dictionary("en_GB")
```


If the files are not in one of the standard paths you can also specify the full path to either or both the dic and aff file:

```{r, eval = FALSE}
dutch <- dictionary("~/workspace/Dictionaries/Dutch.dic")
print(dutch)
```
```
<hunspell dictionary>
 affix: /Users/jeroen/workspace/Dictionaries/Dutch.aff 
 dictionary: /Users/jeroen/workspace/Dictionaries/Dutch.dic 
 encoding: UTF-8 
 wordchars: '-./0123456789\ĳ’ 
```
### Setting a Language

The hunspell R package includes dictionaries for `en_US` and `en_GB`. So if you you don't speak `en_US` you can always switch to the British English:

```{r}
hunspell("My favourite colour to visualise is grey")
hunspell("My favourite colour to visualise is grey", dict = 'en_GB')
```

If you want to use another language you need to make sure that the dictionary is available from your system. The `dictionary` function is used to read in dictionary. 

```{r, eval = FALSE}
dutch <- dictionary("~/workspace/Dictionaries/Dutch.dic")
hunspell("Hij heeft de klok wel horen luiden, maar weet niet waar de klepel hangt", dict = dutch)
```

Note that if the `dict` argument is a string, it will be passed on to the `dictionary` function.

### Installing Dictionaries in RStudio

RStudio users can install various dictionaries via the "Global Options" menu of the IDE. Once these dictionaries are installed they become available to the `hunspell` and `spelling` package as well.

<img data-external="1" src="https://jeroen.github.io/images/rs-hunspell.png" alt="rstudio screenshot">

### Dictionaries on Linux

The best way to install dictionaries on __Linux__ is via the system package manager. For example on if you would like to install the Austrian-German dictionary on **Debian** or **Ubuntu** you either need the [`hunspell-de-at`](https://packages.debian.org/testing/hunspell-de-at) or [`myspell-de-at`](https://packages.debian.org/testing/myspell-de-at) package:

```sh
sudo apt-get install hunspell-de-at
```

On **Fedora** and **CentOS** / **RHEL** all German dialects are included with the [`hunspell-de`](https://src.fedoraproject.org/rpms/hunspell-de) package

```sh
sudo yum install hunspell-de
```

After installing this you should be able to load the dictionary:

```r
dict <- dictionary('de_AT')
```

If that didn't work, verify that the dictionary files were installed in one of the system directories (usually `/usr/share/myspell` or `/usr/share/hunspell`).


### Custom Dictionaries

If your system does not provide standard dictionaries you need to download them yourself. There are a lot of places that provide quality dictionaries. 
 
 - [SCOWL](http://wordlist.aspell.net/dicts/)
 - [OpenOffice](http://ftp.snt.utwente.nl/pub/software/openoffice/contrib/dictionaries/)
 - [LibreOffice](http://archive.ubuntu.com/ubuntu/pool/main/libr/libreoffice-dictionaries/?C=S;O=D)
 - [titoBouzout](https://github.com/titoBouzout/Dictionaries)
 - [wooorm](https://github.com/wooorm/dictionaries)

On OS-X it is [recommended](https://github.com/Homebrew/homebrew-core/blob/master/Formula/hunspell.rb#L38-L47) to put the files in `~/Library/Spelling/` or `/Library/Spelling/`. However you can also put them in your project working directory, or any of the other standard locations. If you wish to store your dictionaries somewhere else, you can make hunspell find them by setting the `DICPATH` environment variable. The `hunspell:::dicpath()` shows which locations your system searches:

```{r}
Sys.setenv(DICPATH = "/my/custom/hunspell/dir")
hunspell:::dicpath()
```

```{r, echo = FALSE, message = FALSE}
unlink(c("1406.4806v1.tar.gz", "content.tex"))
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hunspell.R
\name{hunspell}
\alias{hunspell}
\alias{hunspell_find}
\alias{en_stats}
\alias{dicpath}
\alias{hunspell_parse}
\alias{hunspell_check}
\alias{hunspell_suggest}
\alias{hunspell_analyze}
\alias{hunspell_stem}
\alias{hunspell_info}
\alias{dictionary}
\alias{list_dictionaries}
\title{Hunspell Spell Checking and Morphological Analysis}
\usage{
hunspell(text, format = c("text", "man", "latex", "html", "xml"),
  dict = dictionary("en_US"), ignore = en_stats)

hunspell_parse(text, format = c("text", "man", "latex", "html", "xml"),
  dict = dictionary("en_US"))

hunspell_check(words, dict = dictionary("en_US"))

hunspell_suggest(words, dict = dictionary("en_US"))

hunspell_analyze(words, dict = dictionary("en_US"))

hunspell_stem(words, dict = dictionary("en_US"))

hunspell_info(dict = dictionary("en_US"))

dictionary(lang = "en_US", affix = NULL, add_words = NULL,
  cache = TRUE)

list_dictionaries()
}
\arguments{
\item{text}{character vector with arbitrary input text}

\item{format}{input format; supported parsers are \code{text}, \code{latex}, \code{man},
\code{xml} and \code{html}.}

\item{dict}{a dictionary object or string which can be passed to \code{\link{dictionary}}.}

\item{ignore}{character vector with additional approved words added to the dictionary}

\item{words}{character vector with individual words to spell check}

\item{lang}{dictionary file or language, see details}

\item{affix}{file path to corresponding affix file. If \code{NULL} it is
is assumed to be the same path as \code{dict} with extension \code{.aff}.}

\item{add_words}{a character vector of additional words to add to the dictionary}

\item{cache}{speed up loading of dictionaries by caching}
}
\description{
The \code{\link{hunspell}} function is a high-level wrapper for finding spelling
errors within a text document. It takes a character vector with text (\code{plain},
\code{latex}, \code{man}, \code{html} or \code{xml} format), parses out the words
and returns a list with incorrect words for each line. It effectively combines
\code{\link{hunspell_parse}} with \code{\link{hunspell_check}} in a single step.
Other functions in the package operate on individual words, see details.
}
\details{
Hunspell uses a special dictionary format that defines which stems and affixes are
valid in a given language. The \code{\link{hunspell_analyze}} function shows how a
word breaks down into a valid stem plus affix. The \code{\link{hunspell_stem}}
function is similar but only returns valid stems for a given word. Stemming can be
used to summarize text (e.g in a wordcloud). The \code{\link{hunspell_check}} function
takes a vector of individual words and tests each one for correctness. Finally
\code{\link{hunspell_suggest}} is used to suggest correct alternatives for each
(incorrect) input word.

Because spell checking is usually done on a document, the package includes some
parsers to extract words from various common formats. With \code{\link{hunspell_parse}}
we can parse plain-text, latex and man format. R also has a few built-in parsers
such as \code{\link[tools:RdTextFilter]{RdTextFilter}} and
\code{\link[tools:SweaveTeXFilter]{SweaveTeXFilter}}, see also
\code{\link[utils:aspell]{?aspell}}.

The package searches for dictionaries in the working directory as well as in the
standard system locations. \code{\link{list_dictionaries}} provides a list of all
dictionaries it can find. Additional search paths can be specified by setting
the \code{DICPATH} environment variable. A US English dictionary (\code{en_US}) is
included with the package; other dictionaries need to be installed by the system.
Most operating systems already include compatible dictionaries with names such as
\href{https://packages.debian.org/sid/hunspell-en-gb}{hunspell-en-gb} or
\href{https://packages.debian.org/sid/myspell-en-gb}{myspell-en-gb}.

To manually install dictionaries, copy the corresponding \code{.aff} and \code{.dic}
file to \code{~/Library/Spelling} or a custom directory specified in \code{DICPATH}.
Alternatively you can pass the entire path to the \code{.dic} file as the \code{dict}
parameter. Some popular sources of dictionaries are
\href{http://wordlist.aspell.net/dicts/}{SCOWL},
\href{http://ftp.snt.utwente.nl/pub/software/openoffice/contrib/dictionaries/}{OpenOffice},
\href{http://archive.ubuntu.com/ubuntu/pool/main/libr/libreoffice-dictionaries/?C=S;O=D}{debian},
\href{https://github.com/titoBouzout/Dictionaries}{github/titoBouzout} or
\href{https://github.com/wooorm/dictionaries}{github/wooorm}.

Note that \code{hunspell} uses \code{\link{iconv}} to convert input text to
the encoding used by the dictionary. This will fail if \code{text} contains characters
which are unsupported by that particular encoding. For this reason UTF-8 dictionaries
are preferable over legacy 8-bit dictionaries.
}
\examples{
# Check individual words
words <- c("beer", "wiskey", "wine")
correct <- hunspell_check(words)
print(correct)

# Find suggestions for incorrect words
hunspell_suggest(words[!correct])

# Extract incorrect from a piece of text
bad <- hunspell("spell checkers are not neccessairy for langauge ninja's")
print(bad[[1]])
hunspell_suggest(bad[[1]])

# Stemming
words <- c("love", "loving", "lovingly", "loved", "lover", "lovely", "love")
hunspell_stem(words)
hunspell_analyze(words)

# Check an entire latex document
tmpfile <- file.path(tempdir(), "1406.4806v1.tar.gz")
download.file("https://arxiv.org/e-print/1406.4806v1", tmpfile,  mode = "wb")
untar(tmpfile, exdir = tempdir())
text <- readLines(file.path(tempdir(), "content.tex"), warn = FALSE)
bad_words <- hunspell(text, format = "latex")
sort(unique(unlist(bad_words)))

# Summarize text by stems (e.g. for wordcloud)
allwords <- hunspell_parse(text, format = "latex")
stems <- unlist(hunspell_stem(unlist(allwords)))
words <- head(sort(table(stems), decreasing = TRUE), 200)
}
