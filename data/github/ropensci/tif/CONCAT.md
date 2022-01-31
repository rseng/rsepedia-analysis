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
## tif: Text Interchange Formats

[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/statsmaths/tif?branch=master&svg=true)](https://ci.appveyor.com/project/statsmaths/cleanNLP) [![Travis-CI Build Status](https://travis-ci.org/statsmaths/cleanNLP.svg?branch=master)](https://travis-ci.org/ropensci/tif)

This package describes and validates formats for storing
common object arising in text analysis as native R objects.
Representations of a text corpus, document term matrix, and
tokenized text are included. The tokenized text format is
extensible to include other annotations. There are two versions
of the corpus and tokens objects; packages should accept
both and return or coerce to at least one of these.

## Installation

You can install the development version using devtools:

```{r}
devtools::install_github("ropensci/tif")
```

## Usage

The package can be used to check that a particular object is in a valid 
format. For example, here we see that the object `corpus` is a valid corpus
data frame:

```{r}
library(tif)
corpus <- data.frame(doc_id = c("doc1", "doc2", "doc3"),
                     text = c("Aujourd'hui, maman est morte.",
                      "It was a pleasure to burn.",
                      "All this happened, more or less."),
                     stringsAsFactors = FALSE)

tif_is_corpus_df(corpus)
```
```
TRUE
```

The package also has functions to convert between the list and data frame
formats for corpus and token object. For example:

```{r}
tif_as_corpus_character(corpus)
```
```
                              doc1                               doc2 
   "Aujourd'hui, maman est morte."       "It was a pleasure to burn." 
                              doc3 
"All this happened, more or less." 
```

Note that extra meta data columns will be lost in the conversion from a data
frame to a named character vector.

## Details

This package describes and validates formats for storing
common object arising in text analysis as native R objects.
Representations of a text corpus, document term matrix, and
tokenized text are included. The tokenized text format is
extensible to include other annotations. There are two versions
of the corpus and tokens objects; packages should accept and return
at least one of these.

**corpus** (data frame) - A valid corpus data frame object
is a data frame with at least two columns. The first column
is called doc_id and is a character vector with UTF-8 encoding. Document
ids must be unique. The second column is called text and
must also be a character vector in UTF-8 encoding. Each
individual document is represented by a single row in
the data frame. Addition document-level metadata columns
and corpus level attributes are allowed but not required.

**corpus** (character vector) - A valid character vector corpus
object is an character vector with UTF-8 encoding. If it has
names, this should be a unique character also in UTF-8
encoding. No other attributes should be present.

**dtm** - A valid document term matrix is a sparse matrix with
the row representing documents and columns representing
terms. The row names is a character vector giving the
document ids with no duplicated entries. The column
names is a character vector giving the terms of the
matrix with no duplicated entries. The sparse matrix
should inherit from the Matrix class dgCMatrix.

**tokens** (data frame) - A valid data frame tokens
object is a data frame with at least two columns. There must be
a column called doc_id that is a character vector
with UTF-8 encoding. Document ids must be unique.
There must also be a column called token that must also be a
character vector in UTF-8 encoding.
Each individual token is represented by a single row in
the data frame. Addition token-level metadata columns
are allowed but not required. 

**tokens** (list) - A valid corpus tokens object is (possibly
named) list of character vectors. The character vectors, as
well as names, should be in UTF-8 encoding. No other
attributes should be present in either the list or any of its
elements.
# tif 0.3.0

* Further discussion has lead us to simplify the corpus and token data frame
formats. The doc_id, text, and token columns can be in any position within the
data frame.

# tif 0.2.0

* After a round of input for the initial version of the specification,
we decided to allow two formats for corpus and tokens objects. In addition
to the original data frame variants there is a character vector corpus
object and a list-based tokens object. Converts between the various types
are now included in the package.

### New Functions

* `tif_is_corpus_character` returns TRUE or FALSE for whether the input
is a valid character vector corpus object.

* `tif_is_tokens_list` returns TRUE or FALSE for whether the input
is a valid list-based tokens object.

* `tif_as_corpus_character` takes a valid tif corpus object and returns
a character vector corpus object.

* `tif_as_corpus_df` takes a valid tif corpus object and returns
a data frame corpus object.

* `tif_as_tokens_character` takes a valid tif tokens object and returns
a list-based tokens object.

* `tif_as_tokens_df` takes a valid tif tokens object and returns
a list-based tokens object.

### Renamed Functions

* The old validate functions have been renamed `tif_is_corpus_df`,
`tif_is_dtm` and `tif_is_tokens_df`. This is more in line with base-R
functions and separates the "df" version of the corpus and tokens from
the alternative new forms.

# tif 0.1.0

* This is the initial implementation of the ideas discussed at
the rOpenSci Text Workshop from 21-22 April 2017.

### New Functions

* `tif_corpus_validate` returns TRUE or FALSE for whether the input
is a valid corpus object.

* `tif_dtm_validate` returns TRUE or FALSE for whether the input is
a valid document corpus object.

* `tif_tokens_validate` returns TRUE or FALSE for whether the input is
a valid tokens object.

### Known issues

* do not yet have a test suite for the package

* encoding checkin is not yet working

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validators.R
\name{tif_is_dtm}
\alias{tif_is_dtm}
\title{Validate Document Term Matrix Object}
\usage{
tif_is_dtm(dtm, warn = FALSE)
}
\arguments{
\item{dtm}{a document term matrix object to test
the validity of}

\item{warn}{logical. Should the function produce a
verbose warning for the condition for which
the validation fails. Useful for testing.}
}
\value{
a logical vector of length one indicating
              whether the input is a valid document term
              matrix
}
\description{
A valid document term matrix is a sparse matrix with
the row representing documents and columns representing
terms. The row names is a character vector giving the
document ids with no duplicated entries. The column
names is a character vector giving the terms of the
matrix with no duplicated entries. The spare matrix
should inherit from the Matrix class dgCMatrix.
}
\details{
The tests are run sequentially and the function returns,
with a warning if the warn flag is set, on the first test
that fails. We use this implementation because some tests
may fail entirely or be meaningless if the prior ones are
note passed. For example, if the dtm object is not a matrix
it may not contain row or column names.
}
\examples{
#' @importFrom Matrix Matrix
dtm <- Matrix::Matrix(0, ncol = 26, nrow = 5, sparse = TRUE)
colnames(dtm) <- LETTERS
rownames(dtm) <- sprintf("doc\%d", 1:5)

tif_is_dtm(dtm)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validators.R
\name{tif_is_tokens_list}
\alias{tif_is_tokens_list}
\title{Validate Tokens List Object}
\usage{
tif_is_tokens_list(tokens, warn = FALSE)
}
\arguments{
\item{tokens}{a tokens object to test for validity}

\item{warn}{logical. Should the function produce a
verbose warning for the condition for which
the validation fails. Useful for testing.}
}
\value{
a logical vector of length one indicating
               whether the input is a valid tokens
}
\description{
A valid corpus tokens object is (possibly named) list of
character vectors. The character vectors, as well as
names, should be in UTF-8 encoding. No other attributes
should be present in either the list or any of its elements.
}
\details{
The tests are run sequentially and the function returns,
with a warning if the warn flag is set, on the first test
that fails. We use this implementation because some tests
may fail entirely or be meaningless if the prior ones are
note passed.
}
\examples{
tokens <- list(doc1 = c("aujourd'hui", "maman", "est", "morte"),
               doc2 = c("it", "was", "a", "pleasure", "to", "burn"),
               doc3 = c("all", "this", "happened", "more", "or", "less"))
tif_is_tokens_list(tokens)

names(tokens) <- c("doc1", "doc2", "doc3")
tif_is_tokens_list(tokens)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkg.R
\docType{package}
\name{tif-package}
\alias{tif}
\alias{tif-package}
\title{tif: Text Interchange Formats}
\description{
This package describes and validates formats for storing
common object arising in text analysis as native R objects.
Representations of a text corpus, document term matrix, and
tokenized text are included. The corpus and tokens objects
have multiple valid formats. Packages compliant with the
tif proposal should accept all valid formats and should
directly return, or provide conversion functions, for
converting outputs into at least one of the formats (when
applicable). The tokenized text format is extensible to
include other annotations such as part of speech tags and
named entities.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/tif}
  \item Report bugs at \url{http://github.com/ropensci/tif/issues}
}

}
\author{
\strong{Maintainer}: Taylor Arnold \email{taylor.arnold@acm.org}

Authors:
\itemize{
  \item Ken Benoit \email{k.r.benoit@lse.ac.uk}
  \item Lincoln Mullen \email{lmullen@gmu.edu }
  \item Adam Obeng \email{contact@adamobeng.com}
  \item rOpenSci Text Workshop Participants (2017)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coercion.R
\name{tif_as}
\alias{tif_as}
\alias{tif_as_corpus_character}
\alias{tif_as_corpus_character.default}
\alias{tif_as_corpus_character.character}
\alias{tif_as_corpus_character.data.frame}
\alias{tif_as_corpus_df}
\alias{tif_as_corpus_df.default}
\alias{tif_as_corpus_df.character}
\alias{tif_as_corpus_df.data.frame}
\alias{tif_as_tokens_df}
\alias{tif_as_tokens_df.default}
\alias{tif_as_tokens_df.list}
\alias{tif_as_tokens_df.data.frame}
\alias{tif_as_tokens_list}
\alias{tif_as_tokens_list.default}
\alias{tif_as_tokens_list.list}
\alias{tif_as_tokens_list.data.frame}
\title{Coerce Between tif Object Specifications}
\usage{
tif_as_corpus_character(corpus)

\method{tif_as_corpus_character}{default}(corpus)

\method{tif_as_corpus_character}{character}(corpus)

\method{tif_as_corpus_character}{data.frame}(corpus)

tif_as_corpus_df(corpus)

\method{tif_as_corpus_df}{default}(corpus)

\method{tif_as_corpus_df}{character}(corpus)

\method{tif_as_corpus_df}{data.frame}(corpus)

tif_as_tokens_df(tokens)

\method{tif_as_tokens_df}{default}(tokens)

\method{tif_as_tokens_df}{list}(tokens)

\method{tif_as_tokens_df}{data.frame}(tokens)

tif_as_tokens_list(tokens)

\method{tif_as_tokens_list}{default}(tokens)

\method{tif_as_tokens_list}{list}(tokens)

\method{tif_as_tokens_list}{data.frame}(tokens)
}
\arguments{
\item{corpus}{valid tif corpus object to coerce}

\item{tokens}{valid tif tokens object to coerce}
}
\description{
These functions convert between the various valid
formats for corpus and tokens objects. By using these
in other packages, maintainers need to only handle
whichever specific format they would like to work
with, but gain the freedom to output (or convert
into) the one most suited to their package's paradigm.
}
\details{
No explicit checking is done on the input; the output
is guaranteed to be valid only if the input is a valid
format. In fact, we make an effort to not modify an
object that appears to be in the required format already
due to R's copy on modify semantics.
}
\examples{
# coerce corpus object
corpus <- c("Aujourd'hui, maman est morte.",
            "It was a pleasure to burn.",
            "All this happened, more or less.")
names(corpus) <- c("Camus", "Bradbury", "Vonnegut")

new <- tif_as_corpus_df(corpus)
new
tif_as_corpus_character(new)

# coerce tokens object
tokens <- list(doc1 = c("aujourd'hui", "maman", "est", "morte"),
               doc2 = c("it", "was", "a", "pleasure", "to", "burn"),
               doc3 = c("all", "this", "happened", "more", "or", "less"))

new <- tif_as_tokens_df(tokens)
new
tif_as_tokens_list(new)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validators.R
\name{tif_is_corpus_df}
\alias{tif_is_corpus_df}
\title{Validate Corpus Data Frame Object}
\usage{
tif_is_corpus_df(corpus, warn = FALSE)
}
\arguments{
\item{corpus}{a corpus object to test for validity}

\item{warn}{logical. Should the function produce a
verbose warning for the condition for which
the validation fails. Useful for testing.}
}
\value{
a logical vector of length one indicating
               whether the input is a valid corpus
}
\description{
A valid data frame corpus object is an object that
least two columns. One column must be called doc_id
and be a character vector with UTF-8 encoding. Document
ids must be unique. There must also be a column called text
and must also be a character vector in UTF-8 encoding. Each
individual document is represented by a single row in
the data frame. Addition document-level metadata columns
and corpus level attributes are allowed but not required.
}
\details{
The tests are run sequentially and the function returns,
with a warning if the warn flag is set, on the first test
that fails. We use this implementation because some tests
may fail entirely or be meaningless if the prior ones are
note passed. For example, if the corpus object does not
have a variable named "text" it does not make sense to
check whether this column is a character vector.
}
\examples{
corpus <- data.frame(doc_id = c("doc1", "doc2", "doc3"),
                     text = c("Aujourd'hui, maman est morte.",
                      "It was a pleasure to burn.",
                      "All this happened, more or less."),
                     stringsAsFactors = FALSE)

tif_is_corpus_df(corpus)

corpus$author <- c("Camus", "Bradbury", "Vonnegut")
tif_is_corpus_df(corpus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validators.R
\name{tif_is_corpus_character}
\alias{tif_is_corpus_character}
\title{Validate Corpus Character Vector Object}
\usage{
tif_is_corpus_character(corpus, warn = FALSE)
}
\arguments{
\item{corpus}{a corpus object to test for validity}

\item{warn}{logical. Should the function produce a
verbose warning for the condition for which
the validation fails. Useful for testing.}
}
\value{
a logical vector of length one indicating
               whether the input is a valid corpus
}
\description{
A valid character vector corpus object is an character
vector with UTF-8 encoding. If it has names, this should
be a unique character also in UTF-8 encoding. No other
attributes should be present.
}
\details{
The tests are run sequentially and the function returns,
with a warning if the warn flag is set, on the first test
that fails. We use this implementation because some tests
may fail entirely or be meaningless if the prior ones are
note passed.
}
\examples{
corpus <- c("Aujourd'hui, maman est morte.",
            "It was a pleasure to burn.",
            "All this happened, more or less.")

tif_is_corpus_character(corpus)

names(corpus) <- c("Camus", "Bradbury", "Vonnegut")
tif_is_corpus_character(corpus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validators.R
\name{tif_is_tokens_df}
\alias{tif_is_tokens_df}
\title{Validate Tokens Data Frame Object}
\usage{
tif_is_tokens_df(tokens, warn = FALSE)
}
\arguments{
\item{tokens}{a tokens object to test for validity}

\item{warn}{logical. Should the function produce a
verbose warning for the condition for which
the validation fails. Useful for testing.}
}
\value{
a logical vector of length one indicating
               whether the input is a valid tokens object
}
\description{
A valid tokens data frame object is a data frame or an
object that inherits a data frame. It has no row names
and has at least two columns. It must a contain column called
doc_id that is a character vector with UTF-8 encoding.
Document ids must be unique. It must also contain a column called
token that must also be a character vector in UTF-8 encoding.
Each individual token is represented by a single row in
the data frame. Addition token-level metadata columns
are allowed but not required.
}
\details{
The tests are run sequentially and the function returns,
with a warning if the warn flag is set, on the first test
that fails. We use this implementation because some tests
may fail entirely or be meaningless if the prior ones are
note passed. For example, if the tokens object does not
have a variable named "doc_id" it does not make sense to
check whether this column is a character vector.
}
\examples{
tokens <- data.frame(doc_id = c("doc1", "doc1", "doc1", "doc1",
                                "doc2",  "doc2", "doc2", "doc2",
                                "doc2", "doc2", "doc3", "doc3",
                                "doc3", "doc3", "doc3", "doc3"),
                     token = c("aujourd'hui", "maman", "est",
                               "morte", "it", "was", "a", "pleasure",
                               "to", "burn", "all", "this", "happened",
                               "more", "or", "less"),
                     stringsAsFactors = FALSE)

tif_is_tokens_df(tokens)

tokens$pos <- "NOUN"
tokens$NER <- ""
tokens$sentiment <- runif(16L)
tif_is_tokens_df(tokens)
}
