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

