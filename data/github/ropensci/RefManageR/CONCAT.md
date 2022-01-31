---
title: 'RefManageR: Import and Manage BibTeX and BibLaTeX References in R'
tags:
  - R
  - reference management
  - BibLaTeX
  - document generation
authors:
 - name: Mathew W. McLean
   orcid: 0000-0002-7891-9645
   affiliation: 1
affiliations:
 - name: University of Technology Sydney
   index: 1
date: 25 May 2017
bibliography: paper.bib
---

# Summary

RefManageR provides tools for importing and working with bibliographic
references.  It greatly enhances the bibentry class in R by providing
a class BibEntry which stores BibTeX and BibLaTeX references, supports
UTF-8 encoding, and can be easily searched by any field, by date
ranges, and by various formats for name lists (author by last names,
translator by full names, etc.) using R's person class
[@hornik2012]. Entries can be updated, combined, sorted, printed in a
number of styles, and exported. BibTeX [@bibtex] and BibLaTeX
[@biblatex] .bib files can be read into R and converted to BibEntry
objects.  Interfaces to NCBI's Entrez, CrossRef, and Zotero are
provided for importing references and references can be created from
locally stored PDFs using Poppler.  Includes functions for citing and
generating a bibliography with hyperlinks for documents prepared with
RMarkdown or RHTML.  A vignette is available to further demonstrate
all functionality [@refmanager].
 
BibEntry objects can be created directly in R or a .bib file can be
read into R to create the object.  Tools are
provided for importing references from Crossref, Zotero, Google
Scholar, and PDFs and for looking up PubMed ID's and DOIs.

BibEntry objects may be searched and indexed by field values, name
lists, keys, dates, date ranges, etc.  They can be printed in a number
of formats (e.g. text, html) and most of the base bibliography styles
available with BibLaTeX (e.g. alphabetic, numeric, authortitle, and
authoryear).  All sorting methods for bibliographies available in the
BibLaTeX LaTeX package have been implemented.  A function is provided
to convert a BibEntry object to a character vector containing lines of
a BibTeX or BibLaTeX file, converting fields, entry types and
expanding crossreferences to coerce BibLaTeX entries to BibTeX if
requested.  The results can also be written to a file.

Citations can be gerenated in a number of styles using one of the
available functions for citations.  A list of references can be
printed based on the works the user has cited thus far in their
document.  The citations and bibliography can be printed including
hyperlinks using either the R Markdown or R HTML formats.  A function
is provided to open electronic copies of references in a PDF viewer or
web browser.  A simple function is provided for setting default
formatting options throughout the session.

# References
RefManageR
========
[![R-CMD-check](https://github.com/ropensci/refmanager/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/refmanager/actions)
[![AppVeyor Build Status](http://ci.appveyor.com/api/projects/status/github/ropensci/RefManageR?branch=master&svg=true)](http://ci.appveyor.com/project/ropensci/RefManageR)
[![Coverage Status](https://coveralls.io/repos/github/ropensci/RefManageR/badge.svg?branch=master)](https://coveralls.io/github/ropensci/RefManageR?branch=master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/RefManageR)](https://cran.r-project.org/package=RefManageR)
[![CRAN_Download_Badge](http://cranlogs.r-pkg.org/badges/RefManageR)](https://cran.r-project.org/package=RefManageR)
[![](https://badges.ropensci.org/119_status.svg)](https://github.com/ropensci/software-review/issues/119)
[![](https://joss.theoj.org/papers/10.21105/joss.00338/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00338)
`RefManageR` provides tools for importing and working with
bibliographic references.  It greatly enhances the `bibentry` class by
providing a class `BibEntry` which stores `BibTeX` and `BibLaTeX` references,
supports `UTF-8` encoding, and can be easily searched by any field, by date
ranges, and by various formats for name lists (author by last names,
translator by full names, etc.). Entries can be updated, combined, sorted,
printed in a number of styles, and exported. `BibTeX` and `BibLaTeX` `.bib` files
can be read into `R` and converted to `BibEntry` objects.  Interfaces to
`NCBI Entrez`, `CrossRef`, and `Zotero` are provided for importing references and
references can be created from locally stored `PDF` files using `Poppler`.  Includes
functions for citing and generating a bibliography with hyperlinks for
documents prepared with `RMarkdown` or `RHTML`.

Please see the [vignette](https://arxiv.org/pdf/1403.2036v1)
for an introduction and [NEWS](https://github.com/ropensci/RefManageR/blob/master/inst/NEWS.md)
for the latest changes.

To install the latest version from `GitHub`:

```
install.packages("remotes")
remotes::install_github("ropensci/RefManageR")
```
[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# Contributing to RefManageR

The following is modified from the [plotly contributing
guidelines](https://github.com/ropensci/plotly/blob/master/CONTRIBUTING.md).

## Opening issues

When reporting an issue please provide a reproducible example.  You
may find the package [reprex](https://github.com/tidyverse/reprex)
useful for this.

## Development Guidelines

If you'd like to contribute changes to RefManageR, see [the GitHub
flow](https://guides.github.com/introduction/flow/index.html) for
proposing, submitting, reviewing, and accepting changes. If you
haven't done this before, Hadley Wickham provides a nice overview of
git (<http://r-pkgs.had.co.nz/git.html>), as well as best practices
for submitting pull requests
(<http://r-pkgs.had.co.nz/git.html#pr-make>). I also recommend using
his style guide when writing code
(<http://adv-r.had.co.nz/Style.html>).  If your pull request fixes a
bug, or implements a new feature, it's a good idea to write a test
(<http://r-pkgs.had.co.nz/tests.html>) to demonstrate it's working.
Changes in Version 1.3.0 (2019-10-30)
========================================================

* Package `bibtex` moved to Suggests in `DESCRIPTION` due to it
currently being orphaned on CRAN. Work is underway by the ROpenSci
team to rectify this. The package can still be installed from
[GitHub](https://github.com/ROpenSci/bibtex). In the event that
`bibtex` is not installed, the functions `ReadBib()`,
`GetBibEntryWithDOI()`, `ReadCrossRef()`, and `ReadZotero()` 
throw an appropriate message asking the user to install `bibtex` and invisibly
return `NULL`.
* The old CrossRef API can no longer be used. If `use.old.api` is set
  to \code{TRUE} in \code{ReadCrossRef()}, it will be ignored with a warning.
* `GetDOIs()` had to be removed due to changes to the CrossRef API. It
  will hopefully return in the next release.
  

Changes in Version 1.2.13 (2019-04-03)
========================================================

## BUG FIXES

* When working in single-byte locales, the `print` method for
`BibEntry` objects is more robust against accented characters being
converted to incorrect ones when `bib.style = "authoryear"`. 
Additionally, for this style, a period could be removed
from the last initial in the first author's given name when
`first.inits = TRUE`. This has been corrected.

Changes in Version 1.2.12 (2019-04-02)
========================================================

* The serial comma is now used when formatting name lists

## BUG FIXES

* Fix issue that could lead to Unicode characters being
converted to latin1 when printing (h/t joaochenriques #62)
* Fix issue with extracting DOI from CrossRef results
in `GetDOIs`
* Fixes for `ReadCrossRef` when `use.old.api` is `TRUE`. Scores
are sometimes not returned by the API call and when this occurs
the entries will now be added to the output `BibEntry` object with
a message indicating that no score was available.
* A comma no longer appears before "et al." when the `max.names`
options is set to `1` (h/t davidaknowles #56)

Changes in Version 1.2.8 (2018-12-10)
========================================================

## BUG FIXES

* Fix extraction of citation counts in `ReadGS` that was occasionally 
causing errors from some `scholar.id`s (h/t Miao Sun, #59)
* Fix printing of thesis and report entries types in authoryear style
when "type" field missing (h/t Hugo Grunson, #58)
* Fix for `PrintBibliography` for the case of `BibEntry` objects with
  a single entry, an NA value could appear next to the year in the
  output (#60)
* Fix issue that could lead to Unicode characters being
converted to latin1 when printing (h/t joaochenriques #62)

Changes in Version 1.2.2 (2018-05-31)
========================================================

## BUG FIXES

* `GetPubMedByID` is better at extracting years and months
from the results returned by NCBI Entrez (h/t Dale Steele #52)
* The `as.data.frame`method for `BibEntry` objects now correctly
handles the case of a single entry with name list fields containing multiple 
names (h/t Damon Bayer #51)

Changes in Version 1.2.0 (2018-04-24)
========================================================

## NEW FEATURES

* `+.BibEntry` and `merge.BibEntry` gain an argument ignore.case,
which defaults to `BibOptions()$ignore.case` (`TRUE`) so that case is 
ignore when checking for duplicate fields (h/t Justin Calabrese #47)
* Improved warning message when printing entries with unknown 
LaTeX macros (the entry key is now included). (h/t Justin Calabrese #49)
* The entry key is now included in warning messages when entries are 
missing fields and `BibOptions()$check.entries == "warn"` (h/t Justin
Calabrese #48)

## BUG FIXES

* Entries are now only checked once to ensure all required
  fields are present in `ReadBib`

Changes in Version 1.1.0 (2018-04-02)
========================================================

* `PrintBibliography` gains parameters "start" and "end"
to allow for printing only a subset of all cited entries from
a BibEntry object (h/t Joseph Casillas #45, #46)


Changes in Version 1.0.0 (2018-02-19)
========================================================

* Use https for all links (h/t Katrin Leinweber)
* Use preferred DOI resovler (h/t Katrin Leinweber)
* Add support for latex macro ast for asterisks (h/t Melinda Higgins)

Changes in Version 0.14.25 (2017-12-25)
========================================================

## BUG FIXES

* Fix `GetDOIs` to use https
* Fix download of bibliographic info from DOI in `ReadPDFs`

Changes in Version 0.14.23 (2017-11-12)
========================================================

## BUG FIXES

* Fix writing of BibEntry object to stdout in
`WriteBib` (h/t Stephane Plaisance)
* `ReadBib` won't add an attribute "strings" if there
are none present in read bib file (h/t Stephane Plaisance)

Changes in Version 0.14.21 (2017-09-05)
========================================================

## BUG FIXES

* Fix deletion of temporary file if user supplies a DOI to 
`ReadCrossRef` (h/t Ben Raymond)


Changes in Version 0.14.20 (2017-07-27)
========================================================

## NEW FEATURES

* Documentation example improvements
* Improve error handling for API query functions
* Package peer reviewed and accepted by rOpenSci (h/t Noam Ross, 
Carl Boettiger, and Amelia McNamara)

## BUG FIXES

* Remove missing plot from Rhtml vignette
* URL field returned by `GetBibEntryWithDOI` is now decoded properly
* Fix hyperlinks from bibliography to citations in vignettes
* Remove a incorrect message occasionally output from the addition
operator for `BibEntry` objects

Changes in Version 0.14.12 (2017-06-30)
========================================================

## NEW FEATURES

* Package now uses httr, xml2, jsonlite packages instead of RCurl,
XML, RJSONIO for scaffolding
* No more R CMD check NOTE regarding foreign function call to bibtex
(h/t Romain Francois)

## BUG FIXES

* Fix printing when `BibOptions(style = 'yaml)`
* Remove invalid character in inst/Bib/RJC.bib
* Correct parsing of interval dates when creating unique labels for
authoryear style citations
* `c.BibEntry` throws an error if not all objects are `bibentry` objects
* Fix typos in documentation
* Literal ampersands are now printed correctly (not as '\&') (h/t Yue Hu)
* Ensure BibTeX month macros are processed properly by lubridate
in non-English locales (h/t Sergio Oller)

Changes in Version 0.13.4 (2017-04-25)
========================================================

## BUG FIXES

* Unescape special characters in URL fields returned by CrossRef (h/t
  Michael Schubert)
* Remove square brackets from custom entry type names (h/t Hugh
  Parsonage)

Changes in Version 0.13.1 (2016-11-14)
========================================================

## BUG FIXES

* Feature involving `LaTeX` macros added in package version 0.12.0 can
only be used for R 3.3.z and higher; this corrects cause of failed
checks on R 3.2.z


Changes in Version 0.13.0 (2016-11-09)
========================================================

## BUG FIXES

* Updated calls to NCBI Entrez for functions `ReadPubMed`,
`GetPubMedByID`, etc.  to use https as now required by NCBI (h/t Dale
Steele and Anthony Crane)
* Change reference to www.omegahat.org to www.omegahat.net (h/t Kurt Hornik)
* Documentation for `ReadPubMed` is updated to reflect that the
default number of entries returned (controlled by the argument
`retmax`) is 20 (h/t Dale Steele)

Changes in Version 0.12.0 (2016-09-30)
========================================================

## NEW FEATURES

* Some `LaTeX` macros unknown to R are now defined as macros in the
package, and will be parsed using `macros` arg in `tools::parse_Rd`
(assuming `getRversion() >= "3.2.0"` Note: corrected in 0.13.1 to be
`getRversion() >= "3.3.0"`)

## BUG FIXES

* Parse `LaTeX` macro `\textquotesingle` in author names (h/t Bill Denney)
* Avoid "Request-URI too large" errors in GetPubMedByID if requesting a large number
of IDs (h/t Maurits Evers)

Changes in Version 0.11.0 (2016-09-10)
========================================================

## NEW FEATURES

* `ReadCrossRef` now uses the
[newer CrossRef API](https://github.com/CrossRef/rest-api-doc/blob/master/rest_api.md) 
and gains arguments `filter` and `offset` to use with the new API; an
additional argument `use.old.api` is added if the user wishes to use
the old API (h/t Carl Boettiger)
* `ReadCrossRef` now parses the results returned by CrossRef to create
the `BibEntry` object when using the new API; for the old API (and
hence, older versions of the package) the query only returns DOIs and
`ReadCrossRef` would then use the DOIs to request the corresponding
BibTeX entries from CrossRef (i.e. less HTTP requests when using the
new API)

## BUG FIXES

* Fix generation of entry keys when the word used from the title for key
generation contains a non-ascii character (h/t Mark Johnson)
* RefManageR will no longer hang due to a bug in `tools::latexToUtf8`
([PR\#17138](https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=17138)) that
is occasionally encountered when that function processes an unknown
macro (h/t Eric Bryant)
* Entries with no title field can now be printed without error when
`BibOptions()$check.entries` is *not* set to "error" (default is "error")

Changes in Version 0.10.15 (2016-06-06)
========================================================

## BUG FIXES

* Removed unnecessary use of local/parent.frame; fixes execution with
bytecode compiler (h/t Tomas Kalibera)

Changes in Version 0.10.12 (2016-03-25)
========================================================

## BUG FIXES

* Fixed broken test involving `ReadPDFs` due to changed URL (h/t Kurt
  Hornik)
* `as.data.frame.BibEntry` works for length one BibEntry with multiple
authors  (h/t Dale Steele)
* Use `httr::GET` to fix `ReadGS`
* Fixed broken tests in `test-authors.R` owing to changes to `person` class
* Fixed `ReadCrossRef` tests and error message
* Fixed printing for authoryear style (h/t Joseph Casillas)
* Name list fields (author, editor, etc.) provided to the function `BibEntry` are
now properly parsed when specified as they would be in BibTeX/BibLaTeX;
e.g. `author = "Smith, Bob A. and Doe, Jane"`. 

Changes in Version 0.10.5 (2016-01-03)
========================================================

* The 'key' field in `BibEntry` objects is now always enforced to be unique
* `as.data.frame.BibEntry` is faster and now works if duplicate keys are
present; keys in  (h/t Dale Steele)
* Fix for `ReadCrossRef` if downloaded BibTeX had leading whitespace
(h/t Carl Boettiger)
* `useBytes = TRUE` used for all calls to `grep`, `sub`, etc. (h/t
  HI&RH Lord Ripley of England)
* remove use of deprecated function `lubridate::new_interval`
* updated URL for the BibLaTeX manual
* Fix test in test-search.R that broke because of new year (h/t HI&RH
  Lord Ripley of England)
* add additional functions from utils and stats to NAMESPACE

Changes in Version 0.9.0 (2015-06-10)
========================================================

* Use `bibtex >= 0.4.0.9000` function `do_read_bib` to avoid
`.External` call and `R check` note (request of HI&RH Lord Ripley of
England)

Changes in Version 0.8.63 (2015-06-08)
========================================================

## NEW FEATURES

* Improve parsing of dates in `ReadPDFs`
* Citations using `Cite` family of functions can now be `pandoc`
style, e.g. `[@abibkey]` by setting `BibOptions(cite.style =
"pandoc")` (h/t Dale Steele)
* Added note about locales when parsing string 'month' fields to
`ReadBib` help page (h/t Dieter Menne)

## BUG FIXES

* Fixed merging `BibEntry` objects by multiple fields when no duplicates
* `open.BibEntry` fixed to not use partial matching of field names;
e.g. an error would occur if the specified entry had a 'urldate'
field, but no 'url' field
* `open.BibEntry` will `message` and not throw error if entry cannot
  be opened
* Fixes for `ReadPDFs` when argument `use.metadata` is `FALSE`
* Fix for `ReadPDFs` when when reading *one* file which is a JSTOR pdf
* Fix for sorting by volume (`BibOptions(sorting = "anyvt")` and
  `BibOptions(sorting = "nyvt")`
* Fix for sorting by label (`BibOptions()sorting` equal to  "anyvt" or "anyt")
* `GetBibEntryWithDOI.R` will not `stop` if an error occurs
downloading any of the DOIs (e.g., if one entry in the `doi` vector
has a typo and the rest are valid)

Changes in Version 0.8.52 (2015-01-26)
========================================================

## NEW FEATURES

* `GetPubMedByID`: Now returns some additional fields including
'month' and 'issn' for articles; will print a warning if PubMed does
not return the complete list of authors; will use the name of a
collective if one is available and the individual authors are missing
(h/t Dale Steele)

## BUG FIXES

* `ReadBib`: If a name list field in an entry cannot be parsed in the bib file, the
entry will be ignored, but the rest of the file will still be processed and
returned. In the past, this caused an error and no output would be returned.
* 'Book' entries will now be parsed correctly by `GetPubMedByID` (h/t
  Dale Steele)
* Fix error/warning messages when entry is missing required fields (bug introduced in
Version 0.8.45)
* Name lists containing a comma in braces will now be parsed
correctly, e.g. "Buchalter, Louis and {Murder, Inc.} and Anastasia,
Albert"

Changes in Version 0.8.45 (2014-12-29)
========================================================

## BUG FIXES

* `ReadCrossRef` now correctly handles the small number of cases where
BibTeX information cannot be obtained for a particular DOI, which
resulted in 'stack imbalance' warnings and no results being returned
(h/t Norman L Guinasso Jr).
* `ReadGS` fixed to account for changes to Google "API" (h/t Norman L
  Guinasso Jr).
* Improved parsing for BibTeX format names ending with a '}' (h/t Henrik Bengtsson).
* Printing references with `style = "html"` would not always add an
opening <cite> tag when `bib.style = "numeric"` or `bib.style =
"alphabetic"` (h/t Henrik Bengtsson).
* `format.BibEntry` would ignore the `.style` argument if called
directly by the user.  Note, this function should normally not need to
be called directly. (h/t Henrik Bengtsson)

Changes in Version 0.8.40 (2014-10-28)
========================================================

## NEW FEATURES

* Improved formatting of citation given to CrossRef for increased
chances of finding matches with `GetDOIs` function (h/t Erich
Studerus)
* Additional parsing of 'month' field to accomodate days and ranges of
days and months.  Example bib entries that will be parsed correctly
include `month = jun # "/" # jul`, `month = "20~" # jan`, `month =
"20==25~" # dec`, `month = "10~" # jan # "/" # feb` (request of
Stephen Eglen)
* Added argument 'group' to `ReadZotero` for specifying a groupID to
query a group library instead of a user library (h/t Greg Blomquist).

## BUG FIXES

* DOI's hyperlinks in Markdown format are now correct (h/t Stephen Eglen)
* `print.BibEntry` with `BibOptions(style = "Biblatex")` fixed (h/t
  Artem Klevtsov)
* `unlist.BibEntry` and `RelistBibEntry` now retains `@strings` and
`mheader` and `mfooter` attributes (see ?BibEntry) if they are present

Changes in Version 0.8.34 (2014-08-18)
========================================================

## NEW FEATURES

* Added function `GetDOIs` which searches CrossRef for DOIs for the
citations stored in a `BibEntry` object

## BUG FIXES

* `ReadCrossRef` fixed to account for change to CrossRef API
  endpoint. (h/t Carl Boettiger)
* Abstracts returned by NCBI Entrez can be multiple parts.  This is
now handled correctly and the complete abstract will be returned in
the 'abstract' field. (h/t Erich Studerus)
* DOI's were too naively extracted from NCBI Entrez results, resulting
in some entries having 'doi' fields with length greater than one.  Now
fixed. (h/t Erich Studerus)

Changes in Version 0.8.32 (2014-08-14)
========================================================

## NEW FEATURES

* Functions for interacting with NCBI Entrez return abstract of each
  article (request of Erich Studerus)

## BUG FIXES

* `print.BibEntry` with `BibOptions(style = "citation")` now works properly
* Examples calling web resources should no longer upset the check farm
  (h/t HI&RH Lord Ripley of England)


Changes in Version 0.8.3 (2013-07-30)
=========================================================

## NEW FEATURES

* `as.BibEntry` will create entry key if given a `bibentry` object
with no key.  Useful when citing packages with `citation`.
* PrintBibliography and Cite functions (Cite, Citet, etc.) accept
`bibentry` objects in addition to `BibEntry` objects.
* `$<-.BibEntry` will now accept a single person object, so that a
single author in a multi-author entry may be updated.  An example may
be found at `help("$<-.BibEntry")`.  (h/t Carl Boettiger)

## BUG FIXES

* validated html
* changed example for `WriteBib` that occasionally failed check

Changes in Version 0.8.2 (2013-06-01)
=========================================================

## BUG FIXES

* Cite functions work if a specified entry has no key.  Note that keys
  should always be provided for all entries as they are required for
  all entries in a BibLaTeX bib file (h/t Carl Boettiger)
* Entries returned by Crossref that have entry type 'Data' which is
  not supported by default in BibLaTeX are converted to type 'Online'
  (h/t Carl Boettiger)
* Fix for `ReadGS` when argument `check.entries` is FALSE or "warn"
  (h/t Francisco Rodriguez Sanchez)
* Family names from Scholar in all caps are handled correctly in ReadGS

## NEW FEATURES

* Functions for interacting with PubMed return language of each
  article (h/t Dale Steele)
* Added CITATION file
* Updated License to explicitly include GPL-2 and GPL-3

Changes in Version 0.8.1 (2013-03-09)
=========================================================

## BUG FIXES

* Fix for `names<-.BibEntry`
* Fix for `print.BibEntry` when entry has urldate field but no url field
* Corrections for some documentation typos
* Fix pmidrelated field when `batch.mode = FALSE` in `GetPubMedRelated`
* Fix for `LookupPubMedID` when `index` argument specified
* `open.BibEntry` now works properly
* Fix for converting thesis entries in `toBibtex.BibEntry`
* Fix for `WriteBib` with `biblatex` argument

## NEW FEATURES

* Added Vignettes including user manual and Rmd citation examples
* Added NEWS
* Added HTML output of Rmd and RHTML citation examples to doc/
Add Citations to an RMarkdown Document and Print Bibliography
========================================================

```{r setup, include = TRUE, cache = FALSE}
library(RefManageR)
bib <- ReadBib(system.file("Bib", "biblatexExamples.bib", 
                           package = "RefManageR"), check = FALSE)
BibOptions(check.entries = FALSE, style = "markdown", bib.style = "alphabetic", cite.style = 'alphabetic')
```
  This is an R Markdown document. This is an example of a citation in the text `r Citet(bib, author = "Kastenholz", .opts = list(longnamesfirst = FALSE))`. Now we cite in parentheses `r AutoCite(bib, "baez/online", before = "e.g., ")`.  You can change the default options in a setup chunk at the start of the document or at any other point using the <code>BibOptions</code> function or by specifying options as a list in the `.opts` argument to the cite functions.

These are reports `r Citet(bib, bibtype = "Report", .opts = list(hyperlink = "to.doc", super = TRUE))`.  Their hyperlinks go to their entry in the bibliography.
The link for `r AutoCite(bib, "markey")` will take you to the document in a new window; this is the default behaviour, if a link is available (see `?open.BibEntry`). The following citation has no hyperlink `r AutoCite(bib, location = "Uppsala", .opts = list(hyperlink = FALSE))`.  You can also embed plots, for example: 
  
```{r fig.width=7, fig.height=6}
plot(cars)
```

`r NoCite(bib = bib, title = "CTAN")`I've added a reference to CTAN without citing it.  Look at all my Aristotle: `r AutoCite(bib, author = "Aristotle")`.  

```{r fig.width=7, fig.height=6}
plot(cars)
```

Some papers on the arXiv are `r TextCite(bib, eprinttype = "arxiv")`.

**References**

```{r results = "asis", echo = FALSE}
PrintBibliography(bib, .opts = list(check.entries = FALSE))
```
Add Citations to an RMarkdown Document and Print Bibliography
========================================================

```{r setup, include = FALSE, cache = FALSE}
library(RefManageR)
bib <- ReadBib(system.file("Bib", "biblatexExamples.bib", 
                           package = "RefManageR"), check = FALSE)
bib2 <- ReadBib(system.file("Bib", "RJC.bib", package = "RefManageR"))[[seq_len(20)]]
BibOptions(check.entries = FALSE, style = "markdown", cite.style = "authoryear",
           bib.style = "numeric")
```
  This is an R Markdown document. This is an example of a citation in the text `r Citet(bib, "loh")`. Now we cite in parentheses `r AutoCite(bib, "baez/online", before = "e.g., ")`.  Notice the useful 'b' beside the year in the citation.  You can change the default options in a setup chunk at the start of the document or at any other point using the <code>BibOptions</code> function or by specifying options as a list in the `.opts` argument to the cite functions.  In this example we mix `"authoyear"` citation style with `"numeric"` bibliography style.

Note that I do not only have to cite by key, and may use all the features of the `SearchBib` function to index into the BibEntry object.  Here are all the entries of type `Report` in my bibliography `r Citet(bib, bibtype = "Report", .opts = list(hyperlink = "to.doc"))`.  The hyperlinks will take you to their entry in the bibliography.  The link for `r TextCite(bib, "markey")` will open the document in a new window; this is the default behaviour, if a link is available (see `?open.BibEntry`). The following citation has no hyperlink `r AutoCite(bib, location = "Uppsala", .opts = list(hyperlink = FALSE))`.  You can also embed plots, to make the page longer: 
  
```{r fig.width=7, fig.height=6}
plot(cars)
```
`r NoCite(bib = bib, title = "CTAN")`I've added a reference to CTAN without citing it using the `NoCite` function.  Now I'm adding a reference from another bibliography (a second `BibEntry` object) `r AutoCite(bib2, title = "binary longitudinal data")`.  Look at all my Aristotle: `r TextCite(bib, author = "Aristotle")`.  

```{r fig.width=7, fig.height=6}
plot(cars)
```

Some papers on the arXiv are `r TextCite(bib, eprinttype = "arxiv")`.

**References**

```{r results = "asis", echo = FALSE}
PrintBibliography(bib, .opts = list(check.entries = FALSE, sorting = "ynt"))
```

**More References**

```{r results = "asis", echo = FALSE}
PrintBibliography(bib2)
```Add Citations to an RMarkdown Document and Print Bibliography
========================================================

```{r setup, include = TRUE, cache = FALSE}
library(RefManageR)
bib <- ReadBib(system.file("Bib", "biblatexExamples.bib", 
                           package = "RefManageR"), check = FALSE)
BibOptions(check.entries = FALSE, style = "markdown", bib.style = "numeric", cite.style = "numeric")
```
  This is an R Markdown document. This is an example of a citation in the text `r Citet(bib, 12, .opts = list(longnamesfirst = FALSE))`. Now we cite in parentheses `r AutoCite(bib, "baez/online", before = "e.g., ")`.  You can change the default options in a setup chunk at the start of the document or at any other point using the <code>BibOptions</code> function or by specifying options as a list in the `.opts` argument to the cite functions.

See what happens when `r AutoCite(bib, author = "kant")` use the shorthand field?

These are reports `r Citet(bib, bibtype = "Report", .opts = list(hyperlink = "to.doc", super = TRUE))`.  Their hyperlinks go to their entry in the bibliography.
The link for `r AutoCite(bib, "markey")` will take you to the document in a new window; this is the default behaviour, if a link is available (see `?open.BibEntry`). The following citation has no hyperlink `r AutoCite(bib, location = "Uppsala", .opts = list(hyperlink = FALSE))`.  You can also embed plots, for example: 
  
```{r fig.width=7, fig.height=6}
plot(cars)
```
`r NoCite(bib = bib, title = "CTAN")`I've added a reference to CTAN without citing it.  Look at all my Aristotle: `r AutoCite(bib, author = "Aristotle")`.  

```{r fig.width=7, fig.height=6}
plot(cars)
```

Some papers on the arXiv are `r TextCite(bib, eprinttype = "arxiv")`.

**References**

```{r results = "asis", echo = FALSE}
PrintBibliography(bib, .opts = list(check.entries = FALSE))
```
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{RMarkdown Citations - Alphabetic Style}
-->
Add Citations to an RMarkdown Document and Print Bibliography
========================================================

```{r setup, include = TRUE, cache = FALSE}
library(RefManageR)
bib <- ReadBib(system.file("Bib", "biblatexExamples.bib", 
                           package = "RefManageR"), check = FALSE)
BibOptions(check.entries = FALSE, style = "markdown", bib.style = "alphabetic", cite.style = 'alphabetic')
```
  This is an R Markdown document. This is an example of a citation in the text `r Citet(bib, 12, .opts = list(longnamesfirst = FALSE))`. Now we cite in parentheses `r AutoCite(bib, "baez/online", before = "e.g., ")`.  You can change the default options in a setup chunk at the start of the document or at any other point using the <code>BibOptions</code> function or by specifying options as a list in the `.opts` argument to the cite functions.

These are reports `r Citet(bib, bibtype = "Report", .opts = list(hyperlink = "to.doc", super = TRUE))`.  Their hyperlinks go to their entry in the bibliography.
The link for `r AutoCite(bib, "markey")` will take you to the document in a new window; this is the default behaviour, if a link is available (see `?open.BibEntry`). The following citation has no hyperlink `r AutoCite(bib, location = "Uppsala", .opts = list(hyperlink = FALSE))`.  You can also embed plots, for example: 
  
```{r fig.width=7, fig.height=6}
plot(cars)
```
`r NoCite(bib = bib, title = "CTAN")`I've added a reference to CTAN without citing it.  Look at all my Aristotle: `r AutoCite(bib, author = "Aristotle")`.  

```{r fig.width=7, fig.height=6}
plot(cars)
```

Some papers on the arXiv are `r TextCite(bib, eprinttype = "arxiv")`.

**References**

```{r results = "asis", echo = FALSE}
PrintBibliography(bib, .opts = list(check.entries = FALSE))
```
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{RMarkdown Citations Example}
-->
Add Citations to an RMarkdown Document and Print Bibliography
========================================================

```{r unload, include = FALSE}
## Needed to clear internal package citation list from previous vignette builds
unloadNamespace("RefManageR")
```

```{r setup, include = FALSE, cache = FALSE}
library(RefManageR)
bib <- ReadBib(system.file("Bib", "biblatexExamples.bib", 
                           package = "RefManageR"), check = FALSE)
bib2 <- ReadBib(system.file("Bib", "RJC.bib", package = "RefManageR"))[[seq_len(20)]]
BibOptions(check.entries = FALSE, style = "markdown", cite.style = "authoryear",
           bib.style = "numeric")
```
  This is an R Markdown document. This is an example of a citation in the text `r Citet(bib, "loh")`. Now we cite in parentheses `r AutoCite(bib, "baez/online", before = "e.g., ")`.  Notice the useful 'b' beside the year in the citation.  You can change the default options in a setup chunk at the start of the document or at any other point using the <code>BibOptions</code> function or by specifying options as a list in the `.opts` argument to the cite functions.  In this example we mix `"authoyear"` citation style with `"numeric"` bibliography style.

Note that I do not only have to cite by key, and may use all the features of the `SearchBib` function to index into the BibEntry object.  Here are all the entries of type `Report` in my bibliography `r Citet(bib, bibtype = "Report", .opts = list(hyperlink = "to.doc"))`.  The hyperlinks will take you to their entry in the bibliography.  The link for `r TextCite(bib, "markey")` will open the document in a new window; this is the default behaviour, if a link is available (see `?open.BibEntry`). The following citation has no hyperlink `r AutoCite(bib, location = "Uppsala", .opts = list(hyperlink = FALSE))`.  You can also embed plots, to make the page longer: 
  
```{r fig.width=7, fig.height=6}
plot(cars)
```
`r NoCite(bib = bib, title = "CTAN")`I've added a reference to CTAN without citing it using the `NoCite` function.  Now I'm adding a reference from another bibliography (a second `BibEntry` object) `r AutoCite(bib2, title = "binary longitudinal data")`.  Look at all my Aristotle: `r TextCite(bib, author = "Aristotle")`.  

```{r fig.width=7, fig.height=6}
plot(cars)
```

Some papers on the arXiv are `r TextCite(bib, eprinttype = "arxiv")`.

**References**

```{r results = "asis", echo = FALSE}
PrintBibliography(bib, .opts = list(check.entries = FALSE, sorting = "ynt"))
```

**More References**

```{r results = "asis", echo = FALSE}
PrintBibliography(bib2)
```% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadBib.R
\name{ReadBib}
\alias{ReadBib}
\title{BibLaTeX/BibTeX .bib file parser}
\usage{
ReadBib(
  file,
  .Encoding = "UTF-8",
  header = if (length(preamble)) paste(preamble, sep = "\\n") else "",
  footer = "",
  check = BibOptions()$check.entries
)
}
\arguments{
\item{file}{string; bib file to parse.}

\item{.Encoding}{encoding}

\item{header}{header of the citation list. By default this is made from the Preamble
entries found in the bib file.}

\item{footer}{footer of the citation list.}

\item{check}{\dQuote{error}, \dQuote{warn}, or logical \code{FALSE}.  What action
should be taken if an entry is missing required fields?  \code{FALSE} means no
checking is done, \dQuote{warn} means entry is added with an error.
\dQuote{error} means the entry will not be added.  See \code{\link{BibOptions}}.}
}
\description{
Parser for bibliography databases in the bib format containing either BibLaTeX or
BibTeX entries.
}
\note{
Date fields are parsed using the locale specified by
\code{Sys.getlocale("LC_TIME")}.  To read a bib file with character \sQuote{month}
fields in a language other than the current locale, \code{Sys.setlocale} should be
used to change \sQuote{LC_TIME}` to match the bib file before calling \code{ReadBib}.

Keys will be made unique by calling \code{\link[base]{make.unique}} with
\code{sep = ":"}.
}
\examples{
if (requireNamespace("bibtex")) {
    file.name <- system.file("Bib", "RJC.bib", package="RefManageR")
    bib <- ReadBib(file.name)
}
}
\seealso{
\code{\link[bibtex]{read.bib}} in package \code{bibtex}
}
\author{
McLean, M. W., based on code in \code{bibtex} package by Francois, R.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NamesAssign.R, R/names.R
\name{names<-.BibEntry}
\alias{names<-.BibEntry}
\alias{names.BibEntry}
\title{Names (keys) of a BibEntry object}
\usage{
\method{names}{BibEntry}(x) <- value

\method{names}{BibEntry}(x)
}
\arguments{
\item{x}{an object of class BibEntry}

\item{value}{character vector of new key values to replace into \code{x}}
}
\value{
\code{names<-} the updated BibEntry object.

\code{names} - character vector of the keys of the BibEntry object.
}
\description{
Functions to get and set the keys of an object of class BibEntry
}
\examples{
if (requireNamespace("bibtex")) {
    bib <- ReadBib(system.file("Bib", "test.bib", package = "RefManageR"))
    names(bib)
    names(bib)[1] <- 'newkey'
}
}
\author{
McLean, M. W. \email{mathew.w.mclean@gmail.com}
}
\keyword{attribute}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BibEntryListExtract.R
\name{[[.BibEntry}
\alias{[[.BibEntry}
\title{Extract entries from a BibEntry object by index}
\usage{
\method{[[}{BibEntry}(x, i, drop = FALSE)
}
\arguments{
\item{x}{a BibEntry object}

\item{i}{numeric indices of entries to extract, or a character vector of keys corresponding to the entries to be
extracted.}

\item{drop}{logical, should attributes besides class be dropped from result?}
}
\value{
an object of class BibEntry.
}
\description{
Operator for extracting BibEntry objects by index.
}
\note{
This method is different than the usual operator \code{[[} for lists in that a vector of indices may be specified.

This method behaves differently than the \code{[} operator for BibEntry objects in that it does not expand
crossreferences when returning, so that a parent entry or xdata entry will be dropped if it is not also indexed
when indexing the child entry.

This method is not affected by the value of \code{BibOptions()$return.ind}.
}
\examples{
if (requireNamespace("bibtex")) {
    file.name <- system.file("Bib", "biblatexExamples.bib", package="RefManageR")
    bib <- suppressMessages(ReadBib(file.name))
    bib[[20:21]]
    bib[c("hyman", "loh")]

    ## Note this is FALSE because [[ does not inherit from the dropped parent entry while [ does.
    identical(bib[1], bib[[1]])
}
}
\seealso{
Other operators: 
\code{\link{$.BibEntry}()},
\code{\link{$<-.BibEntry}()},
\code{\link{+.BibEntry}()},
\code{\link{[.BibEntry}()},
\code{\link{[<-.BibEntry}()},
\code{\link{[[<-.BibEntry}()},
\code{\link{c.BibEntry}()}
}
\concept{operators}
\keyword{database}
\keyword{list}
\keyword{manip}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadCrossRef.R
\name{ReadCrossRef}
\alias{ReadCrossRef}
\title{Search CrossRef for citations.}
\usage{
ReadCrossRef(
  query = "",
  filter = list(),
  limit = 5,
  offset = 0,
  sort = "relevance",
  year = NULL,
  min.relevance = 2,
  temp.file = tempfile(fileext = ".bib"),
  delete.file = TRUE,
  verbose = FALSE,
  use.old.api = FALSE
)
}
\arguments{
\item{query}{string; search term}

\item{filter}{named list of possible filters; see \code{Details}
and \code{References}; ignored if \code{use.old.api = TRUE}}

\item{limit}{numeric; maximum number of entries to return}

\item{offset}{numeric; CrossRef will not return the first
\code{offset} results (default 0); ignored if \code{use.old.api
= TRUE}}

\item{sort}{string; how specifying how the results from CrossRef
should be sorted.  Possible values when \code{use.old.api =
FALSE} are \code{"score"} (default; same as
\code{"relevance"}), \code{"updated"}, \code{"deposited"},
\code{"indexed"}, or \code{"published"}; see the references}

\item{year}{numeric; if specified, only results from this year will
be returned.}

\item{min.relevance}{numeric; only results with a CrossRef-assigned
relevance score at least this high will be returned.}

\item{temp.file}{string; file name to use for storing Bibtex
information returned by CrossRef.}

\item{delete.file}{boolean; should the bib file be deleted on exit?}

\item{verbose}{boolean; if \code{TRUE}, additional messages are
output regarding the results of the query.}

\item{use.old.api}{boolean; should the older CrossRef API be used
for the search? NO LONGER SUPPORTED, all queries need to use
the new API.}
}
\value{
An object of class \code{BibEntry}.
}
\description{
Provides an interface to the CrossRef API, searching for citations
given a string query.  Results are written to a bib file, read back
into \code{R} using \code{\link{WriteBib}}, and returned as a
BibEntry object.
}
\details{
When \code{use.old.api = TRUE}, the query HTTP request only returns DOIs,
which are then used to make HTTP requests for the corresponding BibTeX entries from
CrossRef; when \code{use.old.api = FALSE}, the query HTTP request is parsed to create
the \code{BibEntry} object (i.e. there are less HTTP requests when using the new API).

CrossRef assigns a score between 0 and 100 based on how relevant a
reference seems to be to your query.  The \emph{old} API
documentation warns that while false negatives are unlikely, the
search can be prone to false positives.  Hence, setting
\code{min.revelance} to a high value may be necessary if
\code{use.old.api = TRUE}. In some instances with the old API, no
score is returned, if this happens, the entries are added with a
message indicating that no score was available.

Possible values for the \emph{names} in \code{filter} are \code{"has-funder"},
\code{"funder"}, \code{"prefix"}, \code{"member"}, \code{"from-index-date"},
\code{"until-index-date"},
\code{"from-deposit-date"}, \code{"until-deposit-date"}, \code{"from-update-date"},
\code{"until-update-date"}, \code{"from-created-date"}, \code{"until-created-date"},
\code{"from-pub-date"}, \code{"until-pub-date"}, \code{"has-license"}, \code{"license.url"},
\code{"license.version"}, \code{"license.delay"}, \code{"has-full-text"},
\code{"full-text.version"}, \code{"full-text.type"}, \code{"public-references"},
\code{"has-references"}, \code{"has-archive"}, \code{"archive"}, \code{"has-orcid"},
\code{"orcid"}, \code{"issn"}, \code{"type"}, \code{"directory"}, \code{"doi"},
\code{"updates"}, \code{"is-update"}, \code{"has-update-policy"}, \code{"container-title"},
\code{"publisher-name"}, \code{"category-name"}, \code{"type-name"}, \code{"award.number"},
\code{"award.funder"}, \code{"assertion-group"}, \code{"assertion"}, \code{"affiliation"},
\code{"has-affiliation"}, \code{"alternative-id"}, and \code{"article-number"}.
See the first reference for a description of their meanings.
}
\note{
The entries returned by Crossref are frequently missing
    fields required by BibTeX, if you want the entries to be
    returned anyway, set \code{BibOptions()$check.entries} to
    \code{FALSE} or \code{"warn"}

Fields \code{"score"} (the relevancy score) and \code{"license"} will be
returned when \code{use.old.api = FALSE}.
}
\examples{
if (interactive() && !httr::http_error("https://search.crossref.org/")){
  BibOptions(check.entries = FALSE)
  ## 3 results from the American Statistical Association involving "regression"
  ReadCrossRef("regression", filter = list(prefix="10.1198"), limit = 3)

  ## Some JRSS-B papers published in 2010 or later, note the quotes for filter
  ##   names with hypens
  ReadCrossRef(filter = list(issn = "1467-9868", "from-pub-date" = 2010),
               limit = 2, min.relevance = 0)

  ## Articles published by Institute of Mathematical Statistics
  ReadCrossRef(filter = list(prefix = "10.1214"), limit = 5, min.relevance = 0)

  ## old API
  ReadCrossRef(query = 'rj carroll measurement error', limit = 2, sort = "relevance",
    min.relevance = 80, use.old.api = TRUE)

  ReadCrossRef(query = 'carroll journal of the american statistical association',
    year = 2012, limit = 2, use.old.api = TRUE)
}
}
\references{
Newer API: \url{https://github.com/CrossRef/rest-api-doc/blob/master/rest_api.md},
Older API: \url{https://search.crossref.org/help/api}
}
\seealso{
\code{\link{ReadZotero}}, \code{\link{BibEntry}},
package \code{rcrossref} for larger queries and deep paging

Other pubmed: 
\code{\link{GetPubMedByID}()},
\code{\link{GetPubMedRelated}()},
\code{\link{LookupPubMedID}()},
\code{\link{ReadPubMed}()}
}
\concept{pubmed}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/08asBibEntry.R
\name{as.BibEntry}
\alias{as.BibEntry}
\alias{is.BibEntry}
\title{Coerce to a BibEntry object}
\usage{
as.BibEntry(x)

is.BibEntry(x)
}
\arguments{
\item{x}{any \code{R} object.}
}
\value{
\code{as.BibEntry} - if successful, an object of class BibEntry.

\code{is.BibEntry} - logical; \code{TRUE} if \code{x} is a BibEntry
object.
}
\description{
Functions to check if an object is a BibEntry, or coerce it if possible.
}
\details{
\code{as.BibEntry} is able to coerce suitably formatted character
vectors, \code{\link{bibentry}} objects, lists,
and data.frames to BibEntry objects.  See the examples.
}
\note{
Each entry to be coerced should have a bibtype, key, and all required
fields for the specified bibtype.
}
\examples{
if (requireNamespace("bibtex")) {
    file.name <- system.file("Bib", "biblatexExamples.bib", package="RefManageR")
    bib <- suppressMessages(ReadBib(file.name))[[20:21]]
    identical(as.BibEntry(unlist(bib)), bib)  ## see also RelistBibEntry

    identical(as.BibEntry(unclass(bib)), bib)

    identical(as.BibEntry(as.data.frame(bib)), bib)
 }

bib <- c(bibtype = "article", key = "mclean2014", title = "My New Article",
  author = "Mathew W. McLean", journaltitle = "The Journal", date = "2014-01")
as.BibEntry(bib)

bib <- bibentry(bibtype = "article", key = "mclean2014", title = "My New Article",
journal = "The Journal", year = 2014, author = "Mathew W. McLean")
print(bib, .bibstyle = "JSS")
as.BibEntry(bib)

bib <- list(c(bibtype = "article", key = "mclean2014a", title = "My New Article",
  author = "Mathew W. McLean", journaltitle = "The Journal", date = "2014-01"),
  c(bibtype = "article", key = "mclean2014b", title = "Newer Article",
  author = "Mathew W. McLean", journaltitle = "The Journal", date = "2014-02"))
as.BibEntry(bib)
}
\seealso{
\code{\link{BibEntry}}
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadPubMed.R
\name{ReadPubMed}
\alias{ReadPubMed}
\title{Search NCBI's E-Utilities for citation information}
\usage{
ReadPubMed(query, database = "PubMed", ...)
}
\arguments{
\item{query}{string; search term.}

\item{database}{string; the Entrez database to search.}

\item{...}{additional parameters to use for the search.
See the \emph{Details}.}
}
\value{
an object of class BibEntry.
}
\description{
This function takes a query and searches an Entrez database for
references using NCBI's E-Utilities, returning the results in a BibEntry
object.
}
\details{
Optional additional parameters to pass to the server include
\itemize{
\item \code{retstart} - index of the first retrieved ID that should be
included in the results.
\item \code{retmax} - maximum number of IDs the server will
return (default 20). 
\item \code{field} - limits the query to search only the specified
field (e.g. \dQuote{title}).
\item \code{datetype} - type of date to use when limiting search by
dates. E.g. \dQuote{mdat}
for modification date or \dQuote{pdat} for publication date.
\item \code{reldate} - integer; only items that have (\code{datetype})
date values within \code{reldate} \emph{days}
of the current date will be returned.
\item \code{mindate}, \code{maxdate} - date ranges to restrict search
results.  Possible formats are
\dQuote{YYYY}, \dQuote{YYYY/MM}, and \dQuote{YYYY/MM/DD}.
}
}
\note{
The returned entries will have type either \sQuote{Article} or
\sQuote{Misc} depending on whether journal information was retrieved.
See the Entrez documentation listed in the \emph{References}.

The language of the entry will be returned in the field \dQuote{language}
and the abstract will be returned in the field \dQuote{abstract}, if they
are available.
}
\examples{
if (interactive() && !httr::http_error("https://eutils.ncbi.nlm.nih.gov/"))
  ReadPubMed(query = "raymond carroll measurement error", retmax = 5, mindate = 1990)
}
\references{
\url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch}
}
\seealso{
Other pubmed: 
\code{\link{GetPubMedByID}()},
\code{\link{GetPubMedRelated}()},
\code{\link{LookupPubMedID}()},
\code{\link{ReadCrossRef}()}
}
\concept{pubmed}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadPDFs.R
\name{ReadPDFs}
\alias{ReadPDFs}
\title{Create bibliographic information from PDF Metadata.}
\usage{
ReadPDFs(
  path,
  .enc = "UTF-8",
  recursive = TRUE,
  use.crossref = TRUE,
  use.metadata = TRUE,
  progress = FALSE
)
}
\arguments{
\item{path}{character; path to directory containing pdfs or filename of
one pdf. \code{normalizePath} is used on the specified path}

\item{.enc}{character; text encoding to use for reading pdf and creating
BibEntry object. Available encodings for Poppler can be found using
\code{system("pdfinfo -listenc")}.  The encoding must also be listed
in \code{iconvlist()}.}

\item{recursive}{logical; same as \code{\link{list.files}}.  Should pdfs
in subdirectories of path be used?}

\item{use.crossref}{logical; should an attempt be made to download bibliographic
information from CrossRef if
any Document Object Identifiers (DOIs) are found? This is only supported if the Suggeseted package \code{bibtex} is found.}

\item{use.metadata}{logical; should the PDF metadata also be used to help
create entries?}

\item{progress}{logical; should progress bar be generated when fetching from
CrossRef?}
}
\value{
An object of class BibEntry.
}
\description{
This function creates bibliographic information by reading the Metadata and
text of PDFs stored in a user specified directory using Poppler
(\url{https://poppler.freedesktop.org/}).  IF requested, the function
first searches for DOIs and downloads \code{BibTeX} entries from
\code{\link{ReadCrossRef}} if DOIs are found.  If this is not requested or
a DOI is not found for an entry, an attempt is made to build a BibTeX
entry from the metadata and text.
}
\details{
This function requires that the \code{pdfinfo} utility from Poppler PDF
\url{https://poppler.freedesktop.org/} be installed.

This function will create only \code{Article} or \code{Misc} \code{BibTeX} entries.

The absolute path to each file will be stored in the bib entry in a field
called \sQuote{file}, which is recognized by \code{BibLaTeX} (though not printed by
any standard style) and can be used by the
\code{\link{open.BibEntry}} function to open the PDF in the default viewer.

If the keywords metadata field is available, it will be added to the bib entry
in a field \sQuote{keywords}, which is recognized by \code{BibLaTeX}.
}
\examples{
\dontrun{
path <- system.file("doc", package = "RefManageR")
ReadPDFs(path)
}
}
\references{
\url{https://poppler.freedesktop.org/}
}
\seealso{
\code{\link{ReadCrossRef}}, \code{\link{BibEntry}},
\code{\link{open.BibEntry}}
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06BibEntry.R
\name{BibEntry}
\alias{BibEntry}
\title{Enhanced Bibliographic Entries}
\usage{
BibEntry(
  bibtype,
  textVersion = NULL,
  header = NULL,
  footer = NULL,
  key = NULL,
  ...,
  other = list(),
  mheader = NULL,
  mfooter = NULL
)
}
\arguments{
\item{bibtype}{a character string with a BibTeX entry type. See Entry Types
for details.}

\item{textVersion}{a character string with a text representation of the
reference to optionally be employed for printing.}

\item{header}{a character string with optional header text.}

\item{footer}{a character string with optional footer text.}

\item{key}{a character string giving the citation key for the entry.}

\item{...}{arguments of the form \code{tag = value} giving the fields
of the entry, with \code{tag} and \code{value} the name and value of the
field, respectively. Arguments with empty values are dropped. See
  Entry Fields for details.}

\item{other}{list; additional way to specify fields and their values}

\item{mheader}{string; optional \dQuote{outer} header text}

\item{mfooter}{string; optional \dQuote{outer} footer text}
}
\value{
an object of class BibEntry
}
\description{
Provides a new class \code{BibEntry} which builds on \code{\link{bibentry}}
to provide enhanced functionality for representing, manipulating, importing,
etc. bibliographic information in BibTeX or BibLaTeX style.
}
\details{
The BibEntry objects created by BibEntry can represent an
arbitrary positive number of references, as with \code{bibentry}, but
many additional methods are defined for building and manipulating a database
  of references.
}
\note{
Date fields are parsed using the locale specified by
\code{Sys.getlocale("LC_TIME")} (relevant when specifying a character
\sQuote{month} field, instead of the recommended integer format)

Name list fields (author, editor, etc.) should be specified as they would be for
BibTeX/BibLaTeX; e.g. \code{author = "Doe, Jane and Smith, Bob A."}.
}
\section{Entry Types}{

  bibentry creates "bibentry" objects, which are modeled after BibLaTeX and
BibTeX entries. The entry should
  be a valid BibLaTeX or BibTeX entry type.  For a list of valid BibTeX entry
types, see \code{\link{bibentry}}.  BibLaTeX supports all entry types from
BibTeX for backwards compatibility. BibLaTeX defines following entry types
  '  \itemize{
   \item \emph{article} - An article in a journal, magazine, newspaper, or other periodical which forms a
  self-contained unit with its own title.  Required fields: author, title, journal/journaltitle, year/date.
   \item \emph{book} - A single-volume book with one or more authors where the authors share credit for
  the work as a whole.  Required fields: author, title, year/date. (Also covers BibTeX @inbook).
   \item \emph{mvbook} - A multi-volume \emph{book}. For backwards compatibility, multi-volume books are also
  supported by the entry type @book.  Required fields: author, title, year/date.
   \item \emph{inbook} - A part of a book which forms a self-contained unit with its own title. Note that the
  profile of this entry type is different from standard BibTeX.  Required fields: author, title,
  booktitle, year/date
   \item \emph{bookinbook} This type is similar to \emph{inbook} but intended for works originally published as a
  stand-alone book.
   \item \emph{suppbook} - Supplemental material in a \emph{book}. This type is closely related to the \emph{inbook} entry
  type.
   \item \emph{booklet} - A book-like work without a formal publisher or sponsoring institution.  Required fields:
   author/editor, title, year/date.
   \item \emph{collection} - A single-volume collection with multiple, self-contained contributions by distinct
   authors which have their own title.  Required fields: editor, title, year/date.
   \item \emph{mvcollection} - A multi-volume \emph{collection}.  Also supported by \emph{collection}.
   Required fields: editor, title, year/date.
   \item \emph{incollection} - A contribution to a collection which forms a self-contained unit with a distinct
   author and title.  Required fields: author, editor, title, booktitle, year/date.
   \item \emph{suppcollection} - Supplemental material in a \emph{collection}.
   \item \emph{manual} - Technical or other documentation, not necessarily in printed form. Required fields:
     author/editor, title, year/date
   \item \emph{misc} - A fallback type for entries which do not fit into any other category.  Required fields:
   author/editor, title, year/date
   \item \emph{online} - An online resource.  Required fields: author, title, number, year/date.
   \item \emph{patent} - A patent or patent request.  Required fields: author, title, number, year/date.
   \item \emph{periodical} - A complete issue of a periodical, such as a special issue of a journal.  Required
   fields: editor, title, year/date
   \item \emph{suppperiodical} - Supplemental material in a \emph{periodical}.
   \item \emph{proceedings} - A single-volume conference proceedings.  Required fields: editor, title, year/date.
   \item \emph{mvproceedings} - A multi-volume @proceedings entry. Required fields: editor, title, year/date.
   \item \emph{inproceedings} - An article in a conference proceedings.  Required fields: author, editor, title,
   booktitle, year/date.
   \item \emph{reference} - A single-volume work of reference such as an encyclopedia or a dictionary.  Alias
         for \emph{collection} in standard styles.
   \item \emph{mvreference} - A multi-volume \emph{reference} entry.
   \item \emph{inreference} - An article in a work of reference.  Alias for \emph{incollection} in most styles.
   \item \emph{report} - A technical report, research report, or white paper published by a university or some
   other institution.  Required fields: author, title, type, institution, year/date.
   \item \emph{set} - An entry set. This entry type is special, see BibLaTeX manual.
   \item \emph{thesis} - A thesis written for an educational institution to satisfy the requirements for a degree.
   Use the type field to specify the type of thesis.  Required fields: author, title, type, institution,
   year/date.
   \item \emph{unpublished} A work with an author and a title which has not been formally published, such as a
   manuscript or the script of a talk.  Required fields: author, title, year/date.
   \item \emph{xdata} - This entry type is special. \emph{xdata} entries hold data which may be inherited by other
   entries. (Biber only.)
   \item \emph{custom[a-f]} Custom types (up to five) for special bibliography styles. Not used by the standard
    styles.
 }
}

\examples{
BibEntry(bibtype = "Article", key = "mclean2014", title = "An Article Title",
  author = "McLean, Mathew W. and Wand, Matt P.", journaltitle = "The Journal Title",
  date = "2014-02-06", pubstate = "forthcoming")
bib <- BibEntry(bibtype = "XData", key = "arxiv_data", eprinttype = "arxiv",
eprintclass = "stat.ME", year = 2013, urldate = "2014-02-01", pubstate = "submitted")
bib <- c(bib, BibEntry(bibtype = "Misc", key = "mclean2014b",
  title = "Something On the {arXiv}", author = "Mathew W. McLean", eprint = "1312.9999",
  xdata = "arxiv_data", url = "https://arxiv.org/abs/1310.5811"))
bib
toBiblatex(bib)
}
\references{
BibLaTeX manual \url{https://mirror.pregi.net/tex-archive/macros/latex/contrib/biblatex/doc/biblatex.pdf}
}
\seealso{
\code{\link{bibentry}}
}
\author{
McLean, M. W. \email{mathew.w.mclean@gmail.com}
}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetBibEntryWithDOI.R
\name{GetBibEntryWithDOI}
\alias{GetBibEntryWithDOI}
\title{Lookup a Bibtex entry using a Digital Object Identifier}
\usage{
GetBibEntryWithDOI(
  doi,
  temp.file = tempfile(fileext = ".bib"),
  delete.file = TRUE
)
}
\arguments{
\item{doi}{character vector; DOIs to use to retrieve bibliographic information.}

\item{temp.file}{string; a file to write the Bibtex data returned by the
DOI System to.}

\item{delete.file}{logical; should \code{temp.file} be deleted when the
function exits?}
}
\value{
an object of class BibEntry.
}
\description{
Uses the DOI System API to look up bibliography information given a set of DOIs.
}
\details{
The bibliographic information returned by the search of the \url{https://doi.org/}
API is temporarily
written to a file and then read back into \code{R} and return as a
\code{BibEntry} object.
}
\examples{
if (interactive() && !httr::http_error("https://doi.org/"))
  GetBibEntryWithDOI(c("10.1016/j.iheduc.2003.11.004", "10.3998/3336451.0004.203"))
}
\references{
\url{https://www.doi.org/tools.html}
}
\seealso{
\code{\link{ReadCrossRef}}, \code{\link{BibEntry}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BibEntryDollarExtract.R
\name{$.BibEntry}
\alias{$.BibEntry}
\title{Extract fields from a BibEntry object}
\usage{
\method{$}{BibEntry}(x, name)
}
\arguments{
\item{x}{an object of class BibEntry}

\item{name}{the field to extract}
}
\value{
a named list of values for the field specified by name for each entry; \code{NULL} if the field is not present for
a particular entry.  The names attribute of the returned list contains the entry keys (potentially back-quoted).
}
\description{
used to extract a single field from each entry in a BibEntry object
}
\note{
\code{name} may be \dQuote{bibtype} to extract entry types or \dQuote{key} to extract keys.
}
\examples{
if (requireNamespace("bibtex")) {
    file.name <- system.file("Bib", "biblatexExamples.bib", package="RefManageR")
    bib <- suppressMessages(ReadBib(file.name))
    bib[[50:55]]$author
    bib[[seq_len(5)]]$bibtype
 }
}
\seealso{
Other operators: 
\code{\link{$<-.BibEntry}()},
\code{\link{+.BibEntry}()},
\code{\link{[.BibEntry}()},
\code{\link{[<-.BibEntry}()},
\code{\link{[[.BibEntry}()},
\code{\link{[[<-.BibEntry}()},
\code{\link{c.BibEntry}()}
}
\concept{operators}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BibEntryReplaceBibEntry.R
\name{[[<-.BibEntry}
\alias{[[<-.BibEntry}
\title{Assign a BibEntry entry to another BibEntry object}
\usage{
\method{[[}{BibEntry}(x, i) <- value
}
\arguments{
\item{x}{- a BibEntry object}

\item{i}{- a numeric index or a string entry key}

\item{value}{- a single entry BibEntry object or an object that can be 
coerced to BibEntry using \code{\link{as.BibEntry}}}
}
\value{
an object of class BibEntry
}
\description{
Replace one entry in a BibEntry object with another
}
\details{
This function will replace the specified entry in \code{x} with the entry given
by \code{value}.  To replace multiple entries see \code{\link{[<-.BibEntry}}.
}
\seealso{
Other operators: 
\code{\link{$.BibEntry}()},
\code{\link{$<-.BibEntry}()},
\code{\link{+.BibEntry}()},
\code{\link{[.BibEntry}()},
\code{\link{[<-.BibEntry}()},
\code{\link{[[.BibEntry}()},
\code{\link{c.BibEntry}()}
}
\concept{operators}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadZotero.R
\name{ReadZotero}
\alias{ReadZotero}
\title{Get Bibliography Information From a Zotero Library.}
\usage{
ReadZotero(
  user,
  group,
  .params,
  temp.file = tempfile(fileext = ".bib"),
  delete.file = TRUE
)
}
\arguments{
\item{user}{Zotero userID for use in calls to the Zotero API.  This is not
the same as your Zotero username.  The userID for accessing user-owned
libraries can be found at \verb{https://www.zotero.org/settings/keys}
after logging in.}

\item{group}{Zotero groupID for use in calls to the Zotero API.  Only one
of \code{user} and \code{group} should be specified; \code{group} will be
ignored if both are specified.}

\item{.params}{A \emph{named} list of parameters to use in requests to the
Zotero API with possible values
 \itemize{
   \item q - Search string to use to search the library
   \item qmode - Search mode. Default is "titleCreatorYear".  Use "everything"
         to include full-text content in search.
   \item key - API key.  This must be specified to access non-public libraries.
   \item collection - name of a specific collection within the library to search
   \item itemType - type of entry to search for; e.g., "book" or "journalArticle"
   \item tag - name of tag to search for in library
   \item limit - maximum number of entries to return
   \item start - index of first entry to return
 }}

\item{temp.file}{character; file name where the BibTeX data returned by
Zotero will be temporarily written.}

\item{delete.file}{boolean; should \code{temp.file} be removed on exit?}
}
\value{
An object of class BibEntry
}
\description{
Get Bibliography Information From a Zotero Library.
}
\examples{
\dontrun{
## first two entries in library with bayesian in title
ReadZotero(user = "1648676", .params = list(q = "bayesian",
  key = "7lhgvcwVq60CDi7E68FyE3br", limit=2))

## Search specific collection
## collection key can be found by reading uri when collection is selected in Zotero
ReadZotero(user = "1648676", .params=list(q = "yu", key = "7lhgvcwVq60CDi7E68FyE3br",
  collection = "3STEQRNU"))

## Search by tag
## Notice the issue with how Zotero uses a TechReport entry for arXiv manuscripts
## This is one instance where the added fields of BibLaTeX are useful
ReadZotero(user = "1648676", .params=list(key = "7lhgvcwVq60CDi7E68FyE3br",
  tag = "Statistics - Machine Learning"))

## To read these in you must set check.entries to FALSE or "warn"
old.opts <- BibOptions(check.entries = FALSE)
length(ReadZotero(user = "1648676", .params = list(key = "7lhgvcwVq60CDi7E68FyE3br",
  tag = "Statistics - Machine Learning")))

## Example using groups
ReadZotero(group = "13495", .params = list(q = "Schmidhuber",
  collection = "QU23T27Q"))
BibOptions(old.opts)
}
}
\references{
\verb{https://www.zotero.org/support/dev/server_api/v2/read_requests}
}
\seealso{
\code{\link{BibEntry}}
}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadPubMed.R
\name{GetPubMedByID}
\alias{GetPubMedByID}
\title{Retrieve citation information from NCBI's Entrez for a set of PubMed IDs}
\usage{
GetPubMedByID(id, db = "pubmed", ...)
}
\arguments{
\item{id}{character vector; PubMed ID's for searching NCBI's Entrez.}

\item{db}{string; Entrez database to search.}

\item{...}{additional parameters to use for the search.
See the Entrez documentation listed in the \emph{References}.}
}
\value{
a BibEntry object.
}
\description{
Uses NCBI's E-Utilities to retrieve bibliographic information given a
vector of PubMed ID's and returns the results as a BibEntry object.
}
\note{
Returned entries will have \code{bibtype} \dQuote{Article} or \dQuote{Book},
unless a collection title is present -- in which case the \code{bibtype} will be
\dQuote{InBook} -- or there is no journal information returned for an article -- in
which case the \code{bibtype} will be \dQuote{Misc}.
}
\examples{
if (interactive() && !httr::http_error("https://eutils.ncbi.nlm.nih.gov/"))
  GetPubMedByID(c("11209037", "21245076"))
}
\references{
\url{https://www.ncbi.nlm.nih.gov/books/NBK25500/}
}
\seealso{
Other pubmed: 
\code{\link{GetPubMedRelated}()},
\code{\link{LookupPubMedID}()},
\code{\link{ReadCrossRef}()},
\code{\link{ReadPubMed}()}
}
\concept{pubmed}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BibEntryReplaceOp.R
\name{[<-.BibEntry}
\alias{[<-.BibEntry}
\title{Update Different Fields of Multiple Entries of a BibEntry Object}
\usage{
\method{[}{BibEntry}(x, i, j, ...) <- value
}
\arguments{
\item{x}{- a BibEntry object.}

\item{i}{- see \code{\link{[.BibEntry}}}

\item{j}{- see \code{\link{[.BibEntry}}}

\item{...}{- see \code{\link{[.BibEntry}}}

\item{value}{- values to be assigned to \code{x}.  To update one entry only,
should be a named character vector with names corresponding to fields.  To
update multiple entries, should be a list of named character vectors.  Can
also be an object of class BibEntry.}
}
\value{
an object of class BibEntry.
}
\description{
Assign new values for specified fields in a BibEntry object using a named
character vector or list of named character vectors.
}
\details{
Date and name list fields should be in the format expected
by Biblatex (see \code{\link{BibEntry}}).

To clear a field \sQuote{field_name} from an entry use \code{field_name = ""}.
}
\examples{
if (requireNamespace("bibtex")) {
    file.name <- system.file("Bib", "RJC.bib", package="RefManageR")
    bib <- ReadBib(file.name)
    print(bib[seq_len(3L)], .opts = list(sorting = "none", bib.style = "alphabetic"))
    ## add month to Serban et al., add URL and urldate to Jennings et al., and
    ##   add DOI and correct journal to Garcia et al.
    bib[seq_len(3L)] <- list(c(date="2013-12"),
                            c(url="https://bsb.eurasipjournals.com/content/2013/1/13",
                              urldate = "2014-02-02"),
                            c(doi="10.1093/bioinformatics/btt608",
                              journal = "Bioinformatics"))
    print(bib[seq_len(3L)], .opts = list(sorting = "none", bib.style = "alphabetic"))
    bib2 <- bib[seq_len(3L)]
    bib2[2:3] <- bib[5:6]
    bib2
    bib2[3] <- c(journal='', eprinttype = "arxiv", eprint = "1308.5427",
      eprintclass = "math.ST", pubstate = "submitted", bibtype = "misc")
    bib2
}
}
\seealso{
Other operators: 
\code{\link{$.BibEntry}()},
\code{\link{$<-.BibEntry}()},
\code{\link{+.BibEntry}()},
\code{\link{[.BibEntry}()},
\code{\link{[[.BibEntry}()},
\code{\link{[[<-.BibEntry}()},
\code{\link{c.BibEntry}()}
}
\concept{operators}
\keyword{manip}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadPubMed.R
\name{LookupPubMedID}
\alias{LookupPubMedID}
\title{Retrieve PubMed ID's for a BibEntry object}
\usage{
LookupPubMedID(bib, index)
}
\arguments{
\item{bib}{a bibentry object}

\item{index}{indices specifying which entries of \code{bib} will be
searched for.  If \code{missing}, all entries
are searched for.}
}
\value{
a BibEntry object - \code{bib} with additional eprinttype and eprint
fields when the search is successful
for an entry.
}
\description{
Uses the NCBI E-utilities to to search for PubMed ID's for citations
stored in a BibEntry object.
}
\details{
For each entry a citation string is created using the fields
journaltitle/journal, date/year,
  volume, pages, and author; and these strings are then used to search the
NCBI database for PubMed ID's.

  If an ID is found for an entry, the entry is updated so that the eprinttype
field is assigned the value
  \dQuote{pubmed} and the eprint field is assigned the ID.
}
\examples{
if (interactive() && !httr::http_error("https://eutils.ncbi.nlm.nih.gov/")){
  file.name <- system.file("Bib", "RJC.bib", package = "RefManageR")
  bib <- ReadBib(file.name)
  LookupPubMedID(bib[[101:102]])
}
}
\seealso{
Other pubmed: 
\code{\link{GetPubMedByID}()},
\code{\link{GetPubMedRelated}()},
\code{\link{ReadCrossRef}()},
\code{\link{ReadPubMed}()}
}
\concept{pubmed}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/relist.R, R/unlist.R
\name{RelistBibEntry}
\alias{RelistBibEntry}
\alias{unlist.BibEntry}
\title{Flatten and unflatten BibEntry objects}
\usage{
RelistBibEntry(flesh, skeleton = NULL)

\method{unlist}{BibEntry}(x, recursive = FALSE, use.names = TRUE)
}
\arguments{
\item{flesh}{list; an \code{unlist}ed BibEntry object}

\item{skeleton}{currently ignored}

\item{x}{a BibEntry object to flatten}

\item{recursive}{ignored.}

\item{use.names}{ignored.}
}
\value{
\code{RelistBibEntry} - an object of class BibEntry

For \code{unlist}, a list with bib entries collapsed into a single list.
}
\description{
\code{RelistBibEntry} unflattens a BibEntry object that has been
flattened with \code{unlist}.

\code{unlist} flattens a BibEntry object to a single list where every field (including \code{bibtype} and \code{key})
of every entry is a separate element in the list.
}
\details{
\code{RelistBibEntry} is only intended for use with
\code{unlist}ed BibEntry objects.
}
\note{
The names of the list elements from an unlisted BibEntry object will not be unique.  To do this see \code{\link{make.unique}}.
}
\examples{
bib <- list(c(bibtype = "article", key = "mclean2014a", title = "My New Article",
  author = "Mathew W. McLean", journaltitle = "The Journal", date = "2014-01"),
  c(bibtype = "article", key = "mclean2014b", title = "My Newer Article",
  author = "Mathew W. McLean", journaltitle = "The Journal", date = "2014-02"))
bib <- as.BibEntry(bib)
unlist(bib)
RelistBibEntry(unlist(bib))
}
\seealso{
\code{\link{as.BibEntry}}
}
\keyword{database}
\keyword{list}
\keyword{manip}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BibEntryExtractOp.R, R/SearchBib.R
\name{[.BibEntry}
\alias{[.BibEntry}
\alias{SearchBib}
\title{Search BibEntry objects by field}
\usage{
\method{[}{BibEntry}(x, i, j, ..., drop = FALSE)

SearchBib(x, .opts = list(), ...)
}
\arguments{
\item{x}{an object of class BibEntry}

\item{i}{A named list or character vector of search terms with names
corresponding to the field to search for the
search term.  Alternatively, a vector of entry key values or numeric or
logical indices specifying which entries to extract.}

\item{j}{A named list or character vector, as \code{i}.  Entries matching the
search specified by i \emph{OR} matching
the query specified by \code{j} will be return}

\item{...}{arguments in the form \code{bib.field = search.term}, or as \code{j}
list\emph{s} or character vector\emph{s} for additional searches.  For
\code{SearchBib}, can alternatively have same form as \code{i}.}

\item{drop}{logical, should attributes besides class be dropped from result?}

\item{.opts}{list of search options with \code{name = value} entries.  Any option described 
in \code{\link{BibOptions}} is valid, with the following being the most relevant ones
\itemize{
\item \code{use.regex} - logical; are the search terms regular expressions or should exact matching be used?
\item \code{ignore.case} - logical; should case be ignored when comparing strings?
\item \code{match.date} - how should the date fields date, urldate, eventdate, and origdate.  Default is \dQuote{year.only}, so 
that months and days in dates are ignored when comparing.  Currently, specifying any other value results the full date being
used.  See the Note section.
\item \code{match.author} - character string; how should name fields be searched? If \dQuote{family.only}, only family names are
compared; if \dQuote{family.with.initials}, family name and given name initials are used; if \dQuote{exact}, full 
names  are used.
\item \code{return.ind} - logical; if TRUE the returned object is numeric indices of match locations; otherwise, a BibEntry
object is returned
}}
}
\value{
an object of class BibEntry (the results of the search/indexing),
\emph{or} if \code{BibOptions()$return.ind=TRUE}, the indices in \code{x} that
match the search terms.
}
\description{
Allows for searching and indexing a BibEntry object by fields, including
names and dates.  The extraction operator and the \code{SearchBib} function
simply provide different interfaces to the same search functionality.
}
\note{
The arguments to the SearchBib function that control certain search
features can also be changed for the extraction
operator by changing the corresponding option in the .BibOptions object; see
\code{\link{BibOptions}}.
}
\examples{
if (requireNamespace("bibtex")) {
    file.name <- system.file("Bib", "biblatexExamples.bib", package="RefManageR")
    bib <- suppressMessages(ReadBib(file.name))

    ## author search, default is to use family names only for matching
    bib[author = "aristotle"]

    ## Aristotle references before 1925
    bib[author="aristotle", date = "/1925"]

    ## Aristotle references before 1925 *OR* references with editor Westfahl
    bib[list(author="aristotle", date = "/1925"),list(editor = "westfahl")]

    ## Change some searching and printing options and search for author
    old.opts <- BibOptions(bib.style = "authoryear", match.author = "exact",
      max.names = 99, first.inits = FALSE)
    bib[author="Mart\u00edn, Jacinto and S\u00e1nchez, Alberto"]
    BibOptions(old.opts)  ## reset options

    ## Some works of Raymond J. Carroll's
    file.name <- system.file("Bib", "RJC.bib", package="RefManageR")
    bib <- ReadBib(file.name)
    length(bib)

    ## index by key
    bib[c("chen2013using", "carroll1978distributions")]

    ## Papers with someone with family name Wang
    length(SearchBib(bib, author='Wang', .opts = list(match.author = "family")))

    ## Papers with Wang, N.
    length(SearchBib(bib, author='Wang, N.', .opts = list(match.author = "family.with.initials")))

    ## tech reports with Ruppert
    length(bib[author='ruppert',bibtype="report"])

    ##Carroll and Ruppert tech reports at UNC
    length(bib[author='ruppert',bibtype="report",institution="north carolina"])

    ## Carroll and Ruppert papers since leaving UNC
    length(SearchBib(bib, author='ruppert', date="1987-07/",
       .opts = list(match.date = "exact")))
}

## Carroll and Ruppert papers NOT in the 1990's
\dontrun{
if (requireNamespace("bibtex")) {
    length(SearchBib(bib, author='ruppert', date = "!1990/1999"))
    identical(SearchBib(bib, author='ruppert', date = "!1990/1999"),
      SearchBib(bib, author='ruppert', year = "!1990/1999"))
    table(unlist(SearchBib(bib, author='ruppert', date="!1990/1999")$year))

    ## Carroll + Ruppert + Simpson
    length(bib[author="Carroll, R. J. and Simpson, D. G. and Ruppert, D."])

    ## Carroll + Ruppert OR Carroll + Simpson
    length(bib[author=c("Carroll, R. J. and Ruppert, D.", "Carroll, R. J. and Simpson, D. G.")])

    ## Carroll + Ruppert tech reports at UNC "OR" Carroll and Ruppert JASA papers
    length(bib[list(author='ruppert',bibtype="report",institution="north carolina"),
      list(author="ruppert",journal="journal of the american statistical association")])
}
}
}
\seealso{
Other operators: 
\code{\link{$.BibEntry}()},
\code{\link{$<-.BibEntry}()},
\code{\link{+.BibEntry}()},
\code{\link{[<-.BibEntry}()},
\code{\link{[[.BibEntry}()},
\code{\link{[[<-.BibEntry}()},
\code{\link{c.BibEntry}()}
}
\concept{operators}
\keyword{database}
\keyword{list}
\keyword{manip}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/09sort.R
\name{sort.BibEntry}
\alias{sort.BibEntry}
\title{Sort a BibEntry Object}
\usage{
\method{sort}{BibEntry}(
  x,
  decreasing = FALSE,
  sorting = BibOptions()$sorting,
  .bibstyle = BibOptions()$bib.style,
  ...
)
}
\arguments{
\item{x}{an object of class BibEntry}

\item{decreasing}{logical; should the sort be increasing or decreasing?}

\item{sorting}{sort method to use, see \bold{Details}.}

\item{.bibstyle}{bibliography style; used when \code{sort} is called by
\code{\link{print.BibEntry}}}

\item{...}{internal use only}
}
\value{
the sorted BibEntry object
}
\description{
Sorts a \code{BibEntry} object by specified fields.  The possible fields used
for sorting and the order they are used in correspond with the options
available in BibLaTeX.
}
\details{
The possible values for argument \code{sorting} are
\itemize{
\item nty - sort by name, then by title, then by year
\item nyt - sort by name, then by year, then title
\item nyvt - sort by name, year, volume, title
\item anyt - sort by alphabetic label, name, year, title
\item anyvt - sort by alphabetic label, name, year, volume, title
\item ynt - sort by year, name, title
\item ydnt - sort by year (descending), name, title
\item debug - sort by keys
\item none - no sorting is performed
}

All sorting methods first consider the field presort, if available.
Entries with no presort field are assigned presort
value \dQuote{mm}. Next the sortkey field is used.

When sorting by name, the sortname field is used first.  If it is not present,
the author field is used,
if that is not present editor is used, and if that is not present translator is
used.  All of these fields are affected
by the value of \code{max.names} in .BibOptions()$max.names.

When sorting by title, first the field sorttitle is considered.  Similarly,
when sorting by year, the field sortyear is
first considered.

When sorting by volume, if the field is present it is padded to four digits
with leading zeros; otherwise, the string \dQuote{0000} is used.

When sorting by alphabetic label, the labels that would be generating with
the \dQuote{alphabetic} bibstyle are used.  First the shorthand field is
considered, then label, then shortauthor, shorteditor, author, editor,
and translator.  Refer to the BibLaTeX manual Sections 3.1.2.1 and 3.5 and
Appendix C.2 for more information.
}
\examples{
if (requireNamespace("bibtex")) {
    file.name <- system.file("Bib", "biblatexExamples.bib", package="RefManageR")
    bib <- suppressMessages(ReadBib(file.name)[[70:73]])
    BibOptions(sorting = "none")
    bib
    sort(bib, sorting = "nyt")
    sort(bib, sorting = "ynt")
    BibOptions(restore.defaults = TRUE)
}
}
\references{
Lehman, Philipp and Kime, Philip and Boruvka, Audrey and
Wright, J. (2013). The biblatex Package.
\url{https://mirror.pregi.net/tex-archive/macros/latex/contrib/biblatex/doc/biblatex.pdf}.
}
\seealso{
\code{\link{BibEntry}}, \code{\link{print.BibEntry}}, \code{\link{order}}
}
\keyword{manip}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/levels.R
\name{levels.BibEntry}
\alias{levels.BibEntry}
\alias{fields}
\title{Extract all fields present in a BibEntry object}
\usage{
\method{levels}{BibEntry}(x)

fields(x)
}
\arguments{
\item{x}{a BibEntry object.}
}
\value{
a list with the same length as \code{x} of character vectors giving the 
fields present in each entry of a BibEntry object.
}
\description{
These functions return a list of all fields present in a BibEntry object.
}
\note{
The only difference between \code{fields} and \code{levels} is that
\code{levels} returns a list with element names corresponding to entry keys.
}
\examples{
bib <- as.BibEntry(list(c(bibtype = "Article", key = "mclean2014a", title = "My New Article", 
  author = "Mathew W. McLean", 
  journaltitle = "The Journal", date = "2014-01"), c(bibtype = "Book", key = "mclean2014b", 
  title = "My New Book", editor = "Mathew W. McLean", ISBN = "247123837", date = "2014-02")))       
fields(bib)
levels(bib)
}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadPubMed.R
\name{GetPubMedRelated}
\alias{GetPubMedRelated}
\title{Retrieve related articles from PubMed using PubMed ID's}
\usage{
GetPubMedRelated(
  id,
  database = "pubmed",
  batch.mode = TRUE,
  max.results = 10,
  return.sim.scores = FALSE,
  return.related.ids = FALSE
)
}
\arguments{
\item{id}{either a character vector of PubMed ID's or a BibEntry object,
which is expected to have at least some entries with
\code{eprinttype = "pubmed"} and eprint field specifying a PubMed ID.}

\item{database}{string; the Entrez database to search}

\item{batch.mode}{logical; if \code{TRUE}, the PubMed IDs in \code{id}
are combined by Entrez when searching for linked
IDs so that only one set of linked IDs is returned.  If \code{FALSE}, a
set of linked IDs is obtained for each ID
in \code{id}.
will be returned}

\item{max.results}{numeric vector; the maximum number of results to
return if \code{batch.mode} \code{TRUE}; or if \code{batch.mode} is
\code{FALSE}, this should have the same length
as \code{id} with each element giving the maximum number of results to
return for the corresponding ID.}

\item{return.sim.scores}{logical; Entrez returns a similarity score with
each returned citation giving a measure of how similar the returned entry
is to the ones specified by the query.  If \code{TRUE} these scores are added
to the returned BibEntry object in a field called \sQuote{score}.}

\item{return.related.ids}{logical; should the original PubMed ID(s) that a
returned entry is related to be stored in a field called \sQuote{PMIDrelated}.}
}
\value{
an object of class BibEntry.
}
\description{
Searches PubMed for articles related to a set of PubMed ID's using
NCBI's E-Utilities.
}
\examples{
if (interactive() && !httr::http_error("https://eutils.ncbi.nlm.nih.gov/")){
  file.name <- system.file("Bib", "RJC.bib", package="RefManageR")
  bib <- ReadBib(file.name)
  bib <- LookupPubMedID(bib[[101:102]])
  toBiblatex(GetPubMedRelated(bib, batch.mode = TRUE, max.results = 2,
  return.sim.scores = TRUE, return.related.ids = TRUE))
  GetPubMedRelated(bib, batch.mode = FALSE, max.results = c(2, 2))
}
}
\references{
\url{https://www.ncbi.nlm.nih.gov/books/NBK25500/}
}
\seealso{
Other pubmed: 
\code{\link{GetPubMedByID}()},
\code{\link{LookupPubMedID}()},
\code{\link{ReadCrossRef}()},
\code{\link{ReadPubMed}()}
}
\concept{pubmed}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asdataframe.R
\name{as.data.frame.BibEntry}
\alias{as.data.frame.BibEntry}
\title{Coerce to a Data Frame}
\usage{
\method{as.data.frame}{BibEntry}(x, row.names = NULL, optional = FALSE, ...)
}
\arguments{
\item{x}{- a BibEntry object}

\item{row.names}{- ignored}

\item{optional}{- ignored}

\item{...}{- ignored}
}
\value{
a data.frame object with row names giving the keys, and first column giving entry type.
}
\description{
Coerces a BibEntry object to a data.frame, with each row of the data frame being a field present in at least one
entry in the BibEntry object being coerced.
}
\examples{
bib <- list(c(bibtype = "article", key = "mclean2014a", title = "My New Article",
  author = "Mathew W. McLean", journaltitle = "The Journal", date = "2014-01"),
  c(bibtype = "article", key = "mclean2014b", volume = 10, title = "My Newer Article",
  author = "Mathew W. McLean", journaltitle = "The Journal", date = "2014-02"))
bib <- as.BibEntry(bib)
as.data.frame(bib)
}
\seealso{
\code{\link{BibEntry}}, \code{\link{as.BibEntry}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BibEntryAddOp.R
\name{+.BibEntry}
\alias{+.BibEntry}
\alias{merge.BibEntry}
\title{Merge two BibEntry objects while discarding duplicates}
\usage{
\method{+}{BibEntry}(e1, e2)

\method{merge}{BibEntry}(
  x,
  y,
  fields.to.check = BibOptions()$merge.fields.to.check,
  ignore.case = BibOptions()$ignore.case,
  ...
)
}
\arguments{
\item{e1}{BibEntry object}

\item{e2}{BibEntry object to be merged with e1}

\item{x}{BibEntry object}

\item{y}{BibEntry object}

\item{fields.to.check}{character vector; which BibLaTeX fields should be
checked to determine if an entry
is a duplicate?  Can include \code{"bibtype"} to check entry type and
\code{"key"} to check entry keys. Specifying \code{"all"} checks all fields
using \code{\link{duplicated}}.}

\item{ignore.case}{logical; if \code{TRUE}, case is ignored when determining
if fields are duplicates.}

\item{...}{ignored}
}
\value{
an object of class BibEntry
}
\description{
Merges two BibEntry objects comparing only the specified fields to detect
duplicates, thus it is can be made less strict
than using \code{duplicated}, \code{unique}, etc.  Attributes are also merged
and keys are ensured to be unique.
\code{merge} and \code{+} simply provide different interfaces for merging.
}
\examples{
if (requireNamespace("bibtex")) {
    file.name <- system.file("Bib", "biblatexExamples.bib", package="RefManageR")
    bib <- suppressMessages(ReadBib(file.name))
    bib1 <- bib[seq_len(44)]
    bib2 <- bib[45:length(bib)]

    ## The following is FALSE because the parent entry of one entry in bib1
    ##   is in bib2, so the child entry is expanded in the BibEntry object
    ##   returned by `[` to include the fields inherited from the dropped parent
    identical(merge(bib1, bib2, 'all'), bib)
    toBiblatex(bib1[[1L]])
    toBiblatex(bib[[1L]])

    ## Alternatively, the operator `[[` for BibEntry objects does not expand
    ##   cross references
    bib1 <- bib[[seq_len(44)]]
    bib2 <- bib[[45:length(bib)]]
    identical(merge(bib1, bib2, 'all'), bib)

    ## Not strict enough
    invisible(merge(bib1, bib2, c('title', 'date')))
 }

## New publications of R.J. Carroll from Google Scholar and Crossref
\dontrun{
if (requireNamespace("bibtex")) {
    bib1 <- ReadGS(scholar.id = "CJOHNoQAAAAJ", limit = '10', sort.by.date = TRUE)
    bib2 <- ReadCrossRef(query = "rj carroll", limit = 10, sort = "relevance",
      min.relevance = 80)
    oldopt <- BibOptions(merge.fields.to.check = "title")
    rjc.new.pubs <- bib1 + bib2
    BibOptions(oldopt)
}
}
}
\seealso{
\code{\link{duplicated}}, \code{\link{unique}}

Other operators: 
\code{\link{$.BibEntry}()},
\code{\link{$<-.BibEntry}()},
\code{\link{[.BibEntry}()},
\code{\link{[<-.BibEntry}()},
\code{\link{[[.BibEntry}()},
\code{\link{[[<-.BibEntry}()},
\code{\link{c.BibEntry}()}
}
\author{
McLean, M. W. \email{mathew.w.mclean@gmail.com}
}
\concept{operators}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BibEntryAssignOp.R
\name{$<-.BibEntry}
\alias{$<-.BibEntry}
\title{Replace values for a particular field in a BibEntry object}
\usage{
\method{$}{BibEntry}(x, name) <- value
}
\arguments{
\item{x}{a BibEntry object}

\item{name}{string; the field to assign the new values to.}

\item{value}{character vector; the replacement field values to be assigned.}
}
\value{
an object of class BibEntry with the updated fields.
}
\description{
Used to replace the values stored for a specified field in a BibEntry object.
}
\note{
The method expects date and name list fields to be in the format
expected by Biblatex.  The 
field specified by \code{name} does not have to be one currently in \code{x}.
}
\examples{
bib <- BibEntry(bibtype = "misc", key = "mclean", author = "Mathew W. McLean", 
  title = "My Work", year = "2012")
bib$year <- 2014
bib$author <- "McLean, M. W. and Carroll, R. J." 
bib$url <- "https://example.com"
bib

bib <- c(bib, as.BibEntry(citation()))
bib[1]$author[2] <- person(c("Raymond", "J."), "Carroll")
bib$author
}
\seealso{
Other operators: 
\code{\link{$.BibEntry}()},
\code{\link{+.BibEntry}()},
\code{\link{[.BibEntry}()},
\code{\link{[<-.BibEntry}()},
\code{\link{[[.BibEntry}()},
\code{\link{[[<-.BibEntry}()},
\code{\link{c.BibEntry}()}
}
\concept{operators}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WriteBib.R
\name{WriteBib}
\alias{WriteBib}
\title{Create a BibTeX File from a BibEntry Object
e
Creates a Bibtex File from a BibEntry object for use with either BibTeX
or BibLaTex.}
\usage{
WriteBib(
  bib,
  file = "references.bib",
  biblatex = TRUE,
  append = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{bib}{a BibEntry object to be written to file}

\item{file}{character string naming a file, should; end in \dQuote{.bib}.
Can be \code{NULL}, in which case the BibEntry object will be written
to \code{\link{stdout}}.}

\item{biblatex}{boolean; if \code{TRUE}, \code{\link{toBiblatex}} is used
and no conversions of the BibEntry object
are done; if \code{FALSE} entries will be converted as described in
\code{\link{toBibtex.BibEntry}}.}

\item{append}{as in \code{write.bib} in package \code{bibtex}}

\item{verbose}{as in \code{write.bib} in package \code{bibtex}}

\item{...}{additional arguments passed to \code{\link{writeLines}}}
}
\value{
\code{bib} - invisibly
}
\description{
Create a BibTeX File from a BibEntry Object
e
Creates a Bibtex File from a BibEntry object for use with either BibTeX
or BibLaTex.
}
\note{
To write the contents of \code{bib} \dQuote{as is}, the argument
\code{biblatex} should be \code{TRUE}, otherwise
conversion is done as in \code{\link{toBibtex.BibEntry}}.
}
\examples{
if (requireNamespace("bibtex")){
    bib <- BibEntry("Article", key = "Carroll_2012",
                    doi = "10.1080/01621459.2012.699793",
                    year = "2012", month = "sep",
                    volume = 107, number = 499,
                    pages = {1166--1177},
      author = "R. Carroll and A. Delaigle and P. Hall",
      title = "Deconvolution When Classifying Noisy Data ...",
      journal = "Journal of the American Statistical Association")

  ## Write bib if no server error and bibtex available
  if (length(bib)){
    tfile <- tempfile(fileext = ".bib")
    WriteBib(bib, tfile, biblatex = TRUE)
    identical(ReadBib(tfile), bib)
    unlink(tfile)
  }
}
}
\seealso{
\code{write.bib} in package \code{bibtex}, \code{\link{ReadBib}},
\code{\link{toBibtex.BibEntry}}, \code{\link{toBiblatex}},
\code{\link{BibEntry}}
}
\author{
McLean, M. W. based on \code{write.bib} by Gaujoux, R.
in package \code{bibtex}.
}
\keyword{IO}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadGS.R
\name{ReadGS}
\alias{ReadGS}
\title{Import book and article references from a public Google Scholar profile by ID.}
\usage{
ReadGS(
  scholar.id,
  start = 0,
  limit = 100,
  sort.by.date = FALSE,
  .Encoding = "UTF-8",
  check.entries = BibOptions()$check.entries
)
}
\arguments{
\item{scholar.id}{character; the Google Scholar ID from which citations will
be imported.  The ID can by found by
visiting an author's Google Scholar profile and noting the value in the uri
for the \dQuote{user} parameter.}

\item{start}{numeric; index of first citation to include.}

\item{limit}{numeric; maximum number of results to return.  Cannot exceed 100.}

\item{sort.by.date}{boolean; if true, newest citations are imported first;
otherwise, most cited works are imported first.}

\item{.Encoding}{character; text encoding to use for importing the results
and creating the bib entries.}

\item{check.entries}{What should be done with incomplete entries (those
containing \dQuote{...} due to long fields)?
Either \code{FALSE} to add them anyway, \code{"warn"} to add with a warning,
or any other value to drop the entry
with a message and continue processing the remaining entries.}
}
\value{
An object of class BibEntry.  If the entry has any citations, the number of
citations is stored in a field \sQuote{cites}.
}
\description{
This function will create a BibEntry object for up to 100 references from a
provided Google Scholar ID,
if the profile is public.  The number of citations for each entry will
also be imported.
}
\details{
This function creates \code{BibTeX} entries from an author's
Google Scholar page.
If the function finds numbers corresponding to volume/number/pages of a journal
article, an \sQuote{Article} entry
is created.  If an arXiv identifier is found, a \sQuote{Misc} entry is created
with \code{eprint}, \code{eprinttype}, and \code{url} fields.  Otherwise, a
\sQuote{TechReport} entry is created; unless the entry has more than ten citations,
in which case a \sQuote{Book} entry is created.

Long author lists, long titles, and long journal/publisher names can all lead to
these fields being incomplete for
a particular entry.  When this occurs, these entries are either dropped or added
with a warning depending on the value of the \code{check.entries} argument.
}
\note{
Read Google's Terms of Service before using.

It is not possible to automatically import BibTeX entries directly from Google
Scholar as no API is available and this violates their Terms of Service.
}
\examples{
if (interactive() && !httr::http_error("https://scholar.google.com")){
  ## R. J. Carroll's ten newest publications
  ReadGS(scholar.id = "CJOHNoQAAAAJ", limit = 10, sort.by.date = TRUE)

  ## Matthias Katzfu\ss
  BibOptions(check.entries = "warn")
  kat.bib <- ReadGS(scholar.id = "vqW0UqUAAAAJ")

  ## retrieve GS citation counts stored in field 'cites'
  kat.bib$cites
}
}
\seealso{
\code{\link{BibEntry}}
}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UpdateFieldName.R
\name{UpdateFieldName}
\alias{UpdateFieldName}
\title{Rename a field in a BibEntry object.}
\usage{
UpdateFieldName(x, old.field, new.field)
}
\arguments{
\item{x}{- a BibEntry object}

\item{old.field}{- string; the current name of the field to be renamed}

\item{new.field}{- string; the new name to replace \code{old.field}}
}
\value{
\code{x}, with the renamed field.
}
\description{
This function will rename a field, in every entry where it is present, in a 
BibEntry object.
}
\examples{
bib <- as.BibEntry(list(c(bibtype = "article", key = "mclean2014a", title = "My New Article", 
  author = "Mathew W. McLean", journal = "The Journal", date = "2014-01"), 
  c(bibtype = "article", key = "mclean2014b", title = "My Newer Article", 
    author = "Mathew W. McLean", journal = "The Journal", date = "2014-02")))       
bib <- UpdateFieldName(bib, "journal", "journaltitle")
toBiblatex(bib)   
}
\keyword{manip}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/open.BibEntry.R
\name{open.BibEntry}
\alias{open.BibEntry}
\title{Open BibEntry in PDF viewer or web browser.}
\usage{
\method{open}{BibEntry}(
  con,
  entry = 1L,
  open.field = c("file", "url", "eprint", "doi"),
  viewer,
  ...
)
}
\arguments{
\item{con}{BibEntry object to extract connections from.}

\item{entry}{numeric index or character key of entry in \code{bib} to open.}

\item{open.field}{character vector of fields to use in \code{bib} to open
the BibEntry.  Possible fields are any combination of \dQuote{file},
\dQuote{url}, \dQuote{eprint}, or \dQuote{doi}.  \dQuote{eprint} is
implemented for \code{eprinttype=} \dQuote{JSTOR}, \dQuote{PubMed},
or \dQuote{arXiv}.  When multiple fields are specified, they are tried in
the order they appear in the vector.}

\item{viewer}{character string giving the name of the program to be used
as hypertext browser.  It should be in the PATH, or a full path specified.
Alternatively, an R function to be called to invoke the browser.  Defaults
to \code{getOptions("pdfviewer")} if \code{open.field = "file"} and
\code{getOptions("browser")}, otherwise.}

\item{...}{not used.}
}
\description{
Attempts to open a connection to an entry in a BibEntry object using fields
such as \sQuote{file}, \sQuote{DOI}, \sQuote{eprint} + \sQuote{eprinttype},
and \sQuote{URL}.
}
\examples{
\dontrun{
if (requireNamespace("bibtex")) {
    testbib <- ReadBib(system.file("REFERENCES.bib", package="bibtex"))
    open(testbib)
    testbib$file <- file.path(R.home("doc/manual"), "R-intro.pdf")
    open(testbib)
}
}
}
\seealso{
\code{\link{browseURL}}
}
\author{
McLean, M. W. \email{mathew.w.mclean@gmail.com}
}
\keyword{connection}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RefManageR-package.R
\docType{package}
\name{RefManageR-package}
\alias{RefManageR-package}
\alias{RefManageR}
\alias{refmanager}
\title{Import and Manage BibTeX and BibLaTeX references with RefManageR}
\description{
RefManageR provides tools for importing and working with
bibliographic references.  It greatly enhances the bibentry class by
providing a class BibEntry which stores BibTeX and BibLaTeX references,
supports UTF-8 encoding, and can be easily searched by any field, by date
ranges, and by various formats for name lists (author by last names,
translator by full names, etc.). Entries can be updated, combined, sorted, printed
in a number of styles, and exported. BibTeX and BibLaTeX .bib files can be
read into R and converted to BibEntry objects.  Interfaces to NCBI's
Entrez, CrossRef, and Zotero are provided for importing references and
references can be created from locally stored PDFs using Poppler.  Includes
functions for citing and generating a bibliography with hyperlinks for
documents prepared with RMarkdown or RHTML.
}
\details{
\bold{Importing and Creating References}

BibEntry objects can be created directly using the \code{\link{BibEntry}} function.  \code{.bib} files can be read into R
using the \code{\link{ReadBib}} function.
Tools are provided for importing references from Crossref, Zotero, Google Scholar, 
and PDFs and looking up PubMed ID's and DOIs.  See \code{\link{ReadPDFs}}, \code{\link{ReadZotero}}, \code{\link{ReadCrossRef}}, \code{\link{ReadGS}},
\code{\link{ReadPubMed}}, \code{\link{GetPubMedByID}}, \code{\link{GetPubMedRelated}}.

\bold{Manipulating BibEntry objects}

BibEntry objects may be searched and indexed by field values, name lists, keys, dates, date ranges, etc.  
See \code{\link{[.BibEntry}}, \code{\link{[<-.BibEntry}}, \code{\link{[[.BibEntry}}, \code{\link{$.BibEntry}}.

\bold{Printing and Exporting Bibliographies}

The \code{\link{print.BibEntry}} function can print in a number of formats (e.g. text, html) and most of the 
base bibliography styles available with BibLaTeX (e.g. alphabetic, numeric, authortitle, and authoryear).  
\code{\link{toBibtex.BibEntry}} will convert a BibEntry object to a character vector containing lines of 
a BibTeX file, converting fields, entry types and expanding crossreferences as needed to coerce BibLaTeX entries to
BibTeX.  \code{\link{toBiblatex}} converts the BibEntry object to a character vector containing lines of 
the corresponding BibLaTeX file.  The results can be written to a file using \code{\link{WriteBib}}.

Citations can be generated in a number of styles using one of the available functions for 
citations.  A list of references can be printed based on the works the user has cited thus far
in their document.  See \code{\link{Cite}}.  The citations and bibliography can be printed 
including hyperlinks using either the R Markdown or R HTML formats.

\bold{Additional features}

All sorting methods for bibliographies available in the BibLaTeX LaTeX package have been implemented see 
\code{\link{sort.BibEntry}} and the references.

Using \code{\link{open.BibEntry}} electronic copies of references can be opened in a PDF viewer or web browser.

The convenience function \code{\link{BibOptions}} is provided for setting defaults for commonly used
functions such as \code{\link{print.BibEntry}}, \code{\link{[.BibEntry}}, and 
\code{\link{Cite}}.  Its interface is similar to \code{\link{options}}.
}
\references{
McLean, M. W. (2014). Straightforward Bibliography Management in R Using the RefManageR Package.
\href{https://arxiv.org/abs/1403.2036}{arXiv: 1403.2036 [cs.DL]}. Submitted.

Lehman, P., P. Kime, A. Boruvka, and J. Wright (2013). The biblatex Package.
\url{https://mirror.pregi.net/tex-archive/macros/latex/contrib/biblatex/doc/biblatex.pdf}.

Hornik, K., D. Murdoch, and A. Zeileis (2012). 
Who Did What? The Roles of R Package Authors and How to Refer to Them. The R Journal \bold{4}, 1.
\url{https://journal.r-project.org/archive/2012-1/RJournal_2012-1_Hornik~et~al.pdf}

Patashnik, O (1988). Bibtexing. \url{https://mirror.pregi.net/tex-archive/biblio/bibtex/contrib/doc/btxdoc.pdf}.
}
\author{
McLean, M. W. \email{mathew.w.mclean@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02BibOptions.R
\name{BibOptions}
\alias{BibOptions}
\title{Set options/hooks for RefManageR}
\usage{
BibOptions(..., restore.defaults = FALSE)
}
\arguments{
\item{...}{a character vector or strings specifying option names to access; or to set options values, 
a named list or vector of option values or options specified in name=value pairs.}

\item{restore.defaults}{logical; if TRUE, \code{...}'s are ignored and all package options are restored to their
defaults.}
}
\value{
if a vector of option names is supplied, the current value of the requested options, or if \code{...} is missing,
all current option values; otherwise, when setting options the old values of the changed options are (invisibly) 
returned as a list.
}
\description{
This function is used to access and set package options for RefManageR, similar to \code{\link{options}}.  
The options are listed in the details
}
\details{
The following are valid package options. 

\bold{Options for searching/indexing a BibEntry object.  See \code{\link{[.BibEntry}} and 
\code{\link{[<-.BibEntry}}}
\enumerate{
\item \code{match.author} - string; controls how name list fields (author, editor, translator, etc.) are matched 
when searching for names.
\dQuote{family.with.initials} require family names and given name initials to match, \dQuote{exact} requires names to match
exactly, and any other value results in only family names being compared (the default).
\item \code{match.date} - string; controls how date fields are matched when searching.  If \dQuote{year.only} (the default),
only years are checked for equality when comparing dates, otherwise months and days will also be compared,
if they are available.
\item \code{use.regex} - logical; if \code{TRUE}, regular expressions are used when searching non-date fields; otherwise, exact
matching is used.
\item \code{ignore.case} - logical; if \code{TRUE}, case is ignored when searching.
\item \code{return.ind} - logical; if \code{TRUE} the return value of \code{\link{SearchBib}} and the operators 
\code{\link{[.BibEntry}}, will be the indices of any matches; otherwise, a \code{BibEntry}
object is returned.
}

\bold{Options for Printing with \code{\link{print.BibEntry}} and \code{\link{PrintBibliography}}}
\enumerate{
\item \code{bib.style} - string; Biblatex bibliography style to use when printing and formatting a BibEntry object.  Possible
values are \dQuote{numeric} (default), \dQuote{authoryear}, \dQuote{authortitle}, \dQuote{alphabetic}, \dQuote{draft}.
\item \code{first.inits} - logical; if \code{TRUE}, only given name initials are displayed when printing; otherwise, full names
are used.
\item \code{dashed} - logical; if \code{TRUE} and \code{bib.style = "authoryear"} or \code{bib.style = "authortitle"},
recurring author and editor names are replaced with \dQuote{---} when printing.
\item \code{sorting} - string; controls how BibEntry objects are sorted.  Possible values are \dQuote{nty}, \dQuote{nyt}, 
\dQuote{nyvt}, \dQuote{anyt}, \dQuote{anyvt}, \dQuote{ynt}, \dQuote{ydnt}, \dQuote{none}, \dQuote{debug};  see 
\code{\link{sort.BibEntry}}
\item \code{max.names} - numeric; maximum number of names to display before using
\dQuote{et al.} when formatting and printing name list fields. This is also the minimum
number of names that will be displayed if \dQuote{et al.} is used (corresponding
to the \sQuote{minnames} package option in Biblatex). See below and option
\code{longnamesfirst} when using this argument with the citation
functions (\code{\link{Citet}}, etc.). 
\item \code{no.print.fields} character vector; fields that should not be printed, 
e.g., doi, url, isbn, etc.
\item \code{style} - character string naming the printing style.  Possible values are 
plain text (style \dQuote{text}), BibTeX (\dQuote{Bibtex}), BibLaTeX (\dQuote{Biblatex}),
a mixture of plain text and BibTeX as 
traditionally used for citations (\dQuote{citation}), HTML (\dQuote{html}), 
LaTeX (\dQuote{latex}), \dQuote{markdown}, \dQuote{yaml}, 
R code (\dQuote{R}), and a simple copy of the textVersion elements 
(style \dQuote{textVersion}, see \code{\link{BibEntry}})
}

\bold{Options for the \code{\link{Cite}} functions}
\enumerate{
\item \code{cite.style} - character string; bibliography style to use to generate citations.  
\item \code{style} - as above, but used to format the citations.  
\item \code{hyperlink} - character string or logical; for use with \code{style = "markdown"}
and \code{style = "html"} (ignored otherwise).  If \code{FALSE}, no hyperlink
will be generated for the citation or in the bibliography when printing.  If
set equal to \code{"to.bib"}, then hyperlinks will be
generated linking the citation and bibliography.  The default value, \code{"to.doc"},
will try to create the hyperlink using the \code{url}, \code{doi}, or \code{eprint} fields of 
entry.  If these fields are not available, the hyperlink will point to the bibliography.  See
also \code{\link{open.BibEntry}}.
\item \code{super} - logical; should superscripts be used for numeric citations?  Ignored if
 \code{cite.style != "numeric"}.
\item \code{max.names} - numeric; same as above, except for citations. Note, that the
first time a reference is cited, this option will be ignored if \code{longfirstnames}
is \code{TRUE}.
\item \code{longnamesfirst} logical; should the first time a citation appears in the text
not be truncated at \code{max.names}?
\item \code{bibpunct} - character vector; punctuation to use in a citation.  The entries
in \code{bibpunct} are as follows
\enumerate{
\item The left delimiter for non-alphabetic and non-numeric citation styles
\item The right delimiter for non-alphabetic and non-numeric citation styles
\item The left delimiter for alphabetic and numeric citation styles
\item The right delimiter for alphabetic and numeric citation styles 
\item The separator between references in a citation.
\item Punctuation to go between the author and year.
}
}

\bold{Other}
\enumerate{
\item \code{check.entries} - string or \code{FALSE}; if \code{FALSE} entries are not checked to ensure that they have all the 
required fields for the type of entry; if \dQuote{warn} then entries are checked, but only a warning is issued and the 
entry is processed anyway; otherwise an error is produced if an entry does not have the required fields (default).  Note that
the majority of fields listed as required for a particular entry type in the Biblatex manual are not actually required for
Biblatex to produce an entry.
\item \code{merge.fields.to.check} - character vector; for \code{\link{merge.BibEntry}} and the operator \code{\link{+.BibEntry}},
the fields that should be checked when comparing entries for equality when merging BibEntry objects.  Specifying 
\dQuote{all} results in all fields be checked with \code{\link{duplicated}}.  The default is \dQuote{key} to only check for
duplicated keys.
}
}
\note{
If \code{...} is missing and \code{restore.defaults = FALSE}, all options and their current values will be returned
as a list.
}
\examples{
BibOptions()
BibOptions("first.inits", "bib.style")

oldopts <- BibOptions(first.inits = FALSE, bib.style = "authoryear")
oldopts
BibOptions(oldopts)

BibOptions(restore.defaults = TRUE)
}
\seealso{
\code{\link{print.BibEntry}}, \code{\link{BibEntry}}, \code{\link{options}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmdCite.R
\name{Cite}
\alias{Cite}
\alias{PrintBibliography}
\alias{TextCite}
\alias{AutoCite}
\alias{Citep}
\alias{Citet}
\alias{NoCite}
\title{Cite a BibEntry object in text and print all citations}
\usage{
Cite(bib, ..., textual = FALSE, before = NULL, after = NULL, .opts = list())

PrintBibliography(bib, .opts = list(), start = 1, end = length(bib))

Citep(bib, ..., before = NULL, after = NULL, .opts = list())

AutoCite(bib, ..., before = NULL, after = NULL, .opts = list())

Citet(bib, ..., before = NULL, after = NULL, .opts = list())

TextCite(bib, ..., before = NULL, after = NULL, .opts = list())

NoCite(bib, ..., .opts = list())
}
\arguments{
\item{bib}{a \code{BibEntry} or \code{bibentry} object}

\item{...}{passed to \code{\link{SearchBib}} for indexing into bib.  A character
vector of keys, for example.}

\item{textual}{logical; if TRUE, a \dQuote{textual} citation is produced, i.e.
what is produced by \\citet in \code{natbib} and \\textcite in \code{BibLaTeX};
otherwise, a parenthetical citation as \\citep and \\autocite.}

\item{before}{string; optional text to display before the citation.}

\item{after}{string; optional text to display after the citation.}

\item{.opts}{list; See the relevant section in \code{\link{BibOptions}} for a
description of all valid options for these functions.}

\item{start}{Integer; specifying the index of the first citation to
print. Useful for printing long bibliographies on multiple
pages/slides.}

\item{end}{Integer; specifying the index of the last citation to
print. Useful for printing long bibliographies on multiple
pages/slides.}
}
\value{
For the cite functions: a character string containing the citation

PrintBibliography: The formatted list of references.

NoCite: no return value; invoked for its side-effect.
}
\description{
The \code{Cite} functions allow for citing a \code{BibEntry} object in text.  The
\code{PrintBibliography} function allows for printing the bibliography of all
the cited entries.  The \code{NoCite} function adds references to the bibliography
without including a citation.  These functions are most useful when used in,
e.g., a RMarkdown or RHTML document.
}
\details{
See the package vignettes and execute the examples below.

If \code{bib.style = "alphabetic"} or \code{bib.style =
    "numeric"}, then sorting needs to be done at the start of the
    document prior to using a cite function as sorting is not done
    by the \code{PrintBibliography} function for those styles (specifying
    \code{sorting} in \code{.opts} is ignored in this case).  If no
    sorting is done, the references are listed in the order they
    were cited in for those two styles.

If the \code{...} argument to NoCite is identical to \dQuote{*}, then all
references in \code{bib} are added to the bibliography without citations.
}
\examples{
if (requireNamespace("bibtex")) {
    file <- system.file("Bib", "biblatexExamples.bib", package = "RefManageR")
    BibOptions(check.entries = FALSE)
    bib <- ReadBib(file)
    Citet(bib, 12)
    NoCite(bib, title = "Alkanethiolate")
    PrintBibliography(bib, .opts = list(style = "latex",
                      bib.style = "authoryear"))
}
\dontrun{
  if (requireNamespace("bibtex")){
    Citep(bib, c("loh", "geer"), .opts = list(cite.style = "numeric"),
          before = "see e.g., ")
    Citet(bib, "loh", .opts = list(cite.style = "numeric", super = TRUE))
    AutoCite(bib, eprinttype = "arxiv", .opts = list(cite.style = "authoryear"))
    AutoCite(bib, eprinttype = "arxiv", .opts = list(cite.style = "pandoc"))
    Citep(bib, author = "kant")
    ## shorthand field in both entries gets used for numeric and alphabetic labels
    TextCite(bib, author = "kant", .opts = list(cite.style = "alphabetic"))
    TextCite(bib, author = "kant", .opts = list(cite.style = "numeric"))
    TextCite(bib, author = "kant", .opts = list(cite.style = "alphabetic",
             style = "html"))
    punct <- unlist(BibOptions("bibpunct"))
    punct[3:4] <- c("(", ")")
    TextCite(bib, 33, .opts = list(bibpunct = punct, cite.style = "alphabetic"))

    BibOptions(restore.defaults = TRUE)
  }
}
\dontrun{
library(knitr)
## See also TestNumeric.Rmd and TestAlphabetic.Rmd for more examples
old.dir <- setwd(tdir <- tempdir())
doc <- system.file("Rmd", "TestRmd.Rmd", package = "RefManageR")
file.show(doc)
tmpfile <- tempfile(fileext = ".html", tmpdir = tdir)
knit2html(doc, tmpfile)
browseURL(tmpfile)

doc <- system.file("Rhtml", "TestAuthorYear.Rhtml", package = "RefManageR")
file.show(doc)
tmpfile <- tempfile(fileext = ".html", tmpdir = tdir)
knit2html(doc, tmpfile)
browseURL(tmpfile)
setwd(old.dir)
unlink(tdir)
}
}
\seealso{
\code{\link{print.BibEntry}}, \code{\link{BibOptions}},
\code{\link[utils]{citeNatbib}}, the package vignettes
bib <-
}
\keyword{methods}
\keyword{print}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toBiblatex.R, R/toBibtex.R
\name{toBiblatex}
\alias{toBiblatex}
\alias{toBibtex.BibEntry}
\alias{toBibtex}
\title{Convert BibEntry objects to BibTeX or BibLaTeX}
\usage{
toBiblatex(object, ...)

\method{toBibtex}{BibEntry}(
  object,
  note.replace.field = c("urldate", "pubsate", "addendum"),
  extra.fields = NULL,
  ...
)
}
\arguments{
\item{object}{an object of class BibEntry to be converted}

\item{...}{ignored}

\item{note.replace.field}{a character vector of BibLaTeX fields.  When converting an entry to BibTeX, the first field in the
entry that matches one specified in this vector will be added to the note field, \emph{if} the note field is not already
present}

\item{extra.fields}{character vector; fields that are not supported in standard BibTeX styles are by default dropped
in the result return by the toBibtex function.
Any fields specified in extra.fields will \emph{not} be dropped if present in an entry.}
}
\value{
an object of class \dQuote{Bibtex} - character vectors where each element holds one line of a BibTeX or BibLaTeX file
}
\description{
toBiblatex converts a BibEntry object to character vectors with BibLaTeX markup.  toBibtex will convert a BibEntry object
to character vectors with BibTeX markup, converting some BibLaTeX fields and all entry types that are not supported
by BibTeX to ones that are supported.
}
\details{
toBiblatex converts the BibEntry object to a vector containing the corresponding BibLaTeX file, it ensures the name
list fields (e.g. author and editor) are formatted properly to be read by bibtex and biber and otherwise prints all fields
as is, thus it is similar to \code{\link{toBibtex}}.

toBibtex will attempt to convert BibLaTeX entries to a format that can be read by bibtex.  Any fields not supported by
bibtex are dropped unless they are specified in \code{extra.fields}.  The fields below, if they are present, are converted
as described and added to a bibtex supported field, unless that field is already present.
\itemize{
\item date - The \code{date} field, if present will be truncated
to a year and added to the \code{year} field, if it is not already present. If a month is specified with the date, it will
be added to the \code{month} field.
\item journaltitle - Will be changed to journal, if it is not already present
\item location - Will be changed to address
\item institution - Converted to \code{school} for thesis entries
\item sortkey - Converted to \code{key}
\item maintitle - Converted to \code{series}
\item issuetitle - Converted to \code{booktitle}
\item eventtitle - Converted to \code{booktitle}
\item eprinttype - Converted to \code{archiveprefix} (for arXiv references)
\item eprintclass - Converted to \code{primaryclass} (for arXiv references)
}

If no \code{note} field is present, the note.replace.field can be used to specified BibLaTeX fields that can be looked for
and added to the note field if they are present.

BibLaTeX entry types that are not supported by bibtex are converted by toBibtex as follows
"mvbook" = "Book", "bookinbook" = "InBook", "suppbook" = "InBook",
\itemize{
\item MvBook,Collection,MvCollection,Reference,MvReference,Proceedings,MvProceedings,Periodical - to Book
\item BookInBook,SuppBook,InReference,SuppPeriodical - to InBook
\item report,patent - to TechReport
\item SuppCollection - to InCollection
\item thesis - to MastersThesis if \code{type = mathesis}, else to PhdThesis
\item \emph{rest} - to Misc
}
}
\examples{
if (requireNamespace("bibtex")) {
    file.name <- system.file("Bib", "biblatexExamples.bib", package="RefManageR")
    bib <- suppressMessages(ReadBib(file.name))
    toBiblatex(bib[70:72])
    toBibtex(bib[70:72])
}
}
\seealso{
\code{\link{toBibtex}}, \code{\link{BibEntry}}, \code{\link{print.BibEntry}}
}
\author{
McLean, M. W. \email{mathew.w.mclean@gmail.com}
}
\keyword{IO}
\keyword{database}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/head.BibEntry.R
\name{head.BibEntry}
\alias{head.BibEntry}
\alias{tail.BibEntry}
\title{Return the first or last part of a BibEntry object}
\usage{
\method{head}{BibEntry}(x, n = 6L, suppress.messages = TRUE, ...)

\method{tail}{BibEntry}(x, n = 6L, suppress.messages = TRUE, ...)
}
\arguments{
\item{x}{an object of class BibEntry.}

\item{n}{a single integer. If positive, size for the resulting object: number of elements for a vector
(including lists), rows for a matrix or data frame or lines for a function. If negative, all but the
n last/first number of elements of x.}

\item{suppress.messages}{boolean; should the head/tail entries be printed via \code{\link{message}}?}

\item{...}{arguments to be passed to or from other methods.}
}
\value{
an object of class BibEntry.
}
\description{
Prints the first or last entries of a BibEntry object (via \code{\link{message}}) and returns them \emph{invisibly}
  (via \code{\link{invisible}}).
}
\details{
If \code{suppress.messages} is \code{FALSE}, the head/tail entries are output to the console along
  with some additional formatting for the \sQuote{bibtype} and \sQuote{key}, in addition to
invisibly returning the entries.
}
\examples{
if (requireNamespace("bibtex")) {
    file <- system.file("Bib", "biblatexExamples.bib", package = "RefManageR")
    BibOptions(check.entries = FALSE)
    bib <- ReadBib(file)
    tail(bib, 2, suppress.messages = FALSE)
    bib <- head(bib, 1, suppress.messages = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.BibEntry}
\alias{print.BibEntry}
\title{Print BibLaTeX bibliography Entries}
\usage{
\method{print}{BibEntry}(x, .opts = list(), ...)
}
\arguments{
\item{x}{a BibEntry object}

\item{.opts}{a list of formatting options from \code{\link{BibOptions}}.  Possible
options are
\itemize{
\item \code{style} - character string naming the printing style.  Possible
values are plain text (style \dQuote{text}), BibTeX (\dQuote{Bibtex}), BibLaTeX
(\dQuote{Biblatex}), a mixture of plain text and BibTeX as
traditionally used for citations (\dQuote{citation}), HTML (\dQuote{html}),
LaTeX (\dQuote{latex}), \dQuote{markdown}, \dQuote{yaml},
R code (\dQuote{R}), and a simple copy of the textVersion elements
(style \dQuote{textVersion}, see \code{\link{BibEntry}})
\item \code{bib.style} - character string specifying BibLaTeX style to use for
formatting references.  Possible values are \dQuote{numeric} (default),
\dQuote{authoryear}, \dQuote{authortitle}, \dQuote{alphabetic}, \dQuote{draft}.
See section 3.3.2 of the BibLaTeX manual.
\item \code{sorting} - how should the entries in \code{x} be sorted?  See
\code{\link{sort.BibEntry}}.
\item \code{max.names} - maximum number of names to display for name list fields before
truncation with \dQuote{et al.}.
\item \code{first.inits} - logical; if true only initials of given names are printed,
otherwise full names are used.
\item \code{dashed} - logical; for \code{.bibstyle = "authoryear"} or
\code{.bibstyle = "authoryear"} only,
if \code{TRUE} duplicate author and editor lists are replaced with \dQuote{---}
when printed.
\item \code{no.print.fields} character vector; fields that should not be printed,
e.g., doi, url, isbn, etc.
}}

\item{...}{extra parameters to pass to the renderer.}
}
\description{
Prints bibliographic information stored in BibEntry objects in BibLaTeX style
}
\note{
setting max.names to \code{value} is equivalent to setting \code{maxnames=value} and
\code{minnames=value} in BibLaTeX.

Custom BibLaTeX styles may be defined using the function \code{\link{bibstyle}}.  To fully
support BibLaTeX, the created
environment must have functions for formatting each of the entry types described in
\code{\link{BibEntry}}.
}
\examples{
if (requireNamespace("bibtex")) {
    file.name <- system.file("Bib", "biblatexExamples.bib", package="RefManageR")
    bib <- suppressMessages(ReadBib(file.name))
    print(bib[author="aristotle"], .opts = list(bib.style = "numeric"))
    print(bib[55:57], .opts = list(bib.style = "authortitle", first.inits = FALSE))
    print(bib[80:88], .opts = list(bib.style = "alphabetic", max.names = 1,
          no.print.fields = "issn"))
    print(bib[32:36], .opts = list(bib.style = "draft"))
    oldopts <- BibOptions(bib.style = "authoryear", dashed = TRUE, sorting = "ydnt")
    bib[editor = "westfahl"]
    BibOptions(oldopts)
}
}
\references{
Lehman, Philipp and Kime, Philip and Boruvka, Audrey and Wright, J. (2013). The
biblatex Package.
\url{https://mirror.pregi.net/tex-archive/macros/latex/contrib/biblatex/doc/biblatex.pdf}.
}
\seealso{
\code{\link{BibEntry}}, \code{\link{ReadBib}}, \code{\link{sort.BibEntry}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BibEntryCombineOp.R
\name{c.BibEntry}
\alias{c.BibEntry}
\title{Combine BibEntry objects.}
\usage{
\method{c}{BibEntry}(..., recursive = FALSE)
}
\arguments{
\item{...}{- BibEntry objects to be concatenated.}

\item{recursive}{- logical; ignored.}
}
\value{
a single BibEntry object.
}
\description{
Combines mutliple BibEntry objects into a single one.
}
\note{
\code{c} will remove all attributes besides \code{class}.

No checking for duplicate entries is performed though keys will be made unique.
}
\examples{
bib <- c(BibEntry(bibtype = "article", key = "mclean2014a", title = "My New Article",
  author = "Mathew W. McLean", journaltitle = "The Journal", date = "2014-01"),
  BibEntry(bibtype = "article", key = "mclean2014b",
  title = "My Newer Article", author = "Mathew W. McLean", journaltitle = "The Journal",
  date = "2014-02"))
}
\seealso{
Other operators: 
\code{\link{$.BibEntry}()},
\code{\link{$<-.BibEntry}()},
\code{\link{+.BibEntry}()},
\code{\link{[.BibEntry}()},
\code{\link{[<-.BibEntry}()},
\code{\link{[[.BibEntry}()},
\code{\link{[[<-.BibEntry}()}
}
\concept{operators}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BibEntryExtractOp.R
\name{FindBibEntry}
\alias{FindBibEntry}
\title{Find a search term in the specified field of a BibEntry object}
\usage{
FindBibEntry(bib, term, field)
}
\description{
Workhorse function for SearchBib
}
\keyword{internal}
\newcommand{\textquotesingle}{'}
\newcommand{\hypen}{-}
\newcommand{\ast}{*}
\newcommand{\mkbibquote}{\dQuote{#1}}
\newcommand{\mkbibemph}{\emph{#1}}
\newcommand{\mkbibbold}{\bold{#1}}
\newcommand{\'I}{\u00cd}
\newcommand{\'i}{\u00ed}
\newcommand{\"I}{\u00cf}
\newcommand{\"i}{\u00ef}
\newcommand{\^I}{\u00ce}
\newcommand{\^i}{\u00ee}
\newcommand{\`I}{\u00cc}
\newcommand{\`i}{\u00ec}
\newcommand{\`\i}{\u00ec}
\newcommand{\a`}{\u00e0}
\newcommand{\a'}{\u00e1}
