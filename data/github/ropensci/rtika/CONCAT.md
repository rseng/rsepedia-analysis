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

# rtika

***Extract text or metadata from over a thousand file types.***

[![R-CMD-check](https://github.com/ropensci/rtika/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rtika/actions/)
[![ROpenSci](https://badges.ropensci.org/191_status.svg)](https://github.com/ropensci/software-review/issues/191/)
[![Coverage
status](https://codecov.io/gh/ropensci/rtika/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/rtika?branch=master)
[![Cranlogs
Downloads](https://cranlogs.r-pkg.org/badges/rtika)](https://CRAN.R-project.org/package=rtika)

> Apache Tika is a content detection and analysis framework, written in
> Java, stewarded at the Apache Software Foundation. It detects and
> extracts metadata and text from over a thousand different file types,
> and as well as providing a Java library, has server and command-line
> editions suitable for use from other programming languages …

> For most of the more common and popular formats, Tika then provides
> content extraction, metadata extraction and language identification
> capabilities. (From <https://en.wikipedia.org/wiki/Apache_Tika>,
> accessed Jan 18, 2018)

This is an R interface to the Tika software.

## Installation

To start, you need R and `Java 8` or `OpenJDK 1.8`. Higher versions
work. To check your version, run the command `java -version` from a
terminal. Get Java installation tips at
<https://www.java.com/en/download/> or
<http://openjdk.java.net/install/>. Because the `rJava` package is
***not*** required, installation is simple. You can cut and paste the
following snippet:

``` r
install.packages('rtika', repos = 'https://cloud.r-project.org')

library('rtika')

# You need to install the Apache Tika .jar once.
install_tika()
```

Read an [introductory
article](https://docs.ropensci.org/rtika/articles/rtika_introduction.html)
at <https://docs.ropensci.org/rtika/articles/rtika_introduction.html>.

## Key Features

-   `tika_text()` to extract plain text.
-   `tika_xml()` and `tika_html()` to get a structured XHMTL rendition.
-   `tika_json()` to get metadata as `.json`, with XHMTL content.
-   `tika_json_text()` to get metadata as `.json`, with plain text
    content.
-   `tika()` is the main function the others above inherit from.
-   `tika_fetch()` to download files with a file extension matching the
    Content-Type.

## Supported File Types

Tika parses and extracts text or metadata from over one thousand digital
formats, including:

-   Portable Document Format (`.pdf`)
-   Microsoft Office document formats (Word, PowerPoint, Excel, etc.)
-   Rich Text Format (`.rtf`)
-   Electronic Publication Format (`.epub`)
-   Image formats (`.jpeg`, `.png`, etc.)
-   Mail formats (`.mbox`, Outlook)
-   HyperText Markup Language (`.html`)
-   XML and derived formats (`.xml`, etc.)
-   Compression and packaging formats (`.gzip`, `.rar`, etc.)
-   OpenDocument Format
-   iWorks document formats
-   WordPerfect document formats
-   Text formats
-   Feed and Syndication formats
-   Help formats
-   Audio formats
-   Video formats
-   Java class files and archives
-   Source code
-   CAD formats
-   Font formats
-   Scientific formats
-   Executable programs and libraries
-   Crypto formats

For a list of MIME types, look for the “Supported Formats” page here:
<https://tika.apache.org/>

## Get Plain Text

**The `rtika` package processes batches of documents efficiently**, so I
recommend batches. Currently, the `tika()` parsers take a tiny bit of
time to spin up, and that will get annoying with hundreds of separate
calls to the functions.

``` r
# Test files
batch <- c(
  system.file("extdata", "jsonlite.pdf", package = "rtika"),
  system.file("extdata", "curl.pdf", package = "rtika"),
  system.file("extdata", "table.docx", package = "rtika"),
  system.file("extdata", "xml2.pdf", package = "rtika"),
  system.file("extdata", "R-FAQ.html", package = "rtika"),
  system.file("extdata", "calculator.jpg", package = "rtika"),
  system.file("extdata", "tika.apache.org.zip", package = "rtika")
)

# batches are best, and can also be piped with magrittr.
text <- tika_text(batch)

# text has one string for each document:
length(text)
#> [1] 7

# A snippet:
cat(substr(text[1], 54, 190)) 
#> ckage ‘jsonlite’
#> June 1, 2017
#> 
#> Version 1.5
#> 
#> Title A Robust, High Performance JSON Parser and Generator for R
#> 
#> License MIT + file LICENSE
```

To learn more and find out how to extract structured text and metadata,
read the vignette:
<https://docs.ropensci.org/rtika/articles/rtika_introduction.html>.

## Enhancements

Tika also can interact with the Tesseract OCR program on some Linux
variants, to extract plain text from images of text. If `tesseract-ocr`
is installed, Tika should automatically locate and use it for images and
PDFs that contain images of text. However, this does not seem to work on
OS X or Windows. To try on Linux, first follow the [Tesseract
installation
instructions](https://github.com/tesseract-ocr/tesseract/wiki). The next
time Tika is run, it should work. For a different approach, I suggest
[`tesseract`](https://github.com/ropensci/tesseract) package by @jeroen,
which is a specialized R interface.

The Apache Tika community welcomes your feedback. Issues regarding the R
interface should be raised at the [`rTika` Github Issue
Tracker](https://github.com/ropensci/tesseract). If you are confident
the issue concerns Tika or one of its underlying parsers, use the [Tika
Bugtracking
System](https://issues.apache.org/jira/projects/TIKA/issues).

## Using the Tika App Directly

If your project or package needs to use the Tika App `.jar`, you can
include `rTika` as a dependency and call the `rtika::tika_jar()`
function to get the path to the Tika app installed on the system.

## Similar R Packages

The are a number of specialized parsers that overlap in functionality.
For example, the [`pdftools`](https://github.com/ropensci/pdftools)
package extracts metadata and text from PDF files, the
[`antiword`](https://github.com/ropensci/antiword) package extracts text
from recent versions of Word, and the
[`epubr`](https://github.com/ropensci/epubr) package by @leonawicz
processes `epub` files. These packages do not depend on Java, while
`rTika` does.

The big difference between Tika and a specialized parser is that Tika
integrates dozens of specialist libraries maintained by the Apache
Foundation. Apache Tika processes over a thousand file types and
multiple versions of each. This eases the processing of digital archives
that contain unpredictable files. For example, researchers use Tika to
process archives from court cases, governments, or the Internet Archive
that span multiple years. These archives frequently contain diverse
formats and multiple versions of each format. Because Tika finds the
matching parser for each individual file, is well suited to diverse sets
of documents. In general, the parsing quality is good and consistently
so. In contrast, specialized parsers may only work with a particular
version of a file, or require extra tinkering.

On the other hand, a specialized library can offer more control and
features when it comes to structured data and formatting. For example,
the [`tabulizer`](https://github.com/ropensci/tabulizer) package by
@leeper and @tpaskhalis includes bindings to the ‘Tabula PDF Table
Extractor Library’. Because PDF files store tables as a series of
positions with no obvious boundaries between data cells, extracting a
`data.frame` or `matrix` requires heuristics and customization which
that package provides. To be fair to Tika, there are some formats where
`rtika` will extract data as table-like XML. For example, with Word and
Excel documents, Tika extracts simple tables as XHTML data that can be
turned into a tabular `data.frame` using the `rvest::html_table()`
function.

## History

In September 2017, github.com user *kyusque* released `tikaR`, which
uses the `rJava` package to interact with Tika (See:
<https://github.com/kyusque/tikaR>). As of writing, it provided similar
text and metadata extraction, but only `xml` output.

Back in March 2012, I started a similar project to interface with Apache
Tika. My code also used low-level functions from the `rJava` package. I
halted development after discovering that the Tika command line
interface (CLI) was easier to use. My empty repository is at
<https://r-forge.r-project.org/projects/r-tika/>.

I chose to finally develop this package after getting excited by Tika’s
new ‘batch processor’ module, written in Java. The batch processor has
very good efficiency when processing tens of thousands of documents.
Further, it is not too slow for a single document either, and handles
errors gracefully. Connecting `R` to the Tika batch processor turned out
to be relatively simple, because the `R` code is simple. It uses the CLI
to point Tika to the files. Simplicity, along with continuous testing,
should ease integration. I anticipate that some researchers will need
plain text output, while others will want `json` output. Some will want
multiple processing threads to speed things up. These features are now
implemented in `rtika`, although apparently not in `tikaR` yet.

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/rtika/blob/master/CONDUCT.md). By
participating in this project you agree to abide by its terms.

[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
rtika 2.0.0 (2021-08-05)
========================= 

    * Updated Tika to 2.0.0. Details are found at https://tika.apache.org/2.0.0/index.html .

rtika 1.23 (2020-04-24)
========================= 

### MINOR IMPROVEMENTS

    * Updated Tika to 1.24.1. Details are found at https://tika.apache.org/1.24.1/index.html .

rtika 1.23 (2019-12-12)
========================= 

### MINOR IMPROVEMENTS

    * Updated Tika to 1.23. Details are found at https://tika.apache.org/1.23/index.html .
    
rtika 1.22 (2019-08-01)
========================= 

### MINOR IMPROVEMENTS

    * Updated Tika to 1.22
    
rtika 1.20 (2019-02-26)
========================= 

### MINOR IMPROVEMENTS

    * Updated Tika to 1.20
    * Includes two config files to either turn on or off OCR. This is only relevant on Linux variants that have the Tesseract OCR engine installed.
    
### BUG FIX

    * Created a workaround because normalizePath() on Windows produced inconsistent results. 

rtika 1.19.1 (2018-07-08)
========================= 

### MINOR IMPROVEMENTS

    * Updated Tika to 1.19.1.
    * Updated the 'sys' package integration, and Jeroen informed the Tika team of an unexpected method of interprocess communication in the batch processor.


 rtika 1.1.19 (2018-07-08)
========================= 

### NEW FEATURES

  * The new java() function is used get the command to invoke Java for all tika() functions, and allows the option of changing its value across sessions. If you want to use a particular installation of Java, set the JAVA_HOME variable using the Sys.setenv(JAVA_HOME = 'my path'). The java() function will check for this variable, and if found return it instead of the default 'java' invocation. 
  * Updated to Tika version 1.19. 

### MINOR IMPROVEMENTS

  * The tika_check function now uses the more advance SHA512 checksum instead of the MD5. To implement this, the 'digest' package is now a dependency.

 rtika 0.1.8 (2018-04-25)
========================= 

### MINOR IMPROVEMENTS

  * The install_tika() function now gets the Tika 1.18 release that came out 2018-04-24, instead of the 1.18 development version.
  
  rtika 0.1.7 (2018-03-08)
========================= 

### MINOR IMPROVEMENTS

  * The new install_tika() function allows this package to be distributed on CRAN. The Tika App jar was too large to go on CRAN directly. The .jar is installed in the directory determined by the rappdirs::user_data_dir() function. 
  * The .onLoad() function now gives various installation advice when starting up. 
  

### DEPRECATED AND DEFUNCT
  * Removed the .jar in favor of the install_tika() function.

 rtika 0.1.6 (2018-03-01)
========================= 

### NEW FEATURES

  *  tika(), tika_xml(), tika_json(), tika_text(), and tika_html() have a new downloader, which preserve the server's content-type encoding as a file extension when possible. This should help Tika identify and parse downloaded files more reliably. It depends on the 'curl' package.
  * Added tika_fetch(), which is a stand alone function to download files and append a file extension matching the content type declared by the server. Additional features for this function include specifying the number of download retries. The output of tika_fetch() can be piped directly into other tika functions.
  * New introductory vignette covers how to use the functions and surveys several applications.
  * tika(), tika_xml(), tika_json(), tika_text(), and tika_html() can now be set to return=FALSE, which does not return any R character vector but invisibly returns NULL. This would be most useful in massive file conversion jobs with hundreds of thousands of files.
  * Used pkgdown to create a website for github pages.
  * New tika_json_text() function gets metadata in .json with plain text content.

### DEPRECATED AND DEFUNCT

  * Previous vignette has been removed in favor of new one.
  * The tikajar package is not required in this version.  Moved the .jar file back into this package to ease installation until I hear from CRAN. 
  
rtika 0.1.5 (2018-02-15)
=========================

### MINOR IMPROVEMENTS

  * Added dependency on 'sys' package because the 'system2' function was causing intermittent errors by ending tika in mid process.
  * Added startup check of the java version, using .onLoad() call to 'java -version'
  * Removed redundant conversion to UTF-8, because the Tika batch routine is already outputting UTF-8. 
  * Increased the speed of building packages (fewer downloads needed for testing, and the examples do not run).
  * Added Code of Conduct to CONDUCT.md file
  * Set default 'cleanup' attribute to TRUE.
  
rtika 0.1.4 (2018-02-15)
=========================

### MINOR IMPROVEMENTS

  * Because it is too big for CRAN, removed the Tika .jar file.
  * Added the Tika .jar to a new tikajar package on github.
  * Put the ropensci review badge on the tikajar package also, since its an essential component of this package.
  * Updated DESCRIPTION, documentation and .travis.yml to reflect the new installation routine.

rtika 0.1.3 (2018-02-04)
=========================

### NEW FEATURES

  * added convenience functions that advertise output format: tika_xml(), tika_json(), tika_text(), tika_html().
  
### MINOR IMPROVEMENTS 

  * README examples use magrittr pipe.
  
rtika 0.1.2 (2018-01-30)
=========================

### NEW FEATURES

  * added many tests to increase code coverage dramatically.
  * integrated the covr package.
  
### MINOR IMPROVEMENTS 

  * for Windows users, the curl package is recommended to prevent base R download.file from corrupting files.

### DEPRECATED AND DEFUNCT

  * removed the n_chars parameter in favor of using  file.size() internally.
  * removed onLoad.R that checked java version in order to speed package loading and simplify code coverage testing.
  
rtika 0.1.1 (2018-01-23)
=========================

### NEW FEATURES

  * allows the user to input the URLs and file paths of documents. URLs will be downloaded first to a temporary directory. The previous interface has been changed.
  
### MINOR IMPROVEMENTS

  * the Tika License file is included in the source.
  * added a vignette on basic text processing using the library.

rtika 0.1.0 (2018-01-19)
=========================

### NEW FEATURES

  * Initial release.
  * R interface to Apache Tika batch processing CLI, found to be the most efficient CLI option.
  * tika function returns processing results as a character vector.
  * includes the Tika App .jar. Tika source is available at: https://github.com/apache/tika






---
output: github_document
---
```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
if(!requireNamespace('rtika')){
    install.packages('rtika', 
        repos = 'https://cloud.r-project.org')
}

library('rtika')

if(is.na(tika_jar())){
 install_tika()
}
```

# rtika
***Extract text or metadata from over a thousand file types.***

[![R-CMD-check](https://github.com/ropensci/rtika/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rtika/actions/)
[![ROpenSci](https://badges.ropensci.org/191_status.svg)](https://github.com/ropensci/software-review/issues/191/)
[![Coverage status](https://codecov.io/gh/ropensci/rtika/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/rtika?branch=master)
[![Cranlogs Downloads](https://cranlogs.r-pkg.org/badges/rtika)](https://CRAN.R-project.org/package=rtika)

>Apache Tika is a content detection and analysis framework, written in Java, stewarded at the Apache Software Foundation. It detects and extracts metadata and text from over a thousand different file types, and as well as providing a Java library, has server and command-line editions suitable for use from other programming languages ...

>For most of the more common and popular formats, Tika then provides content extraction, metadata extraction and language identification capabilities. (From <https://en.wikipedia.org/wiki/Apache_Tika>, accessed Jan 18, 2018)

This is an R interface to the Tika software. 

## Installation
To start, you need R and `Java 8` or `OpenJDK 1.8`. Higher versions work. To check your version, run the command `java -version` from a terminal. 
Get Java installation tips at <https://www.java.com/en/download/> or <http://openjdk.java.net/install/>. 
Because the `rJava` package is ***not*** required, installation is simple. You can cut and paste the following snippet:

```{r, eval=FALSE}
install.packages('rtika', repos = 'https://cloud.r-project.org')

library('rtika')

# You need to install the Apache Tika .jar once.
install_tika()

```

Read an [introductory article](https://docs.ropensci.org/rtika/articles/rtika_introduction.html) at https://docs.ropensci.org/rtika/articles/rtika_introduction.html.

## Key Features

* `tika_text()` to extract plain text.
* `tika_xml()` and `tika_html()` to get a structured XHMTL rendition.
* `tika_json()` to get metadata as `.json`, with XHMTL content. 
* `tika_json_text()` to get metadata as `.json`, with plain text content.
* `tika()` is the main function the others above inherit from. 
* `tika_fetch()` to download files with a file extension matching the Content-Type.

## Supported File Types

Tika parses and extracts text or metadata from over one thousand digital formats, including:

* Portable Document Format (`.pdf`)
* Microsoft Office document formats (Word, PowerPoint, Excel, etc.)
* Rich Text Format (`.rtf`)
* Electronic Publication Format (`.epub`)
* Image formats (`.jpeg`, `.png`, etc.)
* Mail formats (`.mbox`, Outlook)
* HyperText Markup Language (`.html`)
* XML and derived formats (`.xml`, etc.)
* Compression and packaging formats (`.gzip`, `.rar`, etc.)
* OpenDocument Format
* iWorks document formats
* WordPerfect document formats
* Text formats
* Feed and Syndication formats
* Help formats
* Audio formats
* Video formats
* Java class files and archives
* Source code
* CAD formats
* Font formats
* Scientific formats
* Executable programs and libraries
* Crypto formats


For a list of MIME types, look for the "Supported Formats" page here: https://tika.apache.org/

## Get Plain Text

**The `rtika` package processes batches of documents efficiently**, so I recommend batches. 
Currently, the `tika()` parsers take a tiny bit of time to spin up, and that will get annoying with hundreds of separate calls to the functions.
```{r}
# Test files
batch <- c(
  system.file("extdata", "jsonlite.pdf", package = "rtika"),
  system.file("extdata", "curl.pdf", package = "rtika"),
  system.file("extdata", "table.docx", package = "rtika"),
  system.file("extdata", "xml2.pdf", package = "rtika"),
  system.file("extdata", "R-FAQ.html", package = "rtika"),
  system.file("extdata", "calculator.jpg", package = "rtika"),
  system.file("extdata", "tika.apache.org.zip", package = "rtika")
)

# batches are best, and can also be piped with magrittr.
text <- tika_text(batch)

# text has one string for each document:
length(text)

# A snippet:
cat(substr(text[1], 54, 190)) 

```

To learn more and find out how to extract structured text and metadata, read the vignette: https://docs.ropensci.org/rtika/articles/rtika_introduction.html.

## Enhancements

Tika also can interact with the Tesseract OCR program on some Linux variants, to extract plain text from images of text. If `tesseract-ocr` is installed, Tika should automatically locate and use it for images and PDFs that contain images of text. However, this does not seem to work on OS X or Windows. To try on Linux, first follow the [Tesseract installation instructions](https://github.com/tesseract-ocr/tesseract/wiki). The next time Tika is run, it should work. For a different approach, I suggest [`tesseract`](https://github.com/ropensci/tesseract) package by @jeroen, which is a specialized R interface.

The Apache Tika community welcomes your feedback. Issues regarding the R interface should be raised at the [`rTika` Github Issue Tracker](https://github.com/ropensci/tesseract). If you are confident the issue concerns Tika or one of its underlying parsers, use the [Tika Bugtracking System](https://issues.apache.org/jira/projects/TIKA/issues). 

## Using the Tika App Directly

If your project or package needs to use the Tika App `.jar`, you can include `rTika` as a dependency and call the `rtika::tika_jar()` function to get the path to the Tika app installed on the system. 


## Similar R Packages

The are a number of specialized parsers that overlap in functionality. For example, the [`pdftools`](https://github.com/ropensci/pdftools) package extracts metadata and text from PDF files, the [`antiword`](https://github.com/ropensci/antiword) package extracts text from recent versions of Word, and the [`epubr`](https://github.com/ropensci/epubr) package by @leonawicz processes `epub` files. These packages do not depend on Java, while `rTika` does. 

The big difference between Tika and a specialized parser is that Tika integrates dozens of specialist libraries maintained by the Apache Foundation.  Apache Tika processes over a thousand file types and multiple versions of each. This eases the processing of digital archives that contain unpredictable files. For example, researchers use Tika to process archives from court cases, governments, or the Internet Archive that span multiple years. These archives frequently contain diverse formats and multiple versions of each format. Because Tika finds the matching parser for each individual file, is well suited to diverse sets of documents. In general, the parsing quality is good and consistently so. In contrast, specialized parsers may only work with a particular version of a file, or require extra tinkering. 

On the other hand, a specialized library can offer more control and features when it comes to structured data and formatting. For example, the [`tabulizer`](https://github.com/ropensci/tabulizer) package by @leeper and @tpaskhalis includes bindings to the 'Tabula PDF Table Extractor Library'. Because PDF files store tables as a series of positions with no obvious boundaries between data cells, extracting a `data.frame` or `matrix` requires heuristics and customization which that package provides. To be fair to Tika, there are some formats where `rtika` will extract data as table-like XML. For example, with Word and Excel documents, Tika extracts simple tables as XHTML data that can be turned into a tabular `data.frame` using the `rvest::html_table()` function. 


## History

In September 2017, github.com user *kyusque* released `tikaR`, which uses the `rJava` package to interact with Tika (See: <https://github.com/kyusque/tikaR>). 
As of writing, it provided similar text and metadata extraction, but only `xml` output. 

Back in March 2012, I started a similar project to interface with Apache Tika.
My code also used low-level functions from the `rJava` package.
I halted development after discovering that the Tika command line interface (CLI) was easier to use.
My empty repository is at <https://r-forge.r-project.org/projects/r-tika/>.

I chose to finally develop this package after getting excited by Tika's new 'batch processor' module, written in Java.
The batch processor has very good efficiency when processing tens of thousands of documents. 
Further, it is not too slow for a single document either, and handles errors gracefully. 
Connecting `R` to the Tika batch processor turned out to be relatively simple, because the `R` code is simple.
It uses the CLI to point Tika to the files. Simplicity, along with continuous testing, should ease integration.
I anticipate that some researchers will need plain text output, while others will want `json` output.
Some will want multiple processing threads to speed things up.
These features are now implemented in `rtika`, although apparently not in `tikaR` yet.


## Code of Conduct
Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/rtika/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "Introduction to rtika"
author: "Sasha Goodman"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to rtika}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: bibliography.bib
---

```{r setup, include = FALSE}
# only evaluate code if "NOT_CRAN"
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

if(NOT_CRAN){
  if(is.na(rtika::tika_jar())){ rtika::install_tika() }  
}

```
# A Digital Babel Fish

```
                    .----.      
           ____    __\\\\\\__                 
           \___'--"          .-.          
           /___>------rtika  '0'          
          /____,--.        \----B      
                  "".____  ___-"    
                  //    / /                   
                        ]/               
```
Apache Tika is similar to the Babel fish in Douglas Adam's book, "The Hitchhikers' Guide to the Galaxy" [@mattmann2011tika p. 3]. The Babel fish translates any natural language to any other. While Apache Tika does not yet translate natural languages, it starts to tame the tower of babel of digital document formats. As the Babel fish allowed a person to understand Vogon poetry, Tika allows a computer to extract text and objects from Microsoft Word.  

The world of digital file formats is like a place where each community has their own language. Academic, business, government, and online communities use anywhere from a few file types to thousands. Unfortunately, attempts to unify groups around a single format are often fruitless [@mattmann2013computing]. 

This plethora of document formats has become a common concern. Tika is a common library to address this issue. Starting in Apache Nutch in 2005, Tika became its own project in 2007 and then a component of other Apache projects including Lucene, Jackrabbit, Mahout, and Solr [@mattmann2011tika p. 17]. 

With the increased volume of data in digital archives, and terabyte sized data becoming common, Tika's design goals include keeping complexity at bay, low memory consumption, and fast processing [@mattmann2011tika p. 18].  The `rtika` package is an interface to Apache Tika that leverages Tika's batch processor module to parse many documents fairly efficiently. Therefore, I recommend using batches whenever possible. 

# Extract Plain Text 

Video, sound and images are important, and yet much meaningful data remains numeric or textual. Tika can parse many formats and extract alpha-numeric characters, along with a few characters to control the arrangement of text, like line breaks. 

I recommend an analyst start with a directory on the computer and get a vector of paths to each file using `base::list.files()`. The commented code below has a recipe. Here, I use test files that are included with the package.

```{r, eval=NOT_CRAN} 

library('rtika')
library('magrittr')

# Code to get ALL the files in my_path:

# my_path <- "~"
# batch <- file.path(my_path,
#                 list.files(path = my_path,
#                 recursive = TRUE))

# pipe the batch into tika_text() 
# to get plain text

# test files
batch <- c(
  system.file("extdata", "jsonlite.pdf", package = "rtika"),
  system.file("extdata", "curl.pdf", package = "rtika"),
  system.file("extdata", "table.docx", package = "rtika"),
  system.file("extdata", "xml2.pdf", package = "rtika"),
  system.file("extdata", "R-FAQ.html", package = "rtika"),
  system.file("extdata", "calculator.jpg", package = "rtika"),
  system.file("extdata", "tika.apache.org.zip", package = "rtika")
)

text <-  
    batch %>%
    tika_text() 

# normal syntax also works:
# text <- tika_text(batch)

```

The output is a R character vector of the same length and order as the input files.

In the example above, there are several seconds of overhead to start up the Tika batch processor and then process the output. The most costly file was the first one. Large batches are parsed more quickly. For example, when parsing thousands of 1-5 page Word documents, I've measured 1/100th of a second per document on average.

Occasionally, files are not parsable and the returned value for the file will be `NA`. The reasons include corrupt files, disk input/output issues, empty files, password protection, a unhandled format, the document structure is broken, or the document has an unexpected variation. 

These issues should be rare. Tika works well on most documents, but if an archive is very large there may be a small percentage of unparsable files, and you might want to handle those.
```{r, eval=NOT_CRAN}
# Find which files had an issue
# Handle them if needed
batch[which(is.na(text))]
```

Plain text is easy to search using `base::grep()`.

```{r, eval=NOT_CRAN}
length(text)

search <-
    text[grep(pattern = ' is ', x = text)]

length(search)
```

With plain text, a variety of interesting analyses are possible, ranging from word counting to constructing matrices for deep learning. Much of this text processing is handled easily with the well documented `tidytext` package [@silge2017text]. Among other things, it handles tokenization and creating term-document matrices.

# Preserve Content-Type when Downloading 

A general suggestion is to use `tika_fetch()` when downloading files from the Internet, to preserve the server Content-Type information in a file extension. 

Tika's Content-Type detection is improved with file extensions (Tika also relies on other features such as Magic bytes, which are unique control bytes in the file header). The `tika_fetch()` function tries to preserves Content-Type information from the download server by finding the matching extension in Tika's database.

```{r, eval=NOT_CRAN}
download_directory <- tempfile('rtika_')

dir.create(download_directory)

urls <- c('https://tika.apache.org/',
          'https://cran.rstudio.com/web/packages/keras/keras.pdf')

downloaded <- 
    urls %>% 
    tika_fetch(download_directory)

# it will add the appropriate file extension to the downloads
downloaded

```
This `tika_fetch()` function is used internally by the `tika()` functions when processing URLs. By using `tika_fetch()` explicitly with a specified directory, you can also save the files and return to them later. 

# Settings for Big Datasets

Large jobs are possible with `rtika`. However, with hundreds of thousands of documents, the R object returned by the `tika()` functions can be too big for RAM. In such cases, it is good to use the computer's disk more, since running out of RAM slows the computer.

I suggest changing two parameters in any of the `tika()` parsers. First, set `return = FALSE` to prevent returning a big R character vector of text. Second, specify an existing directory on the file system using `output_dir`, pointing to where the processed files will be saved. The files can be dealt with in smaller batches later on. 

Another option is to increase the number of threads, setting `threads` to something like the number of processors minus one. 

```{r, eval=NOT_CRAN}
# create a directory not already in use.
my_directory <-
   tempfile('rtika_')
                  
dir.create(my_directory)

# pipe the batch to tika_text()
batch %>%
tika_text(threads = 4,
          return = FALSE,
          output_dir = my_directory) 

# list all the file locations 
processed_files <- file.path(
                normalizePath(my_directory),
                list.files(path = my_directory,
                recursive = TRUE)
                )

```
 The location of each file in `output_dir` follows a convention from the Apache Tika batch processor: the full path to each file mirrors the original file's path, only within the `output_dir`. 
```{r, eval=NOT_CRAN}
processed_files
```
Note that `tika_text()` produces `.txt` files, `tika_xml()` produces `.xml` files, `tika_html()` produces `.html` files, and both `tika_json()` and `tika_json_text()` produce `.json` files.


# Get a Structured XHTML Rendition
 
Plain text falls short for some purposes. For example, pagination might be important for selecting a particular page in a PDF.  The Tika authors chose HTML as a universal format because it offers semantic elements that are common or familiar. For example, the hyperlink is represented in HTML as the anchor element `<a>` with the attribute `href`. The HTML in Tika preserves this metadata:

```{r , eval=NOT_CRAN}
library('xml2')

# get XHTML text
html <- 
    batch %>%
    tika_html() %>%
    lapply(xml2::read_html)

# parse links from documents
links <-
    html %>%
    lapply(xml2::xml_find_all, '//a') %>%
    lapply(xml2::xml_attr, 'href')

sample(links[[1]],10)
```


Each type of file has different information preserved by Tika's internal parsers. The particular aspects vary. Some notes:


* PDF files retain pagination, with each page starting with the XHTML element `<div class="page">`. 
* PDFs retain hyperlinks in the anchor element `<a>` with the attribute `href`.
* Word and Excel documents retain tabular data as a `<table>` element. The `rvest` package has a function to get tables of data  with `rvest::html_table()`.
* Multiple Excel sheets are preserved as multiple XHTML tables. Ragged tables, where rows have differing numbers of cells, are not supported.


Note that `tika_html()` and `tika_xml()` both produce the same strict form of HTML called XHTML, and either works essentially the same for all the documents I've tried. 

# Access Metadata in the XHTML
The `tika_html()` and `tika_xml()` functions are focused on extracting strict, structured HTML as XHTML. In addition, metadata can be accessed in the `meta` tags of the XHTML. Common metadata fields include `Content-Type`, `Content-Length`, `Creation-Date`, and `Content-Encoding`.

```{r , eval=NOT_CRAN}
# Content-Type
html %>%
lapply(xml2::xml_find_first, '//meta[@name="Content-Type"]') %>%
lapply(xml2::xml_attr, 'content') %>%
unlist()

# Creation-Date
html %>%
lapply(xml2::xml_find_first, '//meta[@name="Creation-Date"]') %>%
lapply(xml2::xml_attr, 'content') %>%
unlist()

```



# Get Metadata in JSON Format

Metadata can also accessed with `tika_json()` and `tika_json_text()`. Consider all that can be found from a single image:


```{r, eval=NOT_CRAN}
library('jsonlite')
# batch <- system.file("extdata", "calculator.jpg", package = "rtika")

# a list of data.frames
metadata <-
    batch %>% 
    tika_json() %>%
    lapply(jsonlite::fromJSON)

# look at metadata for an image
str(metadata[[6]])

```


In addition, each specific format can have its own specialized metadata fields. For example, photos sometimes store latitude and longitude:

```{r, eval=NOT_CRAN}
metadata[[6]]$'geo:lat'
metadata[[6]]$'geo:long'
```


# Get Metadata from "Container" Documents

Some types of documents can have multiple objects within them. For example, a `.gzip` file may contain many other files. The `tika_json()` and `tika_json_text()` functions have a special ability that others do not. They will recurse into a container and examine each file within. The Tika authors call the format `jsonRecursive` for this reason.

In the following example, I created a compressed archive of the Apache Tika homepage, using the command line programs `wget` and `zip`. The small archive includes the HTML page, its images, and required files. 

```{r, eval=NOT_CRAN}
# wget gets a webpage and other files. 
# sys::exec_wait('wget', c('--page-requisites', 'https://tika.apache.org/'))
# Put it all into a .zip file 
# sys::exec_wait('zip', c('-r', 'tika.apache.org.zip' ,'tika.apache.org'))
batch <- system.file("extdata", "tika.apache.org.zip", package = "rtika")

# a list of data.frames
metadata <-
    batch %>% 
    tika_json() %>%
    lapply(jsonlite::fromJSON)

# The structure is very long. See it on your own with: str(metadata)

```

Here are some of the main metadata fields of the recursive `json` output:

```{r, eval=NOT_CRAN}
# the 'X-TIKA:embedded_resource_path' field
embedded_resource_path <- 
    metadata %>%
    lapply(function(x){ x$'X-TIKA:embedded_resource_path' }) 

embedded_resource_path
```
The `X-TIKA:embedded_resource_path` field tells you where in the document hierarchy each object resides. The first item in the character vector is the root, which is the container itself. The other items are embedded one layer down, as indicated by the forward slash `/`. In the context of the `X-TIKA:embedded_resource_path` field,  paths are not literally directory paths like in a file system. In reality, the image `icon_info_sml.gif` is within a folder called `images`. Rather, the number of forward slashes  indicates the level of recursion within the document. One slash `/` reveals a first set of embedded documents. Additional slashes `/` indicate that the parser has recursed into an embedded document within an embedded document. 

```{r, eval=NOT_CRAN}
content_type <-
    metadata %>%
    lapply(function(x){ x$'Content-Type' }) 

content_type
```
The `Content-Type` metadata reveals the first item is the container and has the type `application/zip`. The items after that are deeper and include web formats such as `application/xhtml+xml`, `image/png`, and `text/css`.

```{r, eval=NOT_CRAN}
content <- 
     metadata %>%
    lapply(function(x){ x$'X-TIKA:content' })

str(content)

```

The `X-TIKA:content` field includes the XHTML rendition of an object. It is possible to extract plain text in the `X-TIKA:content` field by calling `tika_json_text()` instead. That is the only difference between `tika_json()` and `tika_json_text()`.

It may be surprising to learn that Word documents are containers (at least the modern `.docx` variety are). By parsing them with `tika_json()` or `tika_json_text()`, the various images and embedded objects can be analyzed. However, there is an added complexity, because each document may  produce a long vector of `Content-Types` for each embedded file, instead of a single `Content-Type` for the container like `tika_xml()` and `tika_html()`.


# Extending rtika

Out of the box, `rtika` uses all the available Tika Detectors and Parsers and runs with sensible defaults. For most, this will work well.

In future versions, Tika uses a configuration file to customize parsing. This config file option is on hold in `rtika`, because Tika's batch module is still new and the config file format will likely change and be backward incompatible. Please stay tuned.

There is also room for improvement with the document formats common in the R community, especially Latex and Markdown. Tika currently reads and writes these formats just fine, captures metadata and recognizes the MIME type when downloading with `tika_fetch()`. However, Tika does not have parsers to fully understand the Latex or Markdown document structure, render it to XHTML, and extract the plain text while ignoring markup. For these cases, Pandoc will be more useful (See: https://pandoc.org/demos.html ). 

You may  find these resources useful:

* Current Tika issues and progress can be seen here: https://issues.apache.org/jira/projects/TIKA
* The Tika Wiki is here: https://cwiki.apache.org/confluence/display/tika/
* Tika sourcecode: https://github.com/apache/tika

# References

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tika_json_text.R
\name{tika_json_text}
\alias{tika_json_text}
\title{Get json Metadata and Plain Text Content}
\usage{
tika_json_text(input, ...)
}
\arguments{
\item{input}{Character vector describing the paths and/or urls to the input documents.}

\item{...}{Other parameters to be sent to \code{tika()}.}
}
\value{
A character vector in the same order and with the same length as \code{input}, of unparsed \code{json}. Unprocessed files are \code{as.character(NA)}.
}
\description{
Tika can parse and extract text from almost anything, including zip, tar, tar.bz2, and other archives that contain documents.
 If you have a zip file with 100 text files in it, you can get the text and metadata for each file nested inside of the zip file.
 This recursive output is currently used for the jsonified mode. See:  https://wiki.apache.org/tika/RecursiveMetadata
 
 The document contents are plain text in the "X-TIKA:content" field.
  
  
  If \code{output_dir} is specified, files will have the \code{.json} file extension.
}
\examples{
\donttest{
batch <- c(
 system.file("extdata", "jsonlite.pdf", package = "rtika"),
 system.file("extdata", "curl.pdf", package = "rtika"),
 system.file("extdata", "table.docx", package = "rtika"),
 system.file("extdata", "xml2.pdf", package = "rtika"),
 system.file("extdata", "R-FAQ.html", package = "rtika"),
 system.file("extdata", "calculator.jpg", package = "rtika"),
 system.file("extdata", "tika.apache.org.zip", package = "rtika")
)
json <- tika_json_text(batch)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tika_jar.R
\name{tika_jar}
\alias{tika_jar}
\title{Path to Apache Tika}
\usage{
tika_jar()
}
\value{
A string describing the file path to the Tika App \code{.jar} file. If not found, \code{NA}.
}
\description{
Gets the path to the Tika App \code{.jar} installed by \code{tika_install()}.
}
\section{Details}{

The \code{tika_jar()} function also checks if the \code{.jar} is actually on the file system.

The file path is used by all of the \code{tika()} functions by default.
}

\section{Alternative Uses}{

You can call Apache Tika directly,
as shown in the examples here.

It is better to use the \code{sys} package and avoid \code{system2()},
which has caused erratic, intermittent errors with Tika.
}

\examples{
\donttest{
jar <- tika_jar()
# see help
sys::exec_wait('java',c('-jar',jar, '--help'))
# detect language of web page
sys::exec_wait('java',c('-jar',jar, '--language','https://tika.apache.org/'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tika_html.R
\name{tika_html}
\alias{tika_html}
\title{Get Structured XHTML}
\usage{
tika_html(input, ...)
}
\arguments{
\item{input}{Character vector describing the paths and/or urls to the input documents.}

\item{...}{Other parameters to be sent to \code{tika()}.}
}
\value{
A character vector in the same order and with the same length as \code{input}, of unparsed \code{XHTML}. Unprocessed files are \code{as.character(NA)}.
}
\description{
If \code{output_dir} is specified, files will have the \code{.html} file extension.
}
\examples{
\donttest{
batch <- c(
 system.file("extdata", "jsonlite.pdf", package = "rtika"),
 system.file("extdata", "curl.pdf", package = "rtika"),
 system.file("extdata", "table.docx", package = "rtika"),
 system.file("extdata", "xml2.pdf", package = "rtika"),
 system.file("extdata", "R-FAQ.html", package = "rtika"),
 system.file("extdata", "calculator.jpg", package = "rtika"),
 system.file("extdata", "tika.apache.org.zip", package = "rtika")
)
html <- tika_html(batch)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tika_check.R
\name{tika_check}
\alias{tika_check}
\title{Check Tika against a checksum}
\usage{
tika_check(digest, jar = tika_jar(), algo = "sha512")
}
\arguments{
\item{digest}{Character vector of length one with the target checksum.}

\item{jar}{Optional alternative path to a Tika \code{jar} file.}

\item{algo}{Optional algorithm used to create checksum. Defaults to SHA512.}
}
\value{
logical if the \code{jar} checksum matches \code{digest}.
}
\description{
This is used by \code{install_tika()} internally, 
or can be called directly on a \code{jar} file.
The latest \code{jar} files and checksums are at https://tika.apache.org/download.html.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install_tika.R
\name{install_tika}
\alias{install_tika}
\title{Install or Update the Apache Tika \code{jar}}
\usage{
install_tika(
  version = "2.0.0",
  digest = paste0("8479aa41611add4f6ffa4917991ee0a74d9f1ece12a",
    "65133368567e35ef91648142279e6981f09820e1658",
    "04d051aa8406d376ed382b17cfe2f691d30e365b9f"),
  mirrors = c("https://ftp.wayne.edu/apache/tika/",
    "http://mirrors.ocf.berkeley.edu/apache/tika/", "http://apache.cs.utah.edu/tika/",
    "http://mirror.cc.columbia.edu/pub/software/apache/tika/"),
  retries = 2,
  url = character()
)
}
\arguments{
\item{version}{The declared Tika version}

\item{digest}{The sha512 checksum. Set to an empty string \code{""} to skip the check.}

\item{mirrors}{A vector of Apache mirror sites. One is picked randomly.}

\item{retries}{The number of times to try the download.}

\item{url}{Optional url to a particular location of the tika app. Setting this to any character string overrides downloading from random mirrors.}
}
\value{
Logical if the installation was successful.
}
\description{
This downloads and installs the Tika App \code{jar} (~60 MB) into a user directory,
and verifies the integrity of the file using a checksum.
The default settings should work fine.
}
\section{Details}{

The default settings of \code{install_tika()} should typically be left as they are.

This function will download the version of the Tika \code{jar} tested to work
with this package, and can verify file integrity using a checksum.

It will normally download from a random Apache mirror.
If the mirror fails,
it tries the archive at \code{http://archive.apache.org/dist/tika/}.
You can also enter a value for \code{url} directly to override this.

It will download into a directory determined
by the \code{rappdirs::user_data_dir()} function,
specific to the operating system.

If \code{tika()} is stopping with an error compalining about the \code{jar},
try running \code{install_tika()} again.
}

\section{Uninstalling}{

If you are uninstalling the entire \code{rtika} package
and want to remove the Tika App \code{jar} also,
run:

 \code{unlink(rappdirs::user_data_dir('rtika'), recursive = TRUE)}

Alternately, navigate to the install folder and delete it manually.
It is the file path returned by \code{rappdirs::user_data_dir('rtika')}.
The path is OS specific, and explained here:
https://github.com/r-lib/rappdirs .
}

\section{Distribution}{

Tika is distributed under the Apache License Version 2.0,
which generally permits distribution of the code "Object" without the "Source".
The master copy of the Apache Tika source code is held in GIT. 
You can fetch (clone) the large source from GitHub ( https://github.com/apache/tika ).
}

\examples{
\donttest{
install_tika()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/java.R
\name{java}
\alias{java}
\title{System Command to Run Java}
\usage{
java()
}
\value{
The system command needed to invoke Java, as a string.
}
\description{
Gets the system command needed to run Java from the command line, as a string.
Typically, this is the string: 'java'. 
However, if the R session has the \code{JAVA_HOME} environmental variable set, it will use that to locate java instead. 
This can be persisted over sessions (see the Details below).
}
\section{Details}{

This function is used by all of the \code{tika()} functions internally as the default value to its \code{java} parameter. 

This function tries to find an environmental variable using \code{Sys.getenv("JAVA_HOME")}.
It looks for the \code{java} executable inside the \code{bin} directory in the \code{JAVA_HOME} directory.

If you want to use a specific version of Java, set the \code{JAVA_HOME} variable using \code{Sys.setenv(JAVA_HOME = 'my path')},
where 'my path' is the path to a folder that has a \code{bin} directory with a \code{java} executable. 

For example, on Windows 10 \code{JAVA_HOME} might be \code{C:/Program Files (x86)/Java/jre1.8.0_171}. 
On Ubuntu and OS X, it might be the \code{/usr} directory.  

The \code{JAVA_HOME} variable can also be set to persist over sessions.
Add the path to the \code{.Rprofile} by adding \code{Sys.setenv(JAVA_HOME = 'my path')}, and it will use that every time R is started.
}

\examples{
\donttest{
# Typically, this function returns the string 'java'.
# If JAVA_HOME is set, it's a path to java in a 'bin' folder.
java()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tika_json.R
\name{tika_json}
\alias{tika_json}
\title{Get json Metadata and XHTML Content}
\usage{
tika_json(input, ...)
}
\arguments{
\item{input}{Character vector describing the paths and/or urls to the input documents.}

\item{...}{Other parameters to be sent to \code{tika()}.}
}
\value{
A character vector in the same order and with the same length as \code{input}, of unparsed \code{json}. Unprocessed files are \code{as.character(NA)}.
}
\description{
Tika can parse and extract text from almost anything, including zip, tar, tar.bz2, and other archives that contain documents.
 If you have a zip file with 100 text files in it, you can get the text and metadata for each file nested inside of the zip file.
 This recursive output is currently used for the jsonified mode. See:  https://wiki.apache.org/tika/RecursiveMetadata
 
 The document content is XHTML in the "X-TIKA:content" field.
  
 If \code{output_dir} is specified, files will have the \code{.json} file extension.
}
\examples{
\donttest{
batch <- c(
 system.file("extdata", "jsonlite.pdf", package = "rtika"),
 system.file("extdata", "curl.pdf", package = "rtika"),
 system.file("extdata", "table.docx", package = "rtika"),
 system.file("extdata", "xml2.pdf", package = "rtika"),
 system.file("extdata", "R-FAQ.html", package = "rtika"),
 system.file("extdata", "calculator.jpg", package = "rtika"),
 system.file("extdata", "tika.apache.org.zip", package = "rtika")
)
json <- tika_json(batch)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tika_text.R
\name{tika_text}
\alias{tika_text}
\title{Get Plain Text}
\usage{
tika_text(input, ...)
}
\arguments{
\item{input}{Character vector describing the paths and/or urls to the input documents.}

\item{...}{Other parameters to be sent to \code{tika()}.}
}
\value{
A character vector in the same order and with the same length as \code{input}, of plain text. Unprocessed files are \code{as.character(NA)}.
}
\description{
If \code{output_dir} is specified, files will have the \code{.txt} file extension.
}
\examples{
\donttest{
batch <- c(
 system.file("extdata", "jsonlite.pdf", package = "rtika"),
 system.file("extdata", "curl.pdf", package = "rtika"),
 system.file("extdata", "table.docx", package = "rtika"),
 system.file("extdata", "xml2.pdf", package = "rtika"),
 system.file("extdata", "R-FAQ.html", package = "rtika"),
 system.file("extdata", "calculator.jpg", package = "rtika"),
 system.file("extdata", "tika.apache.org.zip", package = "rtika")
)
text <- tika_text(batch)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tika_xml.R
\name{tika_xml}
\alias{tika_xml}
\title{Get a Structured XHTML Rendition}
\usage{
tika_xml(input, ...)
}
\arguments{
\item{input}{Character vector describing the paths and/or urls to the input documents.}

\item{...}{Other parameters to be sent to \code{tika()}.}
}
\value{
A character vector in the same order and with the same length as \code{input}, of unparsed \code{XHTML}. Unprocessed files are \code{as.character(NA)}.
}
\description{
If \code{output_dir} is specified, files will have the \code{.xml} file extension.
}
\examples{
\donttest{
batch <- c(
 system.file("extdata", "jsonlite.pdf", package = "rtika"),
 system.file("extdata", "curl.pdf", package = "rtika"),
 system.file("extdata", "table.docx", package = "rtika"),
 system.file("extdata", "xml2.pdf", package = "rtika"),
 system.file("extdata", "R-FAQ.html", package = "rtika"),
 system.file("extdata", "calculator.jpg", package = "rtika"),
 system.file("extdata", "tika.apache.org.zip", package = "rtika")
)
xml <- tika_xml(batch)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtika.R
\docType{package}
\name{rtika}
\alias{rtika}
\title{rtika: R Interface to 'Apache Tika'}
\description{
Extract text or metadata from over a thousand file types. Get either plain text or structured XHTML content.
}
\section{Installing}{


If you have not done so already, finish installing \pkg{rtika} by typing in the R console:

\code{install_tika()}
}

\section{Getting Started}{


The \code{\link{tika_text}} function will extract plain text from many types of documents. It is a good place to start. Please read the Vignette also.
Other main functions include \code{\link{tika_xml}} and \code{\link{tika_html}} that get a structured XHMTL rendition. The \code{\link{tika_json}} function gets metadata as `.json`, with XHMTL content. 

The \code{\link{tika_json_text}} function gets metadata as `.json`, with plain text content.

\code{\link{tika}} is the main function the others above inherit from. 
 
Use \code{\link{tika_fetch}} to download files with a file extension matching the Content-Type.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tika.R
\name{tika}
\alias{tika}
\title{Main R Interface to 'Apache Tika'}
\usage{
tika(
  input,
  output = c("text", "jsonRecursive", "xml", "html")[1],
  output_dir = "",
  return = TRUE,
  java = rtika::java(),
  jar = rtika::tika_jar(),
  threads = 2,
  max_restarts = integer(),
  timeout = 3e+05,
  max_file_size = integer(),
  config = system.file("extdata", "ocr.xml", package = "rtika"),
  args = character(),
  quiet = TRUE,
  cleanup = TRUE,
  lib.loc = .libPaths()
)
}
\arguments{
\item{input}{Character vector describing the paths to the input documents.
Strings starting with 'http://','https://', or 'ftp://' are downloaded to a
temporary directory. On Windows, the local paths cannot span
drives because of a Windows convention.}

\item{output}{Optional character vector of the output format. The default,
\code{"text"}, gets plain text without metadata. \code{"xml"} and
\code{"html"} get \code{XHTML} text with metadata. \code{"jsonRecursive"}
gets \code{XHTML} text and \code{json} metadata.
\code{c("jsonRecursive","text")} or \code{c("J","t")} get plain text and
\code{json} metadata. See the 'Output Details' section.}

\item{output_dir}{Optional directory path to save the converted files in.
Tika may overwrite files so an empty directory is best. See the 'Output
Details' section before using.}

\item{return}{Logical if an R object should be returned. Defaults to
TRUE. If set to FALSE, and output_dir (above) must be specified.}

\item{java}{Optional command to invoke Java. For example, it can be the full
path to a particular Java version. See the Configuration section below.}

\item{jar}{Optional alternative path to a \code{tika-app-X.XX.jar}. Useful
if this package becomes out of date.}

\item{threads}{Integer of the number of file consumer threads Tika uses.
Defaults to 2.}

\item{max_restarts}{Integer of the maximum number of times the watchdog
process will restart the child process. The default is no limit.}

\item{timeout}{Integer of the number of milliseconds allowed to a parse
before the process is killed and restarted. Defaults to 300000.}

\item{max_file_size}{Integer of the maximum bytes allowed.
Do not process files larger than this. The default is unlimited.}

\item{config}{Path to the XML config file. Defaults to \code{system.file("extdata", "ocr.xml", package = "rtika")}'. There is also a \code{no-ocr.xml} file available.}

\item{args}{Optional character vector of additional arguments passed to Tika,
that may not yet be implemented in this R interface, in the pattern of
\code{c('-arg1','setting1','-arg2','setting2')}.}

\item{quiet}{Logical if Tika command line messages and errors are to be
suppressed. Defaults to \code{TRUE}.}

\item{cleanup}{Logical to clean up temporary files after running the command,
which can accumulate. Defaults to \code{TRUE}. They are in \code{tempdir()}. These
files are automatically removed at the end of the R session even if set to
FALSE.}

\item{lib.loc}{Optional character vector describing the library paths.
Normally, it's best to
leave this parameter alone. The parameter is included
mainly for package testing.}
}
\value{
A character vector in the same order and with the same length as
\code{input}. Unprocessed files are \code{as.character(NA)}.
If \code{return = FALSE}, then a \code{NULL} value is invisibly returned.
See the Output Details section below.
}
\description{
Extract text or metadata from over a thousand file types.
Get either plain text or structured \code{XHTML}.
Metadata includes \code{Content-Type}, character encoding, and Exif data from
jpeg or tiff images. See the long list of supported file types, 
click the "Supported Formats" link on this page :
\url{https://tika.apache.org/}.
}
\section{Output Details}{

If an input file did not exist, could not be downloaded, was a directory, or
Tika could not process it, the result will be \code{as.character(NA)} for
that file.

By default, \code{output = "text"} and this produces plain text with no
metadata. Some formatting is preserved in this case using tabs, newlines and
spaces.

Setting \code{output} to either \code{"xml"} or the shortcut \code{"x"} will
produce a strict form of \code{HTML} known as \code{XHTML}, with metadata in
the \code{head} node and formatted text in the \code{body}.
Content retains more formatting with \code{"xml"}. For example, a Word or
Excel table will become a HTML \code{table}, with table data as text in
\code{td} elements. The \code{"html"} option and its shortcut \code{"h"}
seem to produce the same result as \code{"xml"}.
Parse XHTML output with \code{xml2::read_html}.

Setting \code{output} to \code{"jsonRecursive"} or its shortcut \code{"J"}
produces a tree structure in `json`. Metadata fields are at the top level.
The \code{XHTML} or plain text will be found in the \code{X-TIKA:content}
field. By default the text is \code{XHTML}. This can be changed to plain
text like this: \code{output=c("jsonRecursive","text")} or
\code{output=c("J","t")}. This syntax is meant to mirror Tika's. Parse
\code{json} with \code{jsonlite::fromJSON}.

 If \code{output_dir} is specified, then the converted files will also be
 saved to this directory. It's best to use an empty directory because Tika
 may overwrite existing files. Tika seems to add an extra file extension to
 each file to reduce the chance, but it's still best to use an empty
 directory. The file locations within the \code{output_dir} maintain the same
 general path structure as the input files. Downloaded files have a path
 similar to the `tempdir()` that R uses. The original paths are now relative
 to \code{output_dir}.  Files are appended with \code{.txt} for the default
 plain text, but can be \code{.json}, \code{.xml}, or \code{.html} depending
 on the \code{output} setting. One way to get a list of the processed files
 is to use \code{list.files} with \code{recursive=TRUE}.
 If \code{output_dir} is not specified, files are saved to a volatile temp
 directory named by \code{tempdir()} and will be deleted when R shuts down.
 If this function will be run on very large batches repeatedly, these
 temporary files can be cleaned up every time by adding
 \code{cleanup=TRUE}.
}

\section{Background}{

Tika is a foundational library for several Apache projects such as the Apache
Solr search engine. It has been in development since at least 2007. The most
efficient way I've found to process many thousands of documents is Tika's
'batch' mode, which is the only mode used in `rtika`. There are potentially
more things that can be done, given enough time and attention, because
Apache Tika includes many libraries and methods in its .jar file. The source is available at:
\url{https://tika.apache.org/}.
}

\section{Installation}{

 Tika requires Java 8.

 Java installation instructions are at http://openjdk.java.net/install/
or https://www.java.com/en/download/help/download_options.xml.

By default, this R package internally invokes Java by calling the \code{java}
command from the command line. To specify the path to a particular Java
version, set the path in the \code{java} attribute of the \code{tika}
function.
}

\examples{
\donttest{
#extract text
batch <- c(
  system.file("extdata", "jsonlite.pdf", package = "rtika"),
  system.file("extdata", "curl.pdf", package = "rtika"),
  system.file("extdata", "table.docx", package = "rtika"),
  system.file("extdata", "xml2.pdf", package = "rtika"),
  system.file("extdata", "R-FAQ.html", package = "rtika"),
  system.file("extdata", "calculator.jpg", package = "rtika"),
  system.file("extdata", "tika.apache.org.zip", package = "rtika")
)
text = tika(batch)
cat(substr(text[1],45,450))

#more complex metadata
if(requireNamespace('jsonlite')){

  json = tika(batch,c('J','t'))
  # 'J' is shortcut for jsonRecursive
  # 't' for text
  metadata = lapply(json, jsonlite::fromJSON )

  #embedded resources
  lapply(metadata, function(x){ as.character(x$'Content-Type') })

  lapply(metadata, function(x){ as.character(x$'Creation-Date') })

  lapply(metadata, function(x){  as.character(x$'X-TIKA:embedded_resource_path') })
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tika_fetch.R
\name{tika_fetch}
\alias{tika_fetch}
\title{Fetch Files with the Content-Type Preserved in the File Extension}
\usage{
tika_fetch(
  urls,
  download_dir = tempdir(),
  ssl_verifypeer = TRUE,
  retries = 1,
  quiet = TRUE
)
}
\arguments{
\item{urls}{Character vector of one or more URLs to be downloaded.}

\item{download_dir}{Character vector of length one describing the path to the directory to save the results.}

\item{ssl_verifypeer}{Logical, with a default of TRUE. Some server SSL certificates might not be recognized by the host system, and in these rare cases the user can ignore that if they know why.}

\item{retries}{Integer of the number of times to retry each url after a failure to download.}

\item{quiet}{Logical if download warnings should be printed. Defaults to FALSE.}
}
\value{
Character vector of the same length and order as input with the paths describing the locations of the downloaded files. Errors are returned as NA.
}
\description{
On the Internet, Content-Type information is mainly communicated via the server's headers.
This is an issue if a file is saved to disk without examining the headers.
The file can have a missing or incorrect file extension.
For example, a URL ending in a slash (\code{/}) can produce file with the Content-Type of  \code{text/html}.
The same URL might also produce a \code{image/jpeg} or \code{application/pdf} file.
URLs ending in \code{.php}, \code{.cfm} can produce any Content-Type.
The downloaded file will lose the server's declared Content-Type unless its appended as a file extension.
\code{tika_fetch()} gets a file from the URL, examines the server headers,
and appends the matching file extension
from Tika's database.
}
\examples{
\donttest{
tika_fetch('https://tika.apache.org/')
# a unique file name with .html appended to it
}
}
