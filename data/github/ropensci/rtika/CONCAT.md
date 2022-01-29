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






