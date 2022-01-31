# pdftools

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/pdftools?branch=master&svg=true)](https://ci.appveyor.com/project/jeroen/pdftools)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/pdftools)](http://cran.r-project.org/package=pdftools)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/pdftools)](http://cran.r-project.org/web/packages/pdftools/index.html)

## Introduction

Scientific articles are typically locked away in PDF format, a format designed primarily for printing but not so great for searching or indexing. The new pdftools package allows for extracting text and metadata from pdf files in R. From the extracted plain-text one could find articles discussing a particular drug or species name, without having to rely on publishers providing metadata, or pay-walled search engines.

The pdftools slightly overlaps with the Rpoppler package by Kurt Hornik. The main motivation behind developing pdftools was that Rpoppler depends on glib, which does not work well on Mac and Windows. The pdftools package uses the poppler c++ interface together with Rcpp, which results in a lighter and more portable implementation.


## Installation

On Windows and Mac the binary packages can be installed directly from CRAN:

```r
install.packages("pdftools")
```

Installation on Linux requires the poppler development library. On __Ubuntu 16.04 (Xenial)__ and __Ubuntu 18.04 (Bionic)__ we have backports that support the latest `pdf_data()` functionality:

```
sudo add-apt-repository -y ppa:cran/poppler
sudo apt-get update
sudo apt-get install -y libpoppler-cpp-dev
```

On other versions of __Debian__ or __Ubuntu__ simply use::

```
sudo apt-get install libpoppler-cpp-dev
```

If you want to install the package from source on __MacOS__ you need brew:

```
brew install poppler
```

On Fedora:

```
sudo yum install poppler-cpp-devel
```

### Building from source

#### On Ubuntu 

__Update__: Itt is now recommended to use the backport PPA mentioned above. If you really want to build from source, follow the instructions [of this askubuntu.com answer](https://askubuntu.com/a/1112947).

#### On CentOS

On CentOS the `libpoppler-cpp` library is not included with the system so we need to build from source. Note that recent versions of poppler require C++11 which is not available on CentOS, so we build a slightly older version of libpoppler.

```sh
# Build dependencies
yum install wget xz libjpeg-devel openjpeg2-devel

# Download and extract
wget https://poppler.freedesktop.org/poppler-0.47.0.tar.xz
tar -Jxvf poppler-0.47.0.tar.xz
cd poppler-0.47.0

# Build and install
./configure
make
sudo make install
```

By default libraries get installed in `/usr/local/lib` and `/usr/local/include`. On CentOS this is not a default search path so we need to set `PKG_CONFIG_PATH` and  `LD_LIBRARY_PATH` to point R to the right directory:

```sh
export LD_LIBRARY_PATH="/usr/local/lib"
export PKG_CONFIG_PATH="/usr/local/lib/pkgconfig"
```

We can then start R and install `pdftools`.

## Getting started

The `?pdftools` manual page shows a brief overview of the main utilities. The most important function is `pdf_text` which returns a character vector of length equal to the number of pages in the pdf. Each string in the vector contains a plain text version of the text on that page.

```r
library(pdftools)
download.file("http://arxiv.org/pdf/1403.2805.pdf", "1403.2805.pdf", mode = "wb")
txt <- pdf_text("1403.2805.pdf")

# first page text
cat(txt[1])

# second page text
cat(txt[2])
```

In addition, the package has some utilities to extract other data from the PDF file. The `pdf_toc` function shows the table of contents, i.e. the section headers which pdf readers usually display in a menu on the left. It looks pretty in JSON:

```r
# Table of contents
toc <- pdf_toc("1403.2805.pdf")

# Show as JSON
jsonlite::toJSON(toc, auto_unbox = TRUE, pretty = TRUE)
```

Other functions provide information about fonts, attachments and metadata such as the author, creation date or tags.


```r
# Author, version, etc
info <- pdf_info("1403.2805.pdf")

# Table with fonts
fonts <- pdf_fonts("1403.2805.pdf")
```

## Bonus feature: rendering pdf

A bonus feature on most platforms is rendering of PDF files to bitmap arrays. The poppler library provides all functionality to implement a complete PDF reader, including graphical display of the content. In R we can use `pdf_render_page` to render a page of the PDF into a bitmap, which can be stored as e.g. png or jpeg.

```r
# renders pdf to bitmap array
bitmap <- pdf_render_page("1403.2805.pdf", page = 1)

# save bitmap image
png::writePNG(bitmap, "page.png")
jpeg::writeJPEG(bitmap, "page.jpeg")
webp::write_webp(bitmap, "page.webp")
```

This feature is still experimental and currently does not work on Windows.

## Limitations and related packages

### Tables

Data scientists are often interested in data from tables. Unfortunately the pdf format is pretty dumb and does not have notion of a table (unlike for example HTML). Tabular data in a pdf file is nothing more than strategically positioned lines and text, which makes it difficult to extract the raw data with `pdftools`.

```r
txt <- pdf_text("http://arxiv.org/pdf/1406.4806.pdf")

# some tables
cat(txt[18])
cat(txt[19])
```

The [`tabulizer`](https://github.com/ropensci/tabulizer) package is dedicated to extracting tables from PDF, and includes interactive tools for selecting tables. However, `tabulizer` depends on `rJava` and therefore requires additional setup steps or may be impossible to use on systems where Java cannot be installed.

It is possible to use `pdftools` with some creativity to parse tables from PDF documents, which does not require Java to be installed.

### Scanned text

If you want to extract text from scanned text present in a pdf, you'll need to use OCR (optical character recognition). Please refer to the [rOpenSci `tesseract` package](https://github.com/ropensci/tesseract) that provides bindings to the Tesseract OCR engine. In particular read [the section of its vignette about reading from PDF files using `pdftools` and `tesseract`](https://cran.r-project.org/web/packages/tesseract/vignettes/intro.html#read_from_pdf_files).


[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qpdf.R
\docType{import}
\name{qpdf}
\alias{qpdf}
\alias{pdf_combine}
\alias{reexports}
\alias{pdf_compress}
\alias{pdf_length}
\alias{pdf_split}
\alias{pdf_subset}
\title{qpdf utilities}
\seealso{
Other pdftools: 
\code{\link{pdf_ocr_text}()},
\code{\link{pdftools}},
\code{\link{rendering}}
}
\concept{pdftools}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{qpdf}{\code{\link[qpdf:qpdf]{pdf_combine}}, \code{\link[qpdf:qpdf]{pdf_compress}}, \code{\link[qpdf:qpdf]{pdf_length}}, \code{\link[qpdf:qpdf]{pdf_split}}, \code{\link[qpdf:qpdf]{pdf_subset}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/render.R
\name{rendering}
\alias{rendering}
\alias{pdf_render_page}
\alias{render}
\alias{pdf_convert}
\alias{poppler_config}
\title{Render / Convert PDF}
\usage{
pdf_render_page(
  pdf,
  page = 1,
  dpi = 72,
  numeric = FALSE,
  antialias = TRUE,
  opw = "",
  upw = ""
)

pdf_convert(
  pdf,
  format = "png",
  pages = NULL,
  filenames = NULL,
  dpi = 72,
  antialias = TRUE,
  opw = "",
  upw = "",
  verbose = TRUE
)

poppler_config()
}
\arguments{
\item{pdf}{file path or raw vector with pdf data}

\item{page}{which page to render}

\item{dpi}{resolution (dots per inch) to render}

\item{numeric}{convert raw output to (0-1) real values}

\item{antialias}{enable antialiasing. Must be \code{"text"} or \code{"draw"} or \code{TRUE} (both)
or \code{FALSE} (neither).}

\item{opw}{owner password}

\item{upw}{user password}

\item{format}{string with output format such as \code{"png"} or \code{"jpeg"}. Must be equal
to one of \code{poppler_config()$supported_image_formats}.}

\item{pages}{vector with one-based page numbers to render. \code{NULL} means all pages.}

\item{filenames}{vector of equal length to \code{pages} with output filenames. May also be
a format string which is expanded using \code{pages} and \code{format} respectively.}

\item{verbose}{print some progress info to stdout}
}
\description{
High quality conversion of pdf page(s) to png, jpeg or tiff format, or render into a
raw bitmap array for further processing in R.
}
\examples{
# Rendering should be supported on all platforms now
# convert few pages to png
file.copy(file.path(Sys.getenv("R_DOC_DIR"), "NEWS.pdf"), "news.pdf")
pdf_convert("news.pdf", pages = 1:3)

# render into raw bitmap
bitmap <- pdf_render_page("news.pdf")

# save to bitmap formats
png::writePNG(bitmap, "page.png")
jpeg::writeJPEG(bitmap, "page.jpeg")
webp::write_webp(bitmap, "page.webp")

# Higher quality
bitmap <- pdf_render_page("news.pdf", page = 1, dpi = 300)
png::writePNG(bitmap, "page.png")

# slightly more efficient
bitmap_raw <- pdf_render_page("news.pdf", numeric = FALSE)
webp::write_webp(bitmap_raw, "page.webp")

# Cleanup
unlink(c('news.pdf', 'news_1.png', 'news_2.png', 'news_3.png',
 'page.jpeg', 'page.png', 'page.webp'))
}
\seealso{
Other pdftools: 
\code{\link{pdf_ocr_text}()},
\code{\link{pdftools}},
\code{\link{qpdf}}
}
\concept{pdftools}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{pdftools}
\alias{pdftools}
\alias{pdf_info}
\alias{pdf_text}
\alias{pdf_data}
\alias{pdf_fonts}
\alias{pdf_attachments}
\alias{pdf_toc}
\alias{pdf_pagesize}
\title{PDF utilities}
\usage{
pdf_info(pdf, opw = "", upw = "")

pdf_text(pdf, opw = "", upw = "")

pdf_data(pdf, font_info = FALSE, opw = "", upw = "")

pdf_fonts(pdf, opw = "", upw = "")

pdf_attachments(pdf, opw = "", upw = "")

pdf_toc(pdf, opw = "", upw = "")

pdf_pagesize(pdf, opw = "", upw = "")
}
\arguments{
\item{pdf}{file path or raw vector with pdf data}

\item{opw}{string with owner password to open pdf}

\item{upw}{string with user password to open pdf}

\item{font_info}{if TRUE, extract font-data for each box. Be careful, this
requires a very recent version of poppler and will error otherwise.}
}
\description{
Utilities based on libpoppler for extracting text, fonts, attachments
and metadata from a pdf file.
}
\details{
The \code{\link{pdf_text}} function renders all textboxes on a text canvas
and returns a character vector of equal length to the number of pages in the
PDF file. On the other hand, \code{\link{pdf_data}} is more low level and
returns one data frame per page, containing one row for each textbox in the PDF.

Note that \code{\link{pdf_data}} requires a recent version of libpoppler
which might not be available on all Linux systems.
When using \code{\link{pdf_data}} in R packages, condition use on
\code{poppler_config()$has_pdf_data} which shows if this function can be
used on the current system. For Ubuntu 16.04 (Xenial) and 18.04 (Bionic)
you can use \href{https://github.com/ropensci/pdftools#installation}{the PPA}
with backports of Poppler 0.74.0.

Poppler is pretty verbose when encountering minor errors in PDF files,
in especially \code{\link{pdf_text}}. These messages are usually safe
to ignore, use \code{\link{suppressMessages}} to hide them altogether.
}
\examples{
# Just a random pdf file
pdf_file <- file.path(R.home("doc"), "NEWS.pdf")
info <- pdf_info(pdf_file)
text <- pdf_text(pdf_file)
fonts <- pdf_fonts(pdf_file)
files <- pdf_attachments(pdf_file)
}
\seealso{
Other pdftools: 
\code{\link{pdf_ocr_text}()},
\code{\link{qpdf}},
\code{\link{rendering}}
}
\concept{pdftools}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ocr.R
\name{pdf_ocr_text}
\alias{pdf_ocr_text}
\alias{pdf_ocr_data}
\title{OCR text extraction}
\usage{
pdf_ocr_text(
  pdf,
  pages = NULL,
  opw = "",
  upw = "",
  language = "eng",
  dpi = 600
)

pdf_ocr_data(
  pdf,
  pages = NULL,
  opw = "",
  upw = "",
  language = "eng",
  dpi = 600
)
}
\arguments{
\item{pdf}{file path or raw vector with pdf data}

\item{pages}{which pages of the pdf file to extract}

\item{opw}{string with owner password to open pdf}

\item{upw}{string with user password to open pdf}

\item{language}{passed to \link[tesseract:tesseract]{tesseract} to specify the
languge of the engine.}

\item{dpi}{resolution to render image that is passed to \link[tesseract:ocr]{tesseract::ocr}.}
}
\description{
Perform OCR text extraction. This requires you have the \code{tesseract} package.
}
\seealso{
Other pdftools: 
\code{\link{pdftools}},
\code{\link{qpdf}},
\code{\link{rendering}}
}
\concept{pdftools}
