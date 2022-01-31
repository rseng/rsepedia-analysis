# tesseract

> Bindings to [Tesseract-OCR](https://opensource.google/projects/tesseract): 
  a powerful optical character recognition (OCR) engine that supports over 100 languages.
  The engine is highly configurable in order to tune the detection algorithms and
  obtain the best possible results.

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tesseract)](https://cran.r-project.org/package=tesseract)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/tesseract)](https://cran.r-project.org/package=tesseract)

 - Upstream Tesseract-OCR documentation: https://tesseract-ocr.github.io/tessdoc/
 - Introduction: https://docs.ropensci.org/tesseract/articles/intro.html
 - Reference: https://docs.ropensci.org/tesseract/reference/ocr.html

## Hello World

Simple example

```r
# Simple example
text <- ocr("https://jeroen.github.io/images/testocr.png")
cat(text)

# Get XML HOCR output
xml <- ocr("https://jeroen.github.io/images/testocr.png", HOCR = TRUE)
cat(xml)
```

Roundtrip test: render PDF to image and OCR it back to text

```r
# Full roundtrip test: render PDF to image and OCR it back to text
curl::curl_download("https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf", "R-intro.pdf")
orig <- pdftools::pdf_text("R-intro.pdf")[1]

# Render pdf to png image
img_file <- pdftools::pdf_convert("R-intro.pdf", format = 'tiff', pages = 1, dpi = 400)

# Extract text from png image
text <- ocr(img_file)
unlink(img_file)
cat(text)
```

## Installation

On Windows and MacOS the package binary package can be installed from CRAN:

```r
install.packages("tesseract")
```

Installation from source on Linux or OSX requires the `Tesseract` library (see below).

### Install from source

 On __Debian__ or __Ubuntu__ install [libtesseract-dev](https://packages.debian.org/testing/libtesseract-dev) and
[libleptonica-dev](https://packages.debian.org/testing/libleptonica-dev). Also install [tesseract-ocr-eng](https://packages.debian.org/testing/tesseract-ocr-eng) to run examples.

```
sudo apt-get install -y libtesseract-dev libleptonica-dev tesseract-ocr-eng
```

On __Ubuntu__ you can optionally use [this PPA](https://launchpad.net/~alex-p/+archive/ubuntu/tesseract-ocr-devel) to get the latest version of Tesseract:

```
sudo add-apt-repository ppa:alex-p/tesseract-ocr-devel
sudo apt-get install -y libtesseract-dev tesseract-ocr-eng
```

On __Fedora__ we need [tesseract-devel](https://src.fedoraproject.org/rpms/tesseract) and
[leptonica-devel](https://src.fedoraproject.org/rpms/leptonica)

```
sudo yum install tesseract-devel leptonica-devel
````

On __RHEL__ and __CentOS__ we need [tesseract-devel](https://src.fedoraproject.org/rpms/tesseract) and
[leptonica-devel](https://src.fedoraproject.org/rpms/leptonica) from EPEL

```
sudo yum install epel-release
sudo yum install tesseract-devel leptonica-devel
````


On __OS-X__ use [tesseract](https://github.com/Homebrew/homebrew-core/blob/master/Formula/tesseract.rb) from Homebrew:

```
brew install tesseract
```

Tesseract uses training data to perform OCR. Most systems default to English
training data. To improve OCR results for other languages you can to install the
appropriate training data. On Windows and OSX you can do this in R using 
`tesseract_download()`:


```r
tesseract_download('fra')
```

On Linux you need to install the appropriate training data from your distribution. 
For example to install the spanish training data:

  - [tesseract-ocr-spa](https://packages.debian.org/testing/tesseract-ocr-spa) (Debian, Ubuntu)
  - [tesseract-langpack-spa](https://src.fedoraproject.org/rpms/tesseract-langpack) (Fedora, EPEL)

Alternatively you can manually download training data from [github](https://github.com/tesseract-ocr/tessdata)
and store it in a path on disk that you pass in the `datapath` parameter or set a default path via the
`TESSDATA_PREFIX` environment variable. Note that the Tesseract 4 and Tesseract 3 use different 
training data format. Make sure to download training data from the branch that matches your libtesseract version.

---
title: "Using the Tesseract OCR engine in R"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 2
    toc_float: true
    fig_caption: false
vignette: >
  %\VignetteIndexEntry{Using the Tesseract OCR engine in R}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r, echo = FALSE, message = FALSE}
library(tibble)
#knitr::opts_chunk$set(comment = "")
has_nld <- "nld" %in% tesseract::tesseract_info()$available
if(identical(Sys.info()[['user']], 'jeroen')) stopifnot(has_nld)
```

The tesseract package provides R bindings [Tesseract](https://opensource.google/projects/tesseract): a powerful optical character recognition (OCR) engine that supports over 100 languages. The engine is highly configurable in order to tune the detection algorithms and obtain the best possible results.

Keep in mind that OCR (pattern recognition in general) is a very difficult problem for computers. Results will rarely be perfect and the accuracy rapidly decreases with the quality of the input image. But if you can get your input images to reasonable quality, Tesseract can often help to extract most of the text from the image.

## Extract Text from Images

OCR is the process of finding and recognizing text inside images, for example from a screenshot, scanned paper. The image below has some example text:

![test](https://jeroen.github.io/images/testocr.png){data-external=1}

```{r}
library(tesseract)
eng <- tesseract("eng")
text <- tesseract::ocr("http://jeroen.github.io/images/testocr.png", engine = eng)
cat(text)
```

Not bad! The `ocr_data()` function returns all words in the image along with a bounding box and confidence rate.

```{r}
results <- tesseract::ocr_data("http://jeroen.github.io/images/testocr.png", engine = eng)
results
```

## Language Data

The tesseract OCR engine uses language-specific training data in the recognize words. The OCR algorithms bias towards words and sentences that frequently appear together in a given language, just like the human brain does. Therefore the most accurate results will be obtained when using training data in the correct language. 

Use `tesseract_info()` to list the languages that you currently have installed.

```{r}
tesseract_info()
```

By default the R package only includes English training data. Windows and Mac users can install additional training data using `tesseract_download()`. Let's OCR a screenshot from Wikipedia in Dutch (Nederlands) 

[![utrecht](https://jeroen.github.io/images/utrecht2.png)](https://nl.wikipedia.org/wiki/Geschiedenis_van_de_stad_Utrecht)

```{r, eval=FALSE}
# Only need to do download once:
tesseract_download("nld")
```

```{r eval = has_nld}
# Now load the dictionary
(dutch <- tesseract("nld"))
text <- ocr("https://jeroen.github.io/images/utrecht2.png", engine = dutch)
cat(text)
```

As you can see immediately: almost perfect! (OK just take my word). 


## Preprocessing with Magick

The accuracy of the OCR process depends on the quality of the input image. You can often improve results by properly scaling the image, removing noise and artifacts or cropping the area where the text exists. See [tesseract wiki: improve quality](https://github.com/tesseract-ocr/tesseract/wiki/ImproveQuality) for important tips to improve the quality of your input image.

The awesome [magick](https://cran.r-project.org/package=magick/vignettes/intro.html) R package has many useful functions that can be use for enhancing the quality of the image. Some things to try:

 - If your image is skewed, use `image_deskew()` and `image_rotate()` make the text horizontal.
 - `image_trim()` crops out whitespace in the margins. Increase the `fuzz` parameter to make it work for noisy whitespace.
 - Use `image_convert()` to turn the image into greyscale, which can reduce artifacts and enhance actual text.
 - If your image is very large or small resizing with `image_resize()` can help tesseract determine text size.
 - Use `image_modulate()` or `image_contrast()` or `image_contrast()` to tweak brightness / contrast if this is an issue.
 - Try `image_reducenoise()` for automated noise removal. Your mileage may vary.
 - With `image_quantize()` you can reduce the number of colors in the image. This can sometimes help with increasing contrast and reducing artifacts.
 - True imaging ninjas can use `image_convolve()` to use custom [convolution methods](https://ropensci.org/technotes/2017/11/02/image-convolve/). 

Below is an example OCR scan from an online [AI course](https://courses.cs.vt.edu/csonline/AI/Lessons/VisualProcessing/OCRscans.html). The code converts it to black-and-white and resizes + crops the image before feeding it to tesseract to get more accurate OCR results.

![bowers](https://jeroen.github.io/images/bowers.jpg){data-external=1}


```{r}
library(magick)
input <- image_read("https://jeroen.github.io/images/bowers.jpg")

text <- input %>%
  image_resize("2000x") %>%
  image_convert(type = 'Grayscale') %>%
  image_trim(fuzz = 40) %>%
  image_write(format = 'png', density = '300x300') %>%
  tesseract::ocr() 

cat(text)
```


## Read from PDF files

If your images are stored in PDF files they first need to be converted to a proper image format. We can do this in R using the `pdf_convert` function from the pdftools package. Use a high DPI to keep quality of the image.

```{r, eval=require(pdftools)}
pngfile <- pdftools::pdf_convert('https://jeroen.github.io/images/ocrscan.pdf', dpi = 600)
text <- tesseract::ocr(pngfile)
cat(text)
```


## Tesseract Control Parameters

Tesseract supports hundreds of [control parameters](https://tesseract-ocr.github.io/tessdoc/ControlParams) which alter the OCR engine. Use `tesseract_params()` to list all parameters with their default value and a brief description. It also has a handy `filter` argument to quickly find parameters that match a particular string.

```{r}
# List all parameters with *colour* in name or description
tesseract_params('colour')
```

Do note that some of the control parameters have changed between Tesseract engine 3 and 4.

```{r}
tesseract::tesseract_info()['version']
```

### Whitelist / Blacklist characters

One powerful parameter is `tessedit_char_whitelist` which restricts the output to a limited set of characters. This may be useful for reading for example numbers such as a bank account, zip code, or gas meter.

The whitelist parameter works for all versions of Tesseract engine 3 and also engine versions 4.1 and higher, but unfortunately it did not work in Tesseract 4.0.


![receipt](https://jeroen.github.io/images/receipt.png){data-external=1}

```{r}
numbers <- tesseract(options = list(tessedit_char_whitelist = "$.0123456789"))
cat(ocr("https://jeroen.github.io/images/receipt.png", engine = numbers))
```

To test if this actually works, look what happens if we remove the `$` from `tessedit_char_whitelist`:

```{r}
# Do not allow any dollar sign 
numbers2 <- tesseract(options = list(tessedit_char_whitelist = ".0123456789"))
cat(ocr("https://jeroen.github.io/images/receipt.png", engine = numbers2))
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tesseract.R
\name{tesseract}
\alias{tesseract}
\alias{tesseract_params}
\alias{tesseract_info}
\title{Tesseract Engine}
\usage{
tesseract(
  language = "eng",
  datapath = NULL,
  configs = NULL,
  options = NULL,
  cache = TRUE
)

tesseract_params(filter = "")

tesseract_info()
}
\arguments{
\item{language}{string with language for training data. Usually defaults to \code{eng}}

\item{datapath}{path with the training data for this language. Default uses
the system library.}

\item{configs}{character vector with files, each containing one or more parameter
values. These config files can exist in the current directory or one of the standard
tesseract config files that live in the tessdata directory. See details.}

\item{options}{a named list with tesseract parameters. See details.}

\item{cache}{speed things up by caching engines}

\item{filter}{only list parameters containing a particular string}
}
\description{
Create an OCR engine for a given language and control parameters. This can be used by
the \link{ocr} and \link{ocr_data} functions to recognize text.
}
\details{
Tesseract \href{https://tesseract-ocr.github.io/tessdoc/ControlParams}{control parameters}
can be set either via a named list in the
\code{options} parameter, or in a \code{config} file text file which contains the parameter name
followed by a space and then the value, one per line. Use \code{\link[=tesseract_params]{tesseract_params()}} to list
or find parameters. Note that that some parameters are only supported in certain versions
of libtesseract, and that invalid parameters can sometimes cause libtesseract to crash.
}
\examples{
tesseract_params('debug')
}
\references{
\href{https://tesseract-ocr.github.io/tessdoc/ControlParams}{tesseract wiki: control parameters}
}
\seealso{
Other tesseract: 
\code{\link{ocr}()},
\code{\link{tesseract_download}()}
}
\concept{tesseract}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ocr.R
\name{ocr}
\alias{ocr}
\alias{ocr_data}
\title{Tesseract OCR}
\usage{
ocr(image, engine = tesseract("eng"), HOCR = FALSE)

ocr_data(image, engine = tesseract("eng"))
}
\arguments{
\item{image}{file path, url, or raw vector to image (png, tiff, jpeg, etc)}

\item{engine}{a tesseract engine created with \code{\link[=tesseract]{tesseract()}}. Alternatively a
language string which will be passed to \code{\link[=tesseract]{tesseract()}}.}

\item{HOCR}{if \code{TRUE} return results as HOCR xml instead of plain text}
}
\description{
Extract text from an image. Requires that you have training data for the language you
are reading. Works best for images with high contrast, little noise and horizontal text.
See \href{https://github.com/tesseract-ocr/tesseract/wiki/ImproveQuality}{tesseract wiki} and
our package vignette for image preprocessing tips.
}
\details{
The \code{ocr()} function returns plain text by default, or hOCR text if hOCR is set to \code{TRUE}.
The \code{ocr_data()} function returns a data frame with a confidence rate and bounding box for
each word in the text.
}
\examples{
# Simple example
text <- ocr("https://jeroen.github.io/images/testocr.png")
cat(text)

xml <- ocr("https://jeroen.github.io/images/testocr.png", HOCR = TRUE)
cat(xml)

df <- ocr_data("https://jeroen.github.io/images/testocr.png")
print(df)

\donttest{
# Full roundtrip test: render PDF to image and OCR it back to text
curl::curl_download("https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf", "R-intro.pdf")
orig <- pdftools::pdf_text("R-intro.pdf")[1]

# Render pdf to png image
img_file <- pdftools::pdf_convert("R-intro.pdf", format = 'tiff', pages = 1, dpi = 400)
unlink("R-intro.pdf")

# Extract text from png image
text <- ocr(img_file)
unlink(img_file)
cat(text)
}

engine <- tesseract(options = list(tessedit_char_whitelist = "0123456789"))
}
\references{
\href{https://github.com/tesseract-ocr/tesseract/wiki/ImproveQuality}{Tesseract: Improving Quality}
}
\seealso{
Other tesseract: 
\code{\link{tesseract_download}()},
\code{\link{tesseract}()}
}
\concept{tesseract}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tessdata.R
\name{tesseract_download}
\alias{tesseract_download}
\alias{tessdata}
\title{Tesseract Training Data}
\usage{
tesseract_download(lang, datapath = NULL, progress = interactive())
}
\arguments{
\item{lang}{three letter code for language, see \href{https://github.com/tesseract-ocr/tessdata}{tessdata} repository.}

\item{datapath}{destination directory where to download store the file}

\item{progress}{print progress while downloading}
}
\description{
Helper function to download training data from the official
\href{https://tesseract-ocr.github.io/tessdoc/Data-Files}{tessdata} repository. Only use this function on
Windows and OS-X. On Linux, training data can be installed directly with
\href{https://src.fedoraproject.org/rpms/tesseract}{yum} or
\href{https://packages.debian.org/search?suite=stable&section=all&arch=any&searchon=names&keywords=tesseract-ocr-}{apt-get}.
}
\details{
Tesseract uses training data to perform OCR. Most systems default to English
training data. To improve OCR performance for other languages you can to install the
training data from your distribution. For example to install the spanish training data:
\itemize{
\item \href{https://packages.debian.org/testing/tesseract-ocr-spa}{tesseract-ocr-spa} (Debian, Ubuntu)
\item \code{tesseract-langpack-spa} (Fedora, EPEL)
}

On Windows and MacOS you can install languages using the \link{tesseract_download} function
which downloads training data directly from \href{https://github.com/tesseract-ocr/tessdata}{github}
and stores it in a the path on disk given by the \code{TESSDATA_PREFIX} variable.
}
\examples{
\dontrun{
if(is.na(match("fra", tesseract_info()$available)))
  tesseract_download("fra")
french <- tesseract("fra")
text <- ocr("https://jeroen.github.io/images/french_text.png", engine = french)
cat(text)
}
}
\references{
\href{https://tesseract-ocr.github.io/tessdoc/Data-Files}{tesseract wiki: training data}
}
\seealso{
Other tesseract: 
\code{\link{ocr}()},
\code{\link{tesseract}()}
}
\concept{tesseract}
