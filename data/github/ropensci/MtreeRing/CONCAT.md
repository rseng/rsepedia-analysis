
MtreeRing
=======

<!--require(knitr);require(markdown);knit("README.Rmd")-->


**Authors:** [Jingning Shi](https://www.researchgate.net/profile/Jingning-Shi), [Wei Xiang](https://www.researchgate.net/profile/Wei-Xiang-11)<br/>
**License:** [GPL3](https://cran.r-project.org/web/licenses/GPL-3)

<!--pkg badges-->
[![TravisCI Build Status](https://api.travis-ci.org/ropensci/MtreeRing.svg?branch=master)](https://api.travis-ci.org/ropensci/MtreeRing.svg?branch=master)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/cti7i110hecl8kpf?svg=true)](https://ci.appveyor.com/project/JingningShi/MtreeRing)
[![codecov](https://codecov.io/github/ropensci/MtreeRing/coverage.svg?branch=master)](https://codecov.io/github/ropensci/MtreeRing?branch=master)
[![ropensci](https://badges.ropensci.org/287_status.svg)](https://github.com/ropensci/software-review/issues/287)
[![cran checks](https://cranchecks.info/badges/worst/MtreeRing)](https://cranchecks.info/pkgs/MtreeRing)
[![Downloads](https://cranlogs.r-pkg.org/badges/MtreeRing)](https://CRAN.R-project.org/package=MtreeRing)

`MtreeRing` is a tool for automatically measuring tree-ring width using image processing techniques.

## Installation

Install the stable version from CRAN


```r
install.packages("MtreeRing")
```

or the development version from GitHub


```r
# install.packages("devtools")
devtools::install_github("ropensci/MtreeRing")
```

## Ring-width measurement

### 1. Read an image


```r
library(MtreeRing)
## Read and plot a tree ring image
img.name <- system.file("001.png", package = "MtreeRing")
t1 <- ring_read(img = img.name, dpi = 1200, plot = TRUE)
```

`ring_read` supports commonly used image formats, including png, tiff, jpg and bmp.

### 2. Detect ring borders 

After plotting the image, the automatic detection of ring borders can be performed using three alternative methods: (1) watershed algorithm; (2) Canny edge detector; (3) a linear detection algorithm from R package [measuRing](https://CRAN.R-project.org/package=measuRing).


```r
## Split a long core sample into 2 pieces to
## get better display performance and use the
## watershed algorithm to detect ring borders:
t2 <- ring_detect(ring.data = t1, seg = 2, method = 'watershed')
```

<center><img src="inst/README-img001.png" width = "65%" height = "65%" /></center>
<center>Figure 1. The automatic detection of ring borders</center>

### 3. Calculate ring-width series 

If all ring borders are correctly identified, you can generate a ring-width series in data frame format. Use `write.rwl` to export the ring-width series to an rwl file.


```r
rw.df <- ring_calculate(ring.data = t2, seriesID = "940220")
library(dplR) # A dendrochronological analysis package
fn <- tempfile(fileext=".rwl")
write.rwl(rwl.df = rw.df, fname = fn, format = "tucson")
```


## Shiny application

If you are not familiar with R and its command line interface, the shiny-based app is a good alternative.


```r
MtreeRing::ring_app_launch()
```

This command allows to run a Shiny-based application within the system's default web browser. The app provides a beginner-friendly graphical interface and supports more flexible mouse-based interactions, allowing image files to be uploaded up to 150 MB in size.

The dashboard has three components: a header, sidebar and body, like this

<img src="inst/README-img002.png" width = "85%" height = "85%" />

A workflow for the Shiny app can be found at https://ropensci.github.io/MtreeRing/articles/app-MtreeRing.html. Most steps are demonstrated with a gif to make the workflow more understandable.

## Ring width correction

If an increment borer is used to extract samples, it is well known that the auger sometimes fails to traverse the pith of the sampled tree but passes through one side of the pith at a certain distance. Tangent lines of rings close to the pith are therefore not perpendicular to the horizontal path, which may lead to considerable errors in ring widths.

Under such conditions, you can create two paths by setting the argument `incline = TRUE`, or by ticking the checkbox **Inclined tree rings**. See this example.

<img src="inst/RingCorrection.png" width = "80%" height = "80%" /> 

The line segment connecting two dots on the same ring should match the tangent of a tree ring border. The corrected ring width is estimated from the distance between adjacent rings and orientation of ring borders.

## Support

Any feedback, bug reports or suggestions are welcomed. If you have a comment on `MtreeRing`, or you find a bug in the released or beta versions, please submit bugs and/or feature requests at https://github.com/ropensci/MtreeRing/issues. Include the package version, OS, and any command-line required to reproduce the problem.

## Code of conduct

I will try to add new features based on user feedback. It is hoped that others will contribute additional useful features. Please note that the 'MtreeRing' project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
MtreeRing News
=======
### MtreeRing 1.4.5 (2021-04-20)
- Update CSS styles of the graphical user interface

### MtreeRing 1.4.4 (2021-03-15)
- Update links in the README file

### MtreeRing 1.4.3 (2021-03-14)
- Replace ‘spatstat’ with ‘spatstat.geom’ based on feedback from the Spatstat Team

### MtreeRing 1.4.2 (2019-10-03)
- Update GitHub links in `ring_app_launch` function

### MtreeRing 1.4.1 (2019-09-23)
- Update README file

### MtreeRing 1.4.0 (2019-09-01)
- Redesign the interface and measurement process for the Shiny app
- Update documentation for the redesigned Shiny app
- Add a Support section in the README file
- Replace all vignettes with GitHub links

### MtreeRing 1.3.1 (2019-05-29)
- Update README file

### MtreeRing 1.3 (2019-05-20)
- Add a NEWS.md file to track changes
- Add CODE_OF_CONDUCT.md
- Add `plot` parameter to `ring_read`
- Rename all exported functions using snake_case style
- Update documentation for the Shiny app and mouse actions
- Update README file
- Fix a bug where `ring_detect` fails to detect ring borders
- Remove intro vignette and add three new vignettes to demonstrate major functionality

### MtreeRing 1.2 (2019-03-02)
- Add appveyor-CI and travis-CI
- Add README file
- Add tests for all functions and Shiny app
- Add intro vignette
- Update documentation for `autoDetect`
- Update DESCRIPTION URL and bug report

### MtreeRing 1.1 (2018-10-16)
- Initial CRAN release
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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pith_measure.R
\name{pith_measure}
\alias{pith_measure}
\title{Calibrate ring-width series}
\usage{
pith_measure(ring.data, inner.arc = TRUE, last.yr = NULL,
  color = "black", border.type = 16, label.cex = 1.5)
}
\arguments{
\item{ring.data}{A magick image object produced by \code{\link{ring_read}}.}

\item{inner.arc}{A logical value indicating whether to calibrate the 
ring-width series using the arcs of inner rings. See details below.}

\item{last.yr}{\code{NULL} or an integer giving the year of formation 
of the left-most ring. If \code{NULL}, border numbers (starting from 1) 
are used instead of years.}

\item{color}{Color for labels.}

\item{border.type}{Symbol for ring borders. See \code{pch} in 
\code{\link{points}} for possible values and shapes.}

\item{label.cex}{The magnification to be used for years or border numbers.}
}
\value{
A data frame of the calibrated ring-width series. The measurements 
units are millimeters (mm)
}
\description{
This function can calibrate the ring-width series 
using arcs of inner rings.
}
\details{
This function allows the user to create a path, and manually mark 
ring borders by clicking on the graphical window. 

An example demonstrated with pictures can be found in the package vignette. 
Type \code{vignette('pith-MtreeRing')} to see this example.

\itemize{
\item
If \code{inner.arc = TRUE}, the ring-width series is calibrated using arcs 
of inner rings (Duncan, 1989).
 
\bold{Step1}. You can click the left mouse button to add a horizontal path.
The path should traverse an appropriate arc (read the reference below  
for more details).

\bold{Step2}. You can add three points to the selected arc by
left-clicking. The first point should be placed on the left endpoint of 
the arc, and the second point is placed on the right endpoint. 

After adding these two points, a vertical dashed line will be plotted 
automatically according to the (x,y) positions of endpoints you just added. 
The third points should be placed on the intersection of the vertical 
dashed line and the selected arc. 

\bold{Step3}. you are prompted to mark tree rings along the path by 
left-clicking on the image. Every click draws a point.
Note that the left endpoint of the arc will be considered as the last 
ring border without the need to mark it. 

After marking tree rings, the identification process does not automatically 
stop by itself. On the Windows platform, the identification process 
can be terminated by clicking the second button and selecting \bold{Stop} 
from the menu. On the MacOS system, you can press the \bold{Escape} key to 
terminate this process.

The ring-width series are corrected using formulas proposed by Duncan (1989).

\item
If \code{inner.arc = FALSE}, the user can create a path which matches 
the direction of wood growth. 

\bold{Step1}. You can add two points by left-clicking on the image. 
Every click draws a point.
A path passing through these two points will be plotted. The path should 
follow the rays from bark to pith.

\bold{Step2}. You can mark tree rings along the path by left-clicking
on the image. The termination of identification process is similar.
}
}
\examples{
img.path <- system.file("missing_pith.png", package = "MtreeRing")

## Read the image:
t1 <- ring_read(img = img.path, dpi = 1200, plot = FALSE)

## Use the arcs of inner rings to calibrate ring-width series:
\donttest{t2 <- pith_measure(t1, inner.arc = TRUE, last.yr = 2016)}

## Try another method to measure ring widths:
\donttest{t3 <- pith_measure(t1, inner.arc = FALSE, last.yr = 2016)}
}
\references{
Duncan R. (1989) 
An evaluation of errors in tree age estimates based on increment cores 
in Kahikatea (Dacrycarpus dacrydiodes).
\emph{New Zealand Natural Sciences}
\bold{16(4)}, 1-37.
}
\author{
Jingning Shi
}
\name{MtreeRing-package}
\alias{MtreeRing-package}
\alias{MtreeRing}
\docType{package}
\title{A Shiny Application for Automatic Measurements of Tree-Ring Widths 
on Digital Images}
\description{
Use morphological image processing and edge detection algorithms to 
automatically measure tree ring widths on digital images. Users can also 
manually mark tree rings on species with complex anatomical structures. 
The arcs of inner-rings and angles of successive inclined ring boundaries 
are used to correct ring-width series. The package provides a Shiny-based 
application, allowing R beginners to easily analyze tree ring images and 
export ring-width series in standard file formats.
}
\details{

\tabular{ll}{
Package: \tab MtreeRing\cr
Type: \tab Package\cr
License: \tab GPL-3\cr
Maintainer: \tab Jingning Shi <snow940220@bjfu.edu.cn>\cr
LazyData: \tab TRUE\cr
}

}
\author{Jingning Shi, Wei Xiang}
\keyword{ package }% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ring_detect.R
\name{ring_detect}
\alias{ring_detect}
\title{Automatic detection of tree-ring borders}
\usage{
ring_detect(ring.data, seg = 1, auto.path = TRUE, manual = FALSE,
  method = "canny", incline = FALSE, sample.yr = NULL,
  watershed.threshold = "auto", watershed.adjust = 0.8,
  struc.ele1 = NULL, struc.ele2 = NULL, marker.correction = FALSE,
  default.canny = TRUE, canny.t1, canny.t2, canny.smoothing = 2,
  canny.adjust = 1.4, path.dis = 1, origin = 0,
  border.color = "black", border.type = 16, label.color = "black",
  label.cex = 1.2)
}
\arguments{
\item{ring.data}{A magick image object produced by \code{\link{ring_read}}.}

\item{seg}{An integer specifying the number of image segments.}

\item{auto.path}{A logical value. If \code{TRUE}, a path is automatically 
created at the center of the image. If \code{FALSE}, the function allows 
the user to create a sub-image and a path by interactive clickings. 
See details below.}

\item{manual}{A logical value indicating whether to skip the automatic 
detection. If \code{TRUE}, ring borders are visually identified after 
creating the path. See \code{\link{ring_modify}} to learn how to mark
tree rings by clicking on the image.}

\item{method}{A character string specifying how ring borders are detected. 
It requires one of the following characters: \code{"watershed"}, 
\code{"canny"}, or \code{"lineardetect"}. See details below.}

\item{incline}{A logical value indicating whether to correct ring widths. 
If \code{TRUE}, two horizontal paths are added to the image.}

\item{sample.yr}{\code{NULL} or an integer giving the year of formation 
of the left-most ring. If \code{NULL}, use the current year.}

\item{watershed.threshold}{The threshold used for producing the marker 
image, either a numeric from 0 to 1, or the character "auto" (using the 
Otsu algorithm), or a character of the form "XX\%" (e.g., "58\%").}

\item{watershed.adjust}{A numeric used to adjust the Otsu threshold. 
The default is 1 which means that the threshold will not be adjusted. 
The sizes of early-wood regions in the marker image will reduce along 
with the decrease of \code{watershed.adjust}.}

\item{struc.ele1}{\code{NULL} or a vector of length two specifying the 
width and height of the first structuring element. If \code{NULL}, the 
size of the structuring element is determined by the argument \code{dpi}.}

\item{struc.ele2}{\code{NULL} or a vector of length two specifying the 
width and height of the second structuring element. If \code{NULL}, the 
size of the structuring element is determined by the argument \code{dpi}.}

\item{marker.correction}{A logical value indicating whether to relabel 
early-wood regions by comparing the values of their left-side neighbours.}

\item{default.canny}{A logical value. If \code{TRUE}, upper and lower 
Canny thresholds are determined automatically.}

\item{canny.t1}{A numeric giving the threshold for weak edges.}

\item{canny.t2}{A numeric giving the threshold for strong edges.}

\item{canny.smoothing}{An integer specifying the degree of smoothing.}

\item{canny.adjust}{A numeric used as a sensitivity control factor for 
the Canny edge detector. The default is 1 which means that the sensitivity 
will not be adjusted. The number of detected borders will reduce along 
with the increase of this value.}

\item{path.dis}{A numeric specifying the perpendicular distance between 
two paths when the argument \code{incline = TRUE}. The unit is in mm.}

\item{origin}{A numeric specifying the origin in smoothed gray to find 
ring borders. See \code{\link{ringBorders}} from the package 
\code{\link{measuRing}}.}

\item{border.color}{Color for ring borders.}

\item{border.type}{Symbol for ring borders. See \code{pch} in 
\code{\link{points}} for possible values and their shapes.}

\item{label.color}{Color for years and border numbers.}

\item{label.cex}{The magnification to be used for years and border numbers.}
}
\value{
A matrix (grayscale image) or array (color image) 
representing the tree-ring image.
}
\description{
This function is used to automatically detect tree ring 
borders along the user-defined path.
}
\details{
If \code{auto.path = FALSE}, the user can create a rectangular sub-image 
and a horizontal path by interactively clicking on the tree ring image. 
The automatic detection will be performed within this rectangular 
sub-image. 
To create a sub-image and a path, follow these steps.
\itemize{
  \item
  Step 1. Select the left and right edges of the rectangle
  
  The user can point the mouse at any 
  desired locations and click the left mouse button to add each edge. 
  
  \item
  Step 2. Select the top and bottom edges of the rectangle
  
  The user can point the mouse at any desired locations and click the 
  left mouse button to add each edge. The width of the rectangle is 
  defined as the distance between the top and bottom edges, and should 
  not be unnecessarily large to reduce time consumption and memory usage. 
  Creating a long and narrow rectangle if possible.
  
  \item
  Step 3. Create a path
  
  After creating the rectangular sub-image, the user can add a horizontal 
  path by left-clicking on the sub-image (generally at the center of the 
  sub-image, try to choose a clean defect-free area). Ring borders and 
  other markers are plotted along this path. If \code{incline = TRUE}, 
  two paths are added simultaneously.
}

After creating the sub-image and the path, this function will open several 
graphics windows and plot detected ring borders on image segments. The 
number of image segments is controlled by the argument \code{seg}.

Argument \code{method} determines how ring borders are identified. 
\itemize{
  \item
  If \code{method = "watershed"}, this function uses the watershed algorithm 
  to obtain ring borders (Soille and Misson, 2001).
  \item
  If \code{method = "canny"}, this function uses the Canny algorithm 
  to detect borders.
  \item
  If \code{method = "lineardetect"}, a linear detection algorithm from the 
  package \code{\link{measuRing}} is used to identify ring borders (Lara 
  et al., 2015). Note that \code{incline = TRUE} is not supported in this 
  mode, and path will be automatically created at the center of the image. 
}
If the argument \code{method = "watershed"} or \code{"canny"}, the original 
image is processed by morphological openings and closings using rectangular 
structuring elements of increasing size before detecting borders. The first 
small structuring element is used to remove smaller dark spots in early 
wood regions, and the second large structuring element is used to remove 
light strips in late wood regions. More details about morphological 
processing can be found at Soille and Misson (2001).
}
\note{
This function uses \code{\link{locator}} to record mouse 
positions so it only works on "X11", "windows" and "quartz" devices.
}
\examples{
img.path <- system.file("001.png", package = "MtreeRing")

## Read a tree ring image:
t1 <- ring_read(img = img.path, dpi = 1200, plot = FALSE)

## Split a long core sample into 3 pieces to
## get better display performance and use the
## watershed algorithm to detect ring borders:
t2 <- ring_detect(t1, seg = 3, method = 'watershed', border.color = 'green')

}
\references{
Soille, P., Misson, L. (2001)
Tree ring area measurements using morphological image analysis.
\emph{Canadian Journal of Forest Research}
\bold{31}, 1074-1083. {doi: 10.1139/cjfr-31-6-1074}

Lara, W., Bravo, F., Sierra, C.A. (2015)
measuRing: An R package to measure tree-ring widths from scanned images.
\emph{Dendrochronologia}
\bold{34}, 43-50. {doi: 10.1016/j.dendro.2015.04.002}
}
\author{
Jingning Shi
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ring_read.R
\name{ring_read}
\alias{ring_read}
\title{Read and plot a tree-ring image file}
\usage{
ring_read(img, dpi = NULL, RGB = c(0.299, 0.587, 0.114),
  plot = FALSE, rotate = 0, magick = TRUE)
}
\arguments{
\item{img}{A character string indicating the path of the image file. 
Supported formats include png, tiff, jpg and bmp.}

\item{dpi}{An integer specifying the dpi of the image file. A minimum of 
300 dpi is required when running automatic detection.}

\item{RGB}{A numeric vector of length 3 giving the weight of RGB channels.}

\item{plot}{A logical value indicating whether to plot the tree ring image 
when reading it. If \code{FALSE}, the image is not plotted until
function \code{\link{ring_detect}} or \code{\link{pith_measure}} is called.}

\item{rotate}{An integer specifying how many degrees to rotate (clockwise). 
It requires one of the following values:
\code{0}, \code{90}, \code{180} or \code{270}.}

\item{magick}{A logical value. If \code{TRUE}, \code{magick} is used to
read the tree ring image. If \code{FALSE},
packages \code{png}, \code{jpg} and \code{tiff} are used instead.
See details below.}
}
\value{
A magick image object containing the image data.
}
\description{
This function can read an image file from the hard disk and 
plot it in a newly-opened graphics device.
}
\details{
Proper image preparation has a great influence on the measurement of 
ring widths. A tree-ring image should not contain irrelevant or redundant 
features, such as wooden mounts where cores are glued. The larger the file 
size of an image, the slower the image processing operation will be.

\bold{Pith side} of a wood sample should be placed on the \bold{right side} 
of a graphics window. Use \code{rotate} to change its position.

It is highly recommended to use the default value \code{magick = TRUE}, 
because \code{magick} can significantly reduce the memory usage
when reading a large file.
If image data is stored in a non-standard format, image reading may fail.
In that case you can set \code{magick = FALSE} to 
avoid the use of \code{magick}.
}
\examples{
img.path <- system.file("001.png", package = "MtreeRing")

## Read and plot the image:
t1 <- ring_read(img = img.path, dpi = 1200, plot = TRUE)
}
\author{
Jingning Shi
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ring_modify.R
\name{ring_modify}
\alias{ring_modify}
\title{Edit ring borders visually}
\usage{
ring_modify(ring.data, del = NULL, del.u = NULL, del.l = NULL,
  add = FALSE)
}
\arguments{
\item{ring.data}{A matrix or array produced by \code{\link{ring_detect}}.}

\item{del}{A numeric vector giving the border numbers to be removed.}

\item{del.u}{A numeric vector giving the border numbers to be removed 
on the upper path.}

\item{del.l}{A numeric vector giving the border numbers to be removed 
on the lower path.}

\item{add}{A logical value indicating whether to add new ring borders.}
}
\value{
A matrix (grayscale image) or array (color image)
representing the tree ring image.
}
\description{
This function can remove existing ring borders 
or add new borders.
}
\details{
This function is used to remove existing ring borders, or to add new 
borders by interactively clicking on the image segments.

If the user creates one path (\code{incline = FALSE}), the argument 
\code{del} is used to remove ring borders. If the user creates two paths 
(\code{incline = TRUE}), arguments \code{del.u} and \code{del.l} are used 
to remove ring borders.

If \code{add = TRUE}, graphics windows opened by \code{\link{ring_detect}}
will be activated sequentially. When a graphics window is activated, 
the user can add new borders by left-clicking the mouse along the path.
Every click draws a point representing the ring border.
Type \code{vignette('detection-MtreeRing')} to see 
an example of adding ring borders.

The identification process does not automatically stop by itself.

\itemize{
  \item
  On the Windows system, the identification process can be terminated by 
  pressing the right mouse button and selecting \bold{Stop} from the menu.
  \item 
  On the MacOS system, for a X11 device the identification process is 
  terminated by pressing any mouse button other than the first, and for a 
  quartz device this process is terminated by pressing the \bold{ESC} key.
}

Once the user terminates the identification process, the current 
graphics window will be closed automatically, and the graphics window of
the following segment is activated. When all graphics windows are closed,
\code{ring_modify} will re-open graphics windows and plot new borders.

This function can perform both deletion and addition in one call.
The removal of ring borders takes precedence over addition.
}
\examples{
img.path <- system.file("001.png", package = "MtreeRing")

## Read a tree ring image:
t1 <- ring_read(img = img.path, dpi = 1200)

## Split a long core sample into 3 pieces to
## get better display performance and use the
## watershed algorithm to detect ring borders:
t2 <- ring_detect(ring.data = t1, seg = 3, method = 'watershed')

## Do not modify t2, but create a new array object t3. 
## Remove some borders without adding new borders:
t3 <- ring_modify(ring.data = t2, del = c(1, 3, 5, 19:21), add = FALSE)

}
\author{
Jingning Shi
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ring_app_launch.R
\name{ring_app_launch}
\alias{ring_app_launch}
\title{Run Shiny-based Application}
\usage{
ring_app_launch(launch.browser = TRUE)
}
\arguments{
\item{launch.browser}{A logical value. 
If \code{FALSE}, a built-in browser will be launched automatically 
after the app is started. If \code{TRUE}, the system's default 
web browser is used instead. This argument only works for RStudio.
See details below.}
}
\description{
Run a Shiny-based application within the system's default 
web browser.
}
\details{
\code{launch.browser = FALSE} is not recommended, as the file renaming
does not work on the RStudio built-in browser when saving the data.

A workflow for the Shiny app can be found here:
\url{https://ropensci.github.io/MtreeRing/articles/app-MtreeRing.html}. 
Most steps are demonstrated with a gif
to make the workflow more understandable.

To stop the app, go to the R console and press the Escape key. 
You can also click the stop sign icon in the
upper right corner of the RStudio console.
}
\author{
Jingning Shi, Wei Xiang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ring_calculate.R
\name{ring_calculate}
\alias{ring_calculate}
\title{Generate a ring-width series}
\usage{
ring_calculate(ring.data, seriesID)
}
\arguments{
\item{ring.data}{A matrix or array produced by \code{\link{ring_detect}} 
or \code{\link{ring_modify}}.}

\item{seriesID}{A character string specifying the column name of 
the ring-width series.}
}
\value{
A data frame. The series ID is the column name 
and years are row names. The measurements units are millimeters (mm).
}
\description{
This function can calculate ring widths according to 
coordinates of detected ring borders.
}
\examples{
img.path <- system.file("001.png", package = "MtreeRing")

## Read a tree ring image:
t1 <- ring_read(img = img.path, dpi = 1200)

## Split a long core sample into 3 pieces to
## get better display performance and use the
## watershed algorithm to detect ring borders:
t2 <- ring_detect(ring.data = t1, seg = 3, method = 'watershed')

## Calculate ring widths from the attribute list of t2:
rw.df <- ring_calculate(ring.data = t2, seriesID = "940220")
}
\author{
Jingning Shi
}
