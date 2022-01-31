
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rpolyhedra

[![Downloads](http://cranlogs.r-pkg.org/badges/Rpolyhedra?color=brightgreen)](http://www.r-pkg.org/pkg/Rpolyhedra)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/Rpolyhedra?color=brightgreen)](http://www.r-pkg.org/pkg/Rpolyhedra)

<!-- Polyhedra database scraped from publically available sources using R6 objects and 'rgl' visualizing capabilities. -->

This package is a curation made based on the poly package found on
<http://www.netlib.org/polyhedra/> ([Original Help
message](poly_original_help_message.html)), and the polyhedra database
found on <http://dmccooey.com/polyhedra>, both of which provide
polyhedra databases on its own format. As such, Rpolyhedra provides with
the following:

1.  A module to scrape the polyhedra for the different sources found
    with features for incremental correction of issues found and to be
    found in scraping process.
2.  A database of the scraped polyhedra.
3.  An R6 polyhedron representation with ‘rgl’ package visualizing
    capabilities.

| Release                                                                                                  | Usage                                                                                                    | Development                                                                                                                                                                                            |
| :------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [![](https://badges.ropensci.org/157_status.svg)](https://github.com/ropensci/onboarding/issues/157)     | [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.0-blue.svg)](https://cran.r-project.org/) | [![Travis](https://travis-ci.org/ropensci/Rpolyhedra.svg?branch=master)](https://travis-ci.org/ropensci/Rpolyhedra)                                                                                    |
| [![CRAN](http://www.r-pkg.org/badges/version/Rpolyhedra)](https://cran.r-project.org/package=Rpolyhedra) |                                                                                                          | [![codecov](https://codecov.io/gh/ropensci/Rpolyhedra/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/Rpolyhedra)                                                                       |
|                                                                                                          |                                                                                                          | [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) |

# Preview

Through
[Rpolyhedraexplorer](https://qbotics.shinyapps.io/rpolyhedra-explorer/)
you can navigate the polyhedra database without actually installing R
environment.

# How to get started

``` r
install.packages("Rpolyhedra")
```

# How to get started (Development version)

Install the R package using the following commands on the R console:

``` r
devtools::install_github("ropensci/Rpolyhedra", build_opts = NULL)
```

# Loading the database

``` r
library(Rpolyhedra)
```

``` r
# if want to switch to fullDB in user filespace, it will ask you for downloading the full database to your home directory
switchToFullDatabase()
```

# A simple example of 5 regular polyhedra

To get started execute the following commands:

``` r
# 0.  Load libraries
library(knitr)
library(rgl)
# For forarding webgl output to knitr
knit_hooks$set(webgl = hook_webgl)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

``` r

# 1.  Obtain 5 regular solids
polyhedra.2.draw <- getAvailablePolyhedra(source = "netlib")
polyhedra.2.draw <- polyhedra.2.draw %>%
                        filter(scraped.name %in%
                            c("tetrahedron", "octahedron", "cube",
                               "icosahedron", "dodecahedron"))

# 2. Setup colors and scales
n <- nrow(polyhedra.2.draw)
polyhedron.colors <- rainbow(n)

polyhedron.scale <- 5
```

```` r
# For interactive RGL window
#```{r, render, webgl=TRUE}
# 3. open and setup RGL window
open3d()
#> glX 
#>   1
par3d(FOV = 1)
rgl.bg( sphere =FALSE, fogtype = "none", color=c("black"))
rgl.viewpoint(theta = 0, phi=0, zoom=0.8, fov=1)
# 4. for each polyhedron, setup rotation, position and render
for (i in seq_len(n)) {
  # Obtain polyhedron
  polyhedron.row <- polyhedra.2.draw[i,]
  polyhedron.name <- polyhedron.row$scraped.name
  polyhedron <- getPolyhedron(source = polyhedron.row$source, polyhedron.name)

  # Setup angles, position into transformationMatrix
  current.angle <- i/n * 2 * pi
  tm <- rotationMatrix(current.angle, 1, 0, 0)
  x.pos <- round(polyhedron.scale * sin(current.angle), 2)
  y.pos <- round(polyhedron.scale * cos(current.angle), 2)
  tm <- tm %*% translationMatrix(x.pos, y.pos, 0)

  # Render
  print(paste("Drawing ", polyhedron.name, " rotated ", round(current.angle, 2),
              " in (1,0,0) axis. Translated to (", x.pos, ",", y.pos, ",0)",
              " with color ", polyhedron.colors[i], sep = ""))
  shape.rgl <- polyhedron$getRGLModel(transformation.matrix = tm)
  shade3d(shape.rgl, color = polyhedron.colors[i])
}
#> [1] "Drawing tetrahedron rotated 1.26 in (1,0,0) axis. Translated to (4.76,1.55,0) with color #FF0000"
#> [1] "Drawing octahedron rotated 2.51 in (1,0,0) axis. Translated to (2.94,-4.05,0) with color #CCFF00"
#> [1] "Drawing cube rotated 3.77 in (1,0,0) axis. Translated to (-2.94,-4.05,0) with color #00FF66"
#> [1] "Drawing icosahedron rotated 5.03 in (1,0,0) axis. Translated to (-4.76,1.55,0) with color #0066FF"
#> [1] "Drawing dodecahedron rotated 6.28 in (1,0,0) axis. Translated to (0,5,0) with color #CC00FF"

#rgl::rglwidget()
#rgl::rgl.snapshot("man/figures/README-5-polyhedra.png")
````

![5-polyhedra](man/figures/README-5-polyhedra.png)

## sources

### netlib

Includes 142 polyhedra definitions. The PHD format was created to
describe the geometric polyhedron definitions derived mathematically by
Andrew Hume and by the Kaleido program of Zvi Har’El.

PHD files were generated using
[poly2](http://www.netlib.org/poly2/readme) library (no longer
maintained). Although the code is available, specific programming skills
are required to run it.

PHD files can be found in `extdata/www.netlib.org/polyhedra/index.html`

### dmccooey

Includes 767 polyhedra definitions. The [polyhedra
database](http://dmccooey.com/polyhedra/) built by David Mccooey has an
open format which has been scraped to feed Rpolyhedra database

dmccooey files can be found in `extdata/dmccooey.com/polyhedra/`

# Troubleshooting

## devtools

### Ubuntu

``` bash
apt-get install libcurl4-openssl-dev
```

### Windows

run end user CRAN version

### macOS brew

``` bash
brew install openssl
```

After, in R:

``` r
install.packages("devtools")
```

## rgl

### Ubuntu

``` bash
sudo apt-get install r-cran-rgl
```

Please note that the ‘Rpolyhedra’ project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to
this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
Rpolyhedra 0.4.4
============

### MINOR IMPROVEMENTS
* Roxygen changed the way R6 classes are documented
* Problems in testhat with some tests

Rpolyhedra 0.4.2
============

First version published using devtools::release() on [rOpenSci.org](http://www.rOpenSci.org). 

### MINOR IMPROVEMENTS

* Calculates and normalizes polyhedra size using geometry::convhulln instead of bounding box

### BUG FIXES

* A polyhedron now applies internal transformation matrix


Rpolyhedra 0.4.1
============

First version published on [rOpenSci.org](http://www.rOpenSci.org). 

### MINOR IMPROVEMENTS

* Complies with all the prerequisites of rOpenSci and applies the suggestions made by rOpenSci reviewers. 
* Fixes a test that writes on user space.
* Integrates with codecov.io, which allows for better test coverage. 
* Updated examples.

Rpolyhedra 0.4.0
============

### NEW FEATURES

* `Rpolyhedra` can export polyhedra definitions as XML.


Rpolyhedra 0.3.0
============

### NEW FEATURES

* `Rpolyhedra` now has a new database format based on ascii RDSs, which are meant to use less memory, for example when used in a Shiny App.
* `Rpolyhedra` now uses a transformation matrix for general polyhedra manipulation.

### MINOR IMPROVEMENTS

* Applied suggestions from rOpenSci onboarding process. 

# Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at https://github.com/qbotics/Rpolyhedra. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.


## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change. 

Please note we have a code of conduct, please follow it in all your interactions with the project.

## Polyhedra Scraping

In vignettes there is documentation and examples for public functions. However the functionality of the packages goes beyond final user. One major feature has to do with reproducing the scraping process. To that extent, the original polyhedra definitions are shipped (A minimal package version named min-pkg within the package, full-db in a secondary repository accesable with switchToFullDB).

All code is documented within R files. Several functions and classes (for developers) are not included in vignettes to avoid final user confusion, but it should be possible for a developer to get insights of the code following test cases, and extend funcionality or make contributions to the project.

# Reproducibility

The project was built from the ground up with reproducibility in mind. To accommodate for that, each run of the scraping functionality stores information about the run in a Ledger that can be later queried for analytical purposes. 

## Pull Request Process

1. Ensure any install or build dependencies are removed before the end of the layer when doing a 
   build.
2. Update the README.md with details of changes to the interface, this includes new environment 
   variables, exposed ports, useful file locations and container parameters.

## Code style

The code style under use is the recommended in the [google style code (GSC)](https://google.github.io/styleguide/Rguide.xml). 
R6 classes has no code style in GSC. The code style defined for a R6 class-object named FooBar is FooBar.class .




Please ignore previous upload of release 0.4.4 as it doesn't include the sample database. I'm submitting again with the database (It generates automatically but we release a sample already included in the package).

There was some problems running check_rhub() with rgl (because the lack of screen in docker severs, I think), and in some environments it required string which is not a dependence of Rpolyhedra. check_win_devel() built ok in Windows.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  dpi=200,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```



# Rpolyhedra

[![Downloads](http://cranlogs.r-pkg.org/badges/Rpolyhedra?color=brightgreen)](http://www.r-pkg.org/pkg/Rpolyhedra)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/Rpolyhedra?color=brightgreen)](http://www.r-pkg.org/pkg/Rpolyhedra)


 <!-- Polyhedra database scraped from publically available sources using R6 objects and 'rgl' visualizing capabilities. -->


 This package is a curation made based on the poly package found on http://www.netlib.org/polyhedra/ ([Original Help message](poly_original_help_message.html)), and the polyhedra database found on http://dmccooey.com/polyhedra, both of which provide polyhedra databases on its own format. As such, Rpolyhedra provides with the following:

 1. A module to scrape the polyhedra for the different sources found with features for incremental correction of issues found and to be found in scraping process.
 1. A database of the scraped polyhedra.
 1. An R6 polyhedron representation with 'rgl' package visualizing capabilities.


| Release | Usage | Development |
|:--------|:------|:------------|
| [![](https://badges.ropensci.org/157_status.svg)](https://github.com/ropensci/onboarding/issues/157)| [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.0-blue.svg)](https://cran.r-project.org/) | [![Travis](https://travis-ci.org/ropensci/Rpolyhedra.svg?branch=master)](https://travis-ci.org/ropensci/Rpolyhedra) |
| [![CRAN](http://www.r-pkg.org/badges/version/Rpolyhedra)](https://cran.r-project.org/package=Rpolyhedra) | | [![codecov](https://codecov.io/gh/ropensci/Rpolyhedra/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/Rpolyhedra) |
|||[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)|

# Preview
Through [Rpolyhedraexplorer](https://qbotics.shinyapps.io/rpolyhedra-explorer/) you can navigate the polyhedra database without actually installing R environment.



# How to get started
```R
install.packages("Rpolyhedra")
```

# How to get started (Development version)

Install the R package using the following commands on the R console:

```R
devtools::install_github("ropensci/Rpolyhedra", build_opts = NULL)
```

# Loading the database
```{r, rpolyhedra}
library(Rpolyhedra)
```

```R
# if want to switch to fullDB in user filespace, it will ask you for downloading the full database to your home directory
switchToFullDatabase()
```

# A simple example of 5 regular polyhedra

To get started execute the following commands:

```{r, libraries}
# 0.  Load libraries
library(knitr)
library(rgl)
# For forarding webgl output to knitr
knit_hooks$set(webgl = hook_webgl)
library(dplyr)
```

```{r, 1-retrieve}

# 1.  Obtain 5 regular solids
polyhedra.2.draw <- getAvailablePolyhedra(source = "netlib")
polyhedra.2.draw <- polyhedra.2.draw %>%
                        filter(scraped.name %in%
                            c("tetrahedron", "octahedron", "cube",
                               "icosahedron", "dodecahedron"))

# 2. Setup colors and scales
n <- nrow(polyhedra.2.draw)
polyhedron.colors <- rainbow(n)

polyhedron.scale <- 5
```

```{r, render}
# For interactive RGL window
#```{r, render, webgl=TRUE}
# 3. open and setup RGL window
open3d()
par3d(FOV = 1)
rgl.bg( sphere =FALSE, fogtype = "none", color=c("black"))
rgl.viewpoint(theta = 0, phi=0, zoom=0.8, fov=1)
# 4. for each polyhedron, setup rotation, position and render
for (i in seq_len(n)) {
  # Obtain polyhedron
  polyhedron.row <- polyhedra.2.draw[i,]
  polyhedron.name <- polyhedron.row$scraped.name
  polyhedron <- getPolyhedron(source = polyhedron.row$source, polyhedron.name)

  # Setup angles, position into transformationMatrix
  current.angle <- i/n * 2 * pi
  tm <- rotationMatrix(current.angle, 1, 0, 0)
  x.pos <- round(polyhedron.scale * sin(current.angle), 2)
  y.pos <- round(polyhedron.scale * cos(current.angle), 2)
  tm <- tm %*% translationMatrix(x.pos, y.pos, 0)

  # Render
  print(paste("Drawing ", polyhedron.name, " rotated ", round(current.angle, 2),
              " in (1,0,0) axis. Translated to (", x.pos, ",", y.pos, ",0)",
              " with color ", polyhedron.colors[i], sep = ""))
  shape.rgl <- polyhedron$getRGLModel(transformation.matrix = tm)
  shade3d(shape.rgl, color = polyhedron.colors[i])
}

#rgl::rglwidget()
#rgl::rgl.snapshot("man/figures/README-5-polyhedra.png")
```

![5-polyhedra](man/figures/README-5-polyhedra.png)

## sources
### netlib
 Includes 142 polyhedra definitions.
 The PHD format was created to describe the geometric polyhedron definitions derived mathematically by Andrew Hume and by the Kaleido program of Zvi Har'El.

 PHD files were generated using [poly2](http://www.netlib.org/poly2/readme) library (no longer maintained). Although the code is available, specific programming skills are required to run it.

PHD files can be found in `extdata/www.netlib.org/polyhedra/index.html`

### dmccooey
Includes 767 polyhedra definitions.
The [polyhedra database](http://dmccooey.com/polyhedra/) built by David Mccooey has an open format which has been scraped to feed Rpolyhedra database

dmccooey files can be found in `extdata/dmccooey.com/polyhedra/`

# Troubleshooting

## devtools

### Ubuntu

```bash
apt-get install libcurl4-openssl-dev
```

### Windows

run end user CRAN version

### macOS brew

```bash
brew install openssl
```
After, in R:

```R
install.packages("devtools")
```

## rgl

### Ubuntu
```bash
sudo apt-get install r-cran-rgl
```


Please note that the 'Rpolyhedra' project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

---
title: "Rpolyhedra"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Rpolyhedra}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---
# Introduction
This package is a curation made based on the poly package found on http://www.netlib.org/polyhedra/ ([Original Help message](poly_original_help_message.html)), and the polyhedra database found on http://dmccooey.com/polyhedra, both of which provide polyhedra databases on its own format. As such, Rpolyhedra provides with the following:

1. A module to scrape the polyhedra for the different sources found with features for incremental correction of issues found and to be found in scraping process.
1. A database of the scraped polyhedra.
1. An R6 polyhedron representation with 'rgl' package visualizing capabilities. 

# Usage

For final users, the package provides a common interface for accessing public polyhedra databases, analyze properties, compare and visualize them with RGL. 

For advanced users, the package provides a simplified set of R6 objects to scrape and compare polyhedra databases. 
```{r setup, include=FALSE}
library(rgl)
library(dplyr)
library(Rpolyhedra)
setupKnitr()
```


## Get available polyhedra
Once the original files had been processed, a simple call to ```getAvailablePolyhedra()``` retrieves a list of the available polyhedra with properties and status in the polyhedra database:

```{r availablePolyhedra}
#show only the first 10 polyhedra.
head(getAvailablePolyhedra(), n = 10)
```

## Retrieve a polyhedron
The access to a particular polyhedron can be done with a call to ```getPolyhedron(<<source>>, <<polyhedron.name>>)```, which returns a Polyhedron object. For example, to retrieve a cube from the netlib database, the call would be:

```{r getPolyhedron}
cube <- getPolyhedron(source = "netlib", polyhedron.name = "cube")
```

# A demo
To try package functionality, a simple demo can be executed which shows the 5 regular polyhedra.

```{r demo, webgl=TRUE}
# 1.  Obtain 5 regular solids
polyhedra.2.draw <- getAvailablePolyhedra(source = "netlib")
polyhedra.2.draw <- polyhedra.2.draw %>%
                        filter(scraped.name %in%
                            c("tetrahedron", "octahedron", "cube",
                               "icosahedron", "dodecahedron"))

# 2. Setup colors and scales
n <- nrow(polyhedra.2.draw)
polyhedron.colors <- rainbow(n)
polyhedron.scale <- 5

# 3. Open and setup RGL window
open3d()
par3d(FOV = 1)
rgl.bg( sphere =FALSE, fogtype = "none", color=c("black"))
rgl.viewpoint(theta = 0, phi=0, zoom=0.8, fov=1)

# 4. For each polyhedron, setup rotation, position and render
for (i in seq_len(n)) {
  # Obtain polyhedron
  polyhedron.row <- polyhedra.2.draw[i,]
  polyhedron.name <- polyhedron.row$scraped.name
  polyhedron <- getPolyhedron(source = polyhedron.row$source, polyhedron.name)

  # Setup angles, position into transformationMatrix
  current.angle <- i/n * 2 * pi
  tm <- rotationMatrix(current.angle, 1, 0, 0)
  x.pos <- round(polyhedron.scale * sin(current.angle), 2)
  y.pos <- round(polyhedron.scale * cos(current.angle), 2)
  tm <- tm %*% translationMatrix(x.pos, y.pos, 0)

  # Render
  print(paste("Drawing ", polyhedron.name, " rotated ", round(current.angle, 2),
              " in (1,0,0) axis. Translated to (", x.pos, ",", y.pos, ",0)",
              " with color ", polyhedron.colors[i], sep = ""))
  shape.rgl <- polyhedron$getRGLModel(transformation.matrix = tm)
  shade3d(shape.rgl, color = polyhedron.colors[i])
}
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-lib.R
\name{scrapePolyhedra}
\alias{scrapePolyhedra}
\title{Scrape polyhedra objects}
\usage{
scrapePolyhedra(
  scrape.config,
  source.filenames = NULL,
  sources.config = getUserEnvir(".available.sources")
)
}
\arguments{
\item{scrape.config}{predefined configuration for scraping}

\item{source.filenames}{if not null specify which source filenames to scrape}

\item{sources.config}{the sources that will be used by the function}
}
\value{
polyhedra db object
}
\description{
Gets polyhedra objects from text files of
different sources, scheduling and scraping using
predefined configurations.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-lib.R
\name{scrapePolyhedraSources}
\alias{scrapePolyhedraSources}
\title{Scrape polyhedra sources}
\usage{
scrapePolyhedraSources(sources.config =
         getUserEnvir(".available.sources"),
    max.quant.config.schedule = 0,
    max.quant.scrape = 0, time2scrape.source = 30,
    source.filenames = NULL, retry.scrape = FALSE)
}
\arguments{
\item{sources.config}{the sources that will be used by the function}

\item{max.quant.config.schedule}{number of files to schedule}

\item{max.quant.scrape}{number of files scrape}

\item{time2scrape.source}{time applied to scrape source}

\item{source.filenames}{if not null specify which source filenames to scrape}

\item{retry.scrape}{should it retry scrape?}
}
\value{
polyhedra db object
}
\description{
Scrapes polyhedra objects from text files of
different sources, in order to make them available to the
package.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyhedra-lib.R
\docType{class}
\name{PolyhedronStateDeserializer.class}
\alias{PolyhedronStateDeserializer.class}
\title{Polyhedron State Deserializer}
\description{
Polyhedron state for deserialize from database
}
\section{Super class}{
\code{\link[Rpolyhedra:PolyhedronState]{Rpolyhedra::PolyhedronState}} -> \code{PolyhedronStateDeserializer}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{serialized.polyhedron}}{polyhedron definition serialized}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PolyhedronStateDeserializer.class$new()}}
\item \href{#method-scrape}{\code{PolyhedronStateDeserializer.class$scrape()}}
\item \href{#method-clone}{\code{PolyhedronStateDeserializer.class$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="addError">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-addError}{\code{Rpolyhedra::PolyhedronState$addError()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="applyTransformationMatrix">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-applyTransformationMatrix}{\code{Rpolyhedra::PolyhedronState$applyTransformationMatrix()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="buildRGL">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-buildRGL}{\code{Rpolyhedra::PolyhedronState$buildRGL()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="checkEdgesConsistency">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-checkEdgesConsistency}{\code{Rpolyhedra::PolyhedronState$checkEdgesConsistency()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="exportToXML">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-exportToXML}{\code{Rpolyhedra::PolyhedronState$exportToXML()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="getName">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-getName}{\code{Rpolyhedra::PolyhedronState$getName()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="getSolid">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-getSolid}{\code{Rpolyhedra::PolyhedronState$getSolid()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize PolyhedronStateDeserializer object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDeserializer.class$new(serialized.polyhedron)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{serialized.polyhedron}}{a serialized polyhedron}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new  PolyhedronStateDeserializer object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrape"></a>}}
\if{latex}{\out{\hypertarget{method-scrape}{}}}
\subsection{Method \code{scrape()}}{
Generates a PolyhedronStateDefined from a serialized polyhedron
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDeserializer.class$scrape()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A new  PolyhedronStateDefined object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDeserializer.class$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyhedra-lib.R
\name{PolyhedronState.class}
\alias{PolyhedronState.class}
\title{Polyhedron State}
\description{
Polyhedron State

Polyhedron State
}
\details{
This abstract class provide the basis from which every polyhedron state class derivate.
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{source}}{polyhedron definition source}

\item{\code{file.id}}{polyhedron file id}

\item{\code{errors}}{Errors string}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PolyhedronState.class$new()}}
\item \href{#method-addError}{\code{PolyhedronState.class$addError()}}
\item \href{#method-scrape}{\code{PolyhedronState.class$scrape()}}
\item \href{#method-getName}{\code{PolyhedronState.class$getName()}}
\item \href{#method-getSolid}{\code{PolyhedronState.class$getSolid()}}
\item \href{#method-checkEdgesConsistency}{\code{PolyhedronState.class$checkEdgesConsistency()}}
\item \href{#method-applyTransformationMatrix}{\code{PolyhedronState.class$applyTransformationMatrix()}}
\item \href{#method-buildRGL}{\code{PolyhedronState.class$buildRGL()}}
\item \href{#method-exportToXML}{\code{PolyhedronState.class$exportToXML()}}
\item \href{#method-clone}{\code{PolyhedronState.class$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a polyhedronState object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronState.class$new(source, file.id)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{source}}{the source file}

\item{\code{file.id}}{the file id}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new  PolyhedronState object.
'@description
Adds an error to the error string and log it as info
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-addError"></a>}}
\if{latex}{\out{\hypertarget{method-addError}{}}}
\subsection{Method \code{addError()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronState.class$addError(current.error)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{current.error}}{the error to add}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrape"></a>}}
\if{latex}{\out{\hypertarget{method-scrape}{}}}
\subsection{Method \code{scrape()}}{
Scrapes the polyhedra folder files
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronState.class$scrape()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getName"></a>}}
\if{latex}{\out{\hypertarget{method-getName}{}}}
\subsection{Method \code{getName()}}{
Get Polyhedron name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronState.class$getName()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
string with polyhedron name
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getSolid"></a>}}
\if{latex}{\out{\hypertarget{method-getSolid}{}}}
\subsection{Method \code{getSolid()}}{
Returns the object corresponding to the solid
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronState.class$getSolid()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-checkEdgesConsistency"></a>}}
\if{latex}{\out{\hypertarget{method-checkEdgesConsistency}{}}}
\subsection{Method \code{checkEdgesConsistency()}}{
Checks edge consistency
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronState.class$checkEdgesConsistency()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-applyTransformationMatrix"></a>}}
\if{latex}{\out{\hypertarget{method-applyTransformationMatrix}{}}}
\subsection{Method \code{applyTransformationMatrix()}}{
Apply transformation matrix to polyhedron
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronState.class$applyTransformationMatrix(transformation.matrix)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{transformation.matrix}}{the transformation matrix to apply to the polyhedron}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-buildRGL"></a>}}
\if{latex}{\out{\hypertarget{method-buildRGL}{}}}
\subsection{Method \code{buildRGL()}}{
Creates a 'rgl' representation of the object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronState.class$buildRGL(transformation.matrix)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{transformation.matrix}}{the transformation matrix to apply to the polyhedron}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-exportToXML"></a>}}
\if{latex}{\out{\hypertarget{method-exportToXML}{}}}
\subsection{Method \code{exportToXML()}}{
Gets an XML representation out of the polyhedron object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronState.class$exportToXML()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronState.class$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyhedra-lib.R
\docType{class}
\name{Polyhedron.class}
\alias{Polyhedron.class}
\title{Polyhedron}
\description{
Polyhedron container class, which is accessible by the final users upon call
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{file.id}}{Polyhedron file.id}

\item{\code{state}}{Polyhedron state}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Polyhedron.class$new()}}
\item \href{#method-scrapeNetlib}{\code{Polyhedron.class$scrapeNetlib()}}
\item \href{#method-scrapeDmccooey}{\code{Polyhedron.class$scrapeDmccooey()}}
\item \href{#method-deserialize}{\code{Polyhedron.class$deserialize()}}
\item \href{#method-getName}{\code{Polyhedron.class$getName()}}
\item \href{#method-getState}{\code{Polyhedron.class$getState()}}
\item \href{#method-getSolid}{\code{Polyhedron.class$getSolid()}}
\item \href{#method-isChecked}{\code{Polyhedron.class$isChecked()}}
\item \href{#method-getRGLModel}{\code{Polyhedron.class$getRGLModel()}}
\item \href{#method-exportToXML}{\code{Polyhedron.class$exportToXML()}}
\item \href{#method-getErrors}{\code{Polyhedron.class$getErrors()}}
\item \href{#method-checkProperties}{\code{Polyhedron.class$checkProperties()}}
\item \href{#method-clone}{\code{Polyhedron.class$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a polyhedronState object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$new(file.id, state = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{file.id}}{the file id}

\item{\code{state}}{polyhedron state object}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new  Polyhedron object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrapeNetlib"></a>}}
\if{latex}{\out{\hypertarget{method-scrapeNetlib}{}}}
\subsection{Method \code{scrapeNetlib()}}{
scrape Netlib polyhedron definition
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$scrapeNetlib(netlib.p3.lines)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{netlib.p3.lines}}{vector with netlib definition lines}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new  PolyhedronStateDefined object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrapeDmccooey"></a>}}
\if{latex}{\out{\hypertarget{method-scrapeDmccooey}{}}}
\subsection{Method \code{scrapeDmccooey()}}{
scrape Dmccooey polyhedron definition
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$scrapeDmccooey(polyhedra.dmccooey.lines)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{polyhedra.dmccooey.lines}}{vector with Dmccooey definition lines}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new  PolyhedronStateDefined object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize}{}}}
\subsection{Method \code{deserialize()}}{
deserialize a polyhedron state definition
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$deserialize(serialized.polyhedron)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{serialized.polyhedron}}{a serialized version of a polyhedron state}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new  PolyhedronStateDefined object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getName"></a>}}
\if{latex}{\out{\hypertarget{method-getName}{}}}
\subsection{Method \code{getName()}}{
get Polyhedron name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$getName()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
string with polyhedron name
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getState"></a>}}
\if{latex}{\out{\hypertarget{method-getState}{}}}
\subsection{Method \code{getState()}}{
Gets polyhedron state
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$getState()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A new  PolyhedronState object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getSolid"></a>}}
\if{latex}{\out{\hypertarget{method-getSolid}{}}}
\subsection{Method \code{getSolid()}}{
Gets a solid definition
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$getSolid()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A list of vertex vectors composing polyhedron faces.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-isChecked"></a>}}
\if{latex}{\out{\hypertarget{method-isChecked}{}}}
\subsection{Method \code{isChecked()}}{
checks Edges consistency
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$isChecked()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A boolean value
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getRGLModel"></a>}}
\if{latex}{\out{\hypertarget{method-getRGLModel}{}}}
\subsection{Method \code{getRGLModel()}}{
Return an 'rgl' model with an optional transformation described by transformation.matrix parameter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$getRGLModel(transformation.matrix = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{transformation.matrix}}{transformation matrix parameter}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An tmesh3d object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-exportToXML"></a>}}
\if{latex}{\out{\hypertarget{method-exportToXML}{}}}
\subsection{Method \code{exportToXML()}}{
exports an XML definition of current polyhedron
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$exportToXML()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A character object with the XML definition
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getErrors"></a>}}
\if{latex}{\out{\hypertarget{method-getErrors}{}}}
\subsection{Method \code{getErrors()}}{
returns the errors found when processing current polyhedron
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$getErrors()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a data.frame with polyhedron errors
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-checkProperties"></a>}}
\if{latex}{\out{\hypertarget{method-checkProperties}{}}}
\subsection{Method \code{checkProperties()}}{
check properties of current polyhedron
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$checkProperties(expected.vertices, expected.faces)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{expected.vertices}}{expected vertices number}

\item{\code{expected.faces}}{expected faces number}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Unmodified polyhedron object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Polyhedron.class$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyhedra-lib.R
\docType{class}
\name{PolyhedronStateDmccooeyScraper.class}
\alias{PolyhedronStateDmccooeyScraper.class}
\title{Polyhedron State Dmccooey Scraper}
\description{
Polyhedron State Dmccooey Scraper

Polyhedron State Dmccooey Scraper
}
\details{
Scrapes polyhedra from a dmccooey file format
}
\section{Super class}{
\code{\link[Rpolyhedra:PolyhedronState]{Rpolyhedra::PolyhedronState}} -> \code{PolyhedronStateDmccooeyScraper}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{regexp.values.names}}{regexp for scraping values names}

\item{\code{regexp.rn}}{regexp for scraping real numbers}

\item{\code{regexp.values}}{regexp for scraping values}

\item{\code{regexp.vertex}}{regexp for scraping vertices}

\item{\code{regexp.faces}}{regexp for scraping faces}

\item{\code{polyhedra.dmccooey.lines}}{dmccooey polyhedra definition lines}

\item{\code{labels.map}}{labels map where values are}

\item{\code{values}}{labels map where values are}

\item{\code{vertices}}{specification}

\item{\code{vertices.replaced}}{3D values}

\item{\code{faces}}{definition}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PolyhedronStateDmccooeyScraper.class$new()}}
\item \href{#method-setupRegexp}{\code{PolyhedronStateDmccooeyScraper.class$setupRegexp()}}
\item \href{#method-scrapeValues}{\code{PolyhedronStateDmccooeyScraper.class$scrapeValues()}}
\item \href{#method-scrapeVertices}{\code{PolyhedronStateDmccooeyScraper.class$scrapeVertices()}}
\item \href{#method-scrapeFaces}{\code{PolyhedronStateDmccooeyScraper.class$scrapeFaces()}}
\item \href{#method-scrape}{\code{PolyhedronStateDmccooeyScraper.class$scrape()}}
\item \href{#method-getName}{\code{PolyhedronStateDmccooeyScraper.class$getName()}}
\item \href{#method-applyTransformationMatrix}{\code{PolyhedronStateDmccooeyScraper.class$applyTransformationMatrix()}}
\item \href{#method-buildRGL}{\code{PolyhedronStateDmccooeyScraper.class$buildRGL()}}
\item \href{#method-exportToXML}{\code{PolyhedronStateDmccooeyScraper.class$exportToXML()}}
\item \href{#method-clone}{\code{PolyhedronStateDmccooeyScraper.class$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="addError">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-addError}{\code{Rpolyhedra::PolyhedronState$addError()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="checkEdgesConsistency">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-checkEdgesConsistency}{\code{Rpolyhedra::PolyhedronState$checkEdgesConsistency()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="getSolid">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-getSolid}{\code{Rpolyhedra::PolyhedronState$getSolid()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize Dmccooey scraper
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$new(file.id, polyhedra.dmccooey.lines)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{file.id}}{identifier of the definition file.}

\item{\code{polyhedra.dmccooey.lines}}{raw Dmccooey definition file lines}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new  PolyhedronStateDmccooeyScraper object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-setupRegexp"></a>}}
\if{latex}{\out{\hypertarget{method-setupRegexp}{}}}
\subsection{Method \code{setupRegexp()}}{
setupRegexp for Dmccooey definition
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$setupRegexp()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
This PolyhedronStateDmccooeyScraper object with regexp defined.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrapeValues"></a>}}
\if{latex}{\out{\hypertarget{method-scrapeValues}{}}}
\subsection{Method \code{scrapeValues()}}{
scrape values from Dmccooey definition
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$scrapeValues(values.lines)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{values.lines}}{values definitions in Dmccooey source}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
This PolyhedronStateDmccooeyScraper object with values defined.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrapeVertices"></a>}}
\if{latex}{\out{\hypertarget{method-scrapeVertices}{}}}
\subsection{Method \code{scrapeVertices()}}{
scrape polyhedron vertices from definition
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$scrapeVertices(vertices.lines)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{vertices.lines}}{vertices definitions in Dmccooey source}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
This PolyhedronStateDmccooeyScraper object with faces defined.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrapeFaces"></a>}}
\if{latex}{\out{\hypertarget{method-scrapeFaces}{}}}
\subsection{Method \code{scrapeFaces()}}{
scrape polyhedron faces from definition
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$scrapeFaces(faces.lines)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{faces.lines}}{face}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
This PolyhedronStateDmccooeyScraper object with faces defined.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrape"></a>}}
\if{latex}{\out{\hypertarget{method-scrape}{}}}
\subsection{Method \code{scrape()}}{
scrape Dmccooey polyhedron definition
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$scrape()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A new PolyhedronStateDefined object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getName"></a>}}
\if{latex}{\out{\hypertarget{method-getName}{}}}
\subsection{Method \code{getName()}}{
get Polyhedron name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$getName()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
string with polyhedron name
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-applyTransformationMatrix"></a>}}
\if{latex}{\out{\hypertarget{method-applyTransformationMatrix}{}}}
\subsection{Method \code{applyTransformationMatrix()}}{
Apply transformation matrix to polyhedron
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$applyTransformationMatrix(
  transformation.matrix
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{transformation.matrix}}{the transformation matrix to apply to the polyhedron}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-buildRGL"></a>}}
\if{latex}{\out{\hypertarget{method-buildRGL}{}}}
\subsection{Method \code{buildRGL()}}{
Creates a 'rgl' representation of the object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$buildRGL(transformation.matrix)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{transformation.matrix}}{the transformation matrix to apply to the polyhedron}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-exportToXML"></a>}}
\if{latex}{\out{\hypertarget{method-exportToXML}{}}}
\subsection{Method \code{exportToXML()}}{
serializes object in XML
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$exportToXML()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDmccooeyScraper.class$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyhedra-lib.R
\name{PolyhedronStateDefined.class}
\alias{PolyhedronStateDefined.class}
\title{Polyhedron State scraped and defined}
\description{
Polyhedron State scraped and defined

Polyhedron State scraped and defined
}
\section{Super class}{
\code{\link[Rpolyhedra:PolyhedronState]{Rpolyhedra::PolyhedronState}} -> \code{PolyhedronStateDefined}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{file.id}}{polyhedron filename in original}

\item{\code{source}}{polyhedron definition source (netlib|dmccooey)}

\item{\code{name}}{polyhedron name (netlib|dmccooey)}

\item{\code{symbol}}{the eqn(1) input for two symbols separated by a tab;
the Johnson symbol, and the Schlafli symbol (netlib)}

\item{\code{dual}}{the name of the dual polyhedron optionally followed
by a horizontal tab and the number of the dual (netlib)}

\item{\code{sfaces}}{polyhedron solid face list (netlib)}

\item{\code{svertices}}{polyhedron solid vertices list (netlib)}

\item{\code{vertices}}{Polyhedron vertices list (netlib|dmccooey)}

\item{\code{vertices.centered}}{centered vertices for applying
transformation matrices}

\item{\code{net}}{polyhedron 2D net model with vertices defined for
a planar representation (netlib)}

\item{\code{solid}}{polyhedron list of edges which generate a
solid (netlib|dmccooey)}

\item{\code{hinges}}{Polyhedron hinge list (netlib)}

\item{\code{dih}}{Dih attribute (netlib)}

\item{\code{edges}}{polyhedron edges (netlib|dmccooey)}

\item{\code{transformation.matrix}}{transformation matrix for
calculations and visualizing polyhedron}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PolyhedronStateDefined.class$new()}}
\item \href{#method-scrape}{\code{PolyhedronStateDefined.class$scrape()}}
\item \href{#method-getName}{\code{PolyhedronStateDefined.class$getName()}}
\item \href{#method-getSymbol}{\code{PolyhedronStateDefined.class$getSymbol()}}
\item \href{#method-adjustVertices}{\code{PolyhedronStateDefined.class$adjustVertices()}}
\item \href{#method-getVertices}{\code{PolyhedronStateDefined.class$getVertices()}}
\item \href{#method-getNet}{\code{PolyhedronStateDefined.class$getNet()}}
\item \href{#method-getSolid}{\code{PolyhedronStateDefined.class$getSolid()}}
\item \href{#method-inferEdges}{\code{PolyhedronStateDefined.class$inferEdges()}}
\item \href{#method-checkEdgesConsistency}{\code{PolyhedronStateDefined.class$checkEdgesConsistency()}}
\item \href{#method-triangulate}{\code{PolyhedronStateDefined.class$triangulate()}}
\item \href{#method-getConvHull}{\code{PolyhedronStateDefined.class$getConvHull()}}
\item \href{#method-calculateMassCenter}{\code{PolyhedronStateDefined.class$calculateMassCenter()}}
\item \href{#method-getNormalizedSize}{\code{PolyhedronStateDefined.class$getNormalizedSize()}}
\item \href{#method-getTransformedVertices}{\code{PolyhedronStateDefined.class$getTransformedVertices()}}
\item \href{#method-resetTransformationMatrix}{\code{PolyhedronStateDefined.class$resetTransformationMatrix()}}
\item \href{#method-applyTransformationMatrix}{\code{PolyhedronStateDefined.class$applyTransformationMatrix()}}
\item \href{#method-buildRGL}{\code{PolyhedronStateDefined.class$buildRGL()}}
\item \href{#method-exportToXML}{\code{PolyhedronStateDefined.class$exportToXML()}}
\item \href{#method-expectEqual}{\code{PolyhedronStateDefined.class$expectEqual()}}
\item \href{#method-serialize}{\code{PolyhedronStateDefined.class$serialize()}}
\item \href{#method-clone}{\code{PolyhedronStateDefined.class$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="addError">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-addError}{\code{Rpolyhedra::PolyhedronState$addError()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
object initialization routine
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$new(
  source,
  file.id,
  name,
  vertices,
  solid,
  net = NULL,
  symbol = "",
  dual = NULL,
  sfaces = NULL,
  svertices = NULL,
  hinges = NULL,
  dih = NULL,
  normalize.size = TRUE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{source}}{the library to use}

\item{\code{file.id}}{identifier of the definition file.}

\item{\code{name}}{the polyhedron name}

\item{\code{vertices}}{the vertices}

\item{\code{solid}}{the solid object}

\item{\code{net}}{the net}

\item{\code{symbol}}{the symbol}

\item{\code{dual}}{whether it is dual or not}

\item{\code{sfaces}}{the solid faces}

\item{\code{svertices}}{the solid vertices}

\item{\code{hinges}}{the hinges}

\item{\code{dih}}{the dih}

\item{\code{normalize.size}}{whether it has to normalize the size or not}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new  PolyhedronStateDefined object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrape"></a>}}
\if{latex}{\out{\hypertarget{method-scrape}{}}}
\subsection{Method \code{scrape()}}{
scrape polyhedron.
As the state is defined this functions do nothing
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$scrape()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
current object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getName"></a>}}
\if{latex}{\out{\hypertarget{method-getName}{}}}
\subsection{Method \code{getName()}}{
get Polyhedron name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$getName()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
string with polyhedron name
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getSymbol"></a>}}
\if{latex}{\out{\hypertarget{method-getSymbol}{}}}
\subsection{Method \code{getSymbol()}}{
get Polyhedron symbol
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$getSymbol()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
string with polyhedron symbol
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-adjustVertices"></a>}}
\if{latex}{\out{\hypertarget{method-adjustVertices}{}}}
\subsection{Method \code{adjustVertices()}}{
adjust polyhedron Vertices
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$adjustVertices(normalize.size = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{normalize.size}}{whether it has to normalize the size or not}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
modified  PolyhedronStateDefined object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getVertices"></a>}}
\if{latex}{\out{\hypertarget{method-getVertices}{}}}
\subsection{Method \code{getVertices()}}{
Get the polyhedron state
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$getVertices(solid = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{solid}}{toggles the production of solid vertices.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getNet"></a>}}
\if{latex}{\out{\hypertarget{method-getNet}{}}}
\subsection{Method \code{getNet()}}{
Gets the net property
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$getNet()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getSolid"></a>}}
\if{latex}{\out{\hypertarget{method-getSolid}{}}}
\subsection{Method \code{getSolid()}}{
Gets the solid property
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$getSolid()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-inferEdges"></a>}}
\if{latex}{\out{\hypertarget{method-inferEdges}{}}}
\subsection{Method \code{inferEdges()}}{
Infer edges
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$inferEdges(force.recalculation = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{force.recalculation}}{forces the recalculation of the edges}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-checkEdgesConsistency"></a>}}
\if{latex}{\out{\hypertarget{method-checkEdgesConsistency}{}}}
\subsection{Method \code{checkEdgesConsistency()}}{
Checks edges consistency
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$checkEdgesConsistency()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-triangulate"></a>}}
\if{latex}{\out{\hypertarget{method-triangulate}{}}}
\subsection{Method \code{triangulate()}}{
Triangulates the polyhedron
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$triangulate(force = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{force}}{forces the triangulation.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getConvHull"></a>}}
\if{latex}{\out{\hypertarget{method-getConvHull}{}}}
\subsection{Method \code{getConvHull()}}{
Gets the convex hull
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$getConvHull(
  transformation.matrix = self$transformation.matrix,
  vertices.id.3d = private$vertices.id.3d
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{transformation.matrix}}{the transformation matrix}

\item{\code{vertices.id.3d}}{the vertices ids}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
the convex hull
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-calculateMassCenter"></a>}}
\if{latex}{\out{\hypertarget{method-calculateMassCenter}{}}}
\subsection{Method \code{calculateMassCenter()}}{
Calculates the center of mass.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$calculateMassCenter(
  vertices.id.3d = private$vertices.id.3d,
  applyTransformation = TRUE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{vertices.id.3d}}{the vertices ids}

\item{\code{applyTransformation}}{does it need to apply transformations?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getNormalizedSize"></a>}}
\if{latex}{\out{\hypertarget{method-getNormalizedSize}{}}}
\subsection{Method \code{getNormalizedSize()}}{
Gets the normalized size
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$getNormalizedSize(size)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{size}}{the object's size}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getTransformedVertices"></a>}}
\if{latex}{\out{\hypertarget{method-getTransformedVertices}{}}}
\subsection{Method \code{getTransformedVertices()}}{
Gets the transformed vertices
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$getTransformedVertices(
  vertices = self$vertices.centered,
  transformation.matrix = self$transformation.matrix
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{vertices}}{input vertices}

\item{\code{transformation.matrix}}{the transformation matrix}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-resetTransformationMatrix"></a>}}
\if{latex}{\out{\hypertarget{method-resetTransformationMatrix}{}}}
\subsection{Method \code{resetTransformationMatrix()}}{
Resets the transformation matrix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$resetTransformationMatrix()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-applyTransformationMatrix"></a>}}
\if{latex}{\out{\hypertarget{method-applyTransformationMatrix}{}}}
\subsection{Method \code{applyTransformationMatrix()}}{
Apply transformation matrix to polyhedron
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$applyTransformationMatrix(transformation.matrix)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{transformation.matrix}}{the transformation matrix to apply to the polyhedron}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
an applied transformation.matrix
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-buildRGL"></a>}}
\if{latex}{\out{\hypertarget{method-buildRGL}{}}}
\subsection{Method \code{buildRGL()}}{
Build 'rgl'
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$buildRGL(transformation.matrix = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{transformation.matrix}}{the transformation matrix}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-exportToXML"></a>}}
\if{latex}{\out{\hypertarget{method-exportToXML}{}}}
\subsection{Method \code{exportToXML()}}{
Exports the object to XML format
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$exportToXML()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-expectEqual"></a>}}
\if{latex}{\out{\hypertarget{method-expectEqual}{}}}
\subsection{Method \code{expectEqual()}}{
Determines if a polyhedron is equal to this one.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$expectEqual(polyhedron)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{polyhedron}}{the polyhedron to compare to.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize"></a>}}
\if{latex}{\out{\hypertarget{method-serialize}{}}}
\subsection{Method \code{serialize()}}{
Serialize the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$serialize()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateDefined.class$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-lib.R
\name{getAvailablePolyhedra}
\alias{getAvailablePolyhedra}
\title{Get available polyhedra}
\usage{
getAvailablePolyhedra(sources, search.string)
}
\arguments{
\item{sources}{A string vector containing the source, which can be obtained from getAvailableSources().}

\item{search.string}{A search string}
}
\value{
polyhedra names vector
}
\description{
Gets the list of names of available polyhedra and its status in
the polyhedra database, which can be later called with getPolyhedron
}
\examples{

#gets all polyhedra in the database
available.polyhedra <- getAvailablePolyhedra()

#returns all polyhedra from a given source, in this case, netlib
available.netlib.polyhedra <- getAvailablePolyhedra(sources="netlib")

#search within the polyhedron names

cube <- getAvailablePolyhedra(sources="netlib",search.string="cube")
cube
}
\seealso{
getAvailableSources
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-lib.R
\name{getPolyhedraObject}
\alias{getPolyhedraObject}
\title{Get a polyhedra object}
\usage{
getPolyhedraObject()
}
\value{
.polyhedra
}
\description{
Return the polyhedra database handler.
}
\seealso{
PolyhedraDatabase.class
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rpolyhedra-package.R
\docType{package}
\name{Rpolyhedra-package}
\alias{Rpolyhedra-package}
\alias{_PACKAGE}
\alias{Rpolyhedra}
\title{Rpolyhedra: Polyhedra Database}
\description{
A polyhedra database scraped from various sources as R6 objects and 'rgl' visualizing capabilities.
}
\details{
A polyhedra database scraped from:
\itemize{
   \item http://paulbourke.net/dataformats/phd/: PHD files as R6 objects and 'rgl'
         visualizing capabilities. The PHD format was created to describe the geometric
         polyhedra definitions derived mathematically <http://www.netlib.org/polyhedra/>
         by Andrew Hume and by the Kaleido program of Zvi Har'El.
   \item http://dmccooey.com/Polyhedra: Polyhedra text datafiles.
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/Rpolyhedra}
  \item \url{https://github.com/ropensci/Rpolyhedra}
  \item Report bugs at \url{https://github.com/ropensci/Rpolyhedra/issues}
}

}
\author{
\strong{Maintainer}: Alejandro Baranek \email{abaranek@dc.uba.ar} [compiler, copyright holder]

Authors:
\itemize{
  \item Leonardo Belen \email{leobelen@gmail.com} [compiler, copyright holder]
}

Other contributors:
\itemize{
  \item  qbotics \email{qbotics6@gmail.com} [copyright holder]
  \item Barret Schloerke \email{schloerke@gmail.com} [reviewer]
  \item Lijia Yu \email{yu@lijiayu.net} [reviewer]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyhedra-lib.R
\name{PolyhedronStateNetlibScraper.class}
\alias{PolyhedronStateNetlibScraper.class}
\title{Polyhedron State Netlib Scraper}
\description{
Polyhedron State Netlib Scraper

Polyhedron State Netlib Scraper
}
\details{
Scrapes polyhedra from a PHD file format.
}
\section{Super class}{
\code{\link[Rpolyhedra:PolyhedronState]{Rpolyhedra::PolyhedronState}} -> \code{PolyhedronStateNetlibScraper}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{netlib.p3.lines}}{The path to the PHD files}

\item{\code{labels.rows}}{Labels - row of appearance}

\item{\code{labels.map}}{Labels - Map of content}

\item{\code{errors}}{the errors found}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PolyhedronStateNetlibScraper.class$new()}}
\item \href{#method-extractRowsFromLabel}{\code{PolyhedronStateNetlibScraper.class$extractRowsFromLabel()}}
\item \href{#method-getLabels}{\code{PolyhedronStateNetlibScraper.class$getLabels()}}
\item \href{#method-scrapeNet}{\code{PolyhedronStateNetlibScraper.class$scrapeNet()}}
\item \href{#method-extractCFOutBrackets}{\code{PolyhedronStateNetlibScraper.class$extractCFOutBrackets()}}
\item \href{#method-scrapeVertices}{\code{PolyhedronStateNetlibScraper.class$scrapeVertices()}}
\item \href{#method-setupLabelsOrder}{\code{PolyhedronStateNetlibScraper.class$setupLabelsOrder()}}
\item \href{#method-getDataFromLabel}{\code{PolyhedronStateNetlibScraper.class$getDataFromLabel()}}
\item \href{#method-getName}{\code{PolyhedronStateNetlibScraper.class$getName()}}
\item \href{#method-scrape}{\code{PolyhedronStateNetlibScraper.class$scrape()}}
\item \href{#method-applyTransformationMatrix}{\code{PolyhedronStateNetlibScraper.class$applyTransformationMatrix()}}
\item \href{#method-buildRGL}{\code{PolyhedronStateNetlibScraper.class$buildRGL()}}
\item \href{#method-exportToXML}{\code{PolyhedronStateNetlibScraper.class$exportToXML()}}
\item \href{#method-clone}{\code{PolyhedronStateNetlibScraper.class$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="addError">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-addError}{\code{Rpolyhedra::PolyhedronState$addError()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="checkEdgesConsistency">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-checkEdgesConsistency}{\code{Rpolyhedra::PolyhedronState$checkEdgesConsistency()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="Rpolyhedra" data-topic="PolyhedronState" data-id="getSolid">}\href{../../Rpolyhedra/html/PolyhedronState.html#method-getSolid}{\code{Rpolyhedra::PolyhedronState$getSolid()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initializes the object, taking the file.id and PDH file as parameters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$new(file.id, netlib.p3.lines)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{file.id}}{the file id}

\item{\code{netlib.p3.lines}}{the lines to add}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new  PolyhedronStateNetlibScraper object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-extractRowsFromLabel"></a>}}
\if{latex}{\out{\hypertarget{method-extractRowsFromLabel}{}}}
\subsection{Method \code{extractRowsFromLabel()}}{
Extracts data from the label, taking the label number and the
  expected label as parameters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$extractRowsFromLabel(
  label.number,
  expected.label
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{label.number}}{the label number}

\item{\code{expected.label}}{the expected label}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getLabels"></a>}}
\if{latex}{\out{\hypertarget{method-getLabels}{}}}
\subsection{Method \code{getLabels()}}{
get Labels from current netlib file description
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$getLabels()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a list containing labels from netlib file description
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrapeNet"></a>}}
\if{latex}{\out{\hypertarget{method-scrapeNet}{}}}
\subsection{Method \code{scrapeNet()}}{
scrape Net Model from netlib format
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$scrapeNet(net.txt, offset = 0)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{net.txt}}{a vector containing net model in netlib format}

\item{\code{offset}}{in numbering vertices}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a list containing a net model
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-extractCFOutBrackets"></a>}}
\if{latex}{\out{\hypertarget{method-extractCFOutBrackets}{}}}
\subsection{Method \code{extractCFOutBrackets()}}{
Remove brackets for current field content
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$extractCFOutBrackets(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{a string containing brackets}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
value
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrapeVertices"></a>}}
\if{latex}{\out{\hypertarget{method-scrapeVertices}{}}}
\subsection{Method \code{scrapeVertices()}}{
scrape vertices described in netlib format
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$scrapeVertices(vertices.txt)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{vertices.txt}}{vector containing netlib format vertices}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
data.frame containing netlib vertices
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-setupLabelsOrder"></a>}}
\if{latex}{\out{\hypertarget{method-setupLabelsOrder}{}}}
\subsection{Method \code{setupLabelsOrder()}}{
setupLabelsOrder
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$setupLabelsOrder()}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{vertices.txt}}{vector containing netlib format vertices}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
data.frame containing netlib vertices
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getDataFromLabel"></a>}}
\if{latex}{\out{\hypertarget{method-getDataFromLabel}{}}}
\subsection{Method \code{getDataFromLabel()}}{
Get data from label specified as parameter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$getDataFromLabel(label)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{label}}{the label to get data from}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
value
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getName"></a>}}
\if{latex}{\out{\hypertarget{method-getName}{}}}
\subsection{Method \code{getName()}}{
get Polyhedron name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$getName()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
string with polyhedron name
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrape"></a>}}
\if{latex}{\out{\hypertarget{method-scrape}{}}}
\subsection{Method \code{scrape()}}{
scrape Netlib polyhedron definition
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$scrape()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A new PolyhedronStateDefined object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-applyTransformationMatrix"></a>}}
\if{latex}{\out{\hypertarget{method-applyTransformationMatrix}{}}}
\subsection{Method \code{applyTransformationMatrix()}}{
Apply transformation matrix to polyhedron
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$applyTransformationMatrix(
  transformation.matrix
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{transformation.matrix}}{the transformation matrix to apply to the polyhedron}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-buildRGL"></a>}}
\if{latex}{\out{\hypertarget{method-buildRGL}{}}}
\subsection{Method \code{buildRGL()}}{
Creates a 'rgl' representation of the object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$buildRGL(transformation.matrix)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{transformation.matrix}}{the transformation matrix to apply to the polyhedron}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-exportToXML"></a>}}
\if{latex}{\out{\hypertarget{method-exportToXML}{}}}
\subsection{Method \code{exportToXML()}}{
serializes object in XML
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$exportToXML()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedronStateNetlibScraper.class$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serialization-lib.R
\name{polyhedronToXML}
\alias{polyhedronToXML}
\title{Polyhedron to XML}
\usage{
polyhedronToXML(polyhedron.state.defined, is.transformed.vertices = TRUE)
}
\arguments{
\item{polyhedron.state.defined}{the polyhedron to get a representation from}

\item{is.transformed.vertices}{flag which states if vertices are in original position or transformationMatrix applied}
}
\value{
an XML document, ready to be converted to String with XML::saveXML()
}
\description{
Gets an XML representation out of the polyhedron object
}
\examples{
#get the representation of a cube (netlib library)
XML::saveXML(polyhedronToXML(getPolyhedron("netlib", "cube")$state))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-lib.R
\name{getPolyhedron}
\alias{getPolyhedron}
\title{Get polyhedron}
\usage{
getPolyhedron(source = "netlib", polyhedron.name)
}
\arguments{
\item{source}{string vector, which can be obtained from getAvailableSources()}

\item{polyhedron.name}{a valid name of a polyhedron in
the database. Current names can be found with getAvailablePolyhedra()}
}
\value{
polyhedron R6 object
}
\description{
Gets a polyhedron from the database. It returns an R6 Class
with all its characteristics and functions.
The object returned, of type Polyhedron.class, allows to the
user to get access to all the functionality provided.
}
\examples{
tetrahedron <- getPolyhedron(source = 'netlib',
       polyhedron.name = 'tetrahedron')

# returns name of polyhedra
tetrahedron$getName()

# polyhedron state
tetrahedron.state <- tetrahedron$getState()

# Johnson symbol and Schlafli symbol
tetrahedron.state$getSymbol()

# vertex data.frame
tetrahedron.state$getVertices()

# List of faces of solid representation (3D)
tetrahedron.state$getSolid()

# List of faces of net representation (2D)
tetrahedron.state$getNet()
}
\seealso{
getAvailablePolyhedra, getAvailableSources
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-lib.R
\name{switchToFullDatabase}
\alias{switchToFullDatabase}
\title{Switch to full database}
\usage{
switchToFullDatabase(env=NA)
}
\arguments{
\item{env}{The environment to run on, can be PACKAGE,
HOME or NA. If NA, it asks the user for a an Environment.}
}
\value{
.data.env
}
\description{
Prompts user for changing database to fulldb in
user filespace. Also, allows the user to switch back to
the package database, which is a minimal one for testing purposes.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db-lib.R
\docType{class}
\name{PolyhedraDatabase.class}
\alias{PolyhedraDatabase.class}
\title{Polyhedra database}
\description{
Scrapes all polyhedra in data folder to save a representation which
is accessible by the final users upon call to \code{getPolyhedron()}.
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{version}}{version of database file}

\item{\code{polyhedra.rds.file}}{path of rds database file}

\item{\code{sources.config}}{Sources configuration for scraping different sources}

\item{\code{ledger}}{rr ledger of scraping process}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PolyhedraDatabase.class$new()}}
\item \href{#method-getVersion}{\code{PolyhedraDatabase.class$getVersion()}}
\item \href{#method-configPolyhedraRDSPath}{\code{PolyhedraDatabase.class$configPolyhedraRDSPath()}}
\item \href{#method-existsSource}{\code{PolyhedraDatabase.class$existsSource()}}
\item \href{#method-addSourceConfig}{\code{PolyhedraDatabase.class$addSourceConfig()}}
\item \href{#method-existsPolyhedron}{\code{PolyhedraDatabase.class$existsPolyhedron()}}
\item \href{#method-getPolyhedraSourceDir}{\code{PolyhedraDatabase.class$getPolyhedraSourceDir()}}
\item \href{#method-getPolyhedronFilename}{\code{PolyhedraDatabase.class$getPolyhedronFilename()}}
\item \href{#method-getPolyhedron}{\code{PolyhedraDatabase.class$getPolyhedron()}}
\item \href{#method-addPolyhedron}{\code{PolyhedraDatabase.class$addPolyhedron()}}
\item \href{#method-configPolyhedraSource}{\code{PolyhedraDatabase.class$configPolyhedraSource()}}
\item \href{#method-saveRDS}{\code{PolyhedraDatabase.class$saveRDS()}}
\item \href{#method-cover}{\code{PolyhedraDatabase.class$cover()}}
\item \href{#method-scrape}{\code{PolyhedraDatabase.class$scrape()}}
\item \href{#method-testRR}{\code{PolyhedraDatabase.class$testRR()}}
\item \href{#method-generateTestTasks}{\code{PolyhedraDatabase.class$generateTestTasks()}}
\item \href{#method-schedulePolyhedraSources}{\code{PolyhedraDatabase.class$schedulePolyhedraSources()}}
\item \href{#method-getAvailableSources}{\code{PolyhedraDatabase.class$getAvailableSources()}}
\item \href{#method-getAvailablePolyhedra}{\code{PolyhedraDatabase.class$getAvailablePolyhedra()}}
\item \href{#method-clone}{\code{PolyhedraDatabase.class$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new PolyhedraDatabase object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$new()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A new `PolyhedraDatabase` object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getVersion"></a>}}
\if{latex}{\out{\hypertarget{method-getVersion}{}}}
\subsection{Method \code{getVersion()}}{
get the version of the current object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$getVersion()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Database version
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-configPolyhedraRDSPath"></a>}}
\if{latex}{\out{\hypertarget{method-configPolyhedraRDSPath}{}}}
\subsection{Method \code{configPolyhedraRDSPath()}}{
sets the path of the RDS object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$configPolyhedraRDSPath()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Database version
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-existsSource"></a>}}
\if{latex}{\out{\hypertarget{method-existsSource}{}}}
\subsection{Method \code{existsSource()}}{
Determines if the source exists on
  the database
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$existsSource(source)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{source}}{source description}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
boolean value
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-addSourceConfig"></a>}}
\if{latex}{\out{\hypertarget{method-addSourceConfig}{}}}
\subsection{Method \code{addSourceConfig()}}{
add  source.config to the database
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$addSourceConfig(source.config)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{source.config}}{SourceConfig object able to scrape source polyhedra definitions}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
PolyhedraDatabase.class object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-existsPolyhedron"></a>}}
\if{latex}{\out{\hypertarget{method-existsPolyhedron}{}}}
\subsection{Method \code{existsPolyhedron()}}{
Determines if the database includes a polyhedron which name
matches the parameter value
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$existsPolyhedron(source = "netlib", polyhedron.name)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{source}}{source description}

\item{\code{polyhedron.name}}{polyhedron description}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
boolean value
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getPolyhedraSourceDir"></a>}}
\if{latex}{\out{\hypertarget{method-getPolyhedraSourceDir}{}}}
\subsection{Method \code{getPolyhedraSourceDir()}}{
gets polyhedra sources folder
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$getPolyhedraSourceDir(source, create.dir = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{source}}{source description}

\item{\code{create.dir}}{if dir does not exists, create it}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
string with polyhedra sources path
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getPolyhedronFilename"></a>}}
\if{latex}{\out{\hypertarget{method-getPolyhedronFilename}{}}}
\subsection{Method \code{getPolyhedronFilename()}}{
gets the filename of the polyhedron matching parameter.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$getPolyhedronFilename(
  source,
  polyhedron.name,
  extension
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{source}}{source description}

\item{\code{polyhedron.name}}{polyhedron description}

\item{\code{extension}}{extension of the polyhedron filename}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
string with polyhedron filename
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getPolyhedron"></a>}}
\if{latex}{\out{\hypertarget{method-getPolyhedron}{}}}
\subsection{Method \code{getPolyhedron()}}{
gets polyhedron object which name
matches the parameter value
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$getPolyhedron(
  source = "netlib",
  polyhedron.name,
  strict = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{source}}{source description}

\item{\code{polyhedron.name}}{polyhedron description}

\item{\code{strict}}{halts execution if polyhedron not found}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Polyhedron.class object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-addPolyhedron"></a>}}
\if{latex}{\out{\hypertarget{method-addPolyhedron}{}}}
\subsection{Method \code{addPolyhedron()}}{
add polyhedron object to the database
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$addPolyhedron(
  source = "netlib",
  source.filename,
  polyhedron,
  overwrite = FALSE,
  save.on.change = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{source}}{source description}

\item{\code{source.filename}}{filename of the polyhedron source definition}

\item{\code{polyhedron}}{polyhedron object}

\item{\code{overwrite}}{overwrite exiting definition}

\item{\code{save.on.change}}{saves Database state after operation}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Polyhedron.class object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-configPolyhedraSource"></a>}}
\if{latex}{\out{\hypertarget{method-configPolyhedraSource}{}}}
\subsection{Method \code{configPolyhedraSource()}}{
Process parameter filenames using source.config parameter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$configPolyhedraSource(
  source.config,
  source.filenames = NULL,
  max.quant = 0,
  save.on.change = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{source.config}}{source configuration for scraping files}

\item{\code{source.filenames}}{filenames of the polyhedron source definition}

\item{\code{max.quant}}{maximum filenames to process}

\item{\code{save.on.change}}{saves Database state after operation}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Modified `PolyhedraDatabase` object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-saveRDS"></a>}}
\if{latex}{\out{\hypertarget{method-saveRDS}{}}}
\subsection{Method \code{saveRDS()}}{
saveRDS
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$saveRDS(save.on.change = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{save.on.change}}{saves Database state after operation}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
saveRDS return status
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cover"></a>}}
\if{latex}{\out{\hypertarget{method-cover}{}}}
\subsection{Method \code{cover()}}{
Cover objects and applies covering.code parameter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$cover(
  mode,
  sources = names(self$sources.config),
  covering.code,
  polyhedra.names = NULL,
  max.quant = 0,
  save.on.change = FALSE,
  seed = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mode}}{covering mode. Available values are "scrape.queued", "scrape.retry","skipped",  "test"}

\item{\code{sources}}{sources names}

\item{\code{covering.code}}{code for applying in covering}

\item{\code{polyhedra.names}}{polyhedra names to cover (optional)}

\item{\code{max.quant}}{maximum numbers of polyhedra to cover}

\item{\code{save.on.change}}{saves Database state after operation}

\item{\code{seed}}{seed for deterministic random generator}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list with resulting objects covered
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-scrape"></a>}}
\if{latex}{\out{\hypertarget{method-scrape}{}}}
\subsection{Method \code{scrape()}}{
Scrape polyhedra queued sources
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$scrape(
  mode = "scrape.queued",
  sources = names(self$sources.config),
  max.quant = 0,
  time2scrape.source = 30,
  save.on.change = FALSE,
  skip.still.queued = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mode}}{covering mode. Available values are "scrape.queued", "scrape.retry","skipped",  "test"}

\item{\code{sources}}{sources names}

\item{\code{max.quant}}{maximum numbers of polyhedra to cover}

\item{\code{time2scrape.source}}{maximum time to spend scraping each source}

\item{\code{save.on.change}}{saves Database state after operation}

\item{\code{skip.still.queued}}{Flag unscraped files with status `skipped``}

\item{\code{covering.code}}{code for applying in covering}

\item{\code{polyhedra.names}}{polyhedra names to cover (optional)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list with resulting objects covered
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-testRR"></a>}}
\if{latex}{\out{\hypertarget{method-testRR}{}}}
\subsection{Method \code{testRR()}}{
testRR
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$testRR(
  sources = names(self$sources.config),
  max.quant = 0
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{sources}}{sources names}

\item{\code{max.quant}}{maximum numbers of polyhedra to cover}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list with resulting objects tested
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-generateTestTasks"></a>}}
\if{latex}{\out{\hypertarget{method-generateTestTasks}{}}}
\subsection{Method \code{generateTestTasks()}}{
generate Test tasks for selected polyhedra
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$generateTestTasks(
  sources = names(self$sources.config),
  polyhedra.names = NULL,
  TestTaskClass,
  max.quant = 0
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{sources}}{sources names}

\item{\code{polyhedra.names}}{polyhedra names to cover (optional)}

\item{\code{TestTaskClass}}{an R6 TestTaskClass class}

\item{\code{max.quant}}{maximum numbers of polyhedra to cover}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list with resulting TestTasks generated
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-schedulePolyhedraSources"></a>}}
\if{latex}{\out{\hypertarget{method-schedulePolyhedraSources}{}}}
\subsection{Method \code{schedulePolyhedraSources()}}{
Schedules polyhedra sources for scraping
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$schedulePolyhedraSources(
  sources.config = getPackageEnvir(".available.sources"),
  source.filenames = NULL,
  max.quant = 0,
  save.on.change = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{sources.config}}{sources configurations for scraping files}

\item{\code{source.filenames}}{filenames of the polyhedron source definition}

\item{\code{max.quant}}{maximum filenames to process}

\item{\code{save.on.change}}{saves Database state after operation}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Modified `PolyhedraDatabase` object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getAvailableSources"></a>}}
\if{latex}{\out{\hypertarget{method-getAvailableSources}{}}}
\subsection{Method \code{getAvailableSources()}}{
Returns available sources in current database
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$getAvailableSources()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A vector with names of available sources
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getAvailablePolyhedra"></a>}}
\if{latex}{\out{\hypertarget{method-getAvailablePolyhedra}{}}}
\subsection{Method \code{getAvailablePolyhedra()}}{
Retrieves all polyhedron within the source those names match with search.string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$getAvailablePolyhedra(
  sources = self$getAvailableSources(),
  search.string = NULL,
  ignore.case = TRUE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{sources}}{sources names}

\item{\code{search.string}}{string for matching polyhedron names}

\item{\code{ignore.case}}{ignore case in search string}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list with resulting objects covered
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolyhedraDatabase.class$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-lib.R
\name{getAvailableSources}
\alias{getAvailableSources}
\title{Get available sources}
\usage{
getAvailableSources()
}
\value{
sources string vector, which can be obtained from getAvailableSources()
}
\description{
Gets the list of names of available sources in database to be used later as
references to the package.
}
\examples{
#gets all sources in the database
available.sources <- getAvailableSources()

#returns all polyhedra from all sources
available.polyhedra <- getAvailablePolyhedra(sources=available.sources)

#search within the polyhedron names from all sources
cubes <- getAvailablePolyhedra(sources=available.sources,
        search.string="cube")
cubes
}
\seealso{
getAvailablePolyhedra, getPolyhedron
}
