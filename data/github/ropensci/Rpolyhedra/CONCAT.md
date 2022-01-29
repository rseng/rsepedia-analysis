
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
