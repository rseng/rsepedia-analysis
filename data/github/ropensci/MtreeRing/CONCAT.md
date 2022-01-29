
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
