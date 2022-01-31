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
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Vizumap <img src='man/figures/Vizumap_Hex.png' align="right" height="138.5" />

## An R package for visualizing uncertainty in spatial data.

[![R build
status](https://github.com/lydialucchesi/Vizumap/workflows/R-CMD-check/badge.svg)](https://github.com/lydialucchesi/Vizumap/actions)

There is a [`Vizumap` pkgdown
site](https://lydialucchesi.github.io/Vizumap/) with a vignette.

A [`Vizumap` paper](https://doi.org/10.21105/joss.02409) is available in
the Journal of Open Source Software (JOSS). If you use `Vizumap`, please
cite this paper.

## Installation

You can install a development version of the `Vizumap` package using the
command below.

    remotes::install_github(repo = "lydialucchesi/Vizumap", build_vignettes = TRUE, force = TRUE)

## Authors

Lydia Lucchesi, Australian National University & CSIRO Data61, Email:
<Lydia.Lucchesi@anu.edu.au>

Petra Kuhnert, CSIRO Data61, Email: <Petra.Kuhnert@data61.csiro.au>

## About the package

Approaches for visualising uncertainty in spatial data are presented in
this package. These include the three approaches developed in [Lucchesi
and Wikle
(2017)](http://faculty.missouri.edu/~wiklec/LucchesiWikle2017Stat) and a
fourth approach presented in [Kuhnert et
al. (2018)](https://publications.csiro.au/publications/#publication/PIcsiro:EP168206).
The package is outlined in [Lucchesi et
al. (2021)](https://doi.org/10.21105/joss.02409).

#### Bivariate Maps

In these bivariate choropleth maps, two colour schemes, one representing
the estimates and another representing the margins of error, are blended
so that an estimate and its error can be conveyed on a map using a
single colour.

#### Map Pixelation

In this approach, each map region is pixelated. Pixels are filled with
colours representing values within an estimate’s margin of error.
Regions that appear as a solid colour reflect smaller margins of error,
while more pixelated regions indicate greater uncertainty. These maps
can be animated to provide a novel uncertainty visualisation experience.

#### Glyph Rotation

In this method, glyphs located at region centroids are rotated to
represent uncertainty. The colour filling each glyph corresponds to the
estimate.

#### Exceedance Probability Maps

The final map-based exploration is through exceedance probabilities,
which are visualised on a map to highlight regions that exhibit varying
levels of departure from a threshold of concern or target.

## Examples

A vignette for the `Vizumap` package is available and contains examples
relating to each of the visualisation methods.

    vignette("Vizumap")

## Testing

If you would like to install and run the unit tests interactively,
include `INSTALL_opts = "--install-tests"` in the installation code.

    remotes::install_github(repo = "lydialucchesi/Vizumap", build_vignettes = TRUE, force = TRUE, INSTALL_opts = "--install-tests")
    
    testthat::test_package("Vizumap", reporter = "stop")

## Contribute

To contribute to `Vizumap`, please follow these
[guidelines](CONTRIBUTING.md).

Please note that the `Vizumap` project is released with a [Contributor
Code of Conduct](CONDUCT.md). By contributing to this project, you agree
to abide by its terms.

## License

`Vizumap` version 1.2.0 is licensed under [GPLv3](LICENSE.md).

## Citation

Lucchesi et al., (2021). Vizumap: an R package for visualising
uncertainty in spatial data. Journal of Open Source Software, 6(59),
2409, <https://doi.org/10.21105/joss.02409>

    @article{lucchesi2021vizumap,
      title={Vizumap: an R package for visualising uncertainty in spatial data},
      author={Lucchesi, Lydia R and Kuhnert, Petra M and Wikle, Christopher K},
      journal={Journal of Open Source Software},
      volume={6},
      number={59},
      pages={2409},
      year={2021}
    }

## History of Vizumap

Vizumap began as a visualisation project at the University of Missouri
in 2016. Chris Wikle, professor of statistics, posed an interesting
research question to Lydia Lucchesi, a student curious about data
visualisation and R.

How do you include uncertainty on a map displaying areal data estimates?

Over the course of a year, they put together three methods for
visualising uncertainty in spatial statistics: the bivariate choropleth
map, the pixel map, and the glyph map. By mid-2017, there were maps, and
there was a lot of R code, but there was not a tool that others could
use to easily make these types of maps, too. That’s when statistician
Petra Kuhnert recommended developing an R package. Over the course of a
month, Petra and Lydia developed Vizumap (originally named VizU) at
CSIRO Data61 in Canberra, Australia. Since then, the package has been
expanded to include exceedance probability maps, an uncertainty
visualisation method developed by Petra while working on a Great Barrier
Reef (GBR) project.

Vizumap has been used to visualise the uncertainty of American Community
Survey estimates, the prediction errors of sediment estimates in a GBR
catchment, and most recently the [uncertainty of estimated locust
densities in
Australia](https://www.nature.com/articles/s41598-020-73897-1/figures/4).
We would like to assemble a Vizumap gallery that showcases different
applications of the package’s mapping methods. If you use Vizumap to
visualise uncertainty, please feel free to send the map our way. We
would like to see it\!

## References

Kuhnert, P.M., Pagendam, D.E., Bartley, R., Gladish, D.W., Lewis, S.E.
and Bainbridge, Z.T. (2018) [Making management decisions in face of
uncertainty: a case study using the Burdekin catchment in the Great
Barrier Reef, Marine and Freshwater
Research](https://publications.csiro.au/publications/#publication/PIcsiro:EP168206),
69, 1187-1200, <https://doi.org/10.1071/MF17237>.

Lucchesi, L.R. and Wikle C.K. (2017) [Visualizing uncertainty in areal
data with bivariate choropleth maps, map pixelation and glyph
rotation](http://faculty.missouri.edu/~wiklec/LucchesiWikle2017Stat),
Stat, <https://doi.org/10.1002/sta4.150>.

Lucchesi, L.R., Kuhnert, P.M. and Wikle, C.K. (2021) [Vizumap: an R
package for visualising uncertainty in spatial
data](https://doi.org/10.21105/joss.02409), Journal of Open Source
Software, <https://doi.org/10.21105/joss.02409>.
## Vizumap v1.2.0

---

- JOSS revisions

## Vizumap v1.1.0

---

- updates since first release in 2017
- package renamed in November 2018

## VizU v1.0.0

---

- initial package release

# Contributing

Thank you for your cooperation and contribution! We appreciate your help making Vizumap the best it can be.

## How to Propose a Change to Vizumap

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it is a problem. If you have found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

Guidelines for contributors:

* Create a named branch (not master) which merges cleanly with master.
* Limit the scope of the pull request to a clearly defined improvement or fix.
* Follow the tidyverse [style guide](https://style.tidyverse.org). Please do not restyle code unrelated to your PR.
* Use [roxygen2](https://cran.r-project.org/package=roxygen2) for documentation.
* Use [testthat](https://cran.r-project.org/package=testthat). Contributions with test cases included are easier to accept.
* For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes you made followed by your
GitHub username and links to relevant issue(s)/PR(s).
* After you merge a pull request, delete the associated branch if you can. You will see an option for this on the GitHub pull request page.

Every pull request will be reviewed by @lydialucchesi or @pkuhnert. They will either merge it or provide feedback on issues to be addressed. They may add additional commits to the pull request. @lydialucchesi and @pkuhnert can make minor edits to a pull request (e.g., fixing typos) before merging it.

### Code of Conduct

Please note that the Vizumap project is released with a [Contributor Code of Conduct](CONDUCT.md). By contributing to this project, you agree to abide by its terms.

---
title: 'Vizumap: an R package for visualising uncertainty in spatial data'
tags:
  - R
  - spatial statistics
  - visualisation
  - uncertainty
  - maps
authors:
  - name: Lydia R. Lucchesi
    affiliation: "1, 2"
  - name: Petra M. Kuhnert
    affiliation: 2
  - name: Christopher K. Wikle
    affiliation: 3
affiliations:
 - name: Australian National University, Canberra, Australia
   index: 1
 - name: CSIRO Data61, Canberra, Australia
   index: 2
 - name: University of Missouri, Columbia, USA
   index: 3 
date: 3 June 2020
bibliography: paper.bib
---

# Summary

To make a sound data-driven decision, it is necessary to consider the quality and strength of the evidence. Therefore, it is important that the uncertainty of statistical estimates be effectively communicated to decision makers to ensure it is properly considered in the decision-making 
process. Generally, this uncertainty information is shared through visualisation (e.g., error bars). However, in spatial applications, finding methods that can communicate additional information, about the spatial estimates, in an understandable and meaningful way can be challenging. To address this visualisation shortcoming in spatial statistics, we developed the `Vizumap` R package. It is a toolkit designed for statisticians, scientists, data journalists, etc., discussing uncertainty in spatial data.

`Vizumap` contains a series of straightforward functions that can be used to create four different types of maps, which are discussed in detail in @vizMethod and @gbrData. The visualisation methods include the bivariate choropleth map, pixel map, glyph map, and exceedance probability map. Following is a brief description of each method.

* For the bivariate choropleth map, two colour scales are merged to create a bivariate grid that can encode the estimates and errors at once, making it possible to fill each region with a single colour that represents both values.

* The pixel map is generated by dividing each region into small pixels, which are then filled with a colour value sampled from the colour gradient range. The method for assigning colours to pixels depends on the type of distribution specified. Areas of greater uncertainty appear pixelated while areas of less uncertainty appear smooth, as though filled by only a single colour. These pixel maps can be animated so that areas of high uncertainty flicker while areas of low uncertainty appear static to the eye.

* The glyph map is intended to address the issue of unequal region sizes on a choropleth map. Glyphs of equal size are placed at the centroid of each region, filled with a colour that represents the estimate, and rotated depending on the degree of uncertainty.

* Exceedance probability maps show the probability of exceeding a nominated threshold of concern. These probabilities can be pre-calculated and passed to the exceedance probability map function or calculated within the function given a prescribed probability distribution function, estimate, and error. For example, let $X$ denote a random variable, and $x$ is a possible value of that random variable, $X$. $F_x(x)$ is the cumulative distribution function and equals the probability that the value of $X$ is less than or equal to a specific value or threshold, $x$, so 
\begin{equation}\label{eq:eq1}
{F_x(x) = Pr[X {\leq} x].}
\end{equation}
The probability that the value $X$ is greater than the specific threshold $x$ can then be written as 
\begin{equation}\label{eq:eq2}
{Pr[X > x] = 1 - F_x(x) = 1 - Pr[X {\leq} x].}
\end{equation}
In R, assuming an exponential distribution, this can be evaluated through the following expressions: `1 - pexp(q, rate)` or `pexp(q, rate, lower.tail = FALSE)`.

Other related R packages include `pixelate` and `biscale`. The `pixelate` package explores the use of pixelation for uncertainty visualisation on isopleth maps [@pixelatePackage], and the `biscale` package can be used to generate bivariate maps for two variables [@biscalePackage].

A comprehensive vignette that demonstrates how to use `Vizumap` is available after package installation. The functions are divided into three categories: formatting, building, and viewing. Formatting functions prepare data for use in the building functions, which are used to build the colour palettes, maps, and map keys. Viewing functions are used to check and combine the different graphical objects designed with the building functions. Previous applications of the `Vizumap` R package include visualising American Community Survey estimates with their corresponding margins of error. We believe `Vizumap` is useful in a wide range of applications, and we will continue to improve the toolkit to enhance its utility.

As an illustration below, we use `Vizumap` to visualise estimated pollutant loads of sediment from the upper Burdekin catchment in Queensland, Australia, into the Great Barrier Reef (GBR). The predictions and uncertainties, published in @gbrData, were developed from a Bayesian Hierarchical Model (BHM) that assimilated estimates of sediment concentration and flow with modelled output from a catchment model developed on the upper Burdekin catchment. Development of this modelling strategy is discussed in @gbrMethod. Here we just focus on total suspended sediment (TSS).  The export of pollutants from coastal catchments within Australia has important implications for the health of the GBR lagoon, and `Vizumap` offers a variety of methods for communicating these predictions and uncertainties to catchment managers and policy makers.

## Bivariate map

Figure \ref{bivMap} is a bivariate map of the upper Burdekin catchment in Queensland, Australia. It visualises predicted total suspended sediment (TSS) concentrations and prediction uncertainty. In this example, a custom colour palette is created and then used to build the bivariate map and colour key. The bivariate bins are defined using terciles.

![Bivariate map and bivariate key for the upper Burdekin catchment in Queensland, Australia, showing the total suspended sediment (TSS) concentration in mg/L.\label{bivMap}](bivariateMap.png){width=4in height=4in}

## Pixel map

Figure \ref{pixMap} visualises the uncertainty of TSS predictions, while giving a general idea of estimated TSS concentrations. For a closer look at the pixelation, a subset of regions is included to the right of the map (Figure \ref{pixMap}(B)). The colours filling the pixels in each region were sampled from the estimate's relative frequency distribution. This pixel map can be animated so that the pixels flicker between sampled values. If the map below were to be animated, the areas that appear most pixelated in the static map would, correspondingly, have the most visible movement among pixels in the animated map. Movement among pixels in areas of low uncertainty would be hard to detect due to the minimal differences between the sampled values of orange. An example of an animated pixel map can be found in @vizMethod.

![Pixel map showing the TSS concentrations for the upper Burdekin catchment in Queensland, Australia. Figure \ref{pixMap}(B) provides a closer look at the pixelation in five regions from Figure \ref{pixMap}(A). \label{pixMap}](pixelMap.png){width=4in height=4in}

## Glyph map

Figure \ref{glyphMap} depicts estimated TSS concentrations and the uncertainty of these predictions. The colour filling each glyph (located at the region centroid) represents the estimate, and the rotation of the glyph represents the uncertainty.

![Glyph map and glyph key showing the TSS predictions and uncertainties for the upper Burdekin catchment in Queensland, Australia.\label{glyphMap}](glyphMap.png){width=4in height=4in}

## Exceedance probability map

In Figure \ref{exceedMap}, the calculated probability of exceeding a sediment concentration greater than 837 mg/L (a threshold of concern discussed in @gbrData) is plotted in order to draw attention to the high-risk regions on the map. These probabilities were calculated from the posterior distributions of a BHM as outlined in @gbrData.

![Exceedance probability map showing the probability of exceeding the nominated guideline of 837 mg/L for TSS in the upper Burdekin catchment in Queensland, Australia.\label{exceedMap}](exceedMap.png){width=4in height=4in}

# Acknowledgements

We would like to acknowledge the support and funding from the CSIRO Digiscape Future Science Platform for the second author.

`Vizumap` was built using the following R packages: `ggplot2`, `animation`, `broom`, `dplyr`, `geoaxe`, `ggmap`, `grDevices`, `gridExtra`, `maps`, `maptools`, `plyr`, `reshape2`, `rgdal`, `rgeos`, `roxygen2`, `sp`, `spbabel`, `testthat`, `usethis`, and `utils`.

# References
# Vizumap JOSS paper revisions

(DN) - Edit based on revision suggested by Daniel Nüst (https://github.com/openjournals/joss-reviews/issues/2409#issuecomment-664007471)

(SD) - Edit based on revision suggested by Shaoqing Dai (https://github.com/openjournals/joss-reviews/issues/2409#issuecomment-736095413)

### README 
* Edited method descriptions
* (DN) Hyperlinked to open-access papers
* (DN) Fixed license description
* (DN) Switched from devtools to remotes installation
* (DN) Removed Zenodo reference
* (DN) Added section about contributing to the package
  * Generated a Code of Conduct file with usethis::use_code_of_conduct()
  * Added contribution guidelines file
* (DN) Added section about running tests
* (DN) Added action badge
* (DN & editor @bstabler) Added section about Vizumap collaboration history
* (SD) Added links to pkgdown site and vignette

### Website
* (SD) Generated a website for Vizumap using pkgdown

### License
* (DN) Removed incorrect file from repo
* (DN) Added correct full text license file (GPLv3)

### Description
* (DN) Corrected package version (it should be v1.0.0)
  * In the JOSS submission form, it was submitted as v1.2.0 - is this going to be an issue?
  * UPDATED PACKAGE VERSION TO 1.1.0
* (DN) Deleted extra maintainer
* (DN) Added a longer description
* (DN & editor @bstabler) Updated author fields
* (DN) Added ORCIDs for all authors
* Added pkgdown and source code urls

### Vignette
* (DN) Added exceedance map example
* (DN) Improved variable/item naming
* (DN) Fixed Mexico border overlap
* (DN) Added all() wrapper
* Added California border to pixel map to improve map interpretability
* Included running of pixelation code (seems fast enough to include)

### Vizumap code
* (DN) Added `\dontrun{}` to examples instead of commenting out
* (DN) Renamed "shapefile" argument to "geoData"
* Reformatted comments to follow Hadley style guide ("# comment" instead of "#comment")
* (DN) Started adding testing for data formatting
* (DN and SD) Addressed pixelate warnings (remove and then return to original projection)
* (DN) Added comment to pixelate documentation about projection changes

### Paper
* (DN) Moved the directory to Rbuildignore
* (DN) Included an R script with code that produces the paper figures
* Added funding acknowledgement
* (DN) Incorporated Daniel’s rephrasing suggestions
  * “finding methods that add additional elements” to “finding methods that can communicate additional information, about the spatial estimates, in an understandable and meaningful way”
  * “functions can be found in the package download” to “is available after package installation”
* (DN) Clarified audience for statement of need requirement:
  * “However, in spatial applications, finding methods that can communicate additional information, about the spatial estimates, in an understandable and meaningful way can be challenging. To address this visualisation shortcoming in spatial statistics, we developed the Vizumap R package. **It is a toolkit designed for statisticians, scientists, data journalists, etc., discussing uncertainty in spatial data.**”
* (DN) Added zoom view of five regions for the pixel map figure and a sentence about this in the map description
* (DN) Commented on related R packages
* Added testthat and usethis to package acknowledgements
* (DN) Fixed text under pixel figure
* (SD) Added figure captions
* (SD) Numbered equations
* (SD) Added sub-figure labels to pixel map
* Changed "TSS loads" to "TSS concentrations"
* Edits to pixel sampling description







---
output: rmarkdown::github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# Vizumap <img src='man/figures/Vizumap_Hex.png' align="right" height="138.5" />
## An R package for visualizing uncertainty in spatial data. 

[![R build status](https://github.com/lydialucchesi/Vizumap/workflows/R-CMD-check/badge.svg)](https://github.com/lydialucchesi/Vizumap/actions)

There is a [`Vizumap` pkgdown site](https://lydialucchesi.github.io/Vizumap/) with a vignette.

A [`Vizumap` paper](https://doi.org/10.21105/joss.02409) is available in the Journal of Open Source Software (JOSS). If you use `Vizumap`, please cite this paper.

## Installation

You can install a development version of the `Vizumap` package using the command below.

```
remotes::install_github(repo = "lydialucchesi/Vizumap", build_vignettes = TRUE, force = TRUE)
```
## Authors 

Lydia Lucchesi, Australian National University & CSIRO Data61,  Email: Lydia.Lucchesi@anu.edu.au

Petra Kuhnert, CSIRO Data61, Email: Petra.Kuhnert@data61.csiro.au 

## About the package

Approaches for visualising uncertainty in spatial data are presented in this package. These include the three approaches developed in [Lucchesi and Wikle (2017)](http://faculty.missouri.edu/~wiklec/LucchesiWikle2017Stat) and a fourth approach presented in [Kuhnert et al. (2018)](https://publications.csiro.au/publications/#publication/PIcsiro:EP168206). The package is outlined in [Lucchesi et al. (2021)](https://doi.org/10.21105/joss.02409).

#### Bivariate Maps

In these bivariate choropleth maps, two colour schemes, one representing the estimates and another representing the margins of error, are blended so that an estimate and its error can be conveyed on a map using a single colour. 

#### Map Pixelation

In this approach, each map region is pixelated. Pixels are filled with colours representing values within an estimate's margin of error. Regions that appear as a solid colour reflect smaller margins of error, while more pixelated regions indicate greater uncertainty. These maps can be animated to provide a novel uncertainty visualisation experience.

#### Glyph Rotation

In this method, glyphs located at region centroids are rotated to represent uncertainty. The colour filling each glyph corresponds to the estimate.

#### Exceedance Probability Maps

The final map-based exploration is through exceedance probabilities, which are visualised on a map to highlight regions that exhibit varying levels of departure from a threshold of concern or target. 

## Examples

A vignette for the `Vizumap` package is available and contains examples relating to each of the visualisation methods.

```
vignette("Vizumap")
```

## Testing

If you would like to install and run the unit tests interactively, include `INSTALL_opts = "--install-tests"` in the installation code.

```
remotes::install_github(repo = "lydialucchesi/Vizumap", build_vignettes = TRUE, force = TRUE, INSTALL_opts = "--install-tests")

testthat::test_package("Vizumap", reporter = "stop")
```

## Contribute

To contribute to `Vizumap`, please follow these [guidelines](CONTRIBUTING.md).

Please note that the `Vizumap` project is released with a [Contributor Code of Conduct](CONDUCT.md). By contributing to this project, you agree to abide by its terms.

## License

`Vizumap` version 1.2.0 is licensed under [GPLv3](LICENSE.md).

## Citation

Lucchesi et al., (2021). Vizumap: an R package for visualising uncertainty in spatial data. Journal of Open Source Software, 6(59), 2409, https://doi.org/10.21105/joss.02409

```
@article{lucchesi2021vizumap,
  title={Vizumap: an R package for visualising uncertainty in spatial data},
  author={Lucchesi, Lydia R and Kuhnert, Petra M and Wikle, Christopher K},
  journal={Journal of Open Source Software},
  volume={6},
  number={59},
  pages={2409},
  year={2021}
}
```

## History of Vizumap

Vizumap began as a visualisation project at the University of Missouri in 2016. Chris Wikle, professor of statistics, posed an interesting research question to Lydia Lucchesi, a student curious about data visualisation and R.

How do you include uncertainty on a map displaying areal data estimates?

Over the course of a year, they put together three methods for visualising uncertainty in spatial statistics: the bivariate choropleth map, the pixel map, and the glyph map. By mid-2017, there were maps, and there was a lot of R code, but there was not a tool that others could use to easily make these types of maps, too. That’s when statistician Petra Kuhnert recommended developing an R package. Over the course of a month, Petra and Lydia developed Vizumap (originally named VizU) at CSIRO Data61 in Canberra, Australia. Since then, the package has been expanded to include exceedance probability maps, an uncertainty visualisation method developed by Petra while working on a Great Barrier Reef (GBR) project.

Vizumap has been used to visualise the uncertainty of American Community Survey estimates, the prediction errors of sediment estimates in a GBR catchment, and most recently the [uncertainty of estimated locust densities in Australia](https://www.nature.com/articles/s41598-020-73897-1/figures/4). We would like to assemble a Vizumap gallery that showcases different applications of the package’s mapping methods. If you use Vizumap to visualise uncertainty, please feel free to send the map our way. We would like to see it!

## References

Kuhnert, P.M., Pagendam, D.E., Bartley, R., Gladish, D.W., Lewis, S.E. and Bainbridge, Z.T. (2018) [Making management decisions in face of uncertainty: a case study using the Burdekin catchment in the Great Barrier Reef, Marine and Freshwater Research](https://publications.csiro.au/publications/#publication/PIcsiro:EP168206), 69, 1187-1200, https://doi.org/10.1071/MF17237.

Lucchesi, L.R. and Wikle C.K. (2017) [Visualizing uncertainty in areal data with bivariate choropleth maps, map pixelation and glyph rotation](http://faculty.missouri.edu/~wiklec/LucchesiWikle2017Stat), Stat, https://doi.org/10.1002/sta4.150.

Lucchesi, L.R., Kuhnert, P.M. and Wikle, C.K. (2021) [Vizumap: an R package for visualising uncertainty in spatial data](https://doi.org/10.21105/joss.02409), Journal of Open Source Software, https://doi.org/10.21105/joss.02409.
---
title: "Using Vizumap"
author: "Lydia Lucchesi and Petra Kuhnert"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using Vizumap}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(ggplot2)
library(dplyr)
```

**Vizumap** is a tool for visualising uncertainty in spatial data. This vignette demonstrates how to use the package to create bivariate maps, pixel maps, glyph maps, and exceedance probability maps. We have tried to generalise the visualisation approaches so that they are applicable for most types of spatial data.

There are three types of functions in **Vizumap**: formatting functions, building functions, and viewing functions. Formatting functions prepare data frames and spatial polygons for use in building functions. Building functions create maps, animations, keys, and colour palettes. Viewing functions plot these objects. Some of these functions are useful across all four map types, such as `read.uv`, and others are map specific, such as `pixelate`. Below is a table displaying the functions that can be used with each map type. Blue functions format. Green functions create. Red functions plot.

```{r echo = FALSE, fig.align='center', fig.height=5, fig.width=8, fig.asp=.6, warning=FALSE}

b <- data.frame(name = rev(c("read.uv", "build_bmap", "build_bkey", "build_palette", "view", "attach_key")),
                x = rep(0, 6), 
                y = seq(0, 5, 1), 
                col = rev(c("midnightblue", "darkgreen", "darkgreen", "darkgreen", "darkred", "darkred")))

p <- data.frame(name = rev(c("read.uv", "pixelate", "build_pmap", "animate", "view", "")), 
                x = rep(5, 6), 
                y = seq(0, 5, 1), 
                col = rev(c("midnightblue", "midnightblue", "darkgreen", "darkgreen", "darkred", "darkred")))


g <- data.frame(name = rev(c("read.uv", "build_gmap", "build_gkey", "view", "attach_key", "")), 
                x = rep(10, 6), 
                y = seq(0, 5, 1),
                col = rev(c("midnightblue", "darkgreen", "darkgreen", "darkred", "darkred", "darkred")))

e <- data.frame(name = rev(c("read.uv", "build_emap", "view", "", "", "")), 
                x = rep(15, 6), 
                y = seq(0, 5, 1),
                col = rev(c("midnightblue", "darkgreen", "darkred", "darkred", "darkred", "darkred")))

functions <- rbind(b, p, g, e)

headings <- data.frame(name = c("Bivariate map", "Pixel map", "Glyph map", "Exceedance map"), x = c(0, 5, 10, 15), y = rep(6, 4))

g <- ggplot() + geom_text(data = functions, aes(x = x, y = y, label = name, colour=col), size = 4, fontface = "bold") + 
  scale_colour_identity() + 
  geom_text(data = headings, aes(x = x, y = y, label = name), size = 5, fontface = "bold") +
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank()) + lims(x = c(-1, 17), y = c(-.5, 7))
plot(g)
```

## Data

One of the included datasets in the **Vizumap** package is `us_data`, which contains the estimated family poverty rates for all US counties in 2015. Also included is a spatial polygons data frame for a US county map (`us_geo`), which can be used to map the estimates and errors in `us_data`. These data come from the US Census Bureau [Fact Finder website](https://factfinder.census.gov/), and documentation can be found under `?us_data` and `?us_geo`.

```{r, echo=TRUE, fig.align = "center", fig.width=4, fig.height=4, message = FALSE}
library(Vizumap)

data(us_data)
head(us_data)
```

```{r, echo=TRUE, eval = TRUE, fig.align = "center", fig.width=4, fig.height=4, message = FALSE}
data(us_geo)
```


## Building bivariate maps

A bivariate map is based on a bivariate colour palette, which is created by blending two single hue colour palettes. One colour palette represents the variable of interest while the other represents the uncertainty. When the palettes are blended, each colour in the new palette becomes representative of both a variable class and uncertainty class.

There are four pre-prepared bivariate colour palettes included in **Vizumap**: `'BlueYellow'`, `'CyanMagenta'`, `'BlueRed'` and `'GreenBlue'`.

```{r, echo = TRUE, fig.align = "center", fig.width=2, fig.height=2}
# use one of four pre-prepared colour palettes
cmBivPal <- build_palette(name = "CyanMagenta")
view(cmBivPal)
```

You can also design your own bivariate colour palette with `build_palette`. Instead of entering one of the four pre-prepared palettes into `name`, you enter `'usr'`. If you choose two light colours when creating your own palette, the palette will lack important differences in colour across the grid. Therefore, it is best to use two darker colours. For example, `'chartreusu4'` and `'darkblue'` work better than `'tan2'` and `'lightskyblue'` (see comparison below). The numeric vector passed to `difC` controls the colour change across each single hue colour gradient.

```{r, echo = TRUE, fig.align = "center", fig.width=2, fig.height=2}
# test two light colours with low difC values
# creates a bad palette
customBivPal1 <- build_palette(name = "usr", colrange = list(colour = c("tan2", "lightskyblue"), difC = c(1, 1)))
view(customBivPal1)

# change colours for a better palette
customBivPal2 <- build_palette(name = "usr", colrange = list(colour = c("chartreuse4", "darkblue"), difC = c(1, 1)))
view(customBivPal2)

# change difC values to increase colour differences
customBivPal3 <- build_palette(name = "usr", colrange = list(colour = c("chartreuse4", "darkblue"), difC = c(3, 4)))
view(customBivPal3)
```

A bivariate map can be used to simultaneously visualise the estimates and errors in `us_data`. Before the data frame of estimates and errors is passed to `build_bmap`, it must be formatted with `read.uv` as the estimates and errors need to be in the first and second columns of the data frame. You can also import a dataset by passing a CSV file pathway to `file`.

```{r, echo = TRUE, fig.align = "center"}
# load data
data(us_data)
data(us_geo)

# current column order of us_data
head(us_data)

# format data frame
poverty <- read.uv(data = us_data, estimate = "pov_rate", error = "pov_moe")

# estimates and errors are now in the first two columns
head(poverty)
```

In this example, terciles are used to define the numerical bounds for the different estimate and error classes because the dataset contains several outliers.
```{r, echo = TRUE, fig.align = "center", fig.width=8, fig.height=6, fig.asp=.7}
# build a bivariate map with the map data
usBivMap <- build_bmap(data = poverty, geoData = us_geo, id = "GEO_ID", terciles = TRUE)
view(usBivMap)
```

Keys for the bivariate maps are not automatically generated with `build_bmap` and must be created separately with `build_bkey`. It is important that the key arguments match the map arguments. For example, if `terciles` was set to `FALSE`, the key would not accurately reflect the county colour assignments for this example US map.

```{r, echo = TRUE, fig.align = "center", fig.width=4, fig.height=4}
# build a key
usBivKey <- build_bkey(data = poverty, terciles = TRUE)
view(usBivKey)
```

Maps and keys can be viewed together with `attach_key`.
```{r, echo=TRUE, fig.align="center", fig.asp=.5, fig.height=10, fig.width=14, message=FALSE, warning=FALSE}
attach_key(usBivMap, usBivKey)
```

Changing the `terciles` and `palette` arguments leads to a map that looks very different than the previous map even though it is displaying the same dataset. Here the custom colour palette is used.
```{r, echo = TRUE, eval = TRUE, fig.align = "center", fig.width=8, fig.height=8, fig.asp=.75}
# make some changes
usBivMapDif <- build_bmap(data = poverty, geoData = us_geo, id = "GEO_ID", terciles = FALSE, palette = customBivPal3)
view(usBivMapDif)
```

## Building pixel maps

Pixel maps are created by pixelating regions and assigning each pixel in a region a value from an estimate's confidence interval or discrete relative frequency distribution. The first step is to format the data frame containing the data with `read.uv` and pixelate the map with `pixelate`.

A California county map will be used to illustrate the method. A subset of `us_data` that contains only California estimates and errors is created.
```{r, echo = TRUE}
data(us_data)
us_data$GEO.id2 <- as.numeric(us_data$GEO.id2)
ca_data <- subset(us_data, us_data$GEO.id2 > 6000 & us_data$GEO.id2 < 7000)
ca_data <- read.uv(data = ca_data, estimate = "pov_rate", error = "pov_moe")
row.names(ca_data) <- seq(1, nrow(ca_data), 1)
```

The county polygons for California are extracted from `us_geo` and pixelated with `pixelate`. It can take several minutes to pixelate a shapefile depending on the size of the shapefile.
```{r, echo = TRUE, eval = TRUE}
data(us_geo)
ca_geo <- subset(us_geo, us_geo@data$STATE == "06")
pix <- pixelate(ca_geo, id = "region")
```

There must be a shared column between `ca_data` and `pix`. Depending on your shapefile, getting this shared column might require some creativity. `pixelate` returns a column that contains the slot IDs for the polygons. Adding these slot IDs to `ca_data` through the `GEO_ID` column present in both `ca_data` and `ca_geo` creates a shared column that can be used for the `id` argument in `build_pmap`. Ultimately, if the shapefile slot IDs are not referenced in your data frame of estimates and errors, you will need to find a way to way to get this information added to your data frame of estimates and errors before using `build_pmap`.
```{r, echo = TRUE, eval = TRUE}
df <- data.frame(region = sapply(slot(ca_geo, "polygons"), function(x) slot(x, "ID")), name = unique(ca_geo@data$GEO_ID))
ca_data$region <- df[match(ca_data$GEO_ID, df$name), 1]
ca_data$region <- as.character(ca_data$region)

# check that values in shared column match
all(ca_data$region %in% pix$region)
```

If `distribution = "uniform"`, the values are sampled uniformly from an interval, where the lower bound is the estimate minus the error and the upper bound is the estimate plus the error. For the estimates and margins of error in `us_data`, this interval corresponds to the estimate's 90% confidence interval.
```{r, echo = TRUE, fig.height=8, fig.width=8, eval = TRUE, message = FALSE}
# uniform distribution
unifPixMap <- build_pmap(data = ca_data, distribution = "uniform", pixelGeo = pix, id = "region", border = ca_geo)
view(unifPixMap)
```

If `distribution = "normal"`, the values assigned to pixels will be drawn from normal distributions parameterised using the estimates and errors (means and standard deviations).
```{r, echo = TRUE, fig.height=8, fig.width=8, eval = TRUE, message = FALSE}
# normal distribution
ca_data$se <- ca_data$pov_moe / 1.645
ca_data <- read.uv(data = ca_data, estimate = "pov_rate", error = "se")

normPixMap <- build_pmap(data = ca_data, distribution = "normal", pixelGeo = pix, id = "region", border = ca_geo)
view(normPixMap)
```

If `distribution = "discrete"`, a data frame of quantiles, which define the relative frequency distributions for the estimates, must be passed to `q`. Below is what the California map would look like if the relative frequency distributions for the estimates followed exponential distributions.
```{r, echo = TRUE, fig.height=8, fig.width=8, eval = TRUE, message = FALSE}
# experiment with discrete distribution
# exponential - example for q argument
ca_data.q <- with(ca_data, data.frame(p0.05 = qexp(0.05, 1/pov_rate), p0.25 = qexp(0.25, 1/pov_rate), p0.5 = qexp(0.5, 1/pov_rate), p0.75 = qexp(0.75, 1/pov_rate), p0.95 = qexp(0.95, 1/pov_rate)))

head(ca_data.q)

discPixMap <- build_pmap(data = ca_data, distribution = "discrete", 
                  pixelGeo = pix, id = "region", q = ca_data.q, border = ca_geo)
view(discPixMap)
```

Pixel maps can be animated with `animate` so that the pixels flicker between a series of assigned values. `view` saves the animation to your computer as an html file and automatically opens a browser to view it. A longer `aniLength` corresponds to a longer animation as well as a longer running time. Generating the animation with `view` can take several minutes.
```{r, echo = TRUE, eval = TRUE}
# animate the normal distribution map
normPixAni <- animate(normPixMap, aniLength = 30)
# view(normPixAni)
```

## Building glyph maps

The process of creating a glyph map is very similar to the process of creating a bivariate map. The data are formatted; the map is created; the key is created; and the two are merged. Glyphs are plotted at either region centroids or specific sites, and the colour and rotation of the glyph represent the variable of interest and the uncertainty. The method is illustrated below with a Colorado county map.

The datasets are loaded and subsetted. `co_data` is formatted with `read.uv` for use in `build_gmap` and `buid_gkey`.
```{r, echo = TRUE}
data(us_data)
data(us_geo)

co_geo <- subset(us_geo, us_geo@data$STATE == "08")

us_data$GEO.id2 <- as.numeric(us_data$GEO.id2)
co_data <- subset(us_data, us_data$GEO.id2 > 8000 & us_data$GEO.id2 < 9000)
co_data <- read.uv(data = co_data, estimate = "pov_rate", error = "pov_moe")
```

Because `geoData` is included in this example, the `build_gmap` function will plot a glyph at each region centroid. The colour of the glyph represents the estimated poverty rate among families, and the rotation of the glyph represents the margin of error for the estimate.
```{r, echo = TRUE, fig.align = "center", fig.width=8, fig.height=8, fig.asp=.75}
# build a glyph map
usGlyphMap <- build_gmap(data = co_data, geoData = co_geo, id = "GEO_ID", size = 80, glyph = "icone", border = "state")
view(usGlyphMap)
```

Keys for glyph maps are not automatically generated with `build_gmap` and must be created separately with `build_gkey`. It is important that the key arguments match the map arguments.
```{r, echo = TRUE, fig.align = "center", fig.width=4, fig.height=4}
# build a glyph key
usGlyphKey <- build_gkey(data = co_data, glyph = "icone")
view(usGlyphKey)
```

Maps and keys can be viewed together with `attach_key`.
```{r, echo = TRUE, fig.align = "center", fig.width=10, fig.height=6, fig.asp=.5, message = FALSE}
attach_key(usGlyphMap, usGlyphKey)
```

You can change the size, shape, and colour of the glyphs as well as add different borders.
```{r, echo = TRUE, fig.align = "center", fig.width=8, fig.height=8, fig.asp=.75}
# build a glyph map
usGlyphMapDif <- build_gmap(data = co_data, geoData = co_geo, id = "GEO_ID", size = 70, border = "county", glyph = "semi", palette = "Reds")
view(usGlyphMapDif)
```

## Building exceedance probability maps

Exceedance maps plot the probability of exceeding some nominated threshold of concern. The `UB` dataset included in the package contains pre-calculated exceedance probabilities. However, the `us_data` data frame does not contain pre-calculated exceedance probabilities, and these will therefore need to be generated. First, a threshold needs to be selected - let's choose a poverty rate of 30%.

```{r, echo = TRUE}
# load data
data(us_data)
data(us_geo)

# format the data
poverty <- read.uv(data = us_data, estimate = "pov_rate", error = "pov_moe")

# check variable quantiles
quantile(us_data$pov_rate)
```

An exponential distribution is an appropriate choice for this application (a normal distribution would not be). A list containing information about the distribution and threshold is passed to the `build_emap` function, where the probabilities are calculated and plotted. Below is an example of how to prepare this list.

```{r, echo = TRUE}
# define probability distribution (exponential distribution)
pd <- quote({ pexp(q, rate, lower.tail = FALSE) })

# define argument listing
args <- quote({ list(rate = 1/estimate) })

# capture distribution and arguments in a single list
pdflist <- list(dist = pd, args = args, th = 30)
```

Finally, we build the exceedance map using the `build_emap` function. We need to supply the formatted data frame, instructions for calculating the exceedance probabilities, and information about map design.
```{r echo=TRUE, fig.align="center", fig.asp=.6, fig.height=8, fig.width=8, message=FALSE, warning=FALSE}
usExcMap <- build_emap(data = poverty, pdflist = pdflist, geoData = us_geo, id = "GEO_ID", key_label = "Pr[X > 30]")
view(usExcMap)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_pmap.R
\name{build_pmap}
\alias{build_pmap}
\title{Build a pixel map}
\usage{
build_pmap(
  data = NULL,
  distribution = NULL,
  pixelGeo,
  id,
  border = NULL,
  palette = "Blues",
  q = NULL,
  limits = NULL
)
}
\arguments{
\item{data}{A data frame.}

\item{distribution}{Name of the distribution that the pixel assignments will
be drawn from. It must be one of \code{discrete}, \code{normal} or
\code{uniform}. If \code{distribution = "discrete"}, a data frame of the
quantiles that define the relative frequency distribution for the
estimate must be entered for \code{q}. If \code{distribution = "normal"},
the values assigned to pixels will be drawn from a normal distribution
parameterised using the estimates and errors (means and standard
deviations). If \code{distribution = "uniform"}, values will be sampled with
equal probability from a sequence of 5 numbers that spans the estimate minus
its error to the estimate plus its error.}

\item{pixelGeo}{An object from \code{\link{pixelate}}.}

\item{id}{Name of the common column shared by the objects passed to
\code{data}, \code{pixelGeo} and \code{q} (if \code{distribution =
"discrete"}).}

\item{border}{Name of geographical borders to be added to the map. It must be
one of \code{\link[maps]{county}}, \code{\link[maps]{france}},
\code{\link[maps]{italy}}, \code{\link[maps]{nz}},
\code{\link[maps]{state}}, \code{\link[maps]{usa}} or
\code{\link[maps]{world}} (see documentation for
\code{\link[ggplot2]{map_data}} for more information). The borders will be
refined to match latitute and longtidue coordinates provided in the data
frame or spatial polygons data frame. Alternatively, you can supply a \code{SpatialPolygonsDataFrame}.}

\item{palette}{Name of colour palette. It must be one of \code{Blues},
\code{Greens}, \code{Greys}, \code{Oranges}, \code{Purples} or \code{Reds}
(see documentation for \code{\link[ggplot2]{scale_fill_distiller}} for more
information).}

\item{q}{A data frame of quantiles which define the distribution for each
estimate. Each row is an estimate, and each column is a quantile. See
examples for an example of \code{q} input.}

\item{limits}{Limits for the legend. Default is NULL, which takes the limits to be the range of the data.}
}
\description{
This function builds a choropleth map that visualises estimates and errors
simultaneously through pixelation and sampling.
}
\examples{
\dontrun{
# This code will produce a pixelated map when run in R
# It is not run here.
data(us_geo)
ca_geo <- subset(us_geo, us_geo@data$STATE == "06")
pix <- pixelate(ca_geo, id = "region")

data(us_data)
us_data$GEO.id2 <- as.numeric(us_data$GEO.id2)
ca_data <- subset(us_data, us_data$GEO.id2 > 6000 & us_data$GEO.id2 < 7000)
ca_data <- read.uv(data = ca_data, estimate = "pov_rate", error = "pov_moe")
row.names(ca_data) <- seq(1, nrow(ca_data), 1)

df <- data.frame(region = sapply(slot(ca_geo, "polygons"),
 function(x) slot(x, "ID")), name = unique(ca_geo@data$GEO_ID))
ca_data$region <- df[match(ca_data$GEO_ID, df$name), 1]
ca_data$region <- as.character(ca_data$region)

# uniform distribution
m <- build_pmap(data = ca_data, distribution = "uniform", pixelGeo = pix, id = "region")
view(m)

# normal distribution
ca_data$se <- ca_data$pov_moe / 1.645
ca_data <- read.uv(data = ca_data, estimate = "pov_rate", error = "se")

m <- build_pmap(data = ca_data, distribution = "normal", pixelGeo = pix, id = "region")
view(m)

# experiment with discrete distribution
# exponential - example for q argument
ca_data.q <- with(ca_data, data.frame(p0.05 = qexp(0.05, 1/pov_rate),
 p0.25 = qexp(0.25, 1/pov_rate), p0.5 = qexp(0.5, 1/pov_rate),
 p0.75 = qexp(0.75, 1/pov_rate), p0.95 = qexp(0.95, 1/pov_rate)))

m <- build_pmap(data = ca_data, distribution = "discrete", pixelGeo = pix,
 id = "region", q = ca_data.q)
view(m)
}
}
\seealso{
\code{\link{animate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/us_geo-data.r
\docType{data}
\name{us_geo}
\alias{us_geo}
\title{USA shapefile}
\format{
A large spatial polygons data frame with 3142 elements and 7 variables:
 \describe{
 \item{GEO_ID}{14-digit code that identifies summary level of data and region}
 \item{STATE}{state FIPS code}
 \item{COUNTY}{county FIPS code}
 \item{NAME}{region name}
 \item{LSAD}{legal statistical area description}
 \item{SHAPE_AREA}{area measurement in square meters}
 \item{SHAPE_LEN}{perimeter measurement in meters}
 }
}
\source{
\url{http://factfinder.census.gov/}
}
\usage{
us_geo
}
\description{
A dataset that produces a USA county map. It can be used to
visualise data in the us_data data frame spatially.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UB_tss.r
\docType{data}
\name{UB_tss}
\alias{UB_tss}
\title{UB_tss}
\format{
A data frame comprising 411 rows (one for each RCA) and 5 columns.
}
\usage{
UB_tss
}
\description{
A dataset that contains estimates of Total Suspended Sediment (TSS),
the error, i.e. standard deviation (TSS_error), the exceedance probabilities that
correspond to thresholds of 837 mg/L and 2204 mg/L, respectively for each polygon
or reach contributing area (RCA). These thresholds were determined from Bartley et al.
(2012). The TSS estimates, error and exceedance probabilities were derived from the
model developed and described in Gladish et al. (2016).
}
\references{
Bartley, R., Speirs, W. J., Ellis, T. W., and Waters, D. K. (2012).
A review of sediment and nutrient concentration data from Australia for use in
catchment water quality models. Marine Pollution Bulletin 65, 101–116.
doi:10.1016/J.MARPOLBUL.2011.08.009
Gladish, D. W., Kuhnert, P. M., Pagendam, D. E., Wikle, C. K., Bartley, R.,
Searle, R. D., Ellis, R. J., Dougall, C., Turner, R. D. R., Lewis, S. E.,
Bainbridge, Z. T., and Brodie, J. E. (2016). Spatio-temporal assimilation
of modelled catchment loads with monitoring data in the Great Barrier
Reef. The Annals of Applied Statistics 10, 1590–1618. doi:10.1214/16-
AOAS950
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_gkey.R
\name{build_gkey}
\alias{build_gkey}
\title{Build a glyph key}
\usage{
build_gkey(data, glyph = "icone")
}
\arguments{
\item{data}{A data frame.}

\item{glyph}{Name of glyph shape. Options include \code{icone} and \code{semi}.}
}
\description{
This function creates a key of rotated glyphs for a map produced with \code{\link{build_gmap}}.
}
\details{
A key for the glyph map is not automatically generated with
\code{\link{build_gmap}} and must be made using \code{\link{build_gkey}}. It is important
that the arguments passed to \code{\link{build_gkey}} match those passed to
\code{\link{build_gmap}}. The map and key can be viewed together using
\code{\link{attach_key}}.
}
\examples{
data(us_data)
data(us_geo)
co_geo <- subset(us_geo, us_geo@data$STATE == "08")
us_data$GEO.id2 <- as.numeric(us_data$GEO.id2)
co_data <- subset(us_data, us_data$GEO.id2 > 8000 & us_data$GEO.id2 < 9000)
co_data <- read.uv(data = co_data, estimate = "pov_rate", error = "pov_moe")

# build a glyph key
key <- build_gkey(data = co_data, glyph = "icone")
view(key)

}
\seealso{
\code{\link{attach_key}}
}
\name{Vizumap-internal}
\alias{calcbin.r}
\alias{findNbounds.r}
\alias{view.bivmap.r}
\alias{view.key.r}
\alias{view.palette.r}



\title{Internal Vizumap functions}
\description{
These are internal \code{Vizumap} functions that are not intended to be directly called by the user
}


\author{Petra Kuhnert and Lydia Lucchesi}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixelate.R
\name{pixelate}
\alias{pixelate}
\title{Pixelate polygons}
\usage{
pixelate(geoData = NULL, file = NULL, layer = NULL, pixelSize = 2, id = "id")
}
\arguments{
\item{geoData}{An object of class "SpatialPolygons" or
"SpatialPolygonsDataFrame".}

\item{file}{A shapefile pathway.}

\item{layer}{Name of geoData layer (see documentation for
\code{\link[rgdal]{readOGR}} for more information).}

\item{pixelSize}{An integer 1, 2 or 3. One corresponds to the smallest pixel
size, and three corresponds to the largest.}

\item{id}{A name which will be given to the new ID column. This ID corresponds to the slot ID in
the spatial data.}
}
\description{
This function divides a shapefile into pixels so it can be used to create a
pixel map with \code{\link{build_pmap}}.
}
\details{
This function can take several minutes to run depending on the size of the shapefile. Within this function, the projection of the spatial object is removed and then returned to the original projection.
}
\examples{
data(us_geo)
ca_geo <- subset(us_geo, us_geo@data$STATE == "06")
pix <- pixelate(ca_geo, id = "region")


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_bmap.R
\name{build_bmap}
\alias{build_bmap}
\title{Build a bivariate colour map}
\usage{
build_bmap(
  data,
  geoData = NULL,
  id = NULL,
  border = NULL,
  palette = "BlueYellow",
  size = NULL,
  terciles = FALSE,
  bound = NULL,
  flipAxis = FALSE
)
}
\arguments{
\item{data}{A data frame.}

\item{geoData}{A spatial polygons data frame.}

\item{id}{Name of the common column shared by the objects passed to
\code{data} and \code{geoData}. The estimates and errors in the data frame
will be matched to the geographical regions of the spatial polygons data
frame through this column.}

\item{border}{Name of geographical borders to be added to the map. It must be
one of \code{\link[maps]{county}}, \code{\link[maps]{france}},
\code{\link[maps]{italy}}, \code{\link[maps]{nz}},
\code{\link[maps]{state}}, \code{\link[maps]{usa}} or
\code{\link[maps]{world}} (see documentation for
\code{\link[ggplot2]{map_data}} for more information). The borders will be
refined to match latitute and longtidue coordinates provided in the data
frame or spatial polygons data frame.}

\item{palette}{Name of colour palette or character vector of hex colour codes
from the \code{\link{build_palette}} function. Colour palette names include
\code{BlueYellow}, \code{CyanMagenta}, \code{BlueRed} and \code{GreenBlue}.}

\item{size}{An integer between 1 and 20. Value controls the size of points
when \code{geoData = NULL}. If \code{size = NULL}, the points will remain
the default size.}

\item{terciles}{A logical value. This provides the option to define numerical
bounds for the colour key grid using terciles instead of equal intervals.}

\item{bound}{Output from the \code{findNbounds} function if a different
set of data is required to bound the map.  This is useful if you are wanting
to create a bivariate map across multiple years and show colours that correspond
to the same key. Default is NULL.}

\item{flipAxis}{A logical value. Whether to place the axis on the opposite sides or not.}
}
\description{
This function builds a map that visualises estimates and errors simultaneously
with a bivariate colour scheme.
}
\details{
If \code{geoData} remains \code{NULL}, the function will produce a map of
plotted points representing specific sites; in this case, the data frame must
include latitude and longitude coordinates in columns \code{"long"} and
\code{"lat"}.
}
\examples{
data(us_data)
data(us_geo)
poverty <- read.uv(data = us_data, estimate = "pov_rate", error = "pov_moe")

# bivariate map with a spatial polygons data frame
map <- build_bmap(data = poverty, geoData = us_geo, id = "GEO_ID",
 border = "state", terciles = TRUE)
view(map)

}
\seealso{
\code{\link{build_bkey}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attach_key.R
\name{attach_key}
\alias{attach_key}
\title{Attach a key}
\usage{
attach_key(map, mapkey)
}
\arguments{
\item{map}{An object generated with \code{\link{build_bmap}} or \code{\link{build_gmap}}.}

\item{mapkey}{An object generated with \code{\link{build_bkey}} or \code{\link{build_gkey}}.}
}
\description{
This function attaches a bivariate colour key to a
bivariate colour map or a glyph key to a glyph map so that they can be viewed together.
}
\examples{
data(us_data)
data(us_geo)
poverty <- read.uv(data = us_data, estimate = "pov_rate", error = "pov_moe")

# bivariate map and key together
map <- build_bmap(data = poverty, geoData = us_geo, id = "GEO_ID",
 border = "state", terciles = TRUE)
key <- build_bkey(data = poverty, terciles = TRUE)
attach_key(map, key)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdist.r
\name{pdist}
\alias{pdist}
\title{pdist}
\usage{
pdist(pname, th, args)
}
\arguments{
\item{pname}{Distribution to be used which is embedded through the \code{quote} command.
This can comprise an existing function from the \code{stats} package e.g. \code{pnorm}.
Alternatively this can be user defined.}

\item{th}{Threshold for calculating the \[Pr(X > th)\]}

\item{args}{Arguments that correspond to the specified distribution (but excluding the
threshold.}
}
\description{
Calculates the probability of exceedance for a given distribution,
threshold and distribution parameters.
}
\examples{
pname <- quote({ pnorm(q, mean, sd, lower.tail = FALSE) })
pdist(pname = pname, th = 3, args = list(mean = 2, sd = 1))


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_emap.r
\name{build_emap}
\alias{build_emap}
\title{Build an exceedance probability colour map}
\usage{
build_emap(
  data,
  pdflist = NULL,
  geoData = NULL,
  id = NULL,
  key_label,
  palette = "YlOrRd",
  size = NULL,
  border = NULL
)
}
\arguments{
\item{data}{A data frame containing columns that house estimates of the mean,
standard deviation (or margin of error) and exceedance probability (optional).
A number of options are considered for calculating the probability of exceeding
a threshold. See below for more information.}

\item{pdflist}{A list capturing the pdf function that defines the distribution function to use to
calculate the probabilities of exeedence. By default this is NULL and assumes
the exceedance probabilities have been calculated outside of the function and passed
as a third column of the data dataframe.  Functions need to conform to the class
of distribution functions available within R through the \code{stats} package.}

\item{geoData}{A spatial polygons data frame.}

\item{id}{Name of the common column shared by the objects passed to
\code{data} and \code{geoData}. The exceedance probability in the data frame
will be matched to the geographical regions of the spatial polygons data
frame through this column.}

\item{key_label}{Label of legend.}

\item{palette}{Name of colour palette. Colour palette names include
\code{BlueYellow}, \code{CyanMagenta}, \code{BlueRed}, \code{GreenBlue} and \code{YellowRed}.}

\item{size}{An integer between 1 and 20. Value controls the size of points
 when \code{geoData = NULL}. If \code{size = NULL}, the points will remain
 the default size.

 @details An exceedance probability map can be produced using:
 (i) precalculated exceedance probabilities, which are provided as a third column
 to the input dataframe; or
 (ii) exceedance probabilities that are calculated within the function using one of
 the standard probability distributions (e.g. \code{pnorm}) provided in the \code{stats}
 package in R, or
 (iii) exceedance probabilities that are calculated through a user defined function that
 is passed to the package which conforms to a similar structure as the suite of
 \code{Distributions} available in R.  Examples are provided below.}

\item{border}{Name of geographical borders to be added to the map. It must be
one of \code{\link[maps]{county}}, \code{\link[maps]{france}},
\code{\link[maps]{italy}}, \code{\link[maps]{nz}},
\code{\link[maps]{state}}, \code{\link[maps]{usa}} or
\code{\link[maps]{world}} (see documentation for
\code{\link[ggplot2]{map_data}} for more information). The borders will be
refined to match latitute and longtidue coordinates provided in the data
frame or spatial polygons data frame.}
}
\description{
This function builds a map that visualises the probability of exceeding some
nominated threshold of concern.
}
\details{
If \code{geoData} remains \code{NULL}, the function will produce a map of
plotted points representing specific sites; in this case, the data frame must
include latitude and longitude coordinates in columns \code{"long"} and
\code{"lat"}.
}
\examples{
data(us_data)
data(us_geo)
poverty <- read.uv(data = us_data, estimate = "pov_rate", error = "pov_moe")

# Exceedance probability map: Pr[X > 30] (Exponential Distribution)
#---- define probability distribution
pd <- quote({ pexp(q, rate, lower.tail = FALSE) })
#---- define argument listing
args <- quote({ list(rate = 1/estimate) })
#---- capture distribution and arguments in a single list
pdflist <- list(dist = pd, args = args, th = 30)
map <- build_emap(data = poverty, pdflist = pdflist, geoData = us_geo, id = "GEO_ID",
            border = "state", key_label = "Pr[X > 30]")
view(map) + ggplot2::ggtitle("Proper use of build_emap (appropriate distribution choice)")

# Example where an inappropriate distributions is tried
# Exceedance probability map: Pr[X>30] (Normal Distribution)

#---- define probability distribution
pd <- quote({ pnorm(q, mean, sd, lower.tail = FALSE) })
#---- define argument listing
args <- quote({ list(mean = estimate, sd = error/1.645) })
#---- capture distribution and arguments in a single list
pdflist <- list(dist = pd, args = args, th = 30)
map <- build_emap(data = poverty, pdflist = pdflist, geoData = us_geo, id = "GEO_ID",
            border = "state", key_label = "Pr[X > 30]")
view(map) + ggplot2::ggtitle("Misuse of build_emap (inappropriate distribution choice)")

# Example where exceedance probabilities have been supplied (GBR Example)
# Load Upper Burdekin Data
data(UB)

# Build Palette
exc_pal <- build_palette(name = "usr", colrange = list(colour = c("yellow", "red"),
                                                       difC = c(1, 1)))
                                                       view(exc_pal)
# Create map and view it
map <- build_emap(data = UB_tss,  geoData = UB_shp, id = "scID",
            key_label = "Pr[TSS > 837mg/L]")
view(map)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/us_data-data.r
\docType{data}
\name{us_data}
\alias{us_data}
\title{USA Poverty Data}
\format{
A data frame with 3142 elements and 5 variables:
 \describe{
 \item{GEO_ID}{14-digit code that identifies summary level of data and region}
 \item{GEO.id2}{FIPS code}
 \item{GEO.display.label}{geographic name for geographic area}
 \item{pov_rate}{percentage of families whose income was below the poverty level}
 \item{pov_moe}{margin of error for poverty estimate}
 }
}
\source{
\url{http://factfinder.census.gov/}
}
\usage{
us_data
}
\description{
A dataset containing estimated poverty rates among American families for
2015. Included is a margin of error for each estimate,
which can be used to calculate the estimate's 90\% confidence interval.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.uv.r
\name{read.uv}
\alias{read.uv}
\title{Format data}
\usage{
read.uv(data = NULL, file = NULL, estimate, error, exceedance)
}
\arguments{
\item{data}{A data frame.}

\item{file}{A CSV file pathway.}

\item{estimate}{Name of estimate column.}

\item{error}{Name of error column.}

\item{exceedance}{Name of exceedance probability column.}
}
\description{
This function formats a data frame so it can be used with
\code{\link{build_bmap}}, \code{\link{build_gmap}}, \code{\link{build_pmap}},
\code{\link{build_bkey}} and \code{\link{build_gkey}}.
}
\details{
Estimates and errors must be in the first and second columns of a data frame
for \code{\link{build_bmap}}, \code{\link{build_gmap}},
\code{\link{build_pmap}}, \code{\link{build_bkey}} and
\code{\link{build_gkey}}. \code{\link{read.uv}} provides an automated way to
format a data frame already loaded in R or to import a CSV file as a data
frame and arrange the columns correctly.
}
\examples{
data(us_data)
poverty <- read.uv(data = us_data, estimate = "pov_rate", error = "pov_moe")
head(poverty)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/view.r
\name{view}
\alias{view}
\title{View a \strong{Vizumap} object}
\usage{
view(obj)
}
\arguments{
\item{obj}{An object from \code{\link{build_bmap}}, \code{\link{build_gmap}},
\code{\link{build_pmap}}, \code{\link{build_bkey}},
\code{\link{build_gkey}}, \code{\link{build_palette}}, or
\code{\link{animate}}.}
}
\description{
This function views maps, keys, palettes and animations created with
\strong{Vizumap}.
}
\details{
For maps, keys and palettes, assigning the output from \code{\link{view}} to a variable
name saves the output as an object of class \code{"ggplot"}. Animations
created with \code{\link{animate}} are viewed through an internet browser
as html objects.
}
\examples{
gb <- build_palette(name = "GreenBlue")
view(gb)

# ggplot object
p <- view(gb)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_palette.r
\name{build_palette}
\alias{build_palette}
\title{Build a palette}
\usage{
build_palette(
  name,
  colrange = list(colour = NULL, difC = NULL),
  flipVertical = FALSE,
  flipHorizontal = FALSE
)
}
\arguments{
\item{name}{Name of colour palette or \code{usr} for option to design a new
palette. Colour palette names include
\code{BlueYellow}, \code{CyanMagenta}, \code{BlueRed} and
\code{GreenBlue}.}

\item{colrange}{List with a character vector of length two and a numeric
vector of length two.}

\item{flipVertical}{Whether the palette should be flipped vertically (ie. replace top portion with bottom portion)}

\item{flipHorizontal}{Whether the palette should be flipped horizontally (ie. replace left portion with right portion)}

\item{colour}{A character vector of two colour names from the colours() range or valid hexadecimal colors.}

\item{difC}{A numeric vector of two integers 1, 2, 3 or 4. Values
control how much a colour changes in value across the grid. One corresponds
with a small change in colour value, and four corresponds with a large
change in colour value.}
}
\description{
This function prepares one of the four included colour palettes
or builds a new colour palette.
}
\details{
Note that \code{colrange} only needs to be specified if \code{name = "usr."}
When choosing colours, it is best to avoid light colours or tints as these
will lead to a colour palette lacking noticeable differences across the 3 x 3
colour grid.
}
\examples{
# use one of four prepared colour palettes
p <- build_palette(name = "CyanMagenta")
view(p)

# design a new palette
p <- build_palette(name = "usr", colrange =
 list(colour = c("darkblue", "chartreuse4"), difC = c(3, 4)))
view(p)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findNbounds.r
\name{findNbounds}
\alias{findNbounds}
\title{findNbounds}
\usage{
findNbounds(
  data,
  estimate,
  error,
  terciles,
  expertR_est = NA,
  expertR_err = NA
)
}
\arguments{
\item{data}{Dataset from which bounds are to be determined.}

\item{estimate}{Name of the estimate column.}

\item{error}{Name of the error column.}

\item{terciles}{Should terciles be calculated? (Default: FALSE).}

\item{expertR_est}{A vector consisting of the range of expert values for the
estimate (Default: NA).}

\item{expertR_err}{A vector consisting of the range of expert values for the
error (Default: NA).}
}
\description{
Find the bins for the colour grid
}
\details{
Expert ranges for the estimate and error can be supplied to assist
  with the comparison of multiple figures. By default these are set to NA,
  which will result in the bounds being directly computed from the data.
  Note, if the expert values do not span the range of the data, this function
  will default to finding the bounds of the data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_gmap.R
\name{build_gmap}
\alias{build_gmap}
\title{Build a glyph map}
\usage{
build_gmap(
  data,
  geoData = NULL,
  id = NULL,
  size = 50,
  border = NULL,
  glyph = "icone",
  palette = "Blues",
  limits = NULL,
  max_error = NULL
)
}
\arguments{
\item{data}{A data frame.}

\item{geoData}{A spatial polygons data frame.}

\item{id}{Name of the common column shared by the objects passed to
\code{data} and \code{geoData}. The estimates and errors in the data frame
will be matched to the geographical regions of the spatial polygons data
frame through this column.}

\item{size}{An integer between 1 and 100. Value controls the size of the glyphs.}

\item{border}{Name of geographical borders to be added to the map. It must be
one of \code{\link[maps]{county}}, \code{\link[maps]{france}},
\code{\link[maps]{italy}}, \code{\link[maps]{nz}},
\code{\link[maps]{state}}, \code{\link[maps]{usa}} or
\code{\link[maps]{world}} (see documentation for
\code{\link[ggplot2]{map_data}} for more information). The borders will be
refined to match latitute and longtidue coordinates provided in the data
frame or spatial polygons data frame.}

\item{glyph}{Name of glyph shape. Options include \code{icone} and \code{semi}.}

\item{palette}{Name of colour palette. It must be one of \code{Blues},
\code{Greens}, \code{Greys}, \code{Oranges}, \code{Purples} or \code{Reds}
(see documentation for \code{\link[ggplot2]{scale_fill_distiller}} for more
information).}

\item{limits}{Limits for the legend. Default is NULL, which takes the limits
to be the range of the data.}

\item{max_error}{maximum error value. Default is NULL, which takes the
maximum from the error given}
}
\description{
This function builds a map that visualises estimates and errors simultaneously
with rotated glyphs.
}
\details{
If a spatial polygons data frame is used, glyphs will be plotted at region
centroids. If \code{geoData} remains \code{NULL}, glyphs will be plotted at
points on the map representing specific sites; in this case, the data frame
must include latitude and longitude coordinates in columns \code{"long"} and
\code{"lat"}.
}
\examples{
data(us_data)
data(us_geo)
co_geo <- subset(us_geo, us_geo@data$STATE == "08")
us_data$GEO.id2 <- as.numeric(us_data$GEO.id2)
co_data <- subset(us_data, us_data$GEO.id2 > 8000 & us_data$GEO.id2 < 9000)
co_data <- read.uv(data = co_data, estimate = "pov_rate", error = "pov_moe")

# build a glyph map
map <- build_gmap(data = co_data, geoData = co_geo, id = "GEO_ID",
 size = 70, border = "state", glyph = "icone")
view(map)

}
\seealso{
\code{\link{build_gkey}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/animate.R
\name{animate}
\alias{animate}
\title{Animate a pixel map}
\usage{
animate(pmap, flickerSpeed = 0.05, aniLength)
}
\arguments{
\item{pmap}{An object from \code{\link{build_pmap}}.}

\item{flickerSpeed}{A numeric value between 0 and 0.5.}

\item{aniLength}{An integer. A larger value corresponds with a longer
animation; a value of at least 20 is recommended.}
}
\description{
This function animates a pixel map created with \code{\link{build_pmap}} so
that the pixels ficker between a series of sampled values.
}
\examples{
\dontrun{
# This code demonstrates how to animate a pixel map
# It is not run here.
data(us_geo)
ca_geo <- subset(us_geo, us_geo@data$STATE == "06")
pix <- pixelate(ca_geo, id = "region")

data(us_data)
us_data$GEO.id2 <- as.numeric(us_data$GEO.id2)
ca_data <- subset(us_data, us_data$GEO.id2 > 6000 & us_data$GEO.id2 < 7000)
ca_data <- read.uv(data = ca_data, estimate = "pov_rate", error = "pov_moe")
row.names(ca_data) <- seq(1, nrow(ca_data), 1)

df <- data.frame(region = sapply(slot(ca_geo, "polygons"),
 function(x) slot(x, "ID")), name = unique(ca_geo@data$GEO_ID))
ca_data$region <- df[match(ca_data$GEO_ID, df$name), 1]
ca_data$region <- as.character(ca_data$region)

m <- build_pmap(data = ca_data, distribution = "uniform", pixelGeo = pix, id = "region")

a <- animate(m, aniLength = 25)
view(a)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UB_shp.r
\docType{data}
\name{UB_shp}
\alias{UB_shp}
\title{UB_shp}
\format{
A shape file stored as a \code{SpatialPolygonsDataFrame}.
}
\usage{
UB_shp
}
\description{
Upper Burdekin Shape file consisting of 411 sub-catchment
boundaries where estimates of Total Suspended Sediment (TSS) has been
quantified.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_bkey.R
\name{build_bkey}
\alias{build_bkey}
\title{Build a bivariate colour key}
\usage{
build_bkey(
  data,
  palette = "BlueYellow",
  terciles = FALSE,
  flipAxis = FALSE,
  expertR_est = NA,
  expertR_err = NA,
  bound = NULL
)
}
\arguments{
\item{data}{A data frame.}

\item{palette}{Name of colour palette or character vector of hex colour codes
created with \code{\link{build_palette}} function. Colour palette names
include \code{BlueYellow}, \code{CyanMagenta}, \code{BlueRed} and
\code{GreenBlue}.}

\item{terciles}{A logical value. This provides the option to define numerical
bounds for the colour key grid using terciles instead of equal intervals.}

\item{flipAxis}{A logical value. Whether to place the axis on the opposite
sides or not.}

\item{expertR_est}{A vector consisting of the range of expert values for the
estimate (Default: NA).}

\item{expertR_err}{A vector consisting of the range of expert values for the
error (Default: NA).}

\item{bound}{A vector of eight elements representing the bounds for the
estimate and error that will be used on the bivariate colour key.  These can
be created offline using the \code{\link{findNbounds}} function.}
}
\description{
This function links data and a colour palette. Numerical bounds are added to
the 3 x 3 colour grid.
}
\details{
A key for the bivariate map is not automatically generated with
\code{\link{build_bmap}} and must be made using \code{\link{build_bkey}}. It
is important that the arguments passed to \code{\link{build_bkey}} match those
passed to \code{\link{build_bmap}}. The map and key can be viewed together
using \code{\link{attach_key}}.

If bound is NULL, the bounds for the legend will be computed from the
 data and expert bounds (if available) using the findNbound function.  The
 argument bound should only be used if you want to make comparison against
 multiple maps.  For that scenario, the user should use the findNbounds
 function to generate a bound for the larger set of data that wish to
 compare.
}
\examples{
data(us_data)
data(us_geo)
poverty <- read.uv(data = us_data, estimate = "pov_rate", error = "pov_moe")

# use a prepared palette and terciles
key <- build_bkey(data = poverty, terciles = TRUE)
view(key)

# use a created palette
p <- build_palette(name = "usr", colrange =
 list(colour = c("darkblue", "chartreuse4"), difC = c(3, 4)))
key <- build_bkey(data = poverty, palette = p)
view(key)
}
\seealso{
\code{\link{attach_key}}
}
