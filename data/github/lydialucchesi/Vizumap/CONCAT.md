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







