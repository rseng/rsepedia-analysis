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
[![Build Status](https://travis-ci.org/ropensci/iheatmapr.svg?branch=master)](https://travis-ci.org/ropensci/iheatmapr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/iheatmapr?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/iheatmapr)
[![codecov](https://codecov.io/gh/ropensci/iheatmapr/branch/master/graph/badge.svg?token=CTupoUlXNI)](https://codecov.io/gh/ropensci/iheatmapr)
![R version](https://img.shields.io/badge/R%20%3E%3D-3.2.0-blue.svg)
[![](https://badges.ropensci.org/107_status.svg)](https://github.com/ropensci/onboarding/issues/107)
[![JOSS](http://joss.theoj.org/papers/10.21105/joss.00359/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00359)
[![CRAN](https://www.r-pkg.org/badges/version/iheatmapr)](https://cran.r-project.org/package=iheatmapr)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# iheatmapr

`iheatmapr` is an R package for building complex, interactive heatmaps using modular building blocks. "Complex" heatmaps are heatmaps in which subplots along the rows or columns of the main heatmap add more information about each row or column. For example, a one column additional heatmap may indicate what group a particular row or column belongs to. Complex heatmaps may also include multiple side by side heatmaps which show different types of data for the same conditions. Interactivity can improve complex heatmaps by providing tooltips with information about each cell and enabling zooming into interesting features. `iheatmapr` uses the [plotly](https://plotly.com) library for interactivity. 

While there are already plenty of awesome R packages for making heatmaps, including several great packages for making relatively simple interactive heatmaps ([heatmaply](https://github.com/talgalili/heatmaply) and [d3heatmap](https://github.com/rstudio/d3heatmap)) or complex static heatmaps ([ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap)), `iheatmapr` seeks to make it easy to make complex interactive heatmaps. 

## Installation

To install the CRAN version of `iheatmapr`:

```r
install.packages("iheatmapr")
```

To install the github version of `iheatmapr`:

```r
devtools::install_github("ropensci/iheatmapr")
```

## Example Complex Heatmap

As an example of a complex heatmap, we can make a version of the famous vaccines plot from the Wall Street Journal that has been recreated in several other heatmap frameworks in R. 

![](https://raw.githubusercontent.com/ropensci/iheatmapr/master/vaccine.gif)

The code to create this heatmap is:

```R
library(iheatmapr)
data(measles, package = "iheatmapr")

main_heatmap(measles, name = "Measles<br>Cases", x_categorical = FALSE,
             layout = list(font = list(size = 8))) %>%
  add_col_groups(ifelse(1930:2001 < 1961,"No","Yes"),
                  side = "bottom", name = "Vaccine<br>Introduced?",
                  title = "Vaccine?",
                  colors = c("lightgray","blue")) %>%
  add_col_labels(ticktext = seq(1930,2000,10),font = list(size = 8)) %>%
  add_row_labels(size = 0.3,font = list(size = 6)) %>% 
  add_col_summary(layout = list(title = "Average<br>across<br>states"),
                  yname = "summary")  %>%                 
  add_col_title("Measles Cases from 1930 to 2001", side= "top") %>%
  add_row_summary(groups = TRUE, 
                  type = "bar",
                  layout = list(title = "Average<br>per<br>year",
                                font = list(size = 8)))
              
```

Modular components of the plot are added in an iterative fashion to the top, right, left, or bottom of the heatmap. `iheatmapr` also contains a function (`iheatmap`) to make a fairly standard heatmap with optional dendrograms and row or column annotation heatmaps (See vignette).  

All the plots aligned with the main heatmap horizontally share the same y axis and thus zooming in the y direction within the heatmap will also zoom in to those subplots. The plots aligned vertically share an x axis with that heatmap and zooming horizontally within those plots will be linked.  

Hovering over the heatmaps yields a tooltip with the name of the row and column as well as the value represented.

# Documentation

See the [vignette](https://docs.ropensci.org/iheatmapr/articles/full_vignettes/iheatmapr.html) for a more thorough introduction to the package.

# Acknowledgements

This package includes the open source Plotly.js library, which does much of the work of making these interactive plots possible! In creating this package, I also drew inspiration & some code from the great plotly R package; in particular, the code for the `iheatmapr` htmlwidget is adapted from an earlier version of the plotly R package. Additionally, numerous people at Genentech helped provide feedback and guidance for this project, including but not limited to Justin Finkle, August Guang, Michael Lawrence, Gabe Becker, Steve Lianoglou, Pete Haverty... thanks to all who helped review code and/or provide feedback!  This package also went through the on-boarding process for rOpensci -- thanks to the reviewers Carl Ganz and Andee Kaplan and editor Maëlle Salmon for all their helpful feedback! 

[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# iheatmapr 0.5.1 (CRAN Release)

* Bug fix for column annotation labels
* Bug fix for issues with merging colorbars

# iheatmapr 0.5.0 (CRAN Release)

* Remove dependencies on S4Vectors (problematic because it is a Bioconductor package) and plyr 
* Bug fix for single row or column inputs for main heatmap
 
# iheatmapr 0.4.12 (CRAN Release)

* Adjust tests to be compatible with newer version of scales package.

# iheatmapr 0.4.11

* Fix issue for some subplots where an NA value prevents any plotting.

# iheatmapr 0.4.10

* Exports to_plotly_list and adds new to_plotly_json to make it easier to convert to plotly spec

# iheatmapr 0.4.9

* Adds an option to main_heatmap to show colorbar or not, by [@fboehm](https://github.com/fboehm)

# iheatmapr 0.4.8

* Update to fix error in tests due to change in R-devel random sample
* Fix issue where add_col_groups uses continuous scale if all groups are equal
* Fix for iheatmapr_event not working on relayout 
* Fix issue where show_colorbar was always set to be TRUE
* Fix example for save_iheatmap

# iheatmapr 0.4.3

* Minor bug fix for issue on R-devel.

# iheatmapr 0.4.2

* Update tests for compatibility with new version of testthat

# iheatmapr 0.4.1

* Bug fix -- actually show titles for add_row_groups, add_col_groups, etc.

# iheatmapr 0.4.0

* Added option to customize text in tooltips for heatmaps. By default, only show
three sig figs for values.  This change will break rendering of Iheatmap objects
created via older versions of iheatmapr.

# iheatmapr 0.3.0

* Removed plotly r package dependency.  iheatmapr now directly interfaces with plotlyjs instead. While this may mean having two copies of the plotly.js library, the benefit is that the htmlwidget for iheatmapr can be customized for plots created by iheatmapr and that iheatmapr now has more control over when the js component gets updated.

# iheatmapr 0.2.1

* Added a `NEWS.md` file to track changes to the package.



# Contributing to iheatmapr

## Opening Issues

Please report bugs and contribute feature requests using the the Github Issues page. For reporting bugs, please provide a reproducible example. For feature requests, please give an example use case.

## Development guidelines

New features or improvements to the existing codebase or documentation are welcomed. Please use GitHub's pull request feature for proposing and submitting changes.  

New features or changes to existing features should mimic the style used for existing features. For example, if adding a new modular component that can be added onto a complex heatmap, the name of the new function should start with "add_". Additionally, iheatmapr uses the S4 OOP system, so additional functions should generally be written as S4 methods. 

If adding a new feature, a test should be added for that new feature, as well as an update to the vignette to document the new feature.  The website documentation should also be updated using pkgdown.  

Testing is done using the testthat package. A `expect_iheatmap` function is included in helper_expectation.R in the tests/testthat directory. This function wraps a few expectations.  In particular, the expectation will create a saved version of the data that is used to create the plotly graphic. After being run once, the test will check that the data is the same as was previously generated. The view_reference.Rmd file will create the plots based on all the saved data created by the `expect_iheatmap` function. Building view_reference.Rmd can be used to visually inpect the expected result of each test.

## Scope

iheatmapr is intended to be a general purpose package for creating complex, interactive heatmaps. The modular building blocks can be adapted to make domain-specific types of heatmaps, but such specialized adaptations are a better fit for complementary, add-on packages that build upon iheatmapr rather than as components of iheatmapr itself. 

## Questions about these guidelines?

Please use the Issues page for questions about these guidelines. You can also submit a preliminary pull request with questions. 

## Code of Conduct

When contributing to iheatmapr, you are expected to follow the [code of conduct](https://github.com/AliciaSchep/iheatmapr/blob/master/CONDUCT.md).
## Release Overview

Bug fixes based on github issues.

Tested on Windows, Mac, and Linux via Travis, Appveyor, and/or Rhub.

## R CMD check results

0 errors, 0 warnings, 0 notes

## Reverse dependencies

No reverse dependencies on CRAN, checked reverse dependency on Bioconductor (lipidr). 
------
  title: 'iheatmapr: Interactive complex heatmaps in R'
  tags:
    - visualization
    - R
    - heatmap
  authors:
   - name: Alicia N Schep
     orcid: 0000-0002-3915-0618
     affiliation: 1,2
   - name: Sarah K Kummerfeld
     orcid: 0000-0002-0089-2358
     affiliation: 2
  affiliations:
   - name: Stanford University
     index: 1
   - name: Genentech
     index: 2
  date: 2 March 2017
  bibliography: paper.bib
---

# Summary

The iheatmapr package is an R package [@team2000r] for creating complex, interactive heatmaps. Heatmaps are commonly used to visualize patterns in high-dimensional data, particularly in genomics research. Adding annotations and summary plots along the rows or columns of a heatmap can enhance a heatmap visualization. Pairing several heatmaps horizontally can enable association of data across multiple assays. Adding interactive features like tooltips and zooming can further enhance these complex heatmaps, enabling information-rich visualizations linking diverse, high-dimensional data sets. 

There are great tools in R for creating simple interactive heatmaps [@heatmaply, @d3heatmap] or creating static complex heatmaps [@ComplexHeatmap]. However, there are no tools facilitating the creation of complex, interactive heatmaps. The iheatmapr package fills this gap, enabling creation of highly customizable, interactive, complex heatmaps using the plotly library [@plotly]. The resulting interactive visualizations can easily be incorporated into reproducible R Markdown reports for sharing with collaborators. 
  
# References