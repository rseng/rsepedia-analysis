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
  
# References---
title: "Reference plots"  
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(iheatmapr)
```


```{r, echo = FALSE,  results='asis'}
out <- htmltools::tagList()
ref_files <- list.files("testthat/reference/","*.rds", full.names = TRUE)
for (i in seq_along(ref_files)){
   out[[(i*2 - 1)]] <- strsplit(basename(ref_files[i]),".rds")[[1]]
   tmp <- readRDS(ref_files[i])
   out[[i*2]] <- tmp
}

out
```
 
```{r}
Sys.time()
```

---
title: "Minimal iheatmpar vignette"
author: "Alicia Schep"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Minimal iheatmpar vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduction to `iheatmapr` package

`iheatmapr` is an R package for building complex, interactive heatmaps using modular building blocks. "Complex" heatmaps are heatmaps in which subplots along the rows or columns of the main heatmap add more information about each row or column. For example, a one column additional heatmap may indicate what group a particular row or column belongs to. Complex heatmaps may also include multiple side by side heatmaps which show different types of data for the same conditions. Interactivity can improve complex heatmaps by providing tooltips with information about each cell and enabling zooming into interesting features. `iheatmapr` uses the [plotly](https://plotly.com) library for interactivity. 

While there are already plenty of awesome R packages for making heatmaps, including several great packages for making relatively simple interactive heatmaps (e.g. [heatmaply](https://github.com/talgalili/heatmaply) and [d3heatmap](https://github.com/rstudio/d3heatmap)) or complex static heatmaps ([ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap)), `iheatmapr` seeks to make it easy to make complex interactive heatmaps. 

# Usage documentation

For an overview of how to use `iheatmapr`, please see the [documentation website](https://docs.ropensci.org/iheatmapr).---
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{pkgdown}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(fig.width = 4.5, fig.height = 3.5)
```


# Introduction

`iheatmapr` is an R package for building complex, interactive heatmaps using modular building blocks. "Complex" heatmaps are heatmaps in which subplots along the rows or columns of the main heatmap add more information about each row or column. For example, a one column additional heatmap may indicate what group a particular row or column belongs to. Complex heatmaps may also include multiple side by side heatmaps which show different types of data for the same conditions. Interactivity can improve complex heatmaps by providing tooltips with information about each cell and enabling zooming into interesting features.  

While there are already plenty of awesome R packages for making heatmaps, including several great packages for making relatively simple interactive heatmaps ([heatmaply](https://github.com/talgalili/heatmaply) and [d3heatmap](https://github.com/rstudio/d3heatmap)) or complex static heatmaps ([ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap)), `iheatmapr` seeks to make it easy to make complex interactive heatmaps.

For this vignette, we will use the Indometh data included in the datasets package. This data set contains data on the pharmacokinetics of Indometacin, specifically the observed plasma concentration of indometacin at various time points after intravenous injection for 6 different patients. We will first cast the data into a matrix of concentrations by patients (rows) and time (columns).  We will also compute correlation matrices for the patients, and the maximum and minimum concentration of Indometacin per patient. Finally, we will assign patients into some groups -- in this dataset no extra patient information is provided and so we'll just make arbitrary groups that are intended to show how one might use actual groupings like gender or ethnicity. 

```{r, message = FALSE}
library(iheatmapr)
library(datasets)
library(reshape2)

Indometh_matrix <- acast(Indometh, Subject ~ time, value.var = "conc")
Indometh_matrix <- Indometh_matrix[as.character(1:6),]
rownames(Indometh_matrix) <- paste("Patient",rownames(Indometh_matrix))
Indometh_patient_cor <- cor(t(Indometh_matrix))

patient_max_conc <- apply(Indometh_matrix,1,max)
patient_min_conc <- apply(Indometh_matrix,1,min)
patient_groups <- c("A","A","B","A","B","A") # Arbitrary groups
```

## Example Complex Heatmap

We will start off by showing an example of the type of complex heatmap we can create using the `iheatmapr` package. In our example complex heatmap, we will plot:

* a heatmap of the correlation matrix 
* a dendrogram for column clustering
* a dendrogram for row clustering as well as annotation heatmap showing 3 clusters
* annotation heatmaps showing the max and min Indometacin concentrations per patient
* annotation heatmap showing whether patients fall in group "A" or "B"
* a heatmap showing the actual Indometacin timecourse data for each patient
* a plot showing the average Indometacin timecourse response
* titles and labels for heatmaps


```{r, fig.width = 7, fig.height = 4}
main_heatmap(Indometh_patient_cor,name = "Correlation") %>%
  add_col_clustering() %>%
  add_row_clustering(k = 3) %>%
  add_row_title("Patients") %>%
  add_col_title("Patients") %>%
  add_row_annotation(data.frame("Max" = patient_max_conc,
                                "Min" = patient_min_conc,
                                "Groups" = patient_groups)) %>%
  add_main_heatmap(Indometh_matrix,
                   name = "Indometacin<br>Concentration") %>%
  add_col_labels() %>%
  add_col_title("Time") %>%
  add_col_summary()

```

All the plots aligned with the correlation heatmap horizontally share the same y axis and thus zooming in the y direction within the heatmap will also zoom in to those subplots. The plots aligned with either the correlation heatmap or the concentration heatmap vertically share an x axis with that heatmap and zooming horizontally within those plots will be linked.  

Hovering over the heatmaps yields a tooltip with the name of the row and column as well as the value represented.

## Modular system for building up complex heatmaps

In our example complex heatmap above, the various components are added iteratively. There are two main types of plots in `iheatmapr`-- "main heatmaps" and subplots. Subplots can be added to the left, right, top, or bottom of a main heatmap. We will go over the basics of adding together different components of a plot in this section, and the [Modular buidling blocks](#modular-building-blocks) section of the vignette will go into detail into each of the modular components available.   

### Starting with basic heatmap

To initialize a complex heatmap, we start off by making a very simple heatmap using the `main_heatmap` function.  

```{r}
example_heatmap <- main_heatmap(Indometh_patient_cor, name = "Correlation")

example_heatmap
```

### Adding subplots to right or left

To this heatmap, we can add subplots to the top, bottom, left, or right. To illustrate how subplots are added we will start with one kind of subplot -- annotation heatmaps. The function `add_row_annotation` adds annotation heatmap(s) to the side of a heatmap. These functions are most easily chained together using the `%>%` operator from the magrittr package.  

```{r}
example_heatmap <- example_heatmap %>% 
  add_row_annotation(data.frame("Groups" = patient_groups))

example_heatmap
```

After adding one subplot to the right of a heatmap, we can easily add another. We can add as many subplots to any side of the heatmap as we want.    

```{r}
example_heatmap <- example_heatmap %>% 
  add_row_annotation(data.frame("Groups" = patient_groups))

example_heatmap
```

We can also add the subplot to the left instead of the right by specifying `side = left`.

```{r}
example_heatmap <- example_heatmap %>% 
  add_row_annotation(data.frame("Groups" = patient_groups), side = "left")

example_heatmap
```

### Adding subplot to top or bottom 

The analogous function add_col_annotation can be used for adding annotation heatmap(s) to the top or bottom of a heatmap. We will first add one to the top (default side).

```{r}
example_heatmap <- example_heatmap %>% 
  add_col_annotation(data.frame("Groups" = patient_groups))

example_heatmap
```

Now we add one to the bottom:

```{r}
example_heatmap <- example_heatmap %>% 
  add_col_annotation(data.frame("Groups" = patient_groups), side = "bottom")

example_heatmap
```

The axis titles seem to overlap a bit -- to avoid this we can potentially alter the spacing between the plots -- see [Altering plot sizes and spacing](#altering-plot-sizes-and-spacing).  

### Additional main heatmaps

We can add another "main heatmap" using the `add_main_heatmap` function.

```{r}
example_heatmap <- example_heatmap %>% 
  add_main_heatmap(Indometh_patient_cor, name = "Correlation")

example_heatmap
```

Once we've added another main heatmap, when we add a subplot to the top or bottom, by default it will go on top/bottom of the most recently added heatmap.

```{r}
example_heatmap <- example_heatmap %>% 
  add_col_annotation(data.frame("Groups" = patient_groups))

example_heatmap
```

To specify that we want to add the subplot to another main heatmap, the same "xname" argument can be passed both to the desired main heatmap and the subplot to be added. 

### Row and Column ordering

When we initialize a heatmap using `main_heatmap` we can specify the row and/ or column order using the `row_order` and `col_order` arguments respectively. The functions `add_row_clustering` or `add_col_clustering` (See [Add clustering](#add-clustering)) can also be used to order the rows or columns. When the row or column order is specified (either initially or through a later function), this ordering is carried through to all plots along the same axis. 

## The iheatmap and add_iheatmap functions

The `iheatmapr` package includes two functions, `iheatmap` and `add_iheatmap`, that wrap together several of the modular building blocks that are often used in conjunction to make common types of heatmaps. Rather than calling many of the modular functions individually, additional arguments can be passed to `iheatmap` or `add_iheatmap`. Not all modular subcomponents can be specified this way, but such components can still be added directly through their function calls (e.g. `add_col_summary`). The following code generates the same plot as our first example heatmap above:

```{r, fig.width = 7, fig.height = 4}
iheatmap(Indometh_patient_cor, 
         cluster_cols = "hclust", 
         cluster_rows = "hclust",
         col_title = "Patients",
         row_title = "Patients",
         name = "Correlation",
         row_k = 3,
         row_annotation = data.frame("Max" = patient_max_conc,
                                     "Min" = patient_min_conc,
                                     "Groups" = patient_groups)) %>%
  add_iheatmap(Indometh_matrix,
               name = "Indometacin<br>Concentration", # html <br> tag used to split lines
               col_title = "Time",
               col_labels = TRUE) %>% 
  add_col_summary()
```

These functions can be useful for quickly creating a fairly standard plot. However, fewer options are available with `iheatmap` than if using each of the modular functions by themselves and there is much less flexibility in the ordering and sizing of the different subcomponents. 

## Saving plots

The `save_iheatmap` function can be used to save an iheatmap object as a standalone html file or as a static pdf/png/jpeg.

```{r, eval = FALSE}
myplot <- iheatmap(Indometh_patient_cor, 
         cluster_cols = "hclust", 
         cluster_rows = "hclust",
         col_title = "Patients",
         name = "Correlation") 

myplot %>% save_iheatmap("myplot.html") # Save interactive HTML
myplot %>% save_iheatmap("myplot.pdf") # Save static plot (pdf, png, or jpeg)
```


# Modular building blocks

In this section of the vignette, we will explore each of the modular components included in the `iheatmapr` package.  

## Main heatmap

The `main_heatmap` function creates a basic heatmap. By default, the column and row names of the input matrix are used for determining the "Row" and "Column" labels in the tooltip above each cell.   

```{r}
main_heatmap(Indometh_patient_cor, 
             name = "Correlation")
```

Alternatively one can provide the "y" or "x" arguments to use different labels in the tooltip. Providing these arguments will also affect the row and column labels if those are added later (See [Add labels](#add-labels))


```{r}
patient_names <- paste0("Patient ",seq_len(ncol(Indometh_patient_cor)))

main_heatmap(Indometh_patient_cor, 
             x = patient_names,
             y = patient_names,
             name = "Correlation")
```

To change the color use the "colors" argument with either the name of an RColorBrewer palette or a vector of colors:

```{r}
main_heatmap(Indometh_patient_cor,
             colors = "Blues",
             name = "Correlation")
```


## Add main heatmap

The `add_main_heatmap` function is very similar to the main heatmap function, it just adds a new main heatmap to the left or right of the existing plots.  

```{r}
main_heatmap(Indometh_patient_cor, name = "Correlation") %>%
  add_main_heatmap(Indometh_matrix, name = "Indometacin<br>Concentration")
```

Note that if the same "name" argument is given to both the first main heatmap and the additional main heatmap, the two heatmaps will share the same colorscale and colorbar. 

See [Orientation](#orientation) for adding main heataps vertically instead of horizontally!

## Add labels

Row and column labels are not shown by default. To add them, use the `add_row_labels` or `add_col_labels`.  By default, these will use the row and column names provided to the `main_heatmap` or `add_main_heatmap` either directly or indirectly via the row and column names of the matrix. 

```{r}
main_heatmap(Indometh_matrix, name = "Correlation") %>%
  add_row_labels() %>% 
  add_col_labels() %>% 
  add_row_title("Patients") %>% 
  add_col_title("Patients")
```

## Add clustering 

The functions `add_row_clustering` and `add_col_clustering` will cluster the rows or columns. By default, the clustering is via hierarchical clustering and a dendrogram is added. 

```{r}
main_heatmap(Indometh_patient_cor) %>% 
  add_row_clustering() %>% 
  add_col_clustering()
```

Specifying a value of k will result in cluster assignments being made based on the dendrogram and an annotation heatmap gets added showing the clustering. 

```{r}
main_heatmap(Indometh_patient_cor) %>% 
  add_row_clustering(k = 3) %>% 
  add_col_clustering()
```

Alternatively, clustering can be done using kmeans by setting `method = "kmeans"`. In that case, a value of k must be given!

```{r}
main_heatmap(Indometh_matrix) %>% 
  add_row_clustering(k = 3, method = "kmeans") 
```

Another way of clustering the results is by giving group assignments. 

```{r}
main_heatmap(Indometh_matrix) %>% 
  add_row_clustering(method = "groups", groups = c("A","A","B","B","A","A")) 
```

### Using your own clustering

The functions `add_row_clustering` and `add_col_clustering` will perform clustering based on the options given. However, you might in some cases want more control over the clustering or be able to ensure the clustering matches that from some other analysis or visualization. To give you more control over the clustering, there are two lower-level functions that take in a clustering result and simply add them to the plot.  

For hierarchical clustering, you can use the `add_row_dendro` or `add_col_dendro` functions to add the dendrogram. The clustering result should be provided as an "hclust" object.   

```{r}
clust_res <- hclust(as.dist(1 - Indometh_patient_cor))

main_heatmap(Indometh_patient_cor) %>% 
  add_row_dendro(clust_res) %>%
  add_col_dendro(clust_res)
```

If you want to add a dendrogram but have already ordered your matrix, you can set the "reorder" argument to `FALSE`.

For adding cluster assignments from a method like K-means, the `add_col_clusters` and `add_row_clusters` functions can be used. These functions will take in a vector of cluster assignments and both add an annotation heatmap showing the assignments and also re-order the rows or columns based on the clustering.  

```{r}
clust_assign <- kmeans(Indometh_matrix, 3)$cluster

main_heatmap(Indometh_patient_cor) %>% 
  add_row_clusters(clust_assign) %>%
  add_col_clusters(clust_assign)
```

The `add_row_clusters` and `add_col_clusters` methods are very similar to `add_col_groups` and `add_col_groups`; the main difference is that `add_row_clusters` and `dd_col_clusters` will also reorder the rows and columns in addition to adding the annotation heatmap.

## Add annotations

The functions `add_row_annotation` and `add_col_annotation` add one or more annotation heatmaps. Annotations should be provided as a data.frame or something that can be coerced into a data.frame.  
```{r, fig.width = 6}

main_heatmap(Indometh_patient_cor) %>% 
  add_row_annotation(data.frame("Max" = patient_max_conc,
                                "Min" = patient_min_conc,
                                "Groups" = c("A","A","B","B","A","A")))

```

By default, colors will be chosen for each annotation. To assign colors yourself, provide a list of colors, with the names of the list matching the column names of the annotation. Colors can either be the name of an RColorBrewer palette or a vector of colors.  

```{r, fig.width = 6}

main_heatmap(Indometh_patient_cor) %>% 
  add_row_annotation(data.frame("Max" = patient_max_conc,
                                "Min" = patient_min_conc,
                                "Groups" = c("A","A","B","B","A","A")),
                     colors = list("Max" = "Reds",
                                   "Min" = "Blues",
                                   "Groups" = c("purple","pink")))

```

For more control over annotation heatmaps, use the functions `add_col_signal` and `add_row_signal` for adding a single continuous annotation or `add_col_groups` and `add_row_groups` for adding a single discrete annotation.

```{r, fig.width = 6}
main_heatmap(Indometh_patient_cor) %>% 
  add_row_signal(patient_max_conc, "Max<br>Concentration", title = "Max", colors = "Reds") %>%
  add_row_signal(patient_min_conc, "Min<br>Concentration", title = "Min", colors = "Reds") %>%
  add_row_groups(c("A","A","B","B","A","A"), "Groups") 

```

## Add summary

One type of subplot that can be added is a summary of the values over the rows or the columns. For example when plotting Indometh concentrations per patient over time we might want to plot the average time response using `add_col_summary`.  

```{r}
main_heatmap(Indometh_matrix) %>% 
  add_col_summary() 
```

`add_row_summary` is similar but adds plot along left summarizing rows.

With `add_col_summary` or `add_row_summary` we can pass groups to divide by the rows or columns respectively when computing the summary:

```{r}
main_heatmap(Indometh_matrix) %>% 
  add_col_summary(groups = c("A","A","B","B","A","A")) 
```

If groups is set to `TRUE` rather than a vector of groups, then the function will try to use an existing set of groups or clusters that have been added to the plot. This should be avoided if more than one set of groups has been added (only one of the existing sets of groups will be used in that case.)

```{r}
main_heatmap(Indometh_matrix) %>% 
  add_row_clustering(k = 3) %>% 
  add_col_summary(groups = TRUE) 
```

By default, summarization is done using `mean`. However, the `summary_function` argument can be set to 'median', 'sd', 'var', 'mad', 'max', 'min', or 'sum' in order to perform any of those alternate summarizations (Note: 'sum' only introduced in v0.4.4). We can also use the `layout` argument to pass a title to the new axis.    

```{r}
main_heatmap(Indometh_matrix) %>% 
  add_row_clustering(k = 3) %>% 
  add_col_summary(summary_function = "sd", 
                  layout = list(title = "sd")) 
```


## Add barplot

The functions `add_col_barplot` and `add_row_barplot` add barplots along the columns or rows of a main heatmap.  

```{r}
main_heatmap(Indometh_matrix) %>% 
  add_col_barplot(y = as.numeric(colnames(Indometh_matrix)),
                  tracename = "time", 
                  layout = list(title = "Time"))
```

## Add line plot

The functions `add_col_barplot` and `add_row_barplot` add a line plot along the columns or rows of a main heatmap.  

```{r}
main_heatmap(Indometh_matrix) %>% 
  add_col_plot(y = as.numeric(colnames(Indometh_matrix)),
                  tracename = "time", 
                  layout = list(title = "Time"))
```

## Add arbitrary plot

The functions for adding barplots or lineplots assume one data point per row or column. To add an arbitrary plot, one can use `add_subplot`.  

```{r}
main_heatmap(Indometh_matrix) %>% 
  add_subplot(x = 1:3, y = 4:6, side = "top")
```

# Customization options

The next section of this vignette will cover some customization options available in `iheatmapr`.

## Color selection

`iheatmapr` tries to choose colors automatically in a sensible way. If a heatmap with the same "name" argument as another heatmap is added, the two heatmaps will share the same colorbar and scale. Functions to make and add heatmaps and subplots have a "colors" argument that can be used to specify colors manually. These arguments take in RColorBrewer palette names or vectors of colors. 

For example, if we simply repeat the same exact main heatmap with the same "name" argument and the same column names for the annotation matrices, the colorbars are not duplicated.  

```{r, fig.width = 6}
main_heatmap(Indometh_patient_cor, name = "Correlation") %>%
  add_row_annotation(data.frame("Max" = patient_max_conc,
                                "Min" = patient_min_conc,
                                "Groups" = c("A","A","B","B","A","A")),
                     colors = list("Max" = "Reds",
                                   "Min" = "Blues",
                                   "Groups" = c("purple","pink"))) %>% add_main_heatmap(Indometh_patient_cor, name = "Correlation") %>%
  add_row_annotation(data.frame("Max" = patient_max_conc,
                                "Min" = patient_min_conc,
                                "Groups" = c("A","A","B","B","A","A")),
                     colors = list("Max" = "Reds",
                                   "Min" = "Blues",
                                   "Groups" = c("purple","pink")))

```


## Colorbars

When additional heatmaps or annotations with colorscales are added to an iheatmap object, the colorbars will be added according a grid pattern.  The layout and spacing of that grid can be specified using the "colorbar_grid" parameter to iheatmap or simple_heatmap. By default, the grid has three rows. To change the layout of the grid to only have 2 rows and adjust the spacing and sizing of the colorbars, we can use the function setup_colorbar_grid to create a ColorbarGridParameters object that can be passed to the colorbar_grid argument to simple_heatmap.

```{r, fig.width = 7, fig.height = 4}
grid_params <- setup_colorbar_grid(nrows = 2, 
                                   y_length = 0.3, 
                                   x_spacing = 0.2,
                                   y_spacing = 0.5, 
                                   x_start = 1.1, 
                                   y_start = 1)


iheatmap(Indometh_patient_cor, 
         cluster_cols = "hclust", 
         cluster_rows = "hclust",
         col_title = "Patients",
         name = "Correlation",
         row_k = 3,
         row_annotation = data.frame("Max" = patient_max_conc,
                                     "Min" = patient_min_conc,
                                     "Groups" = patient_groups),
         colorbar_grid = grid_params) %>%
  add_iheatmap(Indometh_matrix,
               name = "Indometacin<br>Concentration", # html <br> tag used to split lines
               col_title = "Time",
         row_title = "Patients") %>% 
  add_col_summary()
```


## Continuous Axes

In all the plots above, the axes of the main heatmaps have been treated as categorical, with the columns of equal widths and rows of equal widths. In some cases, it may make sense to treat one or both axes as continuous. There are several ways to indicate to `iheatmapr` that the axes should be continuous. The first is to use the `x_categorical` and `y_categorical` arguments, setting them to false.  

```{r}
main_heatmap(Indometh_matrix, x_categorical = FALSE)
```

By default, names of rows and columns are taken from the row and column names of the matrix, respectively. Alternatively, these can be provided to the x and y arguments to `main_heatmap` or `iheatmap`. If the values are given as numeric and do not correspond simply to the row or column number, then the axis is assumed to be continuous. 

```{r}
main_heatmap(Indometh_matrix, x = as.numeric(colnames(Indometh_matrix)))
```

### Labels

Labels will be a bit different for continuous axes -- rather than labelling every row or column,
labels are chosen more similarly to a regular plot.

```{r}
main_heatmap(Indometh_matrix, x_categorical = FALSE) %>%
  add_col_labels()
```

### Clustering

Another major difference with categorical axes is that continuous axes can't be clustered!

```{r, error = TRUE}
main_heatmap(Indometh_matrix, x_categorical = FALSE) %>%
  add_col_clustering() # will give error
```


## Orientation

By default, additional heatmaps are added horizontally to left or right of existing heatmap(s). It is also possible to build up complex heatmaps vertically rather than horizontally by setting `orientation = vertical` to the initial call to `iheatmap` or `main_heatmap`.

```{r, fig.width = 4, fig.height = 7}
iheatmap(Indometh_patient_cor, 
         cluster_cols = "hclust", 
         cluster_rows = "hclust",
         row_title = "Patients",
         name = "Correlation",
         orientation = "vertical") %>%
  add_iheatmap(t(Indometh_matrix),
               name = "Indometacin<br>Concentration", # html <br> tag used to split lines
               row_title = "Time",
         col_title = "Patients") 
```


## Margins and other layout properties

General plotly layout arguments can be passed to the `layout` argument to `main_heatmap` or `iheatmap`. See the [plotly documentation](https://plotly.com/javascript/reference/#layout) for more information about the various options.

```{r}
main_heatmap(Indometh_matrix, layout = list(margin = list(b = 120)))
```

The layout can also be adjusted after creation of the plot using the `modify_layout` function. 

```{r}
main_heatmap(Indometh_matrix) %>% modify_layout(list(margin = list(b = 120)))
```


## Axis properties

Most functions for adding a subplot include a "layout" argument for providing layout parameters for the additional axis that is created.  This layout parameter accepts a list of layout parameters that plotly uses for altering the layout of an axis.  See the [plotly documentation](https://plotly.com/javascript/reference/#layout-xaxis) for more information about the various options.

As an example, here is a plot with default axis options:

```{r}
iheatmap(Indometh_matrix,
               name = "Indometacin<br>Concentration", # html <br> tag used to split lines
               col_title = "Time",
         row_title = "Patients") %>% 
  add_col_summary()
```

Adding `zeroline = FALSE` removes the zero line, and `title = "Average"` adds yaxis label.  

```{r}
iheatmap(Indometh_matrix,
               name = "Indometacin<br>Concentration", # html <br> tag used to split lines
               col_title = "Time",
         row_title = "Patients") %>% 
  add_col_summary(layout = list(zeroline = FALSE, title = "Average"))
```

## Sharing axes

By default subplots stacked vertically share the X axis with a main heatmap, and subplots stacked horizontally share their Y axis with a main heatmap. Subplots can also be made to share the other axis with an existing subplot. 

Here is an example of two subplots with independent axes:

```{r}
main_heatmap(Indometh_matrix) %>%
  add_col_summary() %>%
  add_main_heatmap(Indometh_matrix) %>%
  add_col_summary()
```

By passing the same "yname" argument to both axes, the axes become shared:

```{r}
main_heatmap(Indometh_matrix) %>%
  add_col_summary(yname = "Summary") %>%
  add_main_heatmap(Indometh_matrix) %>%
  add_col_summary(yname = "Summary")
```

## Altering plot sizes and spacing

Each type of subplot has a default relative size and default relative spacing from the previous plots, but both the size and spacing can be altered by passing a 'buffer' or 'size' argument. Both arguments are relative to the size of the first main heatmap along that dimension.  

```{r}
main_heatmap(Indometh_matrix) %>%
  add_col_summary()
```

```{r}
main_heatmap(Indometh_matrix) %>%
  add_col_summary(buffer = 0.2, size = 1)
```

## Sizing in knitr

To change the sizing within knitr, use the fig.width and fig.height chunk options.

# Shiny 

Iheatmap objects can be used in shiny.  The `renderIheatmap` function should be used in the server code, and the `iheatmaprOutput` in the ui code.  

For observing events, the `iheatmapr_event` function should be used.  The function takes as the first argument the Iheatmap object, and as the second argument a type of event -- "click","hover", or "relayout". 

To see the output from `iheatmapr_event`, you can use the `shiny_test` function from `iheatmapr`.  

```{r, eval = FALSE}
hm <- main_heatmap(Indometh_patient_cor,name = "Correlation") %>%
  add_col_clustering() %>%
  add_row_clustering(k = 3) %>%
  add_row_title("Patients") %>%
  add_col_title("Patients") %>%
  add_row_annotation(data.frame("Max" = patient_max_conc,
                                "Min" = patient_min_conc,
                                "Groups" = patient_groups)) %>%
  add_main_heatmap(Indometh_matrix,
                   name = "Indometacin<br>Concentration") %>%
  add_col_labels() %>%
  add_col_title("Time") %>%
  add_col_summary()

## NOT RUN (runs shiny app)
test_iheatmapr_event(hm, "click")
```

To see the code for the demo shiny app:

```{r}
test_iheatmapr_event
```

## Deploying on shinyapps.io

An issue that has arisen in deploying to shinyapps.io is the use case of adding a download button, e.g.

```{r, eval=FALSE}
# Not evaluated -- example output code
output$download_test<-downloadHandler(
            filename ="test_heatmap.png",
            content = function(file){
                    save_iheatmap(heatmap_reactive(),file,vwidth=2000,vheight=1000)
            },
            contentType = "image/png"
)
```

In this case, the app needs to include the 'webshot' package and have phantomjs installed. With shinyapps.io that can be accomplished via 

```{r, eval=FALSE}
# Not evaluated here -- example of what to include in app for shinyapps.io
library(webshot)
install_phantomjs()
```

With earlier versions of iheatmapr, there was also a BioC dependency, so it was necessary to specify the BioC repos in options. This can be done via `setRepositories()`. Verifying what repositories are set can be done via `getOption("repos")`. This should no longer be needed with versions 1.0 and up of this package.

# Session Info

```{r}
sessionInfo()
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/groups.R
\name{add_col_groups}
\alias{add_col_groups}
\alias{add_col_groups,Iheatmap-method}
\title{add_col_groups}
\usage{
\S4method{add_col_groups}{Iheatmap}(
  p,
  groups,
  name = "Column<br>Groups",
  title = "Groups",
  colors = pick_discrete_colors(groups, p),
  colorbar_position = get_colorbar_position(p),
  show_colorbar = TRUE,
  show_title = TRUE,
  side = c("top", "bottom"),
  layout = list(),
  size = 0.05,
  buffer = 0.015,
  tooltip = setup_tooltip_options(),
  xname = current_xaxis(p),
  yname = NULL,
  pname = name
)
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object}

\item{groups}{vector of group names}

\item{name}{name of colorbar}

\item{title}{name of x axis label}

\item{colors}{palette name or vector of colors}

\item{colorbar_position}{colorbar placement}

\item{show_colorbar}{show the colorbar?}

\item{show_title}{show title as axis label}

\item{side}{side of plot on which to groups annotation}

\item{layout}{list of layout parameters for x axis}

\item{size}{relative size of dendrogram (relative to the main heatmap)}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{tooltip}{tooltip options, see \code{\link{setup_tooltip_options}}}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{pname}{internal name of plot}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds annotation to heatmap indicating what group every column of main heatmap
belongs to
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)
col_groups <- c("A","A","B","D","B")
hm <- iheatmap(mat) \%>\% add_col_groups(col_groups, name = "My Groups")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}, \code{\link{add_row_groups}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iheatmapr.R
\docType{data}
\name{measles}
\alias{measles}
\title{measles}
\description{
Data on measles cases for different states from 1930 to 2001
}
\examples{
data(measles)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{RowLabels-class}
\alias{RowLabels-class}
\alias{RowLabels}
\title{RowLabels}
\description{
Class for storing row labels
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{the names of the tick labels}

\item{\code{positions}}{the positions of the tick labels}

\item{\code{side}}{side of plot on which dendrogram is positioned, controls 
orientation}

\item{\code{textangle}}{angle for text}

\item{\code{font}}{list of font attributes}
}}

\seealso{
\code{\link{add_row_labels}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iheatmapr.R
\docType{package}
\name{iheatmapr}
\alias{iheatmapr}
\title{iheatmapr}
\description{
Interactive complex heatmaps in R
}
\details{
iheatmapr is a package for building complex, interactive heatmaps in R that 
can be explored in interactive R sessions or incorporated into rmarkdown 
documents, shiny applications, or standalone html files. 

The package includes a modular system for building up complex heatmaps, where
subplots get iteratively added to the top/left/right/bottom of the main 
heatmap(s). The \code{\link{iheatmap}} function provides a wrapper around
many of the common modular subcomponents to build fairly standard, moderately
complex heatmap.  

See the vignette for detailed instructions for how to use the package.

iheatmapr uses the plotly javascript library (\url{https://plotly.com/}) for making the 
interactive figures and htmlwidgets (http://www.htmlwidgets.org/) for 
rendering them in R.
}
\seealso{
\code{\link{main_heatmap}}, \code{\link{iheatmap}},
\code{\link{Iheatmap-class}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{ColumnAnnotation-class}
\alias{ColumnAnnotation-class}
\alias{ColumnAnnotation}
\title{ColumnAnnotation}
\description{
Class for storing row annotation
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{vector of annotation values}

\item{\code{colorbar}}{name of colorbar}

\item{\code{show_colorbar}}{show the colorbar?}
}}

\seealso{
\code{\link{add_col_annotation}}, \code{\link{add_col_signal}},
\code{\link{add_col_groups}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapMainX-class}
\alias{IheatmapMainX-class}
\alias{IheatmapMainX}
\title{IheatmapMainX}
\description{
Class for storing X axis information for "main" x axis-- x axis for 
\code{\link{MainHeatmap-class}}
}
\section{Slots}{

\describe{
\item{\code{id}}{plotly id for axis}

\item{\code{domain_start}}{start of domain (0 to 1)}

\item{\code{domain_end}}{end of domain (0 to 1)}

\item{\code{anchor}}{anchor for axis}

\item{\code{layout}}{plotly layout parameters}

\item{\code{categorical}}{is axis categorical?}

\item{\code{order}}{ordering of columns}

\item{\code{text}}{text labels for columns}
}}

\seealso{
\code{\link{Iheatmap-class}}, \code{\link{IheatmapAxis}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapAxis-class}
\alias{IheatmapAxis-class}
\alias{IheatmapAxis}
\title{IheatmapAxis}
\description{
Class for storing axis information
}
\section{Slots}{

\describe{
\item{\code{id}}{plotly id for axis}

\item{\code{domain_start}}{start of domain (0 to 1)}

\item{\code{domain_end}}{end of domain (0 to 1)}

\item{\code{anchor}}{anchor for axis}

\item{\code{layout}}{plotly layout parameters}
}}

\section{SubClasses}{

\itemize{
\item \code{\link{IheatmapX-class}}
\item \code{\link{IheatmapY-class}}
\item \code{\link{IheatmapMainX-class}}
\item \code{\link{IheatmapMainY-class}} 
}
}

\seealso{
\code{\link{Iheatmap-class}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/axis_labels.R
\name{add_row_labels}
\alias{add_row_labels}
\alias{add_row_labels,Iheatmap-method}
\title{add_row_labels}
\usage{
\S4method{add_row_labels}{Iheatmap}(
  p,
  tickvals = NULL,
  ticktext = NULL,
  textangle = 0,
  font = get_layout(p)$font,
  side = c("left", "right"),
  size = 0.1,
  buffer = 0.005,
  xname = NULL,
  yname = current_yaxis(p)
)
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object}

\item{tickvals}{row indices at which to place axis tick labels}

\item{ticktext}{text for axis tick labels}

\item{textangle}{angle for ticktext}

\item{font}{list of plotly font attributes, see 
\url{https://plotly.com/javascript/reference/#layout-font}}

\item{side}{side of plot on which to add subplot}

\item{size}{relative size of subplot relative to main heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal name for xaxis}

\item{yname}{internal name for yaxis}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to 
generate an interactive graphic
}
\description{
Add y axis labels to plot
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm1 <- iheatmap(mat) \%>\% add_row_labels()
hm2 <- iheatmap(mat) \%>\% add_row_labels(ticktext = letters[23:26])


# Print heatmaps if interactive session 
if (interactive()) hm1
if (interactive()) hm2 
}
\seealso{
\code{\link{add_row_title}}, \code{\link{iheatmap}}, 
\code{\link{add_col_labels}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/groups.R
\name{add_row_groups}
\alias{add_row_groups}
\alias{add_row_groups,Iheatmap-method}
\title{add_row_groups}
\usage{
\S4method{add_row_groups}{Iheatmap}(
  p,
  groups,
  name = "Row<br>Groups",
  title = "Groups",
  colors = pick_discrete_colors(groups, p),
  colorbar_position = get_colorbar_position(p),
  show_colorbar = TRUE,
  show_title = TRUE,
  side = c("right", "left"),
  layout = list(),
  size = 0.05,
  buffer = 0.005,
  tooltip = setup_tooltip_options(),
  xname = NULL,
  yname = current_yaxis(p),
  pname = name
)
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object}

\item{groups}{vector of group names}

\item{name}{name of colorbar}

\item{title}{name of x axis label}

\item{colors}{palette name or vector of colors}

\item{colorbar_position}{colorbar placement}

\item{show_colorbar}{show the colorbar?}

\item{show_title}{show title as axis label}

\item{side}{side of plot on which to groups annotation}

\item{layout}{list of layout parameters for x axis}

\item{size}{relative size of dendrogram (relative to the main heatmap)}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{tooltip}{tooltip options, see \code{\link{setup_tooltip_options}}}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{pname}{internal name of plot}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds annotation to heatmap indicating what group every row of main heatmap
belongs to
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)
row_groups <- c("A","A","B","D")
hm <- iheatmap(mat) \%>\% add_row_groups(row_groups, name = "My Groups")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}, \code{\link{add_col_groups}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapY-class}
\alias{IheatmapY-class}
\alias{IheatmapY}
\title{IheatmapY}
\description{
Class for storing Y axis information
}
\section{Slots}{

\describe{
\item{\code{id}}{plotly id for axis}

\item{\code{domain_start}}{start of domain (0 to 1)}

\item{\code{domain_end}}{end of domain (0 to 1)}

\item{\code{anchor}}{anchor for axis}

\item{\code{layout}}{plotly layout parameters}
}}

\seealso{
\code{\link{Iheatmap-class}}, \code{\link{IheatmapAxis}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusters.R
\name{add_col_clusters}
\alias{add_col_clusters}
\alias{add_col_clusters,Iheatmap-method}
\title{add_col_clusters}
\usage{
\S4method{add_col_clusters}{Iheatmap}(
  p,
  clusters,
  name = "Col<br>Clusters",
  reorder = TRUE,
  side = c("top", "bottom"),
  xname = current_xaxis(p),
  ...
)
}
\arguments{
\item{p}{iheatmap object}

\item{clusters}{cluster assignments, should be vector of integers, 
characters, or factors}

\item{name}{name of colorbar indicating cluster membership}

\item{reorder}{reorder rows based on clusters? default is TRUE}

\item{side}{side of plot on which to add subplot}

\item{xname}{name of xaxis}

\item{...}{additional arguments to pass to \code{\link{add_col_groups}} for 
creation of annotation heatmap indicating cluster membership}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Add column groups and order columns based on groups
}
\details{
This function is very similar to \code{\link{add_col_groups}}; the 
main difference is that with this function column will get reordered based on 
the groups.
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)
clusters <- c("A","B","A","B","A")

hm <- iheatmap(mat) \%>\% add_col_clusters(clusters)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_row_clusters}}, \code{\link{add_col_clustering}}, 
\code{\link{iheatmap}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_widget.R
\name{knit_print.Iheatmap}
\alias{knit_print.Iheatmap}
\title{knit_print.Iheatmap}
\usage{
\method{knit_print}{Iheatmap}(x, options)
}
\arguments{
\item{x}{Iheatmap object}

\item{options}{knitr options}
}
\description{
knit_print.Iheatmap
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{ColumnTitle-class}
\alias{ColumnTitle-class}
\alias{ColumnTitle}
\title{ColumnTitle}
\description{
Class for storing column title
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{the title (character)}

\item{\code{side}}{side of plot on which dendrogram is positioned, controls 
orientation}

\item{\code{textangle}}{angle for text}

\item{\code{font}}{list of font attributes}
}}

\seealso{
\code{\link{add_col_title}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapShapes-class}
\alias{IheatmapShapes-class}
\alias{IheatmapShapes}
\title{IheatmapShapes}
\description{
Class for storing \code{\link{IheatmapShapes}} objects
}
\seealso{
\code{\link{IheatmapShape-class}}, \code{\link{Iheatmap-class}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapX-class}
\alias{IheatmapX-class}
\alias{IheatmapX}
\title{IheatmapX}
\description{
Class for storing X axis information
}
\section{Slots}{

\describe{
\item{\code{id}}{plotly id for axis}

\item{\code{domain_start}}{start of domain (0 to 1)}

\item{\code{domain_end}}{end of domain (0 to 1)}

\item{\code{anchor}}{anchor for axis}

\item{\code{layout}}{plotly layout parameters}
}}

\seealso{
\code{\link{Iheatmap-class}}, \code{\link{IheatmapAxis}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{GenericPlot-class}
\alias{GenericPlot-class}
\alias{GenericPlot}
\title{GenericPlot}
\description{
Class for storing an arbitrary subplot
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{list of plotly parameters}
}}

\seealso{
\code{\link{add_subplot}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapAnnotations-class}
\alias{IheatmapAnnotations-class}
\alias{IheatmapAnnotations}
\title{IheatmapAnnotations}
\description{
Class for storing \code{\link{IheatmapAnnotation}} objects
}
\seealso{
\code{\link{IheatmapAnnotation-class}}, \code{\link{Iheatmap-class}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapAxes-class}
\alias{IheatmapAxes-class}
\alias{IheatmapAxes}
\title{IheatmapAxes}
\description{
Class for storing \code{\link{IheatmapAxis}} objects
}
\section{Slots}{

\describe{
\item{\code{axis}}{x or y?}
}}

\seealso{
\code{\link{IheatmapAxis-class}}, \code{\link{Iheatmap-class}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/axis_titles.R
\name{add_row_title}
\alias{add_row_title}
\alias{add_row_title,Iheatmap-method}
\title{add_row_title}
\usage{
\S4method{add_row_title}{Iheatmap}(
  p,
  title,
  textangle = ifelse(side == "left", -90, 90),
  font = get_layout(p)$font,
  side = c("left", "right"),
  size = 0.1,
  buffer = 0.01,
  xname = NULL,
  yname = current_yaxis(p)
)
}
\arguments{
\item{p}{iheatmap object}

\item{title}{title of axis}

\item{textangle}{angle of text}

\item{font}{list of plotly font attributes, see 
\url{https://plotly.com/javascript/reference/#layout-font}}

\item{side}{side of plot on which to add subplot}

\item{size}{relative size of subplot relative to main heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal name for xaxis}

\item{yname}{internal name for yaxis}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Add y axis title to plot
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm <- iheatmap(mat) \%>\% add_row_title("Samples")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_col_title}}, \code{\link{iheatmap}}, 
\code{\link{add_row_labels}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subplot.R
\name{add_subplot}
\alias{add_subplot}
\alias{add_subplot,Iheatmap-method}
\title{add_subplot}
\usage{
\S4method{add_subplot}{Iheatmap}(
  p,
  ...,
  side = c("top", "bottom", "right", "left"),
  layout = list(),
  size = 1,
  buffer = 0.1,
  xname = if (side \%in\% c("top", "bottom")) current_xaxis(p) else NULL,
  yname = if (side \%in\% c("left", "right")) current_yaxis(p) else NULL,
  pname = "subplot"
)
}
\arguments{
\item{p}{iheatmap object}

\item{...}{arguments to pass to plotly trace, see plotly.js documentation at
\url{https://plotly.com/javascript/reference/}}

\item{side}{which side of the current plot to add this heatmap? "right", 
"left","top", or "bottom"}

\item{layout}{axis layout parameters (list)}

\item{size}{relative size of plot.  size relative to first heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{pname}{internal name of plot}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds an arbitrary subplot to iheatmap
}
\examples{
mat <- matrix(rnorm(24), ncol = 6)
hm <- iheatmap(mat) \%>\% add_subplot(x = 1:5, y=1:5, side = "top")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{RowAnnotation-class}
\alias{RowAnnotation-class}
\alias{RowAnnotation}
\title{RowAnnotation}
\description{
Class for storing row annotation
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{vector of annotation values}

\item{\code{colorbar}}{name of colorbar}

\item{\code{show_colorbar}}{show the colorbar?}
}}

\seealso{
\code{\link{add_row_annotation}}, \code{\link{add_row_signal}},
\code{\link{add_row_groups}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_heatmap.R
\name{setup_tooltip_options}
\alias{setup_tooltip_options}
\title{Tooltip Options}
\usage{
setup_tooltip_options(
  row = TRUE,
  col = TRUE,
  value = TRUE,
  prepend_row = "Row: ",
  prepend_col = "Col: ",
  prepend_value = "Value: "
)
}
\arguments{
\item{row}{logical, include row name in tooltip?}

\item{col}{logical, include column name in tooltip?}

\item{value}{logical, include value in tooltip?}

\item{prepend_row}{text to prepend to row name}

\item{prepend_col}{text to prepend to column name}

\item{prepend_value}{text to prepend to value}
}
\value{
a HeatmapTooltipOptions object which stores these options and can be
passed to 'tooltip' argument to main_heatmap and other functions.
}
\description{
This function setups tooltip options for heatmap components of iheatmapr 
complex heatmaps.
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm1 <- main_heatmap(mat, 
   tooltip = setup_tooltip_options(row = FALSE, col = FALSE,
                                   prepend_value = "Value is ")) 

# Print heatmap if interactive session 
if (interactive()) hm1 
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapMainY-class}
\alias{IheatmapMainY-class}
\alias{IheatmapMainY}
\title{IheatmapMainY}
\description{
Class for storing Y axis information for "main" y axis-- y axis for 
\code{\link{MainHeatmap-class}}
}
\section{Slots}{

\describe{
\item{\code{id}}{plotly id for axis}

\item{\code{domain_start}}{start of domain (0 to 1)}

\item{\code{domain_end}}{end of domain (0 to 1)}

\item{\code{anchor}}{anchor for axis}

\item{\code{layout}}{plotly layout parameters}

\item{\code{categorical}}{is axis categorical?}

\item{\code{order}}{ordering of rows}

\item{\code{text}}{text labels for rows}
}}

\seealso{
\code{\link{Iheatmap-class}}, \code{\link{IheatmapAxis}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic_subplots.R
\name{add_col_plot}
\alias{add_col_plot}
\alias{add_col_plot,Iheatmap-method}
\title{add_col_plot}
\usage{
\S4method{add_col_plot}{Iheatmap}(
  p,
  y,
  ...,
  mode = c("lines+markers", "lines", "markers"),
  color = NULL,
  tracename = NA_character_,
  showlegend = !is.na(tracename),
  side = c("top", "bottom"),
  layout = list(),
  size = 0.2,
  buffer = 0.02,
  xname = current_xaxis(p),
  yname = NULL,
  pname = if (!is.na(tracename)) tracename else "col_plot"
)
}
\arguments{
\item{p}{iheatmap object}

\item{y}{y axis values}

\item{...}{additional arguments to add to plotly scatter trace, see
\url{https://plotly.com/javascript/reference/#scatter}}

\item{mode}{mode of plot -- one of "lines+markers","lines", or "markers"}

\item{color}{color of bars}

\item{tracename}{name of trace (for legend and hover)}

\item{showlegend}{show in legend?}

\item{side}{side of plot on which to add subplot}

\item{layout}{yaxis layout list}

\item{size}{relative size of subplot relative to main heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{pname}{internal name of plot}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Add a scatter or line plot with one point per column of the main heatmap
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm <- iheatmap(mat) \%>\% add_col_plot(y = 1:5, tracename = "Strength")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_col_signal}}, \code{\link{iheatmap}}, 
\code{\link{add_col_barplot}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/axis_labels.R
\name{add_col_labels}
\alias{add_col_labels}
\alias{add_col_labels,Iheatmap-method}
\title{add_col_labels}
\usage{
\S4method{add_col_labels}{Iheatmap}(
  p,
  tickvals = NULL,
  ticktext = NULL,
  textangle = -90,
  font = get_layout(p)$font,
  side = c("bottom", "top"),
  size = 0.1,
  buffer = 0.005,
  xname = current_xaxis(p),
  yname = NULL
)
}
\arguments{
\item{p}{\code{link{Iheatmap-class}} object}

\item{tickvals}{column indices at which to place axis tick labels}

\item{ticktext}{text for axis tick labels}

\item{textangle}{angle for ticktext}

\item{font}{list of plotly font attributes, see 
\url{https://plotly.com/javascript/reference/#layout-font}}

\item{side}{side of plot on which to add subplot}

\item{size}{relative size of subplot relative to main heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{name for xaxis}

\item{yname}{name for yaxis}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Add x axis labels to plot
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm1 <- iheatmap(mat) \%>\% add_col_labels()
hm2 <- iheatmap(mat) \%>\% add_col_labels(ticktext = letters[22:26])

# Print heatmap if interactive session 
if (interactive()) hm1
if (interactive()) hm2
}
\seealso{
\code{\link{add_row_title}}, \code{\link{iheatmap}}, 
\code{\link{add_col_labels}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colorbars.R
\name{setup_colorbar_grid}
\alias{setup_colorbar_grid}
\title{setup_colorbar_grid}
\usage{
setup_colorbar_grid(
  nrows = 3,
  y_length = y_spacing * 0.9,
  x_spacing = 0.16,
  y_spacing = y_start/nrows,
  x_start = 1.05,
  y_start = 0.9
)
}
\arguments{
\item{nrows}{number of rows in colorbar grid}

\item{y_length}{length of colorbar}

\item{x_spacing}{spacing along horizonatal axis between colorbars}

\item{y_spacing}{spacing along vertical axis between colorbars}

\item{x_start}{left most position of colorbar grid}

\item{y_start}{top most position of colorbar grid}
}
\value{
\code{\link{IheatmapColorbarGrid-class}} object
}
\description{
function to set parameters controlling colorbar placement in Iheatmap object
}
\examples{

cb_grid <- setup_colorbar_grid(nrows = 2, x_spacing = 0.2)
mat <- matrix(rnorm(24), nrow = 6)
hm <- iheatmap(mat, colorbar_grid = cb_grid, cluster_rows = "kmeans",
         cluster_cols = "kmeans", row_k = 3, col_k = 2)

# Print heatmap if interactive session 
if (interactive()) hm 
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{ContinuousColorbar-class}
\alias{ContinuousColorbar-class}
\alias{ContinuousColorbar}
\title{ContinuousColorbar}
\description{
Class for storing continuous colorbar information
}
\section{Slots}{

\describe{
\item{\code{title}}{title for colorbar}

\item{\code{position}}{integer indicating relative position of colorbar}

\item{\code{colors}}{name of color palette or vector of colors}

\item{\code{zmid}}{midpoint of colorbar}

\item{\code{zmin}}{min of colorbar}

\item{\code{zmax}}{max of colorbar}
}}

\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapPlot-class}
\alias{IheatmapPlot-class}
\alias{IheatmapPlot}
\title{IheatmapPlot}
\description{
Virtual class for storing plot objects
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{main plot data}
}}

\section{SubClasses}{

\itemize{
\item \code{\link{MainHeatmap-class}}
\item \code{\link{RowAnnotation-class}}
\item \code{\link{ColumnAnnotation-class}}
\item \code{\link{RowPlot-class}}
\item \code{\link{ColumnPlot-class}}
\item \code{\link{GenericPlot-class}}
}
}

\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signal.R
\name{add_row_signal}
\alias{add_row_signal}
\alias{add_row_signal,Iheatmap-method}
\title{add_row_signal}
\usage{
\S4method{add_row_signal}{Iheatmap}(
  p,
  signal,
  name,
  title = name,
  xname = NULL,
  yname = current_yaxis(p),
  pname = name,
  colorbar_position = get_colorbar_position(p),
  colors = pick_continuous_colors(zmid, zmin, zmax, p = p),
  zmin = min(signal, na.rm = TRUE),
  zmax = max(signal, na.rm = TRUE),
  zmid = 0,
  side = c("right", "left"),
  size = 0.05,
  buffer = 0.015,
  text = signif(signal, digits = 3),
  tooltip = setup_tooltip_options(),
  show_colorbar = TRUE,
  show_title = TRUE,
  layout = list()
)
}
\arguments{
\item{p}{iheatmap object}

\item{signal}{vector of signal}

\item{name}{name of colorbar}

\item{title}{label for x axis}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{pname}{internal name of plot}

\item{colorbar_position}{colorbar placement}

\item{colors}{color palette or vector of colors}

\item{zmin}{minimum for colorscale}

\item{zmax}{maximum for colorscale}

\item{zmid}{midpoint for colorscale}

\item{side}{side of plot on which to add dendro}

\item{size}{relative size of dendrogram (relative to the main heatmap)}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{text}{text of value to display for data}

\item{tooltip}{tooltip options, see \code{\link{setup_tooltip_options}}}

\item{show_colorbar}{show the colorbar?}

\item{show_title}{show title as axis label}

\item{layout}{list of x axis layout parameters}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds single column heatmap to iheatmap object
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm <- iheatmap(mat) \%>\% add_row_signal(signal = 1:4, name = "Strength")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}, \code{\link{add_col_groups}}

\code{\link{add_col_signal}}, \code{\link{iheatmap}}, 
\code{\link{add_row_annotation}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{add_col_clustering}
\alias{add_col_clustering}
\alias{add_col_clustering,Iheatmap-method}
\title{add_col_clustering}
\usage{
\S4method{add_col_clustering}{Iheatmap}(
  p,
  method = c("hclust", "kmeans", "groups"),
  name = "Col<br>Clusters",
  k = NULL,
  groups = NULL,
  clust_dist = stats::dist,
  colors = NULL,
  show_colorbar = TRUE,
  side = c("top", "bottom"),
  yname = NULL,
  xname = current_xaxis(p)
)
}
\arguments{
\item{p}{iheatmap object}

\item{method}{"hclust" or "kmeans" for hierarchical or k-means clustering, 
respectively}

\item{name}{name of colorbar indicating cluster membership}

\item{k}{number of clusters for rows, needed if order is kmeans or optional 
if hclust}

\item{groups}{vector of group assignments}

\item{clust_dist}{distance function to use for clustering if hierarchical 
clustering}

\item{colors}{colors to use for annotation of grouping, can be RColorBrewer 
palette name or
vector of colors}

\item{show_colorbar}{show the colorbar for the heatmap indicating cluster 
membership}

\item{side}{side of plot on which to add subplot}

\item{yname}{name of yaxis}

\item{xname}{name of xaxis}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
add_col_clustering
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)
hm <- iheatmap(mat) \%>\% add_col_clustering(method = "hclust", k = 2)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_row_clustering}}, \code{\link{iheatmap}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iheatmap.R
\name{add_iheatmap}
\alias{add_iheatmap}
\alias{add_iheatmap,IheatmapHorizontal,matrix-method}
\alias{add_iheatmap,IheatmapVertical,matrix-method}
\title{add_iheatmap}
\usage{
\S4method{add_iheatmap}{IheatmapHorizontal,matrix}(
  p,
  data,
  x = default_x(data),
  cluster_cols = c("none", "hclust", "kmeans", "groups"),
  col_clusters = NULL,
  col_k = NULL,
  col_clust_dist = stats::dist,
  name = "Signal",
  scale = c("none", "rows", "cols"),
  scale_method = c("standardize", "center", "normalize"),
  colors = NULL,
  col_clusters_colors = NULL,
  col_clusters_name = "Col<br>Clusters",
  show_col_clusters_colorbar = TRUE,
  row_annotation = NULL,
  col_annotation = NULL,
  row_annotation_colors = NULL,
  col_annotation_colors = NULL,
  row_labels = NULL,
  col_labels = NULL,
  row_title = NULL,
  col_title = NULL,
  buffer = 0.2,
  ...
)

\S4method{add_iheatmap}{IheatmapVertical,matrix}(
  p,
  data,
  y = default_y(data),
  cluster_rows = c("none", "hclust", "kmeans", "groups"),
  row_clusters = NULL,
  row_k = NULL,
  row_clust_dist = stats::dist,
  name = "Signal",
  scale = c("none", "rows", "cols"),
  scale_method = c("standardize", "center", "normalize"),
  colors = NULL,
  row_clusters_colors = NULL,
  row_clusters_name = "Col<br>Clusters",
  show_row_clusters_colorbar = TRUE,
  row_annotation = NULL,
  col_annotation = NULL,
  row_annotation_colors = NULL,
  col_annotation_colors = NULL,
  row_labels = NULL,
  col_labels = NULL,
  row_title = NULL,
  col_title = NULL,
  buffer = 0.2,
  ...
)
}
\arguments{
\item{p}{iheatmap object}

\item{data}{matrix of values to be plotted as heatmap}

\item{x}{x xaxis labels, by default colnames of data}

\item{cluster_cols}{"none","hclust", or "k-means" for no clustering, 
hierarchical clustering, and k-means clustering of columnsrespectively}

\item{col_clusters}{vector of pre-determined column cluster assignment}

\item{col_k}{number of clusters for columns, needed if cluster_rows is kmeans 
or optional if hclust}

\item{col_clust_dist}{distance function to use for column clustering if 
hierarchical clustering}

\item{name}{Name for colorbar}

\item{scale}{scale matrix by rows, cols or none}

\item{scale_method}{what method to use for scaling, either standardize, 
center, normalize}

\item{colors}{name of RColorBrewer palette or vector of colors for main heatmap}

\item{col_clusters_colors}{colors for col clusters annotation heatmap}

\item{col_clusters_name}{name for col clusters colorbar}

\item{show_col_clusters_colorbar}{show the colorbar for column clusters?}

\item{row_annotation}{row annotation data.frame}

\item{col_annotation}{column annotation data.frame}

\item{row_annotation_colors}{list of colors for row annotations heatmap}

\item{col_annotation_colors}{list of colors for col annotations heatmap}

\item{row_labels}{axis labels for y axis}

\item{col_labels}{axis labels for x axis}

\item{row_title}{x axis title}

\item{col_title}{y axis title}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{...}{additional argument to add_iheatmap}

\item{y}{y axis labels, by default rownames of data}

\item{cluster_rows}{"none","hclust", or "k-means" for no clustering, 
hierarchical clustering, and k-means clustering of rows respectively}

\item{row_clusters}{vector of pre-determined row cluster assignment}

\item{row_k}{number of clusters for rows, needed if cluster_rows is kmeans or 
optional if hclust}

\item{row_clust_dist}{distance function to use for row clustering if 
hierarchical clustering}

\item{row_clusters_colors}{colors for row clusters annotation heatmap}

\item{row_clusters_name}{name for row clusters colorbar}

\item{show_row_clusters_colorbar}{show the colorbar for row clusters?}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
add_iheatmap
}
\details{
By default, no scaling is done of rows or columns. This can be changed by 
specifying the 'scale' argument.  There are three options for scaling 
methods. "standardize" subtracts the mean and divides by standard deviation,
"center" just subtracts the mean, and "normalize" divides by the sum of the 
values.  "normalize" should only be used for data that is all positive!  If 
alternative scaling is desired, the scaling should be done prior to calling 
the iheatmap function.
}
\examples{

mat <- matrix(rnorm(24), nrow = 6)
mat2 <- matrix(rnorm(24), nrow = 6)
annotation = data.frame(gender = c(rep("M", 3),rep("F",3)),
                        age = c(20,34,27,19,23,30))
hm <- iheatmap(mat, 
 cluster_rows = "hclust", 
 cluster_cols = "hclust", 
 col_k = 3) \%>\%
add_iheatmap(mat2, 
 cluster_cols = "hclust", 
 col_k = 3, 
 row_annotation = annotation)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}, \code{\link{main_heatmap}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{modify_layout}
\alias{modify_layout}
\alias{modify_layout,Iheatmap-method}
\title{modify_layout}
\usage{
\S4method{modify_layout}{Iheatmap}(x, new_layout)
}
\arguments{
\item{x}{Iheatmap}

\item{new_layout}{list of new layout parameter}
}
\value{
modified Iheatmap object
}
\description{
modify_layout
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm <- main_heatmap(mat) \%>\% modify_layout(list(margin = list(b = 120))) 

# Print heatmap if interactive session 
if (interactive()) hm 
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{ColumnPlot-class}
\alias{ColumnPlot-class}
\alias{ColumnPlot}
\title{ColumnPlot}
\description{
Class for storing generic column plot
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{vector of values}

\item{\code{name}}{tracename}

\item{\code{type}}{trace type}

\item{\code{showlegend}}{show the legend?}

\item{\code{additional}}{additional plotly parameters}
}}

\seealso{
\code{\link{add_col_plot}}, \code{\link{add_col_barplot}},
\code{\link{add_col_summary}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapColorbars-class}
\alias{IheatmapColorbars-class}
\alias{IheatmapColorbars}
\title{IheatmapColorbars}
\description{
Class for storing \code{\link{IheatmapColorbar}} objects
}
\seealso{
\code{\link{IheatmapColorbar-class}}, \code{\link{Iheatmap-class}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\docType{methods}
\name{add_component}
\alias{add_component}
\alias{add_annotation,Iheatmap,IheatmapAnnotation-method}
\alias{add_axis,IheatmapHorizontal,IheatmapX-method}
\alias{add_axis,IheatmapHorizontal,IheatmapY-method}
\alias{add_axis,IheatmapVertical,IheatmapX-method}
\alias{add_axis,IheatmapVertical,IheatmapY-method}
\alias{add_colorbar,Iheatmap,ContinuousColorbar-method}
\alias{add_colorbar,Iheatmap,DiscreteColorbar-method}
\alias{add_plot,Iheatmap,IheatmapPlot-method}
\alias{add_shape,Iheatmap,IheatmapShape-method}
\alias{add_axis}
\alias{add_colorbar}
\alias{add_plot}
\alias{add_shape}
\alias{add_annotation}
\title{Adding plot components to iheatmapr}
\usage{
add_axis(p, new_axis, ...)

add_colorbar(p, new_colorbar, ...)

add_plot(p, new_plot, ...)

add_shape(p, new_shape, ...)

add_annotation(p, new_anno, ...)
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object}

\item{new_axis}{new \code{\link{IheatmapAxis-class}} object}

\item{new_colorbar}{new \code{\link{IheatmapColorbar-class}} object}

\item{new_plot}{new \code{\link{IheatmapPlot-class}} object}

\item{new_shape}{new \code{\link{IheatmapShape-class}} object}

\item{new_anno}{new \code{\link{IheatmapAnnotation-class}} object}

\item{name}{internal name}
}
\description{
These are generic methods for adding new plot components to an 
\code{link{Iheatmap-class}} object. Not intended for end users; exported for
developers seeking to create new Iheatmap subplots.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colorbars.R, R/components.R
\docType{methods}
\name{colorbars,Iheatmap-method}
\alias{colorbars,Iheatmap-method}
\alias{access_component}
\alias{yaxes,Iheatmap-method}
\alias{xaxes,Iheatmap-method}
\alias{plots,Iheatmap-method}
\alias{shapes,Iheatmap-method,}
\alias{annotations,Iheatmap-method}
\alias{yaxes}
\alias{xaxes}
\alias{plots}
\alias{shapes}
\alias{annotations}
\alias{colorbars}
\alias{colorbars,IheatmapColorbars-method}
\alias{colorbars,IheatmapPlots-method}
\alias{shapes,Iheatmap-method}
\title{Access subcomponents of Iheatmap object}
\usage{
\S4method{colorbars}{Iheatmap}(x, what = c("all", "continuous", "discrete"))

\S4method{yaxes}{Iheatmap}(p, xaxis = NULL)

\S4method{xaxes}{Iheatmap}(p, yaxis = NULL)

\S4method{plots}{Iheatmap}(x)

\S4method{shapes}{Iheatmap}(x)

\S4method{annotations}{Iheatmap}(x)
}
\arguments{
\item{x}{\code{\link{Iheatmap-class}} object}

\item{xaxis}{name of xaxis}

\item{yaxis}{name of yaxis}
}
\description{
These are methods for accessing subcomponents of the Iheatmap object
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{MainHeatmap-class}
\alias{MainHeatmap-class}
\alias{MainHeatmap}
\title{MainHeatmap}
\description{
Class for storing main heatmap
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{matrix of heatmap values}

\item{\code{colorbar}}{name of colorbar}

\item{\code{show_colorbar}}{show the colorbar?}
}}

\seealso{
\code{\link{main_heatmap}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/axes.R
\name{reorder_cols}
\alias{reorder_cols}
\alias{reorder_cols,IheatmapHorizontal,integer-method}
\alias{reorder_cols,IheatmapVertical,integer-method}
\title{reorder_cols}
\usage{
\S4method{reorder_cols}{IheatmapHorizontal,integer}(p, col_order, xname = current_xaxis(p))

\S4method{reorder_cols}{IheatmapVertical,integer}(p, col_order)
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object}

\item{col_order}{integer vector}

\item{xname}{name of xaxis to reorder, only applicable if object is oriented
horizontally}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Reorder the columns of an \code{\link{Iheatmap-class}} object
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
dend <- hclust(dist(t(mat)))
hm <- iheatmap(mat) \%>\% reorder_cols(dend$order)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_row_clustering}}, \code{\link{reorder_cols}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barplot.R
\name{add_col_barplot}
\alias{add_col_barplot}
\alias{add_col_barplot,Iheatmap-method}
\title{add_col_barplot}
\usage{
\S4method{add_col_barplot}{Iheatmap}(
  p,
  y,
  ...,
  color = NULL,
  tracename = NA_character_,
  showlegend = !is.na(tracename),
  side = c("top", "bottom"),
  layout = list(),
  size = 0.2,
  buffer = 0.02,
  xname = current_xaxis(p),
  yname = NULL,
  pname = if (!is.na(tracename)) tracename else "col_barplot"
)
}
\arguments{
\item{p}{iheatmap object}

\item{y}{y axis values}

\item{...}{additional arguments to add to plotly scatter trace, see
\url{https://plotly.com/javascript/reference/#scatter}}

\item{color}{color of bars}

\item{tracename}{name of trace (for legend and hover)}

\item{showlegend}{show in legend?}

\item{side}{side of plot on which to add subplot}

\item{layout}{yaxis layout list}

\item{size}{relative size of subplot relative to main heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{pname}{internal name of plot}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Add bar plot with one bar per column above or below a main heatmap
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm <- iheatmap(mat) \%>\% add_col_barplot(y = 1:5, tracename = "Strength")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_col_signal}}, \code{\link{iheatmap}}, 
\code{\link{add_col_plot}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{RowPlot-class}
\alias{RowPlot-class}
\alias{RowPlot}
\title{RowPlot}
\description{
Class for storing generic row plot
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{vector of values}

\item{\code{name}}{tracename}

\item{\code{type}}{trace type}

\item{\code{showlegend}}{show the legend?}

\item{\code{additional}}{additional plotly parameters}
}}

\seealso{
\code{\link{add_row_plot}}, \code{\link{add_row_barplot}},
\code{\link{add_row_summary}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/axis_titles.R
\name{add_col_title}
\alias{add_col_title}
\alias{add_col_title,Iheatmap-method}
\title{add_col_title}
\usage{
\S4method{add_col_title}{Iheatmap}(
  p,
  title,
  textangle = 0,
  font = get_layout(p)$font,
  side = c("bottom", "top"),
  size = 0.1,
  buffer = 0.01,
  xname = current_xaxis(p),
  yname = NULL
)
}
\arguments{
\item{p}{iheatmap object}

\item{title}{title of axis}

\item{textangle}{angle of text}

\item{font}{list of plotly font attributes, see 
\url{https://plotly.com/javascript/reference/#layout-font}}

\item{side}{side of plot on which to add subplot}

\item{size}{relative size of subplot relative to main heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{name for xaxis}

\item{yname}{name for yaxis}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Add x axis title to plot
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm <- iheatmap(mat) \%>\% add_col_title("My x-axis")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_col_labels}}, \code{\link{iheatmap}}, 
\code{\link{add_row_title}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_heatmap.R
\name{main_heatmap}
\alias{main_heatmap}
\alias{main_heatmap,matrix-method}
\title{main_heatmap}
\usage{
\S4method{main_heatmap}{matrix}(
  data,
  name = "Signal",
  x = default_x(data),
  y = default_y(data),
  colors = pick_continuous_colors(zmid, zmin, zmax),
  colorbar_grid = setup_colorbar_grid(),
  colorbar_position = 1,
  zmid = 0,
  zmin = min(data, na.rm = TRUE),
  zmax = max(data, na.rm = TRUE),
  orientation = c("horizontal", "vertical"),
  x_categorical = NULL,
  y_categorical = NULL,
  row_order = seq_len(nrow(data)),
  col_order = seq_len(ncol(data)),
  text = signif(data, digits = 3),
  tooltip = setup_tooltip_options(),
  xname = "x",
  yname = "y",
  pname = name,
  source = "iheatmapr",
  show_colorbar = TRUE,
  layout = list()
)
}
\arguments{
\item{data}{matrix}

\item{name}{name of colorbar}

\item{x}{x axis labels (by default rownames of data)}

\item{y}{y axis labels (by default colnames of data)}

\item{colors}{color palette or vector of colors}

\item{colorbar_grid}{colorbar grid parameters, should be result from 
\code{\link{setup_colorbar_grid}}}

\item{colorbar_position}{colorbar placement, should be positive integer}

\item{zmid}{midpoint for colorscale}

\item{zmin}{minimum for colorscale}

\item{zmax}{maximum for colorscale}

\item{orientation}{should new main plots be added horizontally or vertically?}

\item{x_categorical}{is x categorical?  will guess if not provided}

\item{y_categorical}{is y categorical?  will guess if not provided}

\item{row_order}{row ordering for this heatmap-- will be used for all 
subsequent elements sharing y axis}

\item{col_order}{column ordering for this heatmap-- will be used for all 
subsequent elements sharing x axis}

\item{text}{text of value to display for data}

\item{tooltip}{tooltip options, see \code{\link{setup_tooltip_options}}}

\item{xname}{internal name for xaxis}

\item{yname}{internal name for yaxis}

\item{pname}{internal plot name}

\item{source}{source name for use with shiny}

\item{show_colorbar}{logical to indicate whether to show colorbar}

\item{layout}{list of layout attributes to pass to plotly, 
eg. list(font = list(size = 15))}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Plots initial heatmap, creates Iheatmap object
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm <- main_heatmap(mat) 

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_iheatmap}}, \code{\link{to_widget}},  
\code{\link{iheatmap}}, \code{\link{Iheatmap-class}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{DiscreteColorbar-class}
\alias{DiscreteColorbar-class}
\alias{DiscreteColorbar}
\title{DiscreteColorbar}
\description{
Class for storing discrete colorbar information
}
\section{Slots}{

\describe{
\item{\code{title}}{title for colorbar}

\item{\code{position}}{integer indicating relative position of colorbar}

\item{\code{colors}}{name of color palette or vector of colors}

\item{\code{ticktext}}{labels for categories}

\item{\code{tickvals}}{integer values for categories}
}}

\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_test.R
\name{test_iheatmapr_event}
\alias{test_iheatmapr_event}
\title{test_iheatmapr_event}
\usage{
test_iheatmapr_event(ihm, event = c("click", "hover", "relayout"))
}
\arguments{
\item{ihm}{Iheatmap object}

\item{event}{name of event, either "click","hover", or "relayout"}
}
\value{
shiny app
}
\description{
test_iheatmapr_event
}
\examples{

\dontrun{
  mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
  hm <- main_heatmap(mat) 
  test_iheatmapr_event(hm, "click")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{RowTitle-class}
\alias{RowTitle-class}
\alias{RowTitle}
\title{RowTitle}
\description{
Class for storing row title
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{the title (character)}

\item{\code{side}}{side of plot on which dendrogram is positioned, controls 
orientation}

\item{\code{textangle}}{angle for text}

\item{\code{font}}{list of font attributes}
}}

\seealso{
\code{\link{add_row_title}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signal.R
\name{add_col_signal}
\alias{add_col_signal}
\alias{add_col_signal,Iheatmap-method}
\title{add_col_signal}
\usage{
\S4method{add_col_signal}{Iheatmap}(
  p,
  signal,
  name,
  title = name,
  yname = NULL,
  xname = current_xaxis(p),
  pname = name,
  colorbar_position = get_colorbar_position(p),
  colors = pick_continuous_colors(zmid, zmin, zmax, p = p),
  zmin = min(signal, na.rm = TRUE),
  zmax = max(signal, na.rm = TRUE),
  zmid = 0,
  side = c("top", "bottom"),
  size = 0.05,
  buffer = 0.015,
  text = signif(signal, digits = 3),
  tooltip = setup_tooltip_options(),
  show_colorbar = TRUE,
  show_title = TRUE,
  layout = list()
)
}
\arguments{
\item{p}{iheatmap object}

\item{signal}{vector of signal}

\item{name}{name of colorbar}

\item{title}{label for y axis}

\item{yname}{internal name of yaxis}

\item{xname}{internal name of xaxis}

\item{pname}{internal name of plot}

\item{colorbar_position}{colorbar placement}

\item{colors}{palette or vector of colors to use}

\item{zmin}{minimum for colorscale}

\item{zmax}{maximum for colorscale}

\item{zmid}{midpoint for colorscale}

\item{side}{side of plot on which to add groups}

\item{size}{relative size of dendrogram (relative to the main heatmap)}

\item{buffer}{amount of space to leave empty before this plot, relative to size 
of first heatmap}

\item{text}{text of value to display for data}

\item{tooltip}{tooltip options, see \code{\link{setup_tooltip_options}}}

\item{show_colorbar}{show the colorbar?}

\item{show_title}{show title as axis label}

\item{layout}{y axis layout parameters to use}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds column signal to iheatmap object
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm <- iheatmap(mat) \%>\% add_col_signal(signal = 1:5, name = "Strength")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}, \code{\link{add_row_groups}}

\code{\link{add_row_signal}}, \code{\link{iheatmap}}, 
\code{\link{add_col_annotation}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotations.R
\name{add_row_annotation}
\alias{add_row_annotation}
\alias{add_row_annotation,Iheatmap-method}
\title{add_row_annotation}
\usage{
\S4method{add_row_annotation}{Iheatmap}(
  p,
  annotation,
  colors = NULL,
  side = c("right", "left"),
  size = 0.05,
  buffer = 0.015,
  inner_buffer = buffer/2,
  layout = list(),
  show_colorbar = TRUE
)
}
\arguments{
\item{p}{\code{link{Iheatmap-class}} object}

\item{annotation}{data.frame or object that can be converted to data frame}

\item{colors}{list of color palettes, with one color per annotation column 
name}

\item{side}{side of plot on which to add row annotation}

\item{size}{relative size of each row annotation}

\item{buffer}{relative size of buffer between previous subplot and row 
annotation}

\item{inner_buffer}{relative size of buffer between each annotation}

\item{layout}{layout properties for new x axis}

\item{show_colorbar}{logical indicator to show or hide colorbar}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds annotation heatmaps for one or more qualitative or quantitative 
annotations for each row of a main heatmap.
}
\examples{

mat <- matrix(rnorm(24), nrow = 6)
annotation <- data.frame(gender = c(rep("M", 3),rep("F",3)),
                        age = c(20,34,27,19,23,30))
hm <- iheatmap(mat) \%>\% add_row_annotation(annotation)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}, \code{\link{add_row_annotation}}, 
\code{\link{add_col_signal}}, \code{\link{add_col_groups}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/axes.R
\name{reorder_rows}
\alias{reorder_rows}
\alias{reorder_rows,IheatmapHorizontal,integer-method}
\alias{reorder_rows,IheatmapVertical,integer-method}
\title{reorder_rows}
\usage{
\S4method{reorder_rows}{IheatmapHorizontal,integer}(p, row_order)

\S4method{reorder_rows}{IheatmapVertical,integer}(p, row_order, yname = current_yaxis(p))
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object}

\item{row_order}{integer vector}

\item{yname}{name of yaxis to reorder, only applicable if object is oriented
vertically}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Reorder the rows of an \code{\link{Iheatmap-class}} object
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
dend <- hclust(dist(mat))
hm <- iheatmap(mat) \%>\% reorder_rows(dend$order)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_row_clustering}}, \code{\link{reorder_cols}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapPlots-class}
\alias{IheatmapPlots-class}
\alias{IheatmapPlots}
\title{IheatmapPlots}
\description{
Class for storing \code{\link{IheatmapPlot}} objects
}
\seealso{
\code{\link{IheatmapPlot-class}}, \code{\link{Iheatmap-class}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_widget.R
\name{to_plotly}
\alias{to_plotly}
\alias{to_plotly_list}
\alias{to_plotly_json}
\title{Convert Iheatmap to plotly spec}
\usage{
to_plotly_list(p)

to_plotly_json(p)
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object to convert}
}
\value{
Returns a JSON for a plotly spec for to_plotly_spec and
as a list of same plotly object for to_plotly_list.
}
\description{
Function  to convert \code{link{Iheatmap-class}} object to a plotly spec 
either as a list or json
}
\examples{

mat <- matrix(rnorm(24), nrow = 6)
hm_json <- iheatmap(mat) \%>\% to_plotly_json()
hm_list <- iheatmap(mat) \%>\% to_plotly_list()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapShape-class}
\alias{IheatmapShape-class}
\alias{IheatmapShape}
\title{IheatmapShape}
\description{
Virtual class for storing shape objects
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{main shape data}
}}

\section{SubClasses}{

\itemize{
\item \code{\link{Dendrogram-class}}
}
}

\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic_subplots.R
\name{add_row_plot}
\alias{add_row_plot}
\alias{add_row_plot,Iheatmap-method}
\title{add_row_plot}
\usage{
\S4method{add_row_plot}{Iheatmap}(
  p,
  x,
  ...,
  mode = c("lines+markers", "lines", "markers"),
  color = NULL,
  tracename = NA_character_,
  showlegend = !is.na(tracename),
  side = c("right", "left"),
  layout = list(),
  size = 0.2,
  buffer = 0.02,
  xname = NULL,
  yname = current_yaxis(p),
  pname = if (!is.na(tracename)) tracename else "row_plot"
)
}
\arguments{
\item{p}{iheatmap object}

\item{x}{x axis values}

\item{...}{additional arguments to add to plotly scatter trace, see
\url{https://plotly.com/javascript/reference/#scatter}}

\item{mode}{mode of plot -- one of "lines+markers","lines", or "markers"}

\item{color}{color of bars}

\item{tracename}{name of trace (for legend and hover)}

\item{showlegend}{show in legend?}

\item{side}{side of plot on which to add subplot}

\item{layout}{yaxis layout list}

\item{size}{relative size of subplot relative to main heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{pname}{internal name of plot}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Add a scatter or line plot with one point per row of the main heatmap
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm <- iheatmap(mat) \%>\% add_row_plot(x = 1:4, tracename = "Strength")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_row_signal}}, \code{\link{iheatmap}}, 
\code{\link{add_row_barplot}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{Dendrogram-class}
\alias{Dendrogram-class}
\alias{Dendrogram}
\title{Dendrogram}
\description{
Class for storing dendrogram
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{hclust object}

\item{\code{side}}{side of plot on which dendrogram is positioned, controls 
orientation}
}}

\seealso{
\code{\link{add_row_dendro}}, \code{\link{add_col_dendro}},
\code{\link{add_row_clustering}}, \code{\link{add_col_clustering}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapAnnotation-class}
\alias{IheatmapAnnotation-class}
\alias{IheatmapAnnotation}
\title{IheatmapAnnotation}
\description{
Virtual class for storing annotation objects
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{main annotation data}
}}

\section{SubClasses}{

\itemize{
\item \code{\link{RowTitle-class}}
\item \code{\link{ColumnTitle-class}}
\item \code{\link{RowLabels-class}}
\item \code{\link{ColumnLabels-class}}
}
}

\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barplot.R
\name{add_row_barplot}
\alias{add_row_barplot}
\alias{add_row_barplot,Iheatmap-method}
\title{add_row_barplot}
\usage{
\S4method{add_row_barplot}{Iheatmap}(
  p,
  x,
  ...,
  color = NULL,
  tracename = NA_character_,
  showlegend = !is.na(tracename),
  side = c("right", "left"),
  layout = list(),
  size = 0.2,
  buffer = 0.02,
  xname = NULL,
  yname = current_yaxis(p),
  pname = if (!is.na(tracename)) tracename else "row_barplot"
)
}
\arguments{
\item{p}{iheatmap object}

\item{x}{x axis values}

\item{...}{additional arguments to add to plotly scatter trace, see
\url{https://plotly.com/javascript/reference/#scatter}}

\item{color}{color of bars}

\item{tracename}{name of trace (for legend and hover)}

\item{showlegend}{show in legend?}

\item{side}{side of plot on which to add subplot}

\item{layout}{yaxis layout list}

\item{size}{relative size of subplot relative to main heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{pname}{internal name of plot}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
add_row_barplot
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm <- iheatmap(mat) \%>\% add_row_barplot(x = 1:4, tracename = "Strength")

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_row_signal}}, \code{\link{iheatmap}}, 
\code{\link{add_row_plot}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotations.R
\name{add_col_annotation}
\alias{add_col_annotation}
\alias{add_col_annotation,Iheatmap-method}
\title{add_col_annotation}
\usage{
\S4method{add_col_annotation}{Iheatmap}(
  p,
  annotation,
  colors = NULL,
  side = c("top", "bottom"),
  size = 0.05,
  buffer = 0.015,
  inner_buffer = buffer/2,
  layout = list(),
  show_colorbar = TRUE
)
}
\arguments{
\item{p}{\code{link{Iheatmap-class}} object}

\item{annotation}{data.frame or object that can be converted to data frame}

\item{colors}{list of color palettes, with one color per annotation column 
name}

\item{side}{side of plot on which to add column annotation}

\item{size}{relative size of each row annotation}

\item{buffer}{relative size of buffer between previous subplot and column 
annotation}

\item{inner_buffer}{relative size of buffer between each annotation}

\item{layout}{layout properties for new y axis}

\item{show_colorbar}{logical indicator to show or hide colorbar}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds annotation heatmaps for one or more qualitative or quantitative 
annotations for each column of a main heatmap.
}
\examples{

mat <- matrix(rnorm(24), ncol = 6)
annotation <- data.frame(gender = c(rep("M", 3),rep("F",3)),
                        age = c(20,34,27,19,23,30))
hm <- iheatmap(mat) \%>\% add_col_annotation(annotation)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}, \code{\link{add_row_annotation}}, 
\code{\link{add_col_signal}}, \code{\link{add_col_groups}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_widget.R
\name{save_iheatmap}
\alias{save_iheatmap}
\alias{save_iheatmap,Iheatmap,character-method}
\title{save_iheatmap}
\usage{
\S4method{save_iheatmap}{Iheatmap,character}(p, filename, ...)
}
\arguments{
\item{p}{\code{link{Iheatmap-class}} object}

\item{filename}{name of file}

\item{...}{additional arguments to \code{\link[htmlwidgets]{saveWidget}} for
saving as html or \code{\link[webshot]{webshot}} for saving as pdf/png/jpeg}
}
\description{
save an \code{link{Iheatmap-class}} object, either as standalone HTML or as static
pdf/png/jpeg
}
\details{
Note that this function requires the webshot package. If deploying
a shiny app that calls this function in shinyapps.io, loading the webshot
library and calling \code{webshot::install_phantomjs()} is needed for the the save
functionality to work.
}
\examples{
mat <- matrix(rnorm(24), nrow = 6)
hm <- iheatmap(mat)
\dontrun{
save_iheatmap(hm, "example_iheatmap.png")
}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dendogram.R
\name{add_row_dendro}
\alias{add_row_dendro}
\alias{add_row_dendro,Iheatmap,hclust-method}
\title{add_row_dendro}
\usage{
\S4method{add_row_dendro}{Iheatmap,hclust}(
  p,
  dendro,
  reorder = TRUE,
  side = c("left", "right"),
  size = 0.15,
  buffer = 0.005,
  xname = NULL,
  yname = current_yaxis(p),
  sname = "row_dendro"
)
}
\arguments{
\item{p}{iheatmap object}

\item{dendro}{hclust object}

\item{reorder}{reorder rows based on dendrogram order?}

\item{side}{side of plot on which to add dendrogram}

\item{size}{relative size of dendrogram (relative to the main heatmap)}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal  name of xaxis}

\item{yname}{internal name of yaxis}

\item{sname}{internal name of shapes}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds row dendrogram to iheatmap object
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
dend <- hclust(dist(mat))
hm <- iheatmap(mat) \%>\% add_row_dendro(dend)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_row_clustering}}, \code{\link{iheatmap}}, 
\code{\link{add_col_dendro}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_utils.R
\docType{methods}
\name{iheatmap_list_utils}
\alias{iheatmap_list_utils}
\alias{length,IheatmapList-method}
\alias{as.list,IheatmapList-method}
\alias{names,IheatmapList-method}
\alias{`names<-`,IheatmapList-method}
\alias{`$`,IheatmapList-method}
\alias{`$<-`,IheatmapList-method}
\alias{`[`,IheatmapList-method}
\alias{`[<-`,IheatmapList-method}
\alias{`[[`,IheatmapList-method}
\alias{`[[<-`,IheatmapList-method}
\alias{lapply,IheatmapList-method}
\alias{vapply,IheatmapList-method}
\alias{[,IheatmapList-method}
\alias{[<-,IheatmapList-method}
\alias{[[,IheatmapList-method}
\alias{[[<-,IheatmapList-method}
\alias{$,IheatmapList-method}
\alias{$<-,IheatmapList-method}
\alias{names<-,IheatmapList-method}
\title{S4 List Utils for Iheatmap classes}
\usage{
\S4method{length}{IheatmapList}(x)

\S4method{as.list}{IheatmapList}(x)

\S4method{[}{IheatmapList}(x, i)

\S4method{[}{IheatmapList}(x, i) <- value

\S4method{[[}{IheatmapList}(x, i)

\S4method{[[}{IheatmapList}(x, i) <- value

\S4method{$}{IheatmapList}(x, name)

\S4method{$}{IheatmapList}(x, name) <- value

\S4method{names}{IheatmapList}(x)

\S4method{names}{IheatmapList}(x) <- value

\S4method{lapply}{IheatmapList}(X, FUN, ...)

\S4method{vapply}{IheatmapList}(X, FUN, FUN.VALUE, ..., USE.NAMES = TRUE)
}
\arguments{
\item{x}{input}

\item{FUN}{function to apply to each element of x}

\item{...}{additional arguments}

\item{FUN.VALUE}{template for return value from FUN}

\item{USE.NAMES}{logical, use names?}
}
\description{
These are utility methods for list-like classes in the package.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapColorbar-class}
\alias{IheatmapColorbar-class}
\alias{IheatmapColorbar}
\title{IheatmapColorbar}
\description{
Virtual class for storing colorbar objects
}
\section{Slots}{

\describe{
\item{\code{title}}{title for colorbar}

\item{\code{position}}{integer indicating relative position of colorbar}

\item{\code{colors}}{name of color palette or vector of colors}
}}

\section{SubClasses}{

\itemize{
\item \code{\link{ContinuousColorbar-class}}
\item \code{\link{DiscreteColorbar-class}}
}
}

\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny.R
\name{iheatmapr_event}
\alias{iheatmapr_event}
\title{Access iheatmapr user input event data in shiny}
\usage{
iheatmapr_event(
  object,
  event = c("hover", "click", "relayout"),
  session = shiny::getDefaultReactiveDomain()
)
}
\arguments{
\item{object}{\code{\link{Iheatmap-class}} object}

\item{event}{The type of plotly event. Currently 'plotly_hover',
'plotly_click', 'plotly_selected', and 'plotly_relayout' are supported.}

\item{session}{a shiny session object (the default should almost always be used).}
}
\description{
This function must be called within a reactive shiny context.
}
\examples{
\dontrun{
shiny::runApp(system.file("examples", "shiny_example", package = "iheatmapr"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iheatmap.R
\name{iheatmap}
\alias{iheatmap}
\alias{iheatmap,matrix-method}
\title{iheatmap}
\usage{
\S4method{iheatmap}{matrix}(
  data,
  x = default_x(data),
  y = default_y(data),
  cluster_rows = c("none", "hclust", "kmeans"),
  cluster_cols = c("none", "hclust", "kmeans"),
  row_clusters = NULL,
  col_clusters = NULL,
  row_k = NULL,
  col_k = NULL,
  row_clust_dist = stats::dist,
  col_clust_dist = stats::dist,
  name = "Signal",
  scale = c("none", "rows", "cols"),
  scale_method = c("standardize", "center", "normalize"),
  colors = NULL,
  col_clusters_colors = NULL,
  col_clusters_name = "Col<br>Clusters",
  row_clusters_colors = NULL,
  row_clusters_name = "Row<br>Clusters",
  show_row_clusters_colorbar = TRUE,
  show_col_clusters_colorbar = TRUE,
  row_annotation = NULL,
  col_annotation = NULL,
  row_annotation_colors = NULL,
  col_annotation_colors = NULL,
  row_labels = NULL,
  col_labels = NULL,
  row_title = NULL,
  col_title = NULL,
  colorbar_grid = setup_colorbar_grid(),
  layout = list(),
  source = "iheatmapr",
  ...
)
}
\arguments{
\item{data}{matrix of values to be plotted as heatmap}

\item{x}{x xaxis labels, by default colnames of data}

\item{y}{y axis labels, by default rownames of data}

\item{cluster_rows}{"none","hclust", or "k-means" for no clustering, 
hierarchical clustering, and k-means clustering of rows respectively}

\item{cluster_cols}{"none","hclust", or "k-means" for no clustering, 
hierarchical clustering, and k-means clustering of columnsrespectively}

\item{row_clusters}{vector of pre-determined row cluster assignment}

\item{col_clusters}{vector of pre-determined column cluster assignment}

\item{row_k}{number of clusters for rows, needed if cluster_rows is kmeans or 
optional if hclust}

\item{col_k}{number of clusters for columns, needed if cluster_rows is kmeans
or optional if hclust}

\item{row_clust_dist}{distance function to use for row clustering if 
hierarchical clustering}

\item{col_clust_dist}{distance function to use for column clustering 
if hierarchical clustering}

\item{name}{Name for colorbar}

\item{scale}{scale matrix by rows, cols or none}

\item{scale_method}{what method to use for scaling, either none, standardize, 
center, normalize}

\item{colors}{name of RColorBrewer palette or vector of colors for main 
heatmap}

\item{col_clusters_colors}{colors for col clusters annotation heatmap}

\item{col_clusters_name}{name for col clusters colorbar}

\item{row_clusters_colors}{colors for row clusters annotation heatmap}

\item{row_clusters_name}{name for row clusters colorbar}

\item{show_row_clusters_colorbar}{show the colorbar for row clusters?}

\item{show_col_clusters_colorbar}{show the colorbar for column clusters?}

\item{row_annotation}{row annotation data.frame}

\item{col_annotation}{column annotation data.frame}

\item{row_annotation_colors}{list of colors for row annotations heatmap}

\item{col_annotation_colors}{list of colors for col annotations heatmap}

\item{row_labels}{axis labels for y axis}

\item{col_labels}{axis labels for x axis}

\item{row_title}{x axis title}

\item{col_title}{y axis title}

\item{colorbar_grid}{colorbar grid parameters, should be result from 
\code{\link{setup_colorbar_grid}}}

\item{layout}{list of layout attributes to pass to plotly, 
eg. list(font = list(size = 15))}

\item{source}{source name for use with shiny}

\item{...}{additional argument to iheatmap}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Make a farily standard interactive heatmap with optional clustering and 
row and column annotations.  For more flexibility and options, see the
\code{\link{main_heatmap}} function and other modular functions as described
 in vignette.
}
\details{
By default, no scaling is done of rows or columns. This can be changed by 
specifying the 'scale' argument.  There are three options for scaling 
methods. "standardize" subtracts the mean and divides by standard deviation,
"center" just subtracts the mean, and "normalize" divides by the sum of the 
values.  "normalize" should only be used for data that is all positive!  If 
alternative scaling is desired, the scaling should be done prior to calling 
the iheatmap function.
}
\examples{
mat <- matrix(rnorm(24), nrow = 6)
annotation = data.frame(gender = c(rep("M", 3),rep("F",3)),
 age = c(20,34,27,19,23,30))
hm <- iheatmap(mat, 
 cluster_rows = "hclust",
 cluster_cols = "kmeans", 
 col_k = 3, 
 row_annotation = annotation)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}, \code{\link{add_iheatmap}}, 
\code{\link{to_widget}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{ColumnLabels-class}
\alias{ColumnLabels-class}
\alias{ColumnLabels}
\title{ColumnLabels}
\description{
Class for storing column labels
}
\section{Slots}{

\describe{
\item{\code{xaxis}}{name of xaxis}

\item{\code{yaxis}}{name of yaxis}

\item{\code{data}}{the names of the tick labels}

\item{\code{positions}}{the positions of the tick labels}

\item{\code{side}}{side of plot on which dendrogram is positioned, controls 
orientation}

\item{\code{textangle}}{angle for text}

\item{\code{font}}{list of font attributes}
}}

\seealso{
\code{\link{add_col_labels}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_widget.R
\name{to_widget}
\alias{to_widget}
\alias{to_widget,Iheatmap-method}
\title{to_widget}
\usage{
\S4method{to_widget}{Iheatmap}(p)
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object to convert}
}
\value{
htmlwidgets object
}
\description{
Function to convert \code{link{Iheatmap-class}} object to widget object
}
\examples{

mat <- matrix(rnorm(24), nrow = 6)
hm <- iheatmap(mat) \%>\% to_widget()
class(hm)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}, \code{\link{main_heatmap}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{add_col_summary}
\alias{add_col_summary}
\alias{add_col_summary,Iheatmap-method}
\title{add_col_summary}
\usage{
\S4method{add_col_summary}{Iheatmap}(
  p,
  groups = NULL,
  heatmap_name = NULL,
  colors = NULL,
  tracename = "Col Summary",
  showlegend = FALSE,
  side = c("top", "bottom"),
  layout = list(),
  size = 0.3,
  buffer = 0.02,
  xname = current_xaxis(p),
  yname = NULL,
  type = c("scatter", "bar"),
  summary_function = c("mean", "median", "sd", "var", "mad", "max", "min", "sum"),
  ...
)
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object}

\item{groups}{vector of group labels, name of groups colorbar, or TRUE -- 
see Details}

\item{heatmap_name}{name of a heatmap within the plot}

\item{colors}{vector of colors or RColorBrewer palette name}

\item{tracename}{name of trace}

\item{showlegend}{show legend?}

\item{side}{side of plot on which to add subplot}

\item{layout}{xaxis layout list}

\item{size}{relative size of subplot relative to main heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{type}{scatter or bar?}

\item{summary_function}{summary function to use, default is mean, options 
are mean, median, sd, var, mad, max, min, and sum}

\item{...}{additional arguments to \code{\link{add_col_plot}} or 
\code{\link{add_col_barplot}}}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds a line plot summarizing the values across columns
}
\details{
If adding the column summary to a vertically oriented heatmap, the summary 
will be based on the topmost heatmap if side is "top" and based on the bottom
heatmap if side is "bottom" unless a "heatmap_name" is specified. The 
heatmap_name should match the "pname" argument given to a previously added
heatmap.

The column summary is based on specific rows if a "groups" argument
is given. The groups argument can either be a vector of group assignments for 
each row, the "pname" for an existing set of groups incorporated into the 
plot using \code{\link{add_row_groups}}, \code{\link{add_row_annotation}}, 
\code{\link{add_row_clusters}}, or \code{\link{add_row_clustering}}.  If 
groups is set to TRUE, then the function will use an existing set of row 
groups added to the plot.
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm1 <- iheatmap(mat) \%>\% add_col_summary()
hm2 <- iheatmap(mat) \%>\% add_col_summary(groups = c("A","A","B","B"))

# Print heatmap if interactive session 
if (interactive()) hm1
if (interactive()) hm2
}
\seealso{
\code{\link{add_row_summary}}, \code{\link{iheatmap}}, 
\code{\link{add_col_plot}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{add_row_clustering}
\alias{add_row_clustering}
\alias{add_row_clustering,Iheatmap-method}
\title{add_row_clustering}
\usage{
\S4method{add_row_clustering}{Iheatmap}(
  p,
  method = c("hclust", "kmeans", "groups"),
  name = "Row<br>Clusters",
  k = NULL,
  groups = NULL,
  clust_dist = stats::dist,
  colors = NULL,
  show_colorbar = TRUE,
  side = c("left", "right"),
  xname = NULL,
  yname = current_yaxis(p)
)
}
\arguments{
\item{p}{iheatmap object}

\item{method}{"hclust" or "kmeans" for hierarchical or k-means clustering, 
respectively}

\item{name}{name of colorbar indicating cluster membership}

\item{k}{number of clusters for rows, needed if order is kmeans or optional 
if hclust}

\item{groups}{vector of group assignments}

\item{clust_dist}{distance function to use for clustering if hierarchical 
clustering}

\item{colors}{colors to use for annotation of grouping, can be RColorBrewer 
palette name or
vector of colors}

\item{show_colorbar}{show the colorbar for the heatmap indicating cluster 
membership}

\item{side}{side of plot on which to add subplot}

\item{xname}{name of xaxis}

\item{yname}{name of yaxis}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
add_row_clustering
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)
hm <- iheatmap(mat) \%>\% add_row_clustering(method = "hclust", k = 2)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_col_clustering}}, \code{\link{iheatmap}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{Iheatmap-class}
\alias{Iheatmap-class}
\alias{IheatmapHorizontal}
\alias{IheatmapVertical}
\alias{IheatmapHorizontal-class}
\alias{IheatmapVertical-class}
\title{Iheatmap-class}
\description{
Class to store complex interactive heatmap objects from iheatmapr package
}
\details{
This is a virtual class with two children classes, 
 IheatmapHorizontal and IheatmapVertical. For IheatmapHorizontal additional
 main heatmaps are added horizontally, and for IheatmapVertical additional
 main heatmaps are added vertically. For details on accessing certain slots 
 of this class, see \code{\link{access_component}} documentation.
}
\section{Slots}{

\describe{
\item{\code{plots}}{list of plot element in \code{\link{IheatmapPlots}} format}

\item{\code{shapes}}{list of shape element in \code{\link{IheatmapShapes}} format}

\item{\code{annotations}}{list of annotation elements in 
\code{\link{IheatmapAnnotations}} format}

\item{\code{xaxes}}{list of x axes in \code{\link{IheatmapAxes}} format}

\item{\code{yaxes}}{list of y axes in \code{\link{IheatmapAxes}} format}

\item{\code{colorbars}}{list of colorbars in \code{\link{IheatmapColorbars}} format}

\item{\code{colorbar_grid}}{colorbar grid parameters in
\code{\link{IheatmapColorbarGrid}} format}

\item{\code{current_xaxis}}{name of current x axis}

\item{\code{current_yaxis}}{name of current y axis}

\item{\code{layout}}{list of plotly layout parameters}

\item{\code{source}}{source name, for use with shiny}
}}

\seealso{
\code{\link{iheatmap}}, \code{\link{main_heatmap}}, 
\code{\link{access_component}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{add_row_summary}
\alias{add_row_summary}
\alias{add_row_summary,Iheatmap-method}
\title{add_row_summary}
\usage{
\S4method{add_row_summary}{Iheatmap}(
  p,
  groups = NULL,
  heatmap_name = NULL,
  colors = NULL,
  tracename = "Row Summary",
  showlegend = FALSE,
  side = c("right", "left"),
  layout = list(),
  size = 0.3,
  buffer = 0.02,
  xname = NULL,
  yname = current_yaxis(p),
  type = c("scatter", "bar"),
  summary_function = c("mean", "median", "sd", "var", "mad", "max", "min", "sum"),
  ...
)
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object}

\item{groups}{vector of group labels, name of groups colorbar, or TRUE -- 
see Details}

\item{heatmap_name}{name of a heatmap within the \code{\link{Iheatmap-class}}
object}

\item{colors}{vector of colors or RColorBrewer palette name}

\item{tracename}{name of trace}

\item{showlegend}{show legend?}

\item{side}{side of plot on which to add subplot}

\item{layout}{xaxis layout list}

\item{size}{relative size of subplot relative to main heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{type}{scatter or bar?}

\item{summary_function}{summary function to use, default is mean, options 
are mean, median, sd, var, mad, max, min, and sum}

\item{...}{additional arguments to \code{\link{add_row_plot}} or 
\code{\link{add_row_barplot}}}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds a line plot summarizing the values across rows
}
\details{
If adding the row summary to a horizontally oriented heatmap, the summary 
will be based on the right-most heatmap if side is "right" and based on the 
left heatmap if side is "left" unless a "heatmap_name" is specified. The 
heatmap_name should match the "pname" argument given to a previously added
heatmap.

The row summary is based on specific columns if a "groups" argument
is given. The groups argument can either be a vector of group assignments for 
each row, the "pname" for an existing set of groups incorporated into the 
plot using \code{\link{add_col_groups}}, \code{\link{add_col_annotation}}, 
\code{\link{add_col_clusters}}, or \code{\link{add_col_clustering}}.  If 
groups is set to TRUE, then the function will use an existing set of column 
groups added to the plot.
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
hm1 <- iheatmap(mat) \%>\% add_row_summary()
hm2 <- iheatmap(mat) \%>\% add_row_summary(groups = c("A","A","B","B","B"))

# Print heatmap if interactive session 
if (interactive()) hm1
if (interactive()) hm2
}
\seealso{
\code{\link{add_col_summary}}, \code{\link{iheatmap}}, 
\code{\link{add_row_plot}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\docType{methods}
\name{make_component}
\alias{make_component}
\alias{make_trace,MainHeatmap-method}
\alias{make_trace,RowAnnotation-method}
\alias{make_trace,ColumnAnnotation-method}
\alias{make_trace,RowPlot-method}
\alias{make_trace,ColumnPlot-method}
\alias{make_trace,GenericPlot-method}
\alias{make_shapes,Dendrogram-method}
\alias{make_annotations,RowTitle-method}
\alias{make_annotations,ColumnTitle-method}
\alias{make_annotations,RowLabels-method}
\alias{make_annotations,ColumnLabels-method}
\alias{make_colorbar,ContinuousColorbar,IheatmapColorbarGrid-method}
\alias{make_colorbar,DiscreteColorbar,IheatmapColorbarGrid-method}
\alias{make_trace}
\alias{make_shapes}
\alias{make_annotations}
\alias{make_colorbar}
\title{Convert iheatmapr subcomponents to plotly format}
\usage{
make_trace(x, ...)

make_shapes(x, ...)

make_annotations(x, ...)

make_colorbar(cb, grid)
}
\arguments{
\item{x}{\code{\link{IheatmapPlot-class}}, \code{\link{IheatmapShape-class}},
or \code{\link{IheatmapAnnotation-class}} object}

\item{...}{additional arguments specific to component}
}
\description{
These are generic methods for converting \code{link{Iheatmap-class}}  plot 
components to plotly lists. Not intended for end users; exported for
developers seeking to create new Iheatmap subplots. Any new 
\code{link{IheatmapPlot}}, \code{link{IheatmapShape}},
 \code{link{IheatmapAnnotation}}, or \code{link{IheatmapColorbar}} child class
 should have one of these methods.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{IheatmapColorbarGrid-class}
\alias{IheatmapColorbarGrid-class}
\alias{IheatmapColorbarGrid}
\title{IheatmapColorbarGrid}
\description{
Parameters for setting up colorbars
}
\section{Slots}{

\describe{
\item{\code{nrows}}{number of rows}

\item{\code{x_spacing}}{spacing between colorbars horizontally}

\item{\code{y_spacing}}{spacing between colorbars vertically}

\item{\code{y_length}}{length of colorbars vertically}

\item{\code{x_start}}{start position horizontally}

\item{\code{y_start}}{start position vertically}
}}

\seealso{
\code{\link{setup_colorbar_grid}}, \code{\link{Iheatmap-class}}
}
\author{
Alicia Schep
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusters.R
\name{add_row_clusters}
\alias{add_row_clusters}
\alias{add_row_clusters,Iheatmap-method}
\title{add_row_clusters}
\usage{
\S4method{add_row_clusters}{Iheatmap}(
  p,
  clusters,
  name = "Row<br>Clusters",
  reorder = TRUE,
  side = c("left", "right"),
  yname = current_yaxis(p),
  ...
)
}
\arguments{
\item{p}{iheatmap object}

\item{clusters}{cluster assignments, should be vector of integers, 
characters, or factors}

\item{name}{name of colorbar indicating cluster membership}

\item{reorder}{reorder rows based on clusters? default is TRUE}

\item{side}{side of plot on which to add subplot}

\item{yname}{name of yaxis}

\item{...}{additional arguments to pass to \code{\link{add_row_groups}} for
 creation of annotation
heatmap indicating cluster membership}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Add row groups and order rows based on groups
}
\details{
This function is very similar to \code{\link{add_row_groups}}; the 
main difference is that with this function rows will get reordered based on 
the groups.
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)
clusters <- c("A","B","A","B")

hm <- iheatmap(mat) \%>\% add_row_clusters(clusters)

# Print heatmap if interactive session 
if (interactive()) hm
}
\seealso{
\code{\link{add_row_clustering}}, \code{\link{add_col_clusters}}, 
\code{\link{iheatmap}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny.R
\name{iheatmapr-shiny}
\alias{iheatmapr-shiny}
\alias{iheatmaprOutput}
\alias{renderIheatmap}
\title{Shiny bindings for iheatmapr}
\usage{
iheatmaprOutput(outputId, width = "100\%", height = "400px")

renderIheatmap(expr, env = parent.frame(), quoted = FALSE)
}
\arguments{
\item{outputId}{output variable to read from}

\item{width, height}{Must be a valid CSS unit (like \code{"100\%"},
\code{"400px"}, \code{"auto"}) or a number, which will be coerced to a
string and have \code{"px"} appended.}

\item{expr}{An expression that generates an Iheatmap object}

\item{env}{The environment in which to evaluate \code{expr}.}

\item{quoted}{Is \code{expr} a quoted expression (with \code{quote()})? This 
is useful if you want to save an expression in a variable.}
}
\description{
Output and render functions for using iheatmapr within Shiny
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_heatmap.R
\name{add_main_heatmap}
\alias{add_main_heatmap}
\alias{add_main_heatmap,IheatmapHorizontal,matrix-method}
\alias{add_main_heatmap,IheatmapVertical,matrix-method}
\title{add_main_heatmap}
\usage{
\S4method{add_main_heatmap}{IheatmapHorizontal,matrix}(
  p,
  data,
  name = "Signal",
  x = default_x(data),
  colors = pick_continuous_colors(zmid, zmin, zmax, p),
  colorbar_position = get_colorbar_position(p),
  show_colorbar = TRUE,
  zmin = min(data, na.rm = TRUE),
  zmax = max(data, na.rm = TRUE),
  zmid = 0,
  col_order = NULL,
  x_categorical = NULL,
  side = c("right", "left"),
  size = 1,
  buffer = 0.04,
  text = signif(data, digits = 3),
  tooltip = setup_tooltip_options(),
  xname = NULL,
  pname = name,
  ...
)

\S4method{add_main_heatmap}{IheatmapVertical,matrix}(
  p,
  data,
  name = "Signal",
  y = default_y(data),
  colors = pick_continuous_colors(zmid, zmin, zmax, p),
  colorbar_position = get_colorbar_position(p),
  show_colorbar = TRUE,
  zmin = min(data, na.rm = TRUE),
  zmax = max(data, na.rm = TRUE),
  zmid = 0,
  row_order = NULL,
  y_categorical = NULL,
  side = c("bottom", "top"),
  size = 1,
  buffer = 0.04,
  text = signif(data, digits = 3),
  tooltip = setup_tooltip_options(),
  yname = NULL,
  pname = name,
  ...
)
}
\arguments{
\item{p}{\code{\link{Iheatmap-class}} object}

\item{data}{matrix}

\item{name}{name of colorbar, will determine if colorbar is shared with 
existing plot}

\item{x}{x axis labels (by default rownames of data); only used if 
orientation is horizontal}

\item{colors}{color palette name or vector of colors}

\item{colorbar_position}{colorbar placement}

\item{show_colorbar}{display the colorbar?}

\item{zmin}{minimum for colorscale}

\item{zmax}{maximum for colorscale}

\item{zmid}{midpoint for scale}

\item{col_order}{column ordering for this heatmap; only used if orientation 
is horizontal}

\item{x_categorical}{is x categorical?  will guess if not provided}

\item{side}{which side of the current plot to add this heatmap?}

\item{size}{relative size of plot.  size relative to first heatmap}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{text}{text of value to display for data}

\item{tooltip}{tooltip options, see \code{\link{setup_tooltip_options}}}

\item{xname}{internal name for x axis}

\item{pname}{internal name for plot}

\item{...}{additional arguments (ignored)}

\item{y}{y axis labels (by default colnames of data); only used if 
orientation is vertical}

\item{row_order}{row ordering for this heatmap; only used if orientation is 
vertical}

\item{y_categorical}{is y categorical?  will guess if not provided}

\item{yname}{internal name for y axis}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds an additional main heatmap to an iheatmap object
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4) 
mat2 <-  matrix(rnorm(24), ncol = 6, nrow = 4) 
hm <- iheatmap(mat) \%>\% add_main_heatmap(mat2)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{iheatmap}}, \code{\link{main_heatmap}}
}
\author{
Alicia Schep
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dendogram.R
\name{add_col_dendro}
\alias{add_col_dendro}
\alias{add_col_dendro,Iheatmap,hclust-method}
\title{add_col_dendro}
\usage{
\S4method{add_col_dendro}{Iheatmap,hclust}(
  p,
  dendro,
  reorder = TRUE,
  side = c("top", "bottom"),
  size = 0.15,
  buffer = 0.005,
  xname = current_xaxis(p),
  yname = NULL,
  sname = "col_dendro"
)
}
\arguments{
\item{p}{iheatmap object}

\item{dendro}{hclust object}

\item{reorder}{reorder rows based on dendrogram order?}

\item{side}{side of plot on which to add dendro}

\item{size}{relative size of dendrogram (relative to the main heatmap)}

\item{buffer}{amount of space to leave empty before this plot, relative to 
size of first heatmap}

\item{xname}{internal name of xaxis}

\item{yname}{internal name of yaxis}

\item{sname}{internal name of shape}
}
\value{
\code{\link{Iheatmap-class}} object, which can be printed to generate 
an interactive graphic
}
\description{
Adds column dendrogram to iheatmap object
}
\examples{

mat <- matrix(rnorm(20), ncol = 5, nrow = 4)  
dend <- hclust(dist(t(mat)))
hm <- iheatmap(mat) \%>\% add_col_dendro(dend)

# Print heatmap if interactive session 
if (interactive()) hm 
}
\seealso{
\code{\link{add_col_clustering}}, \code{\link{iheatmap}}, 
\code{\link{add_row_dendro}}
}
\author{
Alicia Schep
}
