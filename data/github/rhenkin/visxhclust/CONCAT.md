## Resubmission

This is a resubmission. I have fixed:

* Added \value sections to correlation_heatmap.Rd and run_app.Rd
* Replaced \dontrun with if (interactive()) in run_app.Rd example

## Test environments
* macOS Big Sur 10.16, local R installation, R 4.0.5
* ubuntu 16.04 (on travis-ci), R 4.0.5
* win-builder (devel)
* win-builder (release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
# Suggestions

If you want to suggest the inclusion or change of a feature, you can either 1) [fork](https://github.com/rhenkin/visxhclust/fork) the repository to work on it and create an [issue](https://github.com/rhenkin/visxhclust/issues/new) to discuss it before proceeding with a pull request, or 2) create an [issue](https://github.com/rhenkin/visxhclust/issues/new) with your suggestion for others to discuss and potentially work on it.

## Assessing suggestions

Major feature changes or suggestions will be considered based on the following criteria:

1. Scope of the tool: do they significantly expand the scope of the tool and/or are out of place relative to the other features?
1. Visual clutter: do the changes or new inclusions add clutter to the tool or are difficult be visualized with accessible plots?
1. Computation time: are the changes compatible with a fast response time or a progress-bar approach to compensate waiting time?

Minor changes will be discussed following the same criteria but in a less stringent manner.

# Reporting bugs or asking for help

Please report any bugs or ask for help by creating a new [issue](https://github.com/rhenkin/visxhclust/issues/new).

<!-- README.md is generated from README.Rmd. Please edit that file -->

# visxhclust: visual exploration of hierarchical clustering

<!-- badges: start -->

[![R-CMD-check](https://github.com/rhenkin/visxhclust/workflows/R-CMD-check/badge.svg)](https://github.com/rhenkin/visxhclust/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/visxhclust)](https://CRAN.R-project.org/package=visxhclust)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04074/status.svg)](https://doi.org/10.21105/joss.04074)
<!-- badges: end -->

visxhclust is a package that includes a Shiny application for **vis**ual
e**x**ploration of **h**ierarchical **clust**ering. It is aimed at
facilitating iterative workflows of hierarchical clustering on numeric
data. For that, the app allows users to quickly change parameters and
analyse and evaluate results with typical heatmaps with dendrograms and
other charts. Additionally, it includes lightweight data overview plots
such as correlation heatmaps, annotated MDS and PCA plots. On the
evaluation side, it builds on existing packages to compute internal
validation scores and Gap statistic, as well as Dunn’s test to evaluate
significant differences between clusters. Many of the functions are also
exported to facilitate documenting a complete analysis cycle.

## Installation

The latest release can be installed from CRAN:

``` r
install.packages("visxhclust")
```

The latest development version can be installed from GitHub:

``` r
remotes::install_github("rhenkin/visxhclust")
```

Most dependencies are found in CRAN. However, the heatmap drawing
package is part of [Bioconductor](http://www.bioconductor.org/) and may
require a separate installation:

``` r
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

## Getting started

To run the app once the package is installed, use the following
commands:

``` r
library(visxhclust)
# Increases max file size to 30 MB
options(shiny.maxRequestSize = 30*1024^2)
run_app()
```

The app includes multiple help points in the interface (look for the
question marks), and there are also three guides on how to use tool:

-   An [animated
    guide](https://rhenkin.github.io/visxhclust/articles/visxhclust.html)
    on loading data and the basic clustering loop. It’s also accessible
    in R by using the command `vignette("visxhclust")`.
-   An example of how to [reproduce an
    analysis](https://rhenkin.github.io/visxhclust/articles/clusterworkflow.html)
    an analysis using the functions exported by the package. See with
    `vignette("clusterworkflow")` in R.
-   An example of how to [reproduce the evaluation
    workflow](https://rhenkin.github.io/visxhclust/articles/clusterevaluation.html)
    using the functions exported by the package. See with
    `vignette("clusterevaluation")` in R.

## Usage tips and data requirements

To use your data with the tool, you can save a data frame or tibble in
an RDS file, or use comma or tab-delimited files, with .csv, .tsv or
.txt extensions. The clustering method supported by the tool works only
on numeric values; columns containing text will be set aside to annotate
the heatmap if so desired. If a column named `ID` exists, it will be
used as an internal identifier for rows.

Clustering requires complete datasets with no missing values, NULLs or
NAs. If any column contains missing values, it will be set aside to be
used as a heatmap annotation. Badly formatted data will also lead to
unexpected results in the tool. As an alternative, imputation packages
can be used to fill missing data and faulty rows (e.g. text in numeric
columns) should be removed before loading the file into the tool. The
tool provides limited abilities to help with diagnosing issues and
preprocessing data.

# Contributing

Please see the
[guide](https://github.com/rhenkin/visxhclust/blob/master/CONTRIBUTING.md)
for code contribution and suggestions.
# visxhclust 1.0.0

- First version
### Sidebar

This sidebar controls the main clustering workflow: 

1. Load an RDS, CSV, TSV or tab-delimited text file
1. Select a scaling method, distance metric and linkage method
1. Select the desired number of clusters
1. Select the variables/features to use in clustering

Most tabs will then be updated without any additional inputs. The sidebar contains two other controls to help with the analysis:

*Heatmap features:* variables selected from this list will be added to the top of the heatmaps annotations. This list will automatically include any non-numeric variable loaded with the data (e.g. text) or numeric variables that contain missing values, which cannot be used for clustering.

*Correlation-based feature selection:* a one-sided t-test based on the chosen threshold will be use to automatically remove highly correlated variables. You can inspect the relationship between the variables that were deselected in the *Data overview* tab, under *Highly correlated variables*.

### Example data

The following datasets can be downloaded, inspected and loaded into the tool to understand the required format:  

1. [Normally-distributed data](data/sim_normal.csv)  
1. [Normally-distributed data with missing values](data/sim_normal_missing.csv)  
1. [Normally-distributed data with annotation](data/sim_normalannot.csv)  
1. [Binary data](data/sim_binary.csv)  
1. [Log-scaled data](data/logscaled.csv)  
### Internal validation

This tab show the results of a selected internal validation measures for *k* clusters ranging from 2 to 14. The clusterCrit package is used to compute the measures/indices. For more information, please check the documentation of the package at [https://cran.r-project.org/web/packages/clusterCrit/index.html]. To facilitate interpretation, it is possible to highlight first or global minimum or maximum values.

The following is a short summary about which values to look for in each metric:  
**Maximum:** Calinski-Harabasz, Dunn, GDI, Gamma, PBM, Point biserial, Ratwosky-Lance, Silhouette, Tau, Wenmert-Gancarski  
**Minimum:** Banfeld-Raftery, C-index, Davies-Bouldin, G-plus, McClain-Rao, Ray-Turi, Scott-Symons, SD, Xie-Beni  

The remaining measures are the maximum or minimum difference between consecutive k's:  
**Max diff:** Ball-Hall, Ksq, Trace  
**Min diff:** Det ratio, log Det, log SS

### Gap statistic

This tab enables computing the Gap statistic with a user-defined number of boostrap samples. You can also choose a method to highlight the *optimal* number of clusters k. For more information about the underlying functions, check the documentation of the [cluster package](https://cran.r-project.org/package=cluster). For Gap statistic in general, check the original [published article](https://doi.org/10.1111/1467-9868.00293) by Tibshirani et al., 2001.
### Boxplots

This tab contains three views for analysing the clustered data in more detail:

1. **Boxplots/violin plots**: for the variables selected for clustering, this view shows the distribution of the original loaded data (that is, unscaled) across clusters. For small datasets, boxplots overlaid with points are used. For larger datasets, violin plots are used instead. 
1. **Textual summary**: this section shows, for all selected variables, the median, lower and upper quartiles for each cluster.
1. **Annotation distribution**: for variables selected as annotations, either numeric or categorical, this section shows how they are distributed across the clusters
### Managing parameters

Here you can save and load the current settings in the sidebar (distance, linkage and number of clusters), as well as clearing all previously saved settings.
### Data overview

This tab contains three functions:

* **Correlation heatmap:** the heatmap is aimed to support feature selection and displays the correlation among the currently selected features.
* **Highly correlated variables:** a table with variables and correlation that match the threshold set at the bottom of the sidebar, with the corresponding significance level. For pairs of variables that are *not* significantly lower than the threshold, one of the pair is excluded (the second in alphabetical order).
* **Histograms:** to visualize each feature separately and see the effect the of the selected scaling.
### Significance testing

This tab shows the results of Dunn's multiple comparisons test for a specific variable. The selection dropdown will only display variables that were *not* used for computing clusters. Additionally, the test will only work for two or more clusters. The boxplots show, across the clusters, the distributions of the feature chosen for comparison.

The results displayed are the output of the dunn.test function from the [dunn.test](https://cran.r-project.org/web/packages/dunn.test/index.html) package.
### Principal component analysis

This tab has two plots to help understand a bit more the structure of the data used for clustering:

* **2D projection**: this is a projection of two selected PCs annotated with cluster colors.
* **Drivers plot**: this plot shows the strength of the absolute correlation between the original data and the first 8 principal components. The visual part of this plot is inspired by [David Watson's bioplotr package](https://github.com/dswatson/bioplotr).
### Filtering objects

In this tab you can select one cluster to inspect its objects in more detail. You can use Shift or Control/Cmd click to select multiple rows to remove. Note that when you do that, the clusters will be automatically recomputed with the updated dataset. You can also click on "Keep only this cluster" to remove all other objects in dataset that were assigned to the other clusters. **Note: none of these operations can be undone. Reload the file to start again.**
### Distance matrix

This tab contains two views:

* **Projection**: the distance matrix is projected using the classic multidimensional scaling algorithm (stats::cmdscale() in R), and points are annotated with the cluster colors.
* **Heatmap**: the distance matrix is projected as a 2D heatmap, showing pairwise distances. Note that for larger datasets this will be very slow.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# visxhclust: visual exploration of hierarchical clustering
<!-- badges: start -->
[![R-CMD-check](https://github.com/rhenkin/visxhclust/workflows/R-CMD-check/badge.svg)](https://github.com/rhenkin/visxhclust/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/visxhclust)](https://CRAN.R-project.org/package=visxhclust)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04074/status.svg)](https://doi.org/10.21105/joss.04074)
<!-- badges: end -->

visxhclust is a package that includes a Shiny application for **vis**ual e**x**ploration of **h**ierarchical **clust**ering. It is aimed at facilitating iterative workflows of hierarchical clustering on numeric data. For that, the app allows users to quickly change parameters and analyse and evaluate results with typical heatmaps with dendrograms and other charts. Additionally, it includes lightweight data overview plots such as correlation heatmaps, annotated MDS and PCA plots. On the evaluation side, it builds on existing packages to compute internal validation scores and Gap statistic, as well as Dunn's test to evaluate significant differences between clusters. Many of the functions are also exported to facilitate documenting a complete analysis cycle.

## Installation

The latest release can be installed from CRAN:

```{r cran, eval = FALSE}
install.packages("visxhclust")
```

The latest development version can be installed from GitHub:

```{r installation, eval = FALSE}
remotes::install_github("rhenkin/visxhclust")
```

Most dependencies are found in CRAN. However, the heatmap drawing package is part of [Bioconductor](http://www.bioconductor.org/) and may require a separate installation:

```{r bioconductor, eval = FALSE}
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

## Getting started

To run the app once the package is installed, use the following commands:

```{r example, eval = FALSE}
library(visxhclust)
# Increases max file size to 30 MB
options(shiny.maxRequestSize = 30*1024^2)
run_app()
```

The app includes multiple help points in the interface (look for the question marks), and there are also three guides on how to use tool:

- An [animated guide](https://rhenkin.github.io/visxhclust/articles/visxhclust.html) on loading data and the basic clustering loop. It's also accessible in R by using the command `vignette("visxhclust")`.
- An example of how to [reproduce an analysis](https://rhenkin.github.io/visxhclust/articles/clusterworkflow.html) an analysis using the functions exported by the package. See with `vignette("clusterworkflow")` in R.
- An example of how to [reproduce the evaluation workflow](https://rhenkin.github.io/visxhclust/articles/clusterevaluation.html) using the functions exported by the package. See with `vignette("clusterevaluation")` in R.

## Usage tips and data requirements

To use your data with the tool, you can save a data frame or tibble in an RDS file, or use comma or tab-delimited files, with .csv, .tsv or .txt extensions. The clustering method supported by the tool works only on numeric values; columns containing text will be set aside to annotate the heatmap if so desired. If a column named `ID` exists, it will be used as an internal identifier for rows.

Clustering requires complete datasets with no missing values, NULLs or NAs. If any column contains missing values, it will be set aside to be used as a heatmap annotation. Badly formatted data will also lead to unexpected results in the tool. As an alternative, imputation packages can be used to fill missing data and faulty rows (e.g. text in numeric columns) should be removed before loading the file into the tool. The tool provides limited abilities to help with diagnosing issues and preprocessing data.

# Contributing

Please see the [guide](https://github.com/rhenkin/visxhclust/blob/master/CONTRIBUTING.md) for code contribution and suggestions.
---
title: "Documenting a workflow"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Documenting a workflow}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This document shows how to use some functions included in the package to document and reproduce a clustering workflow from the app. Although these functions do not cover all the steps such as selecting features, they allow users of the Shiny app to show the resulting heatmaps and boxplots. The example dataset used here is the `iris` dataset.


```{r setup, message=FALSE, warning=FALSE}
library(visxhclust)
library(dplyr)
library(ggplot2)
```

### Preparing

First we split the numeric and categorical variables and scale the data.

```{r first_step}
numeric_data <- iris %>% select(where(is.numeric))
annotation_data <- iris %>% select(where(is.factor))
```

Let's check the dataset for highly correlated variables that will likely skew the clusters with redundant information:

```{r correlation}
correlation_heatmap(numeric_data)
```

As seen above, petal length and width are highly correlated, so we keep only one of them:

```{r subset}
subset_data <- numeric_data %>% select(Sepal.Length, Sepal.Width, Petal.Width)
```

### Clustering

The clustering itself takes three steps: computing a distance matrix, computing the hierarchical clusters and cutting the tree to find the desired number of clusters. In the app, each of these steps has matching parameters: apply scaling and distance/similarity metric, linkage method and the number of clusters. 

```{r params}
scaling <- TRUE
distance_method <- "euclidean"
linkage_method <- "ward.D2"
# this assumes that, in the app, we identified 3 as the optimal number of clusters
k <- 3 
```

These parameters are used in three functions that the app also uses: `compute_dmat`, `compute_clusters` and `cut_clusters`. You can check the documentation for each function in the package website, or interactively through `?compute_dmat`.

```{r computation}
dmat <- compute_dmat(subset_data, distance_method, TRUE)
clusters <- compute_clusters(dmat, linkage_method)
cluster_labels <- cut_clusters(clusters, k)
```

### Results

Now we can check both the heatmap+dendrogram and boxplots. A function that covers most steps to produce the heatmap is included in the package, with the name: `cluster_heatmaps()`. It plots the dendrogram, the annotation layer, the clustered data heatmap and the heatmap with the rest of the data not used for clustering. In the Shiny app this is done automatically, but outside, plotting the annotation and the unselected data are optional steps; the annotations require an extra step with the function `create_annotations()`. The colors used in the app are also exported by the package as the variable `cluster_colors`.

```{r heatmap}
species_annotation <- create_annotations(iris, "Species")
cluster_heatmaps(scale(subset_data), 
                 clusters,
                 k,
                 cluster_colors,
                 annotation = species_annotation)
```
In addition to the heatmap, the boxplots in the app are also available through functions. There are two steps required to show data through box plots: annotating the original data with the cluster and plotting it.

```{r boxplots}
annotated_data <- annotate_clusters(subset_data, cluster_labels, TRUE)
cluster_boxplots(annotated_data)
```
---
title: "Documenting evaluation"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Documenting evaluation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The other vignette focuses on reproducing a single clustering workflow that assumes that the number of clusters has been decided. As the app includes a few options for evaluating clusters, some of the functions are also made available in the package. The output of the clustering functions can also be used with other packages.

```{r setup, message=FALSE, warning=FALSE}
library(visxhclust)
library(dplyr)
```

### Preprocessing and clustering

```{r prep}
numeric_data <- iris %>% select(Sepal.Length, Sepal.Width, Petal.Width)
dmat <- compute_dmat(numeric_data, "euclidean", TRUE)
clusters <- compute_clusters(dmat, "complete")
```

### Gap statistic

For Gap statistic, the optimal number of clusters depends on the method use to compare cluster solutions. The package cluster includes the function `cluster::maxSE()` to help with that.

```{r gapstat}
gap_results <- compute_gapstat(scale(numeric_data), clusters)
optimal_k <- cluster::maxSE(gap_results$gap, gap_results$SE.sim)
line_plot(gap_results, "k", "gap", xintercept = optimal_k)
```

### Other measures

The Shiny app also includes various other measures computed by [clusterCrit::intCriteria()]. The function `compute_metric` works similarly to `compute_gapstat`, whereas `optimal_score` is similar to maxSE. However, `optimal_score` varies only between first and global minimum and maximum.

```{r dunn}
res <- compute_metric(scale(numeric_data), clusters, "Dunn")
optimal_k <- optimal_score(res$score)
line_plot(res, "k", "score", optimal_k)
```
---
title: "visxhclust Shiny app tutorial"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{visxhclust Shiny app tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This is a simple visual tutorial on the basic loop of the visxhclust Shiny app -- note that almost every tab in the application has a corresponding help icon with more information and tips. To open the Shiny app you need to run:

```{r run, eval=FALSE}
library(visxhclust)
run_app()
```

## Steps for simple iteration

### Step 1: Loading data and setting parameters

![Step 1](load_data.gif)

### Step 2: View clustering results

![Step 2](view_results.gif)

### Step 3: Evaluate

![Step 3](evaluate.gif)

### Step 4: Change parameters again

![Step 4](iterate.gif)

### Step 5: Review results

![Step 5](review_results.gif)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_clustering.R
\name{compute_metric}
\alias{compute_metric}
\title{Compute an internal evaluation metric for clustered data}
\usage{
compute_metric(df, clusters, metric_name, max_k = 14)
}
\arguments{
\item{df}{data frame used to compute clusters}

\item{clusters}{output of \code{\link[=compute_clusters]{compute_clusters()}} or \code{\link[fastcluster:hclust]{fastcluster::hclust()}}}

\item{metric_name}{valid metric name from \code{\link[clusterCrit:getCriteriaNames]{clusterCrit::getCriteriaNames()}} (with TRUE argument)}

\item{max_k}{maximum number of clusters to cut using \code{\link[dendextend:cutree-methods]{dendextend::cutree()}}. Default is 14.}
}
\value{
a data frame with columns \code{k} and \code{score}
}
\description{
Metric will be computed from 2 to max_k clusters. Note that the row number in results will be different from k.
}
\examples{
data_to_cluster <- iris[c("Petal.Length", "Sepal.Length")]
dmat <- compute_dmat(data_to_cluster, "euclidean", TRUE)
clusters <- compute_clusters(dmat, "complete")
compute_metric(scale(data_to_cluster), clusters, "Dunn")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{cluster_boxplots}
\alias{cluster_boxplots}
\title{Plot boxplots with clusters}
\usage{
cluster_boxplots(annotated_data, ...)
}
\arguments{
\item{annotated_data}{data frame returned by \code{annotate_clusters()}}

\item{...}{arguments passed to \code{facet_boxplot()}}
}
\value{
boxplots faceted by clusters
}
\description{
This is a convenience wrapper function for \code{facet_boxplot()}.
Combined with \code{annotate_clusters()}, it
doesn't require specifying axes in \code{facet_boxplot()}.
}
\examples{
dmat <- compute_dmat(iris, "euclidean", TRUE, c("Petal.Length", "Sepal.Length"))
clusters <- compute_clusters(dmat, "complete")
cluster_labels <- cut_clusters(clusters, 2)
annotated_data <- annotate_clusters(iris[, c("Petal.Length", "Sepal.Length")], cluster_labels)
cluster_boxplots(annotated_data, boxplot_colors = visxhclust::cluster_colors)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_annotation_dist}
\alias{plot_annotation_dist}
\title{Plot distribution of annotation data across clusters}
\usage{
plot_annotation_dist(annotations_df, cluster_labels, selected_clusters = NULL)
}
\arguments{
\item{annotations_df}{data frame with variables not used in clustering}

\item{cluster_labels}{output from \code{\link[=cut_clusters]{cut_clusters()}}}

\item{selected_clusters}{optional vector of cluster labels to include in plots}
}
\value{
a \code{patchwork} object
}
\description{
Plot distribution of annotation data across clusters
}
\examples{
dmat <- compute_dmat(iris, "euclidean", TRUE, c("Petal.Length", "Sepal.Length"))
clusters <- compute_clusters(dmat, "complete")
cluster_labels <- cut_clusters(clusters, 2)
plot_annotation_dist(iris["Species"], cluster_labels)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmaps.R
\docType{data}
\name{cluster_colors}
\alias{cluster_colors}
\title{List of colors used in the Shiny app for clusters}
\format{
An object of class \code{character} of length 39.
}
\usage{
cluster_colors
}
\description{
List of colors used in the Shiny app for clusters
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmaps.R
\name{correlation_heatmap}
\alias{correlation_heatmap}
\title{Plot a correlation heatmap}
\usage{
correlation_heatmap(df)
}
\arguments{
\item{df}{numeric data frame to compute correlations}
}
\value{
a \link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap}
}
\description{
Computes pairwise Pearson correlation; if there are fewer than 15 columns, prints
the value of the correlation coefficient inside each tile.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_app.R
\name{run_app}
\alias{run_app}
\title{Runs the Shiny app}
\usage{
run_app()
}
\value{
No return value, runs the app by passing it to print
}
\description{
Runs the Shiny app
}
\examples{
## Only run this example in interactive R sessions
if (interactive()) {
library(visxhclust)
run_app()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_clustering.R
\name{compute_gapstat}
\alias{compute_gapstat}
\title{Compute Gap statistic for clustered data}
\usage{
compute_gapstat(df, clusters, gap_B = 50, max_k = 14)
}
\arguments{
\item{df}{the data used to compute clusters}

\item{clusters}{output of \code{\link[=compute_clusters]{compute_clusters()}} or \code{\link[fastcluster:hclust]{fastcluster::hclust()}}}

\item{gap_B}{number of bootstrap samples for \code{\link[cluster:clusGap]{cluster::clusGap()}} function. Default is 50.}

\item{max_k}{maximum number of clusters to compute the statistic. Default is 14.}
}
\value{
a data frame with the Tab component of \code{\link[cluster:clusGap]{cluster::clusGap()}} results
}
\description{
Compute Gap statistic for clustered data
}
\examples{
data_to_cluster <- iris[c("Petal.Length", "Sepal.Length")]
dmat <- compute_dmat(data_to_cluster, "euclidean", TRUE)
clusters <- compute_clusters(dmat, "complete")
gap_results <- compute_gapstat(scale(data_to_cluster), clusters)
head(gap_results)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{facet_boxplot}
\alias{facet_boxplot}
\title{Faceted boxplots with points or violin plots}
\usage{
facet_boxplot(
  df,
  x,
  y,
  facet_var = NULL,
  boxplot_colors = NULL,
  shape = c("boxplot", "violin"),
  plot_points = TRUE
)
}
\arguments{
\item{df}{a data frame containing all the variables matching the remaining arguments}

\item{x}{categorical variable}

\item{y}{continuous variable}

\item{facet_var}{optional variable to facet data}

\item{boxplot_colors}{list of colors to use as fill for boxplots}

\item{shape}{either "boxplot" or "violin"}

\item{plot_points}{boolean variable to overlay jittered points or not. Default is \code{TRUE}}
}
\value{
a \link[ggplot2:ggplot]{ggplot2::ggplot} object
}
\description{
Faceted boxplots with points or violin plots
}
\examples{
facet_boxplot(iris, x = "Species", y = "Sepal.Length", facet_var = "Species")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_clustering.R
\name{compute_dmat}
\alias{compute_dmat}
\title{Compute a distance matrix from scaled data}
\usage{
compute_dmat(
  x,
  dist_method = "euclidean",
  apply_scaling = FALSE,
  subset_cols = NULL
)
}
\arguments{
\item{x}{a numeric data frame or matrix}

\item{dist_method}{a distance measure to apply to the scaled data. Must be those supported by \code{\link[stats:dist]{stats::dist()}}, plus \code{"mahalanobis"} and \code{"cosine"}. Default is \code{"euclidean"}.}

\item{apply_scaling}{use TRUE to apply \code{\link[base:scale]{base::scale()}}. By default does not scale data.}

\item{subset_cols}{(optional) a list of columns to subset the data}
}
\value{
an object of class "dist" (see \code{\link[stats:dist]{stats::dist()}})
}
\description{
This function applies scaling to the columns of a data frame and
computes and returns a distance matrix from a chosen distance measure.
}
\examples{
dmat <- compute_dmat(iris, "euclidean", TRUE, c("Petal.Length", "Sepal.Length"))
print(class(dmat))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmaps.R
\name{cluster_heatmaps}
\alias{cluster_heatmaps}
\title{Plot heatmap with cluster results and dendrogram}
\usage{
cluster_heatmaps(
  scaled_selected_data,
  clusters,
  k,
  cluster_colors,
  scaled_unselected_data = NULL,
  annotation = NULL
)
}
\arguments{
\item{scaled_selected_data}{scaled matrix or data frame with variables used for clustering}

\item{clusters}{hierarchical cluster results produced by \code{\link[fastcluster:hclust]{fastcluster::hclust()}}}

\item{k}{targeted number of clusters}

\item{cluster_colors}{list of cluster colors to match with boxplots}

\item{scaled_unselected_data}{(optional) scaled matrix or data frame with variables not used for clustering}

\item{annotation}{(optional) \link[ComplexHeatmap:columnAnnotation]{ComplexHeatmap::columnAnnotation} object}
}
\value{
a \link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap}
}
\description{
Plot heatmap with cluster results and dendrogram
}
\examples{
dmat <- compute_dmat(iris, "euclidean", TRUE, c("Petal.Length", "Sepal.Length"))
clusters <- compute_clusters(dmat, "complete")
species_annotation <- create_annotations(iris, "Species")
cluster_heatmaps(scale(iris[c("Petal.Length", "Sepal.Length")]),
                 clusters,
                 3,
                 visxhclust::cluster_colors,
                 annotation = species_annotation)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{logscaled_df}
\alias{logscaled_df}
\title{Simulated logscaled data}
\format{
A data frame with 200 rows and 10 variables:
\describe{
\item{a}{variable a}
\item{b}{variable b}
\item{c}{variable c}
\item{d}{variable d}
\item{e}{variable e}
\item{f}{variable f}
\item{g}{variable g}
\item{h}{variable h}
\item{i}{variable i}
\item{j}{variable j}
}
}
\source{
package author
}
\usage{
logscaled_df
}
\description{
Simulated logscaled data
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_clustering.R
\name{annotate_clusters}
\alias{annotate_clusters}
\title{Annotate data frame with clusters}
\usage{
annotate_clusters(df, cluster_labels, long = TRUE, selected_clusters = NULL)
}
\arguments{
\item{df}{a data frame}

\item{cluster_labels}{list of cluster labels, automatically converted to factor.}

\item{long}{if \code{TRUE}, returned data frame will be in long format. See details for spec. Default is \code{TRUE}.}

\item{selected_clusters}{optional cluster labels to filter}
}
\value{
a wide or long data frame
}
\description{
Annotate data frame with clusters
}
\details{
Long data frame will have columns: \code{Cluster}, \code{Measurement} and \code{Value}.
}
\examples{
dmat <- compute_dmat(iris, "euclidean", TRUE, c("Petal.Length", "Sepal.Length"))
res <- compute_clusters(dmat, "complete")
cluster_labels <- cut_clusters(res, 2)
annotated_data <- annotate_clusters(iris[, c("Petal.Length", "Sepal.Length")], cluster_labels)
head(annotated_data)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmaps.R
\name{plot_cluster_heatmaps}
\alias{plot_cluster_heatmaps}
\title{Draw two heatmaps}
\usage{
plot_cluster_heatmaps(
  top_matrix,
  bottom_matrix,
  dendrograms,
  clusters_set,
  annotation = NULL,
  scaled = FALSE,
  distance_method = NULL,
  cluster_features = TRUE,
  show_col_names = TRUE
)
}
\arguments{
\item{top_matrix}{matrix with selected variables}

\item{bottom_matrix}{matrix with unselected variables}

\item{dendrograms}{to draw above top matrix}

\item{clusters_set}{list of cluster indices}

\item{annotation}{(optional) any kind of annotation object to draw as top_annotation}

\item{scaled}{(optional) boolean to modify colour scale if data has already been scaled}

\item{distance_method}{(optional) if "Binary", use discrete colors for heatmap}

\item{cluster_features}{(optional) If FALSE, row order does not change}

\item{show_col_names}{(optional) If FALSE, does not show column names at base of heatmap}
}
\value{
two concatenated heatmaps drawn with ComplexHeatmap::draw
}
\description{
Draw two heatmaps
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{normal_missing}
\alias{normal_missing}
\title{Simulated normal data with missing values}
\format{
A data frame with 200 rows and 10 variables:
\describe{
\item{a}{variable a}
\item{b}{variable b}
\item{c}{variable c}
\item{d}{variable d}
\item{e}{variable e}
\item{f}{variable f}
\item{g}{variable g}
\item{h}{variable h}
\item{i}{variable i}
\item{j}{variable with randomly missing values}
}
}
\source{
package author
}
\usage{
normal_missing
}
\description{
Simulated normal data with missing values
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_clustering.R
\name{cut_clusters}
\alias{cut_clusters}
\title{Cut a hierarchical tree targeting k clusters}
\usage{
cut_clusters(clusters, k)
}
\arguments{
\item{clusters}{cluster results, produced by e.g. \code{\link[fastcluster:hclust]{fastcluster::hclust()}}}

\item{k}{target number of clusters}
}
\value{
cluster labels
}
\description{
Cut a hierarchical tree targeting k clusters
}
\examples{
dmat <- compute_dmat(iris, "euclidean", TRUE, c("Petal.Length", "Sepal.Length"))
clusters <- compute_clusters(dmat, "complete")
cluster_labels <- cut_clusters(clusters, 2)
head(cluster_labels)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_clustering.R
\name{optimal_score}
\alias{optimal_score}
\title{Find minimum or maximum score in a vector}
\usage{
optimal_score(x, method = c("firstmax", "globalmax", "firstmin", "globalmin"))
}
\arguments{
\item{x}{a numeric vector}

\item{method}{one of "firstmax", "globalmax", "firstmin" or "globalmin"}
}
\value{
the index (not k) of the identified maximum or minimum score
}
\description{
This function is meant to be used with compute_metric. For Gap statistic,
use \code{\link[cluster:clusGap]{cluster::maxSE()}}.
}
\examples{
data_to_cluster <- iris[c("Petal.Length", "Sepal.Length")]
dmat <- compute_dmat(data_to_cluster, "euclidean", TRUE)
clusters <- compute_clusters(dmat, "complete")
res <- compute_metric(scale(data_to_cluster), clusters, "Dunn")
optimal_score(res$score, method = "firstmax")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{normal_annotated}
\alias{normal_annotated}
\title{Simulated normal data with annotations}
\format{
A data frame with 200 rows and 10 variables:
\describe{
\item{a}{variable a}
\item{b}{variable b}
\item{c}{variable c}
\item{d}{variable d}
\item{e}{variable e}
\item{f}{variable f}
\item{g}{variable g}
\item{h}{variable h}
\item{i}{variable i}
\item{j}{variable j}
\item{annot}{annotation column}
}
}
\source{
package author
}
\usage{
normal_annotated
}
\description{
Simulated normal data with annotations
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{line_plot}
\alias{line_plot}
\title{A custom line plot with optional vertical line}
\usage{
line_plot(df, x, y, xintercept = NULL)
}
\arguments{
\item{df}{data source}

\item{x}{variable for horizontal axis}

\item{y}{variable for vertical axis}

\item{xintercept}{optional value in horizontal axis to highlight}
}
\value{
a \link[ggplot2:ggplot]{ggplot2::ggplot} object
}
\description{
A custom line plot with optional vertical line
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{normal_df}
\alias{normal_df}
\title{Simulated normal data}
\format{
A data frame with 200 rows and 10 variables:
\describe{
\item{a}{variable a}
\item{b}{variable b}
\item{c}{variable c}
\item{d}{variable d}
\item{e}{variable e}
\item{f}{variable f}
\item{g}{variable g}
\item{h}{variable h}
\item{i}{variable i}
\item{j}{variable j}
}
}
\source{
package author
}
\usage{
normal_df
}
\description{
Simulated normal data
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmaps.R
\name{create_annotations}
\alias{create_annotations}
\title{Create heatmap annotations from selected variables}
\usage{
create_annotations(df, selected_variables)
}
\arguments{
\item{df}{a data frame. It can be an original unscaled data, or a scaled one}

\item{selected_variables}{list of columns in the data frame to create annotations for}
}
\value{
a \link[ComplexHeatmap:columnAnnotation]{ComplexHeatmap::columnAnnotation} object
}
\description{
This function will create a \link[ComplexHeatmap:columnAnnotation]{ComplexHeatmap::columnAnnotation} object with rows
for each variable passed as argument. Character columns will be coerced into factors.
For factors, the ColorBrewer palette \code{Set3} will be used. For non-negative numeric, the
\code{PuBu} palette will be used, and for columns with negative values, the reversed \code{RdBu} will be used.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{bin_df}
\alias{bin_df}
\title{Simulated binary data}
\format{
A data frame with 200 rows and 10 variables:
\describe{
\item{a}{variable a}
\item{b}{variable b}
\item{c}{variable c}
\item{d}{variable d}
\item{e}{variable e}
\item{f}{variable f}
\item{g}{variable g}
\item{h}{variable h}
\item{i}{variable i}
\item{j}{variable j}
}
}
\source{
package author
}
\usage{
bin_df
}
\description{
Simulated binary data
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_clustering.R
\name{compute_clusters}
\alias{compute_clusters}
\title{Compute clusters hierarchically from distance matrix}
\usage{
compute_clusters(dmat, linkage_method)
}
\arguments{
\item{dmat}{a distance matrix}

\item{linkage_method}{a linkage method supported by \code{\link[fastcluster:hclust]{fastcluster::hclust()}}}
}
\value{
clusters computed by \code{\link[fastcluster:hclust]{fastcluster::hclust()}}
}
\description{
Compute clusters hierarchically from distance matrix
}
\examples{
dmat <- compute_dmat(iris, "euclidean", TRUE, c("Petal.Length", "Sepal.Length"))
res <- compute_clusters(dmat, "complete")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{dmat_projection}
\alias{dmat_projection}
\title{Plot a 2D MDS projection of a distance matrix}
\usage{
dmat_projection(dmat, point_colors = NULL, point_palette = NULL)
}
\arguments{
\item{dmat}{distance matrix}

\item{point_colors}{optional list of labels to color points (will be coerced to factor)}

\item{point_palette}{optional palette used with \code{\link[ggplot2:scale_manual]{ggplot2::scale_colour_manual()}}}
}
\value{
a ggplot object
}
\description{
Plot a 2D MDS projection of a distance matrix
}
\examples{
dmat <- dist(iris[, c("Sepal.Width", "Sepal.Length")])
dmat_projection(dmat)
}
