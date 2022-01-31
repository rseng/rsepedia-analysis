taxview
=======



[![Project Status: Suspended – Initial development has started, but there has not yet been a stable, usable release; work has been stopped for the time being but the author(s) intend on resuming work.](https://www.repostatus.org/badges/latest/suspended.svg)](https://www.repostatus.org/#suspended)
[![R-check](https://github.com/ropensci/taxview/workflows/R-check/badge.svg)](https://github.com/ropensci/taxview/actions)

Summarise and visualize data sets taxonomically

The proposed workflow with `taxview`:

- input:
  - data.frame with data, indicate which column has names or ids
  - vector/list of names or ids (not associated with data)
- gather taxonomic classification data for each taxon
- from previously collected data, compute summary statistics/etc.
- visualize data among taxonomic groups, etc.

## install


```r
remotes::install_github("ropensci/taxview")
```


```r
library(taxview)
```

## use

get some data


```r
x <- system.file("examples/plant_spp.csv", package = "taxview")
```

prepare data: clean, etc.


```r
dat <- tibble::as_tibble(
 data.table::fread(x, stringsAsFactors = FALSE, 
   data.table = FALSE))
out <- tv_prep_ids(x, ids = dat$id, db = "ncbi")
```

Prepare summary. The output of `tv_summarise()` is an S3 class, with a summary of the groupings.


```r
res <- tv_summarise(out)
res
#> <tv_summary>
#>  no. taxa: 602
#>  by rank: N (22)
#>  by rank name: N (602)
#>  within ranks: N (21)
```

The `$summary` slot has the number of taxa in the dataset


```r
res$summary
#> $spp
#> [1] 602
```

The `$by_rank` slot has the breakdown of taxa within each rank category, as a count and percentage.


```r
res$by_rank
#> # A tibble: 22 x 3
#>    rank      count percent
#>    <chr>     <int>   <dbl>
#>  1 genus       124      21
#>  2 species     112      19
#>  3 family       73      12
#>  4 clade        72      12
#>  5 tribe        58      10
#>  6 subfamily    45       7
#>  7 order        41       7
#>  8 subtribe     18       3
#>  9 class        12       2
#> 10 suborder      9       1
#> # … with 12 more rows
```

The `$by_rank_name` slot has the breakdown of taxa ...


```r
res$by_rank_name
#> # A tibble: 602 x 4
#>    name                     rank      count percent
#>    <chr>                    <chr>     <int>   <dbl>
#>  1 50 kb inversion clade    clade         1       0
#>  2 Abrodictyum              genus         1       0
#>  3 Abrodictyum asae-grayi   species       1       0
#>  4 Acacia                   genus         1       0
#>  5 Acacia jonesii           species       1       0
#>  6 Acacieae                 tribe         1       0
#>  7 Acalypheae               tribe         1       0
#>  8 Acalyphoideae            subfamily     1       0
#>  9 Acridocarpus             genus         1       0
#> 10 Acridocarpus spectabilis species       1       0
#> # … with 592 more rows
```

The `$by_within_rank` slot has the breakdown of number of records within each taxon within each rank grouping.


```r
res$by_within_rank[1:2]
#> $clade
#> # A tibble: 72 x 3
#>    name                         count percent
#>    <chr>                        <int>   <dbl>
#>  1 50 kb inversion clade            1       1
#>  2 Amniota                          1       1
#>  3 apioid superclade                1       1
#>  4 Archelosauria                    1       1
#>  5 Archosauria                      1       1
#>  6 asterids                         1       1
#>  7 Bacteroidetes/Chlorobi group     1       1
#>  8 Bifurcata                        1       1
#>  9 Bilateria                        1       1
#> 10 BOP clade                        1       1
#> # … with 62 more rows
#> 
#> $genus
#> # A tibble: 124 x 3
#>    name          count percent
#>    <chr>         <int>   <dbl>
#>  1 Abrodictyum       1       1
#>  2 Acacia            1       1
#>  3 Acridocarpus      1       1
#>  4 Adiantum          1       1
#>  5 Agapetes          1       1
#>  6 Aglaia            1       1
#>  7 Aiouea            1       1
#>  8 Alepidea          1       1
#>  9 Allium            1       1
#> 10 Allocasuarina     1       1
#> # … with 114 more rows
```

visualize (NOT WORKING YET)


```r
tv_viz(res)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/taxview/issues).
* License: MIT
* Get citation information for `taxview` in R doing `citation(package = 'taxview')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/taxview/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/taxview.git`
* Make sure to track progress upstream (i.e., on our version of `taxview` at `ropensci/taxview`) by doing `git remote add upstream https://github.com/ropensci/taxview.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/taxview`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Prefer to Email? Get in touch: [myrmecocystus@gmail.com](mailto:myrmecocystus@gmail.com)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

```

```
taxview
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE
)
```

[![Project Status: Suspended – Initial development has started, but there has not yet been a stable, usable release; work has been stopped for the time being but the author(s) intend on resuming work.](https://www.repostatus.org/badges/latest/suspended.svg)](https://www.repostatus.org/#suspended)
[![R-check](https://github.com/ropensci/taxview/workflows/R-check/badge.svg)](https://github.com/ropensci/taxview/actions)

Summarise and visualize data sets taxonomically

The proposed workflow with `taxview`:

- input:
  - data.frame with data, indicate which column has names or ids
  - vector/list of names or ids (not associated with data)
- gather taxonomic classification data for each taxon
- from previously collected data, compute summary statistics/etc.
- visualize data among taxonomic groups, etc.

## install

```{r eval=FALSE}
remotes::install_github("ropensci/taxview")
```

```{r}
library(taxview)
```

## use

get some data

```{r}
x <- system.file("examples/plant_spp.csv", package = "taxview")
```

prepare data: clean, etc.

```{r}
dat <- tibble::as_tibble(
 data.table::fread(x, stringsAsFactors = FALSE, 
   data.table = FALSE))
out <- tv_prep_ids(x, ids = dat$id, db = "ncbi")
```

Prepare summary. The output of `tv_summarise()` is an S3 class, with a summary of the groupings.

```{r}
res <- tv_summarise(out)
res
```

The `$summary` slot has the number of taxa in the dataset

```{r}
res$summary
```

The `$by_rank` slot has the breakdown of taxa within each rank category, as a count and percentage.

```{r}
res$by_rank
```

The `$by_rank_name` slot has the breakdown of taxa ...

```{r}
res$by_rank_name
```

The `$by_within_rank` slot has the breakdown of number of records within each taxon within each rank grouping.

```{r}
res$by_within_rank[1:2]
```

visualize (NOT WORKING YET)

```{r eval=FALSE}
tv_viz(res)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/taxview/issues).
* License: MIT
* Get citation information for `taxview` in R doing `citation(package = 'taxview')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxview-package.R
\docType{package}
\name{taxview-package}
\alias{taxview-package}
\alias{taxview}
\title{taxview}
\description{
Tools for Vizualizing Data Taxonomically
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tv_viz.R
\name{tv_viz}
\alias{tv_viz}
\title{Visualize data with respect to taxonomy}
\usage{
tv_viz(x)
}
\arguments{
\item{x}{An object of class \code{tv_summary}}
}
\value{
opens interactive vis
}
\description{
Visualize data with respect to taxonomy
}
\examples{
\dontrun{
x <- system.file("examples/plant_spp.csv", package = "taxview")
dat <- tibble::as_tibble(
 data.table::fread(x, stringsAsFactors = FALSE,
   data.table = FALSE))
out <- tv_prep_ids(x, ids = dat$id, db = "ncbi")
(res <- tv_summarise(out))
tv_viz(res)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tv_prep.R
\name{tv_prep}
\alias{tv_prep}
\alias{tv_prep_names}
\alias{tv_prep_ids}
\title{Prepare taxonomic data}
\usage{
tv_prep_names(x, col = NULL, names = NULL, db = NULL)

tv_prep_ids(x, col = NULL, ids = NULL, db = NULL)
}
\arguments{
\item{x}{(data.frame) input data.frame or file path}

\item{col}{(character) column name holding taxonomic names
or taxonomic ids to use}

\item{names}{(character) column name holding taxonomic names
to use. if given, \code{db} is required}

\item{db}{(character) database IDs from. see below for options}

\item{ids}{(character) column name holding taxonomic IDs
to use. if given, \code{db} is required}
}
\value{
an object of class data.frame
}
\description{
Prepare taxonomic data
}
\section{db options}{

\itemize{
\item \code{bold}: Barcode of Life
\item \code{col}: Catalogue of Life
\item \code{eol}: Encyclopedia of Life
\item \code{gbif}: Global Biodiversity Information Facility
\item \code{iucn}: International Union for Conservation of Nature Red List
\item \code{natserv}: Nature Serve
\item \code{nbn}: National Biodiversity Network (UK)
\item \code{tol}: Tree of Life
\item \code{tropicos}: Tropicos
\item \code{itis}: Integrated Taxonomic Information Service
\item \code{ncbi}: National Center for Biotechnology Information
\item \code{worms}: World Register of Marine Species
}
}

\examples{
\dontrun{
x <- system.file("examples/plant_spp.csv", package = "taxview")

# assuming you only have taxonomic names
# tv_prep_names(x, names = "name")

# if you have taxonomic IDs (from set of allowed databases, see above)
## if a column name
# tv_prep_ids(x, ids = "id", db = "eol")
## if a vector of IDs
dat <- tibble::as_tibble(
 data.table::fread(x, stringsAsFactors = FALSE, 
   data.table = FALSE))
out <- tv_prep_ids(x, ids = dat$id, db = "ncbi")
head(out)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tv_summarise.R
\name{tv_summarise}
\alias{tv_summarise}
\title{Summarise data with respect to taxonomy}
\usage{
tv_summarise(x)
}
\arguments{
\item{x}{(data.frame) input data}
}
\value{
an object of class tv_summary
}
\description{
Summarise data with respect to taxonomy
}
\examples{
\dontrun{
x <- system.file("examples/plant_spp.csv", package = "taxview")
dat <- tibble::as_tibble(
 data.table::fread(x, stringsAsFactors = FALSE, 
   data.table = FALSE))
out <- tv_prep_ids(x, ids = dat$id, db = "ncbi")
(res <- tv_summarise(out))
res$summary
res$by_rank
res$by_rank_name
res$by_within_rank
}
}
