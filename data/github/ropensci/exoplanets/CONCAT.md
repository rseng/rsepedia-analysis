
<!-- README.md is generated from README.Rmd. Please edit that file -->

# exoplanets <img src="man/figures/logo.png" align="right" height=150/>

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/exoplanets/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/exoplanets/actions/workflows/check-pak.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/exoplanets/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/exoplanets?branch=master)
[![Peer
review](https://badges.ropensci.org/309_status.svg)](https://github.com/ropensci/software-review/issues/309)
[![CRAN
status](https://www.r-pkg.org/badges/version/exoplanets)](https://CRAN.R-project.org/package=exoplanets)
[![CRAN\_Download\_Badge](https://cranlogs.r-pkg.org/badges/exoplanets)](https://cran.r-project.org/package=exoplanets)
<!-- badges: end -->

The goal of exoplanets is to provide access to [NASA’s Exoplanet Archive
TAP
Service](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html).
The functionality of this package is minimal and is simply an R
interface to access exoplanet data.

<img src="man/figures/README-unnamed-chunk-2-1.png" title="Exoplanets color coded by discovery method" alt="Exoplanets color coded by discovery method" width="100%" />

## Installation

Install the released version of `exoplanets` from CRAN:

``` r
install.packages("exoplanets")
```

Or you can install from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/exoplanets")
```

## Example

This is a basic example which shows you how to access data from the
[k2names](https://exoplanetarchive.ipac.caltech.edu/docs/API_k2names_columns.html)
table:

``` r
library(exoplanets)

options(
  exoplanets.progress = FALSE, # hide progress
  readr.show_types = FALSE     # hide col spec, requires readr 2.0.0 >=
)

exoplanets("k2names")
#> • https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+k2names&format=csv
#> # A tibble: 449 x 3
#>    epic_id        k2_name  pl_name     
#>    <chr>          <chr>    <chr>       
#>  1 EPIC 246199087 K2-112 f TRAPPIST-1 f
#>  2 EPIC 246199087 K2-112 h TRAPPIST-1 h
#>  3 EPIC 211331236 K2-117 c K2-117 c    
#>  4 EPIC 212398486 K2-125 b K2-125 b    
#>  5 EPIC 217941732 K2-130 b K2-130 b    
#>  6 EPIC 228754001 K2-132 b K2-132 b    
#>  7 EPIC 247887989 K2-133 d K2-133 d    
#>  8 EPIC 247589423 K2-136 b K2-136 b    
#>  9 EPIC 247589423 K2-136 d K2-136 d    
#> 10 EPIC 201912552 K2-18 c  K2-18 c     
#> # … with 439 more rows
```

If you wish, you can select only the columns you need:

``` r
exoplanets("ps", c("pl_name", "hostname"))
#> • https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,hostname+from+ps&format=csv
#> # A tibble: 29,683 x 2
#>    pl_name      hostname  
#>    <chr>        <chr>     
#>  1 Kepler-11 c  Kepler-11 
#>  2 Kepler-11 f  Kepler-11 
#>  3 HAT-P-1 b    HAT-P-1   
#>  4 OGLE-TR-10 b OGLE-TR-10
#>  5 TrES-2 b     TrES-2    
#>  6 WASP-3 b     WASP-3    
#>  7 HD 210702 b  HD 210702 
#>  8 BD-08 2823 b BD-08 2823
#>  9 BD-08 2823 c BD-08 2823
#> 10 HAT-P-30 b   HAT-P-30  
#> # … with 29,673 more rows
```

You can also specify the number of rows returned using `limit`:

``` r
exoplanets("keplernames", limit = 5)
#> • https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+keplernames+top+5&format=csv
#> # A tibble: 5 x 4
#>     kepid koi_name  kepler_name   pl_name      
#>     <dbl> <chr>     <chr>         <chr>        
#> 1 7515212 K00679.02 Kepler-212 b  Kepler-212 b 
#> 2 8210018 K02762.01 Kepler-1341 b Kepler-1341 b
#> 3 9008737 K02768.01 Kepler-404 b  Kepler-404 b 
#> 4 4833421 K00232.05 Kepler-122 f  Kepler-122 f 
#> 5 9963524 K00720.02 Kepler-221 d  Kepler-221 d
```

Information on the tables and columns available can be found with:

``` r
tableinfo
#> # A tibble: 546 x 13
#>    table database_column_… table_label description  in_ps_table in_ps_comp_pars…
#>    <chr> <chr>             <chr>       <chr>        <lgl>       <lgl>           
#>  1 ps    default_flag      Default Pa… Boolean fla… TRUE        FALSE           
#>  2 ps    soltype           Solution T… Disposition… TRUE        FALSE           
#>  3 ps    pl_controv_flag   Controvers… Flag indica… TRUE        TRUE            
#>  4 ps    pl_name           Planet Name Planet name… TRUE        TRUE            
#>  5 ps    hostname          Host Name   Stellar nam… TRUE        TRUE            
#>  6 ps    pl_letter         Planet Let… Letter assi… TRUE        TRUE            
#>  7 ps    hd_name           HD ID       Name of the… TRUE        TRUE            
#>  8 ps    hip_name          HIP ID      Name of the… TRUE        TRUE            
#>  9 ps    tic_id            TIC ID      Name of the… TRUE        TRUE            
#> 10 ps    gaia_id           GAIA ID     Name of the… TRUE        TRUE            
#> # … with 536 more rows, and 7 more variables:
#> #   uncertainties_column_positive_negative <chr>, limit_column <chr>,
#> #   default <lgl>, notes <chr>, displayed_string_name <chr>, flag_column <lgl>,
#> #   number_of_measurements <lgl>
```

## Capabilities

At one time, this package used the *Exoplanet Archive Application
Programming Interface (API)*. Since then, a handful of tables have been
transitioned to the *Table Access Protocol (TAP) service*. More tables
will be transitioned to TAP and as such, this package only supports
queries from TAP. For more information, you can read
<https://exoplanetarchive.ipac.caltech.edu/docs/exonews_archive.html#29April2021>.

## Contributing

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
# exoplanets (development version)

# exoplanets 0.2.2

* `tableinfo` has been updated to include additional tables.

# exoplanets 0.2.1

## Breaking changes

* The `progress` parameter in `exoplanets` is no longer available. It has been replaced as an option that can be set with `options`.

## Minor improvements

* Added `forget_exoplanets` to clear the `exoplanets` cache.
* Added package level documentation with `usethis::use_package_doc`.
* Added `limit` parameter to `exoplanets`.
* Added `quiet` option to suppress progress and query message.
* Cleaned up `tableinfo` to remove trailing spaces and other small improvements.

# exoplanets 0.2.0

* Complete rewrite to work with the new TAP service: https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html
* No longer supporting the older API: https://exoplanetarchive.ipac.caltech.edu/docs/program_interfaces.html
* Use of `httr` for progress indicators and better request handling.
* Added `tableinfo` dataset for table and column information.
* Removed all `exo_` functions, added `exoplanets` as a single function for extracting data.
* Added memoization to `exoplanets`.
* Using `httptest` to avoid potential CRAN errors when checking the package.
* Precompiled vignettes to avoid potential CRAN errors when checking the package.
* Transferred ownership to ropensci organization.

# exoplanets 0.1.0

* Added a `NEWS.md` file to track changes to the package.
## Test environments

* macOS-latest, 'release', github actions
* windows-latest, 'release', github actions
* windows-latest, '3.6', github actions
* ubuntu-18.04, 'devel', github actions
* ubuntu-18.04, 'release', github actions
* ubuntu-18.04, 'oldrel', github actions
* ubuntu-18.04, '3.5', github actions
* ubuntu-18.04, '3.4', github actions
* ubuntu-18.04, '3.3', github actions

## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release. This release updates the `tableinfo` dataset.
# parse_url works

    argument "format" is missing, with no default

---

    argument "format" is missing, with no default

---

    argument "format" is missing, with no default

---

    Expected formats are:
    * csv
    * tsv
    * json

<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

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

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the exoplanets project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if

* you have a question, an use case, or otherwise not a bug or feature request for the software itself.
* you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
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
# exoplanets <img src="man/figures/logo.png" align="right" height=150/>

<!-- badges: start -->
[![R build status](https://github.com/ropensci/exoplanets/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/exoplanets/actions/workflows/check-pak.yaml)
[![Codecov test coverage](https://codecov.io/gh/ropensci/exoplanets/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/exoplanets?branch=master)
[![Peer review](https://badges.ropensci.org/309_status.svg)](https://github.com/ropensci/software-review/issues/309)
[![CRAN status](https://www.r-pkg.org/badges/version/exoplanets)](https://CRAN.R-project.org/package=exoplanets)
[![CRAN_Download_Badge](https://cranlogs.r-pkg.org/badges/exoplanets)](https://cran.r-project.org/package=exoplanets)
<!-- badges: end -->

The goal of exoplanets is to provide access to [NASA's Exoplanet Archive TAP Service](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html). The functionality of this package is minimal and is simply an R interface to access exoplanet data.

```{r graphic, echo=FALSE, message=FALSE, dpi=300, fig.height=4, fig.alt="Exoplanets color coded by discovery method"}
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# library(RColorBrewer)
# 
# exoplanets <- exo::exo() %>% 
#   as_tibble()
# 
# ggexoplanets <- exoplanets %>% 
#   select(pl_bmassj, pl_orbper, pl_discmethod) %>% 
#   drop_na() 
# 
# # cols_vec <- randomcoloR::distinctColorPalette(length(unique(ggexoplanets$pl_discmethod)))
# 
# exoplanets <- exo::exo() %>% 
#   as_tibble()
# 
# cols_vec <- c("#86E57B", "#77AFD7", "#DCD955", "#CFD49C", "#D5CCD1", "#84DECA", "#E15AB9", "#D4816C", "#BA91D1", "#A152DF")
# 
# exoplanets %>% 
#   select(pl_name, pl_bmassj, pl_orbper, pl_discmethod) %>% 
#   drop_na(pl_bmassj, pl_orbper, pl_discmethod) %>% 
#   arrange(pl_name) %>% 
#   ggplot(aes(pl_orbper, pl_bmassj)) +
#   geom_point(aes(fill = pl_discmethod), color = "black", shape = 21, size = 2) +
#   scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) +
#   scale_y_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) +
#   labs(
#     x = "Period (days)",
#     y = "Mass (Jupiter Masses)",
#     fill = "Discovery Method"
#   ) +
#   annotation_logticks() +
#   scale_fill_manual(values = cols_vec) +
#   guides(fill = guide_legend(override.aes = list(size = 4))) +
#   geom_curve(aes(x = 100, y = 0.00070, xend = 6.0996151 + 0.6, yend = 0.00194 - 0.0002), 
#              colour = "#555555", 
#              size = 0.5, 
#              curvature = -0.2,
#              arrow = arrow(length = unit(0.03, "npc"))) + 
#   geom_label(aes(x = 100, y = 0.00070, label = "TRAPPIST-1 e"), 
#              hjust = 0, 
#              vjust = 0.5, 
#              label.size = NA, 
#              family = "SFProText-Regular", 
#              size = 5) +
#   theme_bw() +
#   theme(
#     text = element_text(family = "SFProText-Regular"),
#     plot.caption = element_text(family = "SFProText-RegularItalic"),
#     panel.grid = element_blank()
#   )

knitr::include_graphics("man/figures/README-unnamed-chunk-2-1.png")
```

## Installation

Install the released version of `exoplanets` from CRAN:

```r
install.packages("exoplanets")
```

Or you can install from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/exoplanets")
```

## Example

This is a basic example which shows you how to access data from the [k2names](https://exoplanetarchive.ipac.caltech.edu/docs/API_k2names_columns.html) table:

```{r example, message=FALSE}
library(exoplanets)

options(
  exoplanets.progress = FALSE, # hide progress
  readr.show_types = FALSE     # hide col spec, requires readr 2.0.0 >=
)

exoplanets("k2names")
```

If you wish, you can select only the columns you need: 

```{r ps example}
exoplanets("ps", c("pl_name", "hostname"))
```

You can also specify the number of rows returned using `limit`:

```{r keplernames example}
exoplanets("keplernames", limit = 5)
```

Information on the tables and columns available can be found with:

```{r tableinfo examples}
tableinfo
```

## Capabilities

At one time, this package used the *Exoplanet Archive Application Programming Interface (API)*. Since then, a handful of tables have been transitioned to the *Table Access Protocol (TAP) service*. More tables will be transitioned to TAP and as such, this package only supports queries from TAP. For more information, you can read https://exoplanetarchive.ipac.caltech.edu/docs/exonews_archive.html#29April2021. 

## Contributing

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "exoplanets"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{exoplanets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



In this vignette, we will cover all functionality in `exoplanets` by recreating the discovery plot shown in the packages README file. First, we'll need to load the package.


```r
library(exoplanets)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
```

We will need a couple of columns to recreate the plot:

* The planet name
* The planets mass (jupiter)
* The orbital period (in days)
* The discovery method

To identify the column names that represent these things, we have two options:

* Consult the documentation at https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html
* Review the table information with `tableinfo`

I'd encourage you to read the documentation if you're interested in all the details, it is especially useful to familiarize yourself with what you need. In this example, we are interested in the Planetary Systems (PS) table. Per the documentation you will find:

> PS provides a single table view of the all of the ingested planetary systems for each known exoplanet with each row containing a self-contained set of parameters (planet + stellar + system) for each reference. The PS Table contains one row per planet per reference.

The columns we actually need are:

* `pl_name`
* `pl_orbper`
* `pl_massj`
* `discoverymethod`

We can find these columns and their description in the `tableinfo` dataset.


```r
tableinfo %>%
  filter(table == "ps") %>%
  select(database_column_name, description) %>%
  filter(database_column_name %in% c(
    "pl_name",
    "pl_orbper",
    "pl_massj",
    "discoverymethod"
  ))
#> # A tibble: 4 x 2
#>   database_column_name description                                                          
#>   <chr>                <chr>                                                                
#> 1 pl_name              Planet name most commonly used in the literature                     
#> 2 discoverymethod      Method by which the planet was first identified                      
#> 3 pl_orbper            Time the planet takes to make a complete orbit around the host star …
#> 4 pl_massj             Amount of matter contained in the planet, measured in units of masse…
```

In general, if you have a specific set of columns you need, requesting those columns will be quicker than requesting all columns (default behavior). Let's request only what we need.


```r
discovery <- exoplanets(
  table = "ps",
  columns = c(
    "pl_name",
    "pl_orbper",
    "pl_massj",
    "discoverymethod"
  )
)
#> • https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,pl_orbper,pl_massj,discoverymethod+from+ps&format=csv
#> Rows: 29683 Columns: 4
#> ── Column specification ────────────────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (2): pl_name, discoverymethod
#> dbl (2): pl_orbper, pl_massj
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

discovery
#> # A tibble: 29,683 x 4
#>    pl_name      pl_orbper pl_massj discoverymethod
#>    <chr>            <dbl>    <dbl> <chr>          
#>  1 Kepler-11 c      13.0     0.042 Transit        
#>  2 Kepler-11 f      46.7     0.007 Transit        
#>  3 HAT-P-1 b         4.47    0.532 Transit        
#>  4 OGLE-TR-10 b      3.10    0.62  Transit        
#>  5 TrES-2 b          2.47    1.20  Transit        
#>  6 WASP-3 b          1.85    1.76  Transit        
#>  7 HD 210702 b     354.     NA     Radial Velocity
#>  8 BD-08 2823 b      5.6    NA     Radial Velocity
#>  9 BD-08 2823 c    238.     NA     Radial Velocity
#> 10 HAT-P-30 b        2.81    0.711 Transit        
#> # … with 29,673 more rows
```

Finally, we can recreate the plot.


```r
cols_vec <- c("#86E57B", "#77AFD7", "#DCD955", "#CFD49C", "#D5CCD1", "#84DECA", "#E15AB9", "#D4816C", "#BA91D1", "#A152DF")

discovery %>%
  filter(
    !is.na(pl_massj),
    !is.na(pl_orbper),
    !is.na(discoverymethod)
  ) %>%
  ggplot(aes(pl_orbper, pl_massj)) +
  geom_point(aes(fill = discoverymethod), color = "black", shape = 21, size = 1) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(
    x = "Period (days)",
    y = "Mass (Jupiter Masses)",
    fill = "Discovery Method"
  ) +
  annotation_logticks() +
  scale_fill_manual(values = cols_vec) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  geom_curve(
    aes(x = 100, y = 0.00070, xend = 6.0996151 + 0.6, yend = 0.00194 - 0.0002),
    colour = "#555555",
    size = 0.5,
    curvature = -0.2,
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  geom_label(
    aes(x = 100, y = 0.00070, label = "TRAPPIST-1 e"),
    hjust = 0,
    vjust = 0.5,
    label.size = NA,
    size = 5
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())
```

<img src="exoplanets_vignette_img-discovery-1.png" title="plot of chunk discovery" alt="plot of chunk discovery" width="100%" style="display: block; margin: auto;" />

And in case you were wondering, we've highlighted TRAPPIST-1e, an exoplanet considered to be one of the most potentially habitable exoplanets discovered so far.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exoplanets-package.R
\docType{package}
\name{exoplanets-package}
\alias{exoplanets-package}
\alias{_PACKAGE}
\title{exoplanets: Access NASA's Exoplanet Archive Data}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

The goal of exoplanets is to provide access to
    NASA's Exoplanet Archive TAP Service. For more information regarding
    the API please read the documentation
    <https://exoplanetarchive.ipac.caltech.edu/index.html>.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/exoplanets/}
  \item \url{https://github.com/ropensci/exoplanets}
  \item Report bugs at \url{https://github.com/ropensci/exoplanets/issues}
}

}
\author{
\strong{Maintainer}: Tyler Littlefield \email{tylerlittlefield@hey.com} (\href{https://orcid.org/0000-0002-6020-1125}{ORCID})

Other contributors:
\itemize{
  \item Maëlle Salmon [contributor, reviewer]
  \item Mathida Chuk [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exoplanets.R
\name{exoplanets}
\alias{exoplanets}
\title{Retrieve Data from NASAs Exoplanet Archive}
\source{
\url{https://exoplanetarchive.ipac.caltech.edu/}
}
\usage{
exoplanets(table, columns = NULL, limit = NULL, format = "csv")
}
\arguments{
\item{table}{A table name, see `tableinfo`.}

\item{columns}{A vector of valid column names, by default will return all default columns, see `tableinfo`.}

\item{limit}{Number of rows to return. If NULL, returns all data in the table.}

\item{format}{Desired format, either csv, tsv, or json.}
}
\value{
A \code{data.frame} if \code{format="csv"} or \code{format="tsv"}.
A \code{list} if \code{format="json"}.
}
\description{
A simple interface for accessing exoplanet data. At the bare minimum, a table
name is required. Tables names are documented in the `tableinfo` dataset.
}
\details{
At one time, this package used the Exoplanet Archive Application Programming
Interface (API). Since then, a handful of tables have been transitioned to
the Table Access Protocol (TAP) service. More tables will be transitioned to
TAP and as such, this package only supports queries from TAP. For more
information, you can read
\url{https://exoplanetarchive.ipac.caltech.edu/docs/exonews_archive.html#29April2021.}
}
\examples{
if (interactive()) {
  # request all default columns from the `ps` table
  exoplanets("ps")

  # request the planet name and discovery method from the `ps` table
  exoplanets("ps", c("pl_name", "discoverymethod"))

  # request the first 5 rows from the `keplernames` table
  exoplanets("keplernames", limit = 5)

  # request in json format (returns list)
  exoplanets("ps", c("pl_name", "discoverymethod"), format = "json")

}

}
\seealso{
tableinfo
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{tableinfo}
\alias{tableinfo}
\title{Table Information}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 546 rows and 13 columns.
}
\source{
\url{https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html}
}
\usage{
tableinfo
}
\description{
This dataset provides table information for NASA's Exoplanet
Archive TAP service. In particular, the table name, columns, and column
descriptions are provided.
}
\examples{
tableinfo

}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{forget_exoplanets}
\alias{forget_exoplanets}
\title{Clear the exoplanets cache}
\usage{
forget_exoplanets()
}
\description{
Forget past results and reset the \code{exoplanets} cache.
}
\examples{
if (interactive()) {
  system.time(exoplanets("k2names"))
  system.time(exoplanets("k2names"))
  forget_exoplanets()
  system.time(exoplanets("k2names"))
}

}
