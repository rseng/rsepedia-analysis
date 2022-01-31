
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eia <img src="man/figures/logo.png" style="margin-left:10px;margin-bottom:5px;" width="120" align="right">

**Author:** [Matthew Leonawicz](https://github.com/leonawicz)
<a href="https://orcid.org/0000-0001-9452-2771" target="orcid.widget">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
<br/> **License:** [MIT](https://opensource.org/licenses/MIT)<br/>

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis build
status](https://travis-ci.org/ropensci/eia.svg?branch=master)](https://travis-ci.org/ropensci/eia)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ropensci/eia?branch=master&svg=true)](https://ci.appveyor.com/project/leonawicz/eia)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/eia/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/eia?branch=master)

[![](https://badges.ropensci.org/342_status.svg)](https://github.com/ropensci/software-review/issues/342)
[![CRAN
status](https://www.r-pkg.org/badges/version/eia)](https://cran.r-project.org/package=eia)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/eia)](https://cran.r-project.org/package=eia)
[![Github
Stars](https://img.shields.io/github/stars/ropensci/eia.svg?style=social&label=Github)](https://github.com/ropensci/eia/)

The `eia` package provides API access to data from the US [Energy
Information Administration](https://www.eia.gov/) (EIA).

Pulling data from the US Energy Information Administration (EIA) API
requires a registered API key. A key can be obtained at no cost
[here](https://www.eia.gov/opendata/register.php). A valid email and
agreement to the API Terms of Service is required to obtain a key.

`eia` includes functions for searching EIA API data categories and
importing time series and geoset time series datasets. Datasets returned
by these functions are provided in a tidy format or alternatively in
more raw form. It also offers helper functions for working with EIA API
date strings and time formats and for inspecting different summaries of
series metadata. The package also provides control over API key storage
and caching of API request results.

## Installation

Install the CRAN release of `eia` with

``` r
install.packages("eia")
```

To install the development version from GitHub use

``` r
# install.packages("remotes")
remotes::install_github("ropensci/eia")
```

## Example

To begin, store your API key. You can place it somewhere like your
`.Renviron` file and never have to do anything with the key when you use
the package. You can set it with `eia_set_key` in your R session. You
can always pass it explicitly to the `key` argument of a function.

``` r
library(eia)

# not run
eia_set_key("yourkey") # set API key if not already set globally
```

Load a time series of net electricity generation.

``` r
id <- "ELEC.GEN.ALL-AK-99.A"
(d <- eia_series(id, n = 10))
#> # A tibble: 1 x 13
#>   series_id      name                       units        f     description                         copyright source                iso3166 geography start end   updated       data     
#>   <chr>          <chr>                      <chr>        <chr> <chr>                               <chr>     <chr>                 <chr>   <chr>     <chr> <chr> <chr>         <list>   
#> 1 ELEC.GEN.ALL-~ Net generation : all fuel~ thousand me~ A     "Summation of all fuels used for e~ None      EIA, U.S. Energy Inf~ USA-AK  USA-AK    2001  2019  2020-10-27T1~ <tibble ~

d$data[[1]]
#> # A tibble: 10 x 3
#>    value date        year
#>    <dbl> <date>     <int>
#>  1 6071. 2019-01-01  2019
#>  2 6247. 2018-01-01  2018
#>  3 6497. 2017-01-01  2017
#>  4 6335. 2016-01-01  2016
#>  5 6285. 2015-01-01  2015
#>  6 6043. 2014-01-01  2014
#>  7 6497. 2013-01-01  2013
#>  8 6946. 2012-01-01  2012
#>  9 6871. 2011-01-01  2011
#> 10 6760. 2010-01-01  2010

library(ggplot2)
library(tidyr)
unnest(d, cols = data) %>% ggplot(aes(factor(year), value)) + geom_col() + 
  labs(x = "Year", y = d$units, title = d$name, caption = d$description)
```

<img src="man/figures/README-example-1.png" width="100%" />

## References

See the collection of vignette tutorials and examples as well as
complete package documentation available at the `eia` package
[website](https://docs.ropensci.org/eia/).

-----

Please note that the `eia` project is released with a [Contributor Code
of
Conduct](https://github.com/ropensci/eia/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# eia 0.3.7

* Updated URL for `eia_updates`.
* Minor fixes.
* Updated vignettes.

# eia 0.3.6

* Minor fix to canned report.
* Switch from http to https.
* Documentation updates.

# eia 0.3.5

* Minor documentation updates.

# eia 0.3.4

* Added initial report function.
* Minor code improvements.
* Minor documentation updates.

# eia 0.3.3

* Added a wrapper for the series category endpoint.
* Updated documentation, vignette and unit tests.

# eia 0.3.2

* Updated formatting for CRAN.

# eia 0.3.1

* Updated package metadata and improved Travis testing suite configuration.

# eia 0.3.0

* Added convenient key store methods with getter and setter helpers, optionally making it easy to avoid having to provide the key in every function call.
* Moved `key` argument from first to last among relevant function arguments and updated all examples accordingly.
* Added support for hourly time series requests and date format handling.

# eia 0.2.0

* Added optional memoization to API functions, adding a new `cache` argument.
* Added anti-DOS measures, which can be adjusted using `options()`.
* Added more vignettes and documentation.
* Added helper functions for clearing cached results.
* Added helper functions for working with EIA date strings.
* Added helper functions for time series metadata.
* More output formats and consistency between functions.
* Added unit tests.
* Minor updates to functions, documentation.

# eia 0.1.0

* Added package scaffolding.
* Added initial package functions for working with data categories, time series, and geosets.
* Added function documentation, unit tests, vignettes.
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
## Test environments

* local Windows 10 install, R 4.0.3
* Windows 10 (AppVeyor), R 4.0.3
* Ubuntu 16.04 (Travis CI), R-devel, R-release, R-oldrel
* Mac OSX (Travis CI) R-release
* win-builder (devel and release)
* R-hub (various)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is an update release.
* This update includes a maintainer email address update.

Special note: This package is an API wrapper. The particular API requires users to use their own API key. I cannot run function examples or unit tests on CRAN, but all examples and unit tests run successfully in multiple other environments, on local and remote systems. Full test suite also runs on Travis-CI where I am able to import an encrypted key. API key-dependent vignettes are precompiled for CRAN.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", fig.path = "man/figures/README-", 
  fig.width = 7, fig.height = 4, dev = "CairoPNG", dpi = 150, out.width = "100%",
  message = FALSE, warning = FALSE, error = FALSE
)
library(eia)
```
# eia <img src="man/figures/logo.png" style="margin-left:10px;margin-bottom:5px;" width="120" align="right">
**Author:** [Matthew Leonawicz](https://github.com/leonawicz) <a href="https://orcid.org/0000-0001-9452-2771" target="orcid.widget">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
<br/>
**License:** [MIT](https://opensource.org/licenses/MIT)<br/>

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis build status](https://travis-ci.org/ropensci/eia.svg?branch=master)](https://travis-ci.org/ropensci/eia)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/eia?branch=master&svg=true)](https://ci.appveyor.com/project/leonawicz/eia)
[![Codecov test coverage](https://codecov.io/gh/ropensci/eia/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/eia?branch=master)

[![](https://badges.ropensci.org/342_status.svg)](https://github.com/ropensci/software-review/issues/342)
[![CRAN status](https://www.r-pkg.org/badges/version/eia)](https://cran.r-project.org/package=eia)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/eia)](https://cran.r-project.org/package=eia)
[![Github Stars](https://img.shields.io/github/stars/ropensci/eia.svg?style=social&label=Github)](https://github.com/ropensci/eia/)

The `eia` package provides API access to data from the US [Energy Information Administration](https://www.eia.gov/) (EIA).

Pulling data from the US Energy Information Administration (EIA) API requires a registered API key. A key can be obtained at no cost [here](https://www.eia.gov/opendata/register.php). A valid email and agreement to the API Terms of Service is required to obtain a key.

`eia` includes functions for searching EIA API data categories and importing time series and geoset time series datasets. Datasets returned by these functions are provided in a tidy format or alternatively in more raw form. It also offers helper functions for working with EIA API date strings and time formats and for inspecting different summaries of series metadata. The package also provides control over API key storage and caching of API request results.

## Installation

Install the CRAN release of `eia` with

``` r
install.packages("eia")
```

To install the development version from GitHub use

``` r
# install.packages("remotes")
remotes::install_github("ropensci/eia")
```

## Example

To begin, store your API key. You can place it somewhere like your `.Renviron` file and never have to do anything with the key when you use the package. You can set it with `eia_set_key` in your R session. You can always pass it explicitly to the `key` argument of a function.

```{r example1, eval=FALSE}
library(eia)

# not run
eia_set_key("yourkey") # set API key if not already set globally
```

Load a time series of net electricity generation.

```{r example}
id <- "ELEC.GEN.ALL-AK-99.A"
(d <- eia_series(id, n = 10))

d$data[[1]]

library(ggplot2)
library(tidyr)
unnest(d, cols = data) %>% ggplot(aes(factor(year), value)) + geom_col() + 
  labs(x = "Year", y = d$units, title = d$name, caption = d$description)
```

## References

See the collection of vignette tutorials and examples as well as complete package documentation available at the `eia` package [website](https://docs.ropensci.org/eia/).

---

Please note that the `eia` project is released with a [Contributor Code of Conduct](https://github.com/ropensci/eia/blob/master/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "EIA time series data"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{EIA time series data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Finding the series ID

Like category information, time series data is obtained based on its ID. A complete example includes finding the ID for the series if you do not already know it. Chances are you may already know the series IDs you need after using the API explorer on the EIA website. Below, look for total electricity consumption.


```r
library(eia)
# eia_set_key("yourkey") # set API key if not already set globally
eia_cats()
#> $category
#> # A tibble: 1 x 3
#>   category_id name          notes
#>   <chr>       <chr>         <chr>
#> 1 371         EIA Data Sets ""   
#> 
#> $childcategories
#> # A tibble: 14 x 2
#>    category_id name                               
#>          <int> <chr>                              
#>  1           0 Electricity                        
#>  2       40203 State Energy Data System (SEDS)    
#>  3      714755 Petroleum                          
#>  4      714804 Natural Gas                        
#>  5      711224 Total Energy                       
#>  6      717234 Coal                               
#>  7      829714 Short-Term Energy Outlook          
#>  8      964164 Annual Energy Outlook              
#>  9     1292190 Crude Oil Imports                  
#> 10     2123635 U.S. Electric System Operating Data
#> 11     2134384 International Energy Data          
#> 12     2251604 CO2 Emissions                      
#> 13     2631064 International Energy Outlook       
#> 14     2889994 U.S. Nuclear Outages
```

Electricity has category ID 0. Take a closer look there.


```r
eia_cats(0)
#> $category
#> # A tibble: 1 x 4
#>   category_id parent_category_id name        notes
#>   <chr>       <chr>              <chr>       <chr>
#> 1 0           371                Electricity ""   
#> 
#> $childcategories
#> # A tibble: 19 x 2
#>    category_id name                                                              
#>          <int> <chr>                                                             
#>  1           1 Net generation                                                    
#>  2          35 Total consumption                                                 
#>  3          32 Total consumption (Btu)                                           
#>  4          36 Consumption for electricity generation                            
#>  5          33 Consumption for electricity generation (Btu)                      
#>  6          37 Consumption for useful thermal output                             
#>  7          34 Consumption for useful thermal output (Btu)                       
#>  8        1017 Plant level data                                                  
#>  9          38 Retail sales of electricity                                       
#> 10          39 Revenue from retail sales of electricity                          
#> 11          40 Average retail price of electricity                               
#> 12     1718389 Number of customer accounts                                       
#> 13       41137 Fossil-fuel stocks for electricity generation                     
#> 14       41138 Receipts of fossil fuels by electricity plants                    
#> 15       41139 Receipts of fossil fuels by electricity plants (Btu)              
#> 16       41140 Average cost of fossil fuels for electricity generation           
#> 17       41141 Average cost of fossil fuels for electricity generation (per Btu) 
#> 18       41142 Quality of fossil fuels in electricity generation : sulfur content
#> 19       41143 Quality of fossil fuels in electricity generation : ash content
```

There are two categories referring to total consumption. Take category ID 32 as an example and step deeper into the category hierarchy.


```r
eia_cats(32)
#> $category
#> # A tibble: 1 x 4
#>   category_id parent_category_id name                    notes
#>   <chr>       <chr>              <chr>                   <chr>
#> 1 32          0                  Total consumption (Btu) ""   
#> 
#> $childcategories
#> # A tibble: 2 x 2
#>   category_id name        
#>         <int> <chr>       
#> 1         373 By fuel type
#> 2         372 By sector
```

At this point you have a choice between total consumption by sector or fuel type. Select by sector.


```r
eia_cats(372)
#> $category
#> # A tibble: 1 x 4
#>   category_id parent_category_id name      notes
#>   <chr>       <chr>              <chr>     <chr>
#> 1 372         32                 By sector ""   
#> 
#> $childcategories
#> # A tibble: 11 x 2
#>    category_id name                               
#>          <int> <chr>                              
#>  1         388 Electric power (total)             
#>  2         389 Electric utility                   
#>  3         390 Independent power producers (total)
#>  4         391 Electric utility non-cogen         
#>  5         392 Electric utility cogen             
#>  6         393 All commercial (total)             
#>  7         394 Commercial non-cogen               
#>  8         395 Commercial cogen                   
#>  9         396 All industrial (total)             
#> 10         397 Industrial non-cogen               
#> 11         398 Industrial cogen
```

Then total electrical power.


```r
eia_cats(388)
#> $category
#> # A tibble: 1 x 4
#>   category_id parent_category_id name                   notes
#>   <chr>       <chr>              <chr>                  <chr>
#> 1 388         372                Electric power (total) ""   
#> 
#> $childcategories
#> # A tibble: 4 x 2
#>   category_id name             
#>         <int> <chr>            
#> 1         738 Coal             
#> 2         739 Petroleum liquids
#> 3         740 Petroleum coke   
#> 4         741 Natural gas
```

And finally coal.


```r
(x <- eia_cats(738))
#> $category
#> # A tibble: 1 x 4
#>   category_id parent_category_id name  notes
#>   <chr>       <chr>              <chr> <chr>
#> 1 738         388                Coal  ""   
#> 
#> $childseries
#> # A tibble: 174 x 5
#>    series_id                     name                                                                           f     units         updated              
#>    <chr>                         <chr>                                                                          <chr> <chr>         <chr>                
#>  1 ELEC.CONS_TOT_BTU.COW-AK-98.A Total consumption (Btu) : coal : Alaska : electric power (total) : annual      A     million MMBtu 27-OCT-20 06.46.54 PM
#>  2 ELEC.CONS_TOT_BTU.COW-AK-98.M Total consumption (Btu) : coal : Alaska : electric power (total) : monthly     M     million MMBtu 27-OCT-20 06.46.54 PM
#>  3 ELEC.CONS_TOT_BTU.COW-AK-98.Q Total consumption (Btu) : coal : Alaska : electric power (total) : quarterly   Q     million MMBtu 27-OCT-20 06.46.54 PM
#>  4 ELEC.CONS_TOT_BTU.COW-AL-98.A Total consumption (Btu) : coal : Alabama : electric power (total) : annual     A     million MMBtu 27-OCT-20 06.46.54 PM
#>  5 ELEC.CONS_TOT_BTU.COW-AL-98.M Total consumption (Btu) : coal : Alabama : electric power (total) : monthly    M     million MMBtu 27-OCT-20 06.46.54 PM
#>  6 ELEC.CONS_TOT_BTU.COW-AL-98.Q Total consumption (Btu) : coal : Alabama : electric power (total) : quarterly  Q     million MMBtu 27-OCT-20 06.46.54 PM
#>  7 ELEC.CONS_TOT_BTU.COW-AR-98.A Total consumption (Btu) : coal : Arkansas : electric power (total) : annual    A     million MMBtu 27-OCT-20 06.46.54 PM
#>  8 ELEC.CONS_TOT_BTU.COW-AR-98.M Total consumption (Btu) : coal : Arkansas : electric power (total) : monthly   M     million MMBtu 27-OCT-20 06.46.54 PM
#>  9 ELEC.CONS_TOT_BTU.COW-AR-98.Q Total consumption (Btu) : coal : Arkansas : electric power (total) : quarterly Q     million MMBtu 27-OCT-20 06.46.54 PM
#> 10 ELEC.CONS_TOT_BTU.COW-AZ-98.A Total consumption (Btu) : coal : Arizona : electric power (total) : annual     A     million MMBtu 27-OCT-20 06.46.54 PM
#> # ... with 164 more rows
```

## Available child series

At this point you have reached a terminal node of the category tree. Instead of another table of child category IDs and names in the result, there is a `childseries` table.

This table contains:

*    time series IDs
*    names that describe the nested position of the data in the overall category hierarchy
*    time format
*    units
*    the time stamp of the most recent data update

Each row in this table represents a unique time series dataset; in this case for different states and in annual, quarterly and monthly time steps. TO obtain the time series data, make a request using `eia_series` and provide a series ID.

## Time formats

To see how the different time formats are parsed, take the first three IDs for Alaska. Request only the three most recent results for each series.


```r
id <- x$childseries$series_id[1:3]
x1 <- eia_series(id[1], n = 3)
x2 <- eia_series(id[2], n = 3)
x3 <- eia_series(id[3], n = 3)
```

The format of each result is the same. Inspect the first one. It is a data frame with one row. All but the final column, `data`, give metadata about the series. `data` is a list column (in this case of length one) that can be extracted directly or unnested using `tidyr::unnest`.


```r
library(dplyr)
library(tidyr)
library(ggplot2)

x1$data[[1]]
#> # A tibble: 3 x 3
#>   value date        year
#>   <dbl> <date>     <int>
#> 1 11.0  2019-01-01  2019
#> 2 10.4  2018-01-01  2018
#> 3  9.21 2017-01-01  2017

select(x1, series_id, data) %>% unnest(cols = data)
#> # A tibble: 3 x 4
#>   series_id                     value date        year
#>   <chr>                         <dbl> <date>     <int>
#> 1 ELEC.CONS_TOT_BTU.COW-AK-98.A 11.0  2019-01-01  2019
#> 2 ELEC.CONS_TOT_BTU.COW-AK-98.A 10.4  2018-01-01  2018
#> 3 ELEC.CONS_TOT_BTU.COW-AK-98.A  9.21 2017-01-01  2017

unnest(x1, cols = data) %>%
  ggplot(aes(factor(year), value)) + geom_col() +
  labs(x = "Year", y = x1$units, title = "Total Alaska electricity consumption",
       subtitle = "From coal sources", caption = x1$description)
```

<img src="series-series2-1.png" title="plot of chunk series2" alt="plot of chunk series2" width="100%" />

Results are similarly structured for the other series, but the columns containing date information differ.


```r
x1$data[[1]]
#> # A tibble: 3 x 3
#>   value date        year
#>   <dbl> <date>     <int>
#> 1 11.0  2019-01-01  2019
#> 2 10.4  2018-01-01  2018
#> 3  9.21 2017-01-01  2017
x2$data[[1]]
#> # A tibble: 3 x 4
#>    value date        year month
#>    <dbl> <date>     <int> <int>
#> 1 NA     2020-08-01  2020     8
#> 2 NA     2020-07-01  2020     7
#> 3  0.770 2020-06-01  2020     6
x3$data[[1]]
#> # A tibble: 3 x 4
#>   value date        year   qtr
#>   <dbl> <date>     <int> <int>
#> 1 NA    2020-04-01  2020     2
#> 2 NA    2020-01-01  2020     1
#> 3  3.05 2019-10-01  2019     4
```

## Multiple series

The EIA API allows multiple series to be requested in a single API call. You should do this whenever possible to reduce the number of requests you make. To request the same data as above, just provide the `id` vector. Other arguments like `n` are not vectorized.


```r
x <- eia_series(id, n = 3)
x
#> # A tibble: 3 x 13
#>   series_id       name                      units    f     description                            copyright source           iso3166 geography start end   updated     data    
#>   <chr>           <chr>                     <chr>    <chr> <chr>                                  <chr>     <chr>            <chr>   <chr>     <chr> <chr> <chr>       <list>  
#> 1 ELEC.CONS_TOT_~ Total consumption (Btu) ~ million~ A     "Summation of all types of coal; Powe~ None      EIA, U.S. Energ~ USA-AK  USA-AK    2001  2019  2020-10-27~ <tibble~
#> 2 ELEC.CONS_TOT_~ Total consumption (Btu) ~ million~ M     "Summation of all types of coal; Powe~ None      EIA, U.S. Energ~ USA-AK  USA-AK    2001~ 2020~ 2020-10-27~ <tibble~
#> 3 ELEC.CONS_TOT_~ Total consumption (Btu) ~ million~ Q     "Summation of all types of coal; Powe~ None      EIA, U.S. Energ~ USA-AK  USA-AK    2001~ 2020~ 2020-10-27~ <tibble~
```

There are now three rows in the table containing the same data as before. The `data` list column also contains the same structures as before.


```r
x$data
#> [[1]]
#> # A tibble: 3 x 3
#>   value date        year
#>   <dbl> <date>     <int>
#> 1 11.0  2019-01-01  2019
#> 2 10.4  2018-01-01  2018
#> 3  9.21 2017-01-01  2017
#> 
#> [[2]]
#> # A tibble: 3 x 4
#>    value date        year month
#>    <dbl> <date>     <int> <int>
#> 1 NA     2020-08-01  2020     8
#> 2 NA     2020-07-01  2020     7
#> 3  0.770 2020-06-01  2020     6
#> 
#> [[3]]
#> # A tibble: 3 x 4
#>   value date        year   qtr
#>   <dbl> <date>     <int> <int>
#> 1 NA    2020-04-01  2020     2
#> 2 NA    2020-01-01  2020     1
#> 3  3.05 2019-10-01  2019     4
```

These can be unnested and filled in with `NA` as needed.


```r
select(x, series_id, data) %>% unnest(cols = data)
#> # A tibble: 9 x 6
#>   series_id                      value date        year month   qtr
#>   <chr>                          <dbl> <date>     <int> <int> <int>
#> 1 ELEC.CONS_TOT_BTU.COW-AK-98.A 11.0   2019-01-01  2019    NA    NA
#> 2 ELEC.CONS_TOT_BTU.COW-AK-98.A 10.4   2018-01-01  2018    NA    NA
#> 3 ELEC.CONS_TOT_BTU.COW-AK-98.A  9.21  2017-01-01  2017    NA    NA
#> 4 ELEC.CONS_TOT_BTU.COW-AK-98.M NA     2020-08-01  2020     8    NA
#> 5 ELEC.CONS_TOT_BTU.COW-AK-98.M NA     2020-07-01  2020     7    NA
#> 6 ELEC.CONS_TOT_BTU.COW-AK-98.M  0.770 2020-06-01  2020     6    NA
#> 7 ELEC.CONS_TOT_BTU.COW-AK-98.Q NA     2020-04-01  2020    NA     2
#> 8 ELEC.CONS_TOT_BTU.COW-AK-98.Q NA     2020-01-01  2020    NA     1
#> 9 ELEC.CONS_TOT_BTU.COW-AK-98.Q  3.05  2019-10-01  2019    NA     4
```

## Time period and number of results

Here are some things to keep in mind about `eia_series` arguments.

* The function only makes one API call, combining multiple series IDs if provided.
* If there are fewer results than you request, everything available is returned.
* Results are from the most recent unless using `start` alone.
* If you provide `start` and `end`, `n` is ignored.
* If you do not provide a closed period or `n`, you will receive all relevant data available, subject to any API limits.
* Depending on the nature of the request, you may have to construct more than one call.


```r
eia_series(id, n = 10)
eia_series(id, end = 2016, n = 5)
eia_series(id, start = 2000, end = 2016)
```

## Output format

The default is to return tidy data in a tibble data frame. You can set `tidy = FALSE` to return the list returned by `jsonlite::fromJSON` without any further processing.


```r
eia_series(id, n = 3, tidy = FALSE)
#> $request
#> $request$command
#> [1] "series"
#> 
#> $request$series_id
#> [1] "ELEC.CONS_TOT_BTU.COW-AK-98.A;ELEC.CONS_TOT_BTU.COW-AK-98.M;ELEC.CONS_TOT_BTU.COW-AK-98.Q"
#> 
#> 
#> $series
#>                       series_id                                                                         name         units f
#> 1 ELEC.CONS_TOT_BTU.COW-AK-98.A    Total consumption (Btu) : coal : Alaska : electric power (total) : annual million MMBtu A
#> 2 ELEC.CONS_TOT_BTU.COW-AK-98.M   Total consumption (Btu) : coal : Alaska : electric power (total) : monthly million MMBtu M
#> 3 ELEC.CONS_TOT_BTU.COW-AK-98.Q Total consumption (Btu) : coal : Alaska : electric power (total) : quarterly million MMBtu Q
#>                                                                                                   description copyright                                      source iso3166
#> 1 Summation of all types of coal; Power plants owned by companies whose primary purpose is to produce power;       None EIA, U.S. Energy Information Administration  USA-AK
#> 2 Summation of all types of coal; Power plants owned by companies whose primary purpose is to produce power;       None EIA, U.S. Energy Information Administration  USA-AK
#> 3 Summation of all types of coal; Power plants owned by companies whose primary purpose is to produce power;       None EIA, U.S. Energy Information Administration  USA-AK
#>   geography  start    end                  updated                                          data
#> 1    USA-AK   2001   2019 2020-10-27T18:46:54-0400 2019, 2018, 2017, 11.03525, 10.37727, 9.21412
#> 2    USA-AK 200101 202008 2020-10-27T18:46:54-0400       202008, 202007, 202006, NA, NA, 0.77003
#> 3    USA-AK 2001Q1 2020Q2 2020-10-27T18:46:54-0400       2020Q2, 2020Q1, 2019Q4, NA, NA, 3.04527
```

You can also return the raw JSON data in a character string if you need to process this directly with other code.


```r
cat(eia_series(id, n = 3, tidy = NA))
#> {"request":{"command":"series","series_id":"ELEC.CONS_TOT_BTU.COW-AK-98.A;ELEC.CONS_TOT_BTU.COW-AK-98.M;ELEC.CONS_TOT_BTU.COW-AK-98.Q"},"series":[{"series_id":"ELEC.CONS_TOT_BTU.COW-AK-98.A","name":"Total consumption (Btu) : coal : Alaska : electric power (total) : annual","units":"million MMBtu","f":"A","description":"Summation of all types of coal; Power plants owned by companies whose primary purpose is to produce power; ","copyright":"None","source":"EIA, U.S. Energy Information Administration","iso3166":"USA-AK","geography":"USA-AK","start":"2001","end":"2019","updated":"2020-10-27T18:46:54-0400","data":[["2019",11.03525],["2018",10.37727],["2017",9.21412]]},{"series_id":"ELEC.CONS_TOT_BTU.COW-AK-98.M","name":"Total consumption (Btu) : coal : Alaska : electric power (total) : monthly","units":"million MMBtu","f":"M","description":"Summation of all types of coal; Power plants owned by companies whose primary purpose is to produce power; ","copyright":"None","source":"EIA, U.S. Energy Information Administration","iso3166":"USA-AK","geography":"USA-AK","start":"200101","end":"202008","updated":"2020-10-27T18:46:54-0400","data":[["202008",null],["202007",null],["202006",0.77003]]},{"series_id":"ELEC.CONS_TOT_BTU.COW-AK-98.Q","name":"Total consumption (Btu) : coal : Alaska : electric power (total) : quarterly","units":"million MMBtu","f":"Q","description":"Summation of all types of coal; Power plants owned by companies whose primary purpose is to produce power; ","copyright":"None","source":"EIA, U.S. Energy Information Administration","iso3166":"USA-AK","geography":"USA-AK","start":"2001Q1","end":"2020Q2","updated":"2020-10-27T18:46:54-0400","data":[["2020Q2",null],["2020Q1",null],["2019Q4",3.04527]]}]}
```

This allows you to use the returned results with existing code you may have that requires data in one of these less processed structures.

## Helpers functions

### Time series metadata

There are some functions available that make small API calls and return only metadata associated with a time series dataset.


```r
eia_series_metadata(id)
#> # A tibble: 3 x 12
#>   series_id        name                         units    f     description                               copyright source            iso3166 geography start end   updated     
#>   <chr>            <chr>                        <chr>    <chr> <chr>                                     <chr>     <chr>             <chr>   <chr>     <chr> <chr> <chr>       
#> 1 ELEC.CONS_TOT_B~ Total consumption (Btu) : c~ million~ A     "Summation of all types of coal; Power p~ None      EIA, U.S. Energy~ USA-AK  USA-AK    2001  2019  2020-10-27T~
#> 2 ELEC.CONS_TOT_B~ Total consumption (Btu) : c~ million~ M     "Summation of all types of coal; Power p~ None      EIA, U.S. Energy~ USA-AK  USA-AK    2001~ 2020~ 2020-10-27T~
#> 3 ELEC.CONS_TOT_B~ Total consumption (Btu) : c~ million~ Q     "Summation of all types of coal; Power p~ None      EIA, U.S. Energy~ USA-AK  USA-AK    2001~ 2020~ 2020-10-27T~
eia_series_updates(id)
#> # A tibble: 3 x 2
#>   series_id                     updated                 
#>   <chr>                         <chr>                   
#> 1 ELEC.CONS_TOT_BTU.COW-AK-98.A 2020-10-27T18:46:54-0400
#> 2 ELEC.CONS_TOT_BTU.COW-AK-98.M 2020-10-27T18:46:54-0400
#> 3 ELEC.CONS_TOT_BTU.COW-AK-98.Q 2020-10-27T18:46:54-0400
eia_series_dates(id)
#> # A tibble: 333 x 4
#>    series_id                     date       eiadate date_format
#>    <chr>                         <date>     <chr>   <chr>      
#>  1 ELEC.CONS_TOT_BTU.COW-AK-98.A 2001-01-01 2001    A          
#>  2 ELEC.CONS_TOT_BTU.COW-AK-98.A 2002-01-01 2002    A          
#>  3 ELEC.CONS_TOT_BTU.COW-AK-98.A 2003-01-01 2003    A          
#>  4 ELEC.CONS_TOT_BTU.COW-AK-98.A 2004-01-01 2004    A          
#>  5 ELEC.CONS_TOT_BTU.COW-AK-98.A 2005-01-01 2005    A          
#>  6 ELEC.CONS_TOT_BTU.COW-AK-98.A 2006-01-01 2006    A          
#>  7 ELEC.CONS_TOT_BTU.COW-AK-98.A 2007-01-01 2007    A          
#>  8 ELEC.CONS_TOT_BTU.COW-AK-98.A 2008-01-01 2008    A          
#>  9 ELEC.CONS_TOT_BTU.COW-AK-98.A 2009-01-01 2009    A          
#> 10 ELEC.CONS_TOT_BTU.COW-AK-98.A 2010-01-01 2010    A          
#> # ... with 323 more rows
eia_series_range(id)
#> # A tibble: 3 x 7
#>   series_id                     start_date end_date   start  end    date_format     n
#>   <chr>                         <date>     <date>     <chr>  <chr>  <chr>       <int>
#> 1 ELEC.CONS_TOT_BTU.COW-AK-98.A 2001-01-01 2019-01-01 2001   2019   A              19
#> 2 ELEC.CONS_TOT_BTU.COW-AK-98.M 2001-01-01 2020-08-01 200101 202008 M             236
#> 3 ELEC.CONS_TOT_BTU.COW-AK-98.Q 2001-01-01 2020-04-01 2001Q1 2020Q2 Q              78
eia_series_cats(id)
#> # A tibble: 6 x 3
#>   series_id                     category_id name                  
#>   <chr>                               <int> <chr>                 
#> 1 ELEC.CONS_TOT_BTU.COW-AK-98.A         474 Electric power (total)
#> 2 ELEC.CONS_TOT_BTU.COW-AK-98.A         738 Coal                  
#> 3 ELEC.CONS_TOT_BTU.COW-AK-98.M         474 Electric power (total)
#> 4 ELEC.CONS_TOT_BTU.COW-AK-98.M         738 Coal                  
#> 5 ELEC.CONS_TOT_BTU.COW-AK-98.Q         474 Electric power (total)
#> 6 ELEC.CONS_TOT_BTU.COW-AK-98.Q         738 Coal
```

Like `eia_seires`, these functions accept an `id` vector. They always return a tibble data frame. `eia_series_cats` uses the `series categories` endpoint and accepts the `tidy` argument so that output from the endpoint may be a JSON string, list or the tibble data frame.

### EIA date strings

EIA date strings used to specify start and end dates for time series requests are character strings that are not in any standard date formats. There are several functions that assist with moving between these strings and standard dates.

You can convert EIA date strings to dates.


```r
eiadate_to_date(c("201803", "201804"))
#> [1] "2018-03-01" "2018-04-01"
```

or dates to EIA format; here are examples using annual, quarterly and monthly time formats.


```r
date_to_eiadate("2018-05-14", "A")
#> [1] "2018"
date_to_eiadate("2018-05-14", "Q")
#> [1] "2018Q2"
date_to_eiadate("2018-05-14", "M")
#> [1] "201805"
```

It is also easy to create a date sequence from two EIA time stamps. The format is parsed from the first value (they are intended to always be consistent).


```r
(x <- eiadate_to_date_seq("2018Q1", "2018Q4"))
#> [1] "2018-01-01" "2018-04-01" "2018-07-01" "2018-10-01"
date_to_eiadate(x)
#> [1] "2018" "2018" "2018" "2018"
```

## Checking for data updates

It is good practice to minimize the number of API calls you make wherever possible. One way to do this is to not request data that has not changed since you last requested it. You can make an API call to the EIA `updates` endpoint to check update times on data series. If any series have not been updated since you last obtained the data, then you know you do not need to request the data again.

The `eia_series_updates` function shown above is handy for checking the most recent data update times of a specific set of series IDs. However, it must make one or more API calls to do so. For general checks on a potentially large number of series without having to query them all, you should use `eia_updates`. This function takes a category ID and can return the last update times for all series under that category.

If a category level has no series directly associated with it, an empty data frame is returned.


```r
eia_updates(389)
#> # A tibble: 0 x 2
#> # ... with 2 variables: series_id <chr>, updated <chr>
```

If a category has series available, paginated results are returned in a data frame. You can use `n` and `start` (together, unlike for general data requests with `eia_series`) to indicate how many rows to return and where to start. This helps you to cycle through pages of results. By default, `n = 50` and `start = 1`. The EIA API `updates` endpoint allows a maximum of `n = 10000`.


```r
eia_updates(742, n = 5)
#> # A tibble: 5 x 2
#>   series_id                    updated                 
#>   <chr>                        <chr>                   
#> 1 ELEC.CONS_TOT_BTU.COW-AK-1.A 2020-10-27T18:46:54-0400
#> 2 ELEC.CONS_TOT_BTU.COW-AK-1.M 2020-10-27T18:46:54-0400
#> 3 ELEC.CONS_TOT_BTU.COW-AK-1.Q 2020-10-27T18:46:54-0400
#> 4 ELEC.CONS_TOT_BTU.COW-AL-1.A 2020-10-27T18:46:54-0400
#> 5 ELEC.CONS_TOT_BTU.COW-AL-1.M 2020-10-27T18:46:54-0400
```

Set `deep = TRUE` to obtain series associated with child categories. Category 389 above did not have series, but some of the child categories do.


```r
eia_updates(389, n = 5, deep = TRUE)
#> # A tibble: 5 x 2
#>   series_id                    updated                 
#>   <chr>                        <chr>                   
#> 1 ELEC.CONS_TOT_BTU.COW-AK-1.A 2020-10-27T18:46:54-0400
#> 2 ELEC.CONS_TOT_BTU.COW-AK-1.M 2020-10-27T18:46:54-0400
#> 3 ELEC.CONS_TOT_BTU.COW-AK-1.Q 2020-10-27T18:46:54-0400
#> 4 ELEC.CONS_TOT_BTU.COW-AL-1.A 2020-10-27T18:46:54-0400
#> 5 ELEC.CONS_TOT_BTU.COW-AL-1.M 2020-10-27T18:46:54-0400
```

The above example works to show this while using `n = 5` because there were no results for the top parent category. if there had been, it would likely be necessary to request more results in order to see that child series were included.
---
title: "Package overview"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Package overview}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



This vignette provides a brief overview of the most important functions in `eia`. Other vignettes go into greater depth on specific topics and API endpoints.

## API key

### Register a key with EIA

Obtaining an API key is easy and free.

Pulling data from the US Energy Information Administration (EIA) API requires a registered API key. A key can be obtained at no cost [here](https://www.eia.gov/opendata/register.php). A valid email and agreement to the API Terms of Service is required to obtain a key.

It is important to store your API key somewhere secure. Do not commit it to a repository or otherwise share it. For example, you could store it in your `.Renviron` file.

### Key storage and retrieval

You can always provide the `key` argument to every API function call, but you do not have to. There are getter and setter helpers available to make using `eia` functions a more seamless experience.

`eia_set_key` gives you the option of storing your key for the duration of your R session.


```r
library(eia)
# eia_set_key("yourkey")
# eia_get_key() # retrieve it
```

If the key already exists in the system environment and you plan to pass `key` to functions explicitly, you could start as follows.


```r
key <- Sys.getenv("EIA_KEY")

# or:
key <- eia_get_key()
```

In general, however, if your key is set globally such as in `.Renviron`, you do not need to do anything regarding the key when you use the package. See the vignette on API details for more information about all the options you have for key storage.

## EIA categories

Once you have your EIA registered API key and have it in place for your R session by whichever method you prefer, you are ready to begin accessing data from the EIA API.

It is helpful to be aware of various categories of data that are available. Each category has a unique ID number that is required to access associated data. The first call below does not include an ID. The result is a list of two data frames. The first is metadata associated with that position in the category hierarchy. The second is the child category information.

Here is the top-level category information.


```r
eia_cats()
#> $category
#> # A tibble: 1 x 3
#>   category_id name          notes
#>   <chr>       <chr>         <chr>
#> 1 371         EIA Data Sets ""   
#> 
#> $childcategories
#> # A tibble: 14 x 2
#>    category_id name                               
#>          <int> <chr>                              
#>  1           0 Electricity                        
#>  2       40203 State Energy Data System (SEDS)    
#>  3      714755 Petroleum                          
#>  4      714804 Natural Gas                        
#>  5      711224 Total Energy                       
#>  6      717234 Coal                               
#>  7      829714 Short-Term Energy Outlook          
#>  8      964164 Annual Energy Outlook              
#>  9     1292190 Crude Oil Imports                  
#> 10     2123635 U.S. Electric System Operating Data
#> 11     2134384 International Energy Data          
#> 12     2251604 CO2 Emissions                      
#> 13     2631064 International Energy Outlook       
#> 14     2889994 U.S. Nuclear Outages
```

The child category IDs can be used to query results from that category.


```r
eia_cats(0)
#> $category
#> # A tibble: 1 x 4
#>   category_id parent_category_id name        notes
#>   <chr>       <chr>              <chr>       <chr>
#> 1 0           371                Electricity ""   
#> 
#> $childcategories
#> # A tibble: 19 x 2
#>    category_id name                                                              
#>          <int> <chr>                                                             
#>  1           1 Net generation                                                    
#>  2          35 Total consumption                                                 
#>  3          32 Total consumption (Btu)                                           
#>  4          36 Consumption for electricity generation                            
#>  5          33 Consumption for electricity generation (Btu)                      
#>  6          37 Consumption for useful thermal output                             
#>  7          34 Consumption for useful thermal output (Btu)                       
#>  8        1017 Plant level data                                                  
#>  9          38 Retail sales of electricity                                       
#> 10          39 Revenue from retail sales of electricity                          
#> 11          40 Average retail price of electricity                               
#> 12     1718389 Number of customer accounts                                       
#> 13       41137 Fossil-fuel stocks for electricity generation                     
#> 14       41138 Receipts of fossil fuels by electricity plants                    
#> 15       41139 Receipts of fossil fuels by electricity plants (Btu)              
#> 16       41140 Average cost of fossil fuels for electricity generation           
#> 17       41141 Average cost of fossil fuels for electricity generation (per Btu) 
#> 18       41142 Quality of fossil fuels in electricity generation : sulfur content
#> 19       41143 Quality of fossil fuels in electricity generation : ash content
```

## EIA time series data

Time series data is obtained by series ID. Most columns contain metadata. The `data` column contains the time series data.


```r
library(dplyr)
library(tidyr)
library(ggplot2)

id <- "ELEC.GEN.ALL-AK-99.A"
(x <- eia_series(id, start = 2010, end = 2019))
#> # A tibble: 1 x 13
#>   series_id     name                     units       f     description                       copyright source              iso3166 geography start end   updated      data     
#>   <chr>         <chr>                    <chr>       <chr> <chr>                             <chr>     <chr>               <chr>   <chr>     <chr> <chr> <chr>        <list>   
#> 1 ELEC.GEN.ALL~ Net generation : all fu~ thousand m~ A     "Summation of all fuels used for~ None      EIA, U.S. Energy I~ USA-AK  USA-AK    2001  2019  2020-10-27T~ <tibble ~

x$data[[1]]
#> # A tibble: 10 x 3
#>    value date        year
#>    <dbl> <date>     <int>
#>  1 6071. 2019-01-01  2019
#>  2 6247. 2018-01-01  2018
#>  3 6497. 2017-01-01  2017
#>  4 6335. 2016-01-01  2016
#>  5 6285. 2015-01-01  2015
#>  6 6043. 2014-01-01  2014
#>  7 6497. 2013-01-01  2013
#>  8 6946. 2012-01-01  2012
#>  9 6871. 2011-01-01  2011
#> 10 6760. 2010-01-01  2010

select(x, units, data) %>% unnest(cols = data) %>%
  ggplot(aes(factor(year), value)) + geom_col() +
  labs(x = "Year", y = x$units[1], title = "Net Alaska electricity generation",
       subtitle = "From all fuels", caption = x$description[1])
```

<img src="eia-series1-1.png" title="plot of chunk series1" alt="plot of chunk series1" width="100%" />

You can provide arguments like the following:


```r
eia_series(id) # max results
eia_series(id, n = 5) # most recent five
eia_series(id, end = 2016, n = 5) # ending in 2016
eia_series(id, start = 2000, end = 2016) # specific period
```

As with `eia_cats`, the output format does not need to be tidy:


```r
eia_series(id, n = 5, tidy = FALSE) # results of jsonlite::fromJSON
eia_series(id, n = 5, tidy = NA) # origina JSON as character string
```

This allows you to use the returned results with existing code you may have that requires data in one of these less processed structures.

## EIA geosets

Geosets are metadata structures organizing time series datasets that can be mapped. Arguments to `eia_geoset` are the same as `eia_series` with the addition of `region`. Like `id`, `region` can be a vector. Most of the details are the same as before.

In the example below using total electricity generation, get the last two data points for each of and two US states. `dplyr` and `tidyr` are used here to clean up the result a bit for purposes of display. `gpplot2` is used to graph the data after it has been unnested for each state.


```r
id <- c("ELEC.GEN.ALL-99.M") # monthly
region <- c("USA-CA", "USA-NY")
(x <- eia_geoset(id, region, start = "201801", end = "201812"))
#> # A tibble: 2 x 11
#>   geoset_id       setname                             f     units            series_id        name                                       region latlon start end   data        
#>   <chr>           <chr>                               <chr> <chr>            <chr>            <chr>                                      <chr>  <chr>  <chr> <chr> <list>      
#> 1 ELEC.GEN.ALL-9~ Net generation : all fuels : all s~ M     thousand megawa~ ELEC.GEN.ALL-CA~ Net generation : all fuels : California :~ USA-CA <NA>   2001~ 2020~ <tibble [12~
#> 2 ELEC.GEN.ALL-9~ Net generation : all fuels : all s~ M     thousand megawa~ ELEC.GEN.ALL-NY~ Net generation : all fuels : New York : a~ USA-NY <NA>   2001~ 2020~ <tibble [12~

select(x, region, data) %>% unnest(cols = data)
#> # A tibble: 24 x 5
#>    region  value date        year month
#>    <chr>   <dbl> <date>     <int> <int>
#>  1 USA-CA 14364. 2018-12-01  2018    12
#>  2 USA-CA 14605. 2018-11-01  2018    11
#>  3 USA-CA 16864. 2018-10-01  2018    10
#>  4 USA-CA 17592. 2018-09-01  2018     9
#>  5 USA-CA 20994. 2018-08-01  2018     8
#>  6 USA-CA 22356. 2018-07-01  2018     7
#>  7 USA-CA 17514. 2018-06-01  2018     6
#>  8 USA-CA 15932. 2018-05-01  2018     5
#>  9 USA-CA 14811. 2018-04-01  2018     4
#> 10 USA-CA 14039. 2018-03-01  2018     3
#> # ... with 14 more rows

unnest(x, cols = data) %>%
  ggplot(aes(date, value, color = region)) + geom_line() +
  labs(y = x$units[1], title = "Net electricity generation", subtitle = "From all fuels")
```

<img src="eia-geoset1-1.png" title="plot of chunk geoset1" alt="plot of chunk geoset1" width="100%" />

Another convenience of `eia_geoset` is the ability to provide regions in the following forms.

* 2-character US state abbreviations
* State names
* US Census region names
* US census division names

These shortcuts make it easier to construct an API call involving several states.


```r
region <- c("AK", "New England")
x <- eia_geoset(id, region, n = 2)
select(x, region, data) %>% unnest(cols = data)
#> # A tibble: 14 x 5
#>    region value date        year month
#>    <chr>  <dbl> <date>     <int> <int>
#>  1 USA-AK  534. 2020-08-01  2020     8
#>  2 USA-AK  626. 2020-07-01  2020     7
#>  3 USA-CT 4001. 2020-08-01  2020     8
#>  4 USA-CT 4443. 2020-07-01  2020     7
#>  5 USA-MA 2200. 2020-08-01  2020     8
#>  6 USA-MA 2754. 2020-07-01  2020     7
#>  7 USA-ME  861. 2020-08-01  2020     8
#>  8 USA-ME  962. 2020-07-01  2020     7
#>  9 USA-NH 1632. 2020-08-01  2020     8
#> 10 USA-NH 1794. 2020-07-01  2020     7
#> 11 USA-RI  889. 2020-08-01  2020     8
#> 12 USA-RI  845. 2020-07-01  2020     7
#> 13 USA-VT  194. 2020-08-01  2020     8
#> 14 USA-VT  203. 2020-07-01  2020     7

region <- "Middle Atlantic"
x <- eia_geoset(id, region, n = 12)
select(x, region, data) %>% unnest(cols = data)
#> # A tibble: 36 x 5
#>    region value date        year month
#>    <chr>  <dbl> <date>     <int> <int>
#>  1 USA-NJ 6639. 2020-08-01  2020     8
#>  2 USA-NJ 7385. 2020-07-01  2020     7
#>  3 USA-NJ 5721. 2020-06-01  2020     6
#>  4 USA-NJ 4169. 2020-05-01  2020     5
#>  5 USA-NJ 4334. 2020-04-01  2020     4
#>  6 USA-NJ 4469. 2020-03-01  2020     3
#>  7 USA-NJ 4912. 2020-02-01  2020     2
#>  8 USA-NJ 5488. 2020-01-01  2020     1
#>  9 USA-NJ 6193. 2019-12-01  2019    12
#> 10 USA-NJ 5510. 2019-11-01  2019    11
#> # ... with 26 more rows

unnest(x, cols = data) %>%
  ggplot(aes(date, value, color = region)) + geom_line() +
  labs(y = x$units[1], title = "Net electricity generation", subtitle = "From all fuels")
```

<img src="eia-geoset2-1.png" title="plot of chunk geoset2" alt="plot of chunk geoset2" width="100%" />

Even more convenient is that these names are available in R. See the `datasets::state.*` functions and the geoset vignette.
---
title: "EIA categories"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{EIA categories}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Top level categories

It is helpful to be aware of various categories of data that are available. Each category has a unique ID number that is required to access associated data. The first call below does not include an ID. The result is a list of two data frames. The first is metadata associated with that position in the category hierarchy. The second is the child category information.


```r
library(eia)
# eia_set_key("yourkey") # set API key if not already set globally
```

Here is the top-level category information.


```r
eia_cats()
#> $category
#> # A tibble: 1 x 3
#>   category_id name          notes
#>   <chr>       <chr>         <chr>
#> 1 371         EIA Data Sets ""   
#> 
#> $childcategories
#> # A tibble: 14 x 2
#>    category_id name                               
#>          <int> <chr>                              
#>  1           0 Electricity                        
#>  2       40203 State Energy Data System (SEDS)    
#>  3      714755 Petroleum                          
#>  4      714804 Natural Gas                        
#>  5      711224 Total Energy                       
#>  6      717234 Coal                               
#>  7      829714 Short-Term Energy Outlook          
#>  8      964164 Annual Energy Outlook              
#>  9     1292190 Crude Oil Imports                  
#> 10     2123635 U.S. Electric System Operating Data
#> 11     2134384 International Energy Data          
#> 12     2251604 CO2 Emissions                      
#> 13     2631064 International Energy Outlook       
#> 14     2889994 U.S. Nuclear Outages
```

## Child categories

The child category IDs can be used to query results from that category.


```r
eia_cats(0)
#> $category
#> # A tibble: 1 x 4
#>   category_id parent_category_id name        notes
#>   <chr>       <chr>              <chr>       <chr>
#> 1 0           371                Electricity ""   
#> 
#> $childcategories
#> # A tibble: 19 x 2
#>    category_id name                                                              
#>          <int> <chr>                                                             
#>  1           1 Net generation                                                    
#>  2          35 Total consumption                                                 
#>  3          32 Total consumption (Btu)                                           
#>  4          36 Consumption for electricity generation                            
#>  5          33 Consumption for electricity generation (Btu)                      
#>  6          37 Consumption for useful thermal output                             
#>  7          34 Consumption for useful thermal output (Btu)                       
#>  8        1017 Plant level data                                                  
#>  9          38 Retail sales of electricity                                       
#> 10          39 Revenue from retail sales of electricity                          
#> 11          40 Average retail price of electricity                               
#> 12     1718389 Number of customer accounts                                       
#> 13       41137 Fossil-fuel stocks for electricity generation                     
#> 14       41138 Receipts of fossil fuels by electricity plants                    
#> 15       41139 Receipts of fossil fuels by electricity plants (Btu)              
#> 16       41140 Average cost of fossil fuels for electricity generation           
#> 17       41141 Average cost of fossil fuels for electricity generation (per Btu) 
#> 18       41142 Quality of fossil fuels in electricity generation : sulfur content
#> 19       41143 Quality of fossil fuels in electricity generation : ash content
```

View the immediate child categories for a given parent category.


```r
eia_child_cats(389)
#> # A tibble: 4 x 2
#>   category_id name             
#>         <int> <chr>            
#> 1         742 Coal             
#> 2         743 Petroleum liquids
#> 3         744 Petroleum coke   
#> 4         745 Natural gas
```

## Parent categories

View all parent categories for a given child category.


```r
eia_parent_cats(742)
#> # A tibble: 6 x 4
#>   category_id name                    notes parent_category_id
#>   <chr>       <chr>                   <chr> <chr>             
#> 1 371         EIA Data Sets           ""    <NA>              
#> 2 0           Electricity             ""    371               
#> 3 32          Total consumption (Btu) ""    0                 
#> 4 372         By sector               ""    32                
#> 5 389         Electric utility        ""    372               
#> 6 742         Coal                    ""    389
```

## Output format

The default is to return tidy data in a tibble data frame. For `eia_cats` you can set `tidy = FALSE` to return the list returned by `jsonlite::fromJSON` without any further processing or `tidy = NA` to return the raw JSON data as a character string.


```r
eia_cats(tidy = FALSE)
#> $request
#> $request$category_id
#> [1] 371
#> 
#> $request$command
#> [1] "category"
#> 
#> 
#> $category
#> $category$category_id
#> [1] "371"
#> 
#> $category$parent_category_id
#> NULL
#> 
#> $category$name
#> [1] "EIA Data Sets"
#> 
#> $category$notes
#> [1] ""
#> 
#> $category$childcategories
#>    category_id                                name
#> 1            0                         Electricity
#> 2        40203     State Energy Data System (SEDS)
#> 3       714755                           Petroleum
#> 4       714804                         Natural Gas
#> 5       711224                        Total Energy
#> 6       717234                                Coal
#> 7       829714           Short-Term Energy Outlook
#> 8       964164               Annual Energy Outlook
#> 9      1292190                   Crude Oil Imports
#> 10     2123635 U.S. Electric System Operating Data
#> 11     2134384           International Energy Data
#> 12     2251604                       CO2 Emissions
#> 13     2631064        International Energy Outlook
#> 14     2889994                U.S. Nuclear Outages
#> 
#> $category$childseries
#> list()
cat(eia_cats(tidy = NA))
#> {"request":{"category_id":371,"command":"category"},"category":{"category_id":"371","parent_category_id":null,"name":"EIA Data Sets","notes":"","childcategories":[{"category_id":0,"name":"Electricity"},{"category_id":40203,"name":"State Energy Data System (SEDS)"},{"category_id":714755,"name":"Petroleum"},{"category_id":714804,"name":"Natural Gas"},{"category_id":711224,"name":"Total Energy"},{"category_id":717234,"name":"Coal"},{"category_id":829714,"name":"Short-Term Energy Outlook"},{"category_id":964164,"name":"Annual Energy Outlook"},{"category_id":1292190,"name":"Crude Oil Imports"},{"category_id":2123635,"name":"U.S. Electric System Operating Data"},{"category_id":2134384,"name":"International Energy Data"},{"category_id":2251604,"name":"CO2 Emissions"},{"category_id":2631064,"name":"International Energy Outlook"},{"category_id":2889994,"name":"U.S. Nuclear Outages"}],"childseries":[]}}
```
---
title: "EIA geosets"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{EIA geosets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Time series by region

Geosets are metadata structures relating series together. Requesting series using the `geoset` API endpoint is much the same as requesting series data using the `series` endpoint. Instead of `eia_series`, use `eia_geoset`.

The main difference is that you must provide a geoset ID to `id` and a `region` argument. Both may be vectors of multiple series and regions. The function returns the combination of time series datasets and regions that exist. The API will not return all geographic entities associated with a geoset. You are required to specify which region(s) you want and they must be associated with the given geoset ID.


```r
library(eia)
library(dplyr)
library(tidyr)
library(ggplot2)

# eia_set_key("yourkey") # set API key if not already set globally
id <- "ELEC.GEN.ALL-99.A"
region <- c("USA-CA", "USA-NY")
(x <- eia_geoset(id, region[1], n = 3))
#> # A tibble: 1 x 11
#>   geoset_id       setname                             f     units            series_id        name                                       region latlon start end   data        
#>   <chr>           <chr>                               <chr> <chr>            <chr>            <chr>                                      <chr>  <chr>  <chr> <chr> <list>      
#> 1 ELEC.GEN.ALL-9~ Net generation : all fuels : all s~ A     thousand megawa~ ELEC.GEN.ALL-CA~ Net generation : all fuels : California :~ USA-CA <NA>   2001  2019  <tibble [3 ~
```

## Groups of regions

If you want data for all fifty states for example, you can set `region = "USA"`; you do not need to make a vector of all fifty state IDs. However, if you want certain region subsets, there are other options besides making a vector of these values. For example, the `eia_geoset` accepts shorthand descriptions of specific, popular subsets of states; so popular in fact that their labels ship with R itself in the `datasts` package.

You can provide simple state abbreviations (without the `USA-` prefix), state names, and more to the point, US Census regions and divisions. These are two hierarchical sets of US states.


```r
tibble(state.abb, state.name, state.region, state.division)
#> # A tibble: 50 x 4
#>    state.abb state.name  state.region state.division    
#>    <chr>     <chr>       <fct>        <fct>             
#>  1 AL        Alabama     South        East South Central
#>  2 AK        Alaska      West         Pacific           
#>  3 AZ        Arizona     West         Mountain          
#>  4 AR        Arkansas    South        West South Central
#>  5 CA        California  West         Pacific           
#>  6 CO        Colorado    West         Mountain          
#>  7 CT        Connecticut Northeast    New England       
#>  8 DE        Delaware    South        South Atlantic    
#>  9 FL        Florida     South        South Atlantic    
#> 10 GA        Georgia     South        South Atlantic    
#> # ... with 40 more rows
```

Provide the associated label and `eia_geoset` recognizes subsets of US states. Even this can be a vector. In the example below, `region` consists of Alaska plus the states belonging to the New England census division.


```r

(x <- eia_geoset(id, c("AK", "New England"), n = 1))
#> # A tibble: 7 x 11
#>   geoset_id       setname                             f     units           series_id        name                                        region latlon start end   data        
#>   <chr>           <chr>                               <chr> <chr>           <chr>            <chr>                                       <chr>  <chr>  <chr> <chr> <list>      
#> 1 ELEC.GEN.ALL-9~ Net generation : all fuels : all s~ A     thousand megaw~ ELEC.GEN.ALL-AK~ Net generation : all fuels : Alaska : all ~ USA-AK <NA>   2001  2019  <tibble [1 ~
#> 2 ELEC.GEN.ALL-9~ Net generation : all fuels : all s~ A     thousand megaw~ ELEC.GEN.ALL-CT~ Net generation : all fuels : Connecticut :~ USA-CT <NA>   2001  2019  <tibble [1 ~
#> 3 ELEC.GEN.ALL-9~ Net generation : all fuels : all s~ A     thousand megaw~ ELEC.GEN.ALL-MA~ Net generation : all fuels : Massachusetts~ USA-MA <NA>   2001  2019  <tibble [1 ~
#> 4 ELEC.GEN.ALL-9~ Net generation : all fuels : all s~ A     thousand megaw~ ELEC.GEN.ALL-ME~ Net generation : all fuels : Maine : all s~ USA-ME <NA>   2001  2019  <tibble [1 ~
#> 5 ELEC.GEN.ALL-9~ Net generation : all fuels : all s~ A     thousand megaw~ ELEC.GEN.ALL-NH~ Net generation : all fuels : New Hampshire~ USA-NH <NA>   2001  2019  <tibble [1 ~
#> 6 ELEC.GEN.ALL-9~ Net generation : all fuels : all s~ A     thousand megaw~ ELEC.GEN.ALL-RI~ Net generation : all fuels : Rhode Island ~ USA-RI <NA>   2001  2019  <tibble [1 ~
#> 7 ELEC.GEN.ALL-9~ Net generation : all fuels : all s~ A     thousand megaw~ ELEC.GEN.ALL-VT~ Net generation : all fuels : Vermont : all~ USA-VT <NA>   2001  2019  <tibble [1 ~
x$data[[1]]
#> # A tibble: 1 x 3
#>   value date        year
#>   <dbl> <date>     <int>
#> 1 6071. 2019-01-01  2019

region <- "Middle Atlantic"
x <- eia_geoset(id, region, n = 12)
select(x, region, data) %>% unnest(cols = data)
#> # A tibble: 36 x 4
#>    region  value date        year
#>    <chr>   <dbl> <date>     <int>
#>  1 USA-NJ 71019. 2019-01-01  2019
#>  2 USA-NJ 75034. 2018-01-01  2018
#>  3 USA-NJ 75645. 2017-01-01  2017
#>  4 USA-NJ 77611. 2016-01-01  2016
#>  5 USA-NJ 74609. 2015-01-01  2015
#>  6 USA-NJ 68051. 2014-01-01  2014
#>  7 USA-NJ 64751. 2013-01-01  2013
#>  8 USA-NJ 65263. 2012-01-01  2012
#>  9 USA-NJ 64694. 2011-01-01  2011
#> 10 USA-NJ 65682. 2010-01-01  2010
#> # ... with 26 more rows

unnest(x, cols = data) %>%
  ggplot(aes(date, value, color = region)) + geom_line() +
  labs(y = x$units[1], title = "Net electricity generation", subtitle = "From all fuels")
```

<img src="geoset-geoset2-1.png" title="plot of chunk geoset2" alt="plot of chunk geoset2" width="100%" />

## Relations

There is also a `relation` argument that accepts an optional relation ID. If one is provided, `eia_geoset` will switch to the API `relation` endpoint. A relation is another metadata structure that applies to geosets and relates summary statistics associated with geoset IDs to their composite statistics. This makes it easier to obtain variables that facet the data, e.g., by sector or fuel type.

The EIA `relation` API endpoint is officially supported according to the online EIA API documentation, but unfortunately that endpoint does not appear to function at the time of current package release.
---
title: "API details"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{API details}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## API key

### Register a key with EIA

Obtaining an API key is easy and free. 

Pulling data from the US Energy Information Administration (EIA) API requires a registered API key. A key can be obtained at no cost [here](https://www.eia.gov/opendata/register.php). A valid email and agreement to the API Terms of Service is required to obtain a key.

It is important to store your API key somewhere secure. Do not commit it to a repository or otherwise share it. For example, you could store it in your `.Renviron` file.

### Key storage and retrieval

You can always provide the `key` argument to every API function call, but you do not have to. There are getter and setter helpers available to make using `eia` functions a more seamless experience.

`eia_set_key` gives you the option of storing your key in any of three places via the `store` argument:

* `store = "env"`: the package environment that is created when the package is loaded (default method)
* `store = "options"`: in the global `options`
* `store = "sysenv"`: as a system environment variable via `Sys.setenv`.

The last two options require the name-value pair to be named `EIA_KEY = "yourkey"`. These three options also are the order of precedence if you do not specify the `store` argument.

This setup also allows you to store a key that will override another key. This is because `eia_get_key` checks these three storage methods in this order and stops as soon as it finds a key. If you need it to check a specific location, you can specify `store`.

As an example, if the key already exists in the system environment and you plan to pass `key` to functions explicitly, you could start as follows:


```r
library(eia)
key <- Sys.getenv("EIA_KEY") 

# or:
key <- eia_get_key()
```

If you need to set it, you can do so as follows.


```r
# eia_set_key("yourkey")
# eia_get_key() # retrieve it
```

API functions in `eia` use `eia_get_key()` with no arguments as the default value of their `key` argument, checking in the order shown above for an existing key. This way you do not need to repeatedly provide it.

Note that despite the name and behavior, storing an environment variable with `Sys.setenv` (and thus `eia_set_key(key, store = "sysenv")`) is not persistent; the key is lost when the R session terminates, just as it is with the other two session-based options. If you want a persistent key, you must manually add your key somewhere like `.Renviron`. In that cases, you never need `eia_set_key` and `eia_get_key` will retrieve the `EIA_KEY` environment variable. See the package documentation for more details on key options.

In this and subsequent vignettes, you will not see a key being set because it is already an environment variable. You will also not see it used explicitly by any functions because the default behavior is to look up the key in the environment.

## API requests

The EIA API can of course impose its own rate-limiting and other limitations on usage by a given API key. If you use the API improperly or otherwise violate any Terms of Service, the EIA may withdraw your API access. However, the `eia` package also helps prevent accidental overuse by having default settings that limit the potential for making unnecessary API calls. It does this in two ways, both of which allow optional configuration:

*   caching API results in memory using session-based memoization
*   minimum delay between API calls

By default the `eia` package prevents you from accidentally making too many requests too quickly, but it also offers sensible flexibility.

### Memoization

All functions in `eia` that make API calls use memoization by default. They will not make the same API call twice in one R session. A call is made once and the result is cached. Calling the same function with the identical arguments again will only returned the cached result. 

This approach limits the potential for accidentally using the EIA API more than necessary. This is fine for most uses cases. However, if you use your API key to access data that is updated very often, or you have a long-running R process such as a Shiny app on a server that may need to periodically update the data associated with a specific API call, you can set `cache = FALSE`.

Run this example of the same request made with and without memoization. You will notice the cached result by the immediate return.


```r
system.time(eia_cats()) # API call; cache result
system.time(eia_cats()) # read from cache
system.time(eia_cats(cache = FALSE)) # API call
```

Results are cached in memory for the duration of the R session, but you can clear the cache at any time.


```r
eia_clear_cache()
system.time(eia_cats())
```

This allows you to update the cached result. You can reset the cache for only specific endpoints using the following functions.

*    `eia_clear_cats`
*    `eia_clear_series`
*    `eia_clear_geoset`

### Anti-DOS measures

Regardless of overall rate limiting imposed by the EIA API, the `eia` package sets a minimum wait time of one second between successive API calls. In most cases this is an irrelevant safeguard. Most `eia` functions make a single API call and requests for data often take a full second anyway once you factor in the subsequent data manipulation in R.

However, there are cases where you might want to make multiple calls back to back programmatically and perhaps you are initially unsure how many requests will be made or how quickly these requests may execute. The default minimum wait between API calls is a precaution that helps you be a good neighbor.

You can turn this off with `options` if not needed; for example, a case where you know that your API calls will be small in number and you have no reason to be concerned about exceeding the request limits associated with your API key. The default requires you to make an active decision about how to use the API with your own key and API limits in mind.

A call to `eia_parent_cats` is a good example. This function is recursive, but say you know the number of calls is going to be small; it is overkill to impose the additional wait. Note that in order to show this example, it is necessary to turn off memoization to avoid returning a cached result.




```r
system.time(x <- eia_parent_cats(742, cache = FALSE))
#>    user  system elapsed 
#>    0.12    0.00    5.89

options(eia_antidos = 0)
system.time(x <- eia_parent_cats(742, cache = FALSE))
#>    user  system elapsed 
#>    0.08    0.00    0.80

x
#> # A tibble: 6 x 4
#>   category_id name                    notes parent_category_id
#>   <chr>       <chr>                   <chr> <chr>             
#> 1 371         EIA Data Sets           ""    <NA>              
#> 2 0           Electricity             ""    371               
#> 3 32          Total consumption (Btu) ""    0                 
#> 4 372         By sector               ""    32                
#> 5 389         Electric utility        ""    372               
#> 6 742         Coal                    ""    389
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geoset.R
\name{eia_geoset}
\alias{eia_geoset}
\title{EIA geoset data}
\usage{
eia_geoset(
  id,
  region,
  relation = NULL,
  start = NULL,
  end = NULL,
  n = NULL,
  tidy = TRUE,
  cache = TRUE,
  key = eia_get_key()
)
}
\arguments{
\item{id}{character, geoset series ID, may be a vector. See details.}

\item{region}{character, region ID, may be a vector. Data available for the intersection of \code{id} and \code{region} is returned.}

\item{relation}{logical, make a geoset relation query instead of a geoset query. The series \code{id} is the same but is queried differently. Currently not supported, see details.}

\item{start}{start date. Providing only a start date will return up to the maximum 100 results if available.}

\item{end}{end date. Providing only an end date will a single result for that date.}

\item{n}{integer, length of series to return ending at most recent value or at \code{end} date if also provided. Ignored if \code{start} is not \code{NULL}.}

\item{tidy}{logical, return a tidier result. See details.}

\item{cache}{logical, cache result for duration of R session using memoization. See details.}

\item{key}{API key: character if set explicitly; not needed if key is set globally. See \code{\link{eia_set_key}}.}
}
\value{
a tibble data frame (or a list, or character, depending on \code{tidy} value)
}
\description{
Obtain EIA geoset data.
}
\details{
\code{id} may be a vector. This should only be done with \code{tidy = TRUE} if the tidied results can be properly row bound.
The geoset API calls allow multiple regions, but the API expects a single series ID.
This function allows multiple series, but must make one API call per series ID.
There is an expectation of similarly formatted series that can be row bound.
If the IDs are for differently structured data that cannot be tidily row bound,
you may as well make separate requests since each requires a unique API call either way.

By default, additional processing is done to return a tibble data frame.
Set \code{tidy = FALSE} to return only the initial list result of \code{jsonlite::fromJSON}.
Set \code{tidy = NA} to return the original JSON as a character string.

Set to \code{cache = FALSE} to force a new API call for updated data.
Using \code{FALSE} always makes a new API call and returns the result from the server.
\code{TRUE} uses memoization on a per R session basis, caching the result of the function call in memory for the duration of the R session.
You can reset the entire cache by calling \code{eia_clear_cache()}.

The EIA \code{relation} API endpoint is officially supported according to the online EIA API documentation, but that endpoint does not appear to function at the time of current package release.
}
\examples{
\dontrun{
# use eia_set_key() to store stored API key
id <- paste0("ELEC.GEN.ALL-99.", c("A", "Q", "M"))
region <- c("USA-CA", "USA-NY")

eia_geoset(id[1], region[1], start = 2016)
eia_geoset(id[2], region, n = 5)
eia_geoset(id[3], region[2], end = 2016, n = 5)

# multiple series counted as a single API call
x <- eia_geoset(id, region[1], end = 2016, n = 2)
x[, c("region", "data")]

# Use direct US state abbreviations or names;
# Use US Census region and division names.
x <- eia_geoset(id[2], c("AK", "New England"), end = 2016, n = 1)
x[, c("region", "data")]
}
}
\seealso{
\code{\link{eia_clear_cache}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/series.R
\name{eia_series_metadata}
\alias{eia_series_metadata}
\alias{eia_series_updates}
\alias{eia_series_dates}
\alias{eia_series_range}
\alias{eia_series_cats}
\title{EIA series metadata}
\usage{
eia_series_metadata(id, cache = TRUE, key = eia_get_key())

eia_series_updates(id, cache = TRUE, key = eia_get_key())

eia_series_dates(id, cache = TRUE, key = eia_get_key())

eia_series_range(id, cache = TRUE, key = eia_get_key())

eia_series_cats(id, tidy = TRUE, cache = TRUE, key = eia_get_key())
}
\arguments{
\item{id}{character, series ID, may be a vector.}

\item{cache}{logical, cache result for duration of R session using memoization.}

\item{key}{API key: character if set explicitly; not needed if key is set globally. See \code{\link{eia_set_key}}.}

\item{tidy}{logical, return a tidier result. See details.}
}
\value{
a tibble data frame
}
\description{
Make a small request to obtain a data frame containing metadata.
}
\details{
Dates are provided in \code{eia_series_dates} for the convenience of working with the EIA date string format;
for example: maintaining order, generating sequences, computing intervals,
and other operations that work well with dates but would be difficult using arbitrary strings.
Keep in mind that of course these are not real dates, in the sense that you cannot map a year to a specific date.

\code{eia_series_updates} returns a data frame of most recent series update times for \code{id}.
Like the other metadata helpers, this does require an API call to the series to obtain the relevant metadata.
This can be useful if you are only interested in these update times for a specific set of series IDs.
If you need to know the most recent update stamps for a large set of series, you should use \code{\link{eia_updates}}
instead, which makes an API call specifically to the EIA \code{updates} endpoint for specific EIA categories by category ID.

\code{eia_series_cats} differs from the other functions in that it makes an API call directly to the \code{series categories} endpoint.
Like other functions that return endpoint-specific output, it accepts the \code{tidy} argument for control over output structure.
By default, additional processing is done to return a list containing tibble data frames.
Set \code{tidy = FALSE} to return only the initial list result of \code{jsonlite::fromJSON}.
Set \code{tidy = NA} to return the original JSON as a character string.
}
\examples{
\dontrun{
# use eia_set_key() to store stored API key
id <- paste0("ELEC.CONS_TOT_BTU.COW-AK-1.", c("A", "Q", "M"))

eia_series_metadata(id)
eia_series_updates(id)
eia_series_dates(id)
eia_series_range(id)
eia_series_cats(id)
}
}
\seealso{
\code{\link{eia_updates}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/categories.R
\name{eia_cats}
\alias{eia_cats}
\alias{eia_child_cats}
\alias{eia_parent_cats}
\title{EIA categories}
\usage{
eia_cats(id = NULL, tidy = TRUE, cache = TRUE, key = eia_get_key())

eia_child_cats(id, cache = TRUE, key = eia_get_key())

eia_parent_cats(id, cache = TRUE, key = eia_get_key())
}
\arguments{
\item{id}{integer, category ID. If \code{NULL}, the API root category.}

\item{tidy}{logical, return a tidier result. See details.}

\item{cache}{logical, cache result for duration of R session using memoization. See details.}

\item{key}{API key: character if set explicitly; not needed if key is set globally. See \code{\link{eia_set_key}}.}
}
\value{
for \code{eia_cats}, a list of tibble data frames (or a less processed list, or character, depending on \code{tidy} value); others functions return a tibble data frame.
}
\description{
Obtain EIA categories.
}
\details{
By default, additional processing is done to return a list containing tibble data frames.
Set \code{tidy = FALSE} to return only the initial list result of \code{jsonlite::fromJSON}.
Set \code{tidy = NA} to return the original JSON as a character string.

Set to \code{cache = FALSE} to force a new API call for updated data.
Using \code{FALSE} always makes a new API call and returns the result from the server.
\code{TRUE} uses memoization on a per R session basis, caching the result of the function call in memory for the duration of the R session.
You can reset the entire cache by calling \code{eia_clear_cache}.

\code{eia_child_cats} returns only the immediate child categories. \code{eia_parent_cats} returns all parents.
These are wrappers around \code{eia_cats} and always return a tibble data frame.
}
\examples{
\dontrun{
# use eia_set_key() to store stored API key
eia_cats()

eia_child_cats(389) # immedate children
eia_parent_cats(742) # all parents
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/categories.R
\name{eia_updates}
\alias{eia_updates}
\title{EIA data updates}
\usage{
eia_updates(
  id = NULL,
  deep = FALSE,
  n = 50,
  start = 1,
  tidy = TRUE,
  key = eia_get_key()
)
}
\arguments{
\item{id}{integer, category ID, may be a vector. If \code{NULL}, the API root category.}

\item{deep}{logical, if \code{TRUE}, return information on all child series. If \code{FALSE} (default), return only for the category \code{id}.}

\item{n}{integer, maximum number of rows of series to return. Defaults to 50; maximum permitted by the API is 10,000.}

\item{start}{integer, row to start from, defaults to 1.}

\item{tidy}{logical, return a tidier result. See details.}

\item{key}{API key: character if set explicitly; not needed if key is set globally. See \code{\link{eia_set_key}}.}
}
\value{
a tibble data frame (or a list, or character, depending on \code{tidy} value)
}
\description{
Obtain information on EIA data series updates for a given category to avoid having to make requests for data that have not been updated since your last request.
}
\details{
This function returns paginated results of the most recent update dates for data series.
\code{n} and \code{start} help with stepping through chunks.

If you need to know the most recent update stamps for a large set of series, you should use this function,
which makes an API call specifically to the EIA \code{updates} endpoint for specific EIA categories by category ID.
If you are only interested in update times for a specific set of series IDs,
you can use \code{\link{eia_series_updates}}.
Note that while this function accepts a vector of IDs for \code{id}, it must make one API call per ID.

By default, additional processing is done to return a tibble data frame.
Set \code{tidy = FALSE} to return only the initial list result of \code{jsonlite::fromJSON}.
Set \code{tidy = NA} to return the original JSON as a character string.
}
\examples{
\dontrun{
# use eia_set_key() to store stored API key
eia_updates(742, n = 5)
}
}
\seealso{
\code{\link{eia_series_updates}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dates.R
\name{eiadate}
\alias{eiadate}
\alias{eiadate_to_date}
\alias{date_to_eiadate}
\alias{eiadate_to_date_seq}
\title{EIA date parsing}
\usage{
eiadate_to_date(x)

date_to_eiadate(x, date_format = c("A", "Q", "M", "W", "D", "H"))

eiadate_to_date_seq(start, end, weekly = FALSE)
}
\arguments{
\item{x}{character, EIA date string; character or date object for regular dates. See details.}

\item{date_format}{EIA date format: "A", "Q", "M", "W", "D", "H".
These stand for annual, quarterly, monthly, weekly, daily, hourly. See details.}

\item{start}{start EIA date or date.}

\item{end}{end EIA date or date.}

\item{weekly}{logical. See details.}
}
\description{
Helper functions for manipulating and converting between regular year-month-day date strings and EIA date string notation.
}
\details{
There is no reason to mix EIA date formats in this context. Functions that take EIA date strings expect a consistent format.
Also, EIA date formats are parsed automatically from the dates themselves.
However, daily and weekly use the same format. Too avoid ambiguity in \code{eia_date_seq}, daily is assumed; set \code{weekly = TRUE} to treat as weekly.

When providing a real date or date string, such as to \code{date_to_eiadate}, dates should be in \code{yyyy-mm-dd} format,
or at least any format that can be parsed by \code{lubridate::ymd} or \code{lubridate::ymd_hms} for dates and hourly date times, respectively.

\code{"HL"} is not a supported date format. Use \code{"H"}. The API does not
translate the date and time when using \code{"HL"} anyhow; it simply appends
the date string with the number of hours time difference.
}
\examples{
eiadate_to_date(c("201803", "201804"))

date_to_eiadate("2018-05-14", "A")
date_to_eiadate("2018-05-14", "Q")
date_to_eiadate("2018-05-14", "M")

(x <- eiadate_to_date_seq("2018Q1", "2018Q4"))
date_to_eiadate(x, "Q")

(x <- eiadate_to_date("20190102T16Z"))
date_to_eiadate(x, "H")
(x <- eiadate_to_date_seq("20190102T16Z", "20190102T19Z"))
date_to_eiadate(x, "H")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{eia_clear_cache}
\alias{eia_clear_cache}
\alias{eia_clear_cats}
\alias{eia_clear_series}
\alias{eia_clear_geoset}
\title{Clear API results cache}
\usage{
eia_clear_cache()

eia_clear_cats()

eia_clear_series()

eia_clear_geoset()
}
\description{
Reset the results of API calls that are currently cached in memory.
}
\details{
\code{eia_clear_cache} clears the entire cache. The other functions clear the cache associated with specific endpoints.
}
\examples{
\dontrun{
key <- Sys.getenv("EIA_KEY") # your stored API key
system.time(eia_cats(key))
system.time(eia_cats(key))
eia_clear_cache()
system.time(eia_cats(key))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key.R
\name{eia_key}
\alias{eia_key}
\alias{eia_set_key}
\alias{eia_get_key}
\title{Set and get API key}
\usage{
eia_set_key(key, store = c("env", "options", "sysenv"))

eia_get_key(store = c("env", "options", "sysenv"))
}
\arguments{
\item{key}{character, API key.}

\item{store}{character, method for storing API key. See details.}
}
\value{
\code{eia_get_key} returns the key string or \code{NULL} with a warning. \code{eia_set_key} returns a success message or an error.
}
\description{
Set and get API key
}
\details{
Setter and getter helpers allow you to store your EIA API key in one of three ways.
Their use is optional. You can always pass the API key string to the \code{key} argument of any package function that requires it,
but you do not have to.

By default the \code{key} argument for these functions is \code{key = eia_get_key()}.
If your key has been stored in a manner that can be retrieved,
then you can call all the package API functions without having to provide the \code{key} argument repeatedly.
}
\section{Key storage methods}{

If you have already set your key globally somewhere using \code{eia_set_key}, \code{eia_get_key} will retrieve it.
You can add the \code{EIA_KEY = "yourkey"} key-value pair to \code{options()} or as a system environment variable yourself and \code{eia_get_key}
will pick it up as long as you use the name \code{EIA_KEY}. For convenience you can do this in your R session with \code{eia_set_key}.
It gives you three options for how to store the key. The default is to use the \code{eia} package environment that is created when the package is loaded.
}

\section{Precedence}{

Choose one method when setting a key. When getting the key, the three locations are checked in the order:
package environment, \code{options()}, then the system environment. To override the order, specify the method explicitly and the check will only occur there.
This also makes it possible to override a system level key by working with one stored in the package environment or \code{options()}.
}

\section{Persistence}{

Note that none of these three storage methods, including \code{"sysenv"} are persistent; the stored key is lost when the R session is terminated.
A key that is stored outside of R as a system environment variable is retrievable with \code{eia_get_key},
just like those set in an R session with \code{eia_set_key} and \code{store = "sysenv"}.
However, if you truly want the key to persist as an environment variable when R terminates, you must manually add it somewhere like \code{.Renviron};
\code{Sys.setenv} in R cannot achieve this.
}

\examples{
eia_set_key("fake")
eia_get_key()
# eia_get_key("options") returns an error if not set
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eia.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr} package for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reports.R
\name{eia_report}
\alias{eia_report}
\alias{report_drilling_productivity}
\title{Download data for various EIA reports}
\usage{
eia_report(id, ...)

report_drilling_productivity()
}
\arguments{
\item{id}{character, the report ID. See details for the list of available reports.}

\item{...}{arguments passed to individual report data functions.}
}
\value{
a list, typically a list of data frames
}
\description{
These functions download data for various EIA reports found on the EIA website but not necessarily available through the EIA API.
}
\details{
The wrapper function and the individual report functions do not make API calls and do not require an API key.
}
\examples{
\dontrun{
x <- eia_report("drilling productivity")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eia.R
\docType{package}
\name{eia}
\alias{eia}
\title{eia: EIA API wrapper}
\description{
This package provides API access to data from the US \href{https://www.eia.gov/}{Energy Information Administration} (EIA).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/series.R
\name{eia_series}
\alias{eia_series}
\title{EIA data series}
\usage{
eia_series(
  id,
  start = NULL,
  end = NULL,
  n = NULL,
  tidy = TRUE,
  cache = TRUE,
  key = eia_get_key()
)
}
\arguments{
\item{id}{character, series ID, may be a vector.}

\item{start}{start date. Providing only a start date will return up to the maximum 100 results if available.}

\item{end}{end date. Providing only an end date will a single result for that date.}

\item{n}{integer, length of series to return ending at most recent value or at \code{end} date if also provided. Ignored if \code{start} is not \code{NULL}.}

\item{tidy}{logical, return a tidier result. See details.}

\item{cache}{logical, cache result for duration of R session using memoization. See details.}

\item{key}{API key: character if set explicitly; not needed if key is set globally. See \code{\link{eia_set_key}}.}
}
\value{
a tibble data frame (or a list, or character, depending on \code{tidy} value)
}
\description{
Obtain EIA data series.
}
\details{
By default, additional processing is done to return a tibble data frame.
Set \code{tidy = FALSE} to return only the initial list result of \code{jsonlite::fromJSON}.
Set \code{tidy = NA} to return the original JSON as a character string.

Set to \code{cache = FALSE} to force a new API call for updated data.
Using \code{FALSE} always makes a new API call and returns the result from the server.
\code{TRUE} uses memoization on a per R session basis, caching the result of the function call in memory for the duration of the R session.
You can reset the entire cache by calling \code{eia_clear_cache()}.
}
\examples{
\dontrun{
# use eia_set_key() to store stored API key
id <- paste0("ELEC.GEN.ALL-AK-99.", c("A", "Q", "M"))

x1 <- eia_series(id[1], start = 2016)
x2 <- eia_series(id[2], n = 5)
x3 <- eia_series(id[3], end = 2016, n = 5)
x1$data[[1]]
x2$data[[1]]
x3$data[[1]]

# multiple series counted as a single API call
x <- eia_series(id, end = 2016, n = 2)
x$data
}
}
\seealso{
\code{\link{eia_clear_cache}}
}
