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

<!-- README.md is generated from README.Rmd. Please edit that file -->
[![](https://badges.ropensci.org/254_status.svg)](https://github.com/ropensci/onboarding/issues/254) [![Build Status](https://travis-ci.org/ropensci/popler.svg?branch=master)](https://travis-ci.org/ropensci/popler) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/popler?branch=master&svg=true)](https://ci.appveyor.com/project/AldoCompagnoni/popler) [![codecov.io](https://codecov.io/github/ropensci/popler/coverage.svg?branch=master)](https://codecov.io/github/ropensci/popler?branch=master)

Popler
======

`popler` is the R package to browse and query the `popler` data base. `popler` is a PostgreSQL data base that contains population-level datasets from the US long term ecological research (LTER) network. This package is currently only available on GitHub, but our ultimate goal is to submit it to CRAN. A detailed explanation on the structure of the `popler` online database is contained in the [dratf of the manuscript presenting popler](https://github.com/texmiller/popler-ms/blob/master/popler_ms.pdf) and in the [dedicated vignette](https://github.com/AldoCompagnoni/popler/blob/master/vignettes/popler-database-structure.Rmd). The `popler` database is organized around four types of tables:

A. The tables containing information on population abundance. Population abundance can be of five types: count, biomass, cover, density, and at the individual-level.

B. A table containg to link each abundance value to the taxonomic units it refers to. In `popler`, "taxonomic unit" generally refers to a species.

C. The location of each "site". With "site" we refers to the largest spatial replicate available for a dataset. Many datasets provide abundance data from more than one site.

D. Metadata information referring to each separate dataset.

<img src="vignettes/img/schema.png" width="1280" />

Popler rationale
================

The package `popler` aims to facilitate finding, retrieving, and using time-series data of population abundance associated with the US LTER network. To *find datasets*, the functions in `popler` aid in understanding and browsing the metadata information referring to each dataset. To *retrieve* data, the function `pplr_get_data()` downloads single datasets or groups of datasets. These downloads share the same data structure. To *use* downloaded data, the package provides ancillary functions to consult and cite the original data sources, examine the temporal replication of each data set, and methods for a couple of `dplyr` verbs to assist with data manipulation.

### Installation

------------------------------------------------------------------------

``` r
# Install stable version once on CRAN (hopefully soon!)
install.packages('popler')


# Install development version now
if(!require(devtools, quietly = TRUE)) {
  install.packages(devtools)
}

devtools::install_github('ropensci/popler')
```

All exported functions in the `popler` R package use the `pplr_` prefix. Moreover, the functions use lazy and/or tidy evaluation, meaning you do not need to manually quote most inputs.

Finding datasets
================

------------------------------------------------------------------------

##### Dictionary of variables

We suggest to start exploring the metadata variables describing each dataset in popler using `pplr_dictionary()`. This function works in two ways:

1.  It provides a general description of each metadata variable. This happens when this function is called without arguments; for example, when calling `pplr_dictionary()`.
2.  It provides the possible values of its unquoted arguments. For example, when calling `pplr_dictionary( proj_metadata_key )`.

The output of `pplr_dictionary()` is a data frame showing a description of each metadata variable.

``` r
library(popler)
pplr_dictionary()
```

    ##             variable
    ## 1              title
    ## 2  proj_metadata_key
    ## 3             lterid
    ## 4           datatype
    ## 5    structured_data
    ## 6          studytype
    ## 7     duration_years
    ## 8          community
    ## 9          structure
    ## 10         treatment
    ## 11          lat_lter
    ## 12          lng_lter
    ## 13           species
    ## 14           kingdom
    ## 15            phylum
    ## 16             class
    ## 17             order
    ## 18            family
    ## 19             genus
    ##                                                description
    ## 1                                         title of project
    ## 2                                        unique project id
    ## 3                                                lter name
    ## 4              type of abundance data (e.g. count,biomass)
    ## 5  are abundance observations grouped (e.g. based on age)?
    ## 6                     experimental or observational study?
    ## 7                             duration of project in years
    ## 8                     does data set contain multiple taxa?
    ## 9                            types of indidivual structure
    ## 10                                      types of treatment
    ## 11                                      lter site latitude
    ## 12                                     lter site longitude
    ## 13                    specific epithet of a taxonomic unit
    ## 14                                                 kingdom
    ## 15                                                  phylum
    ## 16                                                   class
    ## 17                                                   order
    ## 18                                                  family
    ## 19                                                   genus

However, this function is more powerful when used with an argument. When `pplr_dictionary()` is provided with the name of a metadata variable, it returns the possible unique values of the variable. For example, providing `datatype` shows that popler contains five types of abundance data:

``` r
pplr_dictionary( datatype )
```

    ## $`datatype (NA)`
    ## [1] "individual"  "count"       "cover"       "biomass"     "density"    
    ## [6] "basal_cover"

Additionally, the `pplr_report_dictionary()` function generates an `Rmd` file and renders it into html. This html contains both the meaning of variables, and their unique values.

##### Browsing `popler`

Once you are familiar with the meaning and content of `popler`'s metadata variables, `pplr_browse()` provides the metadata of the studies contained in `popler`. `pplr_browse()` also works with and without an input. Without input, the function produces a data frame including the metadata variables describing every study currently contained in the `popler` database. Note that this data frame is a `tbl` that inherits from the `browse` class. Inputs to `pplr_browse()` allow users to subset this data frame (e.g. `duration_years > 5`). When subsetting, the unique values provided by `pplr_dictionary()` are particularly useful. For more nuanced subsetting of available datasets, the `keyword` argument allows to subset variables using partial matching. Note that `keyword` will act primarily on information contained in the title of studies.

``` r
all_studies <- pplr_browse()

# do not quote logical expressions
long_studies <- pplr_browse(duration_years > 20) 

# keyword is quoted
parasite_studies <- pplr_browse(keyword = 'parasite') 
```

The default settings of both `pplr_browse()` and `pplr_dictionary()` report a subset of the metadata variables contained in popler. To report all variables, set `full_tbl = TRUE`.

``` r
#  vars are quoted
interesting_studies <- pplr_browse(vars = c('duration_years', 'lterid')) 

# Use full_tbl = TRUE to get a table with all possible variables
all_studies_and_vars <- pplr_browse(full_tbl = TRUE)
```

##### Reporting metadata

You can generate a human-readable report on metadata variables of the projects you subset using `pplr_browse` by providing the function with the argument `report = TRUE` . This argument uses `rmarkdown` to render the metadata into an html file, and opens it into your default browser. Alternatively, you can perform the same action described above by providing the `browse` object produced calling `pplr_browse` to the function `pplr_report_metadata()`.

``` r
# generate metadata report for all studies
pplr_browse(report=TRUE)
# alternatively
all_studies <- pplr_browse()
pplr_report_metadata(all_studies)

# generate metadata report for parasite studies
pplr_browse(keyword = 'parasite', report = TRUE)
parasite_studies <- pplr_browse(keyword = 'parasite') 
# alternatively
pplr_report_metadata(parasite_studies)
```

Retrieving data
===============

------------------------------------------------------------------------

Once you explored the metadata and decided which projects interest you, it's time to actually download the data! `pplr_get_data()` connects to the data base via an API, and downloads the raw data based on the criteria supplied. Alternatively, if you're happy with the projects represented in the `browse` object you created earlier, you can simply pass that object to `pplr_get_data()`. Note that if your `browse` object contains 5 rows, `pplr_get_data()` will download 5 separate datasets. All objects created with `pplr_get_data()` inherit from `get_data` and `data.frame` classes.

``` r
# create a browse object and use it to get data

penguins <- pplr_browse(lterid == 'PAL')

# unpack covariates as well

penguin_raw_data <- pplr_get_data(penguins, cov_unpack = TRUE)

# A very specific query

more_raw_data <- pplr_get_data((proj_metadata_key == 43 | 
                                proj_metadata_key == 25) & 
                                year < 1995 )
```

Using data
==========

------------------------------------------------------------------------

We provide three ancillary functions to facilitate the use of the objects downloaded through `pplr_get_data()`.

First, `pplr_metadata_url()` opens up a webpage containing study details. Before doing scientific analyses, we urge the users to review the peculiarities of each dataset by vetting their online documentation. Importantly, `pplr_metadata_url()` also works on objects produced by `pplr_browse`.

Second, `pplr_cov_unpack()` transforms the data contained in the `covariates` column of each downloaded dataset into separate columns. This can or cannot be useful depending on the objectives of the user. Note: you can also transform covariates into a data frame directly through `pplr_get_data()` by providing the function with argument `cov_unpack = TRUE`.

Third, `pplr_citation()` produces a citation for each downloaded dataset.

##### Spatio-temporal replication

The datasets contained in `popler` present many heterogeneities, especially in terms of their spatio-temporal replication. Most studies present at least a few spatial replicates which were not sampled every year. Note that most datasets in `popler` present at least one additional replication level. These spatial replicates are denoted with numbered variables of the form `spatial_replication_level_X`, where `X` refers to the replication level which can go from 1 to 5. The names of these replication levels (e.g. plot, subplot, etc.) are contained in variable `spatial_replication_level_x_label`.

Once you download a dataset, you can examine the temporal replication of the largest spatial replicate (the site, or `spatial_replication_level_1`) using function `pplr_site_rep_plot()`. This function produces a plot showing whether or not a given site was sampled in a year.

``` r
# download and plot yearly spatial replication for dataset 1
kelp_df      <- pplr_get_data( proj_metadata_key == 1)
pplr_site_rep_plot( kelp_df )
```

![](README_files/figure-markdown_github/spatial_rep_plot-1.png)

Additionally, `pplr_site_rep()` produces either a logical vector for subsetting an existing `get_data` object or a summary table of temporal replication for a given spatial resolution. You can control the minimum frequency of sampling and the minimum duration of sampling using the `freq` and `duration` arguments, respectively. Additionally, you can choose the level of spatial replication to filter providing an integer between 1 and 5 to the `rep_level` argument. `return_logical` allows you to control what is returned by the function. `TRUE` returns a logical vector corresponding to rows of the `get_data` that correspond to spatial replicates that meet the criteria of replication specified in the function. `FALSE` returns a summary table describing the number of samples per year at the selected spatial resolution.

``` r
# Example with piping and subsetting w/ the logical vector output

library(dplyr)

SEV_studies <- pplr_get_data( lterid == 'SEV' datatype == 'invidual' )

long_SEV_studies <- SEV_studies %>%
  .[pplr_site_rep(input = .,
                  duration = 12,
                  rep_level = 3), ] %>%
  pplr_site_rep_plot()

# Or, create the summary table

SEV_summary <- SEV_studies %>% 
  pplr_site_rep(duration = 13,
                rep_level = 1,
                return_logical = FALSE)


# Modify the site_rep_plot() by hand using ggplot2 syntax
library(ggplot2)

pplr_site_rep_plot(long_SEV_studies, return_plot = TRUE) +
  ggtitle('Sevilleta LTER Temporal Replication')
```

##### Data manipulation

`popler` supplies methods for a couple of `dplyr` verbs to assist with data manipulation. `filter` and `mutate` methods are available for objects of `browse` and `get_data` classes. Other `dplyr` verbs change the structure of the object too much for those classes to retain their meaning so they are not included in the package, but one can still use them for their own purposes.

``` r
penguins_98 <- filter(penguin_raw_data, year == 1998)

class(penguins_98) # classes are not stripped from objects

penguins_98_true <- mutate(penguins_98, penguins_are = 'Awesome')

class(penguins_98_true)
```

Further information
===================

------------------------------------------------------------------------

Additional information on `popler` is contained in a manuscript, and in the vignettes associated with the R package.

The [manuscript, currently in draft form](https://github.com/texmiller/popler-ms/blob/master/popler_ms.pdf), presents the `popler` database, the R package, and our recommendations on how to use them.

The R package contains three vignettes: one vignette [illustrates the structure of the popler database](https://github.com/AldoCompagnoni/popler/blob/master/vignettes/popler-database-structure.Rmd), and two vignettes provide [an introduction](https://github.com/AldoCompagnoni/popler/blob/master/vignettes/introduction-to-popler.Rmd) and [a more detailed look](https://github.com/AldoCompagnoni/popler/blob/master/vignettes/vetting-popler.Rmd) at the intended workflow of the popler package.

In case these vignettes do not cover your particular use case, you still have questions, or you discover a bug, please don't hesitate to create an [issue](https://github.com/AldoCompagnoni/popler/issues).

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
Popler 0.2.0 (2019-07-23)
=========================

### NEW FEATURES

  * This new version of `Popler` queries the online database via an API

### MINOR IMPROVEMENTS

  * Faster tests through the use of stubs

### BUG FIXES

  * Fix bugs contained in `report_metadata`

### DEPRECATED AND DEFUNCT

  * `hello_world()` now deprecated and will be removed in a
     future version, use `hello_mars()`

### DOCUMENTATION FIXES

  * Created new vignette, `popler-database-structure`, to provide a full documentation the database structure
  * Updated README file to include information on database structure
  * Removed typos and clarified functions documentation


Popler 0.1.0 
=========================

  * Introduces the first stable release of the `popler` _R_ package for querying the PostgreSQL data base of the same name.
## Test environments
* local Windows 10 install, R 3.5.1
* ubuntu 12.04 (on travis-ci), R-oldrelease, R-release, and R-devel
* MacOS Sierra 10.12.6 (on Travis-CI), R-oldrelease, R-release
* Windows Server 2012 (on Appveyor-CI) -- R-oldrelease, R-release
* win-builder (devel and release)

All of the above passed with the following results: 

0 errors | 0 warnings | 0 notes


* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
# Contributing to popler

This outlines how to propose a change to popler. 
### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it's a problem. If you've found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  New code should follow the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). It is highly
unlikely that we will accept changes without associated tests, but there are 
instances where an exception can be made (some code can be very hard to test).
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the current
development version header describing the changes made followed by your GitHub
username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to
abide by its terms.

---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![](https://badges.ropensci.org/254_status.svg)](https://github.com/ropensci/onboarding/issues/254)
[![Build Status](https://travis-ci.org/ropensci/popler.svg?branch=master)](https://travis-ci.org/ropensci/popler)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/popler?branch=master&svg=true)](https://ci.appveyor.com/project/AldoCompagnoni/popler)
[![codecov.io](https://codecov.io/github/ropensci/popler/coverage.svg?branch=master)](https://codecov.io/github/ropensci/popler?branch=master)


# Popler

`popler` is the R package to browse and query the `popler` data base. `popler` is a PostgreSQL data base that contains population-level datasets from the US long term ecological research (LTER) network. This package is currently only available on GitHub, but our ultimate goal is to submit it to CRAN.
A detailed explanation on the structure of the `popler` online database is contained in the [dratf of the manuscript presenting popler ](https://github.com/texmiller/popler-ms/blob/master/popler_ms.pdf) and in the [dedicated vignette](https://github.com/AldoCompagnoni/popler/blob/master/vignettes/popler-database-structure.Rmd). 
The `popler` database is organized around four types of tables:

A. The tables containing information on population abundance. Population abundance can be of five types: count, biomass, cover, density, and at the individual-level.

B. A table containg to link each abundance value to the taxonomic units it refers to. In `popler`, "taxonomic unit" generally refers to a species.

C. The location of each "site". With "site" we refers to the largest spatial replicate available for a dataset. Many datasets provide abundance data from more than one site.

D. Metadata information referring to each separate dataset.

```{r, echo=FALSE}
knitr::include_graphics("vignettes/img/schema.png")
```

# Popler rationale

The package `popler` aims to facilitate finding, retrieving, and using time-series data of population abundance associated with the US LTER network.
To _find datasets_, the functions in `popler` aid in understanding and browsing the metadata information referring to each dataset.
To _retrieve_ data, the function `pplr_get_data()` downloads single datasets or groups of datasets. These downloads share the same data structure.
To _use_ downloaded data, the package provides ancillary functions to consult and cite the original data sources, examine the temporal replication of each data set, and methods for a couple of `dplyr` verbs to assist with data manipulation.


### Installation

------
```{r install popler, eval = FALSE}

# Install stable version once on CRAN (hopefully soon!)
install.packages('popler')


# Install development version now
if(!require(devtools, quietly = TRUE)) {
  install.packages(devtools)
}

devtools::install_github('ropensci/popler')


```

All exported functions in the `popler` R package use the `pplr_` prefix. Moreover, the functions use lazy and/or tidy evaluation, meaning you do not need to manually quote most inputs.


# Finding datasets

------

##### Dictionary of variables

We suggest to start exploring the metadata variables describing each dataset in popler using `pplr_dictionary()`. This function works in two ways:

1. It provides a general description of each metadata variable. This happens when this function is called without arguments; for example, when calling `pplr_dictionary()`.
1. It provides the possible values of its unquoted arguments. For example, when calling `pplr_dictionary( proj_metadata_key )`.

The output of `pplr_dictionary()` is a data frame showing a description of each metadata variable. 


```{r load popler, eval = TRUE, include = FALSE}
library(popler)
```

```{r dictionary funs, eval = TRUE}

library(popler)
pplr_dictionary()

```

However, this function is more powerful when used with an argument. When `pplr_dictionary()` is provided with the name of a metadata variable, it returns the possible unique values of the variable. For example, providing `datatype` shows that popler contains five types of abundance data:

```{r dictionary argument, eval = TRUE}
pplr_dictionary( datatype )

```


Additionally, the `pplr_report_dictionary()` function generates an `Rmd` file and renders it into html. This html contains both the meaning of variables, and their unique values.



##### Browsing `popler`

Once you are familiar with the meaning and content of `popler`'s metadata variables, `pplr_browse()` provides the metadata of the studies contained in `popler`. `pplr_browse()` also works with and without an input. Without input, the function produces a data frame including the metadata variables describing every study currently contained in the `popler` database. Note that this data frame is a `tbl` that inherits from the `browse` class. Inputs to `pplr_browse()` allow users to subset this data frame (e.g. `duration_years > 5`). When subsetting, the unique values provided by `pplr_dictionary()` are particularly useful. For more nuanced subsetting of available datasets, the `keyword` argument allows to subset variables using partial matching. Note that `keyword` will act primarily on information contained in the title of studies.

```{r browse base, eval = FALSE}

all_studies <- pplr_browse()

# do not quote logical expressions
long_studies <- pplr_browse(duration_years > 20) 

# keyword is quoted
parasite_studies <- pplr_browse(keyword = 'parasite') 

```

The default settings of both `pplr_browse()` and `pplr_dictionary()` report a subset of the metadata variables contained in popler. To report all variables, set `full_tbl = TRUE`. 

```{r browse advanced, eval = FALSE}

#  vars are quoted
interesting_studies <- pplr_browse(vars = c('duration_years', 'lterid')) 

# Use full_tbl = TRUE to get a table with all possible variables
all_studies_and_vars <- pplr_browse(full_tbl = TRUE)

```


##### Reporting metadata

You can generate a human-readable report on metadata variables of the projects you subset using `pplr_browse` by providing the function with the argument `report = TRUE` . This argument uses `rmarkdown` to render the metadata into an html file, and opens it into your default browser.
Alternatively, you can perform the same action described above by providing the `browse` object produced calling `pplr_browse` to the function `pplr_report_metadata()`.

``` {r report_metadata, eval = FALSE}

# generate metadata report for all studies
pplr_browse(report=TRUE)
# alternatively
all_studies <- pplr_browse()
pplr_report_metadata(all_studies)

# generate metadata report for parasite studies
pplr_browse(keyword = 'parasite', report = TRUE)
parasite_studies <- pplr_browse(keyword = 'parasite') 
# alternatively
pplr_report_metadata(parasite_studies)

```


# Retrieving data

-----

Once you explored the metadata and decided which projects interest you, it's time to actually download the data! `pplr_get_data()` connects to the data base via an API, and downloads the raw data based on the criteria supplied. Alternatively, if you're happy with the projects represented in the `browse` object you created earlier, you can simply pass that object to `pplr_get_data()`. Note that if your `browse` object contains 5 rows, `pplr_get_data()` will download 5 separate datasets. All objects created with `pplr_get_data()` inherit from `get_data` and `data.frame` classes.

```{r get_data, eval = FALSE}

# create a browse object and use it to get data

penguins <- pplr_browse(lterid == 'PAL')

# unpack covariates as well

penguin_raw_data <- pplr_get_data(penguins, cov_unpack = TRUE)

# A very specific query

more_raw_data <- pplr_get_data((proj_metadata_key == 43 | 
                                proj_metadata_key == 25) & 
                                year < 1995 )

```

# Using data

-----

We provide three ancillary functions to facilitate the use of the objects downloaded through `pplr_get_data()`.

First, `pplr_metadata_url()` opens up a webpage containing study details. Before doing scientific analyses, we urge the users to review the peculiarities of each dataset by vetting their online documentation. Importantly, `pplr_metadata_url()` also works on objects produced by `pplr_browse`.

Second, `pplr_cov_unpack()` transforms the data contained in the `covariates` column of each downloaded dataset into separate columns. This can or cannot be useful depending on the objectives of the user. Note: you can also transform covariates into a data frame directly through `pplr_get_data()` by providing the function with argument `cov_unpack = TRUE`.

Third, `pplr_citation()` produces a citation for each downloaded dataset.


##### Spatio-temporal replication

The datasets contained in `popler` present many heterogeneities, especially in terms of their spatio-temporal replication. Most studies present at least a few spatial replicates which were not sampled every year. Note that most datasets in `popler` present at least one additional  replication level. These spatial replicates are denoted with numbered variables of the form `spatial_replication_level_X`, where `X` refers to the replication level which can go from 1 to 5. The names of these replication levels (e.g. plot, subplot, etc.) are contained in variable `spatial_replication_level_x_label`.

Once you download a dataset, you can examine the temporal replication of the largest spatial replicate (the site, or `spatial_replication_level_1`) using function `pplr_site_rep_plot()`. This function produces a plot showing whether or not a given site was sampled in a year.

```{r spatial_rep_plot, echo=TRUE, results='hide', fig.keep='all', message = FALSE}

# download and plot yearly spatial replication for dataset 1
kelp_df 	 <- pplr_get_data( proj_metadata_key == 1)
pplr_site_rep_plot( kelp_df )

```

Additionally, `pplr_site_rep()` produces either a logical vector for subsetting an existing `get_data` object or a summary table of temporal replication for a given spatial resolution. You can control the minimum frequency of sampling and the minimum duration of sampling using the `freq` and `duration` arguments, respectively. 
Additionally, you can choose the level of spatial replication to filter providing an integer between 1 and 5 to the `rep_level` argument.
`return_logical` allows you to control what is returned by the function. `TRUE` returns a logical vector corresponding to rows of the `get_data` that correspond to spatial replicates that meet the criteria of replication specified in the function. `FALSE` returns a summary table describing the number of samples per year at the selected spatial resolution.


```{r spatial_rep, eval = FALSE}

# Example with piping and subsetting w/ the logical vector output

library(dplyr)

SEV_studies <- pplr_get_data( lterid == 'SEV' datatype == 'invidual' )

long_SEV_studies <- SEV_studies %>%
  .[pplr_site_rep(input = .,
                  duration = 12,
                  rep_level = 3), ] %>%
  pplr_site_rep_plot()

# Or, create the summary table

SEV_summary <- SEV_studies %>% 
  pplr_site_rep(duration = 13,
                rep_level = 1,
                return_logical = FALSE)


# Modify the site_rep_plot() by hand using ggplot2 syntax
library(ggplot2)

pplr_site_rep_plot(long_SEV_studies, return_plot = TRUE) +
  ggtitle('Sevilleta LTER Temporal Replication')


```


##### Data manipulation

`popler` supplies methods for a couple of `dplyr` verbs to assist with data manipulation. `filter` and `mutate` methods are available for objects of `browse` and `get_data` classes. Other `dplyr` verbs change the structure
of the object too much for those classes to retain their meaning so they are not included in the package, but one can still use them for their own purposes.

```{r dplyr-verbs, eval = FALSE}

penguins_98 <- filter(penguin_raw_data, year == 1998)

class(penguins_98) # classes are not stripped from objects

penguins_98_true <- mutate(penguins_98, penguins_are = 'Awesome')

class(penguins_98_true)

```


# Further information

------

Additional information on `popler` is contained in a manuscript, and in the  vignettes associated with the R package.

The [manuscript, currently in  draft form](https://github.com/texmiller/popler-ms/blob/master/popler_ms.pdf), presents the `popler` database, the R package, and our recommendations on how to use them.

The R package contains three vignettes: one vignette [illustrates the structure of the popler database](https://github.com/AldoCompagnoni/popler/blob/master/vignettes/popler-database-structure.Rmd), and two vignettes provide [an introduction](https://github.com/AldoCompagnoni/popler/blob/master/vignettes/introduction-to-popler.Rmd) and [a more detailed look](https://github.com/AldoCompagnoni/popler/blob/master/vignettes/vetting-popler.Rmd) at the intended workflow of the popler package.

In case these vignettes do not cover your particular use case, you still have questions, or you discover a bug, please don't hesitate to create an [issue](https://github.com/AldoCompagnoni/popler/issues).

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Vetting popler"
author: "Aldo Compagnoni, Sam Levin"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vetting popler}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction: identifying groups of data sets

The `popler` R package was built to foster scientific synthesis using LTER long-term population data. The premise of such synthesis is using data from many research projects that share characteristics of scientific interest. To identify projects sharing salient attributes, `popler` uses the metadata information associated with each LTER project. In particular, it is fairly easy to select projects based on one or more of the following features:

 1. Replication, temporal or spatial.
 2. Taxonomic group(s).
 3. Study characteristics. 
 4. Geographic location.

Vetting the database based on these criteria is intuitive. However, `popler` also facilitates identifying data sets in other ways. Below we provide several examples on how to select LTER data based on the four types of features described above. Moreover, in the final section we also show how to carry out more complicated types of searches.

### 1. Replication
#### Temporal replication

If you are interested in long-term data, you will likely want to select projects based on how many years the data was collected for. This is straightforward: 

```{r, warning = FALSE, message = FALSE}

library(popler)
pplr_browse(duration_years > 10)

```

Note that most LTER projects contemplate sampling at a yearly or sub-yearly frequency. Thus, studies longer than 10 years often guarantee a longitudinal series of 10 or more observations. Note that the `duration_years` variable is calculated as `studyendyr - studystartyr`. Thus, an additional variable named `samplefreq` characterizes the approximate sample frequency of each study. 

```{r, warning = FALSE, message = FALSE}

pplr_dictionary(samplefreq)
pplr_browse(samplefreq == "monthly")

```

Note that `samplefreq` is **not** a default variable included in the `pplr_dictionary` or `pplr_browse()` functions. This can be viewed by specifying the `full_tbl = TRUE` argument in either of the above functions.

###1. Spatial replication

#### Before downloading data

If you wish to select data sets based on their spatial replication, you need to consider that `popler` organizes data in nested spatial levels. For example, in many plant studies data is collected at the *plot* level, which can be nested within *block*, which in turn can be nested within *site*. `popler` labels spatial levels using numbers. Spatial level 1 is the coarsest level of replication which contains all other spatial replicates. In the example above, spatial level 1 is *site*, spatial level 2 is *block*, and spatial level 3 is *plot*. `popler` allows for a total of 5 spatial levels. Given the above, you can select studies based on three criteria:

 1. The total number of spatial replicates.
 
 2. The number of replicates within a specific spatial level.
 
 3. The number of nested spatial replicates.

Below we provide three examples for each one of these respective cases.

```{r}
pplr_browse(tot_spat_rep > 100)
pplr_browse(spatial_replication_level_5_number_of_unique_reps > 1)
pplr_browse(n_spat_levs == 3)
```

#### After downloading data

Users can also explore the spatial and temporal replication of the data more explicitly after downloading it with `pplr_get_data()` through two function: `pplr_site_rep()` and `pplr_site_rep_plot()`.

`pplr_site_rep()` provides two options for exploring data that meet temporal replication requirements at a given spatial resolution. The user can choose to filter data by specifying a minimum sampling frequency per year and a minimum number of years that sample with that frequency. Because this function uses the sampling dates to calculate the frequency, it provides additional information that is not contained in the `samplefreq` column of the main metadata table.


```{r, eval = FALSE}
# download some data (note: this download is >100MB)
SEV <- pplr_get_data(proj_metadata_key == 21)

# Create a summary table containing names of replication levels that contain 2 samples per year for 10 years. 
SEV_long_studies <- pplr_site_rep(SEV, 
                                  freq = 2, 
                                  duration = 10, 
                                  return_logical = FALSE)

# you can also subset it directly using the function and specifying it to return a logical vector
subset_vec <- pplr_site_rep(SEV,
                            freq = 2,
                            duration = 10,
                            return_logical = TRUE)
# store subset of data
SEV_long_data <- SEV[subset_vec, ]



```

Users can also visualize the frequency of sampling at the coarsest level of spatial replication using the `pplr_site_rep_plot()` function. This generates a `ggplot` that denotes whether or not a particular site was sampled in a particular year. Note that the coarsest level of spatial replication is called _site_ and it is contained in the variable `spatial_replication_level_1`.

```{r, eval = FALSE}

library(ggplot2)

# return the plot object w/ return_plot = TRUE
pplr_site_rep_plot(SEV_long_data, return_plot = TRUE) +
  ggtitle("Long Term Data from Sevilleta LTER")
  
# or return an invisible copy of the input data and keep piping
library(dplyr)
SEV_long_data %>%
  pplr_site_rep_plot(return_plot = FALSE) %>%
  pplr_report_metadata()

```


###2. Taxonomic group

`popler` is not limited to specific taxonomic groups, and it currently contains mostly data on animals and plants. To select information based on taxonomic groups, simply specify which group and which category you wish to select. The default settings of popler provide seven taxonomic groups: kingdom, phylum, class, order, family, genus, and species in each request. Column `sppcode` provides the identifier, usually an alphanumeric code, associated with each taxonomic entity in the original dataset.
Note that not all LTER studies provide full taxonomic information; hence, browsing studies by taxonomic information will provide partial results (in the example below, not all insects studies will be identified).

```{r, warning = FALSE, message = FALSE}
pplr_dictionary(class)
pplr_browse(class == "Insecta")
```

Note that the taxonomic information returned in `pplr_browse()` is housed in a data structure called _list column_. Each entry of this list column is itself a list that contains a `data.frame` with eight columns. Users can access this information using the following syntax.

```{r, warning = FALSE, message = FALSE}

insects <- pplr_browse(class == 'Insecta')

# access the taxonomic table from the first project in the insects object
insects$taxas[[1]]

# second table (etc.)
insects$taxas[[2]]

```


###3. Study characteristics

Metadata information provides a few variables to select studies based on their design. In particular:

 - `studytype`: indicates whether the study is observational or experimental. Options are `obs` or `exp` for observational and experimental studies, respectively.
 - `treatment_type`: type of treatments, if study is experimental.
 - `community`: indicates whether the project provides data on multiple species. Options are `yes` or `no`.
 - `structured_data`: indicates whether the project provides information on population structure. For example, a population can be sub-divided in age, size, or developmental classes. Options are `yes` or `no`.

Below we show how to use these three fields. 

```{r, warning = FALSE, message = FALSE}
pplr_dictionary(community)
pplr_browse(community == "no") # 20 single-species studies

pplr_dictionary(treatment)
nrow( pplr_browse(treatment == "fire") ) # 21 fire studies

pplr_dictionary(studytype)
nrow( pplr_browse(studytype == "obs") ) # 78 observational studies

```

### 4. Geographic location.

To select studies based on the latitude and longitude of LTER headquarters around which datasets were, or are being collected, simply use the `lat_lter` and `lng_lter` numeric variables: 

```{r, warning = FALSE, message = FALSE}
pplr_dictionary( lat_lter, lng_lter )
pplr_browse( lat_lter > 40 & lng_lter < -100 ) # single-species studies
```

### 5. More complicated searches

Popler allows carrying out more complicated searches by allowing to i)  simultaneously search several types of metadata variables, and ii) search studies matching a string pattern. In the first case, simply provide the function `pplr_browse()` with a logical statement regarding more than one metadata variable. For example, if you want studies on plants with at least 4 nested spatial levels, and 10 years of data:

```{r, eval = FALSE}

pplr_browse(kingdom == "Plantae" & n_spat_levs == 4 & duration_years > 10)

```

In the second case, the keyword argument in function `pplr_browse()` will search for string patterns within the metadata of each study. For example, in case we were interested in studies using traps:

```{r, eval = FALSE}

pplr_browse(keyword = 'trap')

```

Note that the keyword argument works with regular expressions as well:

```{r, eval = FALSE}

# look for studies that include the words "trap" or "spatial"
pplr_browse(keyword = 'trap|spatial')

```
---
title: "Introduction to popler"
author: "Aldo Compagnoni, Sam Levin"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to popler}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The popler R package is an interface that allows browsing and querying population data collected at Long Term Ecological Research (LTER) network sites located across the United States of America. A subset of the population data from the LTER network is contained in an online database called `popler`. The popler R package is an interface to this online database, allowing users to:

- explore *what* type of population data is contained in the `popler` database
- download data contained in the `popler` database
- filter and validate the data once it is downloaded


## Installation

The popler R package is currently in the development phase, and it should be downloaded directly from its [GitHub page](https://github.com/AldoCompagnoni/popler). Before doing this, make sure to install the [devtools R package](https://cran.r-project.org/web/packages/devtools/README.html). Once devtools is installed on your machine, install and load popler:

```{r,warnings = FALSE, message = FALSE}
# devtools::install_github("AldoCompagnoni/popler", build_vignettes = TRUE)
library(popler)
```

## Metadata: *what* type of data is contained in the popler database?

`popler` provides data from hundreds of research projects. The metadata of these projects allow understanding *what* population data are  provided by each project. The `popler` R package provides three functions to explore these metadata.

### pplr_dictionary()

`pplr_dictionary()` shows:

- what the variables contained in the metadata of each project and their meaning are. 
- what data these variables contain.

To see metadata variables and their meaning:  

```{r}
pplr_dictionary()
```

To show what data each variable actually contains, specify one or more variable:

```{r}
pplr_dictionary(lterid, duration_years)
```

Last, but not least, the same information provided by `pplr_dictionary` can be visualized in an html page containing hyperlinks. To open such html page, execute the `pplr_report_dictionary` function. 

```{r eval = FALSE}
pplr_report_dictionary()
```

### pplr_browse()

`pplr_browse()` accesses and subsets the popler metadata table directly. Calling the function returns a table that contains the metadata of all the projects in `popler`:

```{r, eval=FALSE}
all_studies <- pplr_browse()
```

This metadata table can be subset by specifying a logical expression. This is useful to focus on datasets of interest.

```{r}
poa_metadata  <- pplr_browse(genus == "Poa" & species == "fendleriana")
poa_metadata
```

Moreover, akin to `pplr_report_dictionary()`, browse can generate a report and open it as an html page. To do so, set the `report` variable to `TRUE`. Alternatively, you can pass an object created by `pplr_browse()` to `pplr_report_metadata()` to create the same report.  

```{r, eval = FALSE}

pplr_browse(lterid == "SEV", report = TRUE)

SEV <- pplr_browse(lterid == "SEV")

pplr_report_metadata(SEV)

```

#### The keyword argument

`pplr_browse()` can also single out projects based on partial matching across the metadata variables that contain characters. Specify the character string you want to search using the `keyword` argument (note that this function ignores variables that contain numeric values):

```{r eval = FALSE}
pplr_browse(keyword = "parasite", report = TRUE)
```

## Download data

Once you identified one or more datasets of interest, download their raw data using `pplr_get_data()`. You can use this function to download data in three ways:

1. Providing `pplr_get_data()` with an object created through `pplr_browse()`.

2. Providing `pplr_get_data()` with an object created by `pplr_browse()`, and with an additional logical expression to further subset this object of class `browse`.

3. Providing `pplr_get_data()` with a logical expression. This logical expression will typically indicate the specific project(s) the user is interested in downloading.

Below are examples on the three ways to use `pplr_get_data()`:

``` {r eval = FALSE}
# option 1
poa_metadata    <- pplr_browse(genus == "Poa" & species == "fendleriana") 
poa_data        <- pplr_get_data(poa_metadata) 
# option 2
poa_data_11     <- pplr_get_data(poa_metadata, duration_years > 10) 
# option 3
parasite_data   <- pplr_get_data(proj_metadata_key == 25) 

```

Here, we emphasize two important characteristics of `pplr_get_data()`. First, similarly to `pplr_browse()`, the function selects datasets based on the variables described in `pplr_dictionary()`. Second, `pplr_get_data()` will download entire datasets that satisfy user-defined conditions. Hence, for example, in the example above where `genus == "Poa" & species == "fendleriana"`, the function will download three datasets which will include data on _Poa fendleriana_, along with the many other taxa that happen to co-occur with _Poa fendleriana_ in those datasets.


In case you are using a slow internet connection, datasets may take some time to download. Therefore, `popler` provides two utility functions for saving downloaded data locally and efficiently. They are thin wrappers around `saveRDS` and `readRDS` that allow you to store large data sets in highly compressed formats. `.rds` files also have the advantage of rapid read and write times from R, making them optimal for saving data sets for later usage. Note from the examples below: you should *not* specify the file type when specifying the path.

```{r eval = FALSE}

# save the large data set for later usage
pplr_save(poa_data, file = "some/file/path/Poa_Data")

# when you're ready to use it again, pick up where you left off.

poa_data_reloaded <- pplr_load(file = "some/file/path/Poa_Data")

# These will be identical
stopifnot(identical(poa_data, poa_data_reloaded))

```

### Carefully vet the methods of downloaded data sets.

We urge the user to carefully read the documentation of each project before using it for research purposes. Data sets downloaded with `popler` share the same data structure, but each project has its peculiarities. To show the *metadata* of the downloaded data sets, use `pplr_report_metadata` on the data object produced by `pplr_get_data()`. To read the *methods* of each project, click on  the 'metadata link' hyperlink provided in the html page.

```{r eval = FALSE}
pplr_report_metadata(poa_data)

```

### Data structure

In `popler`, datasets are objects produced by `pplr_get_data()` which have the same structure. This structure is documented formally in `vignette('popler-database-structure', package = 'popler')`. Here, we provide a brief description on how spatial replicates and taxonomic information are stored in the database.

Spatial replicates are identified using variables that match the patterns  `spatial_replication_level_X` and `spatial_replication_level_X_label`. Here `X` is a number referring to one of maximum 5 *nested* levels of spatial replication. `X` can vary from 1 to 5, with 1 referring to the largest spatial replication level - the one within which are nested all smaller spatial replicates. So for example, `spatial_replication_level_1` can represent a site, and `spatial_replication_level_2` represents a plot. In this specific case, `spatial_replication_level_1_label` will contain the string 'site', and `spatial_replication_level_2_label` will contain the string 'plot'.

Taxonomic units are identified through species codes in the `sppcode` variable, or through the `genus` and `species` variables. The `sppcode` variable usually contains alphanumeric codes. The `genus` and `species` variables are Latin binomial name. Occasionally, some datasets will contain higher taxonomic classifications (such as `family`, `class`, etc.).

### Spatio-temporal replication

Users can explore the level of temporal replication at each nested level of spatial replication using the `pplr_site_rep_plot()` and `pplr_site_rep()`  functions. 

`pplr_site_rep_plot()` produces a scatterplot that shows which sites (`spatial_replication_level_1`) were sampled in a given year.

`pplr_site_rep()` allows the user to subset datasets downloaded by `pplr_get_data()` based on the _frequency_ and _number of yearly replicates_ contained at a specific level of spatial replication. For example, this function allows to identify the replicates of the second level of spatial replication (e.g. plots within sites) which contain two samples per years (their frequency), for 10 years (the number of yearly replicates).
`pplr_site_rep()` returns a logical vector to subset the `pplr_get_data()` object. For additional examples on how to explore and vet `popler` data, see `vignette('vetting-popler', package = 'popler')`.

## Extra covariates

Most data sets provided by the USA LTER network contain more variables than those accommodated by the schema of `popler`. In order not to loose the original data, `popler` stores all extra information in a character variable named `covariates`. The `popler` package provides two ways to format these covariates into a data frame: the `cov_unpack` argument in `pplr_get data()`, and the `pplr_cov_unpack()` function in `popler`.

Setting the `cov_unpack` argument to `TRUE` returns a data frame that combines the variables of a default query to popler, and the covariates contained in each particular study downloaded through popler:

``` {r}
d_47_cov <- pplr_get_data(proj_metadata_key == 47, cov_unpack = TRUE)
head(d_47_cov)
```

Using the `pplr_cov_unpack()` function on a data frame downloaded using `pplr_get_data()` returns a _separate_ data frame of the covariates contained in the downloaded object.

```{r}
d_47 <- pplr_get_data(proj_metadata_key == 47)
head(pplr_cov_unpack(d_47))
```
---
title: "Popler Database"
author: "Andrew Bibian"
date: "3/16/2019"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Put the title of your vignette here}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

## Popler Overview
Popler is a database of population and individual-level data gathered throughout the Long Term Ecological Research (LTER) stations funded by the USA National Science Foundation.

We define _population data_ datasets as time series on the size or density of a population of a taxonomic unit. The size of a population can be quantified as a count, biomass, or cover class. These measures are always repeated temporally, and they are often repeated spatially as well.

We define _individual-level_ data as information on the attributes of the individuals, or a subset thereof, that make up a population. For example, common attributes of individuals are size, age, and sex.


### Temporal replication
All temporal information within a dataset is stored using 
sampling dates up to, when available, the daily resolution. Since not all datasets have the same temporal resolution, popler stores date information in three separate columns; **day, month, and year**. In any of the date time columns, `NULL` or `-99999` values indicate that the information was not available from the raw data. 

Note that **never do we perform any temporal aggregations of data prior to storing them in Popler.**

### Spatial replication
`popler` subdivides spatial replicates based on their spatial nestedness, because most sampling designs are spatially nested. For example, a study performed at 5 sites, with 3 transects at each site, and 4 quadrats within each transect, contains 3 levels of nested spatial replication:

- Level 1: the site level, 
- Level 2: the transect level
- Level 3: the quadrat level

This referencing scheme allows us to standardize and align datasets collected from a variety of different sampling designs and across different data types.

Metadata on the extent of each spatial sampling unit is recorded when 
available (i.e. km2, m, cm3, etc.). Note, higher levels of spatial replication 
indicate smaller areas of sampling extent. 

Similarly to temporal replication, **never do we aggregate data by levels of spatial replication prior to storing them in Popler.**

#### Spatial Replication Level 1
The first level of spatial replication is also the highest, or coarsest, level of nested spatial replication. We call this level 1 'site'. LTER stations contain permanent 'sites' which are reused across studies. Therefore, "site" labels allow querying data collected at a particular site regardless of the original study a dataset was derived from.

Whenever available in the source metadata, these 'sites' are associated with latitude and longitude coordinates are also recorded if available; `-99999` values exist within these fields for studies that do not record geographic location for their 'site' label. 


## The Popler Database Schema

The following image depicts the relationship structure among the 14 tables of popler. Below we will discuss this schema, table contents, and give definitions. Three tables will not be discussed, because two are related to climate and one is a table for database migration information.

```{r, out.width = "675px", echo=FALSE}
knitr::include_graphics("img/full-schema.png")
```

## 1. LTER Table
Data on the LTER research station.

```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/lter-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='lter')
lter <- read.csv("data/lter.csv")
knitr::kable(lter, booktabs = TRUE)
```

## 2. Study Site Table

This table contains the labels that identify "sites" (spatial replication level 1) used across research projects. Different datasets can use the same "site" label if data were collected at the same sampling location.

```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/study-site-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='study_site')
lter <- read.csv("data/study_site.csv")
knitr::kable(lter, booktabs = TRUE)
```


## 3. Project Table
Metadata related to each separate dataset.

```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/project-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='project')
lter <- read.csv("data/project.csv")
knitr::kable(lter, booktabs = TRUE)
```


## 4. Sites Within Project Table
Site level information regarding starting and ending year of sampling, number of observations, and number of taxonomic units. Note that here, `uniquetaxaunits` refers to the taxonomic units observed within each site, not within each dataset.

```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/site-in-project-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='site_in_project')
lter <- read.csv("data/site_in_project.csv")
knitr::kable(lter, booktabs = TRUE)
```

## 5. Taxonomic Table
Taxonomic information recorded within a project. In this table, taxonomic information refers to each individual site (spatial replication level 1). 

```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/taxa-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='taxa')
lter <- read.csv("data/taxa.csv")
knitr::kable(lter, booktabs = TRUE)
```

## 6. Accepted Taxonomic Table
Table containing "accepted" taxonomic information: this is an attempt to associate taxonomic units in the raw data of popler to taxonomic units accepted in the literature. This taxonomic information also refers to individual sites (spatial replication level 1).

```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/taxa-accepted-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='accepted_taxa')
lter <- read.csv("data/accepted_taxa.csv")
knitr::kable(lter, booktabs = TRUE)
```

## 7. Count Table
Population data where abundance is quantified as number of individuals. Null values filled in.

```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/count-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='count')
lter <- read.csv("data/count.csv")
knitr::kable(lter, booktabs = TRUE)
```

## 8. Biomass Table
Population data where abundance is quantified in terms of biomass.


```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/biomass-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='biomass')
lter <- read.csv("data/biomass.csv")
knitr::kable(lter, booktabs = TRUE)
```

## 9. Density Table
Population data where abundance is quantified in terms of density.


```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/density-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='density')
lter <- read.csv("data/density.csv")
knitr::kable(lter, booktabs = TRUE)
```

## 10. Percent Cover Table
Population data where abundance is quantified in terms of cover.


```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/percent-cover-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='percent_cover')
lter <- read.csv("data/percent_cover.csv")
knitr::kable(lter, booktabs = TRUE)
```

## 11. Individual Table
Individual-level data. The `structure_type` columns refer to the attributes of individuals (e.g. size, age, sex, etc.)


```{r, out.width = "600px", echo=FALSE}
knitr::include_graphics("img/individual-table.png")
```

```{r, echo=FALSE}
# lter <- read_excel("data/metadata.xlsx", sheet='individual')
lter <- read.csv("data/individual.csv")
knitr::kable(lter, booktabs = TRUE)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popler_citation.R
\name{pplr_citation}
\alias{pplr_citation}
\title{Provide citations for a popler object returned by \code{pplr_browse} or
\code{pplr_get_data} object.}
\usage{
pplr_citation(input, bibtex_path = NULL)
}
\arguments{
\item{input}{An object of class \code{browse} or \code{get_data}.}

\item{bibtex_path}{Specify the filename and location for 
the generated markdown file (optional).}
}
\value{
A list of references from \code{input}.
}
\description{
Returns a bibliography, Bibtex citations, and acknowledgment template.
}
\examples{
\dontrun{
# make a browse object
metadata <- pplr_browse(proj_metadata_key \%in\% c(17, 317, 494))

# cite the projects
cite <- pplr_citation(metadata)

# cite$bibliography          # the bibliography
# cite$Bibtex                # Bibtex entries for each dataset
# cite$acknowledgement       # acknowledgement template
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popler.R
\docType{package}
\name{popler}
\alias{popler}
\alias{popler-package}
\title{popler package}
\description{
\code{popler} is a package for interacting with the PostgreSQL data base
with the same name. \code{popler} contains data on long-term population 
dynamics from the LTER network. Every exported function is prefixed with 
\code{pplr_} and then a verb (e.g. \code{pplr_get_data()} or noun (e.g. 
\code{pplr_citation}. Accessing \code{popler} does not require an API key, 
all you need is an internet connection and you are ready to go!
}
\section{Authors}{


Aldo Compagnoni \email{aldo.compagnoni@aggiemail.usu.edu}

Andrew Bibian \email{ajbibian@rice.edu}

Brad Ochocki \email{brad.ochocki@rice.edu}

Sam Levin \email{sam.levin@idiv.de}

Tom Miller \email{tom.miller@rice.edu}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata_url.R
\name{pplr_metadata_url}
\alias{pplr_metadata_url}
\title{Get metadata information from a data object}
\usage{
pplr_metadata_url(input)
}
\arguments{
\item{input}{An object produced by the function \code{pplr_get_data()}.}
}
\description{
Load the webpage containing the metadata of the data sets 
contained in objects produced by \code{pplr_browse} or downloaded through 
\code{pplr_get_data()}. If you downloaded data from multiple projects, 
this function will open multiple webpages. This is a wrapper of function 
\code{browseURL} in \code{base}.
}
\examples{

\dontrun{
# Load the metadata webpages of the projects that contain data from the Poa genus.
fes_d <- pplr_browse(genus == "Festuca")
pplr_metadata_url( fes_d )
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov_unpack.R
\name{pplr_cov_unpack}
\alias{pplr_cov_unpack}
\title{Unpack the covariates contained in the dataset downloaded via 
\code{pplr_get_data()}}
\usage{
pplr_cov_unpack(input)
}
\arguments{
\item{input}{An object of class \code{get_data}.}
}
\value{
A data frame whose columns represent the covariates of the dataset
downloaded via \code{pplr_get_data()}. Note that these covariates are 
contained in the \code{covariates} column datasets downloaded using
\code{pplr_get_data()}.
}
\description{
Create a data frame by "extracting" the \code{covariates} column
contained in an dataset downloaded with \code{pplr_get_data()}.
}
\examples{
\dontrun{
library(dplyr)
demo_d <- pplr_get_data(proj_metadata_key == 8)
as.tbl( pplr_cov_unpack( demo_d ) )
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report_metadata.R
\name{pplr_report_metadata}
\alias{pplr_report_metadata}
\title{Open a report of the metadata of project(s) as an html page}
\usage{
pplr_report_metadata(input, md_file = "./browse.Rmd",
  html_file = "./browse.html")
}
\arguments{
\item{input}{A popler object returned by \code{pplr_browse()} or 
\code{pplr_get_data()}}

\item{md_file}{Specify the filename and location for the generated markdown
file (optional)}

\item{html_file}{Specify the filename and location for
the generated html file (optional)}
}
\value{
An invisible copy of \code{input}.
}
\description{
Generates a readable report of the metadata describing data sets 
contained in popler. The report contains citations, the links to the original 
URL of each data set, and example code to obtain the metadata and data of the
projects represented in the html page.
}
\examples{
\dontrun{
# Full dictionary
one_spp <- pplr_browse(community == "no" & duration_years > 15)
pplr_report_metadata(one_spp)

data <- pplr_get_data(one_spp)
pplr_report_metadata(data) # same as above
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_funs.R
\name{pplr_maps}
\alias{pplr_maps}
\title{Generate maps of LTER sites}
\usage{
pplr_maps(input, return_plot = FALSE)
}
\arguments{
\item{input}{An object created by either \code{pplr_browse()} or
\code{pplr_get_data()}}

\item{return_plot}{logical; if \code{TRUE} function returns the \code{ggplot} 
object for subsequent modification.
If \code{FALSE}, function returns an invisible copy of the \code{input} 
object (useful for piping). Default is \code{FALSE}.}
}
\value{
The \code{input} object (invisibly) or a \code{ggplot2} object.
}
\description{
Generates maps of LTER sites in a given \code{input} object.
Sizes of site markers correspond to the number of studies at a given site.
}
\examples{
\dontrun{

library(dplyr) # make \%>\% available

browse_object <- pplr_browse(proj_metadata_key == 11)

browse_object \%>\%
  pplr_maps() 
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{pplr_summary_table_update}
\alias{pplr_summary_table_update}
\title{Update \code{popler}'s summary table}
\usage{
pplr_summary_table_update()
}
\value{
This function is called for its side effect and does not return 
anything
}
\description{
Automatically retrieve most up to date version of \code{popler}
summary table
}
\note{
The \code{summary_table} is often called internally by popler functions,
 but can also be accessed directly by calling \code{pplr_summary_table_import()}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr_methods.R
\name{mutate.browse}
\alias{mutate.browse}
\alias{mutate.get_data}
\title{Methods for dplyr verbs}
\usage{
\method{mutate}{browse}(.data, ...)

\method{mutate}{get_data}(.data, ...)
}
\arguments{
\item{.data}{A \code{browse} or \code{get_data} object}

\item{...}{Name-value pairs of expressions. Use \code{NULL} to drop a 
variable.}
}
\description{
Add new columns to a \code{browse} or \code{get_data} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/site_replication.R
\name{pplr_site_rep}
\alias{pplr_site_rep}
\alias{pplr_site_rep_plot}
\title{Spatial-temporal replication of data sets}
\usage{
pplr_site_rep(input, freq = 1, duration = 10, rep_level = 1,
  return_logical = TRUE)

pplr_site_rep_plot(input, return_plot = FALSE)
}
\arguments{
\item{input}{An object of produced by \code{pplr_get_data}. Note that this
is not an output from \code{pplr_browse}, as the raw data is required to 
calculate the amount of replication.}

\item{freq}{A number corresponding to the desired annual frequency of 
replicates. Studies that are replicated more frequently will be 
included in the counts and those that replicated less frequently will be 
excluded. 
If \code{return_logical = TRUE}, rows that contain information from sites 
that are replicated at the desired frequency will have a \code{TRUE} value, 
and rows that are not will have a \code{FALSE} value. 
Values greater than 1 will select sampling done multiple times per year. 
For example, \code{freq = 2} indicates a desired sampling frequency of 6
months. Conversely, \code{freq = 0.5} indicates a desired sampling done 
once every 2 years.}

\item{duration}{An integer corresponding to the desired number of yearly
replicates. Rows containing site information from sites with more 
replication will be included, while those with less will be excluded.}

\item{rep_level}{An integer corresponding to the level of spatial 
replication over which verify yearly temporal replication. Values between 1 and 5 
are possible (though higher levels may not be present for some datasets). 
Higher values correspond to higher levels of spatial nestedness. 
The default value of \code{rep_level = 1} corresponds to sites.}

\item{return_logical}{logical; if \code{TRUE}, the function returns a logical
vector. This vector can be used to subset the dataset. If \code{FALSE}, the 
function returns a summary table of class \code{tbl}. This table shows, in 
variable \code{number_of_samples}, how many temporal replicates per year 
are contained by each spatial replicate. Default is \code{TRUE}.}

\item{return_plot}{A logical indicating whether to return a copy of the 
\code{input} data or the \code{ggplot} object created by the function. Use
\code{TRUE} to return the \code{ggplot} object for subsequent modification.
Use \code{FALSE} to return an invisible copy of the \code{input} object 
(useful for piping). Default is \code{FALSE}.}
}
\value{
\code{pplr_site_rep_plot}: \code{input} object (invisibly) or a
\code{ggplot2} object. Use \code{return_plot} to control.

\code{pplr_site_rep}: A \code{tbl} or a logical vector of length 
\code{dim(input)[1]}. Use \code{return_logical} to control.
}
\description{
Functions to examine the number of temporal replicates 
contained within each spatial replication level of a dataset. 
\code{pplr_site_rep_plot} plots the temporal replicates available for 
each site.
\code{pplr_site_rep} produces logical vectors that identify the spatial 
replicates with enough temporal replication, or summary tables.
}
\details{
\code{pplr_site_rep_plot} produces a scatterplot showing the sites 
(\code{spatial_replication_level_1}) and years for which data is available.

\code{pplr_site_rep} works with any level of spatial replication and produces
either a summary table of temporal replication or a logical vector that can be used 
to subset a data set based on the desired frequency and length of time.
}
\examples{
\dontrun{

library(ggplot2)
library(dplyr)

# produce logical vector and subset using it. This can also be piped into a 
# the plotting function for visiualization

good_studies <- pplr_get_data(lterid == 'SEV') \%>\%
                   .[pplr_site_rep(input = .,
                                   duration = 12, 
                                   rep_level = 3), ] \%>\%
                   pplr_site_rep_plot()
                                       

# Or, make a neat summary table and decide where to go from there
SEV <- pplr_get_data(lterid == 'SEV')

rep_table <- pplr_site_rep(input = SEV,
                           freq = 0.5,
                           duration = 12,
                           return_logical = FALSE)
 
# pplr_site_rep_plot ---------------
                                                     
# create an unmodified figure
BNZ <- pplr_get_data(lterid == 'BNZ')

pplr_site_rep_plot(BNZ)

# Return the figure instead of the data for subsequent modification
Antarctica <- pplr_get_data(lterid == 'PAL')

pplr_site_rep_plot(Antarctica,
              return_plot = TRUE) + 
   ggtitle("Penguins Rock!")
   
# Use within pipes. Cannot return and modify the figure this way.
pplr_get_data(lterid == 'SEV') \%>\% 
  pplr_site_rep_plot(return_plot = FALSE) \%>\%
  pplr_report_metadata()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dictionary.R
\name{pplr_dictionary}
\alias{pplr_dictionary}
\title{Dictionary of the popler metadata variables}
\usage{
pplr_dictionary(..., full_tbl = FALSE)
}
\arguments{
\item{...}{A sequence of (unquoted) variables specifying one
or more variables of popler's main table for which dictionary 
information is needed}

\item{full_tbl}{logical; If \code{TRUE}, the function
returns a table describing the variables of the full main table.
If \code{FALSE}, the function returns a table describing the standard 
variables. Default is \code{FALSE}.}
}
\description{
Describes the metadata variables contained
in the popler database, and shows their content.
}
\examples{
\dontrun{
# Column names
column_names <- pplr_dictionary(full_tbl = FALSE)

# Dictionary information
dictionary_lter <- pplr_dictionary(lterid, full_tbl = FALSE)

# multiple columns
dictionary_lter_lat <- pplr_dictionary(lterid,lat_lter, full_tbl = FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr_methods.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{filter}
\alias{mutate}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{dplyr}{\code{\link[dplyr]{filter}}, \code{\link[dplyr]{mutate}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{pplr_summary}
\alias{pplr_summary}
\title{search}
\usage{
pplr_summary(limit = 10, offset = 0, ...)
}
\arguments{
\item{limit}{number of records to return, default: 10}

\item{offset}{record number to start at, default: first record}

\item{...}{curl options passed on to [crul::HttpClient]}
}
\description{
search
}
\examples{
# basic example
pplr_summary()
# pass in curl options for debugging, seeing http request details
pplr_summary(verbose = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{pplr_browse}
\alias{pplr_browse}
\title{Browse the metadata of projects contained in the popler database}
\usage{
pplr_browse(..., full_tbl = FALSE, vars = NULL, view = FALSE,
  keyword = NULL, report = FALSE)
}
\arguments{
\item{...}{A logical expression to subset the table containing the metadata of 
datasets contained in popler}

\item{full_tbl}{logical; Should the function returns the standard 
columns, or the full main table? Default is \code{FALSE}.}

\item{vars}{A vector of characters in case the user want to select 
which variables of popler's main table should be selected?}

\item{view}{If TRUE, opens up a spreadsheet-style data viewer.}

\item{keyword}{A string used to select individual datasets based on
pattern matching. The string is matched to every string element in the 
variables of the metadata table in popler.}

\item{report}{logical; If \code{TRUE}, function produces a markdown 
report about each study's metadata, and opens it as a html page. 
Default is \code{FALSE}.}
}
\value{
A data frame combining the metadata of each project 
and the taxonomic units associated with each project.

This data frame is of class \code{popler}, \code{data.frame},
\code{tbl_df}, and \code{tbl}.
}
\description{
pplr_browse() reports the metadata of LTER studies contained in the popler database. 
The user can subset which datasets, and which metadata variables to visualize.
}
\examples{

\dontrun{
# No arguments return the standard 16 columns of popler's main table
default_vars = pplr_browse()

# full_tbl = TRUE returns the full table
all_vars = pplr_browse(full_tbl = TRUE)

# subset only data from the sevilleta LTER, and open the relative report in a html page
sev_data = pplr_browse(lterid == "SEV", report = TRUE)

# consider only plant data sets 
plant_data = pplr_browse(kingdom == "Plantae")

# Select only the data you need
three_columns = pplr_browse(vars = c("title","proj_metadata_key","genus","species"))

# Select only the data you need
study_21 = pplr_browse( proj_metadata_key == 25)

# Select studies that contain word "parasite"
parasite_studies = pplr_browse( keyword = "parasite")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report_dictionary.R
\name{pplr_report_dictionary}
\alias{pplr_report_dictionary}
\title{A user-friendly dictionary of the popler metadata}
\usage{
pplr_report_dictionary(full_tbl = FALSE, md_file = NULL,
  html_file = NULL)
}
\arguments{
\item{full_tbl}{logical; if \code{TRUE} function returns the variables 
contained in the full main table. If \code{FALSE}, functions returns only the
standard variables. Default is \code{FALSE}.}

\item{md_file}{Specify the filename and location for 
the generated markdown file (optional)}

\item{html_file}{Specify the filename and location for the 
generated html file (optional)}
}
\value{
This function is called for its side effects and does not 
return an R object.
}
\description{
Provides information on the variables of metadata contained in the popler 
database, and the kind of data contained in those variables.
}
\examples{
\dontrun{
# Full dictionary
pplr_report_dictionary(full_tbl = TRUE)

# "Abridged" version
pplr_report_dictionary()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{pplr_search}
\alias{pplr_search}
\title{search}
\usage{
pplr_search(proj_metadata_key, limit = 10, offset = 0, ...)
}
\arguments{
\item{proj_metadata_key}{project metadata key}

\item{limit}{number of records to return, default: 10}

\item{offset}{record number to start at, default: first record}

\item{...}{curl options passed on to [crul::HttpClient]}
}
\description{
search
}
\examples{
# basic example
pplr_search(proj_metadata_key = 13)
# pass in curl options for debugging, seeing http request details
pplr_search(proj_metadata_key = 13, verbose = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr_methods.R
\name{filter.browse}
\alias{filter.browse}
\alias{filter.get_data}
\title{Methods for dplyr verbs}
\usage{
\method{filter}{browse}(.data, ...)

\method{filter}{get_data}(.data, ...)
}
\arguments{
\item{.data}{A \code{browse} or \code{get_data} object}

\item{...}{logical conditions}
}
\description{
Subsets a \code{browse} or \code{get_data} object based on
logical statements.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{pplr_get_data}
\alias{pplr_get_data}
\title{Download data from the popler database}
\usage{
pplr_get_data(..., cov_unpack = FALSE)
}
\arguments{
\item{...}{An object produced by \code{pplr_browse} or a logical expression.}

\item{cov_unpack}{logical; if \code{TRUE}, function \code{pplr_cov_unpack} 
is applied to the variable \code{covariates} of the downloaded dataset in 
order to extract the variables contained in therein and combine the new
columns with the default output. Default is \code{FALSE}.}
}
\value{
This data fame is of class \code{get_data}, and \code{data.frame}.
}
\description{
This function downloads datasets contained in the popler database. 
The user can download data directly, using a logical expression, or indirectly, 
using objects created by \code{pplr_browse}.
}
\details{
. By default, the following variables are included when a user calls
\code{pplr_get_data()}.

\itemize{
  \item{\code{authors}}
  \item{\code{authors_contact}} 
  \item{\code{year}} 
  \item{\code{day}} 
  \item{\code{month}}
  \item{\code{sppcode}} 
  \item{\code{genus}}
  \item{\code{species}}
  \item{\code{datatype}}
  \item{\code{spatial_replication_level_1_label}}
  \item{\code{spatial_replication_level_1}}
  \item{\code{spatial_replication_level_2_label}}
  \item{\code{spatial_replication_level_2}}
  \item{\code{spatial_replication_level_3_label}}
  \item{\code{spatial_replication_level_3}}
  \item{\code{spatial_replication_level_4_label}}
  \item{\code{spatial_replication_level_4}}
  \item{\code{spatial_replication_level_5_label}}
  \item{\code{spatial_replication_level_5}}
  \item{\code{proj_metadata_key}}
  \item{\code{structure_type_1}}
  \item{\code{structure_type_2}}
  \item{\code{structure_type_3}}
  \item{\code{structure_type_4}}
  \item{\code{treatment_type_1}}
  \item{\code{treatment_type_2}}
  \item{\code{treatment_type_3}}
  \item{\code{covariates}}
}
}
\examples{
\dontrun{
# browse a study, then get the data associated with it
parasite = pplr_browse(proj_metadata_key == 25)
gh_data = pplr_get_data(parasite)

# insect data sets from the SEV lter site
insect_sev = pplr_browse(class == "Insecta" & lterid == "SEV")
insect_25_yrs96_99 = pplr_get_data(insect_sev)

insect_21_25 = pplr_get_data( (proj_metadata_key == 43 | 
                               proj_metadata_key == 25) )
}

}
