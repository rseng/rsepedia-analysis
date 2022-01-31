ots
===



[![R-check](https://github.com/ropensci/ots/workflows/R-check/badge.svg)](https://github.com/ropensci/ots/actions?query=workflow%3AR-check)

`ots` is an R client to retrieve data from various ocean time series datasets, including:

* BATS
* HOT
* CALCOFI
* LTER Kelp
* UOPG
* more to come...

Jump over to the issues page to suggest data sets to include or comment on ongoing data source integration progress.

What's the point of getting data from the web in R? This way we only have to solve the problem of how to efficiently get a dataset once, then you can benefit from that. In addition, this should allow you to get any changes to the dataset that appear, or corrections. Last, getting data programatically in R should get you one step closer to a reproducible workflow, one that makes science easier primarily for yourself, and for others using your work.

## Install


```r
install.packages("devtools")
devtools::install_github("ropensci/ots")
```


```r
library("ots")
```

## Easy integration with dplyr


```r
library('dplyr')
tbl_df(bats("zooplankton")$data) %>% 
  filter(sieve_size > 1000) %>% 
  group_by(cruise) %>% 
  summarise(mean_water_vol = mean(water_vol))
```

## BATS - Zooplankton dataset


```r
bats("zooplankton")
```

## BATS - Production dataset


```r
bats("production")
```

## HOT dataset


```r
hot()
```

## Channels Islands National Park kelp data


```r
kelp("benthic_cover")
```

## CALCOFI data


```r
calcofi('hydro_cast')
```

## UOPG data

Various datasets available through this source - in this example getting data from Biowatt, and getting the meteorology data. Note that we still need to fix the column names...


```r
(biowatt_met <- uopg(dataset = 'biowatt', type = "meteorology"))
```

## More coming...

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/ots/issues).
* License: MIT
* Get citation information for `ots` in R doing `citation(package = 'ots')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
ots 0.1.0
===============

### NEW FEATURES

* released to CRAN
I have read and agree to the the CRAN policies at 
http://cran.r-project.org/web/packages/policies.html

R CMD CHECK passed on my local OS X install with R 3.2.2 and
R development version, Ubuntu running on Travis-CI, and Windows
R 3.2.2 and devel on Win-Builder.

This is a new submission.

Thanks! 
Scott Chamberlain
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

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

Please note that the `ots` project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
ots
===

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE, 
  message = FALSE,
  eval = FALSE
)
```

[![R-check](https://github.com/ropensci/ots/workflows/R-check/badge.svg)](https://github.com/ropensci/ots/actions?query=workflow%3AR-check)

`ots` is an R client to retrieve data from various ocean time series datasets, including:

* BATS
* HOT
* CALCOFI
* LTER Kelp
* UOPG
* more to come...

Jump over to the issues page to suggest data sets to include or comment on ongoing data source integration progress.

What's the point of getting data from the web in R? This way we only have to solve the problem of how to efficiently get a dataset once, then you can benefit from that. In addition, this should allow you to get any changes to the dataset that appear, or corrections. Last, getting data programatically in R should get you one step closer to a reproducible workflow, one that makes science easier primarily for yourself, and for others using your work.

## Install

```{r eval=FALSE}
install.packages("devtools")
devtools::install_github("ropensci/ots")
```

```{r}
library("ots")
```

## Easy integration with dplyr

```{r}
library('dplyr')
tbl_df(bats("zooplankton")$data) %>% 
  filter(sieve_size > 1000) %>% 
  group_by(cruise) %>% 
  summarise(mean_water_vol = mean(water_vol))
```

## BATS - Zooplankton dataset

```{r}
bats("zooplankton")
```

## BATS - Production dataset

```{r}
bats("production")
```

## HOT dataset

```{r}
hot()
```

## Channels Islands National Park kelp data

```{r}
kelp("benthic_cover")
```

## CALCOFI data

```{r}
calcofi('hydro_cast')
```

## UOPG data

Various datasets available through this source - in this example getting data from Biowatt, and getting the meteorology data. Note that we still need to fix the column names...

```{r}
(biowatt_met <- uopg(dataset = 'biowatt', type = "meteorology"))
```

## More coming...

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/ots/issues).
* License: MIT
* Get citation information for `ots` in R doing `citation(package = 'ots')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kelp.R
\name{kelp}
\alias{kelp}
\alias{kelp_datasets}
\alias{kelp_metadata}
\title{Get Kelp forest community data.}
\usage{
kelp(dataset = "benthic_cover", path = "~/.ots/kelp", overwrite = TRUE)

kelp_datasets()

kelp_metadata(name = NULL, path = "~/.ots/kelp")
}
\arguments{
\item{dataset}{A dataset code name, see \code{\link{kelp_datasets}}}

\item{path}{A path to store the files, Default: \code{~/.ots/kelp}}

\item{overwrite}{(logical) To overwrite the path to store files in or not, Default: TRUE.}

\item{name}{Metadata table name}
}
\description{
Get Kelp forest community data.
}
\details{
After one download o the dataset, you won't have to download the data again.

The \code{kelp} function is to get datasets, e.g., benthic cover data. The output of each call
to \code{kelp} includes the data, and both the headers and variables metadata tables for that
dataset. In addition, the citation to the data is included. See examples below for how to index
to each of those.

The \code{kelp_datasets} function simply lists the datasets available. You can pass the code
to \code{kelp}.

The \code{kelp_metadata} function is to both list the metadata tables available, and to retrieve
those metadata tables, including: sites, data_updates, metadata_updates, history, and species.
}
\examples{
\dontrun{
# list of datasets
kelp_datasets()

# read in various metadata files
## list metadata tables
kelp_metadata()
## get a table
head( kelp_metadata("sites") )
head( kelp_metadata("data_updates") )
head( kelp_metadata("metadata_updates") )
head( kelp_metadata("history") )
head( kelp_metadata("species") )

# get data
(res <- kelp("benthic_cover"))
head(res$headers)
head(res$vars)
res$citation

(res <- kelp("benthic_density"))
(res <- kelp("fish_density"))
(res <- kelp("fish_size"))
(res <- kelp("invert_size"))
head(res$headers)
(res <- kelp("subtidal"))
(res <- kelp("rdfc"))
(res <- kelp("kelp_size"))
(res <- kelp("kelp_supp_dens"))
(res <- kelp("art_recruit"))
}
}
\references{
\url{http://esapubs.org/archive/ecol/E094/245/metadata.php}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uopg.R
\name{uopg}
\alias{uopg}
\title{Get data from the Upper Ocean Proccesses Group (UOPG).}
\usage{
uopg(
  dataset = "smile",
  type = "meteorology",
  path = "~/.ots/uopg",
  overwrite = TRUE
)
}
\arguments{
\item{dataset}{A dataset code name, one of arabian_sea, asrex_91, asrex_93, biowatt, cmo,
coare, coop, fasinex, lotus, mlml89, mlml91, sesmoor, smile, or subduction. Or their
unique abbreviations.}

\item{type}{A data type, one of meteorology, water_velocity, temperature, or salinity. Or their
unique abbreviations.}

\item{path}{A path to store the files, Default: \code{~/.ots/uopg}}

\item{overwrite}{(logical) To overwrite the path to store files in or not, Default: TRUE.}
}
\description{
Get data from the Upper Ocean Proccesses Group (UOPG).
}
\details{
We download NetCDF files from UOPG, using \pkg{ncdf} to parse the data
to a data.frame.
}
\examples{
\dontrun{
# Smile dataaets
(met <- uopg(dataset = 'smile', type = "meteorology"))
(sal <- uopg(dataset = 'smile', type = "salinity"))
(water <- uopg(dataset = 'smile', type = "water"))
(temp <- uopg(dataset = 'smile', type = "temp"))

# biowatt dataaets
(biowatt_met <- uopg(dataset = 'biowatt', type = "meteorology"))
biowatt_met$data$`1`

# lotus datasets
(lotus_met <- uopg(dataset = 'lotus', type = "meteorology"))
lotus_met$data$lotus
lotus_met$metadata

# coare datasets
(coare_sal <- uopg(dataset = 'coare', type = "salinity"))
}
}
\references{
\url{http://uop.whoi.edu/archives/dataarchives.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcofi.R
\name{calcofi}
\alias{calcofi}
\title{Get California Cooperative Oceanic Fisheries Investigations (CALCOFI) data.}
\usage{
calcofi(dataset = "hydro_cast", path = "~/.ots/calcofi", overwrite = TRUE)
}
\arguments{
\item{dataset}{A dataset code name, one of hydro, macrozoo, or X}

\item{path}{A path to store the files, Default: \code{~/.ots/kelp}}

\item{overwrite}{(logical) To overwrite the path to store files in or not, Default: TRUE.}
}
\description{
FIXME: need to support other data sets in CALCOFI
}
\examples{
\dontrun{
# hydro cast data
(res <- calcofi('hydro_cast'))

# hydro bottle data
(res <- calcofi('hydro_bottle'))

library('dplyr')
res$data \%>\%
 tbl_df \%>\%
 select(-cruz_sta) \%>\%
 group_by(year) \%>\%
 summarise(avg_wind_spd = mean(wind_spd, na.rm = TRUE)) \%>\%
 arrange(desc(avg_wind_spd)) \%>\%
 head
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{type_sum}
\alias{type_sum}
\alias{type_sum.default}
\alias{type_sum.character}
\alias{type_sum.Date}
\alias{type_sum.factor}
\alias{type_sum.integer}
\alias{type_sum.logical}
\alias{type_sum.array}
\alias{type_sum.matrix}
\alias{type_sum.numeric}
\alias{type_sum.POSIXt}
\title{Type summary}
\usage{
type_sum(x)

\method{type_sum}{default}(x)

\method{type_sum}{character}(x)

\method{type_sum}{Date}(x)

\method{type_sum}{factor}(x)

\method{type_sum}{integer}(x)

\method{type_sum}{logical}(x)

\method{type_sum}{array}(x)

\method{type_sum}{matrix}(x)

\method{type_sum}{numeric}(x)

\method{type_sum}{POSIXt}(x)
}
\description{
Type summary
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ots-package.R
\docType{package}
\name{ots-package}
\alias{ots-package}
\alias{ots}
\title{ots}
\description{
ots, an R client for various ocean time series datasets.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bats.R
\name{bats}
\alias{bats}
\title{Get BATS (Bermuda Atlantic Time-series Study) data}
\usage{
bats(dataset = "production")
}
\arguments{
\item{dataset}{(character) One of production, zooplankton, or flux}
}
\description{
Get BATS (Bermuda Atlantic Time-series Study) data
}
\details{
Only these three data sets for now. Others may follow.
}
\examples{
\dontrun{
# Production
(res <- bats("production"))
res$meta
res$vars

# Zooplankton
(res <- bats("zooplankton"))
res$meta
res$vars

# Flux
(res <- bats("flux"))
res$meta
res$vars
}
}
\references{
\url{http://bats.bios.edu/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hot.R
\name{hot}
\alias{hot}
\title{Get Hawaii Ocean Time-series surface CO2 system (HOT) data.}
\usage{
hot()
}
\description{
Get Hawaii Ocean Time-series surface CO2 system (HOT) data.
}
\details{
\itemize{
\item cruise The integer number of the HOT cruise.
\item days The mid-day of the cruise expressed as the number of days from 1 October 1988.
\item date The mid-day of the cruise expressed as a calendar date.
\item temp The mean in situ seawater temperature, in °C.
\item sal The mean seawater salinity, in practical salinity units.
\item phos The mean seawater phosphate concentration, in µmol kg-1.
\item sil The mean seawater silicate concentration, in µmol kg-1.
\item DIC The mean seawater dissolved inorganic carbon (= total CO2) concentration, in µmol kg-1.
\item TA The mean seawater total alkalinity, in µeq kg-1.
\item nDIC The mean seawater salinity-normalized DIC, in µmol kg-1 at salinity = 35.
\item nTA The mean seawater salinity-normalized TA, in µeq kg-1 at salinity = 35.
\item pHmeas_25C The mean measured seawater pH at 25 °C, on the total scale.
\item pHmeas_insitu The mean measured seawater pH, adjusted to in situ temperature, on the total
scale.
\item pHcalc_25C The mean seawater pH, calculated from DIC and TA at 25 °C, on the total scale.
\item pHcalc_insitu The mean seawater pH, calculated from DIC and TA at in situ temperature, on
the total scale.
\item pCO2calc_insitu The mean seawater CO2 partial pressure, in µatm, calculated from DIC and
TA at in situ temperature
\item pCO2calc_20C The mean seawater CO2 partial pressure, in µatm, calculated from DIC and TA
at 20 °C.
\item aragsatcalc_insitu The mean seawater aragonite saturation state (solubility ratio).
\item calcsatcalc_insitu The mean seawater calcite saturation state (solubility ratio).
\item freeCO2_insitu The mean seawater free CO2 concentration, in µmol kg-1.
\item carbonate_insitu The mean seawater carbonate ion concentration, in µmol kg-1.
\item notes This column contains pertinent notes about the data used, especially where sampling
was different from the HOT standard procedures. The below codes are used to indicate
such variations.
}

notes options:

\itemize{
\item a No HOT DIC data; used Keeling DIC
\item b No HOT TA data; used Keeling TA
\item c No pH samples collected
\item d No TA data; assumed nTA of 2303.2
\item e No DIC data in upper 30 dbar; used value from 34.0 dbar
\item f No TA data in upper 30 dbar; used value from 51.3 dbar
\item g Used temp and sal data from downcast before cable broke
\item h No phos data; assumed 0.07
\item i No sil data; assumed 1.04
\item j Rosette lost; samples collected with Go-flo bottle on Kevlar line
\item k No sampling occurred
\item l No rosette; samples collected with Niskin bottles on Kevlar line
\item m No DIC or TA data in upper 30 dbar; used values from 34.9 dbar
\item n No DIC or TA data in upper 30 dbar; used values from 44.9 dbar
\item o No DIC, TA or pH sampling occurred in upper 200 dbar
\item p No DIC, TA or pH sampling occurred
}
}
\examples{
\dontrun{
(res <- hot())
res$meta
res$vars
res$citation
}
}
