
<!-- README.md is generated from README.Rmd. Please edit that file -->

# skynet <img src="man/figures/logo.png" align="right" />

![Build Status](https://travis-ci.org/ropensci/skynet.svg?branch=master)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/skynet)](https://cran.r-project.org/package=skynet)
![](https://cranlogs.r-pkg.org/badges/skynet?color=brightgreen)
[![Coverage
status](https://codecov.io/gh/ropensci/skynet/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/Skynet?branch=master)
[![](https://badges.ropensci.org/214_status.svg)](https://github.com/ropensci/software-review/issues/214)

# Overview

The rationale behind Skynet, is to provide researchers with a unifying
tool overcoming some of the challenges faced when dealing with the
Bureau of Transport Statistics, DB1B and T100 data. The DB1B data
consists of 2 sets of files, Coupon and Ticket. They can be both
downloaded at
<https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=289> and
<https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=272>
respectively while the T100 data can be found here
<https://www.transtats.bts.gov/Tables.asp?DB_ID=111>.

## Note

To comply with R syntax guidelines, we changed to a clearer function
naming from version 1.2.0. Deprecated functions are still present, but
will be removed for the next versions.

## Note on importing from other data sources

We are constantly working on new functions that allow importing data
from different data sources. However, as we can’t cover them all at
least for now, in case you would like to work with a database which is
not covered by skynet, simply create a data.frame with the following
variables:

`itin_id, mkt_id, seq_num, origin_mkt_id, origin, year, quarter,
dest_mkt_id, dest, trip_break, op_carrier, distance, gateway, roundtrip,
itin_yield, passengers, itin_fare, bulk_fare, distance_full`

For more information on the variables, please visit
<https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=289> and
<https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=272>.

Skynet allows that some of this variables have a 0 or NA value, however,
if you’re working with a specific dataset which doesn’t allow an easy
conversion to our format, please feel free to create an issue so we can
look into it. Please make sure to include at least one small example of
a csv file with the data you’re trying to import.

## Installation

You can install skynet from github with:

``` r
# install.packages("devtools")
devtools::install_github("FilipeamTeixeira/skynet")
```

## Import Data

To import data, simply type `import_db1b()` or `import_t100()` including
the path to your desired file.  
**Note**: The Coupon file should take the first argument while the
Ticket file should take the second argument.

``` r
 library(skynet)
 import_db1b("folder/Coupon 2016Q1.csv", "folder/Ticket 2016Q1.csv")
 import_t100("folder/T100_2016.csv")
```

The BTS DB1B data consists of 2 sets of files, `Coupon` and `Ticket`.
They can be both downloaded at
<https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=289> and
<https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=272>
respectively.

Despite being possible to download the complete zipped file, which
includes all variables, due to its size, we recommend selecting the
following set.

| Coupon                     | Ticket             |
| :------------------------- | :----------------- |
| Itinerary ID               | Itinerary ID       |
| Market ID                  | Roundtrip          |
| Sequence Number            | Itinerary Yield    |
| Origin City Market ID      | Passengers         |
| Origin                     | Itinerary Fare     |
| Year                       | Bulkfare Indicator |
| Quarter                    | Distance           |
| Destination City Market ID |                    |
| Destination                |                    |
| Trip Break                 |                    |
| Operating Carrier          |                    |
| Distance                   |                    |
| Gateway                    |                    |

Since version 1.0.2 that the import method changed being the
`netimport()` function no longer available. When importing from the
prezipped DB1B file, just add the argument `zip = TRUE` to the
`import_db1b()` function. This does not apply to the T100 file which can
be simply imported by typing `import_t100()`. In order to save space, it
is possible as well to import the prezipped file, and convert it to a
smaller file with only the necessary variables, with the function
`convert_raw()`.

## Example

To generate a directed network, please type:

    library(skynet)
    # For DB1B data
    import_db1b("folder/Coupon_2011Q1.csv", "folder/Ticket_2011Q1.csv")
    make_net_dir(OD_2011Q1, disp = TRUE, alpha = 0.05)
    
    # For T100 data
    import_t100("folder/T100_2011.csv")
    make_net_dir(T100_2011Q1, disp = TRUE, alpha = 0.05)

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# skynet 0.9.0

* Added a `NEWS.md` file to track changes to the package.

# skynet 0.9.2

* Added tests
* Added Sample Data
* Fixed Imports on Description file

# skynet 0.9.3

* Fixed "no visible binding for global variable" issue
* Replaced disparity filter from semnet package by own package
* Small performace improvements
* Corrected spelling

# skynet 0.9.4

* Added T-100 import
* Corrected issue with general import for international option

# skynet 0.9.7

* Added new map function, now automatically printing different carriers with different colors.
* Improved import functions
* Importing from prezipped file, no longer requires extra function.

# skynet 0.9.7

* Changed way itin_fare was calculated for Directed, Undirected and Metro Networks. Now it uses price per mile and distance between stops to generate that info.

# skynet 0.9.9

* netImport now imports T100 market and segment files.
* netPath airlines renamed to carrier.
* updated vignettes.

# skynet 1.0

* netMetro has been replaced by argument in `netDir()` and `netUnd()`.

# skynet 1.0.1

* Possible to include carriers for undirected networks.
* Possible to filter non-scheduled flights.
* Ground Transport is now included as a carrier.
* Metro Network can be plotted.
* Improved way of calculating airport passenger frequency.
* Minor bug fixes.

# skynet 1.0.2

* New import functions. Now there are separate functions to import csv files from both DB1B and T100 databases.
* New bootnet function to bootstrap networks.

# skynet 1.0.3

* Minor adjustments
* Improved readability

# skynet 1.0.4

* Improved ReadMe file
* Fixed website
* Added extra comments and help information

# skynet 1.1.0

* Changed way files are imported. Now Coupon should take the first argument and Ticket the second.
* Minor adjustmenst to the help files.

# skynet 1.2.0

* Major function naming changes to match syntax etiquette
* Skynet S3 class added

# skynet 1.2.1

* Year and quarter added to skynet object
* Fixed site

# skynet 1.2.2

* Removed convert_raw as it is easier to import using the zip = TRUE argument and select the format to be saved to.

# skynet 1.2.3

* Now it is possible to import directly files from the BTS website with the download_db1b() function.

# skynet 1.3

* We have fixed some bugs and added the download_t100 function, making it finally possible to import both T100 and DB1B datasets without having to navigate to the BTS website.

# skynet 1.3.2

* We've added the possibility to import On-time performance data.

# skynet 1.3.6

* Added dplyr 1.0.0 compatibility

# skynet 1.3.7

* Improved test time

# skynet 1.3.9

* Fixed links
* Added download failed message

# skynet 1.4.1

* Fixed issue where download_t100 requires query to be encoded.
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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
## Release summary

* Complies with CRAN policies on internet resources.
* Replaced error in download_db1b, download_t100, download_t100int and
download_ontime, with message.

## Test environments
* local OS X install, R 4.0.1
* ubuntu 12.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

* I have run R CMD check on the NUMBER downstream dependencies.

* All revdep maintainers were notified of the release on RELEASE DATE.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# skynet <img src="man/figures/logo.png" align="right" />

![Build Status](https://travis-ci.org/ropensci/skynet.svg?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/skynet)](https://cran.r-project.org/package=skynet)
![](https://cranlogs.r-pkg.org/badges/skynet?color=brightgreen)
[![Coverage status](https://codecov.io/gh/ropensci/skynet/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/Skynet?branch=master)
[![](https://badges.ropensci.org/214_status.svg)](https://github.com/ropensci/software-review/issues/214)

# Overview

The rationale behind Skynet, is to provide researchers with a unifying tool overcoming some of the challenges faced when dealing with the Bureau of Transport Statistics, DB1B and T100 data. 
The DB1B data consists of 2 sets of files, Coupon and Ticket. They can be both downloaded at https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=289 and https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=272 respectively while the T100 data can be found here https://www.transtats.bts.gov/Tables.asp?DB_ID=111.

## Note

To comply with R syntax guidelines, we changed to a clearer function naming from version 1.2.0. Deprecated functions are still present, but will be removed for the next versions.

## Note on importing from other data sources

We are constantly working on new functions that allow importing data from different data sources. However, as we can't cover them all at least for now, in case you would like to work with a database which is not covered by skynet, simply create a data.frame with the following variables: 

`itin_id, mkt_id, seq_num, origin_mkt_id, origin, year, quarter, dest_mkt_id, dest, trip_break, op_carrier, distance, gateway, roundtrip, itin_yield, passengers, itin_fare, bulk_fare, distance_full` 

For more information on the variables, please visit https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=289 and https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=272.

Skynet allows that some of this variables have a 0 or NA value, however, if you're working with a specific dataset which doesn't allow an easy conversion to our format, please feel free to create an issue so we can look into it. Please make sure to include at least one small example of a csv file with the data you're trying to import.

## Installation

You can install skynet from github with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("FilipeamTeixeira/skynet")
```

## Import Data

To import data, simply type `import_db1b()` or `import_t100()` including the path to your desired file.  
**Note**: The Coupon file should take the first argument while the Ticket file should take the second argument.
    
```{r, eval=FALSE}
 library(skynet)
 import_db1b("folder/Coupon 2016Q1.csv", "folder/Ticket 2016Q1.csv")
 import_t100("folder/T100_2016.csv")
```

The BTS DB1B data consists of 2 sets of files, `Coupon` and `Ticket`.
They can be both downloaded at https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=289 and https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=272 respectively.

Despite being possible to download the complete zipped file, which includes all variables, due to its size, we recommend selecting the following set. 


```{r, echo=FALSE, message=FALSE, warning=FALSE, results='asis'}
knitr::kable(data.frame(Coupon = c("Itinerary ID", "Market ID", "Sequence Number", "Origin City Market ID",  
"Origin", "Year", "Quarter", "Destination City Market ID", "Destination", "Trip Break", "Operating Carrier", 
"Distance", "Gateway"),
Ticket = c("Itinerary ID", "Roundtrip", "Itinerary Yield", "Passengers",
"Itinerary Fare", "Bulkfare Indicator", "Distance","","","","","","")))
```

Since version 1.0.2 that the import method changed being the `netimport()` function no longer available.
When importing from the prezipped DB1B file, just add the argument `zip = TRUE` to the `import_db1b()` function. This does not apply to the T100 file which can be simply imported by typing `import_t100()`.
In order to save space, it is possible as well to import the prezipped file, and convert it to a smaller file with only the necessary variables, with the function `convert_raw()`. 


## Example

To generate a directed network, please type:

    library(skynet)
    # For DB1B data
    import_db1b("folder/Coupon_2011Q1.csv", "folder/Ticket_2011Q1.csv")
    make_net_dir(OD_2011Q1, disp = TRUE, alpha = 0.05)

    # For T100 data
    import_t100("folder/T100_2011.csv")
    make_net_dir(T100_2011Q1, disp = TRUE, alpha = 0.05)


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

---
title: "Introduction to SKYNET"
author: "Filipe Teixeira"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to SKYNET}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

SKYNET is a flexible R package that allows generating bespoke air transport statistics for urban studies based on publicly available data from the Bureau of Transport Statistics (BTS) in the United States. 

## SKYNET's segments

SKYNET is effectively divided into four segments:

1. Import Data
1. Generate Air Networks
1. Plot Air Networks


## Import Data

To import data, simply type `import_db1b()` or `import_t100()` including the path to your desired file. 
Note: we recommend naming the files with a similar structure as `Ticket 2016Q1.csv` or `Coupon 2016Q1.csv` respectively.
    
```{r, eval=FALSE}
 library(skynet)
 import_db1b("folder/Coupon 2016Q1.csv", "folder/Ticket 2016Q1.csv")
 import_t100("folder/T100_2016.csv")
```

The BTS DB1B data consists of 2 sets of files, `Coupon` and `Ticket`.
They can be both downloaded at https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=289 and https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=272 respectively.

Despite being possible to download the complete zipped file, which includes all variables, due to its size, we recommend selecting the following set. 


```{r, echo=FALSE, message=FALSE, warning=FALSE, results='asis'}
knitr::kable(data.frame(Coupon = c("Itinerary ID", "Market ID", "Sequence Number", "Origin City Market ID",  
"Origin", "Year", "Quarter", "Destination City Market ID", "Destination", "Trip Break", "Operating Carrier", 
"Distance", "Gateway"),
Ticket = c("Itinerary ID", "Roundtrip", "Itinerary Yield", "Passengers",
"Itinerary Fare", "Bulkfare Indicator", "Distance Full","","","","","","")))
```

Since version 1.0.2 that the import method changed being the `netimport()` function no longer available.
When importing from the prezipped DB1B file, just add the argument `zip = TRUE` to the `import_db1b()` function. This does not apply to the T100 file which can be simply imported by typing `import_t100()`.
In order to save space, it is possible as well to import the prezipped file, and convert it to a smaller file with only the necessary variables, with the function `convert_raw()`. 

When importing files from the T100 dataset, we recommend naming the file as `T100 year mkt` for the Market dataset and `T100 year seg` for the Segment dataset.


## Create networks

SKYNET creates three types of networks and an extra option:

1. Directed Network - `make_net_dir()`
1. Undirected Network - `make_net_und()`
1. Path Network - `make_net_path()`
1. Metro Network - To be used as argument in `make_net_dir()` and `make_net_und()`
1. International Option - `make_net_int()`

When generating a network, SKYNET, creates a list which includes:

1. Dataframe from original data (example below)
1. iGraph object
1. Dataframe with nodes(airports) statistics


```{r, echo=FALSE, message=FALSE, warning=FALSE, results='asis'}
library(skynet)
library(dplyr)
library(kableExtra)
options(knitr.table.format = "html") 
data("OD_Sample")
rownames(OD_Sample) <- NULL
knitr::kable(head(OD_Sample, 5)) %>% kable_styling()
```

When generating a network with SKYNET, it is possible use the following arguments:  
1. `carriers` - groups OD data per carrier when `TRUE`  

To extract the backbone of the network:  
1. `cap` (to be used with `pct`) - filters the network based on a given percentage (default percentage = 10%)  
1. `disp` (to be used with `alpha`) - filters the network using the Serrano et all backbone extraction algorithm (default alpha = 0.003)  

## Create Maps

One of SKYNET's advantages is the possibility of plotting maps without having to recur to external software.

Typing `net_map(skynet_object)` plots a ggplot2 based map with OD information. When specifying the group by carrier option when generating a network, `net_map()` distinguishes carriers with different colors. The `pct` argument allow to plot only a percentage of the available data.
It is important to point the path to the dataframe created by SKYNET.

```{r, echo=FALSE, message=FALSE, warning=FALSE,dpi = 300, fig.width = 6, fig.height= 4, out.width=500, }
library(skynet)
data("OD_Sample")
test <- make_net_dir(OD_Sample)
net_map(test, pct = 10)
```

## Extra Functions

SKYNET, allows as well to perform quick searches on both airports and carriers, by their IATA code. `find_airport()`, `find_carrier()`. 

## Bootnet

With version 1.0.2, we included the option to bootstrap networks and retrieve certain network statistics.

```{r}
library(skynet)
test <- make_net_dir(OD_Sample)
boot_network(test$gDir, n = 10)
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMap.R
\name{net_map}
\alias{net_map}
\title{Plot Skynet}
\usage{
net_map(x, pct = 60)
}
\arguments{
\item{x}{Skynet Object
(generated by make_net_dir,make_net_und or make_net_path)}

\item{pct}{percentage of edges to include}
}
\description{
Creates OD ggplot2 generated maps from make.net functions
Shows sample of 60% of flights
}
\examples{
\dontrun{
network <- make.netDir(OD_Sample)
net_map(network, pct = 10)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netPath.R
\name{make_net_path}
\alias{make_net_path}
\title{Path and OD Network}
\usage{
make_net_path(x, leg = FALSE, zero = FALSE, carrier = FALSE)
}
\arguments{
\item{x}{Data frame}

\item{leg}{Generates Leg Count Data frame, based on Path taken.}

\item{zero}{Displays percentage of 0 usd tickets}

\item{carrier}{Groups data per airline

For example, all passengers doing the BOS-ATL-LAX path, are summed by Air Carrier.}
}
\description{
Generates an OD network and a Leg Count data frame(on demand)
}
\examples{
\dontrun{
make_net_path(OD_Sample)

# Generate Leg Count
make_net_path(OD_Sample, leg = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_db1b.R
\name{import_db1b}
\alias{import_db1b}
\title{Import Data from DB1B files}
\usage{
import_db1b(c, t, zip = FALSE, auto = TRUE)
}
\arguments{
\item{c}{Coupon csv file to be imported, in case of DB1B database}

\item{t}{Ticket csv file to be imported, in case of DB1B database}

\item{zip}{Should equal TRUE if original file comes from the BTS prezipped option.}

\item{auto}{Automatically assigns object}
}
\description{
Imports data from BTS/RITA/Transtats files
}
\details{
Coupon files can be found at \url{https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=289}.
Ticket files can be found at \url{https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=272}.
Both files should belong to the same year and same quarter.
\strong{Note}: We do recommend sparklyr to be used for larger sets of data.
More information on variables to select and type of files to use can be found \href{https://github.com/ropensci/skynet}{here}
}
\examples{
\dontrun{

import_db1b(skynet_example("Coupon_2001Q1.csv"), skynet_example("Ticket_2001Q1.csv"))

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disparityfilter.R
\name{disparity_filter}
\alias{disparity_filter}
\title{Disparity Filter}
\usage{
disparity_filter(g, alpha = 0.003)
}
\arguments{
\item{g}{igraph object}

\item{alpha}{Alpha value. Default 0.003}
}
\description{
Uses the Serrano's disparity filter (\url{https://en.wikipedia.org/wiki/Disparity_filter_algorithm_of_weighted_network})
to extract the backbone of the network in "Extracting the multiscale backbone of complex weighted networks"
}
\examples{
\dontrun{
netDir <- make.netDir(OD_Sample)
disparity_filter(netDir$gDir, alpha = 0.003)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netImport.R
\name{netImport}
\alias{netImport}
\title{Import Data}
\usage{
netImport(x = NULL, y = NULL)
}
\arguments{
\item{x}{First csv file to be imported, in case of DB1B database, or in case of using
the T-100 database, the only file to be included.}

\item{y}{Second csv file to be imported.}
}
\description{
Imports data from BTS/RITA/Transtats website
File order doesn't matter, but it is recommended to name the files using the following
syntax: \emph{"Coupon YearQuarter.csv", "Ticket YearQuarter.csv", "T100 Year".}
Note: We do recommend sparklyr to be used for larger sets of data.
}
\examples{
\dontrun{

netImport(skynet_example("Coupon_2001Q1.csv"), skynet_example("Ticket_2001Q1.csv"))

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_t100.R
\name{import_t100}
\alias{import_t100}
\title{Import T-100 Data}
\usage{
import_t100(x, nonsch = FALSE, auto = TRUE)
}
\arguments{
\item{x}{T-100 csv}

\item{nonsch}{Should equal TRUE to include non-scheduled flights}

\item{auto}{Automatically assigns object}
}
\description{
Imports T-100 Data directly from BTS/RITA/Transtats website raw data (prezipped file),
for SKYNET's import function.
}
\details{
Files can be found here \url{https://www.transtats.bts.gov/Tables.asp?DB_ID=111}.
More information on variables to select and type of files to use can be found \href{https://github.com/ropensci/skynet}{here}
}
\examples{
\dontrun{

import_t100(skynet_example("T100_2011_mkt.csv"))

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_t100.R
\name{download_t100}
\alias{download_t100}
\title{Download Data from T100 files}
\usage{
download_t100(y = NULL, type = NULL)
}
\arguments{
\item{y}{year to be imported}

\item{type}{"mkt" for Market, "seg" for Segment databases respectively}
}
\description{
Downloads data from BTS/RITA/Transtats and imports it into R
}
\details{
Note: The BTS often changes the way we can access these files. So please be warned that this is still an experimental feature.
}
\examples{
\dontrun{

download_t100(2010, "mkt")

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_airports.R
\name{airportCodeFull}
\alias{airportCodeFull}
\title{Airport Data - full}
\format{
A dataframe with 6435 observations and 9 variables
}
\description{
USA airport data from the RITA/Transtats database
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createNodes.R
\name{nodeStatsMetro}
\alias{nodeStatsMetro}
\title{Create Metro Nodes}
\usage{
nodeStatsMetro(y)
}
\arguments{
\item{y}{Data Frame}
}
\description{
Create Metro Nodes
}
\examples{
\dontrun{

nodeStatsMetro(OD_Sample)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findCarrier.R
\name{find_carrier}
\alias{find_carrier}
\title{Find Carrier function}
\usage{
find_carrier(x)
}
\arguments{
\item{x}{Carrier}
}
\description{
Searches for airport information based on its IATA code or city name
}
\examples{
\dontrun{
find_carrier("United")

find_carrier("UA")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fromto.R
\name{from_to_stats}
\alias{from_to_stats}
\title{From To function}
\usage{
from_to_stats(x, y, orig)
}
\arguments{
\item{x}{igraph object to query}

\item{y}{origin airport IATA code}

\item{orig}{"from" or "to" options}
}
\description{
Calculate edges weight from IATA Code
}
\examples{
\dontrun{
netDir <- make.netDir(OD_Sample)
from_to_stats(netDir$gDir, "JFK", orig = "from")

from_to_stats(netDir$gDir, "JFK", orig = "to")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_db1b.R
\name{download_db1b}
\alias{download_db1b}
\title{Download Data from DB1B files}
\usage{
download_db1b(y = NULL, q = NULL)
}
\arguments{
\item{y}{year to be imported}

\item{q}{quarter to be imported}
}
\description{
Downloads data from BTS/RITA/Transtats and imports it into R
}
\details{
Coupon files are downloaded from \url{https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=289}.
Ticket files are downloaded from  \url{https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=272}.

Note: The BTS often changes the way we can access these files. So please be warned that this is still an experimental feature.
}
\examples{
\dontrun{

download_db1b(2010, 1)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createNodes.R
\name{createNodes}
\alias{createNodes}
\title{Create Nodes}
\usage{
createNodes(y)
}
\arguments{
\item{y}{Data Frame}
}
\description{
Creates nodes for SKYNET's functions.
Despite being possible to use it individually, it's mainly meant to be used as a complimentary function.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PowerLaw.R
\name{fit_power}
\alias{fit_power}
\title{Power Law}
\usage{
fit_power(graph)
}
\arguments{
\item{graph}{iGraph object}
}
\description{
Plots power law fit
}
\examples{
\dontrun{
netDir <- make.netDir(OD_Sample)
fit_power(netDir$gDir)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_metro.R
\name{MetroFull}
\alias{MetroFull}
\title{Metro (Full) Data}
\format{
A dataframe with 5802 observations and 5 variables
}
\description{
This data comes from the RITA/Transtats database
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netUnd.R
\name{make_net_und}
\alias{make_net_und}
\title{Undirected Network}
\usage{
make_net_und(
  x,
  disp = FALSE,
  alpha = 0.003,
  cap = FALSE,
  pct = 10,
  merge = TRUE,
  carrier = FALSE,
  metro = FALSE
)
}
\arguments{
\item{x}{Data frame}

\item{disp}{Uses the Serrano's disparity filter (\url{https://en.wikipedia.org/wiki/Disparity_filter_algorithm_of_weighted_network})
to extract the backbone of the network.}

\item{alpha}{Argument for disparity filter.}

\item{cap}{Filters original data based on the edge weight.}

\item{pct}{Argument for cap filter. Value should be imput as percentage.}

\item{merge}{When set to FALSE, it keeps parallel edges instead of collapsing them
and summing their weights.}

\item{carrier}{Groups data per carrier and OD}

\item{metro}{Groups data by metropolitan area}
}
\description{
Generates Undirected Network with an iGraph \strong{gUnd} object,
a Data Frame \strong{netUnd} and a Data Frame
with Airport/Nodes statistics \strong{nodes}.
}
\examples{
\dontrun{
make_net_und(OD_Sample)

# Apply Disparity Filter
make_net_und(OD_Sample, disp = TRUE, alpha = 0.05)

# Apply Percentage Cap
make_net_und(OD_Sample, cap = TRUE, pct = 20)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_airports.R
\name{airportCode}
\alias{airportCode}
\title{Airport Data - clean}
\format{
A dataframe with 6435 observations and 5 variables
}
\description{
USA airport data from the RITA/Transtats database
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_ontime.R
\name{download_ontime}
\alias{download_ontime}
\title{Download On-Time}
\usage{
download_ontime(y, m, auto = TRUE)
}
\arguments{
\item{y}{year to be imported}

\item{m}{month to be imported}

\item{auto}{Automatically assigns object}
}
\description{
Download On-Time Performance Data directly from BTS/RITA/Transtats website raw data (prezipped file),
for SKYNET's import function.
}
\examples{
\dontrun{

import_ontime(skynet_example("Ontime.csv"))

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nodeStats.R
\name{node_stats}
\alias{node_stats}
\title{Get node info}
\usage{
node_stats(x)
}
\arguments{
\item{x}{Data Frame to extract information from}
}
\description{
Creates node statistics
Generates Number of Passenger Arrivals, Departures and Transfers
}
\examples{
\dontrun{

node_stats(OD_Sample)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netInt.R
\name{make.netInt}
\alias{make.netInt}
\title{International Data}
\usage{
make.netInt(x = NULL, m = NULL, Q = NULL)
}
\arguments{
\item{x}{T-100 International Segment csv file}

\item{m}{Data set to merge with}

\item{Q}{Desired T-100 Quarter. Should be equal to 1, 2, 3 or 4.}
}
\description{
Imports International data to complement to the DB1B data set.
NOTE: When using this function, certain variables will be skewed as the T100 dataset does not contain
all the data the DB1B dataset contains.
}
\examples{
\dontrun{

make.netInt(skynet_example("T100_2011_int.csv"), OD_Sample, 1)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.skynet}
\alias{summary.skynet}
\title{Displays a summary of a skynet object}
\usage{
\method{summary}{skynet}(object, ...)
}
\arguments{
\item{object}{skynet object to summarise}

\item{...}{other arguments ignored (for compatibility with generic)}
}
\description{
Displays a summary of a skynet object
}
\examples{
net <- make_net_dir(OD_Sample)
summary(net)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_ODsample.R
\name{OD_Sample}
\alias{OD_Sample}
\title{Sample OD data}
\format{
A dataframe with 500.000 observations and 19 variables
}
\description{
Sample data to use with SKYNET functions
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_carriers.R
\name{aircraft_type}
\alias{aircraft_type}
\title{Aircraft type data}
\format{
A dataframe with 422 observations and 2 variables
}
\description{
This data comes from the RITA/Transtats database
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_carriers.R
\name{carriers}
\alias{carriers}
\title{Carrier data}
\format{
A dataframe with 1738 observations and 2 variables
}
\description{
This data comes from the RITA/Transtats database
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findAirport.R
\name{find_airport}
\alias{find_airport}
\title{Find Airport function}
\usage{
find_airport(x)
}
\arguments{
\item{x}{airport IATA code or city name}
}
\description{
Searches for airport information based on its IATA code or city name
It will display multiple airports as it works with partial names.
}
\examples{
\dontrun{
find_airport("Atlanta")

find_airport("ATL")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skynet_example.R
\name{skynet_example}
\alias{skynet_example}
\title{Get path to skynet examples}
\usage{
skynet_example(path = NULL)
}
\arguments{
\item{path}{File name.}
}
\description{
To access csv examples from SKYNET
}
\examples{
\dontrun{
skynet_example()
skynet_example("Coupon 2001Q1.csv")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_t100_int.R
\name{download_t100_int}
\alias{download_t100_int}
\title{Download Data from T100 international files}
\usage{
download_t100_int(y = NULL, type = NULL)
}
\arguments{
\item{y}{year to be imported}

\item{type}{"mkt" for Market, "seg" for Segment databases respectively}
}
\description{
Downloads data from BTS/RITA/Transtats and imports it into R
}
\details{
Note: The BTS often changes the way we can access these files. So please be warned that this is still an experimental feature.
}
\examples{
\dontrun{

download_t100_int(2010, "mkt")

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_net_trip.R
\name{make_net_trip}
\alias{make_net_trip}
\title{Trip directed network}
\usage{
make_net_trip(x, carrier = FALSE)
}
\arguments{
\item{x}{Data frame}

\item{carrier}{Groups data per carrier and OD}
}
\description{
Generates Trip/Route based Directed Network with an iGraph \strong{gDir} object,
a Data Frame \strong{netDir} and a Data Frame
with Airport/Nodes statistics \strong{nodes}.
Returns type of trip:
OD = Origin/Destination pair,
OT = Origin/Transfer pair,
TT = Transfer/Transfer pair,
TD = Transfer/Destination pair
}
\examples{
\dontrun{
make_net_trip(OD_Sample)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootnet.R
\name{boot_network}
\alias{boot_network}
\title{Network bootstrapping}
\usage{
boot_network(g, n = 500, left_ci = 0.005, right_ci = 0.995)
}
\arguments{
\item{g}{iGraph graph or skynet object.}

\item{n}{Number of bootstraps to run. (500 default)}

\item{left_ci}{Confidence interval left limit. (0.005 default)}

\item{right_ci}{Confidence interval left limit (0.995 default)}
}
\description{
Bootstraps a network and returns output containing three network statistics:
Average Path Length, Transitivity, Mean Betweenness.
}
\examples{
\dontrun{
boot_net(g, n = 500)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skynet.R
\docType{package}
\name{skynet-package}
\alias{skynet}
\alias{skynet-package}
\title{skynet: Network analysis for BTS Data}
\description{
Creates networks from the BTS/Transtats data
}
\details{
Given the DB1BCoupon and DB1BTicket, or the T-100 csv's exported
this package allows creating data frames and subsequent igraph graphs.
}
\examples{
NA

}
\references{
NA
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/skynet}
  \item Report bugs at \url{https://github.com/ropensci/skynet/issues}
}

}
\author{
Filipe Teixeira
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_ontime.R
\name{import_ontime}
\alias{import_ontime}
\title{Import on-time Data}
\usage{
import_ontime(x, auto = TRUE)
}
\arguments{
\item{x}{On-time csv (from zipped file)}

\item{auto}{Automatically assigns object}
}
\description{
Imports on-time Data directly from BTS/RITA/Transtats website raw data (prezipped file),
for SKYNET's import function.
}
\details{
Files can be found here \url{https://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=236}.
More information on variables to select and type of files to use can be found \href{https://github.com/ropensci/skynet}{here}
}
\examples{
\dontrun{

import_ontime(skynet_example("Ontime_2011_1.csv"))

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_metro.R
\name{MetroLookup}
\alias{MetroLookup}
\title{Metro Data}
\format{
A dataframe with 5802 observations and 2 variables
}
\description{
This data comes from the RITA/Transtats database
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netDir.R
\name{make_net_dir}
\alias{make_net_dir}
\title{Directed network}
\usage{
make_net_dir(
  x,
  disp = FALSE,
  alpha = 0.003,
  cap = FALSE,
  pct = 10,
  carrier = FALSE,
  metro = FALSE
)
}
\arguments{
\item{x}{Data frame}

\item{disp}{Uses the Serrano's disparity filter (\url{https://en.wikipedia.org/wiki/Disparity_filter_algorithm_of_weighted_network})
to extract the backbone of the network.}

\item{alpha}{Argument for disparity filter.}

\item{cap}{Filters original data based on the edge weight.}

\item{pct}{Argument for cap filter. Value should be imput as percentage.}

\item{carrier}{Groups data per carrier and OD}

\item{metro}{Groups data by metropolitan area}
}
\description{
Generates Directed Network with an iGraph \strong{gDir} object,
a Data Frame \strong{netDir} and a Data Frame
with Airport/Nodes statistics \strong{nodes}.
}
\examples{
\dontrun{
make_net_dir(OD_Sample)

# Apply Disparity Filter
make_net_dir(OD_Sample, disp = TRUE, alpha = 0.05)

# Apply Percentage Cap
make_net_dir(OD_Sample, cap = TRUE, pct = 20)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_airports.R
\name{airportMaster}
\alias{airportMaster}
\title{Airport Data - master}
\format{
A dataframe with 13555 observations and 28 variables
}
\description{
World airport data from the RITA/Transtats database
}
