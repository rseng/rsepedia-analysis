<!-- README.md is generated from README.Rmd. Please edit that file -->
RGPDD
=====

***DEPRECATED*** 

This package is no longer necessary, as the GPDD can now be downloaded directly as `.csv` files from the KNB repostitory: <https://doi.org/10.5063/F1BZ63Z8>


This package originally provided an R interface the [Global Population Dynamics Database](http://www3.imperial.ac.uk/cpb/databases/gpdd), where data was served from a no-longer-maintained API and provided only in the Microsoft Access DataBase format, both of which made access more difficult for R users.  

## Test environments
* local OS X install, R 3.1.2
* ubuntu 12.04 (on travis-ci), R 3.1.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'R6'

  R6 is a build-time dependency.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of httr 
(https://github.com/wch/checkresults/blob/master/httr/r-release). All packages 
that I could install passed except:

* XYZ:...
Data downloaded as MDB file.  
---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# RGPDD

An R interface the [Global Population Dynamics Database](http://www3.imperial.ac.uk/cpb/databases/gpdd)

*Package under development*. Please see the package [issues](/issues) and [milestones](/milestones) for more information.  


## Background and terms

Imperial College, London, Center for Population Biology provides the GPDD in Microsoft Access format available for download without cost after registration with an email address and password. This package aims to streamline the use of this data for R users.  At this time the package primarily provides the 9 data tables of the GPDD directly.  Helper functions for common operations should be added as time permits, and user contributions are welcome. 

Please consult [http://www3.imperial.ac.uk/cpb/databases/gpdd](http://www3.imperial.ac.uk/cpb/databases/gpdd) for more information and official documentation of the GPDD data.  Documentation from the GPDD is being included as part of the R package documentation where appropriate.


At the time of writing, the GPDD database was last updated in 2010 (v2), see `?gpdd_version`. The original (v1) was published in 1999.  Please consider registering for the GPDD and be sure to cite it appropriately:

>  "NERC Centre for Population Biology, Imperial College (2010) The Global Population Dynamics Database Version 2. http://www.sw.ic.ac.uk/cpb/cpb/gpdd.html". 

and notify the database administrator of any publications resulting from this data: cpb-gpdd-dl@imperial.ac.uk

The adminstrator has been notified about the existence of the R package project.


## Quickstart

Install the package:

```r
devtools::install_github("ropensci/rgpdd")
```

Load the data and explore the tables.  

```{r}
library("rgpdd")
library("ggplot2")
library("dplyr")
```

```{r}
ggplot(dplyr::filter(gpdd_data, MainID %in% 1:10)) + geom_line(aes(SeriesStep, Population, col=MainID, group=MainID))
```



---
title: "Vignette Title"
author: "Vignette Author"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

Vignettes are long form documentation commonly included in packages. Because they are part of the distribution of the package, they need to be as compact as possible. The `html_vignette` output type provides a custom style sheet (and tweaks some options) to ensure that the resulting html is as small as possible. The `html_vignette` format:

- Never uses retina figures
- Has a smaller default figure size
- Uses a custom CSS stylesheet instead of the default Twitter Bootstrap style

## Vignette Info

Note the various macros within the `vignette` setion of the metadata block above. These are required in order to instruct R how to build the vignette. Note that you should change the `title` field and the `\VignetteIndexEntry` to match the title of your vignette.

## Styles

The `html_vignette` template includes a basic CSS theme. To override this theme you can specify your own CSS in the document metadata as follows:

    output: 
      rmarkdown::html_vignette:
        css: mystyles.css

## Figures

The figure sizes have been customised so that you can easily put two images side-by-side. 

```{r, fig.show='hold'}
plot(1:10)
plot(10:1)
```

You can enable figure captions by `fig_caption: yes` in YAML:

    output:
      rmarkdown::html_vignette:
        fig_caption: yes

Then you can use the chunk option `fig.cap = "Your figure caption."` in **knitr**.

## More Examples

You can write math expressions, e.g. $Y = X\beta + \epsilon$, footnotes^[A footnote here.], and tables, e.g. using `knitr::kable()`.

```{r, echo=FALSE, results='asis'}
knitr::kable(head(mtcars, 10))
```

Also a quote using `>`:

> "He who gives up [code] safety for [code] speed deserves neither."
([via](https://twitter.com/hadleywickham/status/504368538874703872))
% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gpdd_timeperiod}
\alias{gpdd_timeperiod}
\title{The time period table}
\description{
TimePeriod is a look-up table that provides text descriptions of
the temporal period the sample relates to such as 'January',
'Spring', 'Week 1' and 'Day 1'.
}
\author{
GPDD Administrator \email{cpb-gpdd-dl@imperial.ac.uk}
}
\references{
\url{http://www3.imperial.ac.uk/cpb/databases/gpdd}
}
\keyword{data}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gpdd_main}
\alias{gpdd_main}
\title{main table: metadata for each time series}
\description{
A MAIN record is a 'time series' which is unique
Taxon/Location/LifeCycle combination. Sequential data for multiple
life stages (e.g. eggs, larve and adults) are split into different
Main records and must be amalgamated to create a single time series.
Where more than one adult generation occurs per year generation is
identified in the generation column of the data table
}
\author{
GPDD Administrator \email{cpb-gpdd-dl@imperial.ac.uk}
}
\references{
\url{http://www3.imperial.ac.uk/cpb/databases/gpdd}
}
\keyword{data}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gpdd_biotope}
\alias{gpdd_biotope}
\title{The biotype table}
\description{
The biotype table
}
\author{
GPDD Administrator \email{cpb-gpdd-dl@imperial.ac.uk}
}
\references{
\url{http://www3.imperial.ac.uk/cpb/databases/gpdd}
}
\keyword{data}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gpdd_data}
\alias{gpdd_data}
\title{The data table}
\description{
stores the individual time series abundance records
}
\author{
GPDD Administrator \email{cpb-gpdd-dl@imperial.ac.uk}
}
\references{
\url{http://www3.imperial.ac.uk/cpb/databases/gpdd}
}
\keyword{data}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gpdd_taxon}
\alias{gpdd_taxon}
\title{The taxon table}
\description{
The taxon table stores the taxonomic names relating to Main records.
It is links to the MAIN table with a unique TaxonID.  Most series
are for species.  Some extra information regarding breeding habitats
etc may be found in the notes column
}
\author{
GPDD Administrator \email{cpb-gpdd-dl@imperial.ac.uk}
}
\references{
\url{http://www3.imperial.ac.uk/cpb/databases/gpdd}
}
\keyword{data}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/rgpdd-package.R
\docType{package}
\name{rgpdd}
\alias{rgpdd}
\alias{rgpdd-package}
\title{rgpdd: an R interface to the Global Population Dynamics Database}
\description{
The GPDD was initially compiled by John Prendargast, Ellen Bazeley-White?, Owen Smith, John Lawton and Pablo Inchausti and released in 1999. Version 2.0 was released in 2010 following a substantial restructuring of the database and the addition of 123 new series by David Kidd and Sarah Knight.
}
\details{
Understanding the way in which populations of wild plants and animals behave over long periods of time is crucial to unravelling the way in which communities are assembled and the way in which they respond to disturbance, control or harvesting. The implications for conservation and agriculture are legion. Aside from practicalities, population variation is also intrinsically interesting, and provides a wealth of opportunity for mathematical innovation or exploration, especially when populations have particular cyclic, outbreaking or chaotic properties. For most students of population behaviour, the limiting factor in investigating any of these phenomena, and the development of theory to explain them, is the availability of suitable data. Usually, where analyses are performed and published, authors work on data that they have collected themselves. By definition, the collection of population time series is a lengthy process, and many ecologists have committed themselves to a lifetime of work in order to accumulate detailed information on populations at certain sites over many years.

Studies of population behaviour tend to address a number of themes, each with a typical taxonomic flavour. Thus, students of the chaotic vs cyclic question tend to focus on small mammals, those with an interest in the effects of harvesting or culling generally study fisheries or large mammal data respectively and analyses of insect populations tend to dominate the literature on pest control.

Analysis and subsequent publication can only take place once time series of adequate length have been amassed, but frequently authors will continue to collect data after publication and may follow the first paper with an updated or extended version, or a book or book chapter at a later date. There are examples of data sets that have been assiduously collected but from which no publications have resulted, or from which internal, private or unpublished documents have been generated.

The result of all this fragmentary activity in population dynamics, where data sets are often analysed individually, or in line with certain taxonomic conventions, is that it has been difficult to a) formulate general theory and b) investigate large scale pattern, both spatially and taxonomically. The general unavailability of data has also led to the development of theory through the repeated analysis of the same data sets. The celebrated Canadian lynx/snowshoe hare cycle has been the subject of analyses and publications almost too numerous to count. There is an obvious danger that if individual data sets such as this happen not to be representative of the way in which most populations behave then theoretical understanding may suffer.
}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gpdd_version}
\alias{gpdd_version}
\title{GPDD version information}
\description{
This is Version 2.0 (released 2010).
Version 1.0 (released in 1999) has now been superseded by v2.0 which includes the following enhancements,
}
\details{
-  A consistent definition of a time-series.
- Consistent metadata.
- Units.
- Sampling protocol.
- Consistent temporal coding.
- Addition of missing location information, the spatial bounds of study areas and a spatial accuracy index.
- Abundance data are supplied ‘retro-transformed’ as well as in the published source units.
- Improved documentation.
- 123 additional time-series are included, courtesy of Barry Brook (University of Adelaide).
- Removal of un-cited associated data including body size and biotope information.
}
\author{
GPDD Administrator \email{cpb-gpdd-dl@imperial.ac.uk}
}
\references{
\url{http://www3.imperial.ac.uk/cpb/databases/gpdd}
}
\keyword{data}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gpdd_restricted}
\alias{gpdd_restricted}
\title{Restricted data sets table}
\description{
Due to licensing restrictions 686 series from 6 sources cannot
be distributed without the permission of the owner. These are
data from the British Trust for Ornithology’s Common Bird Census
(97) and Constant Effort Recording Scheme (32), Rothamstead
Experimental Station, UK (9), the National monitoring programme
for wintering wildfowl in Norway 1980 – 93 (T. Nygard, 23),
Phalacrocorax carbo (Great cormorant) and Somateria mollissima
(Common eider) series supplied by N. Rov (2) and data from insect
light trapping supplied by H. Wolda (523).
}
\author{
GPDD Administrator \email{cpb-gpdd-dl@imperial.ac.uk}
}
\references{
\url{http://www3.imperial.ac.uk/cpb/databases/gpdd}
}
\keyword{data}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gpdd_datasource}
\alias{gpdd_datasource}
\title{The data source table}
\description{
Information on where the data was obtained from, relating to the
Main table through a unique DatasourceID.  Sources of data include
published journals, books and unpublished datasets and the references
details are held here.  The table also contains information regarding
access restrictions, contact details and in what medium the data was obtained.
}
\author{
GPDD Administrator \email{cpb-gpdd-dl@imperial.ac.uk}
}
\references{
\url{http://www3.imperial.ac.uk/cpb/databases/gpdd}
}
\keyword{data}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gpdd_location}
\alias{gpdd_location}
\title{Table of locations information for each timeseries}
\description{
Table of locations information for each timeseries
}
\author{
GPDD Administrator \email{cpb-gpdd-dl@imperial.ac.uk}
}
\references{
\url{http://www3.imperial.ac.uk/cpb/databases/gpdd}
}
\keyword{data}

