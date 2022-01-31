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

# roadoi - Use Unpaywall with R

[![R build
status](https://github.com/ropensci/roadoi/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/roadoi/actions)
[![codecov.io](https://codecov.io/github/ropensci/roadoi/coverage.svg?branch=master)](https://codecov.io/github/ropensci/roadoi?branch=master)
[![cran
version](http://www.r-pkg.org/badges/version/roadoi)](https://cran.r-project.org/package=roadoi)
[![rstudio mirror
downloads](http://cranlogs.r-pkg.org/badges/roadoi)](https://github.com/r-hub/cranlogs.app)
[![review](https://badges.ropensci.org/115_status.svg)](https://github.com/ropensci/software-review/issues/115)

roadoi interacts with the [Unpaywall REST
API](https://unpaywall.org/products/api), an openly available
web-interface which returns metadata about open access versions of
scholarly works.

This client supports the most recent API Version 2.

API Documentation: <https://unpaywall.org/products/api>

## How do I use it?

Use the `oadoi_fetch()` function in this package to get open access
status information and full-text links from Unpaywall.

``` r
roadoi::oadoi_fetch(dois = c("10.1038/ng.3260", "10.1093/nar/gkr1047"), 
                    email = "najko.jahn@gmail.com")
#> # A tibble: 2 x 21
#>   doi      best_oa_location  oa_locations oa_locations_emb…
#>   <chr>    <list>            <list>       <list>           
#> 1 10.1038… <tibble [1 × 11]> <tibble [1 … <tibble [0 × 0]> 
#> 2 10.1093… <tibble [1 × 10]> <tibble [6 … <tibble [0 × 0]> 
#> # … with 17 more variables: data_standard <int>,
#> #   is_oa <lgl>, is_paratext <lgl>, genre <chr>,
#> #   oa_status <chr>, has_repository_copy <lgl>,
#> #   journal_is_oa <lgl>, journal_is_in_doaj <lgl>,
#> #   journal_issns <chr>, journal_issn_l <chr>,
#> #   journal_name <chr>, publisher <chr>,
#> #   published_date <chr>, year <chr>, title <chr>,
#> #   updated_resource <chr>, authors <list>
```

There are no API restrictions. However, providing an email address is
required and a rate limit of 100k is suggested. If you need to access
more data, use the [data dump](https://unpaywall.org/products/snapshot)
instead.

### RStudio Addin

This package also has a RStudio Addin for easily finding free full-texts
in RStudio.

![](man/figures/oadoi_addin.gif)

## How do I get it?

Install and load from [CRAN](https://cran.r-project.org/package=roadoi):

``` r
install.packages("roadoi")
library(roadoi)
```

To install the development version, use the [devtools
package](https://cran.r-project.org/package=devtools)

``` r
devtools::install_github("ropensci/roadoi")
library(roadoi)
```

## Documentation

See <https://docs.ropensci.org/roadoi/> to get started.

## Meta

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/roadoi/blob/master/CONDUCT.md). By
participating in this project you agree to abide by its terms.

License: MIT

Please use the [issue
tracker](https://github.com/ropensci/roadoi/issues) for bug reporting
and feature requests.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

## Test environments

- local Mac OS BigSur install (11.4), R version 4.1.0 (2021-05-18)
- GitHub Actions: windows-latest (release), macOS-latest (release), ubuntu-20.04 (release), ubuntu-20-04 (devel)
- CRAN Win Builder

## R CMD check results

On local machine (Mac OS): OK

win-builder (R-release, R-oldrelease, R-devel): OK

## Reverse dependencies

I have run R CMD check on downstream dependencies using revdepcheck::revdep_check(num_workers = 4) and found no problems related to this submission.

---

Dear CRAN maintainers, 

This submission fixes the problems shown on
<https://cran.r-project.org/web/checks/check_results_roadoi.html>.

Thanks for alerting me!

Najko Jahn
# roadoi - Use Unpaywall with R

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![R build status](https://github.com/ropensci/roadoi/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/roadoi/actions)
[![codecov.io](https://codecov.io/github/ropensci/roadoi/coverage.svg?branch=master)](https://codecov.io/github/ropensci/roadoi?branch=master)
[![cran version](http://www.r-pkg.org/badges/version/roadoi)](https://cran.r-project.org/package=roadoi)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/roadoi)](https://github.com/r-hub/cranlogs.app)
[![review](https://badges.ropensci.org/115_status.svg)](https://github.com/ropensci/software-review/issues/115)



roadoi interacts with the [Unpaywall REST API](https://unpaywall.org/products/api), 
an openly available web-interface which returns metadata about open access versions of scholarly works. 

This client supports the most recent API Version 2.

API Documentation: <https://unpaywall.org/products/api>

## How do I use it? 

Use the `oadoi_fetch()` function in this package to get open access status
information and full-text links from Unpaywall.


```{r}
roadoi::oadoi_fetch(dois = c("10.1038/ng.3260", "10.1093/nar/gkr1047"), 
                    email = "najko.jahn@gmail.com")
```

There are no API restrictions. However, providing an email address is required and a rate limit of 100k is suggested. If you need to access more data, use the [data dump](https://unpaywall.org/products/snapshot) instead.

### RStudio Addin

This package also has a RStudio Addin for easily finding free full-texts in RStudio.

![](man/figures/oadoi_addin.gif)

## How do I get it? 

Install and load from [CRAN](https://cran.r-project.org/package=roadoi):

```{r eval=FALSE}
install.packages("roadoi")
library(roadoi)
```

To install the development version, use the [devtools package](https://cran.r-project.org/package=devtools)

```{r eval = FALSE}
devtools::install_github("ropensci/roadoi")
library(roadoi)
```

## Documentation

See <https://docs.ropensci.org/roadoi/> to get started.

## Meta

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/roadoi/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

License: MIT

Please use the [issue tracker](https://github.com/ropensci/roadoi/issues) for bug reporting and feature requests.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Introduction"
author: "Najko Jahn"
date: "`r Sys.Date()`"
vignette: >
  %\VignetteIndexEntry{Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

### About Unpaywall

[Unpaywall](https://unpaywall.org/), developed and maintained by the [team of OurResearch](https://ourresearch.org/team/about), is a non-profit service that finds open access copies of scholarly literature by looking up a DOI (Digital Object Identifier). It not only returns open access full-text links, but also helpful metadata about the open access status of a publication such as licensing or provenance information.

Unpaywall uses different data sources to find open access full-texts including:

- [Crossref](https://www.crossref.org/): a DOI registration agency serving major scholarly publishers.
- [Directory of Open Access Journals (DOAJ)](https://doaj.org/): a registry of open access journals
- Various OAI-PMH metadata sources. OAI-PMH is a protocol often used by open access journals and repositories such as arXiv and PubMed Central.

See [Piwowar et al. (2018)](https://doi.org/10.7717/peerj.4375) for a comprehensive overview of Unpaywall.

### Basic usage

There is one major function to talk with Unpaywall, `oadoi_fetch()`, taking a character vector of DOIs and your email address as required arguments.

```{r}
library(roadoi)
roadoi::oadoi_fetch(dois = c("10.1186/s12864-016-2566-9",
                             "10.1103/physreve.88.012814"), 
                    email = "najko.jahn@gmail.com")
```

#### What's returned?

The client supports API version 2. According to the [Unpaywall Data Format](https://unpaywall.org/data-format), the following variables with the following definitions are returned:

**Column**|**Description**
|:------------|:----------------------------------------------
`doi`|DOI (always in lowercase)
`best_oa_location`|list-column describing the best OA location. Algorithm prioritizes publisher hosted content (e.g. Hybrid or Gold)
`oa_locations`|list-column of all the OA locations. 
`oa_locations_embargoed` | list-column of locations expected to be available in the future based on information like license metadata and journals' delayed OA policies
`data_standard`|Indicates the data collection approaches used for this resource. `1` mostly uses Crossref for hybrid detection. `2` uses more comprehensive hybrid detection methods. 
`is_oa`|Is there an OA copy (logical)? 
`is_paratext`| Is the item an ancillary part of a journal, like a table of contents? See here for more information <https://support.unpaywall.org/support/solutions/articles/44001894783>. 
`genre`|Publication type
`oa_status`|Classifies OA resources by location and license terms as one of: gold, hybrid, bronze, green or closed. See here for more information <https://support.unpaywall.org/support/solutions/articles/44001777288-what-do-the-types-of-oa-status-green-gold-hybrid-and-bronze-mean->.
`has_repository_copy`|Is a full-text available in a repository?
`journal_is_oa`|Is the article published in a fully OA journal? Uses the Directory of Open Access Journals (DOAJ) as source. 
`journal_is_in_doaj`|Is the journal listed in the Directory of Open Access Journals (DOAJ).
`journal_issns`|ISSNs, i.e. unique code to identify journals.
`journal_issn_l`|Linking ISSN.
`journal_name`|Journal title
`publisher`|Publisher
`published_date`|Date published
`year`|Year published. 
`title`|Publication title. 
`updated_resource`|Time when the data for this resource was last updated. 
`authors`|Lists authors (if available)

The columns  `best_oa_location` and  `oa_locations` are list-columns that contain useful metadata about the OA sources found by Unpaywall These are

**Column**|**Description**
|:------------|:----------------------------------------------
`endpoint_id`|Unique repository identifier
`evidence`|How the OA location was found and is characterized by Unpaywall?
`host_type`|OA full-text provided by `publisher` or `repository`. 
`is_best`|Is this location the \code{best_oa_location} for its resource?
`license`|The license under which this copy is published
`oa_date`|When this document first became available at this location
`pmh_id`|OAI-PMH endpoint where we found this location
`repository_institution`|Hosting institution of the repository.
`updated`|Time when the data for this location was last updated
`url`|The URL where you can find this OA copy.
`url_for_landing_page`| The URL for a landing page describing this OA copy.
`url_for_pdf`|The URL with a PDF version of this OA copy.
`versions`|The content version accessible at this location following the DRIVER 2.0 Guidelines  (<https://wiki.surfnet.nl/display/DRIVERguidelines/DRIVER-VERSION+Mappings>)

The Unpaywall schema is also described here: <https://unpaywall.org/data-format>.

The columns  `best_oa_location`. `oa_locations` and `oa_locations_embargoed` are list-columns that contain useful metadata about the OA sources found by Unpaywall. 

If `.flatten = TRUE` the list-column `oa_locations` will be restructured in a long format where each OA fulltext is represented by one row, which allows to take into account all OA locations found by Unpaywall in a data analysis.

```{r}
library(dplyr)
roadoi::oadoi_fetch(dois = c("10.1186/s12864-016-2566-9",
                             "10.1103/physreve.88.012814",
                             "10.1093/reseval/rvaa038",
                             "10.1101/2020.05.22.111294",
                             "10.1093/bioinformatics/btw541"), 
                    email = "najko.jahn@gmail.com", .flatten = TRUE) %>%
  dplyr::count(is_oa, evidence, is_best) 
```


#### Any API restrictions?

There are no API restrictions. However, Unpaywall requires an email address when using its API. If you are too tired to type in your email address every time, you can store the email  in the `.Renviron` file with the option `roadoi_email` 

```
roadoi_email = "najko.jahn@gmail.com"
```

You can open your `.Renviron` file calling 

```r
file.edit("~/.Renviron")`
```

Save the file and restart your R session. To stop sharing the email when using roadoi, delete it from your `.Renviron` file.

#### Keeping track of crawling

To follow your API call, and to estimate the time until completion, use the `.progress` parameter inherited from `plyr` to display a progress bar.

```{r}
roadoi::oadoi_fetch(dois = c("10.1186/s12864-016-2566-9",
                             "10.1103/physreve.88.012814"), 
                    email = "najko.jahn@gmail.com", 
                    .progress = "text")
```


### References 

Piwowar, H., Priem, J., Larivière, V., Alperin, J. P., Matthias, L., Norlander, B., … Haustein, S. (2018). The state of OA: a large-scale analysis of the prevalence and impact of Open Access articles. PeerJ, 6, e4375. <https://doi.org/10.7717/peerj.4375>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oadoiAddins.R
\name{roadoi_addin}
\alias{roadoi_addin}
\title{Find OA copies with RStudio addin}
\usage{
roadoi_addin()
}
\description{
An \href{https://rstudio.github.io/rstudioaddins/}{RStudio addin} to call
\code{oadoi_fetch()}. Shows up as "Find free full-texts" in the RStudio addin
menu.
}
\details{
The addin works as follows:

\enumerate{

\item Copy up to ten line-separated DOIs into the text area
\item Press the button "Run!"
\item Click on the links in the table to download full-text
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roadoi-package.r
\docType{package}
\name{roadoi-package}
\alias{roadoi-package}
\alias{roadoi}
\title{R Client for the Unpaywall-API}
\description{
What is this client for?:
roadoi interacts with the Unpaywall data service, which
links DOIs representing scholarly works with open access versions.
}
\details{
Use the \code{oadoi_fetch()} function in this package to get
open access status information and full-text links from Unpaywall.

You are welcome to contribute to this package. Use
GitHub's issue tracker for bug reporting and feature requests.

More details about the Unpaywall API:
\url{https://unpaywall.org/products/api}.
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oadoi_fetch.r
\name{oadoi_fetch_}
\alias{oadoi_fetch_}
\title{Get open access status information.}
\usage{
oadoi_fetch_(doi = NULL, email = NULL)
}
\arguments{
\item{doi}{character vector,a DOI}

\item{email}{character vector, mandatory!
Unpaywall requires your email address,
  so that they can track usage and notify you when something breaks.
  Set email address in your `.Renviron` file with
  the option
  `roadoi_email` \code{options(roadoi_email = "najko.jahn@gmail.com")}.
  You can open your `.Renviron` file calling `file.edit("~/.Renviron")`.
  Save the file and restart your R session. To stop sharing your email
  when using rcrossref, delete it from your `.Renviron` file.}
}
\value{
A tibble
}
\description{
In general, use \code{\link{oadoi_fetch}} instead. It calls this
method, returning open access status information from all your requests.
}
\examples{
\dontrun{
oadoi_fetch_(doi = c("10.1016/j.jbiotec.2010.07.030"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oadoi_fetch.r
\name{oadoi_fetch}
\alias{oadoi_fetch}
\title{Fetch open access status information and full-text links using Unpaywall}
\usage{
oadoi_fetch(
  dois = NULL,
  email = Sys.getenv("roadoi_email"),
  .progress = "none",
  .flatten = FALSE
)
}
\arguments{
\item{dois}{character vector, search by a single DOI or many DOIs.
A rate limit of 100k requests per day is suggested. If you need to access
more data, request the data dump
\url{https://unpaywall.org/dataset} instead.}

\item{email}{character vector, mandatory!
Unpaywall requires your email address,
so that they can track usage and notify you when something breaks.
Set email address in your `.Renviron` file with
the option
`roadoi_email` \code{options(roadoi_email = "najko.jahn@gmail.com")}.
You can open your `.Renviron` file calling `file.edit("~/.Renviron")`.
Save the file and restart your R session. To stop sharing your email
when using rcrossref, delete it from your `.Renviron` file.}

\item{.progress}{Shows the \code{plyr}-style progress bar.
Options are "none", "text", "tk", "win", and "time".
See \code{\link[plyr]{create_progress_bar}} for details
of each. By default, no progress bar is displayed.}

\item{.flatten}{Simplify open access evidence output. If `TRUE` it
transforms the nested column oa_locations so that each open access
evidence variable has its own column and each row represents a
single full-text.
Following these basic principles of "Tidy Data" makes data analysis
and export as a spreadsheet more straightforward.}
}
\value{
The result is a tibble with each row representing a publication.
  Here are the returned columns and descriptions according to the API docu:


\tabular{ll}{
 \code{doi}              \tab DOI (always in lowercase). \cr
 \code{best_oa_location} \tab list-column describing the best OA location.
 Algorithm prioritizes publisher hosted content (eg Hybrid or Gold),
 then prioritizes versions closer to the  version of record (PublishedVersion
 over AcceptedVersion), then more authoritative  repositories (PubMed Central
 over CiteSeerX). \cr
 \code{oa_locations}     \tab list-column of all the OA locations. \cr
 \code{oa_locations_embargoed}  \tab list-column of
 locations expected to be available in the future based on
  information like license metadata and journals'
  delayed OA policies \cr
 \code{data_standard}    \tab Indicates the data collection approaches used
 for this resource. \code{1} mostly uses Crossref for hybrid detection.
 \code{2} uses a more comprehensive hybrid detection methods. \cr
 \code{is_oa}            \tab Is there an OA copy (logical)? \cr
 \code{is_paratext}      \tab Is the item an ancillary part of a journal,
 like a table of contents? See here for more information
 \url{https://support.unpaywall.org/support/solutions/articles/44001894783}. \cr
 \code{genre}            \tab Publication type \cr
 \code{oa_status}        \tab Classifies OA resources by location
 and license terms as one of: gold, hybrid, bronze, green or closed.
 See here for more information
 \url{https://support.unpaywall.org/support/solutions/articles/44001777288-what-do-the-types-of-oa-status-green-gold-hybrid-and-bronze-mean-}. \cr
 \code{has_repository_copy} \tab Is a full-text
 available in a repository? \cr
 \code{journal_is_oa}    \tab Is the article published in a fully
 OA journal? \cr
 \code{journal_is_in_doaj} \tab Is the journal listed in
  the Directory of Open Access Journals (DOAJ). \cr
 \code{journal_issns}    \tab ISSNs, i.e. unique numbers to identify
 journals. \cr
 \code{journal_issn_l}    \tab Linking ISSN. \cr
 \code{journal_name}     \tab Journal title, not normalized. \cr
 \code{publisher}        \tab Publisher, not normalized. \cr
 \code{published_date}   \tab Date published \cr
 \code{year}             \tab Year published. \cr
 \code{title}            \tab Publication title. \cr
 \code{updated_resource} \tab Time when the data for this resource was last updated. \cr
 \code{authors}          \tab Lists author information (\code{family} name, \code{given}
 name and author role \code{sequence}), if available. \cr
}

The columns  \code{best_oa_location}. \code{oa_locations} and 
\code{oa_locations_embargoed} are list-columns that contain 
useful metadata about the OA sources found by Unpaywall. 

If \code{.flatten = TRUE} the list-column \code{oa_locations} will be
restructured in a long format where each OA fulltext is represented by
one row.

These are:

\tabular{ll}{
 \code{endpoint_id}     \tab Unique repository identifier. \cr
 \code{evidence}        \tab How the OA location was found and is characterized
  by Unpaywall? \cr
 \code{host_type}       \tab OA full-text provided by \code{publisher} or
  \code{repository}. \cr
 \code{is_best}         \tab Is this location the \code{best_oa_location} for its resource? \cr
 \code{license}         \tab The license under which this copy is published,
  e.g. Creative Commons license. \cr
 \code{pmh_id}          \tab OAI-PMH endpoint where we found this location. \cr
 \code{repository institution} \tab Hosting institution of the repository. \cr
 \code{updated}         \tab Time when the data for this location was last updated. \cr
 \code{url}             \tab The \code{url_for_pdf} if there is one; otherwise landing page URL. \cr
 \code{url_for_landing_page} \tab The URL for a landing page describing this OA copy. \cr
 \code{url_for_pdf}     \tab The URL with a PDF version of this OA copy. \cr
 \code{version}        \tab The content version accessible at this location
  following the DRIVER 2.0 Guidelines
 (\url{https://wiki.surfnet.nl/display/DRIVERguidelines/DRIVER-VERSION+Mappings}
 \cr
}

Note that Unpaywall schema is only informally described.
Check also \url{https://unpaywall.org/data-format}.
}
\description{
This is the main function to retrieve comprehensive open access status
information from Unpaywall data service. Please play nice with the API.
For each user, 100k calls per day are suggested. If you need to access
more data, there is also a data dump available.
For more info see \url{https://unpaywall.org/products/snapshot}.
}
\examples{
\dontrun{
oadoi_fetch("10.1038/nature12373", email = "name@example.com")
oadoi_fetch(dois = c("10.1016/j.jbiotec.2010.07.030",
"10.1186/1471-2164-11-245"), email = "name@example.com")
# flatten OA evidence
roadoi::oadoi_fetch(dois = c("10.1186/s12864-016-2566-9",
                            "10.1103/physreve.88.012814",
                            "10.1093/reseval/rvaa038"),
                   email = "najko.jahn@gmail.com", .flatten = TRUE)

}

}
