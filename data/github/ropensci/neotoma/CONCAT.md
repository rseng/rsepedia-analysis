neotoma
========

[![Build Status](https://api.travis-ci.org/ropensci/neotoma.png)](https://travis-ci.org/ropensci/neotoma)
[![Build status](https://ci.appveyor.com/api/projects/status/t2xyqbs0d8h998cb/branch/master)](https://ci.appveyor.com/project/sckott/neotoma/branch/master)
[![codecov.io](https://codecov.io/github/ropensci/neotoma/coverage.svg?branch=master)](https://codecov.io/github/ropensci/neotoma?branch=master)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/neotoma)](https://github.com/metacran/cranlogs.app)
[![cran version](http://www.r-pkg.org/badges/version/neotoma)](https://cran.r-project.org/package=neotoma)
[![NSF-1550707](https://img.shields.io/badge/NSF-1550707-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1550707)
[![NSF-1948926](https://img.shields.io/badge/NSF-1948926-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1948926)

**NOTE** *The Neotoma R package accesses the Windows-based API for the database.  This API is being deprecated, but is still available.  Due to changes in the database servers, the API Help is no longer functioning. Links to API help have been left, but do not resolve.*

The `neotoma` package is a programmatic R interface to the [Neotoma Paleoecological Database](http://www.neotomadb.org/). The package is intended to both allow users to search for sites and to download data for use in analyical workflows of paleoecological research.

`neotoma` is part of the [rOpenSci](http://ropensci.org) project and is also hosted on [Figshare](http://dx.doi.org/10.6084/m9.figshare.677131).  The `neotoma` package has been available on [CRAN](https://cran.r-project.org/package=neotoma) since May 3, 2015.

For more information on the package please refer to:

Goring, S., Dawson, A., Simpson, G. L., Ram, K., Graham, R. W., Grimm, E. C., & Williams, J. W.. (2015). neotoma: A Programmatic Interface to the Neotoma Paleoecological Database. *Open Quaternary*, 1(1), Art. 2. DOI: [10.5334/oq.ab](http://doi.org/10.5334/oq.ab)

For ongoing news, issues or information please join the [Neotoma Slack server](https://bit.ly/2FrZyYD) (and add the #r channel).  [Slack](http://www.slack.com/) is a real time chat application and collaboration hub with mobile and desktop applications.

### Development

*We welcome contributions from any individual, whether code, documentation, or issue tracking.  All participants are expected to follow the [code of conduct](https://github.com/ROpensci/neotoma/blob/master/code_of_conduct.md) for this project.*

+ [Simon Goring](http://goring.org) - University of Wisconsin-Madison, Department of Geography

### Contributors
+ [Gavin Simpson](http://www.fromthebottomoftheheap.net/) - University of Regina, Department of Biology
+ [Jeremiah Marsicek](http://geoweb.uwyo.edu/ggstudent/jmarsice/Site/Home.html) - University of Wyoming, Department of Geology and Geophysics
+ [Karthik Ram](http://nature.berkeley.edu/~kram/) - University of California - Berkely, Berkeley Institue for Data Science.
+ [Luke Sosalla](https://github.com/sosalla) - University of Wisconsin, Department of Geography

Package functions resolve various Neotoma APIs and re-form the data returned by the Neotoma database into R data objects.  The format of the Neotoma data, and the actual API functions can be accessed on the Neotoma API [website](http://api.neotomadb.org/doc/resources/home).

If you have used the package please consider providing us feedback through a [short survey](https://docs.google.com/forms/d/e/1FAIpQLSdRNat6L9grRF0xU5gibkr26xq9jD9wyHgw_AWxhrgn0lWv7w/viewform).

### Install `neotoma`

+ CRAN:
```r
install.packages('neotoma')
```

+ Development version from GitHub:
```r
install.packages("devtools")
library(devtools)
install_github("ropensci/neotoma")
library(neotoma)
```

### Currently implemented in `neotoma`

More functions are available through the package help.  These represent the core functions:

+ `get_site` - obtain information on sites in the Neotoma dataset (which may contain multiple datasets). [API](http://api.neotomadb.org/doc/resources/sites)
+ `get_dataset` - obtain dataset metadata from Neotoma. [API](http://api.neotomadb.org/doc/resources/datasets)
+ `get_download` - obtain full datasets (pollen or mammal) from Neotoma. [API](http://api.neotomadb.org/doc/resources/downloads)
+ `compile_list` - using established pollen-related taxonomies from the literature, take the published taxon list and standardize it to allow cross site analysis.
+ `get_contact` - find contact information for data contributors to Neotoma. [API](http://api.neotomadb.org/doc/resources/contacts)
+ `get_publication` - obtain publication information from Neotoma. [API](http://api.neotomadb.org/doc/resources/publications)
+ `get_table` - return matrices corresponding to one of the Neotoma database tables. [tables](http://api.neotomadb.org/doc/resources/dbtables)
+ `get_taxa` - Get taxon information from Neotoma. [API](http://api.neotomadb.org/doc/resources/taxa)
+ `get_chroncontrol` - Get chronological information used to build the age-depth model for the record. [API](http://api.neotomadb.org/doc/resources/chroncontrol)

### Recent Changes
+ 1.7.6: Updated the API endpoints to correctly point to the new windows API endpoint, in preparation for migration; Fixed the help for `get_closest()`; Changed character encoding for two data tables, `pollen.equiv` and `taxon.list`, made them available in the main package using `data()`.
+ 1.7.4: Bug fix: `get_dataset(gpid=123)` was returning an error, fix corrects the error to allow unassigned `x` variables.  Updated the allowable dataset types for searching to reflect the larger set of dataset types within Neotoma.
+ 1.7.3: Added numeric/integer methods to the `get_site()` and `get_dataset()` functions so that a vector of dataset or siteids can be passed to improve more general workflow methods.
+ 1.7.2: Bugfixes, added the `taxa()` function to easily extract taxa from one or multiple download objects.
+ 1.7.1: Bugfix for `compile_download()`, single sample downloads were failing to compile properly, added the `taxa()` function to extract taxa lists from large download objects.
+ 1.7.0: Added `plot_leaflet()` to allow interactive exploration of downloaded Neotoma data.  Integrates with the Neotoma Explorer.  Minor bugfix for `get_download()` to allow records to be sent to Neotoma and to be filtered.
+ 1.6.2: Improved the basic `plot()` method based on tests against Tilia files in the Neotoma Holding Tank & built more robust interpolation in `read_bacon()` so that age models without interpolated dates can still be imported. `browse()` now opens multiple datastes in the Neotoma Explorer at once.
+ 1.6.1: New `Stratiplot()` method, using the `analogue` package to plot dataset diagrams from `download` and `download_list` objects, bug fixes for `write_agefile()` and a new function, `read_bacon()`, to read in and integrate Bacon chronologies into `download` objects.
+ 1.6.0: Support for vector inputs in the `gpid` selection. Added a `get_closest()` function to find the closest sample site. Mostly clean-up of reported bugs by users. Revised examples for faster check speed.
+ 1.5.1: Minor fix to the `get_dataset()` for site level data to account for some datasets with empty submission data.  Some style changes to code (non-functional changes)
+ 1.5.0: More extensive testing to support multiple dataset types.  Water chemistry datasets still unsupported. Function `read.tilia()` added to read Tilia (http://tiliait.com) style XML files. Moved to using `xml2`, `httr` and `jsonlite` to support parsing.
+ 1.4.1: Small changes to `get_geochron()` to address bug reports and improve object `print` methods.
+ 1.4.0: Added `plot()` method for datasets, sites & downloads.  Fixed a bug with records missing chronologies.

### A few examples:

#### Find the distribution of sites with Mammoth fossils in Neotoma

```r
#  Example requires the mapdata package:
library('mapdata')

#  You may use either '%' or '*' as wildcards for search terms:
test <- get_dataset(taxonname='Mammuthus*')

The API call was successful, you have returned  3273 records.

site.locs <- get_site(test)

# A crude way of making the oceans blue.
plot(1, type = 'n',
     xlim=range(site.locs$long)+c(-10, 10),
     ylim=range(site.locs$lat)+c(-10, 10),
     xlab='Longitude', ylab = 'Latitude')
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightblue")
map('world',
    interior=TRUE,
    fill=TRUE,
    col='gray',
    xlim=range(site.locs$long)+c(-10, 10),
    ylim=range(site.locs$lat)+c(-10, 10),
    add=TRUE)

points(site.locs$long, site.locs$lat, pch=19, cex=0.5, col='red')

```
![thing](inst/img/mammothsites.png)

#### Plot the proportion of publications per year for datasets in Neotoma

```R
# Requires ggplot2
library('ggplot2')
library('plyr')
pubs <- get_publication()

pub.years <- ldply(pubs, "[[", "meta")

ggplot(data=pub.years, aes(x = year)) +
     stat_bin(aes(y=..density..*100, position='dodge'), binwidth=1) +
     theme_bw() +
     ylab('Percent of Publications') +
     xlab('Year of Publication') +
     scale_y_continuous(expand = c(0, 0.1)) +
     scale_x_continuous(breaks = seq(min(pub.years$year, na.rm=TRUE), 2014, by=20))

```

![thing](inst/img/histogramplot.png)

#### Cumulative plot of record uploads to Neotoma since 1998.

Found at [this gist](https://gist.github.com/SimonGoring/718a654f304f2d16ce4b)

![cumulativerecords](https://cloud.githubusercontent.com/assets/1619126/9884174/0026b83a-5ba4-11e5-9f1a-4a6874ceceb6.png)

#### Obtain records & Rebuild Chronologies with Bacon

Found at [this gist](https://gist.github.com/SimonGoring/877dd71cc3ad6bf8531e).  Prepared in part for a Bacon (Blaauw & Christen, 2011) workshop at the 2015 International Limnogeology Conference in Reno-Tahoe, Nevada led by Amy Myrbo (University of Minnesota).

#### Simple paleo-data visualization

Simple paleo-data visualization in R, linking the `rioja`, `neotoma` and `dplyr` packages.  Found at [this gist](https://gist.github.com/SimonGoring/dbb4c8e0087882dc143baa89fa041d2b).

![gif](inst/img/inkspot_neotoma.gif)

#### Find all site elevations in California:

Found at [Simon Goring's gist.](https://gist.github.com/SimonGoring/6a2ba1d55a3a7f78723b37e896b55b70).

#### Match all Neotoma taxa to external databases using `taxize`:

Found at [Simon Goring's gist.](https://gist.github.com/SimonGoring/24fb1228204f768f3f0020f37060db18).

#### Other Resources Using `neotoma`

*  [A simple `neotoma` workshop](http://www.goring.org/resources/Neotoma_Lesson.html)
*  [Data exploration and chronologies with `neotoma`](http://www.goring.org/resources/neotoma_lab_meeting.html)

### `neotoma` Workshops

We have provided a set of educational tools through the [NeotomaDB](http://github.com/neotomaDB) GitHub repository in the [Workshops](http://github.com/NeotomaDB/Workshops) repository.  These are free to share, and can be modified as needed.
## Test Environments
+   [windows passes](https://win-builder.r-project.org/2ssYmFp88U6R/).
+   Linux - [Travis](https://travis-ci.org/github/ropensci/neotoma/builds/728784231) tests are passing.
+   [Appveyor](https://ci.appveyor.com/project/sckott/neotoma/builds/35303211) tests are passing.

## R CMD check results:
+   R CMD check succeeded (0 errors, warnings, 1 note)
+   The single note is for timestamps and error described here: https://stackoverflow.com/questions/63613301/r-cmd-check-note-unable-to-verify-current-time
+   Windows check returns errors for missing API documentation.  Due to a change in the backend services the API documentation is currently offline, but should be restored shortly.

## Downstream Dependencies:
NA
# Contributor Code of Conduct

As contributors and maintainers of this project, and in the interest of fostering an open and welcoming community, we pledge to respect all people who contribute through reporting issues, posting feature requests, updating documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for everyone, regardless of level of experience, gender, gender identity and expression, sexual orientation, disability, personal appearance, body size, race, ethnicity, age, religion, or nationality.

Examples of unacceptable behavior by participants include:
 
  * The use of sexualized language or imagery
  * Personal attacks
  * Trolling or insulting/derogatory comments
  * Public or private harassment
  * Publishing other's private information, such as physical or electronic addresses, without explicit permission
  * Other unethical or unprofessional conduct

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

By adopting this Code of Conduct, project maintainers commit themselves to fairly and consistently applying these principles to every aspect of managing this project. Project maintainers who do not follow or enforce the Code of Conduct may be permanently removed from the project team.

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community.
      
Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting a project maintainer at simon.j.goring@gmail.com. All complaints will be reviewed and investigated and will result in a response that is deemed necessary and appropriate to the circumstances. Maintainers are obligated to maintain confidentiality with regard to the reporter of an incident.

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.3.0, available at [http://contributor-covenant.org/version/1/3/0/][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/3/0/
neotoma
========

[![Build Status](https://api.travis-ci.org/ropensci/neotoma.png)](https://travis-ci.org/ropensci/neotoma)
[![Build status](https://ci.appveyor.com/api/projects/status/t2xyqbs0d8h998cb/branch/master)](https://ci.appveyor.com/project/sckott/neotoma/branch/master)
[![codecov.io](https://codecov.io/github/ropensci/neotoma/coverage.svg?branch=master)](https://codecov.io/github/ropensci/neotoma?branch=master)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/neotoma)](https://github.com/metacran/cranlogs.app)
[![cran version](http://www.r-pkg.org/badges/version/neotoma)](https://cran.r-project.org/package=neotoma)

The `neotoma` package is a programmatic R interface to the [Neotoma Paleoecological Database](http://www.neotomadb.org/). The package is intended to both allow users to search for sites and to download data for use in analyical workflows of paleoecological research.

`neotoma` is part of the [rOpenSci](http://ropensci.org) project and is also hosted on [Figshare](http://dx.doi.org/10.6084/m9.figshare.677131).  The `neotoma` package is also available on [CRAN](https://cran.r-project.org/package=neotoma) as of May 3, 2015.

For more information on the package please refer to: 

Goring, S., Dawson, A., Simpson, G. L., Ram, K., Graham, R. W., Grimm, E. C., & Williams, J. W.. (2015). neotoma: A Programmatic Interface to the Neotoma Paleoecological Database. *Open Quaternary*, 1(1), Art. 2. DOI: [10.5334/oq.ab](http://doi.org/10.5334/oq.ab)

### Development

*We welcome contributions from any individual, whether code, documentation, or issue tracking.  All participants are expected to follow the [code of conduct](https://github.com/ROpensci/neotoma/blob/master/code_of_conduct.md) for this project.*

+ [Simon Goring](http://downwithtime.wordpress.com) - University of Wisconsin-Madison, Department of Geography

### Contributors
+ [Gavin Simpson](http://www.fromthebottomoftheheap.net/) - University of Regina, Department of Biology
+ [Jeremiah Marsicek](http://geoweb.uwyo.edu/ggstudent/jmarsice/Site/Home.html) - University of Wyoming, Department of Geology and Geophysics
+ [Karthik Ram](http://nature.berkeley.edu/~kram/) - University of California - Berkely, Berkeley Institue for Data Science.
+ [Luke Sosalla](https://github.com/sosalla) - University of Wisconsin, Department of Geography

Package functions resolve various Neotoma APIs and re-form the data returned by the Neotoma database into R data objects.  The format of the Neotoma data, and the actual API functions can be accessed on the Neotoma API [website](http://api.neotomadb.org/doc/resources/home).

If you have used the package please consider providing us feedback through a [short survey](https://docs.google.com/forms/d/e/1FAIpQLSdRNat6L9grRF0xU5gibkr26xq9jD9wyHgw_AWxhrgn0lWv7w/viewform).

### Install `neotoma`

+ CRAN:
```r
install.packages('neotoma')
```

+ Development version from GitHub:
```r
install.packages("devtools")
library(devtools)
install_github("ropensci/neotoma")
library(neotoma)
```

### Currently implemented in `neotoma`

More functions are available through the package help.  These represent the core functions:

+ `get_site` - obtain information on sites in the Neotoma dataset (which may contain multiple datasets). [API](http://api.neotomadb.org/doc/resources/sites)
+ `get_dataset` - obtain dataset metadata from Neotoma. [API](http://api.neotomadb.org/doc/resources/datasets)
+ `get_download` - obtain full datasets (pollen or mammal) from Neotoma. [API](http://api.neotomadb.org/doc/resources/downloads)
+ `compile_list` - using established pollen-related taxonomies from the literature, take the published taxon list and standardize it to allow cross site analysis.
+ `get_contact` - find contact information for data contributors to Neotoma. [API](http://api.neotomadb.org/doc/resources/contacts)
+ `get_publication` - obtain publication information from Neotoma. [API](http://api.neotomadb.org/doc/resources/publications)
+ `get_table` - return matrices corresponding to one of the Neotoma database tables. [tables](http://api.neotomadb.org/doc/resources/dbtables)
+ `get_taxa` - Get taxon information from Neotoma. [API](http://api.neotomadb.org/doc/resources/taxa)
+ `get_chroncontrol` - Get chronological information used to build the age-depth model for the record. [API](http://api.neotomadb.org/doc/resources/chroncontrol)

### Recent Changes
+ 1.7.4: Bug fix: `get_dataset(gpid=123)` was returning an error, fix corrects the error to allow unassigned `x` variables.  Updated the allowable dataset types for searching to reflect the larger set of dataset types within Neotoma.
+ 1.7.3: Added numeric/integer methods to the `get_site()` and `get_dataset()` functions so that a vector of dataset or siteids can be passed to improve more general workflow methods.
+ 1.7.2: Bugfixes, added the `taxa()` function to easily extract taxa from one or multiple download objects.
+ 1.7.1: Bugfix for `compile_download()`, single sample downloads were failing to compile properly, added the `taxa()` function to extract taxa lists from large download objects.
+ 1.7.0: Added `plot_leaflet()` to allow interactive exploration of downloaded Neotoma data.  Integrates with the Neotoma Explorer.  Minor bugfix for `get_download()` to allow records to be sent to Neotoma and to be filtered.
+ 1.6.2: Improved the basic `plot()` method based on tests against Tilia files in the Neotoma Holding Tank & built more robust interpolation in `read_bacon()` so that age models without interpolated dates can still be imported. `browse()` now opens multiple datastes in the Neotoma Explorer at once.
+ 1.6.1: New `Stratiplot()` method, using the `analogue` package to plot dataset diagrams from `download` and `download_list` objects, bug fixes for `write_agefile()` and a new function, `read_bacon()`, to read in and integrate Bacon chronologies into `download` objects.
+ 1.6.0: Support for vector inputs in the `gpid` selection. Added a `get_closest()` function to find the closest sample site. Mostly clean-up of reported bugs by users. Revised examples for faster check speed.
+ 1.5.1: Minor fix to the `get_dataset()` for site level data to account for some datasets with empty submission data.  Some style changes to code (non-functional changes)
+ 1.5.0: More extensive testing to support multiple dataset types.  Water chemistry datasets still unsupported. Function `read.tilia()` added to read Tilia (http://tiliait.com) style XML files. Moved to using `xml2`, `httr` and `jsonlite` to support parsing.
+ 1.4.1: Small changes to `get_geochron()` to address bug reports and improve object `print` methods.
+ 1.4.0: Added `plot()` method for datasets, sites & downloads.  Fixed a bug with records missing chronologies.

### A few examples:

#### Find the distribution of sites with Mammoth fossils in Neotoma

```r
#  Example requires the mapdata package:
library('mapdata')

#  You may use either '%' or '*' as wildcards for search terms:
test <- get_dataset(taxonname='Mammuthus*')

The API call was successful, you have returned  3273 records.

site.locs <- get_site(test)

# A crude way of making the oceans blue.
plot(1, type = 'n',
     xlim=range(site.locs$long)+c(-10, 10),
     ylim=range(site.locs$lat)+c(-10, 10),
     xlab='Longitude', ylab = 'Latitude')
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightblue")
map('world',
    interior=TRUE,
    fill=TRUE,
    col='gray',
    xlim=range(site.locs$long)+c(-10, 10),
    ylim=range(site.locs$lat)+c(-10, 10),
    add=TRUE)

points(site.locs$long, site.locs$lat, pch=19, cex=0.5, col='red')

```
![thing](inst/img/mammothsites.png)

#### Plot the proportion of publications per year for datasets in Neotoma

```R
# Requires ggplot2
library('ggplot2')
library('plyr')
pubs <- get_publication()

pub.years <- ldply(pubs, "[[", "meta")

ggplot(data=pub.years, aes(x = year)) +
     stat_bin(aes(y=..density..*100, position='dodge'), binwidth=1) +
     theme_bw() +
     ylab('Percent of Publications') +
     xlab('Year of Publication') +
     scale_y_continuous(expand = c(0, 0.1)) +
     scale_x_continuous(breaks = seq(min(pub.years$year, na.rm=TRUE), 2014, by=20))

```

![thing](inst/img/histogramplot.png)

#### Cumulative plot of record uploads to Neotoma since 1998.

Found at [this gist](https://gist.github.com/SimonGoring/718a654f304f2d16ce4b)

![cumulativerecords](https://cloud.githubusercontent.com/assets/1619126/9884174/0026b83a-5ba4-11e5-9f1a-4a6874ceceb6.png)

#### Obtain records & Rebuild Chronologies with Bacon

Found at [this gist](https://gist.github.com/SimonGoring/877dd71cc3ad6bf8531e).  Prepared in part for a Bacon (Blaauw & Christen, 2011) workshop at the 2015 International Limnogeology Conference in Reno-Tahoe, Nevada led by Amy Myrbo (University of Minnesota).

#### Simple paleo-data visualization

Simple paleo-data visualization in R, linking the `rioja`, `neotoma` and `dplyr` packages.  Found at [this gist](https://gist.github.com/SimonGoring/dbb4c8e0087882dc143baa89fa041d2b).

![gif](inst/img/inkspot_neotoma.gif)

#### Find all site elevations in California:

Found at [Simon Goring's gist.](https://gist.github.com/SimonGoring/6a2ba1d55a3a7f78723b37e896b55b70).

#### Match all Neotoma taxa to external databases using `taxize`:

Found at [Simon Goring's gist.](https://gist.github.com/SimonGoring/24fb1228204f768f3f0020f37060db18).

#### Other Resources Using `neotoma`

*  [A simple `neotoma` workshop](http://www.goring.org/resources/Neotoma_Lesson.html)
*  [Data exploration and chronologies with `neotoma`](http://www.goring.org/resources/neotoma_lab_meeting.html)

### `neotoma` Workshops

We have provided a set of educational tools through the [NeotomaDB](http://github.com/neotomaDB) GitHub repository in the [Workshops](http://github.com/NeotomaDB/Workshops) repository.  These are free to share, and can be modified as needed.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol.download_list}
\alias{get_chroncontrol.download_list}
\title{Function to return chronological control tables from a \code{download_list} object.}
\usage{
\method{get_chroncontrol}{download_list}(x, chronology = 1,
  verbose = TRUE, add = FALSE)
}
\arguments{
\item{x}{A \code{download_list} object.}

\item{chronology}{When \code{download} objects have more than associated chronology, which chronology do you want?  Default is \code{1}.}

\item{verbose}{logical; should messages on API call be printed?}

\item{add}{Should the \code{chroncontrol} be added to the download object (default \code{FALSE})}
}
\description{
Using a \code{download_list}, return the default chron-control table as a \code{data.frame}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/translate.table.R
\docType{data}
\name{translate.table}
\alias{translate.table}
\title{A table to convert the original taxa to standardized lists.}
\format{a \code{data.frame} object}
\source{
The Neotoma database.
}
\usage{
translate.table
}
\description{
A list of standardized (published) taxonomies from the literature to help standardize taxonomies for synthesis work.
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse}
\alias{browse}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
browse(x)
}
\arguments{
\item{x}{A numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object.}
}
\value{
Returns a NULL value, opens a browser.
}
\description{
Using a \code{download} or \code{dataset} object, open up a browser window in the users default browser. Passing a \code{download_list} or \code{dataset_list} will open Neotoma Explorer with the first object and return a warning.

Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
\examples{
\dontrun{
# Where are the XRF data?

xrf.data <- get_dataset(datasettype='X-ray fluorescence (XRF)')
browse(xrf.data)

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/sites
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_methods.R
\name{read.tilia}
\alias{read.tilia}
\title{Read proxy data from Tilia TLX files}
\usage{
read.tilia(file)
}
\arguments{
\item{file}{a string representing a Tilia TLX format file.}
}
\value{
Return a `download` object.
}
\description{
Read proxy data from a Tilia TLX format file.
}
\examples{
\dontrun{
  crystal <- read.tilia('crystal.tlx')
}

}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.numeric}
\alias{get_site.numeric}
\title{Return Site information from a vector of numeric elements.}
\usage{
\method{get_site}{numeric}(sitename, ...)
}
\arguments{
\item{sitename}{A numeric value or vector of numeric elements.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol.default}
\alias{get_chroncontrol.default}
\title{Function to return chronological control tables from a chronologic ID.}
\usage{
\method{get_chroncontrol}{default}(x, chronology = 1, verbose = TRUE,
  add = FALSE)
}
\arguments{
\item{x}{A single numeric chronology ID or a vector of numeric chronology IDs as returned by \code{get_datasets}.}

\item{chronology}{For \code{download} methods, which chronology controls should be used?}

\item{verbose}{logical; should messages on API call be printed?}

\item{add}{logical, should this chron control be added to the download object?}
}
\description{
Using the chronology ID, return the chron control table as a \code{data.frame}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param_check.R
\name{param_check}
\alias{param_check}
\title{Internal function to check passed parameters.}
\usage{
param_check(cl)
}
\arguments{
\item{cl}{Contact ID is a numerical value associated with the Neotoma
Contact table's numerical Contact ID.}
}
\value{
A list with two components:

 \item{flag}{Returns a 0 if everything's fine, a 1 if there's a problem.}
 \item{message}{A list of error messages.}
}
\description{
Functions \code{\link{get_site}}, \code{\link{get_dataset}} and others pass parameters to \code{param_check}, which tells them if there's a problem.
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{internal}
\keyword{misc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_contact.R
\name{get_contact}
\alias{get_contact}
\title{Get contact information.}
\usage{
get_contact(contactid, contactname, contactstatus, familyname)
}
\arguments{
\item{contactid}{Contact ID is a numerical value associated with the Neotoma
Contact table's numerical Contact ID.}

\item{contactname}{A character string indicating the data contributors' project,
organization or personal name.  May be a partial string and can include wildcards.}

\item{contactstatus}{The current status of the contact.  Possible values include:
active, deceased, defunct, extant, inactive, retired, unknown.}

\item{familyname}{A character string.  Full or partial string indicating the
contact's last name.}
}
\value{
The function takes parameters defined by the user and returns a list
   of contact information supplied by the Neotoma Paleoecological Database.
   The user may define all or none of the possible fields.  The function contains
   data checks for each defined parameter.

   The function returns either a single item of class \code{"try-error"} describing
   the reason for failure (either mis-defined parameters or an error from the Neotoma API),
   or a table of contacts, with rows corresponding to the number of individual
   contacts returned by the Neotoma API.  Each row entry includes the following parameters:

 \item{ \code{contact.name} }{Full name of the person, last name first (e.g. \code{"Simpson, George Gaylord"}) or name of organization or project (e.g. \code{"Great Plains Flora Association"}).}
 \item{ \code{contact.status} }{Current status of the person, organization, or project. Field links to the ContactStatuses lookup table.}
 \item{ \code{family.name} }{Family or surname name of a person.}
 \item{ \code{leading.initials} }{Leading initials for given or forenames without spaces (e.g. \code{"G.G."}).}
 \item{ \code{given.names} }{Given or forenames of a person (e.g. \code{"George Gaylord"}). Initials with spaces are used if full given names are not known (e.g. \code{"G. G")}.}
 \item{ \code{suffix} }{Suffix of a person's name (e.g. \code{"Jr."}, \code{"III"}).}
 \item{ \code{title} }{A person's title (e.g. \code{"Dr."}, \code{"Prof."}, \code{"Prof. Dr"}).}
 \item{ \code{phone} }{Telephone number.}
 \item{ \code{fax} }{Fax number.}
 \item{ \code{email} }{Email address.}
 \item{ \code{url} }{Universal Resource Locator, an Internet World Wide Web address.}
 \item{ \code{address} }{Full mailing address.}
 \item{ \code{notes} }{Free form notes or comments about the person, organization, or project.}
 \item{ \code{contact.id} }{Unique database record identifier for the contact.}
 \item{ \code{alias.id} }{The ContactID of a person's current name. If the AliasID is different from the ContactID, the ContactID refers to the person's former name.}
}
\description{
A function to obtain contact information for data contributors from the Neotoma Paleoecological Database.
}
\examples{
\dontrun{
#  To find all data contributors who are active:
active.cont <- get_contact(contactstatus = 'active')

# To find all data contributors who have the last name "Smith"
smith.cont <- get_contact(familyname = 'Smith')
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse.dataset}
\alias{browse.dataset}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
\method{browse}{dataset}(x)
}
\arguments{
\item{x}{A \code{dataset} object.}
}
\description{
Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.site}
\alias{get_dataset.site}
\title{Obtain dataset information from an existing \code{site} object.}
\usage{
\method{get_dataset}{site}(x, ...)
}
\arguments{
\item{x}{An object of class \code{site}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_leaflet.R
\name{plot_leaflet}
\alias{plot_leaflet}
\title{Leaflet plots for neotoma data.}
\usage{
plot_leaflet(x, providerTiles = "Stamen.TerrainBackground", ...)
}
\arguments{
\item{x}{A neotoma data object}

\item{providerTiles}{Default "Stamen.TerrainBackground", a character string indicating the tile background to be used for plotting.}

\item{...}{Other terms to be passed to the function.}
}
\value{
A \code{leaflet} object
}
\description{
A plotting function to provide interactive data investigation using the leaflet tools.
  This package requires a connection to the internet for proper functioning.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_closest.R
\name{get_closest}
\alias{get_closest}
\title{Find the closest dataset records to a site, dataset or long/lat pair in Neotoma}
\usage{
get_closest(x, n, buffer, ...)
}
\arguments{
\item{x}{A vector long/lat pair, or a dataset, site or download.}

\item{n}{The number of records to return.}

\item{buffer}{The size of the buffer for dataset search (in kilometers)}

\item{...}{optional arguments to pass into \code{get_dataset}.}
}
\value{
This command returns a \code{dataset} or \code{dataset_list}, or NULL if no records exist within the bounding box.
}
\description{
Passing in a download object the function outputs a Bacon or Clam formatted file to a
user defined destination for age modelling with existing age-depth modeling software.
}
\examples{
\dontrun{
#  The point of pulling chronology tables is to re-build or examine the chronological 
#  information that was used to build the age-depth model for the core.
# Find the closest records to Madison, WI:
get_closest(x = c(-89.4012, 43.0731), n = 10, buffer = 5, datasettype = "pollen")
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}, Andria Dawson \email{andria.dawson@gmail.com}
}
\keyword{API}
\keyword{Neotoma}
\keyword{Palaeoecology}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_agefile.R
\name{write_agefile}
\alias{write_agefile}
\title{Write age control file to disk formatted for either Bacon or Clam}
\usage{
write_agefile(download, chronology = 1, path, corename,
  cal.prog = "Bacon")
}
\arguments{
\item{download}{A single site returned by \code{get_download}.}

\item{chronology}{Default is \code{1}, the default chronology for the core.  If a core has more than one chronology the user can define a different set of chronological controls.}

\item{path}{The location of the 'Cores' folder & working directory for Bacon.  Do not include "Cores" in the path name.}

\item{corename}{The intended handle for the core, to be used in writing to file.}

\item{cal.prog}{The method intended to build the age model, either \code{'Bacon'} or \code{'Clam'}.}
}
\value{
This command returns a file in location \code{path/Cores} containing all the relevant information required to build either the default or prior chronology for a core.
}
\description{
Passing in a download object the function outputs a Bacon or Clam formatted file to a
user defined destination for age modelling with existing age-depth modeling software.
}
\examples{
\dontrun{
# Find a particular record:

three_pines <- get_download(get_dataset(get_site("Three Pines Bog"), 
                                        datasettype = "pollen"))

# You will need to edit the `path` argument here to point to a directory that 
# contains a `Cores` directory.

write_agefile(download = three_pines[[1]], 
              path = "./inst", 
              corename = "THREEPINES", 
              cal.prog = "Bacon")
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{API}
\keyword{Neotoma}
\keyword{Palaeoecology}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse.dataset_list}
\alias{browse.dataset_list}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
\method{browse}{dataset_list}(x)
}
\arguments{
\item{x}{A \code{dataset_list} object.}
}
\description{
Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.default}
\alias{get_site.default}
\title{Return Site Information.}
\usage{
\method{get_site}{default}(sitename, ...)
}
\arguments{
\item{sitename}{A character string representing the full or partial site name.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.integer}
\alias{get_dataset.integer}
\title{Obtain dataset information from a vector of dataset IDs.}
\usage{
\method{get_dataset}{integer}(x, ...)
}
\arguments{
\item{x}{A single numeric dataset id, or a numeric vector.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_download.R
\name{get_download.default}
\alias{get_download.default}
\title{Function to return full download records using \code{numeric} dataset IDs.}
\usage{
\method{get_download}{default}(x, verbose = TRUE)
}
\arguments{
\item{x}{A single numeric dataset ID or a vector of numeric dataset IDs as returned by \code{get_datasets}.}

\item{verbose}{logical; should messages on API call be printed?}
}
\description{
Using the dataset ID, return all records associated with the data as a \code{download_list}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.numeric}
\alias{get_dataset.numeric}
\title{Obtain dataset information from a vector of dataset IDs.}
\usage{
\method{get_dataset}{numeric}(x, ...)
}
\arguments{
\item{x}{A single numeric dataset id, or a numeric vector.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol.download}
\alias{get_chroncontrol.download}
\title{Function to return chronological control tables from a \code{download} object.}
\usage{
\method{get_chroncontrol}{download}(x, chronology = 1, verbose = TRUE,
  add = FALSE)
}
\arguments{
\item{x}{A single \code{download} object.}

\item{chronology}{For \code{download} methods, which chronology controls should be used?}

\item{verbose}{logical; should messages on API call be printed?}

\item{add}{Should the \code{chroncontrol} be added to the download object (default \code{FALSE})}
}
\description{
Using a \code{download}, return the default chron-control table as a \code{data.frame}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_bacon.R
\name{read_bacon}
\alias{read_bacon}
\title{Function to read in defined Bacon outputs.}
\usage{
read_bacon(x, path = ".", add = FALSE, chron_name = "Bacon",
  as_default = TRUE, download = NULL, sections = NULL,
  age_field = "median", interp = TRUE)
}
\arguments{
\item{x}{A folder path that contains a Bacon \code{age} file.}

\item{path}{The location of the \code{Cores} folder.}

\item{add}{Should the results be added to an existing \code{download}? Defaults to \code{FALSE}.}

\item{chron_name}{The name for the chronology if the Bacon file is being added to a \code{download}.}

\item{as_default}{Should the chronology become the default?}

\item{download}{The target \code{download} if \code{add} is \code{TRUE}.}

\item{sections}{If there are multiple Bacon runs in a folder, identify the file by the number of sections in the run.}

\item{age_field}{Should the age be assigned to the \code{"median"} or the \code{"wmean"}?}

\item{interp}{If the depths don't match up, should we interpolate from the Bacon output? (default \code{TRUE})}
}
\description{
Reads in Bacon output and formats it for inclusion in a download object.
}
\details{
The function expects that you are in a working directory containing a "Cores" which would then contain output files from Bacon runs.  The output can either be added to an existing record (for example, replacing the default age model returned by Neotoma), or it can be loaded on its own.
If the depths for the loaded file do not match with the depths in the `download`'s `sample.meta` then the user can use the `interp` parameter to interpolate between depths.  This method uses linear interpolation.
}
\examples{
\dontrun{
# Download the record for Lake O' Pines:
lake_o_dl <- get_download(15925)

# This assumes that you have Bacon installed in a folder and have
# set it to your working directory.

write_agefile(lake_o_dl[[1]], path = ".", chronology = 1, 
              corename = "LAKEPINES", cal.prog = 'Bacon') 

source("Bacon.R") 

# These defaults just help the core run quickly, they're not 
# neccesarily good parameters.

Bacon("LAKEPINES", acc.mean = 10, 
      thick = 50, depths.file = TRUE, 
      suggest = FALSE, ask = FALSE)

lake_o_dl <- read_bacon("LAKEPINES", add = TRUE, 
                        download = download, sections = 17)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gp.table.R
\docType{data}
\name{gp.table}
\alias{gp.table}
\title{A list of all the geopolitical entities in the Neotoma database.}
\format{a \code{data.frame} object}
\source{
The Neotoma database.
}
\usage{
gp.table
}
\description{
A list of geopolitical entities with associated numeric ID values.
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse.download}
\alias{browse.download}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
\method{browse}{download}(x)
}
\arguments{
\item{x}{A \code{download} object.}
}
\description{
Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bind.R
\name{bind}
\alias{bind}
\title{Function to bind objects together into a longer object.}
\usage{
bind(x, ...)
}
\arguments{
\item{x}{An object returned by one of the \code{get_*} commands for download, site or dataset.}

\item{...}{other objects of the same class.}
}
\value{
This command returns a larger list.
}
\description{
From multiple \code{download*}s, \code{dataset*}s or \code{site}s, join them together into a single object.
}
\details{
To support further synthesis and analysis \code{compile_download} works to transform a list
returned by \code{\link{get_download}} into a large data frame with columns for site and sample attributes
and also with the associated assemblage data at each sample depth.  This function also does the same for
single sites.
}
\examples{
\dontrun{
#  Search for sites with "Thuja" pollen that are older than 8kyr BP and
#  that are on the west coast of North America:
t8kyr.poa <- get_dataset(taxonname="Thuja*", 
                         loc=c(-150, 20, -100, 60), ageyoung = 8000)
t8kyr.canis <- get_dataset(taxonname="Canis*", 
                           loc=c(-150, 20, -100, 60), ageyoung = 8000)

t8kyr.co_site <- bind(t8kyr.poa, t8kyr.canis)
plot(t8kyr.co_site)

####
# We want to look at four different dataset types across a forest-prairie 
# boundary:
dataset_types <- c("ostracode surface sample",
                   "water chemistry",
                   "diatom surface sample",
                   "pollen surface sample")

# Run the `get_dataset` function for each of the different dataset types 
dataset_lists <- lapply(dataset_types, 
                          function(x) { 
                            get_dataset(datasettype=x, 
                                        loc = c(-100,43,-92,48))
                                        })

# Using do.call here to make sure that I don't have to split the list out.
new_datasets <- do.call(bind, dataset_lists)

# And voila!
plot(new_datasets)

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pollen.equiv.R
\docType{data}
\name{pollen.equiv}
\alias{pollen.equiv}
\title{A table to convert the pollen taxa identified by investigators to standardized lists.}
\format{a \code{data.frame} object}
\usage{
translate.table
}
\description{
A list of standardized (published) taxonomies from the literature to help standardize taxonomies for synthesis work.
}
\details{
Taxon conversion table (readable).
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}; Jeremiah Marsicek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.dataset}
\alias{get_site.dataset}
\title{Return Site Information from a numeric list of site ids.}
\usage{
\method{get_site}{dataset}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{dataset}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compile_downloads.R
\name{compile_downloads}
\alias{compile_downloads}
\title{Function to convert multiple downloads into a single large table.}
\usage{
compile_downloads(downloads)
}
\arguments{
\item{downloads}{A download_list as returned by \code{\link{get_download}}, or multiple downloads joined in a list.}
}
\value{
This command returns a data frame.
}
\description{
From the assemblage data for multiple cores, return a single data.frame with columns for site
metadata and assemblage data.
}
\details{
To support further synthesis and analysis \code{compile_download} works to transform a list
returned by \code{\link{get_download}} into a large data frame with columns for site and sample attributes
and also with the associated assemblage data at each sample depth.  This function also does the same for
single sites.
}
\examples{
\dontrun{
#  Search for sites with "Thuja" pollen that are older than 8kyr BP and
#  that are on the west coast of North America:
t8kyr.datasets <- get_dataset(taxonname='Thuja*', 
                              loc=c(-150, 20, -100, 60), 
                              ageyoung = 8000)

#  Returns 3 records (as of 04/04/2013), get dataset for the first record, 
#  Gold Lake Bog.
thuja.sites <- get_download(t8kyr.datasets)

gold.p25 <- compile_taxa(thuja.sites, 'P25')

all.gold <- compile_downloads(gold.p25)

pollen.sums <- rowSums(all.gold[,11:ncol(all.gold)], na.rm=TRUE)

plot(x = all.gold$age, 
     y = all.gold$Cupressaceae.Taxaceae / pollen.sums, 
     col = all.gold$site.name,
     pch = 19)

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org

Gavin DG, Oswald WW, Wahl ER, Williams JW. 2003. A statistical approach to evaluating distance metrics and analog assignments for pollen records. Quaternary Research 60: 356-367.

Whitmore J, Gajewski K, Sawada M, Williams JW, Shuman B, Bartlein PJ, Minckley T, Viau AE, Webb III T, Shafer S, Anderson P, Brubaker L. 2005. Modern pollen data from North America and Greenland for multi-scale paleoenvironmental applications. Quaternary Science Reviews 24: 1828-1848.

Williams J, Shuman B. 2008. Obtaining accurate and precise environmental reconstructions from the modern analog technique and North American surface pollen dataset. Quaternary Science Reviews. 27:669-687.

API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clamodel.R, R/download.R
\name{download}
\alias{download}
\title{A class for download objects.}
\description{
A \code{download} is an object with the full record for a single dataset.

A \code{download} is an object with the full record for a single dataset.
}
\details{
TO DO

TO DO
}
\author{
Simon Goring

Simon Goring
}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication.download}
\alias{get_publication.download}
\title{A function to get publications for downloads in the Neotoma Database using the API.}
\usage{
\method{get_publication}{download}(x, ...)
}
\arguments{
\item{x}{an object of class \code{download}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
The function takes a \code{download} and returns a table with publication information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/counts.R
\name{counts}
\alias{counts}
\alias{counts.download}
\alias{counts.download_list}
\title{Access proxy count data}
\usage{
counts(obj, ...)

\method{counts}{download}(obj, ...)

\method{counts}{download_list}(obj, ...)
}
\arguments{
\item{obj}{an R object from which counts are to be extracted.}

\item{...}{arguments passed to other methods.}
}
\value{
Either a data frame of counts or a list of such objects.
}
\description{
Extract pollen or other proxy counts from data objects and returns them in a useful format.
}
\details{
Methods are available for "download" and "download_list" objects.
}
\examples{
\dontrun{
marion <- get_site('Marion Lake\%')
louise <- get_site('Louise Pond\%')
western.sites <- rbind(marion, louise)
western.data  <- get_dataset(western.sites)

western.dl <- get_download(western.data)
western.cnt <- counts(western.dl)
sapply(western.cnt, dim)
marion.cnt<- counts(western.dl[[1]])
dim(marion.cnt)
}
}
\author{
Gavin Simpson
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Stratiplot.download.R
\name{Stratiplot.download_list}
\alias{Stratiplot.download_list}
\title{Palaeoecological stratigraphic diagrams}
\usage{
\method{Stratiplot}{download_list}(x, yaxis = "age", method = "none",
  group = NULL, ...)
}
\arguments{
\item{x}{A \code{download_list} object.}

\item{yaxis}{One of the columns in \code{sample.meta}, including \code{depth}, \code{age}, \code{age.younger}, or \code{age.older}, default \code{age}.}

\item{method}{An option for axis transformation using \code{tran} from the \code{analogue} package.  \code{"none"} by default.}

\item{group}{An ecological group from the taxon table.}

\item{...}{variables to be passed to \code{Stratiplot}.}
}
\value{
A \code{trellis} object.
}
\description{
Draws paleoecological diagrams from a \code{download_list} object.  Allows control of variable type (using the \code{tran} function from the \code{analogue} package), and taxonomic grouping.  
This function only works for \code{download_list} objects that contain a single object.
}
\details{
A wrapper for the \code{analogue} package's \code{Stratiplot} function.  Allowing the user to plot a stratigraphic diagram directly from a \code{download} object.
}
\examples{
\dontrun{
lake_o_dl <- get_download(15925)
# This works:
Stratiplot(lake_o_dl)

lakes_o_nw <- get_download(get_site(sitename = "Lake B\%"))
# This Fails:
# Stratiplot(lake_o_nw)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.dataset_list}
\alias{get_site.dataset_list}
\title{Return Site Information from a \code{dataset_list}}
\usage{
\method{get_site}{dataset_list}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{dataset_list}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_geochron.R
\name{get_geochron}
\alias{get_geochron}
\title{Function to return geochronological data from records.}
\usage{
get_geochron(x, verbose = TRUE)
}
\arguments{
\item{x}{A numeric dataset ID or a vector of numeric dataset IDs, or an object of class of class \code{site}, \code{dataset}, \code{dataset_list}, \code{download} or \code{download_list} for which geochrons are required.}

\item{verbose}{logical; should messages on API call be printed?}
}
\value{
This command returns either an object of class \code{"try-error"}' (see \code{\link{try}}) defined by the error returned
   from the Neotoma API call, or a \code{geochronologic} object, which is a list with two components, a \code{dataset} and a geochronology table, a \code{data.frame} with the following components:

 \item{ \code{sample.id} }{A unique identifier for the geochronological unit.}
 \item{ \code{age.type} }{String.  The age type, one of calendar years, radiocarbon years, etc.}
 \item{ \code{age} }{Dated age of the material.}
 \item{ \code{e.older} }{The older error limit of the age value.  Commonly 1 standard deviation.}
 \item{ \code{e.young} }{The younger error limit of the age value.}
 \item{ \code{delta13C} }{The measured or assumed delta13C value for radiocarbon dates, if provided.}
 \item{ \code{material.dated} }{A table describing the collection, including dataset information, PI data compatible with \code{\link{get_contact}} and site data compatable with \code{\link{get_site}}.}
 \item{ \code{geo.chron.type} }{Text string, type of geochronological analysis, i.e., Radiocarbon dating, luminesence.}
 \item{ \code{notes} }{Text string}
 \item{ \code{infinite} }{Boolean, does the dated material return an "infinite" date?}

 A full data object containing all the relevant geochronological data available for a dataset.
}
\description{
Using the dataset ID, return all geochronological data associated with the dataID.  At present,
   only returns the dataset in an unparsed format, not as a data table.   This function will only download one dataset at a time.
}
\examples{
\dontrun{
#  Search for the sites around Marion Lake, BC.  I want to find sites within 
#  about 1km.

marion <- get_site(sitename = "Marion Lake*")

marion_close <- get_closest(marion, n = 10, buffer = 1)

#  Returns 116 records (as of 13/07/2015).  These are the pollen records though, 
#  we want the sites:
geochron.records <- get_geochron(marion_close)

#  We want to extract all the radiocarbon ages from the records:

get_ages <- function(x){
  any.ages <- try(x[[2]]$age[x[[2]]$age.type == 'Radiocarbon years BP'])
  if(class(any.ages) == 'try-error') output <- NA
  if(!class(any.ages) == 'try-error') output <- unlist(any.ages)
  output
}

radio.chron <- unlist(sapply(geochron.records, get_ages))

hist(radio.chron[radio.chron<40000], breaks=seq(0, 25000, by = 1000),
     main = 'Radiocarbon dates for Pseudotsuga records',
     xlab = 'Radiocarbon date (14C years before 1950)')
}

}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.download}
\alias{get_site.download}
\title{Return Site Information from a \code{download}}
\usage{
\method{get_site}{download}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{download}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compile_taxa.R
\name{compile_taxa}
\alias{compile_taxa}
\title{Function to convert assemblage taxa to standardized lists.}
\usage{
compile_taxa(object, list.name, alt.table = NULL, cf = TRUE,
  type = TRUE)
}
\arguments{
\item{object}{A pollen object returned by \code{\link{get_download}}.}

\item{list.name}{The taxon compilation list, one of a set of lists from the literature (e.g., \code{"P25"}, \code{"WhitmoreFull"}).  More detail in section Details.}

\item{alt.table}{A user provided table formatted with at least two columns, one called 'taxon' and the other named as in \code{list.name}.}

\item{cf}{Should taxa listed as *cf*s (*e.g.*, *cf*. *Gilia*) be considered highly resolved?}

\item{type}{Should taxa listed as types (*e.g.*, *Iva annua*-type) be considered highly resolved?}
}
\value{
This command returns a list object with the same structure as the parent pollen object returned by \code{\link{get_download}}, or a matrix (or data frame) depending on whether \code{object} is one or the other.  Any pollen taxon not included in the major taxa defined in the pollen gets returned as 'Other'.
}
\description{
From the assemblage data for the core return assemblage data with the assemblage taxa
Currently implemented only for pollen data.
}
\details{
The data object uses the smaller pollen subset.  As this package develops we will add the capacity to summarize data output from the translation. Currently we can return only subsets that have been defined in the literature.  These lists include:
\itemize{
  \item{\code{"P25"} }{ This list is derived from Gavin et al., (2003), and includes 25 pollen taxa.}
  \item{\code{"WS64"} }{  This list is derived from Williams and Shuman (2008).}
  \item{\code{"WhitmoreFull"} }{  This is the full list associated with the Whitmore et al., (2005) North American Modern Pollen Database.}
  \item{\code{"WhitmoreSmall"} }{  As above, but taxa for which both fully resolved and undifferentiated exist these taxa are summed.}
}
}
\examples{
\dontrun{
#  Search for sites with "Thuja" pollen that are older than 8kyr BP and
#  that are on the west coast of North America:
t8kyr.datasets <- get_dataset(taxonname='Thuja*', loc=c(-150, 20, -100, 60), ageyoung = 8000)

#  Returns 3 records (as of 04/04/2013), get dataset for the first record, Gold Lake Bog.
GOLDKBG <- get_download(t8kyr.datasets[[1]])

gold.p25 <- compile_taxa(GOLDKBG, 'P25')

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org

Gavin DG, Oswald WW, Wahl ER, Williams JW. 2003. A statistical approach to 
  evaluating distance metrics and analog assignments for pollen records. 
  Quaternary Research 60: 356-367.

Whitmore J, Gajewski K, Sawada M, Williams JW, Shuman B, Bartlein PJ, Minckley T, 
  Viau AE, Webb III T, Shafer S, Anderson P, Brubaker L. 2005. Modern pollen data 
  from North America and Greenland for multi-scale paleoenvironmental applications. 
  Quaternary Science Reviews 24: 1828-1848.

Williams J, Shuman B. 2008. Obtaining accurate and precise environmental 
  reconstructions from the modern analog technique and North American surface pollen 
  dataset. Quaternary Science Reviews. 27:669-687.

API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication.download_list}
\alias{get_publication.download_list}
\title{A function to get publications for datasets in the Neotoma Database using the API.}
\usage{
\method{get_publication}{download_list}(x, ...)
}
\arguments{
\item{x}{an object of class \code{download_list}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
The function takes a \code{download_list} and returns a table with publication information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol}
\alias{get_chroncontrol}
\title{Function to return chronological control tables used to build age models.}
\usage{
get_chroncontrol(x, chronology = 1, verbose = TRUE, add = FALSE)
}
\arguments{
\item{x}{A single numeric chronology ID, a vector of numeric dataset IDs as returned by \code{\link{get_dataset}} or a \code{download} or \code{download_list} object.}

\item{chronology}{When \code{download} objects have more than associated chronology, which chronology do you want?  Default is \code{1}.}

\item{verbose}{logical, should messages on API call be printed?}

\item{add}{logical, should this chron control be added to the download object?}
}
\value{
This command returns either an object of class  \code{"try-error"} containing the error returned
   from the Neotoma API call, or a full data object containing all the relevant information required to build either the default or prior chronology for a core.
   When \code{download} or \code{download_list} objects are passes, the user can \code{add} the chroncontrol to the
   \code{download} object explicitly, in which case the function will return a download with \code{chroncontrol} embedded.
   
   This is a list comprising the following items:

 \item{ \code{chron.control} }{A table describing the collection, including dataset information, PI data compatable with \code{\link{get_contact}} and site data compatable with \code{\link{get_site}}.}
 \item{ \code{meta} }{Dataset information for the core, primarily the age-depth model and chronology.  In cases where multiple age models exist for a single record the most recent chronology is provided here.}
 
 If Neotoma returns empty content, either the control table or the associated metadata (which happens in approximately 25% of cases) then the data.frames are returned with NA content.
}
\description{
Using the dataset ID, return all records associated with the data.  At present,
   only returns the dataset in an unparsed format, not as a data table.   This function will only download one dataset at a time.
}
\examples{
\dontrun{
#  The point of pulling chronology tables is to re-build or examine the 
#  chronological information that was used to build the age-depth model for 
#  the core.  You can do this by hand, but the `write_agefile` function works 
#  with `download` objects directly.

three_pines <- get_download(get_dataset(get_site("Three Pines Bog"), 
                                        datasettype = "pollen"))
pines_chron <- get_chroncontrol(three_pines)

# Spline interpolation:
model <- smooth.spline(x = pines_chron[[1]]$chron.control$depth,
                       y = pines_chron[[1]]$chron.control$age)
                       
new_ages <- predict(model, x = three_pines[[1]]$sample.meta$depth)

}
}
\references{
+ Neotoma Project Website: http://www.neotomadb.org
+ API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Stratiplot.download.R
\name{Stratiplot.download}
\alias{Stratiplot.download}
\title{Palaeoecological stratigraphic diagrams}
\usage{
\method{Stratiplot}{download}(x, yaxis = "age", method = "none",
  group = NULL, ...)
}
\arguments{
\item{x}{A \code{download} object.}

\item{yaxis}{One of the columns in \code{sample.meta}, including \code{depth}, \code{age}, \code{age.younger}, or \code{age.older}, default \code{age}.}

\item{method}{An option for axis transformation using \code{tran} from the \code{analogue} package.  \code{"none"} by default.}

\item{group}{An ecological group from the taxon table.}

\item{...}{variables to be passed to \code{Stratiplot}.}
}
\value{
A \code{trellis} object.
}
\description{
Draws paleoecological diagrams from a \code{download} object.  Allows control of variable type (using the \code{tran} function from the \code{analogue} package), and taxonomic grouping.
}
\details{
A wrapper for the \code{analogue} package's \code{Stratiplot} function.  Allowing the user to plot a stratigraphic diagram directly from a \code{download} object.
}
\examples{
\dontrun{
lake_o_dl <- get_download(15925)
Stratiplot(lake_o_dl[[1]])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol.dataset_list}
\alias{get_chroncontrol.dataset_list}
\title{Function to return chronological control tables from a \code{dataset_list}.}
\usage{
\method{get_chroncontrol}{dataset_list}(x, chronology = 1,
  verbose = TRUE, add = FALSE)
}
\arguments{
\item{x}{A \code{dataset_list} object.}

\item{chronology}{When \code{download} objects have more than associated chronology, which chronology do you want?  Default is \code{1}.}

\item{verbose}{logical; should messages on API call be printed?}

\item{add}{Should the \code{chroncontrol} be added to the download object (only accepts \code{FALSE})}
}
\description{
Using a \code{dataset_list}, return the default chron-control table.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_download.R
\name{get_download}
\alias{get_download}
\title{Function to return full download records using \code{site}s, \code{dataset}s, or dataset IDs.}
\usage{
get_download(x, verbose = TRUE)
}
\arguments{
\item{x}{A single numeric dataset ID or a vector of numeric dataset IDs as returned by \code{get_datasets}, or a \code{site}, \code{dataset}, or \code{dataset_list}.}

\item{verbose}{logical; should messages on API call be printed?}
}
\value{
This command returns either object of class \code{"try-error"}' (see \code{\link{try}}) defined by the error returned from the Neotoma API call, or an object of class \code{download_list}, containing a set of \code{download} objects, each with relevant assemblage information and metadata:
The \code{download} object is a list of lists and data frames that describe an assemblage, the constituent taxa, the chronology, site and PIs who contributed the data. The following are important components:

 \item{ \code{dataset} }{A table describing the collection, including dataset information, PI data compatible with \code{\link{get_contact}} and site data compatable with \code{\link{get_site}}.}
 \item{ \code{sample.meta} }{Dataset information for the core, primarily the age-depth model and chronology.  In cases where multiple age models exist for a single record the most recent chronology is provided here.}
 \item{ \code{taxon.list} }{The list of taxa contained within the dataset, unordered, including information that can be used in \code{\link{get_taxa}}}
 \item{ \code{counts} }{The assemblage data for the dataset, arranged with each successive depth in rows and the taxa as columns.  All taxa are described in \code{taxon.list}, the chronology is in \code{sample.data}}
 \item{ \code{lab.data} }{A data frame of laboratory data, such as exotic pollen spike, amount of sample counted, charcoal counts, etc.}
 \item{ \code{chronologies} }{A list of existing chronologies.  If only a single chronology exists for a record then this is the same as the age-model in \code{sample.meta}.}
}
\description{
Using the dataset ID, site object or dataset object, return all records associated with the data as a \code{download_list}.
}
\section{Note}{

The function returns a warning in cases where single taxa are defined by multiple taphonomic characteristics, for example grains that are identified separately as crumpled and torn in the same sample and sums these values within a sample.
In the case that a geochronology dataset is passed to \code{get_download} the function returns a message and a NULL object (that is later excised).  Use \code{get_geochron} for these objects.
The chronologies can be augmented using the function \code{get_chroncontrol}, where the individual chronology objects in \code{chronologies} will consist of a table equivalent to \code{sample.meta} and a \code{chroncontrol} object.
}

\examples{
\dontrun{
#  Search for sites with "Pseudotsuga" pollen that are older than 8kyr BP and
#  that are roughly within western British Columbia:
t8kyr.datasets <- get_dataset(taxonname='*Picea*', loc=c(-90, 41, -89, 44),
                              ageold = 20000, ageyoung=10000)

#  Returns 20 records (as of 04/04/2013), get the dataset for all records:
pollen.records <- get_download(t8kyr.datasets)

#  Standardize the taxonomies for the different records using the WS64 taxonomy.
compiled.sites <- compile_taxa(pollen.records, list.name='WS64')

#  Extract the Pseudotsuga curves for the sites:
get.curve <- function(x, taxa) {
               if (taxa \%in\% colnames(x$counts)) {
                 count <- x$counts[,taxa]/rowSums(x$counts, na.rm=TRUE)
               } else {
                 count <- rep(0, nrow(x$count))
               }
               data.frame(site = x$dataset$site.data$site.name,
               age = x$sample.meta$age,
               count = count)
             }

curves <- do.call(rbind.data.frame,
                  lapply(compiled.sites, get.curve, taxa = 'Larix/Pseudotsuga'))

#  For illustration, remove the sites with no Pseudotsuga occurance:
curves <- curves[curves$count > 0, ]

smooth.curve <- predict(loess(sqrt(count)~age, data=curves),
                        data.frame(age=seq(20000, 0, by = -100)))

plot(sqrt(count) ~ age, data = curves,
     ylab = '\% Pseudotsuga/Larix', xlab='Calibrated Years BP', pch=19,
     col=rgb(0.1, 0.1, 0.1, 0.1), xlim=c(0, 20000))
lines(seq(20000, 0, by = -100), smooth.curve, lwd=2, lty=2, col=2)

#  This figure shows us an apparent peak in Larix/Pseudotsuga pollen in the
#  early-Holocene that lends support to a warmer, drier early-Holocene in
#  western North America.
}

}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse.download_list}
\alias{browse.download_list}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
\method{browse}{download_list}(x)
}
\arguments{
\item{x}{A \code{download_list} object.}
}
\description{
Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse.default}
\alias{browse.default}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
\method{browse}{default}(x)
}
\arguments{
\item{x}{A numeric value with the dataset ID.}
}
\description{
Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_table.R
\name{get_table}
\alias{get_table}
\title{Get Neotoma value tables.}
\usage{
get_table(table.name = NULL)
}
\arguments{
\item{table.name}{Call one of the available tables in the Neotoma Database.
A full listing of tables can be found here: \url{http://api.neotomadb.org/doc/resources/dbtables}.
By default it returns all objects in the table.}
}
\description{
Get Neotoma value tables.
}
\details{
A table of values corresponding to the parameter of interest.
}
\examples{
\dontrun{
taxon.table <- get_table('Taxa')

#  Get the frequency of a random taxon in Neotoma.
tax_sample <- sample(nrow(taxon.table), 1)
cat("The taxon", 
    taxon.table$TaxonName[tax_sample], 
    "occurs in Neotoma", 
    length(get_dataset(taxonname = taxon.table$TaxonName[tax_sample])), 
    "times.")
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa.R
\name{taxa}
\alias{taxa}
\alias{taxa.download}
\alias{taxa.download_list}
\title{Access proxy taxonomic data}
\usage{
taxa(obj, ...)

\method{taxa}{download}(obj, ...)

\method{taxa}{download_list}(obj, collapse = TRUE, hierarchy = FALSE,
  ...)
}
\arguments{
\item{obj}{an R object from which counts are to be extracted.}

\item{...}{arguments passed to other methods.}

\item{collapse}{should the results be returned as a list, one for each site (\code{FALSE}), or a single dataframe of all taxa? Default is \code{TRUE}}

\item{hierarchy}{Should the taxonomic hierarchy be included?}
}
\value{
Either a data frame of taxa or a list of such objects.
}
\description{
Extracts taxa from \code{download} objects and returns them in a useful format.
}
\details{
Methods are available for "download" and "download_list" objects.
}
\examples{
\dontrun{
ostracodes <- get_dataset(datasettype = 'ostracode')

ostro.dl <- get_download(ostracodes)
ostro.taxa <- taxa(ostro.dl)
}
}
\author{
Simon Goring
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.download}
\alias{get_dataset.download}
\title{Obtain dataset information from an existing \code{download} object.}
\usage{
\method{get_dataset}{download}(x, ...)
}
\arguments{
\item{x}{An object of class \code{download}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access a \code{dataset} within a \code{download} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.integer}
\alias{get_site.integer}
\title{Return Site Information from a vector of integers.}
\usage{
\method{get_site}{integer}(sitename, ...)
}
\arguments{
\item{sitename}{An integer or vector of integers.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ages.R
\name{ages}
\alias{ages}
\alias{ages.download}
\alias{ages.download_list}
\title{Access proxy age data}
\usage{
ages(obj, ...)

\method{ages}{download}(obj, ...)

\method{ages}{download_list}(obj, ...)
}
\arguments{
\item{obj}{an R object from which counts are to be extracted.}

\item{...}{arguments passed to other methods.}
}
\value{
Either a data frame of ages or a list of such objects.
}
\description{
Extracts age information from objects and returns them in a useful format.
}
\details{
Methods are available for "download" and "download_list" objects.
}
\examples{
\dontrun{
ostracodes <- get_dataset(datasettype = 'ostracode')

ostro.dl <- get_download(ostracodes)
ostro.ages <- ages(ostro.dl)
}
}
\author{
Simon Goring
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.download_list}
\alias{get_site.download_list}
\title{Return Site Information from a \code{download_list}}
\usage{
\method{get_site}{download_list}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{download_list}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_taxa.R
\name{get_taxa}
\alias{get_taxa}
\title{Get taxon information from Neotoma.}
\usage{
get_taxa(taxonid, taxonname, status, taxagroup, ecolgroup)
}
\arguments{
\item{taxonid}{Numeric taxon identifier used in Neotoma}

\item{taxonname}{A character string representing the full or partial name of taxa of interest.}

\item{status}{The current status of the taxon, one of 'extinct', 'extant', 'all'.}

\item{taxagroup}{The taxonomic grouping for the taxa. See \url{http://api.neotomadb.org/doc/resources/taxa} for the list of approved groupings.}

\item{ecolgroup}{The ecological group of the taxa. More detailed than \code{taxagroup}, can be obtained using \code{get_table("EcolGroupTypes")}.}
}
\value{
Returns a data frame with the following components:

 \item{ \code{TaxonID} }{Unique database record identifier for a taxon}
 \item{ \code{TaxonCode} }{Shorthand notation for a taxon identification}
 \item{ \code{TaxonName} }{Name of the taxon}
 \item{ \code{Author} }{Author(s) of the name. Used almost exclusively with beetle taxa}
 \item{ \code{Extinct} }{True if extinct; false if extant}
 \item{ \code{TaxaGroup} }{Code for taxa group to which taxon belongs}
 \item{ \code{EcolGroups} }{Array of ecological group codes to which the taxon belongs}
 \item{ \code{HigherTaxonID} }{TaxonID of the next higher taxonomic rank}
 \item{ \code{PublicationID} }{Publication identification number}
 \item{ \code{Notes} }{Free-form notes or comments about the taxon}
}
\description{
Get taxon information from Neotoma.
}
\examples{
\dontrun{
## Return all species taxa with "Abies" in name - note wildcard
taxa <- get_taxa(taxonname = "Abies*")
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_download.R
\name{get_download.site}
\alias{get_download.site}
\title{Function to return full download records using a \code{site}.}
\usage{
\method{get_download}{site}(x, verbose = TRUE)
}
\arguments{
\item{x}{An object of class \code{site}.}

\item{verbose}{logical; should messages on API call be printed?}
}
\description{
Using a \code{site}, return all records associated with the data as a \code{download_list}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_download.R
\name{get_download.dataset_list}
\alias{get_download.dataset_list}
\title{Function to return full download records using a \code{dataset_list}.}
\usage{
\method{get_download}{dataset_list}(x, verbose = TRUE)
}
\arguments{
\item{x}{An object of class \code{dataset_list}.}

\item{verbose}{logical; should messages on API call be printed?}
}
\description{
Using a \code{dataset_list}, return all records associated with the data as a \code{download_list}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication}
\alias{get_publication}
\title{A function to get publications for sites or datasets in the Neotoma Database using the API.}
\usage{
get_publication(x, contactid, datasetid, author, pubtype, year, search)
}
\arguments{
\item{x}{Numeric Publication ID value, either from \code{\link{get_dataset}} or known.}

\item{contactid}{Numeric Contact ID value, either from \code{\link{get_dataset}} or \code{\link{get_contact}}}

\item{datasetid}{Numeric Dataset ID, known or from \code{\link{get_dataset}}}

\item{author}{Character string for full or partial author's name.  Can include wildcards such as 'Smit*' for all names beginning with 'Smit'.}

\item{pubtype}{Character string, one of eleven allowable types, see \code{\link{get_table}}. For a list of allowed types run \code{get_table("PublicationTypes")}.}

\item{year}{Numeric publication year.}

\item{search}{A character string to search for within the article citation.}
}
\value{
A list is returned with two data frame components:

 \item{ \code{meta} }{A single row with Publication ID, type, year of publication and full citation.}
 \item{ \code{Authors} }{\code{data.frame} of author names, order and IDs, can be of variable length.}
}
\description{
The function takes the parameters, defined by the user, and returns a table with publication information from the Neotoma Paleoecological Database.
}
\examples{
\dontrun{
#  To find all publications from 1998:
year.cont <- get_publication(year = 1998)

# To find all data contributors who have the last name "Smith"
smith.cont <- get_publication(author = 'Smith')
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_download.R
\name{get_download.dataset}
\alias{get_download.dataset}
\title{Function to return full download records using a \code{dataset}.}
\usage{
\method{get_download}{dataset}(x, verbose = TRUE)
}
\arguments{
\item{x}{An object of class \code{dataset}.}

\item{verbose}{logical; should messages on API call be printed?}
}
\description{
Using a \code{dataset}, return all records associated with the data as a \code{download_list}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon.table.R
\docType{data}
\name{taxon.list}
\alias{taxon.list}
\title{Neotoma taxon list}
\format{a \code{data.frame} object}
\source{
The Neotoma database.
}
\usage{
taxon.list
}
\description{
The taxonomy table for datasets in neotoma, as would be returned by \code{\link{get_table}}
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.download_list}
\alias{get_dataset.download_list}
\title{Obtain dataset information from a \code{download_list}.}
\usage{
\method{get_dataset}{download_list}(x, ...)
}
\arguments{
\item{x}{An object of class \code{download_list}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to return datasets corresponding to the objects within a \code{download_list}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol.dataset}
\alias{get_chroncontrol.dataset}
\title{Function to return chronological control tables from a \code{dataset}.}
\usage{
\method{get_chroncontrol}{dataset}(x, chronology = 1, verbose = TRUE,
  add = FALSE)
}
\arguments{
\item{x}{A \code{dataset}.}

\item{chronology}{When \code{download} objects have more than associated chronology, which chronology do you want?  Default is \code{1}.}

\item{verbose}{logical; should messages on API call be printed?}

\item{add}{Should the \code{chroncontrol} be added to the download object (only accepts \code{FALSE})}
}
\description{
Using a \code{dataset}, return the default chron-control table.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/depths.R
\name{depths}
\alias{depths}
\alias{depths.default}
\alias{depths.download}
\alias{depths.download_list}
\title{Extracts the depth values from a `download` object}
\usage{
depths(obj, ...)

\method{depths}{default}(obj, ...)

\method{depths}{download}(obj, ...)

\method{depths}{download_list}(obj, ...)
}
\arguments{
\item{obj}{A \code{download} object.}

\item{...}{arguments passed to other methods.}
}
\value{
Returns a vector of depths.
}
\description{
Using a \code{download} object, return the sample depths (if available).

Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
\examples{
\dontrun{
# Provide a vector of depths to generate a new age model:
# The dataset id 684 is for Devils Lake, a record published by Louis Maher Jr.

pollen.data <- get_download(684)
pollen.chron <- get_chroncontrol(pollen.data)[[1]]

age_sds <- pollen.chron$chron.control$age - focal$chron.control$age.young,
get_curves <- ifelse(regexpr("Radiocarbon",
                             pollen.chron$chron.control$control.type) > -1, 
                     'intcal13', 'normal')

new_chron <- Bchron::Bchronology(ages   = pollen.chron$chron.control$age,
                                 ageSds = age_sds
                                 positions = pollen.chron$chron.control$depth,
                                 calCurves = , 
                                 predictPositions = depths(pollen.data))

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/sites
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset}
\alias{get_dataset}
\title{Obtain dataset information from the Neotoma Paleoecological Database or an existing object.}
\usage{
get_dataset(x, datasettype, piid, altmin, altmax, loc, gpid, taxonids,
  taxonname, ageold, ageyoung, ageof, subdate)
}
\arguments{
\item{x}{An optional value, either a \code{numeric} site ID or object of class \code{download}, \code{download_list} or \code{site}.}

\item{datasettype}{A character string corresponding to one of the allowed dataset types in the Neotoma Database.  Allowed types include: \code{"geochronologic"}, \code{"loss-on-ignition"}, \code{"pollen"}, \code{"plant macrofossils"}, \code{"vertebrate fauna"}, \code{"mollusks"}, and \code{"pollen surface sample"}.}

\item{piid}{Numeric value for the Principle Investigator's ID number.}

\item{altmin}{Numeric value indicating the minimum altitude for the site (can be used alone or with \code{altmax}).}

\item{altmax}{Numeric value indicating the maximum altitude for the site (can be used alone or with \code{altmin}).}

\item{loc}{A numeric vector \code{c(lonW, latS, lonE, latN)} representing the bounding box within which to search for sites.  The convention here is to use negative values for longitudes west of Greenwich or longitudes south of the equator}

\item{gpid}{A character string or numeric value, must correspond to a valid geopolitical identity in the Neotoma Database.  Use get.tables('GeoPoliticalUnits') for a list of acceptable values, or link here: \url{http://api.neotomadb.org/apdx/geopol.htm}}

\item{taxonids}{A numeric identifier for the taxon.  See \code{\link{get_table}} and use \code{get_tables('Taxa')} for a list of acceptable values.}

\item{taxonname}{A character string corresponding to a valid taxon identity in the Neotoma Database.  See \code{\link{get_table}} and use \code{get_table('Taxa')} for a list of acceptable values.}

\item{ageold}{The oldest date acceptable for the search (in years before present).}

\item{ageyoung}{The youngest date acceptable for the search.}

\item{ageof}{If a taxon ID or taxon name is defined this parameter must be set to \code{"taxon"}, otherwise it may refer to \code{"sample"}, in which case the age bounds are for any samples within datasets or \code{"dataset"} if you want only datasets that are within the bounds of ageold and ageyoung.}

\item{subdate}{Date of dataset submission, either YYYY-MM-DD or MM-DD-YYYY.}
}
\value{
More details on the use of these parameters can be obtained from
   \url{http://api.neotomadb.org/doc/resources/datasets}.

   A list of class `dataset_list`, with each item corresponding to an individual record.
   Searches that return no items will result in a NULL value being returned.
   Otherwise each list item (each dataset record) includes the following components:

 \item{ \code{dataset.id} }{Unique database record identifier for the dataset.}
 \item{ \code{dataset.name}  }{Name of the dataset; not commonly used.}
 \item{ \code{CollUnitHandle}  }{Code name of the Collection Unit with which the dataset is associated. This code may be up to 10 characters. Data are frequently distributed by Collection Unit, and the Handle is used for file names.}
 \item{ \code{CollUnitID}  }{Unique database record identifier for the collection unit.}
 \item{ \code{CollType}  }{The collection type. Types include cores, sections, excavations, and animal middens.}
 \item{ \code{DatasetType}  }{The dataset type, such as: geochronologic, loss-on-ignition, pollen, plant macrofossils, vertebrate fauna, etc.}
 \item{ \code{AgeOldest}  }{The oldest of all sample ages (in calendar years before present) in the dataset.}
 \item{ \code{AgeYoungest}  }{The youngest of all sample ages (in calendar years before present) in the dataset.}
 \item{ \code{SubDates}  }{An array of objects that describe dataset submission events.  If multiple submissions occurred then this is a table.}
 \item{ \code{DatasetPIs}  }{An array of objects that describe Principal Investigators associated with a dataset.}
 \item{ \code{Site}  }{An object describing the site where the dataset samples were taken.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
\examples{
\dontrun{
# Search for sites with "Thuja" pollen that are older than 8kyr BP and
# that are on the west coast of North America:
t8kyr.datasets <- get_dataset(taxonname='Thuja*', 
                              loc=c(-150, 20, -100, 60), 
                              ageyoung = 8000)

# Search for vertebrate fossils in Canada (gpid: 756) within the last 2kyr.
gpids <- get_table(table.name='GeoPoliticalUnits')
canID <- gpids[which(gpids$GeoPoliticalName == 'Canada'),1]

v2kyr.datasets <- get_dataset(datasettype='vertebrate fauna', 
                              gpid=canID, 
                              ageold = 2000)
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site}
\alias{get_site}
\title{Return Site Information.}
\usage{
get_site(sitename, altmin, altmax, loc, gpid, ...)
}
\arguments{
\item{sitename}{character string representing the full or partial site name, or an object of class \code{dataset}, \code{dataset_list}, \code{download} or \code{download_list}}

\item{altmin}{Minimum site altitude  (in m).}

\item{altmax}{Maximum site altitude (in m).}

\item{loc}{A numeric vector c(lonW, latS, lonE, latN) representing the bounding box within which to search for sites.  The convention here is to use negative values for longitudes west of Grewnwich or longitudes south of the equator.}

\item{gpid}{A character string or numeric value, must correspond to a valid geopolitical identity in the Neotoma Database.  Use get.tables('GeoPoliticalUnits') for a list of acceptable values, or link here: http://api.neotomadb.org/apdx/geopol.htm}

\item{...}{Optional additional arguments}
}
\value{
A data frame:

 \item{\code{siteid}}{Unique database record identifier for the site.}
 \item{\code{sitename}}{Name of the site.}
 \item{\code{long}}{Mean longitude, in decimal degrees, for a site (-180 to 180).}
 \item{\code{lat}}{Mean latitude, in decimal degrees, for a site (-90 to 90).}
 \item{\code{elev}}{Elevation in meters.}
 \item{\code{description}}{Free form description of a site, including such information as physiography and vegetation around the site.}
 \item{\code{long_acc}}{If the site is described by a bounding box this is the box width.}
 \item{\code{lat_acc}}{If the site is described by a bounding box this is the box height.}
}
\description{
Return site information from the Neotoma Paleoecological Database.

\code{get_site} returns site information from the Neotoma Paleoecological Database
   based on parameters defined by the user.
}
\examples{
\dontrun{
#  What is the distribution of site elevations in Neotoma?
all.sites <- get_site()  #takes a bit of time.

plot(density(all.sites$elev, from = 0, na.rm=TRUE),
main = 'Altitudinal Distribution of Neotoma Sites', xlab = 'Altitude (m)', log='x')

#  Get site information from a dataset:
nw.datasets <- get_dataset(loc = c(-140, 50, -110, 65), 
                           datasettype='pollen',
                           taxonname='Pinus*')
                           
nw.sites <- get_site(nw.datasets)

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://api.neotomadb.org/doc/resources/sites
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.geochronologic_list}
\alias{get_site.geochronologic_list}
\title{Return Site Information from a \code{geochronologic_list}}
\usage{
\method{get_site}{geochronologic_list}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{geochronologic_list}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.geochronologic_list}
\alias{get_dataset.geochronologic_list}
\title{Obtain dataset information from an object of class \code{geochronologic_list}.}
\usage{
\method{get_dataset}{geochronologic_list}(x, ...)
}
\arguments{
\item{x}{An object of class \code{geochronologic_list}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.geochronologic}
\alias{get_site.geochronologic}
\title{Return Site Information from a \code{geochronologic}}
\usage{
\method{get_site}{geochronologic}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{geochronologic}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.geochronologic}
\alias{get_dataset.geochronologic}
\title{Obtain dataset information from an object of class \code{geochronologic}.}
\usage{
\method{get_dataset}{geochronologic}(x, ...)
}
\arguments{
\item{x}{An object of class \code{geochronologic}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication.dataset_list}
\alias{get_publication.dataset_list}
\title{A function to get publications for dataset_lists in the Neotoma Database using the API.}
\usage{
\method{get_publication}{dataset_list}(x, ...)
}
\arguments{
\item{x}{an object of class \code{dataset_list}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
The function takes a \code{dataset_list} and returns a table with publication information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.default}
\alias{get_dataset.default}
\title{Obtain dataset information from the Neotoma Paleoecological Database or an existing object.}
\usage{
\method{get_dataset}{default}(x, datasettype, piid, altmin, altmax, loc,
  gpid, taxonids, taxonname, ageold, ageyoung, ageof, subdate)
}
\arguments{
\item{x}{A numeric value corresponding to the site ID.}

\item{datasettype}{A character string corresponding to one of the allowed dataset types in the Neotoma Database.  Allowed types include: \code{"geochronologic"}, \code{"loss-on-ignition"}, \code{"pollen"}, \code{"plant macrofossils"}, \code{"vertebrate fauna"}, \code{"mollusks"}, and \code{"pollen surface sample"}.}

\item{piid}{Numeric value for the Principle Investigator's ID number.}

\item{altmin}{Numeric value indicating the minimum altitude for the site (can be used alone or with \code{altmax}).}

\item{altmax}{Numeric value indicating the maximum altitude for the site (can be used alone or with \code{altmin}).}

\item{loc}{A numeric vector \code{c(lonW, latS, lonE, latN)} representing the bounding box within which to search for sites.  The convention here is to use negative values for longitudes west of Greenwich or longitudes south of the equator}

\item{gpid}{A character string or numeric value, must correspond to a valid geopolitical identity in the Neotoma Database.  Use get.tables('GeoPoliticalUnits') for a list of acceptable values, or link here: \url{http://api.neotomadb.org/apdx/geopol.htm}}

\item{taxonids}{A numeric identifier for the taxon.  See \code{\link{get_table}} and use \code{get_tables('Taxa')} for a list of acceptable values.}

\item{taxonname}{A character string corresponding to a valid taxon identity in the Neotoma Database.  See \code{\link{get_table}} and use \code{get_table('Taxa')} for a list of acceptable values.}

\item{ageold}{The oldest date acceptable for the search (in years before present).}

\item{ageyoung}{The youngest date acceptable for the search.}

\item{ageof}{If a taxon ID or taxon name is defined this parameter must be set to \code{"taxon"}, otherwise it may refer to \code{"sample"}, in which case the age bounds are for any samples within datasets or \code{"dataset"} if you want only datasets that are within the bounds of ageold and ageyoung.}

\item{subdate}{Date of dataset submission, either YYYY-MM-DD or MM-DD-YYYY.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication.dataset}
\alias{get_publication.dataset}
\title{A function to get publications for datasets in the Neotoma Database using the API.}
\usage{
\method{get_publication}{dataset}(x, ...)
}
\arguments{
\item{x}{an object of class \code{dataset}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
The function takes a \code{dataset} and returns a table with publication information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication.default}
\alias{get_publication.default}
\title{A function to get publications for sites or datasets in the Neotoma Database using the API.}
\usage{
\method{get_publication}{default}(x, contactid, datasetid, author, pubtype,
  year, search)
}
\arguments{
\item{x}{Numeric Publication ID value, either from \code{\link{get_dataset}} or known.}

\item{contactid}{Numeric Contact ID value, either from \code{\link{get_dataset}} or \code{\link{get_contact}}}

\item{datasetid}{Numeric Dataset ID, known or from \code{\link{get_dataset}}}

\item{author}{Character string for full or partial author's name.  Can include wildcards such as 'Smit*' for all names beginning with 'Smit'.}

\item{pubtype}{Character string, one of eleven allowable types, see \code{\link{get_table}}. For a list of allowed types run \code{get_table("PublicationTypes")}.}

\item{year}{Numeric publication year.}

\item{search}{A character string to search for within the article citation.}
}
\description{
The function takes the parameters, defined by the user, and returns a table with
   publication information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol.download_list}
\alias{get_chroncontrol.download_list}
\title{Function to return chronological control tables from a \code{download_list} object.}
\usage{
\method{get_chroncontrol}{download_list}(x, chronology = 1, verbose = TRUE, add = FALSE)
}
\arguments{
\item{x}{A \code{download_list} object.}

\item{chronology}{When \code{download} objects have more than associated chronology, which chronology do you want?  Default is \code{1}.}

\item{verbose}{logical; should messages on API call be printed?}

\item{add}{Should the \code{chroncontrol} be added to the download object (default \code{FALSE})}
}
\description{
Using a \code{download_list}, return the default chron-control table as a \code{data.frame}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse}
\alias{browse}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
browse(x)
}
\arguments{
\item{x}{A numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object.}
}
\value{
Returns a NULL value, opens a browser.
}
\description{
Using a \code{download} or \code{dataset} object, open up a browser window in the users default browser. Passing a \code{download_list} or \code{dataset_list} will open Neotoma Explorer with the first object and return a warning.

Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
\examples{
\dontrun{
# Where are the XRF data?

xrf.data <- get_dataset(datasettype='X-ray fluorescence (XRF)')
browse(xrf.data)

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/sites
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_methods.R
\name{read.tilia}
\alias{read.tilia}
\title{Read proxy data from Tilia TLX files}
\usage{
read.tilia(file)
}
\arguments{
\item{file}{a string representing a Tilia TLX format file.}
}
\value{
Return a `download` object.
}
\description{
Read proxy data from a Tilia TLX format file.
}
\examples{
\dontrun{
  crystal <- read.tilia('crystal.tlx')
}

}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.numeric}
\alias{get_site.numeric}
\title{Return Site information from a vector of numeric elements.}
\usage{
\method{get_site}{numeric}(sitename = NULL, ...)
}
\arguments{
\item{sitename}{A numeric value or vector of numeric elements.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol.default}
\alias{get_chroncontrol.default}
\title{Function to return chronological control tables from a chronologic ID.}
\usage{
\method{get_chroncontrol}{default}(x, chronology = 1, verbose = TRUE, add = FALSE)
}
\arguments{
\item{x}{A single numeric chronology ID or a vector of numeric chronology IDs as returned by \code{get_datasets}.}

\item{chronology}{For \code{download} methods, which chronology controls should be used?}

\item{verbose}{logical; should messages on API call be printed?}

\item{add}{logical, should this chron control be added to the download object?}
}
\description{
Using the chronology ID, return the chron control table as a \code{data.frame}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param_check.R
\name{param_check}
\alias{param_check}
\title{Internal function to check passed parameters.}
\usage{
param_check(cl)
}
\arguments{
\item{cl}{Contact ID is a numerical value associated with the Neotoma
Contact table's numerical Contact ID.}
}
\value{
A list with two components:

 \item{flag}{Returns a 0 if everything's fine, a 1 if there's a problem.}
 \item{message}{A list of error messages.}
}
\description{
Functions \code{\link{get_site}}, \code{\link{get_dataset}} and others pass parameters to \code{param_check}, which tells them if there's a problem.
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{internal}
\keyword{misc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_contact.R
\name{get_contact}
\alias{get_contact}
\title{Get contact information.}
\usage{
get_contact(contactid, contactname, contactstatus, familyname)
}
\arguments{
\item{contactid}{Contact ID is a numerical value associated with the Neotoma
Contact table's numerical Contact ID.}

\item{contactname}{A character string indicating the data contributors' project,
organization or personal name.  May be a partial string and can include wildcards.}

\item{contactstatus}{The current status of the contact.  Possible values include:
active, deceased, defunct, extant, inactive, retired, unknown.}

\item{familyname}{A character string.  Full or partial string indicating the
contact's last name.}
}
\value{
The function takes parameters defined by the user and returns a list
   of contact information supplied by the Neotoma Paleoecological Database.
   The user may define all or none of the possible fields.  The function contains
   data checks for each defined parameter.

   The function returns either a single item of class \code{"try-error"} describing
   the reason for failure (either mis-defined parameters or an error from the Neotoma API),
   or a table of contacts, with rows corresponding to the number of individual
   contacts returned by the Neotoma API.  Each row entry includes the following parameters:

 \item{ \code{contact.name} }{Full name of the person, last name first (e.g. \code{"Simpson, George Gaylord"}) or name of organization or project (e.g. \code{"Great Plains Flora Association"}).}
 \item{ \code{contact.status} }{Current status of the person, organization, or project. Field links to the ContactStatuses lookup table.}
 \item{ \code{family.name} }{Family or surname name of a person.}
 \item{ \code{leading.initials} }{Leading initials for given or forenames without spaces (e.g. \code{"G.G."}).}
 \item{ \code{given.names} }{Given or forenames of a person (e.g. \code{"George Gaylord"}). Initials with spaces are used if full given names are not known (e.g. \code{"G. G")}.}
 \item{ \code{suffix} }{Suffix of a person's name (e.g. \code{"Jr."}, \code{"III"}).}
 \item{ \code{title} }{A person's title (e.g. \code{"Dr."}, \code{"Prof."}, \code{"Prof. Dr"}).}
 \item{ \code{phone} }{Telephone number.}
 \item{ \code{fax} }{Fax number.}
 \item{ \code{email} }{Email address.}
 \item{ \code{url} }{Universal Resource Locator, an Internet World Wide Web address.}
 \item{ \code{address} }{Full mailing address.}
 \item{ \code{notes} }{Free form notes or comments about the person, organization, or project.}
 \item{ \code{contact.id} }{Unique database record identifier for the contact.}
 \item{ \code{alias.id} }{The ContactID of a person's current name. If the AliasID is different from the ContactID, the ContactID refers to the person's former name.}
}
\description{
A function to obtain contact information for data contributors from the Neotoma Paleoecological Database.
}
\examples{
\dontrun{
#  To find all data contributors who are active:
active.cont <- get_contact(contactstatus = 'active')

# To find all data contributors who have the last name "Smith"
smith.cont <- get_contact(familyname = 'Smith')
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse.dataset}
\alias{browse.dataset}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
\method{browse}{dataset}(x)
}
\arguments{
\item{x}{A \code{dataset} object.}
}
\description{
Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.site}
\alias{get_dataset.site}
\title{Obtain dataset information from an existing \code{site} object.}
\usage{
\method{get_dataset}{site}(x, ...)
}
\arguments{
\item{x}{An object of class \code{site}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_leaflet.R
\name{plot_leaflet}
\alias{plot_leaflet}
\title{Leaflet plots for neotoma data.}
\usage{
plot_leaflet(x, providerTiles = "Stamen.TerrainBackground", ...)
}
\arguments{
\item{x}{A neotoma data object}

\item{providerTiles}{Default "Stamen.TerrainBackground", a character string indicating the tile background to be used for plotting.}

\item{...}{Other terms to be passed to the function.}
}
\value{
A \code{leaflet} object
}
\description{
A plotting function to provide interactive data investigation using the leaflet tools.
  This package requires a connection to the internet for proper functioning.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_closest.R
\name{get_closest}
\alias{get_closest}
\title{Find the closest dataset records to a site, dataset or long/lat pair in Neotoma}
\usage{
get_closest(x, n, buffer, ...)
}
\arguments{
\item{x}{A vector long/lat pair, or a dataset, site or download.}

\item{n}{The maximum number of records to return (in the case of ties the return may be larger)}

\item{buffer}{The size of the buffer for dataset search (in meters)}

\item{...}{optional arguments to pass into \code{get_dataset}.}
}
\value{
This command returns a \code{dataset} or \code{dataset_list}, or NULL if no records exist within the bounding box.
}
\description{
Passing in a download object the function outputs a Bacon or Clam formatted file to a
user defined destination for age modelling with existing age-depth modeling software.
}
\details{
The function uses the \code{sf} package to generate a circular buffer around a point of interest.
From there a square bounding box is sent to Neotoma using the \code{get_dataset()} function.  To use the buffering
function we must convert from long/lat to UTM coordinates, which we do by guessing the UTM zone of the point of interest.
Details can be found in the function's R code hosted on GitHub: \url{https://github.com/ropensci/neotoma/blob/master/R/get_closest.R}
}
\examples{
\dontrun{
#  The point of pulling chronology tables is to re-build or examine the chronological
#  information that was used to build the age-depth model for the core.
# Find the closest records to Madison, WI:
get_closest(x = c(-89.4012, 43.0731), n = 10, buffer = 5000, datasettype = "pollen")
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}, Andria Dawson \email{andria.dawson@gmail.com}
}
\keyword{API}
\keyword{Neotoma}
\keyword{Palaeoecology}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_agefile.R
\name{write_agefile}
\alias{write_agefile}
\title{Write age control file to disk formatted for either Bacon or Clam}
\usage{
write_agefile(download, chronology = 1, path, corename, cal.prog = "Bacon")
}
\arguments{
\item{download}{A single site returned by \code{get_download}.}

\item{chronology}{Default is \code{1}, the default chronology for the core.  If a core has more than one chronology the user can define a different set of chronological controls.}

\item{path}{The location of the 'Cores' folder & working directory for Bacon.  Do not include "Cores" in the path name.}

\item{corename}{The intended handle for the core, to be used in writing to file.}

\item{cal.prog}{The method intended to build the age model, either \code{'Bacon'} or \code{'Clam'}.}
}
\value{
This command returns a file in location \code{path/Cores} containing all the relevant information required to build either the default or prior chronology for a core.
}
\description{
Passing in a download object the function outputs a Bacon or Clam formatted file to a
user defined destination for age modelling with existing age-depth modeling software.
}
\examples{
\dontrun{
# Find a particular record:

three_pines <- get_download(get_dataset(get_site("Three Pines Bog"), 
                                        datasettype = "pollen"))

# You will need to edit the `path` argument here to point to a directory that 
# contains a `Cores` directory.

write_agefile(download = three_pines[[1]], 
              path = "./inst", 
              corename = "THREEPINES", 
              cal.prog = "Bacon")
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{API}
\keyword{Neotoma}
\keyword{Palaeoecology}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse.dataset_list}
\alias{browse.dataset_list}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
\method{browse}{dataset_list}(x)
}
\arguments{
\item{x}{A \code{dataset_list} object.}
}
\description{
Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.default}
\alias{get_site.default}
\title{Return Site Information.}
\usage{
\method{get_site}{default}(sitename, ...)
}
\arguments{
\item{sitename}{A character string representing the full or partial site name.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.integer}
\alias{get_dataset.integer}
\title{Obtain dataset information from a vector of dataset IDs.}
\usage{
\method{get_dataset}{integer}(x = NULL, ...)
}
\arguments{
\item{x}{A single numeric dataset id, or a numeric vector.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_download.R
\name{get_download.default}
\alias{get_download.default}
\title{Function to return full download records using \code{numeric} dataset IDs.}
\usage{
\method{get_download}{default}(x, verbose = TRUE)
}
\arguments{
\item{x}{A single numeric dataset ID or a vector of numeric dataset IDs as returned by \code{get_datasets}.}

\item{verbose}{logical; should messages on API call be printed?}
}
\description{
Using the dataset ID, return all records associated with the data as a \code{download_list}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.numeric}
\alias{get_dataset.numeric}
\title{Obtain dataset information from a vector of dataset IDs.}
\usage{
\method{get_dataset}{numeric}(x = NULL, ...)
}
\arguments{
\item{x}{A single numeric dataset id, or a numeric vector.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol.download}
\alias{get_chroncontrol.download}
\title{Function to return chronological control tables from a \code{download} object.}
\usage{
\method{get_chroncontrol}{download}(x, chronology = 1, verbose = TRUE, add = FALSE)
}
\arguments{
\item{x}{A single \code{download} object.}

\item{chronology}{For \code{download} methods, which chronology controls should be used?}

\item{verbose}{logical; should messages on API call be printed?}

\item{add}{Should the \code{chroncontrol} be added to the download object (default \code{FALSE})}
}
\description{
Using a \code{download}, return the default chron-control table as a \code{data.frame}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_bacon.R
\name{read_bacon}
\alias{read_bacon}
\title{Function to read in defined Bacon outputs.}
\usage{
read_bacon(
  x,
  path = ".",
  add = FALSE,
  chron_name = "Bacon",
  as_default = TRUE,
  download = NULL,
  sections = NULL,
  age_field = "median",
  interp = TRUE
)
}
\arguments{
\item{x}{A folder path that contains a Bacon \code{age} file.}

\item{path}{The location of the \code{Cores} folder.}

\item{add}{Should the results be added to an existing \code{download}? Defaults to \code{FALSE}.}

\item{chron_name}{The name for the chronology if the Bacon file is being added to a \code{download}.}

\item{as_default}{Should the chronology become the default?}

\item{download}{The target \code{download} if \code{add} is \code{TRUE}.}

\item{sections}{If there are multiple Bacon runs in a folder, identify the file by the number of sections in the run.}

\item{age_field}{Should the age be assigned to the \code{"median"} or the \code{"wmean"}?}

\item{interp}{If the depths don't match up, should we interpolate from the Bacon output? (default \code{TRUE})}
}
\description{
Reads in Bacon output and formats it for inclusion in a download object.
}
\details{
The function expects that you are in a working directory containing a "Cores" which would then contain output files from Bacon runs.  The output can either be added to an existing record (for example, replacing the default age model returned by Neotoma), or it can be loaded on its own.
If the depths for the loaded file do not match with the depths in the `download`'s `sample.meta` then the user can use the `interp` parameter to interpolate between depths.  This method uses linear interpolation.
}
\examples{
\dontrun{
# Download the record for Lake O' Pines:
lake_o_dl <- get_download(15925)

# This assumes that you have Bacon installed in a folder and have
# set it to your working directory.

write_agefile(lake_o_dl[[1]], path = ".", chronology = 1, 
              corename = "LAKEPINES", cal.prog = 'Bacon') 

source("Bacon.R") 

# These defaults just help the core run quickly, they're not 
# neccesarily good parameters.

Bacon("LAKEPINES", acc.mean = 10, 
      thick = 50, depths.file = TRUE, 
      suggest = FALSE, ask = FALSE)

lake_o_dl <- read_bacon("LAKEPINES", add = TRUE, 
                        download = download, sections = 17)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gp.table.R
\docType{data}
\name{gp.table}
\alias{gp.table}
\title{A list of all the geopolitical entities in the Neotoma database.}
\format{
a \code{data.frame} object
}
\source{
The Neotoma database.
}
\usage{
gp.table
}
\description{
A list of geopolitical entities with associated numeric ID values.
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse.download}
\alias{browse.download}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
\method{browse}{download}(x)
}
\arguments{
\item{x}{A \code{download} object.}
}
\description{
Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bind.R
\name{bind}
\alias{bind}
\title{Function to bind objects together into a longer object.}
\usage{
bind(x, ...)
}
\arguments{
\item{x}{An object returned by one of the \code{get_*} commands for download, site or dataset.}

\item{...}{other objects of the same class.}
}
\value{
This command returns a larger list.
}
\description{
From multiple \code{download*}s, \code{dataset*}s or \code{site}s, join them together into a single object.
}
\details{
To support further synthesis and analysis \code{compile_download} works to transform a list
returned by \code{\link{get_download}} into a large data frame with columns for site and sample attributes
and also with the associated assemblage data at each sample depth.  This function also does the same for
single sites.
}
\examples{
\dontrun{
#  Search for sites with "Thuja" pollen that are older than 8kyr BP and
#  that are on the west coast of North America:
t8kyr.poa <- get_dataset(taxonname="Thuja*", 
                         loc=c(-150, 20, -100, 60), ageyoung = 8000)
t8kyr.canis <- get_dataset(taxonname="Canis*", 
                           loc=c(-150, 20, -100, 60), ageyoung = 8000)

t8kyr.co_site <- bind(t8kyr.poa, t8kyr.canis)
plot(t8kyr.co_site)

####
# We want to look at four different dataset types across a forest-prairie 
# boundary:
dataset_types <- c("ostracode surface sample",
                   "water chemistry",
                   "diatom surface sample",
                   "pollen surface sample")

# Run the `get_dataset` function for each of the different dataset types 
dataset_lists <- lapply(dataset_types, 
                          function(x) { 
                            get_dataset(datasettype=x, 
                                        loc = c(-100,43,-92,48))
                                        })

# Using do.call here to make sure that I don't have to split the list out.
new_datasets <- do.call(bind, dataset_lists)

# And voila!
plot(new_datasets)

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pollen.equiv.R
\docType{data}
\name{pollen.equiv}
\alias{pollen.equiv}
\title{A table to convert the pollen taxa identified by investigators to standardized lists.}
\format{
a \code{data.frame} object
}
\usage{
data(pollen.equiv)
}
\description{
A list of standardized (published) taxonomies from the literature to help standardize taxonomies for synthesis work.
}
\details{
Taxon conversion table (readable).
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}; Jeremiah Marsicek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.dataset}
\alias{get_site.dataset}
\title{Return Site Information from a numeric list of site ids.}
\usage{
\method{get_site}{dataset}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{dataset}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compile_downloads.R
\name{compile_downloads}
\alias{compile_downloads}
\title{Compile download objects}
\usage{
compile_downloads(downloads)
}
\arguments{
\item{downloads}{A download_list as returned by \code{\link{get_download}}, or multiple downloads joined in a list.}
}
\value{
This command returns a data frame.
}
\description{
Function to convert multiple downloads into a single large table.

From the assemblage data for multiple cores, return a single data.frame with columns for site
metadata and assemblage data.

To support further synthesis and analysis \code{compile_download} works to transform a list
returned by \code{\link{get_download}} into a large data frame with columns for site and sample attributes
and also with the associated assemblage data at each sample depth.  This function also does the same for
single sites.
}
\examples{
\dontrun{
#  Search for sites with "Thuja" pollen that are older than 8kyr BP and
#  that are on the west coast of North America:
t8kyr.datasets <- get_dataset(taxonname='Thuja*', 
                              loc=c(-150, 20, -100, 60), 
                              ageyoung = 8000)

#  Returns 3 records (as of 04/04/2013), get dataset for the first record, 
#  Gold Lake Bog.
thuja.sites <- get_download(t8kyr.datasets)

gold.p25 <- compile_taxa(thuja.sites, 'P25')

all.gold <- compile_downloads(gold.p25)

pollen.sums <- rowSums(all.gold[,11:ncol(all.gold)], na.rm=TRUE)

plot(x = all.gold$age, 
     y = all.gold$Cupressaceae.Taxaceae / pollen.sums, 
     col = all.gold$site.name,
     pch = 19)

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org

Gavin DG, Oswald WW, Wahl ER, Williams JW. 2003. A statistical approach to evaluating distance metrics and analog assignments for pollen records. Quaternary Research 60: 356-367.

Whitmore J, Gajewski K, Sawada M, Williams JW, Shuman B, Bartlein PJ, Minckley T, Viau AE, Webb III T, Shafer S, Anderson P, Brubaker L. 2005. Modern pollen data from North America and Greenland for multi-scale paleoenvironmental applications. Quaternary Science Reviews 24: 1828-1848.

Williams J, Shuman B. 2008. Obtaining accurate and precise environmental reconstructions from the modern analog technique and North American surface pollen dataset. Quaternary Science Reviews. 27:669-687.

API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clamodel.R, R/download.R
\name{download}
\alias{download}
\title{A class for download objects.}
\description{
A \code{download} is an object with the full record for a single dataset.

A \code{download} is an object with the full record for a single dataset.
}
\details{
TO DO

TO DO
}
\author{
Simon Goring

Simon Goring
}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication.download}
\alias{get_publication.download}
\title{A function to get publications for downloads in the Neotoma Database using the API.}
\usage{
\method{get_publication}{download}(x, ...)
}
\arguments{
\item{x}{an object of class \code{download}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
The function takes a \code{download} and returns a table with publication information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/counts.R
\name{counts}
\alias{counts}
\alias{counts.download}
\alias{counts.download_list}
\title{Access proxy count data}
\usage{
counts(obj, ...)

\method{counts}{download}(obj, ...)

\method{counts}{download_list}(obj, ...)
}
\arguments{
\item{obj}{an R object from which counts are to be extracted.}

\item{...}{arguments passed to other methods.}
}
\value{
Either a data frame of counts or a list of such objects.
}
\description{
Extract pollen or other proxy counts from data objects and returns them in a useful format.
}
\details{
Methods are available for "download" and "download_list" objects.
}
\examples{
\dontrun{
marion <- get_site('Marion Lake\%')
louise <- get_site('Louise Pond\%')
western.sites <- rbind(marion, louise)
western.data  <- get_dataset(western.sites)

western.dl <- get_download(western.data)
western.cnt <- counts(western.dl)
sapply(western.cnt, dim)
marion.cnt<- counts(western.dl[[1]])
dim(marion.cnt)
}
}
\author{
Gavin Simpson
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Stratiplot.download.R
\name{Stratiplot.download_list}
\alias{Stratiplot.download_list}
\title{Palaeoecological stratigraphic diagrams}
\usage{
\method{Stratiplot}{download_list}(x, yaxis = "age", method = "none", group = NULL, ...)
}
\arguments{
\item{x}{A \code{download_list} object.}

\item{yaxis}{One of the columns in \code{sample.meta}, including \code{depth}, \code{age}, \code{age.younger}, or \code{age.older}, default \code{age}.}

\item{method}{An option for axis transformation using \code{tran} from the \code{analogue} package.  \code{"none"} by default.}

\item{group}{An ecological group from the taxon table.}

\item{...}{variables to be passed to \code{Stratiplot}.}
}
\value{
A \code{trellis} object.
}
\description{
Draws paleoecological diagrams from a \code{download_list} object.  Allows control of variable type (using the \code{tran} function from the \code{analogue} package), and taxonomic grouping.  
This function only works for \code{download_list} objects that contain a single object.
}
\details{
A wrapper for the \code{analogue} package's \code{Stratiplot} function.  Allowing the user to plot a stratigraphic diagram directly from a \code{download} object.
}
\examples{
\dontrun{
lake_o_dl <- get_download(15925)
# This works:
Stratiplot(lake_o_dl)

lakes_o_nw <- get_download(get_site(sitename = "Lake B\%"))
# This Fails:
# Stratiplot(lake_o_nw)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.dataset_list}
\alias{get_site.dataset_list}
\title{Return Site Information from a \code{dataset_list}}
\usage{
\method{get_site}{dataset_list}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{dataset_list}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_geochron.R
\name{get_geochron}
\alias{get_geochron}
\title{Function to return geochronological data from records.}
\usage{
get_geochron(x, verbose = TRUE)
}
\arguments{
\item{x}{A numeric dataset ID or a vector of numeric dataset IDs, or an object of class of class \code{site}, \code{dataset}, \code{dataset_list}, \code{download} or \code{download_list} for which geochrons are required.}

\item{verbose}{logical; should messages on API call be printed?}
}
\value{
This command returns either an object of class \code{"try-error"}' (see \code{\link{try}}) defined by the error returned
   from the Neotoma API call, or a \code{geochronologic} object, which is a list with two components, a \code{dataset} and a geochronology table, a \code{data.frame} with the following components:

 \item{ \code{sample.id} }{A unique identifier for the geochronological unit.}
 \item{ \code{age.type} }{String.  The age type, one of calendar years, radiocarbon years, etc.}
 \item{ \code{age} }{Dated age of the material.}
 \item{ \code{e.older} }{The older error limit of the age value.  Commonly 1 standard deviation.}
 \item{ \code{e.young} }{The younger error limit of the age value.}
 \item{ \code{delta13C} }{The measured or assumed delta13C value for radiocarbon dates, if provided.}
 \item{ \code{material.dated} }{A table describing the collection, including dataset information, PI data compatible with \code{\link{get_contact}} and site data compatable with \code{\link{get_site}}.}
 \item{ \code{geo.chron.type} }{Text string, type of geochronological analysis, i.e., Radiocarbon dating, luminesence.}
 \item{ \code{notes} }{Text string}
 \item{ \code{infinite} }{Boolean, does the dated material return an "infinite" date?}

 A full data object containing all the relevant geochronological data available for a dataset.
}
\description{
Using the dataset ID, return all geochronological data associated with the dataID.  At present,
   only returns the dataset in an unparsed format, not as a data table.   This function will only download one dataset at a time.
}
\examples{
\dontrun{
#  Search for the sites around Marion Lake, BC.  I want to find sites within 
#  about 1km.

marion <- get_site(sitename = "Marion Lake*")

marion_close <- get_closest(marion, n = 10, buffer = 1)

#  Returns 116 records (as of 13/07/2015).  These are the pollen records though, 
#  we want the sites:
geochron.records <- get_geochron(marion_close)

#  We want to extract all the radiocarbon ages from the records:

get_ages <- function(x){
  any.ages <- try(x[[2]]$age[x[[2]]$age.type == 'Radiocarbon years BP'])
  if(class(any.ages) == 'try-error') output <- NA
  if(!class(any.ages) == 'try-error') output <- unlist(any.ages)
  output
}

radio.chron <- unlist(sapply(geochron.records, get_ages))

hist(radio.chron[radio.chron<40000], breaks=seq(0, 25000, by = 1000),
     main = 'Radiocarbon dates for Pseudotsuga records',
     xlab = 'Radiocarbon date (14C years before 1950)')
}

}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.download}
\alias{get_site.download}
\title{Return Site Information from a \code{download}}
\usage{
\method{get_site}{download}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{download}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compile_taxa.R
\name{compile_taxa}
\alias{compile_taxa}
\title{Function to convert assemblage taxa to standardized lists.}
\usage{
compile_taxa(object, list.name, alt.table = NULL, cf = TRUE, type = TRUE)
}
\arguments{
\item{object}{A pollen object returned by \code{\link{get_download}}.}

\item{list.name}{The taxon compilation list, one of a set of lists from the literature (e.g., \code{"P25"}, \code{"WhitmoreFull"}).  More detail in section Details.}

\item{alt.table}{A user provided table formatted with at least two columns, one called 'taxon' and the other named as in \code{list.name}.}

\item{cf}{Should taxa listed as *cf*s (*e.g.*, *cf*. *Gilia*) be considered highly resolved?}

\item{type}{Should taxa listed as types (*e.g.*, *Iva annua*-type) be considered highly resolved?}
}
\value{
This command returns a list object with the same structure as the parent pollen object returned by \code{\link{get_download}}, or a matrix (or data frame) depending on whether \code{object} is one or the other.  Any pollen taxon not included in the major taxa defined in the pollen gets returned as 'Other'.
}
\description{
From the assemblage data for the core return assemblage data with the assemblage taxa
Currently implemented only for pollen data.
}
\details{
The data object uses the smaller pollen subset.  As this package develops we will add the capacity to summarize data output from the translation. Currently we can return only subsets that have been defined in the literature.  These lists include:
\itemize{
  \item{\code{"P25"} }{ This list is derived from Gavin et al., (2003), and includes 25 pollen taxa.}
  \item{\code{"WS64"} }{  This list is derived from Williams and Shuman (2008).}
  \item{\code{"WhitmoreFull"} }{  This is the full list associated with the Whitmore et al., (2005) North American Modern Pollen Database.}
  \item{\code{"WhitmoreSmall"} }{  As above, but taxa for which both fully resolved and undifferentiated exist these taxa are summed.}
}
}
\examples{
\dontrun{
#  Search for sites with "Thuja" pollen that are older than 8kyr BP and
#  that are on the west coast of North America:
t8kyr.datasets <- get_dataset(taxonname='Thuja*', loc=c(-150, 20, -100, 60), ageyoung = 8000)

#  Returns 3 records (as of 04/04/2013), get dataset for the first record, Gold Lake Bog.
GOLDKBG <- get_download(t8kyr.datasets[[1]])

gold.p25 <- compile_taxa(GOLDKBG, 'P25')

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org

Gavin DG, Oswald WW, Wahl ER, Williams JW. 2003. A statistical approach to 
  evaluating distance metrics and analog assignments for pollen records. 
  Quaternary Research 60: 356-367.

Whitmore J, Gajewski K, Sawada M, Williams JW, Shuman B, Bartlein PJ, Minckley T, 
  Viau AE, Webb III T, Shafer S, Anderson P, Brubaker L. 2005. Modern pollen data 
  from North America and Greenland for multi-scale paleoenvironmental applications. 
  Quaternary Science Reviews 24: 1828-1848.

Williams J, Shuman B. 2008. Obtaining accurate and precise environmental 
  reconstructions from the modern analog technique and North American surface pollen 
  dataset. Quaternary Science Reviews. 27:669-687.

API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication.download_list}
\alias{get_publication.download_list}
\title{A function to get publications for datasets in the Neotoma Database using the API.}
\usage{
\method{get_publication}{download_list}(x, ...)
}
\arguments{
\item{x}{an object of class \code{download_list}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
The function takes a \code{download_list} and returns a table with publication information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol}
\alias{get_chroncontrol}
\title{Function to return chronological control tables used to build age models.}
\usage{
get_chroncontrol(x, chronology = 1, verbose = TRUE, add = FALSE)
}
\arguments{
\item{x}{A single numeric chronology ID, a vector of numeric dataset IDs as returned by \code{\link{get_dataset}} or a \code{download} or \code{download_list} object.}

\item{chronology}{When \code{download} objects have more than associated chronology, which chronology do you want?  Default is \code{1}.}

\item{verbose}{logical, should messages on API call be printed?}

\item{add}{logical, should this chron control be added to the download object?}
}
\value{
This command returns either an object of class  \code{"try-error"} containing the error returned
   from the Neotoma API call, or a full data object containing all the relevant information required to build either the default or prior chronology for a core.
   When \code{download} or \code{download_list} objects are passes, the user can \code{add} the chroncontrol to the
   \code{download} object explicitly, in which case the function will return a download with \code{chroncontrol} embedded.
   
   This is a list comprising the following items:

 \item{ \code{chron.control} }{A table describing the collection, including dataset information, PI data compatable with \code{\link{get_contact}} and site data compatable with \code{\link{get_site}}.}
 \item{ \code{meta} }{Dataset information for the core, primarily the age-depth model and chronology.  In cases where multiple age models exist for a single record the most recent chronology is provided here.}
 
 If Neotoma returns empty content, either the control table or the associated metadata (which happens in approximately 25% of cases) then the data.frames are returned with NA content.
}
\description{
Using the dataset ID, return all records associated with the data.  At present,
   only returns the dataset in an unparsed format, not as a data table.   This function will only download one dataset at a time.
}
\examples{
\dontrun{
#  The point of pulling chronology tables is to re-build or examine the 
#  chronological information that was used to build the age-depth model for 
#  the core.  You can do this by hand, but the `write_agefile` function works 
#  with `download` objects directly.

three_pines <- get_download(get_dataset(get_site("Three Pines Bog"), 
                                        datasettype = "pollen"))
pines_chron <- get_chroncontrol(three_pines)

# Spline interpolation:
model <- smooth.spline(x = pines_chron[[1]]$chron.control$depth,
                       y = pines_chron[[1]]$chron.control$age)
                       
new_ages <- predict(model, x = three_pines[[1]]$sample.meta$depth)

}
}
\references{
+ Neotoma Project Website: http://www.neotomadb.org
+ API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Stratiplot.download.R
\name{Stratiplot.download}
\alias{Stratiplot.download}
\title{Palaeoecological stratigraphic diagrams}
\usage{
\method{Stratiplot}{download}(x, yaxis = "age", method = "none", group = NULL, ...)
}
\arguments{
\item{x}{A \code{download} object.}

\item{yaxis}{One of the columns in \code{sample.meta}, including \code{depth}, \code{age}, \code{age.younger}, or \code{age.older}, default \code{age}.}

\item{method}{An option for axis transformation using \code{tran} from the \code{analogue} package.  \code{"none"} by default.}

\item{group}{An ecological group from the taxon table.}

\item{...}{variables to be passed to \code{Stratiplot}.}
}
\value{
A \code{trellis} object.
}
\description{
Draws paleoecological diagrams from a \code{download} object.  Allows control of variable type (using the \code{tran} function from the \code{analogue} package), and taxonomic grouping.
}
\details{
A wrapper for the \code{analogue} package's \code{Stratiplot} function.  Allowing the user to plot a stratigraphic diagram directly from a \code{download} object.
}
\examples{
\dontrun{
lake_o_dl <- get_download(15925)
Stratiplot(lake_o_dl[[1]])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol.dataset_list}
\alias{get_chroncontrol.dataset_list}
\title{Function to return chronological control tables from a \code{dataset_list}.}
\usage{
\method{get_chroncontrol}{dataset_list}(x, chronology = 1, verbose = TRUE, add = FALSE)
}
\arguments{
\item{x}{A \code{dataset_list} object.}

\item{chronology}{When \code{download} objects have more than associated chronology, which chronology do you want?  Default is \code{1}.}

\item{verbose}{logical; should messages on API call be printed?}

\item{add}{Should the \code{chroncontrol} be added to the download object (only accepts \code{FALSE})}
}
\description{
Using a \code{dataset_list}, return the default chron-control table.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_download.R
\name{get_download}
\alias{get_download}
\title{Function to return full download records using \code{site}s, \code{dataset}s, or dataset IDs.}
\usage{
get_download(x, verbose = TRUE)
}
\arguments{
\item{x}{A single numeric dataset ID or a vector of numeric dataset IDs as returned by \code{get_datasets}, or a \code{site}, \code{dataset}, or \code{dataset_list}.}

\item{verbose}{logical; should messages on API call be printed?}
}
\value{
This command returns either object of class \code{"try-error"}' (see \code{\link{try}}) defined by the error returned from the Neotoma API call, or an object of class \code{download_list}, containing a set of \code{download} objects, each with relevant assemblage information and metadata:
The \code{download} object is a list of lists and data frames that describe an assemblage, the constituent taxa, the chronology, site and PIs who contributed the data. The following are important components:

 \item{ \code{dataset} }{A table describing the collection, including dataset information, PI data compatible with \code{\link{get_contact}} and site data compatable with \code{\link{get_site}}.}
 \item{ \code{sample.meta} }{Dataset information for the core, primarily the age-depth model and chronology.  In cases where multiple age models exist for a single record the most recent chronology is provided here.}
 \item{ \code{taxon.list} }{The list of taxa contained within the dataset, unordered, including information that can be used in \code{\link{get_taxa}}}
 \item{ \code{counts} }{The assemblage data for the dataset, arranged with each successive depth in rows and the taxa as columns.  All taxa are described in \code{taxon.list}, the chronology is in \code{sample.data}}
 \item{ \code{lab.data} }{A data frame of laboratory data, such as exotic pollen spike, amount of sample counted, charcoal counts, etc.}
 \item{ \code{chronologies} }{A list of existing chronologies.  If only a single chronology exists for a record then this is the same as the age-model in \code{sample.meta}.}
}
\description{
Using the dataset ID, site object or dataset object, return all records associated with the data as a \code{download_list}.
}
\section{Note}{

The function returns a warning in cases where single taxa are defined by multiple taphonomic characteristics, for example grains that are identified separately as crumpled and torn in the same sample and sums these values within a sample.
In the case that a geochronology dataset is passed to \code{get_download} the function returns a message and a NULL object (that is later excised).  Use \code{get_geochron} for these objects.
The chronologies can be augmented using the function \code{get_chroncontrol}, where the individual chronology objects in \code{chronologies} will consist of a table equivalent to \code{sample.meta} and a \code{chroncontrol} object.
}

\examples{
\dontrun{
#  Search for sites with "Pseudotsuga" pollen that are older than 8kyr BP and
#  that are roughly within western British Columbia:
t8kyr.datasets <- get_dataset(taxonname='*Picea*', loc=c(-90, 41, -89, 44),
                              ageold = 20000, ageyoung=10000)

#  Returns 20 records (as of 04/04/2013), get the dataset for all records:
pollen.records <- get_download(t8kyr.datasets)

#  Standardize the taxonomies for the different records using the WS64 taxonomy.
compiled.sites <- compile_taxa(pollen.records, list.name='WS64')

#  Extract the Pseudotsuga curves for the sites:
get.curve <- function(x, taxa) {
               if (taxa \%in\% colnames(x$counts)) {
                 count <- x$counts[,taxa]/rowSums(x$counts, na.rm=TRUE)
               } else {
                 count <- rep(0, nrow(x$count))
               }
               data.frame(site = x$dataset$site.data$site.name,
               age = x$sample.meta$age,
               count = count)
             }

curves <- do.call(rbind.data.frame,
                  lapply(compiled.sites, get.curve, taxa = 'Larix/Pseudotsuga'))

#  For illustration, remove the sites with no Pseudotsuga occurance:
curves <- curves[curves$count > 0, ]

smooth.curve <- predict(loess(sqrt(count)~age, data=curves),
                        data.frame(age=seq(20000, 0, by = -100)))

plot(sqrt(count) ~ age, data = curves,
     ylab = '\% Pseudotsuga/Larix', xlab='Calibrated Years BP', pch=19,
     col=rgb(0.1, 0.1, 0.1, 0.1), xlim=c(0, 20000))
lines(seq(20000, 0, by = -100), smooth.curve, lwd=2, lty=2, col=2)

#  This figure shows us an apparent peak in Larix/Pseudotsuga pollen in the
#  early-Holocene that lends support to a warmer, drier early-Holocene in
#  western North America.
}

}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse.download_list}
\alias{browse.download_list}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
\method{browse}{download_list}(x)
}
\arguments{
\item{x}{A \code{download_list} object.}
}
\description{
Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse.default}
\alias{browse.default}
\title{Open a browser window to display a Neotoma dataset within the Neotoma Explorer}
\usage{
\method{browse}{default}(x)
}
\arguments{
\item{x}{A numeric value with the dataset ID.}
}
\description{
Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_table.R
\name{get_table}
\alias{get_table}
\title{Get Neotoma value tables.}
\usage{
get_table(table.name = NULL)
}
\arguments{
\item{table.name}{Call one of the available tables in the Neotoma Database.
A full listing of tables can be found here: \url{http://wnapi.neotomadb.org/doc/resources/dbtables}.
By default it returns all objects in the table.}
}
\description{
Get Neotoma value tables.
}
\details{
A table of values corresponding to the parameter of interest.
}
\examples{
\dontrun{
taxon.table <- get_table('Taxa')

#  Get the frequency of a random taxon in Neotoma.
tax_sample <- sample(nrow(taxon.table), 1)
cat("The taxon", 
    taxon.table$TaxonName[tax_sample], 
    "occurs in Neotoma", 
    length(get_dataset(taxonname = taxon.table$TaxonName[tax_sample])), 
    "times.")
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa.R
\name{taxa}
\alias{taxa}
\alias{taxa.download}
\alias{taxa.download_list}
\title{Access proxy taxonomic data}
\usage{
taxa(obj, ...)

\method{taxa}{download}(obj, ...)

\method{taxa}{download_list}(obj, collapse = TRUE, hierarchy = FALSE, ...)
}
\arguments{
\item{obj}{an R object from which counts are to be extracted.}

\item{...}{arguments passed to other methods.}

\item{collapse}{should the results be returned as a list, one for each site (\code{FALSE}), or a single dataframe of all taxa? Default is \code{TRUE}}

\item{hierarchy}{Should the taxonomic hierarchy be included?}
}
\value{
Either a data frame of taxa or a list of such objects.
}
\description{
Extracts taxa from \code{download} objects and returns them in a useful format.
}
\details{
Methods are available for "download" and "download_list" objects.
}
\examples{
\dontrun{
ostracodes <- get_dataset(datasettype = 'ostracode')

ostro.dl <- get_download(ostracodes)
ostro.taxa <- taxa(ostro.dl)
}
}
\author{
Simon Goring
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.download}
\alias{get_dataset.download}
\title{Obtain dataset information from an existing \code{download} object.}
\usage{
\method{get_dataset}{download}(x, ...)
}
\arguments{
\item{x}{An object of class \code{download}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access a \code{dataset} within a \code{download} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.integer}
\alias{get_site.integer}
\title{Return Site Information from a vector of integers.}
\usage{
\method{get_site}{integer}(sitename, ...)
}
\arguments{
\item{sitename}{An integer or vector of integers.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ages.R
\name{ages}
\alias{ages}
\alias{ages.download}
\alias{ages.download_list}
\title{Access proxy age data}
\usage{
ages(obj, ...)

\method{ages}{download}(obj, ...)

\method{ages}{download_list}(obj, ...)
}
\arguments{
\item{obj}{an R object from which counts are to be extracted.}

\item{...}{arguments passed to other methods.}
}
\value{
Either a data frame of ages or a list of such objects.
}
\description{
Extracts age information from objects and returns them in a useful format.
}
\details{
Methods are available for "download" and "download_list" objects.
}
\examples{
\dontrun{
ostracodes <- get_dataset(datasettype = 'ostracode')

ostro.dl <- get_download(ostracodes)
ostro.ages <- ages(ostro.dl)
}
}
\author{
Simon Goring
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.download_list}
\alias{get_site.download_list}
\title{Return Site Information from a \code{download_list}}
\usage{
\method{get_site}{download_list}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{download_list}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_taxa.R
\name{get_taxa}
\alias{get_taxa}
\title{Get taxon information from Neotoma.}
\usage{
get_taxa(taxonid, taxonname, status, taxagroup, ecolgroup)
}
\arguments{
\item{taxonid}{Numeric taxon identifier used in Neotoma}

\item{taxonname}{A character string representing the full or partial name of taxa of interest.}

\item{status}{The current status of the taxon, one of 'extinct', 'extant', 'all'.}

\item{taxagroup}{The taxonomic grouping for the taxa. See \url{http://wnapi.neotomadb.org/doc/resources/taxa} for the list of approved groupings.}

\item{ecolgroup}{The ecological group of the taxa. More detailed than \code{taxagroup}, can be obtained using \code{get_table("EcolGroupTypes")}.}
}
\value{
Returns a data frame with the following components:

 \item{ \code{TaxonID} }{Unique database record identifier for a taxon}
 \item{ \code{TaxonCode} }{Shorthand notation for a taxon identification}
 \item{ \code{TaxonName} }{Name of the taxon}
 \item{ \code{Author} }{Author(s) of the name. Used almost exclusively with beetle taxa}
 \item{ \code{Extinct} }{True if extinct; false if extant}
 \item{ \code{TaxaGroup} }{Code for taxa group to which taxon belongs}
 \item{ \code{EcolGroups} }{Array of ecological group codes to which the taxon belongs}
 \item{ \code{HigherTaxonID} }{TaxonID of the next higher taxonomic rank}
 \item{ \code{PublicationID} }{Publication identification number}
 \item{ \code{Notes} }{Free-form notes or comments about the taxon}
}
\description{
Get taxon information from Neotoma.
}
\examples{
\dontrun{
## Return all species taxa with "Abies" in name - note wildcard
taxa <- get_taxa(taxonname = "Abies*")
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_download.R
\name{get_download.site}
\alias{get_download.site}
\title{Function to return full download records using a \code{site}.}
\usage{
\method{get_download}{site}(x, verbose = TRUE)
}
\arguments{
\item{x}{An object of class \code{site}.}

\item{verbose}{logical; should messages on API call be printed?}
}
\description{
Using a \code{site}, return all records associated with the data as a \code{download_list}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_download.R
\name{get_download.dataset_list}
\alias{get_download.dataset_list}
\title{Function to return full download records using a \code{dataset_list}.}
\usage{
\method{get_download}{dataset_list}(x, verbose = TRUE)
}
\arguments{
\item{x}{An object of class \code{dataset_list}.}

\item{verbose}{logical; should messages on API call be printed?}
}
\description{
Using a \code{dataset_list}, return all records associated with the data as a \code{download_list}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication}
\alias{get_publication}
\title{A function to get publications for sites or datasets in the Neotoma Database using the API.}
\usage{
get_publication(x, contactid, datasetid, author, pubtype, year, search)
}
\arguments{
\item{x}{Numeric Publication ID value, either from \code{\link{get_dataset}} or known.}

\item{contactid}{Numeric Contact ID value, either from \code{\link{get_dataset}} or \code{\link{get_contact}}}

\item{datasetid}{Numeric Dataset ID, known or from \code{\link{get_dataset}}}

\item{author}{Character string for full or partial author's name.  Can include wildcards such as 'Smit*' for all names beginning with 'Smit'.}

\item{pubtype}{Character string, one of eleven allowable types, see \code{\link{get_table}}. For a list of allowed types run \code{get_table("PublicationTypes")}.}

\item{year}{Numeric publication year.}

\item{search}{A character string to search for within the article citation.}
}
\value{
A list is returned with two data frame components:

 \item{ \code{meta} }{A single row with Publication ID, type, year of publication and full citation.}
 \item{ \code{Authors} }{\code{data.frame} of author names, order and IDs, can be of variable length.}
}
\description{
The function takes the parameters, defined by the user, and returns a table with publication information from the Neotoma Paleoecological Database.
}
\examples{
\dontrun{
#  To find all publications from 1998:
year.cont <- get_publication(year = 1998)

# To find all data contributors who have the last name "Smith"
smith.cont <- get_publication(author = 'Smith')
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_download.R
\name{get_download.dataset}
\alias{get_download.dataset}
\title{Function to return full download records using a \code{dataset}.}
\usage{
\method{get_download}{dataset}(x, verbose = TRUE)
}
\arguments{
\item{x}{An object of class \code{dataset}.}

\item{verbose}{logical; should messages on API call be printed?}
}
\description{
Using a \code{dataset}, return all records associated with the data as a \code{download_list}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon.table.R
\docType{data}
\name{taxon.list}
\alias{taxon.list}
\title{Neotoma taxon list}
\format{
a \code{data.frame} object
}
\source{
The Neotoma database.
}
\usage{
data(taxon.list)
}
\description{
The taxonomy table for datasets in neotoma, as would be returned by \code{\link{get_table}}
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.download_list}
\alias{get_dataset.download_list}
\title{Obtain dataset information from a \code{download_list}.}
\usage{
\method{get_dataset}{download_list}(x, ...)
}
\arguments{
\item{x}{An object of class \code{download_list}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to return datasets corresponding to the objects within a \code{download_list}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chroncontrol.R
\name{get_chroncontrol.dataset}
\alias{get_chroncontrol.dataset}
\title{Function to return chronological control tables from a \code{dataset}.}
\usage{
\method{get_chroncontrol}{dataset}(x, chronology = 1, verbose = TRUE, add = FALSE)
}
\arguments{
\item{x}{A \code{dataset}.}

\item{chronology}{When \code{download} objects have more than associated chronology, which chronology do you want?  Default is \code{1}.}

\item{verbose}{logical; should messages on API call be printed?}

\item{add}{Should the \code{chroncontrol} be added to the download object (only accepts \code{FALSE})}
}
\description{
Using a \code{dataset}, return the default chron-control table.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/depths.R
\name{depths}
\alias{depths}
\alias{depths.default}
\alias{depths.download}
\alias{depths.download_list}
\title{Extracts the depth values from a `download` object}
\usage{
depths(obj, ...)

\method{depths}{default}(obj, ...)

\method{depths}{download}(obj, ...)

\method{depths}{download_list}(obj, ...)
}
\arguments{
\item{obj}{A \code{download} object.}

\item{...}{arguments passed to other methods.}
}
\value{
Returns a vector of depths.
}
\description{
Using a \code{download} object, return the sample depths (if available).

Using a numeric value, \code{download}, \code{download_list}, \code{dataset} or \code{dataset_list} object, open up a browser window in the users default browser. Very large objects
}
\examples{
\dontrun{
# Provide a vector of depths to generate a new age model:
# The dataset id 684 is for Devils Lake, a record published by Louis Maher Jr.

pollen.data <- get_download(684)
pollen.chron <- get_chroncontrol(pollen.data)[[1]]

age_sds <- pollen.chron$chron.control$age - focal$chron.control$age.young,
get_curves <- ifelse(regexpr("Radiocarbon",
                             pollen.chron$chron.control$control.type) > -1, 
                     'intcal13', 'normal')

new_chron <- Bchron::Bchronology(ages   = pollen.chron$chron.control$age,
                                 ageSds = age_sds
                                 positions = pollen.chron$chron.control$depth,
                                 calCurves = , 
                                 predictPositions = depths(pollen.data))

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/sites
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset}
\alias{get_dataset}
\title{Obtain dataset information from the Neotoma Paleoecological Database or an existing object.}
\usage{
get_dataset(
  x,
  datasettype,
  piid,
  altmin,
  altmax,
  loc,
  gpid,
  taxonids,
  taxonname,
  ageold,
  ageyoung,
  ageof,
  subdate
)
}
\arguments{
\item{x}{An optional value, either a \code{numeric} site ID or object of class \code{download}, \code{download_list} or \code{site}.}

\item{datasettype}{A character string corresponding to one of the allowed dataset types in the Neotoma Database.  Allowed types include: \code{"geochronologic"}, \code{"loss-on-ignition"}, \code{"pollen"}, \code{"plant macrofossils"}, \code{"vertebrate fauna"}, \code{"mollusks"}, and \code{"pollen surface sample"}.  See note in Details delow.}

\item{piid}{Numeric value for the Principle Investigator's ID number.}

\item{altmin}{Numeric value indicating the minimum altitude for the site (can be used alone or with \code{altmax}).}

\item{altmax}{Numeric value indicating the maximum altitude for the site (can be used alone or with \code{altmin}).}

\item{loc}{A numeric vector \code{c(lonW, latS, lonE, latN)} representing the bounding box within which to search for sites.  The convention here is to use negative values for longitudes west of Greenwich or longitudes south of the equator}

\item{gpid}{A character string or numeric value, must correspond to a valid geopolitical identity in the Neotoma Database.  Use get.tables('GeoPoliticalUnits') for a list of acceptable values, or link here: \url{http://wnapi.neotomadb.org/apdx/geopol.htm}}

\item{taxonids}{A numeric identifier for the taxon.  See \code{\link{get_table}} and use \code{get_table('Taxa')} for a list of acceptable values.}

\item{taxonname}{A character string corresponding to a valid taxon identity in the Neotoma Database.  See \code{\link{get_table}} and use \code{get_table('Taxa')} for a list of acceptable values.}

\item{ageold}{The oldest date acceptable for the search (in years before present).}

\item{ageyoung}{The youngest date acceptable for the search.}

\item{ageof}{If a taxon ID or taxon name is defined this parameter must be set to \code{"taxon"}, otherwise it may refer to \code{"sample"}, in which case the age bounds are for any samples within datasets or \code{"dataset"} if you want only datasets that are within the bounds of ageold and ageyoung.}

\item{subdate}{Date of dataset submission, either YYYY-MM-DD or MM-DD-YYYY.}
}
\value{
More details on the use of these parameters can be obtained from
   \url{http://wnapi.neotomadb.org/doc/resources/datasets}.

   A list of class `dataset_list`, with each item corresponding to an individual record.
   Searches that return no items will result in a NULL value being returned.
   Otherwise each list item (each dataset record) includes the following components:

 \item{ \code{dataset.id} }{Unique database record identifier for the dataset.}
 \item{ \code{dataset.name}  }{Name of the dataset; not commonly used.}
 \item{ \code{CollUnitHandle}  }{Code name of the Collection Unit with which the dataset is associated. This code may be up to 10 characters. Data are frequently distributed by Collection Unit, and the Handle is used for file names.}
 \item{ \code{CollUnitID}  }{Unique database record identifier for the collection unit.}
 \item{ \code{CollType}  }{The collection type. Types include cores, sections, excavations, and animal middens.}
 \item{ \code{DatasetType}  }{The dataset type, such as: geochronologic, loss-on-ignition, pollen, plant macrofossils, vertebrate fauna, etc.}
 \item{ \code{AgeOldest}  }{The oldest of all sample ages (in calendar years before present) in the dataset.}
 \item{ \code{AgeYoungest}  }{The youngest of all sample ages (in calendar years before present) in the dataset.}
 \item{ \code{SubDates}  }{An array of objects that describe dataset submission events.  If multiple submissions occurred then this is a table.}
 \item{ \code{DatasetPIs}  }{An array of objects that describe Principal Investigators associated with a dataset.}
 \item{ \code{Site}  }{An object describing the site where the dataset samples were taken.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
\details{
With regards to \code{datasettypes}, because Neotoma is a "living" database, and new dataset types are being added in an ongoing manner as new research disciplines use the database, you can use \code{get_table("datasettypes")} to see the full list of available dataset types in the database.
}
\examples{
\dontrun{
# Search for sites with "Thuja" pollen that are older than 8kyr BP and
# that are on the west coast of North America:
t8kyr.datasets <- get_dataset(taxonname='Thuja*', 
                              loc=c(-150, 20, -100, 60), 
                              ageyoung = 8000)

# Search for vertebrate fossils in Canada (gpid: 756) within the last 2kyr.
gpids <- get_table(table.name='GeoPoliticalUnits')
canID <- gpids[which(gpids$GeoPoliticalName == 'Canada'),1]

v2kyr.datasets <- get_dataset(datasettype='vertebrate fauna', 
                              gpid=canID, 
                              ageold = 2000)
}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/contacts
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site}
\alias{get_site}
\title{Return Site Information.}
\usage{
get_site(sitename, altmin, altmax, loc, gpid, ...)
}
\arguments{
\item{sitename}{character string representing the full or partial site name, or an object of class \code{dataset}, \code{dataset_list}, \code{download} or \code{download_list}}

\item{altmin}{Minimum site altitude  (in m).}

\item{altmax}{Maximum site altitude (in m).}

\item{loc}{A numeric vector c(lonW, latS, lonE, latN) representing the bounding box within which to search for sites.  The convention here is to use negative values for longitudes west of Grewnwich or longitudes south of the equator.}

\item{gpid}{A character string or numeric value, must correspond to a valid geopolitical identity in the Neotoma Database.  Use get.tables('GeoPoliticalUnits') for a list of acceptable values, or link here: http://wnapi.neotomadb.org/apdx/geopol.htm}

\item{...}{Optional additional arguments}
}
\value{
A data frame:

 \item{\code{siteid}}{Unique database record identifier for the site.}
 \item{\code{sitename}}{Name of the site.}
 \item{\code{long}}{Mean longitude, in decimal degrees, for a site (-180 to 180).}
 \item{\code{lat}}{Mean latitude, in decimal degrees, for a site (-90 to 90).}
 \item{\code{elev}}{Elevation in meters.}
 \item{\code{description}}{Free form description of a site, including such information as physiography and vegetation around the site.}
 \item{\code{long_acc}}{If the site is described by a bounding box this is the box width.}
 \item{\code{lat_acc}}{If the site is described by a bounding box this is the box height.}
}
\description{
Return site information from the Neotoma Paleoecological Database.

\code{get_site} returns site information from the Neotoma Paleoecological Database
   based on parameters defined by the user.
}
\examples{
\dontrun{
#  What is the distribution of site elevations in Neotoma?
all.sites <- get_site()  #takes a bit of time.

plot(density(all.sites$elev, from = 0, na.rm=TRUE),
main = 'Altitudinal Distribution of Neotoma Sites', xlab = 'Altitude (m)', log='x')

#  Get site information from a dataset:
nw.datasets <- get_dataset(loc = c(-140, 50, -110, 65), 
                           datasettype='pollen',
                           taxonname='Pinus*')
                           
nw.sites <- get_site(nw.datasets)

}
}
\references{
Neotoma Project Website: http://www.neotomadb.org
API Reference:  http://wnapi.neotomadb.org/doc/resources/sites
}
\author{
Simon J. Goring \email{simon.j.goring@gmail.com}
}
\keyword{IO}
\keyword{connection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.geochronologic_list}
\alias{get_site.geochronologic_list}
\title{Return Site Information from a \code{geochronologic_list}}
\usage{
\method{get_site}{geochronologic_list}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{geochronologic_list}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.geochronologic_list}
\alias{get_dataset.geochronologic_list}
\title{Obtain dataset information from an object of class \code{geochronologic_list}.}
\usage{
\method{get_dataset}{geochronologic_list}(x, ...)
}
\arguments{
\item{x}{An object of class \code{geochronologic_list}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site.geochronologic}
\alias{get_site.geochronologic}
\title{Return Site Information from a \code{geochronologic}}
\usage{
\method{get_site}{geochronologic}(sitename, ...)
}
\arguments{
\item{sitename}{An object of class \code{geochronologic}.}

\item{...}{Arguments passed from the generic method, not used.}
}
\description{
Return site information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.geochronologic}
\alias{get_dataset.geochronologic}
\title{Obtain dataset information from an object of class \code{geochronologic}.}
\usage{
\method{get_dataset}{geochronologic}(x, ...)
}
\arguments{
\item{x}{An object of class \code{geochronologic}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication.dataset_list}
\alias{get_publication.dataset_list}
\title{A function to get publications for dataset_lists in the Neotoma Database using the API.}
\usage{
\method{get_publication}{dataset_list}(x, ...)
}
\arguments{
\item{x}{an object of class \code{dataset_list}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
The function takes a \code{dataset_list} and returns a table with publication information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset.default}
\alias{get_dataset.default}
\title{Obtain dataset information from the Neotoma Paleoecological Database or an existing object.}
\usage{
\method{get_dataset}{default}(
  x,
  datasettype,
  piid,
  altmin,
  altmax,
  loc,
  gpid,
  taxonids,
  taxonname,
  ageold,
  ageyoung,
  ageof,
  subdate
)
}
\arguments{
\item{x}{A numeric value corresponding to the site ID.}

\item{datasettype}{A character string corresponding to one of the allowed dataset types in the Neotoma Database.  You can find the full list of allowed datasettypes using: \code{get_table("datasettypes")}.}

\item{piid}{Numeric value for the Principle Investigator's ID number.}

\item{altmin}{Numeric value indicating the minimum altitude for the site (can be used alone or with \code{altmax}).}

\item{altmax}{Numeric value indicating the maximum altitude for the site (can be used alone or with \code{altmin}).}

\item{loc}{A numeric vector \code{c(lonW, latS, lonE, latN)} representing the bounding box within which to search for sites.  The convention here is to use negative values for longitudes west of Greenwich or longitudes south of the equator}

\item{gpid}{A character string or numeric value, must correspond to a valid geopolitical identity in the Neotoma Database.  Use get.tables('GeoPoliticalUnits') for a list of acceptable values, or link here: \url{http://wnapi.neotomadb.org/apdx/geopol.htm}}

\item{taxonids}{A numeric identifier for the taxon.  See \code{\link{get_table}} and use \code{get_table('Taxa')} for a list of acceptable values.}

\item{taxonname}{A character string corresponding to a valid taxon identity in the Neotoma Database.  See \code{\link{get_table}} and use \code{get_table('Taxa')} for a list of acceptable values.}

\item{ageold}{The oldest date acceptable for the search (in years before present).}

\item{ageyoung}{The youngest date acceptable for the search.}

\item{ageof}{If a taxon ID or taxon name is defined this parameter must be set to \code{"taxon"}, otherwise it may refer to \code{"sample"}, in which case the age bounds are for any samples within datasets or \code{"dataset"} if you want only datasets that are within the bounds of ageold and ageyoung.}

\item{subdate}{Date of dataset submission, either YYYY-MM-DD or MM-DD-YYYY.}
}
\description{
A function to access the Neotoma API and return datasets corresponding to the parameters defined by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication.dataset}
\alias{get_publication.dataset}
\title{A function to get publications for datasets in the Neotoma Database using the API.}
\usage{
\method{get_publication}{dataset}(x, ...)
}
\arguments{
\item{x}{an object of class \code{dataset}.}

\item{...}{objects passed from the generic.  Not used in the call.}
}
\description{
The function takes a \code{dataset} and returns a table with publication information from the Neotoma Paleoecological Database.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_publication.R
\name{get_publication.default}
\alias{get_publication.default}
\title{A function to get publications for sites or datasets in the Neotoma Database using the API.}
\usage{
\method{get_publication}{default}(x, contactid, datasetid, author, pubtype, year, search)
}
\arguments{
\item{x}{Numeric Publication ID value, either from \code{\link{get_dataset}} or known.}

\item{contactid}{Numeric Contact ID value, either from \code{\link{get_dataset}} or \code{\link{get_contact}}}

\item{datasetid}{Numeric Dataset ID, known or from \code{\link{get_dataset}}}

\item{author}{Character string for full or partial author's name.  Can include wildcards such as 'Smit*' for all names beginning with 'Smit'.}

\item{pubtype}{Character string, one of eleven allowable types, see \code{\link{get_table}}. For a list of allowed types run \code{get_table("PublicationTypes")}.}

\item{year}{Numeric publication year.}

\item{search}{A character string to search for within the article citation.}
}
\description{
The function takes the parameters, defined by the user, and returns a table with
   publication information from the Neotoma Paleoecological Database.
}
