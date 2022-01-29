
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
