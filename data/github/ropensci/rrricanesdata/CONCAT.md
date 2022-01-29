# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and this project adheres to [Semantic Versioning](http://semver.org/).

## [0.1.7] 2018-12-01

### Added
  - EP, AL cyclones for Nov, 2018

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

## [0.1.3] 2018-07-04

### Added
  - EP, AL cyclones for Jul, 2018

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

## [0.1.2] 2018-07-01

### Added
  - EP, AL cyclones for June, 2018

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

## [0.1.1] 2018-06-01

### Added
  - EP, AL cyclones for May, 2018

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

## [0.1] 2018-01-01

### Added
  - Includes all storm data through the 2017 season.

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

## [0.0.1.6] 2017-12-09

### Added
  - Added Oct, Nov cylones

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

## [0.0.1.5] 2017-10-05

### Added
  - Added September cylones

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - `fcst` now has 120 hour data. (See [rrricanes #107](https://github.com/ropensci/rrricanes/issues/107))
  - Updated `storms` (#2)

## [0.0.1.4] 2017-09-01

### Added
  - Updated storm data for month of August

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

## [0.0.1.3] 2017-08-05

### Added
  - `./data-raw/update.R` to update existing datasets.

### Changed
  - All relevant datasets updated through 2017-07-31.

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

## [0.0.1.1] 2017-07-18

### Added
  - NA

### Changed
  - Correct documentation

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

## [0.0.1] 2017-07-16

### Added
  - All datasets added, both basins, 1998 to Jun 13, 2017.

### Changed
  - NA

### Removed
  - NA

### Deprecated
  - NA

### Fixed
  - NA

[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.3.3-6666ff.svg)](https://cran.r-project.org/) [![GitHub (pre-)release](https://img.shields.io/github/release/ropensci/rrricanesdata/all.svg)](https://github.com/ropensci/rrricanesdata/tags) [![](https://badges.ropensci.org/118_status.svg)](https://github.com/ropensci/onboarding/issues/118) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rrricanesdata)](https://cran.r-project.org/package=rrricanesdata) [![Build Status](https://img.shields.io/travis/ropensci/rrricanesdata/master.svg)](https://travis-ci.org/ropensci/rrricanesdata) [![AppVeyor Build Status](https://img.shields.io/appveyor/ci/timtrice/rrricanesdata-on6xt/master.svg)](https://ci.appveyor.com/project/timtrice/rrricanesdata-on6xt)

rrricanesdata
=============

`rrricanesdata` is the complementary dataset for [rrricanes](https://github.com/ropensci/rrricanes). Currently it contains all product datasets for every cyclone that has developed in the north Atlantic and northeast Pacific basins since 1998.

`rrricanesdata` will be updated on the first of every month.

Prerequisites
-------------

No prerequisites.

Installing
----------

`rrricanesdata` can be installed with the following command:

``` r
install.packages("rrricanesdata", 
                 repos = "https://timtrice.github.io/drat/", 
                 type = "source")
```

This package is not nor will be made available in CRAN.

### Updating

As `rrricanesdata` will be updated monthly (provided any advisory was issued the previous month), you will want to ensure you have the latest datasets prior to running any analysis. You need to add the package repository to your `repos` options, such as:

``` r
options(repos = c(getOption("repos"), "https://timtrice.github.io/drat/")
```

You can then run `update.packages` to check for and install any new updates.

Datasets
--------

### adv

-   Key: Unique identifier of cyclone
-   Adv: Advisory number
-   Date: Date and time of advisory
-   Status: Classification of cyclone
-   Name: Name of cyclone
-   Lat: Latitude of cyclone center
-   Lon: Longitude of cyclone center
-   Wind: Maximum sustained one-minute winds in knots
-   Gust: Maximum sustained one-minute gusts in knots
-   Pressure: Minimum central pressure in millibars
-   PosAcc: Position accuracy of cyclone in nautical miles
-   FwdDir: Compass angle of forward motion
-   FwdSpeed: Forward speed in miles per hour
-   Eye: Size of eye in nautical miles
-   SeasNE: Radius of 12ft seas in northeast quadrant
-   SeasSE: Radius of 12ft seas in southeast quadrant
-   SeasSW: Radius of 12ft seas in southwest quadrant
-   SeasNW: Radius of 12ft seas in northwest quadrant

### discus

-   Status: Classification of storm, e.g., Tropical Storm, Hurricane, etc.
-   Name: Name of storm
-   Adv: Advisory Number
-   Date: Date of advisory issuance
-   Key: ID of cyclone
-   Contents: Text content of product

### fcst

-   Key: Unique identifier of cyclone
-   Adv: Advisory number
-   Date: Date and time of advisory
-   FcstDate: Forecast date and time in UTC
-   Lat: Forecast latitude
-   Lon: Forecast Longitude
-   Wind: Forecast wind in knots
-   Gust: Forecast gust in knots

### fcst\_wr

-   Key: Unique identifier of cyclone
-   Adv: Advisory number
-   Date: Date and time of advisory
-   FcstDate: Forecast date and time in UTC
-   WindField: Minimum sustained wind field for quadrants
-   NE: Radius in nautical miles for northeast quadrant
-   SE: Radius in nautical miles for southeast quadrant
-   SW: Radius in nautical miles for southwest quadrant
-   NW: Radius in nautical miles for northwest quadrant

### fstadv

-   Status: Classification of cyclone
-   Name: Name of cyclone
-   Adv: Advisory number
-   Date: Date and time of advisory
-   Key: Unique identifier of cyclone
-   Lat: Latitude of cyclone center
-   Lon: Longitude of cyclone center
-   Wind: Maximum sustained one-minute winds in knots
-   Gust: Maximum sustained one-minute gusts in knots
-   Pressure: Minimum central pressure in millibars
-   PosAcc: Position accuracy of cyclone in nautical miles
-   FwdDir: Compass angle of forward motion
-   FwdSpeed: Forward speed in miles per hour
-   Eye: Size of eye in nautical miles
-   NE64: Radius of &gt;=64kt winds in northeast quadrant
-   SE64: Radius of &gt;=64kt winds in southeast quadrant
-   SW64: Radius of &gt;=64kt winds in southwest quadrant
-   NW64: Radius of &gt;=64kt winds in northwest quadrant
-   NE50: Radius of &gt;=50kt winds in northeast quadrant
-   SE50: Radius of &gt;=50kt winds in southeast quadrant
-   SW50: Radius of &gt;=50kt winds in southwest quadrant
-   NW50: Radius of &gt;=50kt winds in northwest quadrant
-   NE34: Radius of &gt;=34kt winds in northwest quadrant
-   SE34: Radius of &gt;=34kt winds in southeast quadrant
-   SW34: Radius of &gt;=34kt winds in southwest quadrant
-   NW34: Radius of &gt;=34kt winds in northwest quadrant
-   Hr{n}FcstDate: Forecast valid date
-   Hr{n}Lat: Forecast latitude in n hours
-   Hr{n}Lon: Forecast longitude in n hours
-   Hr{n}Wind: Forecast maximum wind in n hours
-   Hr{n}Gust: Forecast maximum gust in n hours
-   Hr{n}NE64: Forecast wind radius in n hours
-   Hr{n}SE64: Forecast wind radius in n hours
-   Hr{n}SW64: Forecast wind radius in n hours
-   Hr{n}NW64: Forecast wind radius in n hours
-   Hr{n}NE50: Forecast wind radius in n hours
-   Hr{n}SE50: Forecast wind radius in n hours
-   Hr{n}SW50: Forecast wind radius in n hours
-   Hr{n}NW50: Forecast wind radius in n hours
-   Hr{n}NE34: Forecast wind radius in n hours
-   Hr{n}SE34: Forecast wind radius in n hours
-   Hr{n}SW34: Forecast wind radius in n hours
-   Hr{n}NW34: Forecast wind radius in n hours
-   SeasNE: Radius of 12ft seas in northeast quadrant
-   SeasSE: Radius of 12ft seas in southeast quadrant
-   SeasSW: Radius of 12ft seas in southwest quadrant
-   SeasNW: Radius of 12ft seas in northwest quadrant

### posest

-   Status: Classification of storm, e.g., Tropical Storm, Hurricane, etc.
-   Name: Name of storm
-   Date: Date of advisory issuance
-   Contents: Text content of product

### prblty

-   Status: Classification of storm, e.g., Tropical Storm, Hurricane, etc.
-   Name: Name of storm
-   Adv: Advisory Number
-   Date: Date of advisory issuance
-   Location: Location for which the probability statistics rely
-   A: Probability of a strike within the next 12 hours
-   B: Probability of a strike between 12 and 24 hours
-   C: Probability of a strike between 24 and 36 hours
-   D: Probability of a strike between 36 and 48 hours
-   E: Probability of a strike between 48 and 72 hours

### public

-   Status: Classification of storm, e.g., Tropical Storm, Hurricane, etc.
-   Name: Name of storm
-   Adv: Advisory Number
-   Date: Date of advisory issuance
-   Key: Unique ID of the cyclone
-   Contents: Text content of product

### storms

-   Key: Storm ID
-   Name: Storm name
-   Wind: Peak wind speed in knots
-   StartDate: Date/time of first advisory
-   EndDate: Date/time of last advisory

### update

-   Status: Classification of storm, e.g., Tropical Storm, Hurricane, etc.
-   Name: Name of storm
-   Date: Date of advisory issuance
-   Key: Unique ID of cyclone
-   Contents: Text content of product

### wndprb

-   Status: Classification of storm, e.g., Tropical Storm, Hurricane, etc.
-   Name: Name of storm
-   Adv: Advisory Number
-   Date: Date of advisory issuance
-   Wind: Minimum wind speed for which probabilities reference
-   Wind12: Probability of sustained Wind within 12 hours
-   Wind24: Probability of sustained Wind within 24 hours
-   Wind24Cum: Cumulative probability through 24 hours
-   Wind36: Probability of sustained Wind within 36 hours
-   Wind36Cum: Cumulative probability through 36 hours
-   Wind48: Probability of sustained Wind within 48 hours
-   Wind48Cum: Cumulative probability through 48 hours
-   Wind72: Probability of sustained Wind within 72 hours
-   Wind72Cum: Cumulative probability through 72 hours
-   Wind96: Probability of sustained Wind within 96 hours
-   Wind96Cum: Cumulative probability through 96 hours
-   Wind120: Probability of sustained Wind within 120 hours
-   Wind120Cum: Cumulative probability through 120 hours

### wr

-   Key: Unique identifier of cyclone
-   Adv: Advisory number
-   Date: Date and time of advisory
-   Windfield: Minimum wind speed expected
-   NE: Radius of Windfield in the northeast quadrant
-   SE: Radius of Windfield in the southeast quadrant
-   SW: Radius of Windfield in the southwest quadrant
-   NW: Radius of Windfield in the northwest quadrant

Built With
----------

-   [R 3.3.3](https://www.r-project.org/) - The R Project for Statistical Computing

Contributing
------------

Please read [CONTRIBUTING.md](https://github.com/ropensci/rrricanesdata/blob/master/.github/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

Versioning
----------

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/ropensci/rrricanesdata/tags).

Authors
-------

-   **Tim Trice** - *Initial work* - [timtrice](https://github.com/timtrice)

See also the list of [contributors](https://github.com/ropensci/rrricanesdata/contributors) who participated in this project.

License
-------

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

Acknowledgments
---------------

-   [Salmon, Maëlle](https://github.com/maelle)
-   [Stachelek, Joseph](https://github.com/jsta)
-   [Robinson, Emily](https://github.com/robinsones)

Known Data Quality Issues
-------------------------

1.  Hurricane Juan (AL152003), Adv 15; no status leads to improper `Status` and `Name` values in some datasets. ([\#82](https://github.com/ropensci/rrricanes/issues/82))
rrricanesdata 0.2.1 (2019-07-01)
==================================

### NEW FEATURES

* Add cyclone data for Jun, 2019

rrricanesdata 0.2.0 (2019-06-17)
==================================

### NEW FEATURES

* Add AL012019 (May, 2019)

rrricanesdata 0.1.7 (2018-12-01)
==================================

### NEW FEATURES

* Added cyclones for Nov., 2018

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.1.6 (2018-11-01)
==================================

### NEW FEATURES

* Added cyclones for Oct., 2018

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.1.5 (2018-10-01)
==================================

### NEW FEATURES

* Added cyclones for Sep., 2018

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.1.4 (2018-09-05)
==================================

### NEW FEATURES

* Added cyclones for Aug., 2018

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.1.3 (2018-08-04)
==================================

### NEW FEATURES

* Added cyclones for Jul, 2018

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.1.2 (2018-07-01)
==================================

### NEW FEATURES

* Added cyclones for June, 2018

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.1.1 (2018-06-01)
==================================

### NEW FEATURES

* Added cyclones for May, 2018

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.1 (2018-01-01)
==================================

### NEW FEATURES

* Includes all data through 2017 season.

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.0.1.6 (2017-12-09)
==================================

### NEW FEATURES

* Datasets updated through 2017-12-08

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.0.1.5 (2017-10-05)
==================================

### NEW FEATURES

* Datasets updated through 2017-09-30

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.0.1.3 (2017-08-05)
==================================

### NEW FEATURES

* Datasets updated through 2017-07-31

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA

rrricanesdata 0.0.1 (2017-07-16)
==================================

### NEW FEATURES

* All product and tidied datasets, both basins, added; covers 1998 to July 13, 2017.

### MINOR IMPROVEMENTS

* NA

### BUG FIXES

* NA

### DEPRECATED AND DEFUNCT

* NA
# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change.

Please note we have a code of conduct, please follow it in all your interactions with the project.

## Pull Request Process

1. Ensure any install or build dependencies are removed before the end of the layer when doing a
   build.
2. Update the README.md with details of changes to the interface, this includes new environment
   variables, exposed ports, useful file locations and container parameters.
3. Increase the version numbers in any examples files and the README.md to the new version that this
   Pull Request would represent. The versioning scheme we use is [SemVer](http://semver.org/).
4. You may merge the Pull Request in once you have the sign-off of two other developers, or if you
   do not have permission to do that, you may request the second reviewer to merge it for you.

## Code of Conduct

### Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

### Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at tim.trice@gmail.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/