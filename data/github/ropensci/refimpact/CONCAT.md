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
refimpact
=========

[![Build Status](https://travis-ci.org/ropensci/refimpact.svg?branch=master)](https://travis-ci.org/ropensci/refimpact) 
[![Build status](https://ci.appveyor.com/api/projects/status/jxj1yela4a6ym6wb/branch/master?svg=true)](https://ci.appveyor.com/project/perrystephenson/refimpact/branch/master) 
[![codecov](https://codecov.io/gh/ropensci/refimpact/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/refimpact) 
[![](https://badges.ropensci.org/78_status.svg)](https://github.com/ropensci/onboarding/issues/78) 
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/refimpact)](https://CRAN.R-project.org/package=refimpact)

**refimpact** provides an API wrapper for the UK Research Excellence Framework 2014 Impact Case Studies Database. You can find more information about this database at <http://impact.ref.ac.uk/CaseStudies/>.

The data may be of interest to you if you are interested in:

-   text mining
-   directed graphs
-   policies for research funding

Case studies in the database are licenced under a CC-BY 4.0 license. The full license can be found [here](https://creativecommons.org/licenses/by/4.0/legalcode) and a more user-friendly version of the license can be be obtained [here](https://creativecommons.org/licenses/by/4.0/).

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

Installation
------------

### Install from CRAN

``` r
install.packages("refimpact")
```

### Install from Github

``` r
install.packages("devtools")
devtools::install_github("perrystephenson/refimpact")
```

Usage
-----

See the vignette:

``` r
vignette("refimpact")
```

More Information
----------------

For more information about a specific function you can use the help commands (for example `?ref_get`).

To raise bug reports and issues, please use the issue tracker in Github.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# refimpact 1.0.0

* Major breaking changes. There is now a single user-facing function `ref_get()`
  which takes the API method as an argument. This standardises a lot of the 
  input validation and error handling, as well as reducing the risk of bugs (as 
  there are less lines of code). Functions from previous version of the package
  are still available, but deprecated.
* A vignette has been added, and the help documentation has been improved.
* The package now uses **httr** when calling the API, which improves reliability
  and provides much better error messages to the end-user when things go wrong.
* The `phrase` parameter to the SearchCaseStudies method now allows text queries
  of any length and complexity
* Bundled a `ref_tags` dataset with the package, to save the end-user from 
  having to iterate through the ListTagValues API method in order to find tags
  to use as parameters when searching the database
* Added a contributor code of conduct
* The entire package was re-written.

# refimpact 0.1.0

* Initial release.

## Test environments

* local OS X install, R 3.4.1
* Ubuntu 12.04 (on travis-ci), R 3.4.0, R 3.3.3, R-devel.
* Windows (on AppVeyor), R 3.4.1, R 3.3.3, R-devel

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTEs

devtools::win_builder() shows three potential spelling errors - each of these 
has been checked and confirmed accurate.

## Downstream dependencies

There are currently no downstream dependencies for this package.
