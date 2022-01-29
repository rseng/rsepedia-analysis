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
helminthR
=======

[![R build status](https://github.com/ropensci/helminthR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/helminthR/actions)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/rmq9euldm5gy9qup?svg=true)](https://ci.appveyor.com/project/taddallas/helminthr)
[![codecov.io](https://codecov.io/github/ropensci/helminthR/coverage.svg?branch=master)](https://codecov.io/github/ropensci/helminthR?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/helminthR)](https://github.com/r-hub/cranlogs.app)


> Programmatically access the London Natural History Museum's [helminth database](https://www.nhm.ac.uk/research-curation/scientific-resources/taxonomy-systematics/host-parasites/index.html).

See software note in _Ecography_ ([available here](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.02131))


### Installation

From GitHub


```r
# install.packages("devtools")
devtools::install_github("rOpenSci/helminthR")
library("helminthR")
```

From CRAN


```r
install.packages("helminthR")
```



### Main functions

#### `findHost()`

Given a host genus and (optionally) species and location, this function returns all host-parasite associations of a given host species. The example below determines all parasite records for helminth infections of _Gorilla gorilla_.


```r
gorillaParasites <- findHost('Gorilla', 'gorilla')
head(gorillaParasites)
```

#### `findParasite()`

Given a helminth parasite genus (and optionally species, and location), this function returns a list of host-parasite records for that parasite. In the example below, I query the database for occurrences of the genus _Strongyloides_.


```r
strongHosts <- findParasite(genus='Strongyloides')
str(strongHosts)
```



#### `listLocations()` and `findLocation()`

List all location names (`listLocations()`). These names can be given to the `findLocation()` function, which finds all host-parasite associations that have occurred in the given location. Below, I look at host-parasite associations recorded in France.



```r
FrenchHostPars <- findLocation(location='France')
str(FrenchHostPars)
```




### Contribute!

Feel free to fork it and contribute some functionality.  



## Meta

* Please [report any issues or bugs](https://github.com/ropensci/helminthR/issues).
* License: GPL-3
* Get citation information for `helminthR` in R doing `citation(package = 'helminthR')`
* Please note that this project is released with a [Contributor Code of Conduct](https://www.contributor-covenant.org/).
By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
helminthR 1.0.9
==============

Fixing some NOTES and WARNINGS in the vignette creation process for the package. Switched from Travis CI to GitHub Actions for CI on GitHub. 



helminthR 1.0.8
==============

Small fix to the cleanData function to now use taxize to get host species taxonomic data.


helminthR 1.0.7
==============

Removed geocoding functionality previously present in listLocations(), as this now requires an API key. A cached version of the geographic coordinates of locations is provided as package data (`data(locations)`). 

Added extra catch in `cleanDat.R` to remove species who are identified as "something spp." instead of just removing those identified as "something sp.". 




helminthR 1.0.6
==============
* bug fix that was causing null results for some location specifications. 





helminthR 1.0.5
==============

* Released to CRAN.
## Test environments

* ubuntu 20.04-release
* ubuntu 20.04-devel
* windows-latest
* macOS-latest


## R CMD check results

R CMD check results
0 errors | 0 warnings | 0 notes


## Reverse dependencies

There are no reverse dependencies


---

I have read and agree to the the CRAN
policies at https://cran.r-project.org/web/packages/policies.html

This is an update to handle some NOTES and a few WARNINGS that were being thrown by certain operating systems when compiling the vignette. 

Thanks!
Tad Dallas
