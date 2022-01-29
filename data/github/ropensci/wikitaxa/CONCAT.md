wikitaxa
========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/wikitaxa)](https://cranchecks.info/pkgs/wikitaxa)
[![R-check](https://github.com/ropensci/wikitaxa/workflows/R-check/badge.svg)](https://github.com/ropensci/wikitaxa/actions/)
[![codecov](https://codecov.io/gh/ropensci/wikitaxa/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/wikitaxa)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/wikitaxa)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/wikitaxa)](https://cran.r-project.org/package=wikitaxa)

`wikitaxa` - taxonomy data from Wikipedia/Wikidata/Wikispecies

`wikitaxa` docs: https://docs.ropensci.org/wikitaxa/

See also the taxize book: https://books.ropensci.org/taxize/


### Low level API

The low level API is meant for power users and gives you more control,
but requires more knowledge.

* `wt_wiki_page()`
* `wt_wiki_page_parse()`
* `wt_wiki_url_build()`
* `wt_wiki_url_parse()`
* `wt_wikispecies_parse()`
* `wt_wikicommons_parse()`
* `wt_wikipedia_parse()`

### High level API

The high level API is meant to be easier and faster to use.

* `wt_data()`
* `wt_data_id()`
* `wt_wikispecies()`
* `wt_wikicommons()`
* `wt_wikipedia()`

Search functions:

* `wt_wikicommons_search()`
* `wt_wikispecies_search()`
* `wt_wikipedia_search()`

## Installation

CRAN version


```r
install.packages("wikitaxa")
```

Dev version


```r
remotes::install_github("ropensci/wikitaxa")
```


```r
library('wikitaxa')
```

## Contributors

* [Ethan Welty](https://github.com/ezwelty)
* [Scott Chamberlain](https://github.com/sckott)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/wikitaxa/issues).
* License: MIT
* Get citation information for `wikitaxa` in R doing `citation(package = 'wikitaxa')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![ropensci](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[coc]: https://github.com/ropensci/wikitaxa/blob/master/CODE_OF_CONDUCT.md
wikitaxa 0.4.0
==============

### MINOR IMPROVEMENTS

* remove docs link to httr that's not used in package (#19)
* remove egs from readme, change vignette name


wikitaxa 0.3.0
==============

### MINOR IMPROVEMENTS

* integration with vcr for test caching for all HTTP requests (#17) (#18)
* link to `taxize` book and `wikitaxa` vignette in readme (#16)

### BUG FIXES

* fix to `wt_wikipedia()` to separate `<br>` tags appropriately (#15)


wikitaxa 0.2.0
==============

### BUG FIXES

* `wt_wikicommons()` fails better now when a page does not exist, and is now consitent with the rest of package (#14)
* `wt_wikicommons()` fixed - classification objects were not working correctly as the data used is a hot mess - tried to improve parsing of that text (#13)
* `wt_data()` fix - was failing due to i think a change in the internal pkg `WikidataR` (#12)


wikitaxa 0.1.4
==============

### NEW FEATURES

* `wt_wikipedia()` and `wt_wikipedia_search()` gain parameter `wiki`
to give the wiki language, which defaults to `en` (#9)

### MINOR IMPROVEMENTS

* move some examples to dontrun (#11)


wikitaxa 0.1.0
==============

### NEW FEATURES

* Released to CRAN
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
(https://contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
## Test environments

* local OS X install, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 2 reverse dependencies.
  (Summary at <https://github.com/ropensci/wikitaxa/blob/master/revdep/README.md>). No problems were found.

---

This version fixes a documentation link to a package that is not in Depends,  Imports, or Suggests in this package.

This is a re-submission of this version, removing Remotes from DESCRIPTION.

Thanks!
Scott Chamberlain
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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/wikitaxa/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/wikitaxa.git`
* Make sure to track progress upstream (i.e., on our version of `wikitaxa` at `ropensci/wikitaxa`) by doing `git remote add upstream https://github.com/ropensci/wikitaxa.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/wikitaxa`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

```

```
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.2 (2020-06-22) |
|os       |macOS Catalina 10.15.5       |
|system   |x86_64, darwin17.0           |
|ui       |X11                          |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |US/Pacific                   |
|date     |2020-06-28                   |

# Dependencies

|package  |old   |new   |Δ  |
|:--------|:-----|:-----|:--|
|wikitaxa |0.3.0 |0.4.0 |*  |
|openssl  |NA    |1.4.1 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*