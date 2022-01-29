addressable
============

![R-CMD-check](https://github.com/ropensci/addressable/workflows/R-CMD-check/badge.svg)
[![codecov](https://codecov.io/gh/ropensci/addressable/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/addressable)

Email Address Validation

## Install


```r
remotes::install_github("ropensci/addressable@main")
```


```r
library("addressable")
```

## Address


```r
x <- Address$new("User+tag@example.com")
x$host$host_name
#> [1] "example.com"
x$local$local
#> [1] "user+tag"
x$valid()
#> [1] TRUE
x$fail()
#> NULL
```


```r
x <- Address$new("user1")
x$valid()
#> [1] FALSE
x$fail()
#> [1] "Invalid Domain Name"
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/addressable/issues).
* License: MIT
* Get citation information for `addressable` in R doing `citation(package = 'addressable')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[coc]: https://github.com/ropensci/addressable/blob/maddressable/CODE_OF_CONDUCT.md
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

* Submit an issue on the [Issues page](https://github.com/ropensci/addressable/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/addressable.git`
* Make sure to track progress upstream (i.e., on our version of `addressable` at `ropensci/addressable`) by doing `git remote add upstream https://github.com/ropensci/addressable.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/addressable`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
