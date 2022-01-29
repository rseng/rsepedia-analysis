pendulum
=====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/pendulum/workflows/R-check/badge.svg)](https://github.com/ropensci/pendulum/actions?query=workflow%3AR-check)

Time classes for R, w/ time mocking capability via [timefuzz][]

**Why?** - Some programming languages support [monkey patching][monkey]. However, R does not officially allow this (there's a way to do it, but it's a hack). Therefore, the only way that I can think of to "mock time" is to use a different set of time classes/objects. That is, we can't just mock time using `Sys.Date()` directly as there's no way to monkey patch the code behind that function.  Some may say it's good that R doesn't officially support monkey patching; there's definitely upsides; but a downside is that mocking things like time becomes difficult.

- `clock()` is designed for setting a specific time
- `now()` is designed for getting the current time AND is mockable
- `sys_time()` and `sys_date()` are drop in replacements for `Sys.time()` and `Sys.Date()` AND are mockable

See the docs to get started: https://docs.ropensci.org/pendulum/

## Installation


```r
remotes::install_github("ropensci/pendulum")
```


```r
library(pendulum)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/pendulum/issues).
* License: MIT
* Get citation information for `pendulum` in R doing `citation(package = 'pendulum')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[timefuzz]: https://github.com/ropensci/timefuzz
[monkey]: https://en.wikipedia.org/wiki/Monkey_patch
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

* Submit an issue on the [Issues page](https://github.com/ropensci/pendulum/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/pendulum.git`
* Make sure to track progress upstream (i.e., on our version of `pendulum` at `ropensci/pendulum`) by doing `git remote add upstream https://github.com/ropensci/pendulum.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/pendulum`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
