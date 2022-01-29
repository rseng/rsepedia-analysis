cchecks
=======



[![R-CMD-check](https://github.com/ropensci/cchecks/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/cchecks/actions/)
[![codecov.io](https://codecov.io/github/ropensci/cchecks/coverage.svg?branch=master)](https://codecov.io/github/ropensci/cchecks?branch=master)

R client for the CRAN checks API at <https://cranchecks.info>

[CRAN checks API docs][docs]

authentication is only needed for the CRAN checks API for the functions that start with `cchn`

See https://docs.ropensci.org/cchecks for full documentation on `cchecks`

## Install


```r
remotes::install_github("ropensci/cchecks")
```


```r
library("cchecks")
```

## heartbeat

- `cch_heartbeat()`

## packages

- current day package check data: `cch_pkgs()` or `cch_pkgs("packagename")`
- historical package check data (30 days back): `cch_pkgs_history()`
- historical package check data for all packages by day: `cch_history()`
- search historical data: `cch_pkgs_search()`

There's an important shortcoming of historical data. The links in the historical
data in the `checks` field are not date specific. If you go to a link in historical
data, for example for April 2nd, 2020, links in that set of data link to whatever
the current check data is for that package. The `check_details` field is
date specific though; the text is scraped from the package checks page each day
and stored, so you can count on that to be date specific. There are sometimes
links to further checks, often of compiled packages on various types of checks
that CRAN runs; we do not have those check results - we could get them but 
have not take the time to sort that out.


## maintainers

- all maintainers: `cch_maintainers()`
- maintainers by email: `cch_maintainers("maelle.salmon_at_yahoo.se")`

## notifications

Functions for working with notifications are all prefixed with `cchn`. 

`cchn` functions are designed to be used from within an R package directory. The functions
look for the package name and maintainer email address. Functions copy heavily from 
https://github.com/r-hub/rhub

The functions

- `cchn_register()`: registration
- `cchn_pkg_rule_list()`/`cchn_rule_list()`: list your own rules
- `cchn_pkg_rule_get()`/`cchn_rule_get()`: get a rule by id
- `cchn_pkg_rule_add()`/`cchn_rule_add()`: create a rule
- `cchn_pkg_rule_delete()`/`cchn_rule_delete()`: delete a rule by id (get id from `cchn_pkg_rule_list`/`cchn_rule_list`)

Functions prefixed with `cchn_pkg_` operate within a package 
directory. That is, your current working directory is an R
package, and is the package for which you want to handle CRAN checks 
notifications. These functions make sure that you are inside of 
an R package, and use the email address and package name based
on the directory you're in.

Functions prefixed with just `cchn_` do not operate within a package.
These functions do not guess package name at all, but require the user
to supply a package name (for those functions that require a package name);
and instead of guessing an email address from your package, we guess email
from the cached cchecks email file (see `?cchn_register`).

The first thing to do is to register an email address. In an R session in a working directory for 
one of your packages that is on CRAN, run `cchn_register()`. This function:

- registers your email address
- a validation email is sent to you right away with your token
- paste the token from the amil into the R session
- the email and token are then saved in a file on your machine
- all `cchn_rule*` functions use this cached token
- you don't need to pass the token in any `cchn` function calls

If you run `cchn_register()` in the same package directory (with the same email address), 
you'll be issued a new token, which will be updated in your cached token file.

It's entirely possible to have more than one email address you use across different R packages. 
If you run `cchn_register()` in a different package directory (with a different email address
from a previous run of `cchn_register()`), you'll be issued a different token associated 
with that new email address.

See `?cchn_rules` for details on how the rules work and many examples of adding rules.

Note that you can only manage your own rules. You can not list, get, or delete rules 
of other users.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/cchecks/issues).
* License: MIT
* Get citation information for `cchecks` in R doing `citation(package = 'cchecks')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[docs]: https://cranchecks.info/docs
[coc]: https://github.com/ropensci/cchecks/blob/master/CODE_OF_CONDUCT.md
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

* Submit an issue on the [Issues page](https://github.com/ropensci/cchecks/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/cchecks.git`
* Make sure to track progress upstream (i.e., on our version of `cchecks` at `ropensci/cchecks`) by doing `git remote add upstream https://github.com/ropensci/cchecks.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/cchecks`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
