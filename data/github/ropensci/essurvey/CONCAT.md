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

## essurvey <img src="man/figures/ess_logo.png" align="right" />

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/essurvey)](https://cran.r-project.org/package=essurvey)
[![R build
status](https://github.com/ropensci/essurvey/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/essurvey/actions)
[![Coverage
status](https://codecov.io/gh/ropensci/essurvey/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/essurvey?branch=master)
[![rOpensci\_Badge](https://badges.ropensci.org/201_status.svg)](https://github.com/ropensci/software-review/issues/201)

## Description

The European Social Survey (ESS) is an academically driven
cross-national survey that has been conducted across Europe since its
establishment in 2001. Every two years, face-to-face interviews are
conducted with newly selected, cross-sectional samples. The survey
measures the attitudes, beliefs and behavior patterns of diverse
populations in more than thirty nations. Taken from the [ESS
website](http://www.europeansocialsurvey.org/about/).

Note: The `essurvey` package was originally called `ess`. Since
`essurvey 1.0.0` all `ess_*` functions have been deprecated in favor of
the `import_*` and `download_*` functions. Also, versions less than and
including `essurvey 1.0.1` returned wrong countries. Please install the
latest CRAN/Github version.

The `essurvey` package is designed to download the ESS data as easily as
possible. It has a few helper functions to download rounds (a term
synonym to waves to denote the same survey in different time points),
rounds for a selected country and to show which rounds/countries are
available. Check out the vignette and other documentation in the
[package’s website](https://docs.ropensci.org/essurvey/) for more
detailed examples of the `essurvey` package.

## Installation

You can install and load the development version with these commands:

``` r
# install.packages("devtools") in case you don't have it
devtools::install_github("ropensci/essurvey")
```

or the stable version with:

``` r
install.packages("essurvey")
```

## Usage

First, you need to register at the ESS website, in case you haven’t.
Please visit the
[register](http://www.europeansocialsurvey.org/user/new) section from
the ESS website. If your email is not registered at their website, an
error will be raised prompting you to go register.

Set your valid email as en environment variable.

``` r
set_email("your@email.com")
```

To explore which rounds/countries are present in the ESS use the
`show_*()` family of functions.

``` r
library(essurvey)
show_countries()
#>  [1] "Albania"            "Austria"            "Belgium"           
#>  [4] "Bulgaria"           "Croatia"            "Cyprus"            
#>  [7] "Czechia"            "Denmark"            "Estonia"           
#> [10] "Finland"            "France"             "Germany"           
#> [13] "Greece"             "Hungary"            "Iceland"           
#> [16] "Ireland"            "Israel"             "Italy"             
#> [19] "Kosovo"             "Latvia"             "Lithuania"         
#> [22] "Luxembourg"         "Montenegro"         "Netherlands"       
#> [25] "Norway"             "Poland"             "Portugal"          
#> [28] "Romania"            "Russian Federation" "Serbia"            
#> [31] "Slovakia"           "Slovenia"           "Spain"             
#> [34] "Sweden"             "Switzerland"        "Turkey"            
#> [37] "Ukraine"            "United Kingdom"
```

To download the first round to use in R:

``` r
one_round <- import_rounds(1)
```

This will return a data frame containing the first round. Typically, the
European Social Survey data files comes with a script that recodes
missing values to `NA` for different programs (Stata, SPSS, SAS).

Use `recode_missings` to recode all values automatically.

``` r
library(tidyverse)

one_round <-
  import_rounds(1) %>%
  recode_missings()
```

See the package vignette for greater detail or see the help page with
`?recode_missings`. You can also download several rounds by supplying
the number of rounds.

``` r
five_rounds <- import_rounds(1:5)
```

This will download all latest versions of rounds 1 through 5 and return
a list of length 5 with each round as a data frame inside the list.

You can check the available rounds with `show_rounds()` because if you
supply a non existent round, the function will return an error.

``` r
two_rounds <- import_rounds(c(1, 22))
#> Error in round_url(rounds) : 
#> ESS round 22 is not a available. Check show_rounds() 
```

Alternatively, you can download all available rounds with
`import_all_rounds()`.

You can also download rounds by country:

``` r
dk_two <- import_country("Denmark", 1:2)
```

Use `show_countries()` to see available countries and
`show_country_rounds("Denmark")` to see available rounds for chosen
country. Alternatively, use `import_all_cntrounds()` to download all
available rounds of a country.

You should be be aware that data from the ESS survey should by analyzed
by taking into consideration the sampling and weights of the survey. A
useful example comes from the work of Anthony Damico and Daniel Oberski
[here](http://asdfree.com/european-social-survey-ess.html).

## Stata, SPSS and SAS users

I’m quite aware that most ESS users don’t know R, that is why the
package also allows to download the data in Stata, SPSS or SAS format
with just one line of code. Instead of the `import_*` functions, use the
`download_*` functions.

``` r
download_rounds(c(1, 2),
                output_dir = "my/new/directory",
                format = 'spss')
```

This will save the ESS rounds into separate folders and unzip them in
the specified directory (if you want to know your current directory,
type `getwd()`). This works the same way for `download_country()`. Be
aware that if you download the files manually you should read them into
R with the `haven` package for all `essurvey` related functions to work.

-----

Please note that this project is released with a [Contributor Code of
Conduct](https://docs.ropensci.org/essurvey/CONDUCT.html). By
participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# essurvey 1.0.8.9999

# essurvey 1.0.8

- Fixes URL creation bug (#53)

# essurvey 1.0.7

- CRAN maintenance release. All vignettes are now precompiled to avoid errors when the ESS website breaks for some reason.
- `ess_email` environmental variable has been renamed to `ESS_EMAIL` to comply with Github Actions standards

# essurvey 1.0.6

- CRAN maintenance release to fix Solaris warnings `Warning in engine$weave(file, quiet = quiet, encoding = enc) : Pandoc (>= 1.12.3) and/or pandoc-citeproc not available. Falling back to R Markdown v1` on CRAN. Tested on Rhub and all passes OK, notifying CRAN.

- Removes automatic citation message when loading package. It's actually annoying.

# essurvey 1.0.5

CRAN maintenance release to add more informative message when the status code of the HTTP request of 'www.europeansocialsurvey.org' is more than 300.

## Minor changes

All tests/examples are now excluded from running on CRAN based on the warning from Brian Ripley:

'Packages which use Internet resources should fail gracefully with an
informative message if the resource is not available (and not give a check
warning nor error).'

They are all forced to run on Travis and Appveyor and this is made clear on the `cran-comments.md`

# essurvey 1.0.4

CRAN maintenance check after release of ESS round 9.

# essurvey 1.0.3

## Breaking changes

* If you don't know which format is available for a round/country, `import_*` and `download_*` functions now accept a NULL argument which runs through `'stata'`, `'spss'` and `'sas'` formats automatically. By default, `import_*` functions have now format set to `NULL` to automatically try the three different formats. This breaks backward dependency but only slightly where it had 'stata' set as default.

## New features

* Users can now download SDDF (weight data) for each country/round combination of files. Functions `show_sddf_cntrounds`, `import_sddf_country` and `download_sddf_country` are now introduced. For technical purposes, `show_sddf_cntrounds` needs for the user to have set their registered ESS email with `set_email`. [#9]

## Minor changes

* Bumps `haven` to minimum package version 2.1.1
* New package website at https://docs.ropensci.org/essurvey

## Internal

* `read_format_data` now tries to read data using `haven` but falls backs to `foreign` in case there's an error. This should only work for SDDF data [#38].
* `read_format_data` and `read_sddf_data` now always return a list. Checking the length of data to return a data frame now happens within each `import_*` function.

## Bug fixes

* Removes an unnecessary if statement in `set_email` that didn't allow to overwrite the email once set.

# essurvey 1.0.2

## Minor changes

* `show_country_rounds` checks if there are missing values and excludes them.

## Breaking changes

`import_all_cntrounds` and `import_country` returned incorrect countries [#31]

# essurvey 1.0.1

## Minor changes

* `ess_email` is now checked that it is not `""` because it wasn't raising an error before.

* Removes the `round` argument from `import_all_cntrounds` because it was a mistake. It already grabs the rounds internally.

# essurvey 1.0.1

Minor release

* Fixes test that checks the number of rounds that each country has. This test was a mistake
because the rounds will change as time passes by and precise country rounds shouldn't be
tested.

# essurvey 1.0.0

The `ess` package has been renamed to `essurvey` for a name conflict with Emacs Speaks Statistics (ESS). See R-pkg mailing list, the post related to the release of ess-0-0-1.

## Breaking changes

* `ess_rounds` and `ess_all_rounds` are deprecated and will be removed in the next release. Use `import_rounds` instead [#22]

* `ess_country` and `ess_all_cntrounds` are deprecated and will be removed in the next release. Use `import_countries` instead [#22]

* The `your_email` argument name of `ess_*` functions has be changed to `ess_email` [#23]

## New features

* `import_rounds`, `import_all_rounds` and `download_rounds` have been introduced as
replacements of `ess_rounds` and `ess_all_rounds`. Same changes were repeated for
`ess_country` and `ess_all_cntrounds` [#22]

* `set_email` to set your email as environmental variable rather than write it in each call [#23]

* All requests to the ESS website are now done through HTTPS rather than HTTP [#24]

* Add package level documentation [#20]

## Minor changes

* `ess_email` had no default value but now has `NULL` as default [#23]

* The `format` argument is now checked through `match.arg` rather than manual check [#25]

# ess 0.1.1 (2018-03-05)

## Breaking changes

* Downloading 1 round both for countries or single rounds now returns a data frame rather than a list. If download is more than two rounds it returns a list. [#8]

## New features

* remove_missings() together with remove_numeric_missings() and remove_character_missings() now allow you to recode the typical categories 'Not applicable', 'Don't know', etc.. into NA's. See the vignette example for more details. [#1]

* Can download files in 'stata', 'spss' and 'sas' formats for all functions (both for downloading to user's directory and for reading data). [#11]

* show_themes() and show_theme_rounds() now available to see which themes have been included in which rounds. [#7]

* show_rounds_country() is now available to see which countries participated in which rounds [#14]

## Bug fixes

* The `ouput_dir` argument is now set to `getwd()` rather than `NULL` as default. [#16]

* When parsing country rounds from the ESS table from the website, shaded dots were being interpreted as valid rounds when in fact they're not. show_* funs new exclude shaded dots until they've been added as valid rounds

* If any `ess_*` function can not connect to the ESS website they will return an explicit R error. [#12]

* `ess_all_cntrounds` and `ess_all_rounds` were returning the directory of each of the files. Now they only return the single directory where the files where saved as a message

# ess 0.0.1 (2017-11-07)

First release
## Test environments
- local Ubuntu 20.0.0 LTS, 4.0.3
- Windows Server 2019 (r-release)
- Mac OS X (r-release)
- Ubuntu 20.04.2 (r-release)
- Ubuntu 20.04.2 (r-devel)


## R CMD check results

- - 0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Reverse dependencies

There are currently no downstream dependencies for this package.

---

Fixed bugs from the library. Standard patch release.

- All tests are run weekly on Github Actions, which are available at https://github.com/ropensci/essurvey/action
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

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

## Best Practices
<!--- Did you remember to include documentation, examples andtests? Unless you're just changing
grammar, please include new tests for your change -->
The following have been updated or added as needed:
[ ] Documentation
[ ] Examples in documentation
[ ] Vignettes
[ ] `testthat` Tests
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/essurvey/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/ropensci/essurvey.git`
* Make sure to track progress upstream (i.e., on our version of `essurvey` at `ropensci/essurvey`) by doing `git remote add upstream https://github.com/ropensci/essurvey.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/essurvey`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!## Expected Behavior
<!--- If you're describing a bug, tell us what should happen -->
<!--- If you're suggesting a change/improvement, tell us how it should work -->

## Current Behavior
<!--- If describing a bug, tell us what happens instead of the expected behavior -->
<!--- If suggesting a change/improvement, explain the difference from current behavior -->

## Steps to Reproduce (for bugs)
<!--- Provide a link to a live example, or an unambiguous set of steps to -->
<!--- reproduce this bug. Include code to reproduce, if relevant -->
1.
2.
3.
4.

## Possible Solution
<!--- Not obligatory, but suggest a fix/reason for the bug, -->
<!--- or ideas how to implement the addition or change -->

## Context
<!--- How has this issue affected you? What are you trying to accomplish? -->
<!--- Providing context helps us come up with a solution that is most useful in the real world -->

## Your Environment
<!--- Include the output of "devtools::session_info()" to help us understand your system environment -->
