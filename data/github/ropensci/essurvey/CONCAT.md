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
---
output:
  github_document:
    html_preview: false
---
```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-",
  eval = FALSE
)
options(tibble.print_min = 5, tibble.print_max = 5)
```

## essurvey <img src="man/figures/ess_logo.png" align="right" />

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/essurvey)](https://cran.r-project.org/package=essurvey)
[![R build status](https://github.com/ropensci/essurvey/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/essurvey/actions)
[![Coverage status](https://codecov.io/gh/ropensci/essurvey/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/essurvey?branch=master)
[![rOpensci_Badge](https://badges.ropensci.org/201_status.svg)](https://github.com/ropensci/software-review/issues/201)


## Description

The European Social Survey (ESS) is an academically driven cross-national survey that has been conducted across Europe since its establishment in 2001. Every two years, face-to-face interviews are conducted with newly selected, cross-sectional samples. The survey measures the attitudes, beliefs and behavior patterns of diverse populations in more than thirty nations. Taken from the [ESS website](http://www.europeansocialsurvey.org/about/).

Note: The `essurvey` package was originally called `ess`. Since `essurvey 1.0.0` all `ess_*` functions have been deprecated in favor of the `import_*` and `download_*` functions. Also, versions less than and including `essurvey 1.0.1` returned wrong countries. Please install the latest CRAN/Github version.

The `essurvey` package is designed to download the ESS data as easily as possible. It has a few helper functions to download rounds (a term synonym to waves to denote the same survey in different time points), rounds for a selected country and to show which rounds/countries are available. Check out the vignette and other documentation in the [package's website](https://docs.ropensci.org/essurvey/) for more detailed examples of the `essurvey` package.

## Installation

You can install and load the development version with these commands:

```{r}
# install.packages("devtools") in case you don't have it
devtools::install_github("ropensci/essurvey")
```

or the stable version with:

```{r}
install.packages("essurvey")
```

## Usage

First, you need to register at the ESS website, in case you haven't. Please visit the [register](http://www.europeansocialsurvey.org/user/new) section from the ESS website. If your email is not registered at their website, an error will be raised prompting you to go register.

Set your valid email as en environment variable.

```{r}
set_email("your@email.com")
```

To explore which rounds/countries are present in the ESS use the `show_*()` family of functions.

```{r, eval = TRUE}
library(essurvey)
show_countries()
```

To download the first round to use in R:

```{r}
one_round <- import_rounds(1)
```

This will return a data frame containing the first round. Typically, the European Social Survey data files comes with a script that recodes missing values to `NA` for different programs (Stata, SPSS, SAS).

Use `recode_missings` to recode all values automatically.

```{r}
library(tidyverse)

one_round <-
  import_rounds(1) %>%
  recode_missings()
```

See the package vignette for greater detail or see the help page with `?recode_missings`. You can also download several rounds by supplying the number of rounds.

```{r}
five_rounds <- import_rounds(1:5)
```

This will download all latest versions of rounds 1 through 5 and return a list of length 5 with each round as a data frame inside the list. 

You can check the available rounds with `show_rounds()` because if you supply a non existent round, the function will return an error.

```{r}
two_rounds <- import_rounds(c(1, 22))
#> Error in round_url(rounds) : 
#> ESS round 22 is not a available. Check show_rounds() 
```

Alternatively, you can download all available rounds with `import_all_rounds()`.

You can also download rounds by country:

```{r}
dk_two <- import_country("Denmark", 1:2)
```

Use `show_countries()` to see available countries and `show_country_rounds("Denmark")` to see available rounds for chosen country. Alternatively, use `import_all_cntrounds()` to download all available rounds of a country.

You should be be aware that data from the ESS survey should by analyzed by taking into consideration the sampling and weights of the survey. A useful example comes from the work of Anthony Damico and Daniel Oberski [here](http://asdfree.com/european-social-survey-ess.html).

## Stata, SPSS and SAS users

I'm quite aware that most ESS users don't know R, that is why the package also allows to download the data in Stata, SPSS or SAS format with just one line of code. Instead of the `import_*` functions, use the `download_*` functions.
```{r}
download_rounds(c(1, 2),
                output_dir = "my/new/directory",
                format = 'spss')
```

This will save the ESS rounds into separate folders and unzip them in the specified directory (if you want to know your current directory, type `getwd()`). This works the same way for `download_country()`. Be aware that if you download the files manually you should read them into R with the `haven` package for all `essurvey` related functions to work.

---

Please note that this project is released with a [Contributor Code of Conduct](https://docs.ropensci.org/essurvey/CONDUCT.html). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Introduction to the essurvey package"
author: "Jorge Cimentada"
date: "2021-03-09"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to the essurvey package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Introduction

Using the `essurvey` package is fairly easy. There are are two main families of functions: `import_*` and `show_*`. They each complement each other and allow the user to almost never have to go to the European Social Survey (ESS) website. The only scenario where you need to enter the ESS website is to validate your email. If you haven't registered, create an account at http://www.europeansocialsurvey.org/user/new. For those unfamiliar with the ESS, this vignette uses the term rounds, here a synonym of waves to denote the same survey in different time points.

Once you register visit your email account to validate the account and you're ready to access the data.

Given that some `essurvey` functions require your email address, this vignette will use a fake email but everything should work accordingly if you registered with the ESS.

## Downloading country specific rounds

Note: versions less than and including `essurvey 1.0.1` returned wrong countries. Please install the latest CRAN/Github version.

To install and load development version of the package use:

```r
# install.packages("devtools")
devtools::install_github("ropensci/essurvey")
```

to install the stable version from CRAN use:


```r
install.packages("essurvey")
```



Downloading the ESS data requires validating your email every time you download data. We can set our email as an environment variable with `set_email`.


```r
set_email("your@email.com")
```

Once that's executed you can delete the previous line and any `import_*` call will look for the email automatically, stored as an environment variable.

Let's suppose you don't know which countries or rounds are available for the ESS. Then the `show_*` family of functions is your friend.

To find out which countries have participated you can use `show_countries()`


```r
show_countries()
```

```
##  [1] "Albania"            "Austria"            "Belgium"           
##  [4] "Bulgaria"           "Croatia"            "Cyprus"            
##  [7] "Czechia"            "Denmark"            "Estonia"           
## [10] "Finland"            "France"             "Germany"           
## [13] "Greece"             "Hungary"            "Iceland"           
## [16] "Ireland"            "Israel"             "Italy"             
## [19] "Kosovo"             "Latvia"             "Lithuania"         
## [22] "Luxembourg"         "Montenegro"         "Netherlands"       
## [25] "Norway"             "Poland"             "Portugal"          
## [28] "Romania"            "Russian Federation" "Serbia"            
## [31] "Slovakia"           "Slovenia"           "Spain"             
## [34] "Sweden"             "Switzerland"        "Turkey"            
## [37] "Ukraine"            "United Kingdom"
```

This function actually looks up the countries in the ESS website. If new countries enter, this will automatically grab those countries as well. Let's check out Turkey. How many rounds has Turkey participated in? We can use `show_country_rounds()`


```r
tk_rnds <- show_country_rounds("Turkey")
tk_rnds
```

```
## [1] 2 4
```
Note that country names are case sensitive. Use the exact name printed out by `show_countries()`

Using this information, we can download those specific rounds easily with `import_country`. Since `essurvey 1.0.0` all `ess_*` functions have been deprecated in favor of the `import_*` and `download_*` functions.


```r
turkey <-
  import_country(
    country = "Turkey",
    rounds = c(2, 4)
    )
```

`turkey` will now be a list of `length(rounds)` containing a data frame for each round. If you only specified one round, then all `import_*` functions return a data frame. `import_country` is useful for when you want to download specific rounds, but not all. To download all rounds for a country automatically you can use `import_all_cntrounds`.


```r
import_all_cntrounds("Turkey")
```

The `import_*` family is  concerned with downloading the data and thus always returns a list containing data frames unless only one round is specified, in which it returns a `tibble`. Conversely, the `show_*` family grabs information from the ESS website and always returns vectors.

## Download complete rounds

Similarly, we can use other functions to download rounds. To see which rounds are currently available, use `show_rounds`.


```r
show_rounds()
```

```
## [1] 1 2 3 4 5 6 7 8 9
```

Similar to `show_countries`, `show_rounds` interactively looks up rounds in the ESS website, so any future rounds will automatically be included.

To download all available rounds, use `import_all_rounds`


```r
all_rounds <- import_all_rounds()
```

Alternatively, use `import_rounds` for selected ones.


```r
selected_rounds <- import_rounds(c(1, 3, 6))
```

## Downloading data for Stata, SPSS and SAS users

All `import_*` functions have an equivalent `download_*` function that allows the user to save the datasets in a specified folder in `'stata'`, `'spss'` or `'sas'` formats.

For example, to save round two from Turkey in a folder called `./my_folder`, we use:


```r
download_country("Turkey", 2,
                 output_dir = "./myfolder/")
```

By default it saves the data as `'stata'` files. Alternatively you can use `'spss'` or `'sas'`.


```r
download_country("Turkey", 2,
                 output_dir = "./myfolder/",
                 format = 'sas')
```

This will save the data to `./myfolder/ESS_Turkey` and inside that folder there will be the `ESS2` folder that contains the data.

## Correcting for missing values

Whenever you download the ESS data, it comes together with a script that recodes the values 6 = 'Not applicable', 7 = 'Refusal', 8 = 'Don't know', 9 = 'No answer' and 9 = 'Not available' as missings. However, that is the case for variables that have a scaling of 1-5. For variables which have a scaling from 1-10 the corresponding missings are 66, 77, and so on. At first glance new users might not know this and start calculating statistics with these variables such as...


```r
sp <- import_country("Spain", 1)
mean(sp$tvtot)
# 4.622406
```

..but that vector contains numbers such as `66`, `77`, that shouldn't be there. `recode_missings()` removes the corresponding missings for numeric variables as well as for character variables. It accepts the complete `tibble` and recodes all variables that should be recoded.


```r
new_coding <- recode_missings(sp)
mean(new_coding$tvtot, na.rm = TRUE)
# 4.527504
```

It also gives you the option of recoding only specific categories. For example...


```r
other_newcoding <- recode_missings(sp, c("Don't know", "Refusal"))
table(other_newcoding$tvpol)
#  0   1   2   3   4   5   6   7  66 
# 167 460 610 252  95  36  26  31  45 
```

...still has missing values but recoded the ones that were specified. I strongly suggest the user not to recode these categories as missing without looking at the data as there might be substantial differences between people who didn't and who did answer questions. If the user is decided to do so, use `recode_missings` to recode everything and the corresponding `recode_*_missings` functions for numeric and character recodings separately. See the documentation of `?recode_missings` for more information.

## Analyzing ESS data

Be aware that for analyzing data from the ESS survey you should take into consideration the sampling and weights of the survey. The [survey](http://r-survey.r-forge.r-project.org/survey/) package provides very good support for this. A useful example comes from the work of Anthony Damico and Daniel Oberski [here](http://asdfree.com/european-social-survey-ess.html).

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_any_rounds.R
\name{show_country_rounds}
\alias{show_country_rounds}
\title{Return available rounds for a country in the European Social Survey}
\usage{
show_country_rounds(country)
}
\arguments{
\item{country}{A character of length 1 with the full name of the country.
Use \code{\link{show_countries}}for a list of available countries.}
}
\value{
numeric vector with available rounds for \code{country}
}
\description{
Return available rounds for a country in the European Social Survey
}
\examples{

\dontrun{

show_country_rounds("Spain")

show_country_rounds("Turkey")

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_sddf_cntrounds.R
\name{show_sddf_cntrounds}
\alias{show_sddf_cntrounds}
\title{Return available SDDF rounds for a country in the European Social Survey}
\usage{
show_sddf_cntrounds(country, ess_email = NULL)
}
\arguments{
\item{country}{A character of length 1 with the full name of the country.
Use \code{\link{show_countries}} for a list of available countries.}

\item{ess_email}{a character vector with your email, such as "your_email@email.com".
If you haven't registered in the ESS website, create an account at 
\url{http://www.europeansocialsurvey.org/user/new}. A preferred method is to login
through \code{\link{set_email}}.}
}
\value{
numeric vector with available rounds for \code{country}
}
\description{
Return available SDDF rounds for a country in the European Social Survey
}
\details{
SDDF data are the equivalent weight data used to analyze the European Social Survey
properly. For more information, see the details section of \code{\link{import_sddf_country}}.
As an exception to the \code{show_*} family of functions, \code{show_sddf rounds}
needs your ESS email to check which rounds are available. Be sure to add it
with \code{\link{set_email}}.
}
\examples{

\dontrun{
set_email("your_email@email.com")

show_sddf_cntrounds("Spain")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_funs.R
\name{show_countries}
\alias{show_countries}
\title{Return available countries in the European Social Survey}
\usage{
show_countries()
}
\value{
character vector with available countries
}
\description{
Return available countries in the European Social Survey
}
\examples{

\dontrun{
show_countries()
}


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_any_rounds.R
\name{show_rounds_country}
\alias{show_rounds_country}
\title{Return countries that participated in \strong{all} of the specified rounds.}
\usage{
show_rounds_country(rounds, participate = TRUE)
}
\arguments{
\item{rounds}{A numeric vector specifying the rounds from which to return the countries.
Use \code{\link{show_rounds}}for a list of available rounds.}

\item{participate}{A logical that controls whether to show participating countries in that/those
rounds or countries that didn't participate. Set to \code{TRUE} by default.}
}
\value{
A character vector with the country names
}
\description{
Return countries that participated in \strong{all} of the specified rounds.
}
\details{
\code{show_rounds_country} returns the countries that participated in
\strong{all} of the specified rounds. That is, \code{show_rounds_country(1:2)}
will return countries that participated both in round 1 and round 2. Conversely,
if \code{participate = FALSE} it will return the countries that did not
participate in \strong{both} round 1 and round 2.
}
\examples{

\dontrun{

# Return countries that participated in round 2

show_rounds_country(2)

# Return countries that participated in all rounds

show_rounds_country(1:8)

# Return countries that didn't participate in the first three rounds

show_rounds_country(1:3, participate = FALSE)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_sddf_country.R
\name{import_sddf_country}
\alias{import_sddf_country}
\alias{import_all_sddf_cntrounds}
\alias{download_sddf_country}
\title{Download SDDF data by round for countries from the European Social Survey}
\usage{
import_sddf_country(country, rounds, ess_email = NULL, format = NULL)

import_all_sddf_cntrounds(country, ess_email = NULL, format = NULL)

download_sddf_country(
  country,
  rounds,
  ess_email = NULL,
  output_dir = getwd(),
  format = "stata"
)
}
\arguments{
\item{country}{a character of length 1 with the full name of the country.
Use \code{\link{show_countries}} for a list of available countries.}

\item{rounds}{a numeric vector with the rounds to download. See \code{\link{show_sddf_cntrounds}}
for all available rounds for any given country.}

\item{ess_email}{a character vector with your email, such as "your_email@email.com".
If you haven't registered in the ESS website, create an account at
\url{http://www.europeansocialsurvey.org/user/new}. A preferred method is to login
through \code{\link{set_email}}.}

\item{format}{the format from which to download the data. By default it is NULL for \code{import_*} functions and tries to read 'stata', 'spss' and 'sas' in the specific order. This can be useful if some countries don't have a particular format available.  Alternatively, the user can specify the format which can either be 'stata', 'spss' or 'sas'.
For the \code{download_*} functions it is set to 'stata' because the format should be
specified before downloading. Setting it to \code{NULL} will iterate over 'stata',
'spss' and 'sas' and download the first that is available. When using \code{import_country}
the data will be downloaded and read in the \code{format} specified. For \code{download_country},
the data is downloaded from the specified \code{format} (only 'spss' and 'stata' supported,
see details).}

\item{output_dir}{a character vector with the output directory in case you want to
only download the files using \code{download_sddf_country}. Defaults to your working
directory. This will be interpreted as a \strong{directory} and not a path with
a file name.}
}
\value{
for \code{import_sddf_country} if \code{length(rounds)} is 1, it returns a tibble with
the latest version of that round. Otherwise it returns a list of \code{length(rounds)}
containing the latest version of each round. For \code{download_sddf_country}, if
\code{output_dir} is a valid directory, it returns the saved directories invisibly and saves
all the rounds in the chosen \code{format} in \code{output_dir}
}
\description{
Download SDDF data by round for countries from the European Social Survey
}
\details{
SDDF data (Sample Design Data Files) are data sets that contain additional columns with the
sample design and weights for a given country in a given round. These additional columns are
required to perform any complex weighted analysis of the ESS data. Users interested in using this data
should read the description of SDDF files \href{http://www.europeansocialsurvey.org/methodology/ess_methodology/sampling.html}{here}
and should read \href{http://www.europeansocialsurvey.org/data/download_sample_data.html}{here} for the
sampling design of the country of analysis for that specific round.

Use \code{import_sddf_country} to download the SDDF data by country into R.
\code{import_all_sddf_cntrounds} will download all available SDDF data for a given country by
default and \code{download_sddf_country} will download SDDF data and save them in a specified
\code{format} in the supplied directory.

The \code{format} argument from \code{import_country} should not matter to the user
because the data is read into R either way. However, different formats might have
different handling of the encoding of some questions. This option was preserved
so that the user can switch between formats if any encoding errors are found in the data. For more
details see the discussion \href{https://github.com/ropensci/essurvey/issues/11}{here}.

Additionally, given that the SDDF data is not very complete, some countries do not have SDDF data
in Stata or SPSS formats. For that reason, the \code{format} argument is not used in \code{import_sddf_country}.
Internally, \code{Stata} is chosen over \code{SPSS} and \code{SPSS} over \code{SAS} in that
order of preference.

For this particular argument, 'sas' is not supported because the data formats have
changed between ESS waves and separate formats require different functions to be
read. To preserve parsimony and format errors between waves, the user should use
'stata' or 'spss'.

Starting from round 7 (including), the ESS switched the layout of SDDF data.
Before the rounds, SDDF data was published separately by wave-country
combination. From round 7 onwards, all SDDF data is released as a single
integrated file with all countries combined for that given round. \code{import_sddf_country}
takes care of this nuance by reading the data and filtering the chosen
country automatically. \code{download_sddf_country} downloads the raw file but also
reads the data into memory to subset the specific country requested. This
process should be transparent to the user but beware that reading/writing the data might delete
some of it's properties such as dropping the labels or label attribute.
}
\examples{
\dontrun{

set_email("your_email@email.com")

sp_three <- import_sddf_country("Spain", 5:6)

show_sddf_cntrounds("Spain")

# Only download the files, this will return nothing

temp_dir <- tempdir()

download_sddf_country(
 "Spain",
 rounds = 5:6,
 output_dir = temp_dir
)

# By default, download_sddf_country downloads 'stata' files but
# you can also download 'spss' or 'sas' files.

download_sddf_country(
 "Spain",
 rounds = 1:8,
 output_dir = temp_dir,
 format = 'spss'
)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_email.R
\name{set_email}
\alias{set_email}
\title{Save your ESS email as an environment variable}
\usage{
set_email(ess_email)
}
\arguments{
\item{ess_email}{a character string with your registered email.}
}
\description{
Save your ESS email as an environment variable
}
\details{
You should only run \code{set_email()} once and every \code{import_} and \code{download_} function
should work fine. Make sure your email is registered at
\url{http://www.europeansocialsurvey.org/} before setting the email.
}
\examples{

\dontrun{
set_email("my_registered@email.com")

import_rounds(1)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_funs.R
\name{show_themes}
\alias{show_themes}
\title{Return available themes in the European Social Survey}
\usage{
show_themes()
}
\value{
character vector with available themes
}
\description{
This function returns the available themes in the European Social Survey.
However, contrary to \code{\link{show_countries}} and \code{\link{show_country_rounds}},
themes can not be downloaded as separate datasets. This and
\code{\link{show_theme_rounds}} serve purely for informative purposes.
}
\examples{

\dontrun{
show_themes()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_funs.R
\name{show_rounds}
\alias{show_rounds}
\title{Return available rounds in the European Social Survey}
\usage{
show_rounds()
}
\value{
numeric vector with available rounds
}
\description{
Return available rounds in the European Social Survey
}
\examples{

\dontrun{
show_rounds()
}


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/essurvey-package.R
\docType{package}
\name{essurvey-package}
\alias{essurvey}
\alias{essurvey-package}
\title{essurvey: Download Data from the European Social Survey on the Fly}
\description{
Download data from the European Social Survey directly from their website <http://www.europeansocialsurvey.org/>. There are two families of functions that allow you to download and interactively check all countries and rounds available.
}
\details{
Note that this package is for downloading data from the ESS survey.
For analyzing the data, the user should consider the weights and sampling
design of each country/round combination as well as sampling between rounds.
For some examples, check out the work of Anthony Damico and Daniel Oberski 
\href{http://asdfree.com/european-social-survey-ess.html}{here}. For detailed
examples on how to explore/download data using this package, visit the package website
at \url{https://docs.ropensci.org/essurvey/}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/essurvey/}
  \item \url{https://github.com/ropensci/essurvey}
  \item Report bugs at \url{https://github.com/ropensci/essurvey/issues}
}

}
\author{
\strong{Maintainer}: Jorge Cimentada \email{cimentadaj@gmail.com}

Other contributors:
\itemize{
  \item Thomas Leeper (Thomas reviewed the package for rOpensci,see https://github.com/ropensci/software-review/issues/201) [reviewer]
  \item Nujcharee Haswell (Nujcharee reviewed the package for rOpensci, see https://github.com/ropensci/software-review/issues/201) [reviewer]
  \item Jorge Lopez \email{jorge@loperez.com} [contributor]
  \item François Briatte \email{f.briatte@gmail.com} [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_rounds.R
\name{import_rounds}
\alias{import_rounds}
\alias{import_all_rounds}
\alias{download_rounds}
\title{Download integrated rounds from the European Social Survey}
\usage{
import_rounds(rounds, ess_email = NULL, format = NULL)

import_all_rounds(ess_email = NULL, format = NULL)

download_rounds(
  rounds,
  ess_email = NULL,
  output_dir = getwd(),
  format = "stata"
)
}
\arguments{
\item{rounds}{a numeric vector with the rounds to download. See \code{\link{show_rounds}}
for all available rounds.}

\item{ess_email}{a character vector with your email, such as "your_email@email.com".
If you haven't registered in the ESS website, create an account at
\url{http://www.europeansocialsurvey.org/user/new}. A preferred method is to login
through \code{\link{set_email}}.}

\item{format}{the format from which to download the data. By default it is NULL for \code{import_*} functions and tries to read 'stata', 'spss' and 'sas' in the specific order. This can be useful if some countries don't have a particular format available.  Alternatively, the user can specify the format which can either be 'stata', 'spss' or 'sas'. For the \code{download_*} functions it is set to 'stata' because the format should be specified before downloading. When using \code{import_country} the data will be downloaded and read in the \code{format} specified. For \code{download_country}, the data is downloaded from the specified \code{format} (only 'spss' and 'stata' supported, see details).}

\item{output_dir}{a character vector with the output directory in case you want to only download the files using
the \code{download_rounds}. Defaults to your working directory. This will be interpreted as
a \strong{directory} and not a path with a file name.}
}
\value{
for \code{import_rounds} if \code{length(rounds)} is 1, it returns a tibble
with the latest version of that round. Otherwise it returns a list of \code{length(rounds)}
containing the latest version of each round. For \code{download_rounds}, if
\code{output_dir} is a valid directory, it returns the saved directories invisibly
and saves all the rounds in the chosen \code{format} in \code{output_dir}
}
\description{
Download integrated rounds from the European Social Survey
}
\details{
Use \code{import_rounds} to download specified rounds and import them to R.
\code{import_all_rounds} will download all rounds by default and \code{download_rounds}
will download rounds and save them in a specified \code{format} in the supplied
directory.

The \code{format} argument from \code{import_rounds} should not matter to the user
because the data is read into R either way. However, different formats might have
different handling of the encoding of some questions. This option was preserved
so that the user
can switch between formats if any encoding errors are found in the data. For more
details see the discussion \href{https://github.com/ropensci/essurvey/issues/11}{here}.
For this particular argument in, 'sas' is not supported because the data formats have
changed between ESS waves and separate formats require different functions to be
read. To preserve parsimony and format errors between waves, the user should use
'spss' or 'stata'.
}
\examples{

\dontrun{

set_email("your_email@email.com")

# Get first three rounds
three_rounds <- import_rounds(1:3)

temp_dir <- tempdir()

# Only download the files to output_dir, this will return nothing.
download_rounds(
 rounds = 1:3,
 output_dir = temp_dir,
)

# By default, download_rounds saves a 'stata' file. You can
# also download 'spss' and 'sas' files.

download_rounds(
 rounds = 1:3,
 output_dir = temp_dir,
 format = 'spss'
)

# If rounds are repeated, will download only unique ones
two_rounds <- import_rounds(c(1, 1))

# If email is not registered at ESS website, error will arise
two_rounds <- import_rounds(c(1, 2), "wrong_email@email.com")

# Error in authenticate(ess_email) :
# The email address you provided is not associated with any registered user.
# Create an account at https://www.europeansocialsurvey.org/user/new

# If selected rounds don't exist, error will arise

two_rounds <- import_rounds(c(1, 22))
# Error in round_url(rounds) :
# ESS round 22 is not a available. Check show_rounds()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_any_rounds.R
\name{show_theme_rounds}
\alias{show_theme_rounds}
\title{Return available rounds for a theme in the European Social Survey}
\usage{
show_theme_rounds(theme)
}
\arguments{
\item{theme}{A character of length 1 with the full name of the theme.
Use \code{\link{show_themes}}for a list of available themes.}
}
\value{
numeric vector with available rounds for \code{country}
}
\description{
This function returns the available rounds for any theme from
\code{\link{show_themes}}. However, contrary to \code{\link{show_country_rounds}}
themes can not be downloaded as separate datasets. This and the 
\code{\link{show_themes}} function serve purely for informative purposes.
}
\examples{

\dontrun{
chosen_theme <- show_themes()[3]

# In which rounds was the topic of 'Democracy' asked?
show_theme_rounds(chosen_theme)

# And politics?
show_theme_rounds("Politics")

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_countries.R
\name{import_country}
\alias{import_country}
\alias{import_all_cntrounds}
\alias{download_country}
\title{Download integrated rounds separately for countries from the European Social
Survey}
\usage{
import_country(country, rounds, ess_email = NULL, format = NULL)

import_all_cntrounds(country, ess_email = NULL, format = NULL)

download_country(
  country,
  rounds,
  ess_email = NULL,
  output_dir = getwd(),
  format = "stata"
)
}
\arguments{
\item{country}{a character of length 1 with the full name of the country. 
Use \code{\link{show_countries}} for a list of available countries.}

\item{rounds}{a numeric vector with the rounds to download. See \code{\link{show_rounds}}
for all available rounds.}

\item{ess_email}{a character vector with your email, such as "your_email@email.com".
If you haven't registered in the ESS website, create an account at 
\url{http://www.europeansocialsurvey.org/user/new}. A preferred method is to login
through \code{\link{set_email}}.}

\item{format}{the format from which to download the data. By default it is NULL for \code{import_*} functions and tries to read 'stata', 'spss' and 'sas' in the specific order. This can be useful if some countries don't have a particular format available.  Alternatively, the user can specify the format which can either be 'stata', 'spss' or 'sas'. For the \code{download_*} functions it is set to 'stata' because the format should be specified before the downloading. When using \code{import_country} the data will be downloaded and read in the \code{format} specified. For \code{download_country}, the data is downloaded from the specified \code{format} (only 'spss' and 'stata' supported, see details).}

\item{output_dir}{a character vector with the output directory in case you want to
only download the files using \code{download_country}. Defaults to your working
directory. This will be interpreted as a \strong{directory} and not a path with
a file name.}
}
\value{
for \code{import_country} if \code{length(rounds)} is 1, it returns a tibble
with the latest version of that round. Otherwise it returns a list of \code{length(rounds)}
containing the latest version of each round. For \code{download_country}, if 
\code{output_dir} is a valid directory, it returns the saved directories invisibly
and saves all the rounds in the chosen \code{format} in \code{output_dir}
}
\description{
Download integrated rounds separately for countries from the European Social
Survey
}
\details{
Use \code{import_country} to download specified rounds for a given country and
import them to R.
\code{import_all_cntrounds} will download all rounds for a given country by default
and \code{download_country} will download rounds and save them in a specified
\code{format} in the supplied directory.

The \code{format} argument from \code{import_country} should not matter to the user
because the data is read into R either way. However, different formats might have
different handling of the encoding of some questions. This option was preserved
so that the user
can switch between formats if any encoding errors are found in the data. For more
details see the discussion \href{https://github.com/ropensci/essurvey/issues/11}{here}.
For this particular argument, 'sas' is not supported because the data formats have
changed between ESS waves and separate formats require different functions to be
read. To preserve parsimony and format errors between waves, the user should use
'spss' or 'stata'.
}
\examples{
\dontrun{

set_email("your_email@email.com")

# Get first three rounds for Denmark
dk_three <- import_country("Denmark", 1:3)

# Only download the files, this will return nothing

temp_dir <- tempdir()

download_country(
 "Turkey",
 rounds = c(2, 4),
 output_dir = temp_dir
)

# By default, download_country downloads 'stata' files but
# you can also download 'spss' or 'sas' files.

download_country(
 "Turkey",
 rounds = c(2, 4),
 output_dir = temp_dir,
 format = 'spss'
)

# If email is not registered at ESS website, error will arise
uk_one <- import_country("United Kingdom", 5, "wrong_email@email.com")
# Error in authenticate(ess_email) : 
# The email address you provided is not associated with any registered user.
# Create an account at http://www.europeansocialsurvey.org/user/new

# If selected rounds don't exist, error will arise

czech_two <- import_country("Czech Republic", c(1, 22))

# Error in country_url(country, rounds) : 
# Only rounds ESS1, ESS2, ESS4, ESS5, ESS6, ESS7, ESS8 available
# for Czech Republic
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recode_missings.R
\name{recode_missings}
\alias{recode_missings}
\alias{recode_numeric_missing}
\alias{recode_strings_missing}
\title{Recode pre-defined missing values as NA}
\usage{
recode_missings(ess_data, missing_codes)

recode_numeric_missing(x, missing_codes)

recode_strings_missing(y, missing_codes)
}
\arguments{
\item{ess_data}{data frame or \code{\link[tibble]{tibble}} with data from the
European Social Survey. This data frame should come either
from \code{\link{import_rounds}}, \code{\link{import_country}} or read with
\code{\link[haven]{read_dta}} or \code{\link[haven]{read_spss}}. This is the case because it
identifies missing values using \code{\link[haven]{labelled}} classes.}

\item{missing_codes}{a character vector with values 'Not applicable',
'Refusal', 'Don't Know', 'No answer' or 'Not available'. By default
all values are chosen. Note that the wording is case sensitive.}

\item{x}{a \code{\link[haven]{labelled}} numeric}

\item{y}{a character vector}
}
\value{
The same data frame or \code{\link[tibble]{tibble}} but with values 'Not applicable',
'Refusal', 'Don't Know', 'No answer' and 'Not available' recoded
as NA.
}
\description{
This function is not needed any more, please see the details section.
}
\details{
Data from the European Social Survey is always accompanied by a script
that recodes the categories 'Not applicable', 'Refusal', 'Don't Know',
'No answer' and 'Not available' to missing. This function recodes
these categories to NA



The European Social Survey now provides these values recoded automatically
in Stata data files. These missing categories are now read as missing values
by \code{\link[haven]{read_dta}}, reading the missing categories correctly from Stata.For an example on how these values are coded, see \href{https://github.com/ropensci/essurvey/issues/35}{here}.

Old details:

When downloading data directly from the European Social Survey's website,
the downloaded .zip file contains a script that recodes some categories
as missings in Stata and SPSS formats. 

For recoding numeric variables \code{recode_numeric_missings}
uses the labels provided by the \code{\link[haven]{labelled}}
class to delete the labels matched in \code{missing_codes}. For the
character variables matching is done with the underlying number assigned to
each category, namely 6, 7, 8, 9 and 9 for 'Not applicable', Refusal',
'Don't Know', No answer' and 'Not available'.

The functions are a direct translation of the Stata script that comes
along when downloading one of the rounds. The Stata script is the same
for all rounds and all countries, meaning that these functions work
for all rounds.
}
\examples{
\dontrun{
seven <- import_rounds(7, your_email)

attr(seven$tvtot, "labels")
mean(seven$tvtot, na.rm = TRUE)

names(table(seven$lnghom1))
# First three are actually missing values

seven_recoded <- recode_missings(seven)

attr(seven_recoded$tvtot, "labels")
# All missings have been removed
mean(seven_recoded$tvtot, na.rm = TRUE)

names(table(seven_recoded$lnghom1))
# All missings have been removed

# If you want to operate on specific variables
# you can use other recode_*_missing 

seven$tvtot <- recode_numeric_missing(seven$tvtot)

# Recode only 'Don't know' and 'No answer' to missing
seven$tvpol <- recode_numeric_missing(seven$tvpol, c("Don't know", "No answer"))


# The same can be done with recode_strings_missing
}

}
