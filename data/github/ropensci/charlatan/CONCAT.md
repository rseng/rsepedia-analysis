charlatan
=========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/charlatan/workflows/R-check/badge.svg)](https://github.com/ropensci/charlatan/actions?query=workflow%3AR-check)
[![cran checks](https://cranchecks.info/badges/worst/charlatan)](https://cranchecks.info/pkgs/charlatan)
[![codecov](https://codecov.io/gh/ropensci/charlatan/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/charlatan)
[![cran version](https://www.r-pkg.org/badges/version/charlatan)](https://cran.r-project.org/package=charlatan)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/charlatan)](https://github.com/r-hub/cranlogs.app)
[![](https://badges.ropensci.org/94_status.svg)](https://github.com/ropensci/onboarding/issues/94)

`charlatan` makes fake data, inspired from and borrowing some code from Python's faker (https://github.com/joke2k/faker)

Make fake data for:

* person names
* jobs
* phone numbers
* colors: names, hex, rgb
* credit cards
* DOIs
* numbers in range and from distributions
* gene sequences
* geographic coordinates
* emails
* URIs, URLs, and their parts
* IP addresses
* more coming ...

Possible use cases for `charlatan`:

* Students in a classroom setting learning any task that needs a dataset.
* People doing simulations/modeling that need some fake data
* Generate fake dataset of users for a database before actual users exist
* Complete missing spots in a dataset
* Generate fake data to replace sensitive real data with before public release
* Create a random set of colors for visualization
* Generate random coordinates for a map
* Get a set of randomly generated DOIs (Digital Object Identifiers) to
assign to fake scholarly artifacts
* Generate fake taxonomic names for a biological dataset
* Get a set of fake sequences to use to test code/software that uses
sequence data

Reasons to use `charlatan`:

* Lite weight, few dependencies
* Relatively comprehensive types of data, and more being added
* Comprehensive set of languages supported, more being added
* Useful R features such as creating entire fake data.frame's

## Installation

cran version


```r
install.packages("charlatan")
```

dev version


```r
remotes::install_github("ropensci/charlatan")
```


```r
library("charlatan")
```

## high level function

... for all fake data operations


```r
x <- fraudster()
x$job()
#> [1] "Toxicologist"
x$name()
#> [1] "Bart Franecki"
x$color_name()
#> [1] "IndianRed"
```

## locale support

Adding more locales through time, e.g.,

Locale support for job data


```r
ch_job(locale = "en_US", n = 3)
#> [1] "Ranger/warden"       "Psychotherapist"     "Immigration officer"
ch_job(locale = "fr_FR", n = 3)
#> [1] "Géotechnicien"                               
#> [2] "Professeur documentaliste"                   
#> [3] "Ingénieur efficacité énergétique du bâtiment"
ch_job(locale = "hr_HR", n = 3)
#> [1] "Policajac"                           "Voditelj projekta"                  
#> [3] "Zdravstveno laboratorijski tehničar"
ch_job(locale = "uk_UA", n = 3)
#> [1] "Фотограф" "Зоолог"   "Мірошник"
ch_job(locale = "zh_TW", n = 3)
#> [1] "CNC電腦程式編排人員" "特用化學工程師"      "財務或會計主管"
```

For colors:


```r
ch_color_name(locale = "en_US", n = 3)
#> [1] "DarkSlateGray" "Indigo"        "NavajoWhite"
ch_color_name(locale = "uk_UA", n = 3)
#> [1] "Червоно-буро-помаранчевий" "Темно-лососевий"          
#> [3] "Блідо-брунатний"
```

More coming soon ...

## generate a dataset


```r
ch_generate()
#> # A tibble: 10 x 3
#>    name                     job                        phone_number      
#>    <chr>                    <chr>                      <chr>             
#>  1 Mr. Posey Stehr III      Immigration officer        +61(2)7879379341  
#>  2 Ms. Henriette Wiegand    Catering manager           1-580-580-8638x830
#>  3 Irena Russel             Retail banker              +04(7)9699546042  
#>  4 Dr. Daniel Bechtelar DDS Architectural technologist 1-834-397-4529x863
#>  5 Dr. Kasey Davis          Designer, jewellery        351.022.9534x24105
#>  6 London Hansen-Hackett    Graphic designer           +06(5)1147537086  
#>  7 Lilyana Runte            Counsellor                 01692508550       
#>  8 Shaquana Herzog          Theme park manager         667.617.8036x99553
#>  9 Maybell Raynor-Hartmann  Writer                     (616)978-2091     
#> 10 Averie Murphy            Community pharmacist       1-111-441-1704
```


```r
ch_generate('job', 'phone_number', n = 30)
#> # A tibble: 30 x 2
#>    job                                         phone_number      
#>    <chr>                                       <chr>             
#>  1 Armed forces training and education officer 1-673-556-2393x997
#>  2 Soil scientist                              1-296-630-3970    
#>  3 Optician, dispensing                        1-678-990-8871    
#>  4 Learning disability nurse                   461.171.6544      
#>  5 Editor, commissioning                       05011328685       
#>  6 Designer, exhibition/display                +26(6)2762788230  
#>  7 Financial risk analyst                      1-636-012-0957x508
#>  8 Scientist, biomedical                       719.524.4489      
#>  9 Teacher, English as a foreign language      +54(0)1232453568  
#> 10 Lecturer, higher education                  (853)580-9291x3186
#> # … with 20 more rows
```


## person name


```r
ch_name()
#> [1] "Kara Boehm"
```


```r
ch_name(10)
#>  [1] "Rebecca Monahan"        "Suzann Franecki"        "Debby Nikolaus"        
#>  [4] "Ama Ullrich"            "Arba Volkman"           "Antony Mueller"        
#>  [7] "Ms. Cinnamon Anderson"  "Iver Hermann"           "Shirleen Mills-Schmidt"
#> [10] "Hadley Little"
```


## phone number


```r
ch_phone_number()
#> [1] "+36(0)2342842531"
```


```r
ch_phone_number(10)
#>  [1] "08296463291"        "970.366.6818"       "01055866557"       
#>  [4] "01717878683"        "785-103-9978"       "1-079-787-2377x619"
#>  [7] "323.362.8212"       "1-303-274-5722"     "493.066.7885x8181" 
#> [10] "610.791.1645x3705"
```

## job


```r
ch_job()
#> [1] "Therapeutic radiographer"
```


```r
ch_job(10)
#>  [1] "Environmental manager"               "Designer, blown glass/stained glass"
#>  [3] "Conservator, furniture"              "Copy"                               
#>  [5] "Administrator, local government"     "Investment analyst"                 
#>  [7] "Public librarian"                    "Engineer, materials"                
#>  [9] "Mechanical engineer"                 "Forest/woodland manager"
```

## credit cards


```r
ch_credit_card_provider()
#> [1] "VISA 16 digit"
ch_credit_card_provider(n = 4)
#> [1] "VISA 16 digit"    "JCB 15 digit"     "JCB 15 digit"     "American Express"
```


```r
ch_credit_card_number()
#> [1] "561223593016571"
ch_credit_card_number(n = 10)
#>  [1] "54998053024724596"   "869968125239286630"  "210063772612064392" 
#>  [4] "4060155369087233"    "501898051709842"     "3712676203745602"   
#>  [7] "3461064670166497"    "3096517555374787348" "3158434698000233509"
#> [10] "3037311974396594"
```


```r
ch_credit_card_security_code()
#> [1] "811"
ch_credit_card_security_code(10)
#>  [1] "598"  "164"  "0297" "083"  "741"  "519"  "948"  "452"  "6641" "286"
```

## Usage in the wild

- eacton/R-Utility-Belt-ggplot2 (https://github.com/eacton/R-Utility-Belt-ggplot2/blob/836a6bd303fbfde4a334d351e0d1c63f71c4ec68/furry_dataset.R)


## Contributors

* Scott Chamberlain (https://github.com/sckott)
* Kyle Voytovich (https://github.com/kylevoyto)
* Martin Pedersen (https://github.com/MartinMSPedersen)

## similar art

* wakefield (https://github.com/trinker/wakefield)
* ids (https://github.com/richfitz/ids)
* rcorpora (https://github.com/gaborcsardi/rcorpora)
* synthpop (https://cran.r-project.org/package=synthpop)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/charlatan/issues).
* License: MIT
* Get citation information for `charlatan` in R doing `citation(package = 'charlatan')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
charlatan 0.4.0
===============

### NEW FEATURES

* gains new vignette "Contributing to charlatan" - given that it can be complicated to contribute, this vignette should make the process easier (#49) (#84)
* `InternetProvider` gains new method `slug` (#67)
* `MiscProvider` gains two new methods `boolean` and `null_boolean` (#70)
* `es_PE` locale support added to `PhoneNumberProvider` (#108)
* `en_NZ` locale support added to `AddressProvider`, `InternetProvider`, and `PersonProvider` (#109)
* main vignette gains examples on using the `MissingDataProvider` thanks to @KKulma (#110)
* `PhoneNumberProvider` gains support for locales: `dk_DK`, `en_NZ`, `id_ID`, `th_TH`, and `tw_GH` (#100)
* each R6 provider gains new method `allowed_locales()` - the exported character vector of allowed locales for each provider has moved inside of the R6 class in `$private` because there's no reason for the user to modify allowed locales - `allowed_locales()` reads this vector for each provider

### MINOR IMPROVEMENTS

* convert all documentation to the new R6 support in roxygen2


charlatan 0.3.0
===============

### NEW FEATURES

* `ch_job()` and `JobsProvider` gains `da_DK` locale support (#94) from @MartinMSPedersen

### MINOR IMPROVEMENTS

* fixes for `PersonProvider` for locale `fr_FR`: fix accents; avoid awkward french names; now can do double first names; removed some duplicate names   (#35) (#83) from @kylevoyto
* remove leading and trailing whitespace in `JobsProvider` and `PersonProvider` where found; and remove some blank suffixes for `fa_IR` `PersonProvider` (#88) (#91) from @kylevoyto
* standardization of locale names to always be `xx_XX` where first two letters are lowercase and second two are uppercase (#90) from @kylevoyto
* change locale for Danish/Denmark from `dk_DK` to `da_DK` to comply with ISO-3166 (#93) from @MartinMSPedersen
* fix Danish phone number formats to match phone numbers actually used there (#93) from @MartinMSPedersen
* remove duplicates and sort names across `PersonProvider` for various locales (#96) from @MartinMSPedersen
* mention similar packages (#72)


charlatan 0.2.2
===============

### BUG FIXES

* run examples conditionally if packages installed for packages in Suggests: `iptools` and `stringi` (#82)


charlatan 0.2.0
===============

### NEW FEATURES

* new package author: <https://github.com/kylevoyto>
* gains `ElementsProvider` and associated methods `ch_element_element()` and `ch_element_symbol()` for getting element names and symbols (#55)
* gains `InternetProvider` with many methods, including for domain names, urls (and their parts), emails, tld's, etc. (#66)
* gains `MiscProvider` with methods for getting locale names and locale codes  (#69)
* gains `UserAgentProvider` for user agent strings (#57)
* gains `FileProvider` with methods for mime type, file extension, file names and paths (#59)
* gains `LoremProvider` with methods for words, sentences and paragraphs (#58)
* `JobProvider` gains Finnish locale (#79)

### MINOR IMPROVEMENTS

* mention usage in the wild in README (#54)
* change behavior when a locale doesn't have a data type from erroring to a zero length string (#64)
* switch to markdown docs (#68)
* fix `PersonProvider` for locale `en_GB` - we were ignoring probabilities of different names (#63) (#75)
* fix `ColorProvider`: generate only the 216 colors in safe web colors (https://en.wikipedia.org/wiki/Web_colors#Web-safe_colors) - and fix method for generating hex colors (#18) (#42) (#76)
* fix to have `safe_color_name` within `ColorProvider` be sensitive to locale (#17) (#77)
* packages `stringi` and `iptools` moved from Imports to Suggests - not required for package use now unless a few specific methods used (#71)
* `AddressProvider` gains methods `street_name`, `street_address`, `postcode`, and `address`. in addition, various fixes to `AddressProvider`  (#62) (#80)


charlatan 0.1.0
===============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 3.6.2 Patched
* ubuntu 14.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

I have checked the 1 reverse dependency, and there were no problems.
See (<https://github.com/ropensci/charlatan/tree/master/revdep>).

---

This version adds support for a variety of new language locales to various data type methods.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/charlatan/issues)

### Code contributions

Check out the [Contributing to charlatan](https://github.com/ropensci/charlatan/blob/master/vignettes/contributing.Rmd) vignette for details if you're adding a new provider or locale.

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/charlatan.git`
* Make sure to track progress upstream (i.e., on our version of `charlatan` at `ropensci/charlatan`) by doing `git remote add upstream https://github.com/ropensci/charlatan.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/charlatan`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 3.6.2 Patched (2019-12-12 r77564) |
|os       |macOS Mojave 10.14.6                        |
|system   |x86_64, darwin15.6.0                        |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-01-23                                  |

# Dependencies

|package   |old   |new   |Δ  |
|:---------|:-----|:-----|:--|
|charlatan |0.3.0 |0.4.0 |*  |
|rlang     |NA    |0.4.3 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*charlatan
=========

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/charlatan/workflows/R-check/badge.svg)](https://github.com/ropensci/charlatan/actions?query=workflow%3AR-check)
[![cran checks](https://cranchecks.info/badges/worst/charlatan)](https://cranchecks.info/pkgs/charlatan)
[![codecov](https://codecov.io/gh/ropensci/charlatan/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/charlatan)
[![cran version](https://www.r-pkg.org/badges/version/charlatan)](https://cran.r-project.org/package=charlatan)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/charlatan)](https://github.com/r-hub/cranlogs.app)
[![](https://badges.ropensci.org/94_status.svg)](https://github.com/ropensci/onboarding/issues/94)

`charlatan` makes fake data, inspired from and borrowing some code from Python's faker (https://github.com/joke2k/faker)

Make fake data for:

* person names
* jobs
* phone numbers
* colors: names, hex, rgb
* credit cards
* DOIs
* numbers in range and from distributions
* gene sequences
* geographic coordinates
* emails
* URIs, URLs, and their parts
* IP addresses
* more coming ...

Possible use cases for `charlatan`:

* Students in a classroom setting learning any task that needs a dataset.
* People doing simulations/modeling that need some fake data
* Generate fake dataset of users for a database before actual users exist
* Complete missing spots in a dataset
* Generate fake data to replace sensitive real data with before public release
* Create a random set of colors for visualization
* Generate random coordinates for a map
* Get a set of randomly generated DOIs (Digital Object Identifiers) to
assign to fake scholarly artifacts
* Generate fake taxonomic names for a biological dataset
* Get a set of fake sequences to use to test code/software that uses
sequence data

Reasons to use `charlatan`:

* Lite weight, few dependencies
* Relatively comprehensive types of data, and more being added
* Comprehensive set of languages supported, more being added
* Useful R features such as creating entire fake data.frame's

## Installation

cran version

```{r eval=FALSE}
install.packages("charlatan")
```

dev version

```{r eval=FALSE}
remotes::install_github("ropensci/charlatan")
```

```{r}
library("charlatan")
```

## high level function

... for all fake data operations

```{r}
x <- fraudster()
x$job()
x$name()
x$color_name()
```

## locale support

Adding more locales through time, e.g.,

Locale support for job data

```{r}
ch_job(locale = "en_US", n = 3)
ch_job(locale = "fr_FR", n = 3)
ch_job(locale = "hr_HR", n = 3)
ch_job(locale = "uk_UA", n = 3)
ch_job(locale = "zh_TW", n = 3)
```

For colors:

```{r}
ch_color_name(locale = "en_US", n = 3)
ch_color_name(locale = "uk_UA", n = 3)
```

More coming soon ...

## generate a dataset

```{r}
ch_generate()
```

```{r}
ch_generate('job', 'phone_number', n = 30)
```


## person name

```{r}
ch_name()
```

```{r}
ch_name(10)
```


## phone number

```{r}
ch_phone_number()
```

```{r}
ch_phone_number(10)
```

## job

```{r}
ch_job()
```

```{r}
ch_job(10)
```

## credit cards

```{r}
ch_credit_card_provider()
ch_credit_card_provider(n = 4)
```

```{r}
ch_credit_card_number()
ch_credit_card_number(n = 10)
```

```{r}
ch_credit_card_security_code()
ch_credit_card_security_code(10)
```

## Usage in the wild

- eacton/R-Utility-Belt-ggplot2 (https://github.com/eacton/R-Utility-Belt-ggplot2/blob/836a6bd303fbfde4a334d351e0d1c63f71c4ec68/furry_dataset.R)


## Contributors

* Scott Chamberlain (https://github.com/sckott)
* Kyle Voytovich (https://github.com/kylevoyto)
* Martin Pedersen (https://github.com/MartinMSPedersen)

## similar art

* wakefield (https://github.com/trinker/wakefield)
* ids (https://github.com/richfitz/ids)
* rcorpora (https://github.com/gaborcsardi/rcorpora)
* synthpop (https://cran.r-project.org/package=synthpop)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/charlatan/issues).
* License: MIT
* Get citation information for `charlatan` in R doing `citation(package = 'charlatan')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "Introduction to the charlatan package"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to the charlatan package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

`charlatan` makes fake data, inspired from and borrowing some code from Python's [faker](https://github.com/joke2k/faker)

Why would you want to make fake data? Here's some possible use cases to
give you a sense for what you can do with this package:

* Students in a classroom setting learning any task that needs a dataset.
* People doing simulations/modeling that need some fake data
* Generate fake dataset of users for a database before actual users exist
* Complete missing spots in a dataset
* Generate fake data to replace sensitive real data with before public release
* Create a random set of colors for visualization
* Generate random coordinates for a map
* Get a set of randomly generated DOIs (Digial Object Identifiers) to
assign to fake scholarly artifacts
* Generate fake taxonomic names for a biological dataset
* Get a set of fake sequences to use to test code/software that uses
sequence data

## Contributing

See the [**Contributing to charlatan**](https://docs.ropensci.org/charlatan/articles/contributing.html) vignette

## Package API

* Low level interfaces: All of these are `R6` objects that a user can
initialize and then call methods on. These contain all the logic that
the below interfaces use.
* High level interfaces: There are high level functions prefixed with
`ch_*()` that wrap low level interfaces, and are meant to be easier
to use and provide an easy way to make many instances of a thing.
* `ch_generate()` - generate a data.frame with fake data, choosing which columns to include from the data types provided in `charlatan`
* `fraudster()` - single interface to all fake data methods, - returns
vectors/lists of data - this function wraps the `ch_*()` functions described above

## Install

Stable version from CRAN

```{r eval=FALSE}
install.packages("charlatan")
```

Development version from Github

```{r eval=FALSE}
devtools::install_github("ropensci/charlatan")
```

```{r}
library("charlatan")
```

## high level function

... for all fake data operations

```{r}
x <- fraudster()
x$job()
x$name()
x$job()
x$color_name()
```

## locale support

Adding more locales through time, e.g.,

Locale support for job data

```{r}
ch_job(locale = "en_US", n = 3)
ch_job(locale = "fr_FR", n = 3)
ch_job(locale = "hr_HR", n = 3)
ch_job(locale = "uk_UA", n = 3)
ch_job(locale = "zh_TW", n = 3)
```

For colors:

```{r}
ch_color_name(locale = "en_US", n = 3)
ch_color_name(locale = "uk_UA", n = 3)
```

More coming soon ...

## generate a dataset

```{r}
ch_generate()
```

```{r}
ch_generate('job', 'phone_number', n = 30)
```

## Data types

### person name

```{r}
ch_name()
```

```{r}
ch_name(10)
```


### phone number

```{r}
ch_phone_number()
```

```{r}
ch_phone_number(10)
```

### job

```{r}
ch_job()
```

```{r}
ch_job(10)
```

### credit cards

```{r}
ch_credit_card_provider()
ch_credit_card_provider(n = 4)
```

```{r}
ch_credit_card_number()
ch_credit_card_number(n = 10)
```

```{r}
ch_credit_card_security_code()
ch_credit_card_security_code(10)
```

## Missing data 
`charlatan` makes it very easy to generate fake data with missing entries. First, you need to run `MissingDataProvider()` and then make an appropriate `make_missing()` call specifying the data type to be generated. This method picks a random number (`N`) of slots in the input `make_missing` vector and then picks `N` random positions that will be replaced with NA matching the input class.

```{r}
testVector <- MissingDataProvider$new()
```

### character strings

```{r}
testVector$make_missing(x = ch_generate()$name) 
```


### numeric data

```{r}
testVector$make_missing(x = ch_integer(10)) 
```


### logicals

```{r}
set.seed(123)
testVector$make_missing(x = sample(c(TRUE, FALSE), 10, replace = TRUE)) 
```

## Messy data

Real data is messy, right?  `charlatan` makes it easy to create
messy data. This is still in the early stages so is not available
across most data types and languages, but we're working on it.

For example, create messy names:

```{r}
ch_name(50, messy = TRUE)
```

Right now only suffixes and prefixes for names in `en_US` locale
are supported. Notice above some variation in prefixes and suffixes.
---
title: "Contributing to charlatan"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Contributing to charlatan}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

`charlatan` is a wee bit complex. This vignette aims to help you contribute
to the package. For a general introduction on contributing to rOpenSci packages
see our [Contributing guide](https://devguide.ropensci.org/contributingguide.html).

Let's start with some definitions.

## Definitions

For the purposes of this package:

* **Provider**: a type of data that can be generated in `charlatan`. For example,
we have providers for phone numbers, addresses and people's names. Adding a provider
may involve a single file, more than one file; and a single R6 class or many
R6 classes.
* **Locale**: a locale for our purposes is a specific spoken language that's
associated with a specific country. You can have more than one locale for a
given language (e.g., `en-US`, `en-GB`). Some fakers won't have any locales,
whereas others can have many.

If you aren't familiar with R6, have a look at the
[R6 website](https://r6.r-lib.org/), in particular the
[introductory vignette](https://r6.r-lib.org/articles/Introduction.html).

## Communication

Open an issue if you want to add a new provider or locale to
an existing provider; it helps make sure there's no duplicated effort and
we can help make sure you have the knowledge you need.

## Adding a new provider

Providers are generally first created by making an R6 class. Let's start with a
heavily simplified base R6 class that defines some utility methods. We
call it `BaseProvider` in `charlatan`, but here we'll call it `MyBaseProvider`
to avoid confusion.

```{r}
library(R6)
MyBaseProvider <- R6::R6Class(
  'MyBaseProvider',
  public = list(
    random_element = function(x) {
      if (length(x) == 0) return('')
      if (inherits(x, "character")) if (!any(nzchar(x))) return('')
      x[sample.int(n = length(x), size = 1)]
    },

    random_int = function(min = 0, max = 9999, size = 1) {
      stopifnot(max >= min)
      num <- max - min + 1
      sample.int(n = num, size = size, replace = TRUE) + (min - 1)
    }
  )
)
```

### Providers without locale support

If you don't need to handle locales it becomes simpler:


```{r}
FooBar <- R6::R6Class(
  'FooBar',
  inherit = charlatan::BaseProvider,
  public = list(
    integer = function(n = 1, min = 1, max = 1000) {
      super$random_int(min, max, n)
    }
  )
)
```

We can create an instance of the `FooBar` class by calling `$new()` on it.
It only has one method `integer()`, which we can call to get a random
integer.

```{r}
x <- FooBar$new()
x
x$integer()
```

### Providers with locale support

If your provider will need to handle different locales, it gets a bit more
complex. In the Python library [faker][] from which this package draws
inspiration, you can create separate folders for each provider within the
Python library.

However, R doesn't allow this, so instead we categorize different locales for
each provider within the file names. For example, for the address provider we
have files in the package:

- [address-provider.R](https://github.com/ropensci/charlatan/blob/master/R/address-provider.R)
- [address-provider-en_US.R](https://github.com/ropensci/charlatan/blob/master/R/address-provider-en_US.R)
- [address-provider-en_GB.R](https://github.com/ropensci/charlatan/blob/master/R/address-provider-en_GB.R)

Where the latter two provides specific data for each locale, and the first
file has the `AddressProvider` class that pulls in the locale specific data.

Here, we'll create a very simplified `AddressProvider` class using an
example locale file.

```{r}
library(charlatan)
file <- system.file("examples", "address-provider-en_US.R", package = "charlatan")
source(file)
MyAddressProvider <- R6::R6Class(
  inherit = MyBaseProvider,
  'MyAddressProvider',
  lock_objects = FALSE,
  public = list(
    locale = NULL,
    city_suffixes = NULL,

    initialize = function() {
      self$locale <- 'en_us'
      self$city_suffixes <-
        eval(parse(text = paste0("city_suffixes_", self$locale)))
    },

    city_suffix = function() {
      super$random_element(self$city_suffixes)
    }
  )
)
```

We can create an instance of the `MyAddressProvider` class by calling `$new()` on it.
It only has one method `city_suffix()`, which we can call to get a random
city suffix.

```{r}
x <- MyAddressProvider$new()
x
x$city_suffix()
```

#### Adding a new locale

When you want to add a new locale to an existing provider, look in the `R/` folder
of the package and the locales that are available are in the file names.

Pick one of the locale files for the provider you're extending, make a duplicate of it
and rename the file with your new locale. Then modify the duplicate, copying the
format but putting in place the appropriate information for the new locale.

Where the data comes from for the new locale may vary. One easy way to start may be
porting over locales in the [faker][] Python library that are not yet in `charlatan`.

If it's a locale for which you can't easily port over from another library, you need
to get the data from a variety of sources. There are some R based packages that
should help:

- [WikidataR][]
- [humaniformat][]
- [WikidataQueryServiceR][]

Keep in mind when using data to look at their license, if any, and any implications
with respect to whether it can be used in this package.


## How locale specific data are used in providers

It's a little tricky how this is done. In the `initialize()` block of each main provider
file (e.g., `address-provider.R`) we pull in the appropriate locale specific data
based on the user input locale. For example, here's an abbreviated `initialize` block from
the `AddressProvider`:

```r
initialize = function(locale = NULL) {
  if (!is.null(locale)) {
    # check global locales
    super$check_locale(locale)
    # check address provider locales
    check_locale_(locale, address_provider_locales)
    self$locale <- locale
  } else {
    self$locale <- 'en_US'
  }

  self$city_prefixes <- parse_eval("city_prefixes_", self$locale)
}
```

A few things to note:

* if no locale is given, we by default use `en_US`
* we check that the locale given is in the allowed set
* for each data type, here just city prefixes shown, use the internal function
`parse_eval()` to pull in the data. Essentially, `parse_eval()` makes the
string `city_prefixes_en_US`, then finds that in the package environment
and `eval()`'s it to bring the data into the R6 object in the `city_prefixes`
slot. We repeat this for each data type. The result is the user initialized
class with locale specific data.


[faker]: https://github.com/joke2k/faker
[humaniformat]: https://github.com/Ironholds/humaniformat
[WikidataQueryServiceR]: https://cran.r-project.org/package=WikidataQueryServiceR
[WikidataR]: https://cran.r-project.org/package=WikidataR
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinate-provider.R
\name{CoordinateProvider}
\alias{CoordinateProvider}
\title{CoordinateProvider}
\description{
coordinates methods
}
\examples{
z <- CoordinateProvider$new()
z$lon()
z$lat()
z$position()
z$position(bbox = c(-120, 30, -110, 60))
}
\keyword{internal}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-lon}{\code{CoordinateProvider$lon()}}
\item \href{#method-lat}{\code{CoordinateProvider$lat()}}
\item \href{#method-position}{\code{CoordinateProvider$position()}}
\item \href{#method-clone}{\code{CoordinateProvider$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-lon"></a>}}
\if{latex}{\out{\hypertarget{method-lon}{}}}
\subsection{Method \code{lon()}}{
a latitude value
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoordinateProvider$lon()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-lat"></a>}}
\if{latex}{\out{\hypertarget{method-lat}{}}}
\subsection{Method \code{lat()}}{
a longitude value
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoordinateProvider$lat()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-position"></a>}}
\if{latex}{\out{\hypertarget{method-position}{}}}
\subsection{Method \code{position()}}{
a position, of form \verb{[longitude,latitude]}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoordinateProvider$position(bbox = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{bbox}}{optionally, specify a bounding box for the position
to be in, of the form \verb{[west,south,east,north]} - checks that the
bbox has valid values for lat and long}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoordinateProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missing.R
\name{ch_missing}
\alias{ch_missing}
\title{Create missing data}
\usage{
ch_missing(x, n = 1)
}
\arguments{
\item{x}{Input vector, can be any class - only 1 vetor}

\item{n}{(integer) number of things to get, any non-negative integer}
}
\description{
Create missing data
}
\examples{
ch_missing(letters)
ch_missing(letters, 10)
ch_missing(letters, 20)
}
\seealso{
\link{MissingDataProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/company.R
\name{ch_company}
\alias{ch_company}
\title{Create fake company names and other company bits}
\usage{
ch_company(n = 1, locale = NULL)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}

\item{locale}{(character) the locale to use. See
\code{CompanyProvider$new()$allowed_locales()} for locales supported.}
}
\description{
Create fake company names and other company bits
}
\examples{
ch_company()
ch_company(10)
ch_company(500)

ch_company(locale = "fr_FR", n = 10)
ch_company(locale = "cs_CZ", n = 10)
ch_company(locale = "es_MX", n = 10)
ch_company(locale = "hr_HR", n = 10)
}
\seealso{
\link{CompanyProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phone_number.R
\name{ch_phone_number}
\alias{ch_phone_number}
\title{Create fake phone numbers}
\usage{
ch_phone_number(n = 1, locale = NULL)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}

\item{locale}{(character) the locale to use. See
\code{PhoneNumberProvider$new()$allowed_locales()} for locales
supported (default: en_US)}
}
\description{
Create fake phone numbers
}
\examples{
ch_phone_number()
ch_phone_number(10)
ch_phone_number(500)

# locales
ch_phone_number(locale = "fr_FR")
ch_phone_number(locale = "uk_UA")
ch_phone_number(locale = "en_CA")
ch_phone_number(locale = "lv_LV")
}
\seealso{
\link{PhoneNumberProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/address-provider.R
\name{AddressProvider}
\alias{AddressProvider}
\title{AddressProvider}
\description{
address methods
}
\examples{
(z <- AddressProvider$new())
z$locale
z$allowed_locales()
z$city_suffix()
z$street_suffix()
z$building_number()
z$city()
z$country()
z$street_name()
z$street_address()
z$address()
z$country()
z$country_code()
z$postcode()

# en_GB
(z <- AddressProvider$new('en_GB'))
z$locale
z$locale_data
z$locale_data$postcode_sets
z$postcode
z$postcode()
z$street_name()

# en_NZ
(z <- AddressProvider$new('en_NZ'))
z$locale
z$street_name()

# es_ES
(z <- AddressProvider$new('es_ES'))
z$locale
z$street_name()

# nl_NL
(z <- AddressProvider$new('nl_NL'))
z$locale
z$street_name()
z$postcode()
z$city()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{AddressProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) xxx}

\item{\code{city_prefixes}}{(character) xxx}

\item{\code{city_suffixes}}{(character) xxx}

\item{\code{street_suffixes}}{(character) xxx}

\item{\code{building_number_formats}}{(character) xxx}

\item{\code{postcode_formats}}{(character) xxx}

\item{\code{states}}{(character) xxx}

\item{\code{states_abbr}}{(character) xxx}

\item{\code{military_state_abbr}}{(character) xxx}

\item{\code{military_ship_prefix}}{(character) xxx}

\item{\code{military_apo_format}}{(character) xxx}

\item{\code{military_dpo_format}}{(character) xxx}

\item{\code{city_formats}}{(character) xxx}

\item{\code{street_name_formats}}{(character) xxx}

\item{\code{street_address_formats}}{(character) xxx}

\item{\code{address_formats}}{(character) xxx}

\item{\code{secondary_address_formats}}{(character) xxx}

\item{\code{countries}}{(character) xxx}

\item{\code{country_codes}}{(character) xxx}

\item{\code{locale_data}}{(character) xxx}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-allowed_locales}{\code{AddressProvider$allowed_locales()}}
\item \href{#method-new}{\code{AddressProvider$new()}}
\item \href{#method-city_suffix}{\code{AddressProvider$city_suffix()}}
\item \href{#method-street_suffix}{\code{AddressProvider$street_suffix()}}
\item \href{#method-building_number}{\code{AddressProvider$building_number()}}
\item \href{#method-city}{\code{AddressProvider$city()}}
\item \href{#method-street_name}{\code{AddressProvider$street_name()}}
\item \href{#method-street_address}{\code{AddressProvider$street_address()}}
\item \href{#method-postcode}{\code{AddressProvider$postcode()}}
\item \href{#method-address}{\code{AddressProvider$address()}}
\item \href{#method-country}{\code{AddressProvider$country()}}
\item \href{#method-country_code}{\code{AddressProvider$country_code()}}
\item \href{#method-clone}{\code{AddressProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-allowed_locales"></a>}}
\if{latex}{\out{\hypertarget{method-allowed_locales}{}}}
\subsection{Method \code{allowed_locales()}}{
fetch the allowed locales for this provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$allowed_locales()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{AddressProvider} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$new(locale = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. See
\verb{$allowed_locales()} for locales supported (default: en_US)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{AddressProvider} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-city_suffix"></a>}}
\if{latex}{\out{\hypertarget{method-city_suffix}{}}}
\subsection{Method \code{city_suffix()}}{
city suffix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$city_suffix()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-street_suffix"></a>}}
\if{latex}{\out{\hypertarget{method-street_suffix}{}}}
\subsection{Method \code{street_suffix()}}{
street suffix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$street_suffix()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-building_number"></a>}}
\if{latex}{\out{\hypertarget{method-building_number}{}}}
\subsection{Method \code{building_number()}}{
building numeber
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$building_number()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-city"></a>}}
\if{latex}{\out{\hypertarget{method-city}{}}}
\subsection{Method \code{city()}}{
city
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$city()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-street_name"></a>}}
\if{latex}{\out{\hypertarget{method-street_name}{}}}
\subsection{Method \code{street_name()}}{
street name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$street_name()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-street_address"></a>}}
\if{latex}{\out{\hypertarget{method-street_address}{}}}
\subsection{Method \code{street_address()}}{
street address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$street_address()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-postcode"></a>}}
\if{latex}{\out{\hypertarget{method-postcode}{}}}
\subsection{Method \code{postcode()}}{
postal code
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$postcode()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-address"></a>}}
\if{latex}{\out{\hypertarget{method-address}{}}}
\subsection{Method \code{address()}}{
address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$address()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-country"></a>}}
\if{latex}{\out{\hypertarget{method-country}{}}}
\subsection{Method \code{country()}}{
country name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$country()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-country_code"></a>}}
\if{latex}{\out{\hypertarget{method-country_code}{}}}
\subsection{Method \code{country_code()}}{
country code
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$country_code()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AddressProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elements.R
\name{elements}
\alias{elements}
\alias{ch_element_symbol}
\alias{ch_element_element}
\title{Get elements}
\usage{
ch_element_symbol(n = 1)

ch_element_element(n = 1)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}
}
\description{
Get elements
}
\examples{
ch_element_symbol()
ch_element_symbol(10)
ch_element_symbol(50)

ch_element_element()
ch_element_element(10)
ch_element_element(50)
}
\seealso{
\link{ElementProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc-provider.R
\name{MiscProvider}
\alias{MiscProvider}
\title{MiscProvider}
\description{
miscellaneous methods
}
\examples{
(x <- MiscProvider$new())
x$language_locale_codes
x$language_code()
x$locale()
x$boolean()
x$null_boolean()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{MiscProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{language_locale_codes}}{(list) locale codes by locale family}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-boolean}{\code{MiscProvider$boolean()}}
\item \href{#method-null_boolean}{\code{MiscProvider$null_boolean()}}
\item \href{#method-locale}{\code{MiscProvider$locale()}}
\item \href{#method-language_code}{\code{MiscProvider$language_code()}}
\item \href{#method-clone}{\code{MiscProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-boolean"></a>}}
\if{latex}{\out{\hypertarget{method-boolean}{}}}
\subsection{Method \code{boolean()}}{
get a random boolean, \code{TRUE} or \code{FALSE}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MiscProvider$boolean(chance_of_getting_true = 50)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{chance_of_getting_true}}{(integer) an integer, default: 50}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-null_boolean"></a>}}
\if{latex}{\out{\hypertarget{method-null_boolean}{}}}
\subsection{Method \code{null_boolean()}}{
get a random boolean, \code{TRUE} or \code{FALSE}, or NULL
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MiscProvider$null_boolean()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-locale"></a>}}
\if{latex}{\out{\hypertarget{method-locale}{}}}
\subsection{Method \code{locale()}}{
get a random locale
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MiscProvider$locale()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-language_code"></a>}}
\if{latex}{\out{\hypertarget{method-language_code}{}}}
\subsection{Method \code{language_code()}}{
random language code
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MiscProvider$language_code()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MiscProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doi.R
\name{ch_doi}
\alias{ch_doi}
\title{Create fake DOIs (Digital Object Identifiers)}
\usage{
ch_doi(n = 1)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}
}
\description{
Create fake DOIs (Digital Object Identifiers)
}
\examples{
ch_doi()
ch_doi(10)
ch_doi(100)
}
\seealso{
\link{DOIProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate.R
\name{ch_generate}
\alias{ch_generate}
\title{Generate a fake dataset}
\usage{
ch_generate(..., n = 10, locale = NULL)
}
\arguments{
\item{...}{columns to include. must be in the allowed set. See
\strong{Allowed column names} below. Three default columns are included
(name, job, phone_number) if nothing is specified - but are overridden
by any input.}

\item{n}{(integer) number of things to get, any non-negative integer}

\item{locale}{(character) the locale to use. options: only supported
for data types that have locale support, See each data provider for
details.}
}
\description{
Generate a fake dataset
}
\section{Allowed column names}{

\itemize{
\item name (default included)
\item job (default included)
\item phone_number (default included)
\item currency
\item color_name
\item rgb_color
\item rgb_css_color
}
}

\examples{
ch_generate()
ch_generate(n = 1)
ch_generate(n = 100)

ch_generate('job')
ch_generate('job', 'name')
ch_generate('job', 'color_name')

# locale
ch_generate(locale = "en_US")
ch_generate(locale = "fr_FR")
ch_generate(locale = "fr_CH")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/person-provider.R
\name{PersonProvider}
\alias{PersonProvider}
\title{PersonProvider}
\description{
person names methods
}
\details{
Note that with the male/female versions if the locale
doesn't provide a male/female version then we fall back to the
generic thing, e.g., if no female first name we give you first
name
}
\examples{
x <- PersonProvider$new()
x$locale
x$render()
x$first_name()
x$first_name_female()
x$first_name_male()
x$last_name()
x$last_name_female()
x$last_name_male()

x <- PersonProvider$new(locale = "en_GB")
x$locale
x$render()
x$first_name()
x$first_name_female()
x$first_name_male()
x$last_name()
x$last_name_female()
x$last_name_male()

z <- PersonProvider$new(locale = "fr_FR")
z$locale
z$render()
z$first_name()
z$first_name_female()
z$first_name_male()
z$last_name()
z$last_name_female()
z$last_name_male()
z$prefix()

z <- PersonProvider$new(locale = "de_AT")
z$locale
z$render()
z$first_name()
z$last_name()
z$prefix()

z <- PersonProvider$new(locale = "cs_CZ")
z$locale
z$render()
z$first_name()
z$first_name_female()
z$first_name_male()
z$last_name()
z$last_name_female()
z$last_name_male()
z$prefix()

z <- PersonProvider$new(locale = "es_MX")
z$locale
z$render()
z$first_name()
z$first_name_female()
z$first_name_male()
z$last_name()
z$prefix()

z <- PersonProvider$new(locale = "en_NZ")
z$locale
z$render()
z$first_name()
z$first_name_female()
z$first_name_male()
z$last_name()

PersonProvider$new(locale = "fr_CH")$render()
PersonProvider$new(locale = "fi_FI")$render()
PersonProvider$new(locale = "fa_IR")$render()
PersonProvider$new(locale = "es_ES")$render()
PersonProvider$new(locale = "de_DE")$render()
PersonProvider$new(locale = "de_AT")$render()
PersonProvider$new(locale = "cs_CZ")$render()
PersonProvider$new(locale = "bg_BG")$render()
PersonProvider$new(locale = "da_DK")$render()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{PersonProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) the locale}

\item{\code{formats}}{(character) person name formats}

\item{\code{person}}{(character) person name data}

\item{\code{messy}}{(logical) the messy setting, \code{TRUE} or \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-allowed_locales}{\code{PersonProvider$allowed_locales()}}
\item \href{#method-new}{\code{PersonProvider$new()}}
\item \href{#method-render}{\code{PersonProvider$render()}}
\item \href{#method-first_name}{\code{PersonProvider$first_name()}}
\item \href{#method-first_name_female}{\code{PersonProvider$first_name_female()}}
\item \href{#method-first_name_male}{\code{PersonProvider$first_name_male()}}
\item \href{#method-last_name}{\code{PersonProvider$last_name()}}
\item \href{#method-last_name_female}{\code{PersonProvider$last_name_female()}}
\item \href{#method-last_name_male}{\code{PersonProvider$last_name_male()}}
\item \href{#method-prefix}{\code{PersonProvider$prefix()}}
\item \href{#method-prefix_female}{\code{PersonProvider$prefix_female()}}
\item \href{#method-prefix_male}{\code{PersonProvider$prefix_male()}}
\item \href{#method-suffix}{\code{PersonProvider$suffix()}}
\item \href{#method-suffix_female}{\code{PersonProvider$suffix_female()}}
\item \href{#method-suffix_male}{\code{PersonProvider$suffix_male()}}
\item \href{#method-clone}{\code{PersonProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-allowed_locales"></a>}}
\if{latex}{\out{\hypertarget{method-allowed_locales}{}}}
\subsection{Method \code{allowed_locales()}}{
fetch the allowed locales for this provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$allowed_locales()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{PersonProvider} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$new(locale = NULL, messy = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. See
\verb{$allowed_locales()} for locales supported (default: en_US)}

\item{\code{messy}}{(logical) make some messy data. Default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{PersonProvider} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-render"></a>}}
\if{latex}{\out{\hypertarget{method-render}{}}}
\subsection{Method \code{render()}}{
Make a person's name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$render(fmt = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{fmt}}{(character) a name format, default: \code{NULL}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-first_name"></a>}}
\if{latex}{\out{\hypertarget{method-first_name}{}}}
\subsection{Method \code{first_name()}}{
make a first name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$first_name()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-first_name_female"></a>}}
\if{latex}{\out{\hypertarget{method-first_name_female}{}}}
\subsection{Method \code{first_name_female()}}{
make a female first name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$first_name_female()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-first_name_male"></a>}}
\if{latex}{\out{\hypertarget{method-first_name_male}{}}}
\subsection{Method \code{first_name_male()}}{
make a male first name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$first_name_male()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-last_name"></a>}}
\if{latex}{\out{\hypertarget{method-last_name}{}}}
\subsection{Method \code{last_name()}}{
make a last name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$last_name()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-last_name_female"></a>}}
\if{latex}{\out{\hypertarget{method-last_name_female}{}}}
\subsection{Method \code{last_name_female()}}{
make a female last name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$last_name_female()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-last_name_male"></a>}}
\if{latex}{\out{\hypertarget{method-last_name_male}{}}}
\subsection{Method \code{last_name_male()}}{
make a male last name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$last_name_male()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-prefix"></a>}}
\if{latex}{\out{\hypertarget{method-prefix}{}}}
\subsection{Method \code{prefix()}}{
make a name prefix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$prefix()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-prefix_female"></a>}}
\if{latex}{\out{\hypertarget{method-prefix_female}{}}}
\subsection{Method \code{prefix_female()}}{
make a female name prefix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$prefix_female()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-prefix_male"></a>}}
\if{latex}{\out{\hypertarget{method-prefix_male}{}}}
\subsection{Method \code{prefix_male()}}{
make a male name prefix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$prefix_male()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-suffix"></a>}}
\if{latex}{\out{\hypertarget{method-suffix}{}}}
\subsection{Method \code{suffix()}}{
make a name suffix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$suffix()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-suffix_female"></a>}}
\if{latex}{\out{\hypertarget{method-suffix_female}{}}}
\subsection{Method \code{suffix_female()}}{
make a female name suffix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$suffix_female()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-suffix_male"></a>}}
\if{latex}{\out{\hypertarget{method-suffix_male}{}}}
\subsection{Method \code{suffix_male()}}{
make a male name suffix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$suffix_male()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PersonProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/currency-provider.R
\name{CurrencyProvider}
\alias{CurrencyProvider}
\title{CurrencyProvider}
\description{
currencies
}
\examples{
z <- CurrencyProvider$new()
z$render()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{CurrencyProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{formats}}{(character) currency formats character vector}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-render}{\code{CurrencyProvider$render()}}
\item \href{#method-clone}{\code{CurrencyProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-render"></a>}}
\if{latex}{\out{\hypertarget{method-render}{}}}
\subsection{Method \code{render()}}{
get a currency
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CurrencyProvider$render()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(character string) of length 3
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CurrencyProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credit_card-provider.R
\name{CreditCardProvider}
\alias{CreditCardProvider}
\title{CreditCardProvider}
\description{
credit card methods
}
\examples{
z <- CreditCardProvider$new()
z$credit_card_provider()
z$credit_card_number()
z$credit_card_security_code()
z$generate_number(13)
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{CreditCardProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{luhn_lookup}}{(list) luhn lookup, named list}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-credit_card_type}{\code{CreditCardProvider$credit_card_type()}}
\item \href{#method-generate_number}{\code{CreditCardProvider$generate_number()}}
\item \href{#method-credit_card_provider}{\code{CreditCardProvider$credit_card_provider()}}
\item \href{#method-credit_card_number}{\code{CreditCardProvider$credit_card_number()}}
\item \href{#method-credit_card_security_code}{\code{CreditCardProvider$credit_card_security_code()}}
\item \href{#method-clone}{\code{CreditCardProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-credit_card_type"></a>}}
\if{latex}{\out{\hypertarget{method-credit_card_type}{}}}
\subsection{Method \code{credit_card_type()}}{
Returns a random credit card type
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CreditCardProvider$credit_card_type(card_type = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{card_type}}{(character) a card type, see \code{credit_card_types}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-generate_number"></a>}}
\if{latex}{\out{\hypertarget{method-generate_number}{}}}
\subsection{Method \code{generate_number()}}{
make a credit card number with specific starting numbers
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CreditCardProvider$generate_number(prefix, length = 13)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{prefix}}{the start of the CC number as a string, any number of digits.}

\item{\code{length}}{the length of the CC number to generate. Typically 13 or 16}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-credit_card_provider"></a>}}
\if{latex}{\out{\hypertarget{method-credit_card_provider}{}}}
\subsection{Method \code{credit_card_provider()}}{
credit card provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CreditCardProvider$credit_card_provider(card_type = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{card_type}}{(character) a card type, see \code{credit_card_types}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-credit_card_number"></a>}}
\if{latex}{\out{\hypertarget{method-credit_card_number}{}}}
\subsection{Method \code{credit_card_number()}}{
credit card number
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CreditCardProvider$credit_card_number(card_type = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{card_type}}{(character) a card type, see \code{credit_card_types}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-credit_card_security_code"></a>}}
\if{latex}{\out{\hypertarget{method-credit_card_security_code}{}}}
\subsection{Method \code{credit_card_security_code()}}{
credit card security code
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CreditCardProvider$credit_card_security_code(card_type = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{card_type}}{(character) a card type, see \code{credit_card_types}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CreditCardProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequences.R
\name{ch_gene_sequence}
\alias{ch_gene_sequence}
\title{Create fake gene sequences}
\usage{
ch_gene_sequence(n = 1, length = 30)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}

\item{length}{(integer) length of sequence to create}
}
\description{
Create fake gene sequences
}
\examples{
ch_gene_sequence()
ch_gene_sequence(10)
ch_gene_sequence(100)

ch_gene_sequence(length = 500)
ch_gene_sequence(10, length = 500)
}
\seealso{
\link{SequenceProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/useragent-provider.R
\name{UserAgentProvider}
\alias{UserAgentProvider}
\title{UserAgentProvider}
\description{
user agent methods
}
\examples{
(x <- UserAgentProvider$new())
x$locale
x$mac_processor()
x$linux_processor()
x$user_agent()
x$chrome()
x$firefox()
x$internet_explorer()
x$opera()
x$safari()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{UserAgentProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) the locale}

\item{\code{user_agents}}{(character) user agent browser specific strings}

\item{\code{windows_platform_tokens}}{(character) windows platform tokens}

\item{\code{linux_processors}}{(character) linux processor options}

\item{\code{mac_processors}}{(character) mac processor options}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-allowed_locales}{\code{UserAgentProvider$allowed_locales()}}
\item \href{#method-new}{\code{UserAgentProvider$new()}}
\item \href{#method-mac_processor}{\code{UserAgentProvider$mac_processor()}}
\item \href{#method-linux_processor}{\code{UserAgentProvider$linux_processor()}}
\item \href{#method-user_agent}{\code{UserAgentProvider$user_agent()}}
\item \href{#method-chrome}{\code{UserAgentProvider$chrome()}}
\item \href{#method-firefox}{\code{UserAgentProvider$firefox()}}
\item \href{#method-safari}{\code{UserAgentProvider$safari()}}
\item \href{#method-opera}{\code{UserAgentProvider$opera()}}
\item \href{#method-internet_explorer}{\code{UserAgentProvider$internet_explorer()}}
\item \href{#method-windows_platform_token}{\code{UserAgentProvider$windows_platform_token()}}
\item \href{#method-linux_platform_token}{\code{UserAgentProvider$linux_platform_token()}}
\item \href{#method-mac_platform_token}{\code{UserAgentProvider$mac_platform_token()}}
\item \href{#method-clone}{\code{UserAgentProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-allowed_locales"></a>}}
\if{latex}{\out{\hypertarget{method-allowed_locales}{}}}
\subsection{Method \code{allowed_locales()}}{
fetch the allowed locales for this provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$allowed_locales()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{UserAgentProvider} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$new(locale = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. See
\verb{$allowed_locales()} for locales supported (default: en_US)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{UserAgentProvider} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mac_processor"></a>}}
\if{latex}{\out{\hypertarget{method-mac_processor}{}}}
\subsection{Method \code{mac_processor()}}{
a mac processor
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$mac_processor()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-linux_processor"></a>}}
\if{latex}{\out{\hypertarget{method-linux_processor}{}}}
\subsection{Method \code{linux_processor()}}{
a linux processor
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$linux_processor()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-user_agent"></a>}}
\if{latex}{\out{\hypertarget{method-user_agent}{}}}
\subsection{Method \code{user_agent()}}{
a random user agent string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$user_agent()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-chrome"></a>}}
\if{latex}{\out{\hypertarget{method-chrome}{}}}
\subsection{Method \code{chrome()}}{
a chrome user agent string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$chrome(
  version_from = 13,
  version_to = 63,
  build_from = 800,
  build_to = 899
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{version_from}}{(integer) minimum version}

\item{\code{version_to}}{(integer) maximum version}

\item{\code{build_from}}{(integer) minimum build}

\item{\code{build_to}}{(integer) maximum build}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-firefox"></a>}}
\if{latex}{\out{\hypertarget{method-firefox}{}}}
\subsection{Method \code{firefox()}}{
a firefox user agent string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$firefox()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-safari"></a>}}
\if{latex}{\out{\hypertarget{method-safari}{}}}
\subsection{Method \code{safari()}}{
a safari user agent string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$safari()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-opera"></a>}}
\if{latex}{\out{\hypertarget{method-opera}{}}}
\subsection{Method \code{opera()}}{
an opera user agent string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$opera()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-internet_explorer"></a>}}
\if{latex}{\out{\hypertarget{method-internet_explorer}{}}}
\subsection{Method \code{internet_explorer()}}{
an internet explorer user agent string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$internet_explorer()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-windows_platform_token"></a>}}
\if{latex}{\out{\hypertarget{method-windows_platform_token}{}}}
\subsection{Method \code{windows_platform_token()}}{
a windows platform token
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$windows_platform_token()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-linux_platform_token"></a>}}
\if{latex}{\out{\hypertarget{method-linux_platform_token}{}}}
\subsection{Method \code{linux_platform_token()}}{
a linux platform token
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$linux_platform_token()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mac_platform_token"></a>}}
\if{latex}{\out{\hypertarget{method-mac_platform_token}{}}}
\subsection{Method \code{mac_platform_token()}}{
a mac platform token
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$mac_platform_token()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UserAgentProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fraudster.R
\name{fraudster}
\alias{fraudster}
\title{Fraudster - catch all client to make all types of fake data}
\usage{
fraudster(locale = NULL)
}
\arguments{
\item{locale}{(character) the locale to use. options: en_US (default),
fr_FR, fr_CH, hr_FR, fa_IR, pl_PL, ru_RU, uk_UA, zh_TW.}
}
\description{
Fraudster - catch all client to make all types of fake data
}
\examples{
# English - the default locale
(x <- fraudster())
x$job()
x$name()
x$color_name()
x$safe_color_name()
x$hex_color()
x$safe_hex_color()
x$rgb_color()
x$rgb_css_color()

# different locales
## French
(y <- fraudster(locale = "fr_FR"))
y$job()

## Croatian
(z <- fraudster(locale = "hr_HR"))
z$job()

## Ukranian
(w <- fraudster(locale = "uk_UA"))
w$job()
w$color_name()

# geospatial
x$lat()
x$lon()
x$position()

# DOIs (Digital Object Identifier)
x$doi()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phonenumbers-provider.R
\name{PhoneNumberProvider}
\alias{PhoneNumberProvider}
\title{PhoneNumberProvider}
\description{
methods for generating phone numbers
}
\examples{
z <- PhoneNumberProvider$new()
z$render()

PhoneNumberProvider$new(locale = "fr_FR")$render()
PhoneNumberProvider$new(locale = "sk_SK")$render()

# locales with area codes
PhoneNumberProvider$new(locale = "en_AU")$render()
PhoneNumberProvider$new(locale = "en_NZ")$render()
PhoneNumberProvider$new(locale = "es_PE")$render()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{PhoneNumberProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) the locale}

\item{\code{formats}}{phone number formats}

\item{\code{area_code_formats}}{area code formats}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-allowed_locales}{\code{PhoneNumberProvider$allowed_locales()}}
\item \href{#method-new}{\code{PhoneNumberProvider$new()}}
\item \href{#method-render}{\code{PhoneNumberProvider$render()}}
\item \href{#method-clone}{\code{PhoneNumberProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-allowed_locales"></a>}}
\if{latex}{\out{\hypertarget{method-allowed_locales}{}}}
\subsection{Method \code{allowed_locales()}}{
fetch the allowed locales for this provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PhoneNumberProvider$allowed_locales()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{PhoneNumberProvider} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PhoneNumberProvider$new(locale = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. See
\verb{$allowed_locales()} for locales supported (default: en_US)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{PhoneNumberProvider} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-render"></a>}}
\if{latex}{\out{\hypertarget{method-render}{}}}
\subsection{Method \code{render()}}{
Make a phone number
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PhoneNumberProvider$render()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PhoneNumberProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/company-provider.R
\name{CompanyProvider}
\alias{CompanyProvider}
\title{CompanyProvider}
\description{
company name/etc. methods
}
\examples{
x <- CompanyProvider$new()
x$locale
x$company()
x$company_suffix()
x$catch_phrase()
x$bs()

x <- CompanyProvider$new(locale = "fr_FR")
x$locale
x$company()
x$company_suffix()
x$siren()

x <- CompanyProvider$new(locale = "hr_HR")
x$locale
x$company()
x$company_suffix()

x <- CompanyProvider$new(locale = "it_IT")
x$locale
x$company()
x$company_suffix()
x$bs()

CompanyProvider$new(locale = "es_MX")$bs()
CompanyProvider$new(locale = "es_MX")$company_prefix()
CompanyProvider$new(locale = "es_MX")$catch_phrase()

CompanyProvider$new(locale = "bg_BG")$company()
CompanyProvider$new(locale = "cs_CZ")$company()
CompanyProvider$new(locale = "de_DE")$company()
CompanyProvider$new(locale = "fa_IR")$company()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{CompanyProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) xxx}

\item{\code{formats}}{(character) xxx}

\item{\code{prefixes}}{(character) xxx}

\item{\code{suffixes}}{(character) xxx}

\item{\code{catch_phrase_words}}{(character) xxx}

\item{\code{bsWords}}{(character) xxx}

\item{\code{siren_format}}{(character) xxx}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-allowed_locales}{\code{CompanyProvider$allowed_locales()}}
\item \href{#method-new}{\code{CompanyProvider$new()}}
\item \href{#method-company}{\code{CompanyProvider$company()}}
\item \href{#method-company_prefix}{\code{CompanyProvider$company_prefix()}}
\item \href{#method-company_suffix}{\code{CompanyProvider$company_suffix()}}
\item \href{#method-catch_phrase}{\code{CompanyProvider$catch_phrase()}}
\item \href{#method-bs}{\code{CompanyProvider$bs()}}
\item \href{#method-siren}{\code{CompanyProvider$siren()}}
\item \href{#method-clone}{\code{CompanyProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-allowed_locales"></a>}}
\if{latex}{\out{\hypertarget{method-allowed_locales}{}}}
\subsection{Method \code{allowed_locales()}}{
fetch the allowed locales for this provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CompanyProvider$allowed_locales()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{CompanyProvider} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CompanyProvider$new(locale = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. See
\verb{$allowed_locales()} for locales supported (default: en_US)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{CompanyProvider} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-company"></a>}}
\if{latex}{\out{\hypertarget{method-company}{}}}
\subsection{Method \code{company()}}{
a company name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CompanyProvider$company()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-company_prefix"></a>}}
\if{latex}{\out{\hypertarget{method-company_prefix}{}}}
\subsection{Method \code{company_prefix()}}{
a company prefix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CompanyProvider$company_prefix()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-company_suffix"></a>}}
\if{latex}{\out{\hypertarget{method-company_suffix}{}}}
\subsection{Method \code{company_suffix()}}{
a company suffix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CompanyProvider$company_suffix()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-catch_phrase"></a>}}
\if{latex}{\out{\hypertarget{method-catch_phrase}{}}}
\subsection{Method \code{catch_phrase()}}{
a catch phrase
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CompanyProvider$catch_phrase()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-bs"></a>}}
\if{latex}{\out{\hypertarget{method-bs}{}}}
\subsection{Method \code{bs()}}{
BS words
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CompanyProvider$bs()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-siren"></a>}}
\if{latex}{\out{\hypertarget{method-siren}{}}}
\subsection{Method \code{siren()}}{
a siren
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CompanyProvider$siren()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CompanyProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/numerics.R
\name{numerics}
\alias{numerics}
\alias{ch_double}
\alias{ch_integer}
\alias{ch_unif}
\alias{ch_norm}
\alias{ch_lnorm}
\alias{ch_beta}
\title{Create numbers}
\usage{
ch_double(n = 1, mean = 0, sd = 1)

ch_integer(n = 1, min = 1, max = 1000)

ch_unif(n = 1, min = 0, max = 9999)

ch_norm(n = 1, mean = 0, sd = 1)

ch_lnorm(n = 1, mean = 0, sd = 1)

ch_beta(n = 1, shape1, shape2, ncp = 0)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}

\item{mean}{mean value}

\item{sd}{standard deviation}

\item{min}{minimum value}

\item{max}{maximum value}

\item{shape1, shape2}{non-negative parameters of the Beta distribution}

\item{ncp}{non-centrality parameter}
}
\description{
Create numbers
}
\examples{
ch_double()
ch_double(10)
ch_double(100)

ch_integer()
ch_integer(10)
ch_integer(100)

ch_unif()
ch_norm()
ch_lnorm()
ch_beta(shape1 = 1, shape2 = 1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/color-provider.R
\name{ColorProvider}
\alias{ColorProvider}
\title{ColorProvider}
\description{
methods for colors
}
\examples{
x <- ColorProvider$new()
x$locale
x$color_name()
x$safe_color_name()
x$hex_color()
x$safe_hex_color()
x$rgb_color()
x$rgb_css_color()

x <- ColorProvider$new(locale = "uk_UA")
x$locale
x$color_name()
x$safe_color_name()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{ColorProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) xxx}

\item{\code{all_colors}}{(character) xxx}

\item{\code{safe_colors}}{(character) xxx}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-allowed_locales}{\code{ColorProvider$allowed_locales()}}
\item \href{#method-new}{\code{ColorProvider$new()}}
\item \href{#method-color_name}{\code{ColorProvider$color_name()}}
\item \href{#method-safe_color_name}{\code{ColorProvider$safe_color_name()}}
\item \href{#method-hex_color}{\code{ColorProvider$hex_color()}}
\item \href{#method-safe_hex_color}{\code{ColorProvider$safe_hex_color()}}
\item \href{#method-rgb_color}{\code{ColorProvider$rgb_color()}}
\item \href{#method-rgb_css_color}{\code{ColorProvider$rgb_css_color()}}
\item \href{#method-clone}{\code{ColorProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-allowed_locales"></a>}}
\if{latex}{\out{\hypertarget{method-allowed_locales}{}}}
\subsection{Method \code{allowed_locales()}}{
fetch the allowed locales for this provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ColorProvider$allowed_locales()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{ColorProvider} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ColorProvider$new(locale = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. See
\verb{$allowed_locales()} for locales supported (default: en_US)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{ColorProvider} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-color_name"></a>}}
\if{latex}{\out{\hypertarget{method-color_name}{}}}
\subsection{Method \code{color_name()}}{
color name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ColorProvider$color_name()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-safe_color_name"></a>}}
\if{latex}{\out{\hypertarget{method-safe_color_name}{}}}
\subsection{Method \code{safe_color_name()}}{
safe color name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ColorProvider$safe_color_name()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-hex_color"></a>}}
\if{latex}{\out{\hypertarget{method-hex_color}{}}}
\subsection{Method \code{hex_color()}}{
hex color
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ColorProvider$hex_color()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-safe_hex_color"></a>}}
\if{latex}{\out{\hypertarget{method-safe_hex_color}{}}}
\subsection{Method \code{safe_hex_color()}}{
safe hex color
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ColorProvider$safe_hex_color()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-rgb_color"></a>}}
\if{latex}{\out{\hypertarget{method-rgb_color}{}}}
\subsection{Method \code{rgb_color()}}{
RGB color
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ColorProvider$rgb_color()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-rgb_css_color"></a>}}
\if{latex}{\out{\hypertarget{method-rgb_css_color}{}}}
\subsection{Method \code{rgb_css_color()}}{
RGB CSS color
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ColorProvider$rgb_css_color()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ColorProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy-provider.R
\name{TaxonomyProvider}
\alias{TaxonomyProvider}
\title{TaxonomyProvider}
\description{
Taxonomy provider
}
\section{Names}{

Names were taken from Theplantlist. 500 genera names and 500
epithets were chosen at random from the set of 10,000 names in the
dataset in the \code{taxize} package. Theplantlist is, as it says on the
tin, composed of plant names - so these fake names are derived from
plant names if that matters to you. These may generate names that match
those of real taxa, but may not as well.
}

\section{Taxonomic authority}{

Randomly, the taxonomic authority is in parentheses - which represents
that the given authority was not the original authority.
}

\examples{
(z <- TaxonomyProvider$new())
z$genus()
z$epithet()
z$species()
z$species(authority = TRUE)
## FIXME - datetimeprovider slow - may be related to unix time problem
# z$species(authority = TRUE, date = TRUE)
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{TaxonomyProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{genera}}{(character) vector of generic names}

\item{\code{epithets}}{(character) vector of eptithet names}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-genus}{\code{TaxonomyProvider$genus()}}
\item \href{#method-epithet}{\code{TaxonomyProvider$epithet()}}
\item \href{#method-species}{\code{TaxonomyProvider$species()}}
\item \href{#method-clone}{\code{TaxonomyProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-genus"></a>}}
\if{latex}{\out{\hypertarget{method-genus}{}}}
\subsection{Method \code{genus()}}{
Get a genus name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TaxonomyProvider$genus()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-epithet"></a>}}
\if{latex}{\out{\hypertarget{method-epithet}{}}}
\subsection{Method \code{epithet()}}{
Get an epithet name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TaxonomyProvider$epithet()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-species"></a>}}
\if{latex}{\out{\hypertarget{method-species}{}}}
\subsection{Method \code{species()}}{
Get a binomial name (genus + epithet)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TaxonomyProvider$species(authority = FALSE, date = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{authority}}{Include authority. default: \code{FALSE}}

\item{\code{date}}{Include authority date. If \code{authority = FALSE},
this is ignored. default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TaxonomyProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{taxonomy}
\alias{taxonomy}
\alias{ch_taxonomic_genus}
\alias{ch_taxonomic_epithet}
\alias{ch_taxonomic_species}
\title{Create fake taxonomic names}
\usage{
ch_taxonomic_genus(n = 1)

ch_taxonomic_epithet(n = 1)

ch_taxonomic_species(n = 1)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}
}
\description{
Create fake taxonomic names
}
\section{Names}{

Names were taken from Theplantlist. 500 genera names and 500
epithets were chosen at random from the set of 10,000 names in the
dataset in the \code{taxize} package. Theplantlist is, as it says on the
tin, composed of plant names - so these fake names are derived from
plant names if that matters to you. These may generate names that match
those of real taxa, but may not as well.
}

\section{Taxonomic authority}{

Randomly, the taxonomic authority is in parentheses - which represents
that the given authority was not the original authority.
}

\examples{
ch_taxonomic_genus()
ch_taxonomic_genus(10)
ch_taxonomic_genus(500)

ch_taxonomic_epithet()
ch_taxonomic_epithet(10)
ch_taxonomic_epithet(500)

ch_taxonomic_species()
ch_taxonomic_species(10)
ch_taxonomic_species(500)
}
\seealso{
\link{TaxonomyProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internet-provider.R
\name{InternetProvider}
\alias{InternetProvider}
\title{InternetProvider}
\description{
internet methods, e.g., email addresses, domain names
}
\note{
Note that if a locale you set doesn't have a locale specific set
of data for \link{PersonProvider} or \link{CompanyProvider} we fall back to
\code{en_US}
}
\examples{
(x <- InternetProvider$new())
x$locale

# uri/url/tld/etc.
x$tld()
x$slug()
x$domain_word()
x$domain_name()
x$domain_name(levels = 2)
x$domain_name(levels = 3)
x$domain_name(levels = 10)
## url's
x$url()
x$url(schemes = c('hbbp', 'hggp'))
x$image_url()
## uri's
x$uri()
x$uri_page()
x$uri_extension()
x$uri_path()
x$uri_path(deep = 1)
x$uri_path(deep = 2)
x$uri_path(deep = 3)
x$uri_path(deep = 4)

# user name
x$user_name()

# emails
x$email()
x$safe_email()
x$free_email()
x$company_email()
x$free_email_domain()
x$ascii_email()
x$ascii_safe_email()
x$ascii_free_email()
x$ascii_company_email()

# addresses, mac, ipv4, ipv6
x$mac_address()
if (requireNamespace("ipaddress", quietly=TRUE)) {
  x$ipv4()
  x$ipv4(network = TRUE)
  x$ipv6()
  x$ipv6(network = TRUE)
}

# different locales
(x <- InternetProvider$new(locale = "en_AU"))
x$locale
x$tld()
x$email()
x$free_email_domain()

(x <- InternetProvider$new(locale = "de_DE"))
x$locale
x$tld()
x$uri()
x$email()
x$ascii_email()

(x <- InternetProvider$new(locale = "bg_BG"))
x$locale
x$tld()
x$uri()
x$url()
x$user_name()
x$email()
x$ascii_email()

(x <- InternetProvider$new(locale = "cs_CZ"))
x$url()
x$user_name()
x$email()

(x <- InternetProvider$new(locale = "en_NZ"))
x$free_email_domain()
x$tld()

(x <- InternetProvider$new(locale = "fa_IR"))
x$url()

(x <- InternetProvider$new(locale = "fr_FR"))
x$url()
x$user_name()
x$email()

(x <- InternetProvider$new(locale = "hr_HR"))
x$url()
x$user_name()
x$email()

# convert a string to ascii with stringi pkg
if (requireNamespace("stringi", quietly=TRUE)) {
  x$to_ascii("anï")
}
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{InternetProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) the locale}

\item{\code{safe_email_tlds}}{(character) safe email tlds}

\item{\code{free_email_domains}}{(character) free email tlds}

\item{\code{tlds}}{(character) tlds}

\item{\code{uri_pages}}{(character) uri pages}

\item{\code{uri_paths}}{(character) uri paths}

\item{\code{uri_extensions}}{(character) uri extensions}

\item{\code{user_name_formats}}{(character) user name formats}

\item{\code{email_formats}}{(character) email formats}

\item{\code{url_formats}}{(character) url formats}

\item{\code{uri_formats}}{(character) uri formats}

\item{\code{image_placeholder_services}}{(character) image uri formats}

\item{\code{replacements}}{(list) a list}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-allowed_locales}{\code{InternetProvider$allowed_locales()}}
\item \href{#method-new}{\code{InternetProvider$new()}}
\item \href{#method-to_ascii}{\code{InternetProvider$to_ascii()}}
\item \href{#method-email}{\code{InternetProvider$email()}}
\item \href{#method-safe_email}{\code{InternetProvider$safe_email()}}
\item \href{#method-free_email}{\code{InternetProvider$free_email()}}
\item \href{#method-company_email}{\code{InternetProvider$company_email()}}
\item \href{#method-ascii_email}{\code{InternetProvider$ascii_email()}}
\item \href{#method-ascii_safe_email}{\code{InternetProvider$ascii_safe_email()}}
\item \href{#method-ascii_free_email}{\code{InternetProvider$ascii_free_email()}}
\item \href{#method-ascii_company_email}{\code{InternetProvider$ascii_company_email()}}
\item \href{#method-user_name}{\code{InternetProvider$user_name()}}
\item \href{#method-tld}{\code{InternetProvider$tld()}}
\item \href{#method-free_email_domain}{\code{InternetProvider$free_email_domain()}}
\item \href{#method-url}{\code{InternetProvider$url()}}
\item \href{#method-domain_name}{\code{InternetProvider$domain_name()}}
\item \href{#method-domain_word}{\code{InternetProvider$domain_word()}}
\item \href{#method-ipv4}{\code{InternetProvider$ipv4()}}
\item \href{#method-ipv6}{\code{InternetProvider$ipv6()}}
\item \href{#method-mac_address}{\code{InternetProvider$mac_address()}}
\item \href{#method-uri_page}{\code{InternetProvider$uri_page()}}
\item \href{#method-uri_path}{\code{InternetProvider$uri_path()}}
\item \href{#method-uri_extension}{\code{InternetProvider$uri_extension()}}
\item \href{#method-uri}{\code{InternetProvider$uri()}}
\item \href{#method-slug}{\code{InternetProvider$slug()}}
\item \href{#method-image_url}{\code{InternetProvider$image_url()}}
\item \href{#method-clone}{\code{InternetProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-allowed_locales"></a>}}
\if{latex}{\out{\hypertarget{method-allowed_locales}{}}}
\subsection{Method \code{allowed_locales()}}{
fetch the allowed locales for this provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$allowed_locales()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{InternetProvider} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$new(locale = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. See
\verb{$allowed_locales()} for locales supported (default: en_US)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{InternetProvider} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_ascii"></a>}}
\if{latex}{\out{\hypertarget{method-to_ascii}{}}}
\subsection{Method \code{to_ascii()}}{
convert to ascii
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$to_ascii(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{the stringn to convert to ascii}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-email"></a>}}
\if{latex}{\out{\hypertarget{method-email}{}}}
\subsection{Method \code{email()}}{
get an email address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$email(domain = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{domain}}{(character) a domain name, if not given, a random
name is chosen}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-safe_email"></a>}}
\if{latex}{\out{\hypertarget{method-safe_email}{}}}
\subsection{Method \code{safe_email()}}{
get a safe email address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$safe_email()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-free_email"></a>}}
\if{latex}{\out{\hypertarget{method-free_email}{}}}
\subsection{Method \code{free_email()}}{
a free email address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$free_email()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-company_email"></a>}}
\if{latex}{\out{\hypertarget{method-company_email}{}}}
\subsection{Method \code{company_email()}}{
company email address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$company_email()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ascii_email"></a>}}
\if{latex}{\out{\hypertarget{method-ascii_email}{}}}
\subsection{Method \code{ascii_email()}}{
ascii email address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$ascii_email()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ascii_safe_email"></a>}}
\if{latex}{\out{\hypertarget{method-ascii_safe_email}{}}}
\subsection{Method \code{ascii_safe_email()}}{
safe ascii email address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$ascii_safe_email()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ascii_free_email"></a>}}
\if{latex}{\out{\hypertarget{method-ascii_free_email}{}}}
\subsection{Method \code{ascii_free_email()}}{
an ascii free email address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$ascii_free_email()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ascii_company_email"></a>}}
\if{latex}{\out{\hypertarget{method-ascii_company_email}{}}}
\subsection{Method \code{ascii_company_email()}}{
ascii company email address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$ascii_company_email()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-user_name"></a>}}
\if{latex}{\out{\hypertarget{method-user_name}{}}}
\subsection{Method \code{user_name()}}{
a user name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$user_name()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-tld"></a>}}
\if{latex}{\out{\hypertarget{method-tld}{}}}
\subsection{Method \code{tld()}}{
a tld
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$tld()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-free_email_domain"></a>}}
\if{latex}{\out{\hypertarget{method-free_email_domain}{}}}
\subsection{Method \code{free_email_domain()}}{
free email domain name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$free_email_domain()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-url"></a>}}
\if{latex}{\out{\hypertarget{method-url}{}}}
\subsection{Method \code{url()}}{
a url
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$url(schemes = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{schemes}}{(character vector) a url scheme, defaults are http and https}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-domain_name"></a>}}
\if{latex}{\out{\hypertarget{method-domain_name}{}}}
\subsection{Method \code{domain_name()}}{
Produce an Internet domain name with the specified
number of subdomain levels
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$domain_name(levels = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{levels}}{(integer) how many levels, must be >1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-domain_word"></a>}}
\if{latex}{\out{\hypertarget{method-domain_word}{}}}
\subsection{Method \code{domain_word()}}{
a domain word
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$domain_word()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ipv4"></a>}}
\if{latex}{\out{\hypertarget{method-ipv4}{}}}
\subsection{Method \code{ipv4()}}{
an ipv4 address or network
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$ipv4(network = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{network}}{(logical) produce a network}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ipv6"></a>}}
\if{latex}{\out{\hypertarget{method-ipv6}{}}}
\subsection{Method \code{ipv6()}}{
an ipv6 address or network
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$ipv6(network = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{network}}{(logical) produce a network}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mac_address"></a>}}
\if{latex}{\out{\hypertarget{method-mac_address}{}}}
\subsection{Method \code{mac_address()}}{
a mac address
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$mac_address()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-uri_page"></a>}}
\if{latex}{\out{\hypertarget{method-uri_page}{}}}
\subsection{Method \code{uri_page()}}{
a uri page
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$uri_page()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-uri_path"></a>}}
\if{latex}{\out{\hypertarget{method-uri_path}{}}}
\subsection{Method \code{uri_path()}}{
a uri path
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$uri_path(deep = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{how deep to go, an integer, if not given an integer
between 1 and 4 (inclusive) is chosen}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-uri_extension"></a>}}
\if{latex}{\out{\hypertarget{method-uri_extension}{}}}
\subsection{Method \code{uri_extension()}}{
a uri extension
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$uri_extension()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-uri"></a>}}
\if{latex}{\out{\hypertarget{method-uri}{}}}
\subsection{Method \code{uri()}}{
a uri
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$uri()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-slug"></a>}}
\if{latex}{\out{\hypertarget{method-slug}{}}}
\subsection{Method \code{slug()}}{
a slug
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$slug(value = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{value}}{(character) a string, if given, returns itself, if not, uses
\link{LoremProvider} to get a random string. default: \code{NULL}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-image_url"></a>}}
\if{latex}{\out{\hypertarget{method-image_url}{}}}
\subsection{Method \code{image_url()}}{
Returns URL to placeholder image -
Example: http://placehold.it/640x480
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$image_url(width = NULL, height = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{width}}{image width, in pixels}

\item{\code{height}}{image height, in pixels}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InternetProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/currency.R
\name{ch_currency}
\alias{ch_currency}
\title{Create fake currencies}
\usage{
ch_currency(n = 1)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}
}
\description{
Create fake currencies
}
\examples{
ch_currency()
ch_currency(10)
ch_currency(500)
}
\seealso{
\link{CurrencyProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name.R
\name{ch_name}
\alias{ch_name}
\title{Create fake person names}
\usage{
ch_name(n = 1, locale = NULL, messy = FALSE)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}

\item{locale}{(character) the locale to use. See
\code{PersonProvider$new()$allowed_locales()} for locales supported
(default: en_US)}

\item{messy}{(logical) make some messy data. Default: \code{FALSE}}
}
\description{
Create fake person names
}
\examples{
ch_name()
ch_name(10)
ch_name(500)

ch_name(locale = "fr_FR", n = 10)
ch_name(locale = "fr_CH", n = 10)
ch_name(locale = "fa_IR", n = 10)
ch_name(locale = "fi_FI", n = 10)
}
\seealso{
\link{PersonProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequence-provider.R
\name{SequenceProvider}
\alias{SequenceProvider}
\title{SequenceProvider}
\description{
genetic sequence generator
}
\examples{
z <- SequenceProvider$new()
z$render()
z$render(10)
z$render(100)
z$render(500)
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{letters}}{(character) nucleotides to include in sequences}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-render}{\code{SequenceProvider$render()}}
\item \href{#method-clone}{\code{SequenceProvider$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-render"></a>}}
\if{latex}{\out{\hypertarget{method-render}{}}}
\subsection{Method \code{render()}}{
Make a sequence
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SequenceProvider$render(length = 30)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{length}}{(integer) length of sequence to create. default: 30}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SequenceProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fraudster.R
\name{FraudsterClient}
\alias{FraudsterClient}
\title{FraudsterClient}
\description{
Fraudster R6 client
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) the locale to use}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{FraudsterClient$new()}}
\item \href{#method-print}{\code{FraudsterClient$print()}}
\item \href{#method-job}{\code{FraudsterClient$job()}}
\item \href{#method-name}{\code{FraudsterClient$name()}}
\item \href{#method-color_name}{\code{FraudsterClient$color_name()}}
\item \href{#method-safe_color_name}{\code{FraudsterClient$safe_color_name()}}
\item \href{#method-hex_color}{\code{FraudsterClient$hex_color()}}
\item \href{#method-safe_hex_color}{\code{FraudsterClient$safe_hex_color()}}
\item \href{#method-rgb_color}{\code{FraudsterClient$rgb_color()}}
\item \href{#method-rgb_css_color}{\code{FraudsterClient$rgb_css_color()}}
\item \href{#method-lat}{\code{FraudsterClient$lat()}}
\item \href{#method-lon}{\code{FraudsterClient$lon()}}
\item \href{#method-position}{\code{FraudsterClient$position()}}
\item \href{#method-doi}{\code{FraudsterClient$doi()}}
\item \href{#method-timezone}{\code{FraudsterClient$timezone()}}
\item \href{#method-unix_time}{\code{FraudsterClient$unix_time()}}
\item \href{#method-date_time}{\code{FraudsterClient$date_time()}}
\item \href{#method-genus}{\code{FraudsterClient$genus()}}
\item \href{#method-epithet}{\code{FraudsterClient$epithet()}}
\item \href{#method-species}{\code{FraudsterClient$species()}}
\item \href{#method-sequence}{\code{FraudsterClient$sequence()}}
\item \href{#method-phone_number}{\code{FraudsterClient$phone_number()}}
\item \href{#method-double}{\code{FraudsterClient$double()}}
\item \href{#method-integer}{\code{FraudsterClient$integer()}}
\item \href{#method-uniform}{\code{FraudsterClient$uniform()}}
\item \href{#method-norm}{\code{FraudsterClient$norm()}}
\item \href{#method-lnorm}{\code{FraudsterClient$lnorm()}}
\item \href{#method-beta}{\code{FraudsterClient$beta()}}
\item \href{#method-currency}{\code{FraudsterClient$currency()}}
\item \href{#method-credit_card_provider}{\code{FraudsterClient$credit_card_provider()}}
\item \href{#method-credit_card_number}{\code{FraudsterClient$credit_card_number()}}
\item \href{#method-credit_card_security_code}{\code{FraudsterClient$credit_card_security_code()}}
\item \href{#method-clone}{\code{FraudsterClient$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{FraudsterClient} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$new(locale = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. options: en_US (default),
fr_FR, fr_CH, hr_FR, fa_IR, pl_PL, ru_RU, uk_UA, zh_TW.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{RequestSignature} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the \code{FraudsterClient} class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-job"></a>}}
\if{latex}{\out{\hypertarget{method-job}{}}}
\subsection{Method \code{job()}}{
jobs
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$job(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-name"></a>}}
\if{latex}{\out{\hypertarget{method-name}{}}}
\subsection{Method \code{name()}}{
names
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$name(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-color_name"></a>}}
\if{latex}{\out{\hypertarget{method-color_name}{}}}
\subsection{Method \code{color_name()}}{
colors
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$color_name(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-safe_color_name"></a>}}
\if{latex}{\out{\hypertarget{method-safe_color_name}{}}}
\subsection{Method \code{safe_color_name()}}{
safe color name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$safe_color_name(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-hex_color"></a>}}
\if{latex}{\out{\hypertarget{method-hex_color}{}}}
\subsection{Method \code{hex_color()}}{
hex color
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$hex_color(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-safe_hex_color"></a>}}
\if{latex}{\out{\hypertarget{method-safe_hex_color}{}}}
\subsection{Method \code{safe_hex_color()}}{
safe hex color
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$safe_hex_color(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-rgb_color"></a>}}
\if{latex}{\out{\hypertarget{method-rgb_color}{}}}
\subsection{Method \code{rgb_color()}}{
rgb color
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$rgb_color(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-rgb_css_color"></a>}}
\if{latex}{\out{\hypertarget{method-rgb_css_color}{}}}
\subsection{Method \code{rgb_css_color()}}{
rgb css color
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$rgb_css_color(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-lat"></a>}}
\if{latex}{\out{\hypertarget{method-lat}{}}}
\subsection{Method \code{lat()}}{
latitude
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$lat(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-lon"></a>}}
\if{latex}{\out{\hypertarget{method-lon}{}}}
\subsection{Method \code{lon()}}{
longitude
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$lon(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-position"></a>}}
\if{latex}{\out{\hypertarget{method-position}{}}}
\subsection{Method \code{position()}}{
long/lat coordinate pair
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$position(n = 1, bbox = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}

\item{\code{bbox}}{a bounding box, see \code{\link[=ch_position]{ch_position()}}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-doi"></a>}}
\if{latex}{\out{\hypertarget{method-doi}{}}}
\subsection{Method \code{doi()}}{
DOIs
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$doi(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-timezone"></a>}}
\if{latex}{\out{\hypertarget{method-timezone}{}}}
\subsection{Method \code{timezone()}}{
date times
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$timezone(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-unix_time"></a>}}
\if{latex}{\out{\hypertarget{method-unix_time}{}}}
\subsection{Method \code{unix_time()}}{
unix time
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$unix_time(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-date_time"></a>}}
\if{latex}{\out{\hypertarget{method-date_time}{}}}
\subsection{Method \code{date_time()}}{
date time
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$date_time(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-genus"></a>}}
\if{latex}{\out{\hypertarget{method-genus}{}}}
\subsection{Method \code{genus()}}{
taxonomic genus
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$genus(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-epithet"></a>}}
\if{latex}{\out{\hypertarget{method-epithet}{}}}
\subsection{Method \code{epithet()}}{
taxonomic epithet
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$epithet(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-species"></a>}}
\if{latex}{\out{\hypertarget{method-species}{}}}
\subsection{Method \code{species()}}{
taxonomic species (genus + epithet)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$species(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-sequence"></a>}}
\if{latex}{\out{\hypertarget{method-sequence}{}}}
\subsection{Method \code{sequence()}}{
random genetic sequence
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$sequence(n = 1, length = 30)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}

\item{\code{length}}{(integer) length of the sequence. default: 30}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-phone_number"></a>}}
\if{latex}{\out{\hypertarget{method-phone_number}{}}}
\subsection{Method \code{phone_number()}}{
phone number
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$phone_number(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-double"></a>}}
\if{latex}{\out{\hypertarget{method-double}{}}}
\subsection{Method \code{double()}}{
a double
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$double(n = 1, mean = 0, sd = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}

\item{\code{mean}}{mean value, default: 0}

\item{\code{sd}}{standard deviation, default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-integer"></a>}}
\if{latex}{\out{\hypertarget{method-integer}{}}}
\subsection{Method \code{integer()}}{
an integer
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$integer(n = 1, min = 1, max = 1000)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}

\item{\code{min}}{minimum value, default: 1}

\item{\code{max}}{maximum value, default: 1000}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-uniform"></a>}}
\if{latex}{\out{\hypertarget{method-uniform}{}}}
\subsection{Method \code{uniform()}}{
an integer from a uniform distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$uniform(n = 1, min = 0, max = 9999)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}

\item{\code{min}}{minimum value, default: 0}

\item{\code{max}}{maximum value, default: 9999}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-norm"></a>}}
\if{latex}{\out{\hypertarget{method-norm}{}}}
\subsection{Method \code{norm()}}{
an integer from a normal distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$norm(n = 1, mean = 0, sd = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}

\item{\code{mean}}{mean value, default: 0}

\item{\code{sd}}{standard deviation, default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-lnorm"></a>}}
\if{latex}{\out{\hypertarget{method-lnorm}{}}}
\subsection{Method \code{lnorm()}}{
an integer from a lognormal distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$lnorm(n = 1, mean = 0, sd = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}

\item{\code{mean}}{mean value, default: 0}

\item{\code{sd}}{standard deviation, default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-beta"></a>}}
\if{latex}{\out{\hypertarget{method-beta}{}}}
\subsection{Method \code{beta()}}{
an integer from a beta distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$beta(n = 1, shape1, shape2, ncp = 0)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}

\item{\code{shape1}}{non-negative parameters of the Beta distribution}

\item{\code{shape2}}{non-negative parameters of the Beta distribution}

\item{\code{ncp}}{non-centrality parameter, default: 0}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-currency"></a>}}
\if{latex}{\out{\hypertarget{method-currency}{}}}
\subsection{Method \code{currency()}}{
currency
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$currency(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-credit_card_provider"></a>}}
\if{latex}{\out{\hypertarget{method-credit_card_provider}{}}}
\subsection{Method \code{credit_card_provider()}}{
credit card provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$credit_card_provider(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-credit_card_number"></a>}}
\if{latex}{\out{\hypertarget{method-credit_card_number}{}}}
\subsection{Method \code{credit_card_number()}}{
credit card number
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$credit_card_number(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-credit_card_security_code"></a>}}
\if{latex}{\out{\hypertarget{method-credit_card_security_code}{}}}
\subsection{Method \code{credit_card_security_code()}}{
credit card security code
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$credit_card_security_code(n = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of random things to generaate. an integer; default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FraudsterClient$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file-provider.R
\name{FileProvider}
\alias{FileProvider}
\title{FileProvider}
\description{
file methods
}
\examples{
(x <- FileProvider$new())
x$locale
x$mime_type()
x$file_extension()
x$file_name()
x$file_path()
x$file_path(depth = 2)
x$file_path(depth = 3)
x$file_path(depth = 6)
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{FileProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) the locale}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-allowed_locales}{\code{FileProvider$allowed_locales()}}
\item \href{#method-new}{\code{FileProvider$new()}}
\item \href{#method-mime_type}{\code{FileProvider$mime_type()}}
\item \href{#method-file_name}{\code{FileProvider$file_name()}}
\item \href{#method-file_extension}{\code{FileProvider$file_extension()}}
\item \href{#method-file_path}{\code{FileProvider$file_path()}}
\item \href{#method-clone}{\code{FileProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-allowed_locales"></a>}}
\if{latex}{\out{\hypertarget{method-allowed_locales}{}}}
\subsection{Method \code{allowed_locales()}}{
fetch the allowed locales for this provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileProvider$allowed_locales()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{FileProvider} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileProvider$new(locale = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. See
\verb{$allowed_locales()} for locales supported (default: en_US)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{FileProvider} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mime_type"></a>}}
\if{latex}{\out{\hypertarget{method-mime_type}{}}}
\subsection{Method \code{mime_type()}}{
a random mime type
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileProvider$mime_type(category = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{category}}{(character) a mime type category of mime types, one
of application, audio, image, message, model, multipart, text or
video. default: \code{NULL}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-file_name"></a>}}
\if{latex}{\out{\hypertarget{method-file_name}{}}}
\subsection{Method \code{file_name()}}{
a random file name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileProvider$file_name(category = NULL, extension = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{category}}{(character) a category of file extension type, one of
audio, image, office, text or video. default: \code{NULL}. If this is
given, \code{extension} is ignored}

\item{\code{extension}}{(character) a file extension. if this is given,
\code{category} is ignored.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-file_extension"></a>}}
\if{latex}{\out{\hypertarget{method-file_extension}{}}}
\subsection{Method \code{file_extension()}}{
a random file extension
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileProvider$file_extension(category = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{category}}{(character) a category of file extension type, one of
audio, image, office, text or video. default: \code{NULL}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-file_path"></a>}}
\if{latex}{\out{\hypertarget{method-file_path}{}}}
\subsection{Method \code{file_path()}}{
a random file path
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileProvider$file_path(depth = 1, category = NULL, extension = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{depth}}{(character) depth of the file (depth >= 0). default: 1}

\item{\code{category}}{(character) a category of file extension type, one of
audio, image, office, text or video. default: \code{NULL}. If this is
given, \code{extension} is ignored}

\item{\code{extension}}{(character) a file extension. if this is given,
\code{category} is ignored.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/numerics-provider.R
\name{NumericsProvider}
\alias{NumericsProvider}
\title{NumericsProvider}
\description{
numeric methods, generate numbers
}
\examples{
z <- NumericsProvider$new()

z$double()
z$double(10)
z$double(10, mean = 100)
z$double(10, mean = 100, sd = 17)

z$integer()
z$integer(10)
z$integer(10, 1, 20)
z$integer(10, 1, 10000000L)

z$unif()
z$norm()
z$lnorm(10)
z$beta(10, 1, 1)
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{NumericsProvider}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-double}{\code{NumericsProvider$double()}}
\item \href{#method-integer}{\code{NumericsProvider$integer()}}
\item \href{#method-unif}{\code{NumericsProvider$unif()}}
\item \href{#method-norm}{\code{NumericsProvider$norm()}}
\item \href{#method-lnorm}{\code{NumericsProvider$lnorm()}}
\item \href{#method-beta}{\code{NumericsProvider$beta()}}
\item \href{#method-clone}{\code{NumericsProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-double"></a>}}
\if{latex}{\out{\hypertarget{method-double}{}}}
\subsection{Method \code{double()}}{
get a double, pulls from normal distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NumericsProvider$double(n = 1, mean = 0, sd = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{(integer) number of values, default: 1}

\item{\code{mean}}{mean value, default: 0}

\item{\code{sd}}{standard deviation, default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-integer"></a>}}
\if{latex}{\out{\hypertarget{method-integer}{}}}
\subsection{Method \code{integer()}}{
get an integer, runs \code{\link[=sample]{sample()}} on range given
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NumericsProvider$integer(n = 1, min = 1, max = 1000)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{(integer) number of values, default: 1}

\item{\code{min}}{minimum value, default: 1}

\item{\code{max}}{maximum value, default: 1000}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-unif"></a>}}
\if{latex}{\out{\hypertarget{method-unif}{}}}
\subsection{Method \code{unif()}}{
get numbers from the uniform distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NumericsProvider$unif(n = 1, min = 0, max = 9999)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{(integer) number of values, default: 1}

\item{\code{min}}{minimum value, default: 1}

\item{\code{max}}{maximum value, default: 1000}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-norm"></a>}}
\if{latex}{\out{\hypertarget{method-norm}{}}}
\subsection{Method \code{norm()}}{
get numbers from the normal distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NumericsProvider$norm(n = 1, mean = 0, sd = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{(integer) number of values, default: 1}

\item{\code{mean}}{mean value, default: 0}

\item{\code{sd}}{standard deviation, default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-lnorm"></a>}}
\if{latex}{\out{\hypertarget{method-lnorm}{}}}
\subsection{Method \code{lnorm()}}{
get numbers from the lognormal distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NumericsProvider$lnorm(n = 1, mean = 0, sd = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{(integer) number of values, default: 1}

\item{\code{mean}}{mean value, default: 0}

\item{\code{sd}}{standard deviation, default: 1}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-beta"></a>}}
\if{latex}{\out{\hypertarget{method-beta}{}}}
\subsection{Method \code{beta()}}{
get numbers from the beta distribution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NumericsProvider$beta(n = 1, shape1, shape2, ncp = 0)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{(integer) number of values, default: 1}

\item{\code{shape1}}{non-negative parameters of the Beta distribution}

\item{\code{shape2}}{non-negative parameters of the Beta distribution}

\item{\code{ncp}}{non-centrality parameter, default: 0}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NumericsProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jobs-provider.R
\name{JobProvider}
\alias{JobProvider}
\title{JobProvider}
\description{
generate jobs
}
\examples{
z <- JobProvider$new()
z$render()

z <- JobProvider$new(locale = "fr_FR")
z$locale
z$render()

z <- JobProvider$new(locale = "hr_HR")
z$locale
z$render()

z <- JobProvider$new(locale = "fa_IR")
z$locale
z$render()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{JobProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) the locale}

\item{\code{formats}}{(character) vector of possible formats}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-allowed_locales}{\code{JobProvider$allowed_locales()}}
\item \href{#method-new}{\code{JobProvider$new()}}
\item \href{#method-render}{\code{JobProvider$render()}}
\item \href{#method-clone}{\code{JobProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-allowed_locales"></a>}}
\if{latex}{\out{\hypertarget{method-allowed_locales}{}}}
\subsection{Method \code{allowed_locales()}}{
fetch the allowed locales for this provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JobProvider$allowed_locales()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{JobProvider} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JobProvider$new(locale = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. See
\verb{$allowed_locales()} for locales supported (default: en_US)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{JobProvider} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-render"></a>}}
\if{latex}{\out{\hypertarget{method-render}{}}}
\subsection{Method \code{render()}}{
Make a job
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JobProvider$render()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JobProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/available_locales.R
\name{charlatan_locales}
\alias{charlatan_locales}
\title{Available locales}
\usage{
charlatan_locales()
}
\value{
a data.frame of the available locales in this package.
See \link{available_locales_df} for structure.

Not all functions support all locales. Check the docs for each one
to see what locales they support.

You can find out more about each locale by running your locale
though \code{stringi::stri_locale_info()}
}
\description{
Available locales
}
\examples{
charlatan_locales()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missing-data-provider.R
\name{MissingDataProvider}
\alias{MissingDataProvider}
\title{MissingDataProvider}
\description{
make missing data
}
\examples{
z <- MissingDataProvider$new()
z$make_missing(x = letters)
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{MissingDataProvider}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-make_missing}{\code{MissingDataProvider$make_missing()}}
\item \href{#method-clone}{\code{MissingDataProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-make_missing"></a>}}
\if{latex}{\out{\hypertarget{method-make_missing}{}}}
\subsection{Method \code{make_missing()}}{
make missing data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MissingDataProvider$make_missing(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{a vector of characters, numeric, integers, logicals, etc}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
This method picks a random number (\code{N}) of slots in
the input vector \code{x} (up to \code{length(x)}). Then picks \code{N} random
positions and replaces them with NA matching the input class.
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MissingDataProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/charlatan_settings.R
\name{charlatan_settings}
\alias{charlatan_settings}
\title{charlatan settings}
\usage{
charlatan_settings(messy = NULL)
}
\arguments{
\item{messy}{(logical) make some messy data. Default: \code{NULL}}
}
\description{
charlatan settings
}
\section{More deets}{

\itemize{
\item \code{messy} - When \code{FALSE}, nothing is different from normal.
When \code{TRUE}, we select incorrect/wrong values with probability X.
Messy mode is only available for \strong{en-US} for now, and only for
some data types. The default setting is \code{NULL}, meaning it is
ignored.
}
}

\examples{
charlatan_settings()
charlatan_settings(messy = TRUE)
charlatan_settings(messy = FALSE)

# with PersonProvider - overrides local messy param in all cases
x <- PersonProvider$new()
x$messy
charlatan_settings(messy = TRUE)
x <- PersonProvider$new()
x$messy
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/date_time.R
\name{date_time}
\alias{date_time}
\alias{ch_timezone}
\alias{ch_unix_time}
\alias{ch_date_time}
\title{Create dates and times}
\usage{
ch_timezone(n = 1)

ch_unix_time(n = 1)

ch_date_time(n = 1)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}
}
\description{
Create dates and times
}
\examples{
ch_timezone()
ch_timezone(10)

ch_unix_time()
ch_unix_time(20)

ch_date_time()
ch_date_time(20)
}
\seealso{
\link{DateTimeProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/charlatan-package.R
\docType{package}
\name{charlatan-package}
\alias{charlatan-package}
\alias{charlatan}
\title{charlatan}
\description{
Make fake data, supporting addresses, person names, dates,
times, colors, coordinates, currencies, digital object identifiers
(DOIs), jobs, phone numbers, DNA sequences, doubles and integers
from distributions and within a range.
}
\section{Package API}{

\itemize{
\item \code{\link[=ch_generate]{ch_generate()}}: generate a data.frame with fake data
\item \code{\link[=fraudster]{fraudster()}}: single interface to all fake data methods
\item High level interfaces: There are high level functions prefixed with
\code{ch_} that wrap low level interfaces, and are meant to be easier
to use and provide easy way to make many instances of a thing.
\item Low level interfaces: All of these are R6 objects that a user can
initialize and then call methods on the them.
}
}

\examples{
# generate individual types of data
ch_name()
ch_phone_number()
ch_job()

# generate a data.frame
ch_generate()

# one interface to all data types - generate the class first
#  reports the locale to be used, can change optionally
(x <- fraudster())
x$job()
x$name()
x$color_name()
x$hex_color()

# low level interfaces to "data providers"
# these are exported by hidden from package man page
# as most users will likely not interact with these
x <- ColorProvider$new()
x$color_name()
x$hex_color()
}
\author{
Scott Chamberlain \email{myrmecocystus+r@gmail.com}

Kyle Voytovich

Martin Pedersen
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/color.R
\name{ch_color}
\alias{ch_color}
\alias{ch_color_name}
\alias{ch_safe_color_name}
\alias{ch_hex_color}
\alias{ch_safe_hex_color}
\alias{ch_rgb_color}
\alias{ch_rgb_css_color}
\title{Create fake colors}
\usage{
ch_color_name(n = 1, locale = NULL)

ch_safe_color_name(n = 1, locale = NULL)

ch_hex_color(n = 1)

ch_safe_hex_color(n = 1)

ch_rgb_color(n = 1)

ch_rgb_css_color(n = 1)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}

\item{locale}{(character) the locale to use. See
\code{ColorProvider$new()$allowed_locales()} for locales supported.
Affects the \code{ch_color_name} and \code{ch_safe_color_name} functions}
}
\description{
Create fake colors
}
\examples{
ch_color_name()
ch_color_name(10)
ch_color_name(500)

ch_safe_color_name()
ch_safe_color_name(10)

ch_hex_color()
ch_hex_color(10)
ch_hex_color(1000)

ch_safe_hex_color()
ch_safe_hex_color(10)

ch_rgb_color()
ch_rgb_color(10)

ch_rgb_css_color()
ch_rgb_css_color(10)

ch_color_name(locale = "uk_UA")
ch_color_name(n = 10, locale = "uk_UA")

ch_safe_color_name(locale = "uk_UA")
ch_safe_color_name(n = 10, locale = "uk_UA")
}
\seealso{
\link{ColorProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lorem-provider.R
\name{LoremProvider}
\alias{LoremProvider}
\title{LoremProvider}
\description{
lorem ipsum methods
}
\examples{
(x <- LoremProvider$new())
x$locale
x$word()
x$words(3)
x$words(6)
x$sentence()
x$sentences(3)
x$sentences(6)
x$paragraph()
x$paragraphs(3)
x$paragraphs(6)
cat(x$paragraphs(6), sep = "\n")
x$text(6)
x$text(10)
x$text(19)
x$text(25)
x$text(50)
x$text(300)
x$text(2000)

# set a different sentence_punctuation or word_connector
(x <- LoremProvider$new(sentence_punctuation = ";"))
x$paragraph(4)
(x <- LoremProvider$new(word_connector = " --- "))
x$paragraph(4)

# different locales
LoremProvider$new(locale = "ar_AA")$word()
LoremProvider$new(locale = "el_GR")$word()
LoremProvider$new(locale = "he_IL")$word()
LoremProvider$new(locale = "ja_JP")$word()
LoremProvider$new(locale = "zh_TW")$word()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{LoremProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{locale}}{(character) the locale}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-allowed_locales}{\code{LoremProvider$allowed_locales()}}
\item \href{#method-new}{\code{LoremProvider$new()}}
\item \href{#method-word}{\code{LoremProvider$word()}}
\item \href{#method-words}{\code{LoremProvider$words()}}
\item \href{#method-sentence}{\code{LoremProvider$sentence()}}
\item \href{#method-sentences}{\code{LoremProvider$sentences()}}
\item \href{#method-paragraph}{\code{LoremProvider$paragraph()}}
\item \href{#method-paragraphs}{\code{LoremProvider$paragraphs()}}
\item \href{#method-text}{\code{LoremProvider$text()}}
\item \href{#method-clone}{\code{LoremProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-allowed_locales"></a>}}
\if{latex}{\out{\hypertarget{method-allowed_locales}{}}}
\subsection{Method \code{allowed_locales()}}{
fetch the allowed locales for this provider
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LoremProvider$allowed_locales()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{LoremProvider} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LoremProvider$new(
  locale = NULL,
  sentence_punctuation = ".",
  word_connector = " "
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locale}}{(character) the locale to use. See
\verb{$allowed_locales()} for locales supported (default: en_US)}

\item{\code{sentence_punctuation}}{(character) End of sentence punctuation}

\item{\code{word_connector}}{(character) Default connector between words}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{LoremProvider} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-word"></a>}}
\if{latex}{\out{\hypertarget{method-word}{}}}
\subsection{Method \code{word()}}{
Generate a random word
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LoremProvider$word(ext_words = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{ext_words}}{a character vector of words you would like to have
instead of 'Lorem ipsum'}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a single word
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-words"></a>}}
\if{latex}{\out{\hypertarget{method-words}{}}}
\subsection{Method \code{words()}}{
Generate a character vector of random words
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LoremProvider$words(nb = 3, ext_words = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nb}}{(integer) how many words to return}

\item{\code{ext_words}}{a character vector of words you would like to have
instead of 'Lorem ipsum'}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
many words
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-sentence"></a>}}
\if{latex}{\out{\hypertarget{method-sentence}{}}}
\subsection{Method \code{sentence()}}{
Generate a random sentence
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LoremProvider$sentence(
  nb_words = 6,
  variable_nb_words = TRUE,
  ext_words = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nb_words}}{(integer) around how many words the sentence should
contain}

\item{\code{variable_nb_words}}{set to \code{FALSE} if you want exactly \code{nb}
words returned, otherwise the result may include a number of words
of \code{nb} +/-40\% (with a minimum of 1)}

\item{\code{ext_words}}{a character vector of words you would like to have
instead of 'Lorem ipsum'}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a single sentence
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-sentences"></a>}}
\if{latex}{\out{\hypertarget{method-sentences}{}}}
\subsection{Method \code{sentences()}}{
Generate a character vector of random sentences
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LoremProvider$sentences(nb = 3, ext_words = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nb}}{(integer) how many sentences to return}

\item{\code{ext_words}}{a character vector of words you would like to have
instead of 'Lorem ipsum'}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
many sentences
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-paragraph"></a>}}
\if{latex}{\out{\hypertarget{method-paragraph}{}}}
\subsection{Method \code{paragraph()}}{
Generate a single paragraph
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LoremProvider$paragraph(
  nb_sentences = 3,
  variable_nb_sentences = TRUE,
  ext_words = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nb_sentences}}{(integer) around how many sentences the paragraph
should contain}

\item{\code{variable_nb_sentences}}{set to \code{FALSE} if you want exactly \code{nb}
sentences returned, otherwise the result may include a number of
sentences of \code{nb} +/-40\% (with a minimum of 1)}

\item{\code{ext_words}}{a character vector of words you would like to have
instead of 'Lorem ipsum'}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a single paragraph
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-paragraphs"></a>}}
\if{latex}{\out{\hypertarget{method-paragraphs}{}}}
\subsection{Method \code{paragraphs()}}{
Generate many paragraphs
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LoremProvider$paragraphs(nb = 3, ext_words = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nb}}{(integer) how many paragraphs to return}

\item{\code{ext_words}}{a character vector of words you would like to have
instead of 'Lorem ipsum'}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
many paragraphs
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-text"></a>}}
\if{latex}{\out{\hypertarget{method-text}{}}}
\subsection{Method \code{text()}}{
Generate a random text string. Depending on the
\code{max_nb_chars}, returns a string made of words, sentences, or
paragraphs.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LoremProvider$text(max_nb_chars = 200, ext_words = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{max_nb_chars}}{Maximum number of characters the text should
contain (minimum 5)}

\item{\code{ext_words}}{a character vector of words you would like to have
instead of 'Lorem ipsum'}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
character string of words
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LoremProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base-provider.R
\name{BaseProvider}
\alias{BaseProvider}
\title{BaseProvider}
\description{
BaseProvider

BaseProvider
}
\examples{
(x <- BaseProvider$new())

x$numerify("#\%\%asdf221?")
x$lexify("#\%\%asdf221?")
x$bothify("#\%\%asdf221?")

z <- PhoneNumberProvider$new()
x$numerify(z$render())

x$random_element(letters)
x$random_int()
x$random_digit()
x$random_digit_not_zero()
x$random_digit_or_empty()
x$random_digit_not_zero_or_empty()
x$random_letter()
x$check_locale("es_ES")
## fails
# x$check_locale("es_EQ")

x$randomize_nb_elements()
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-random_element}{\code{BaseProvider$random_element()}}
\item \href{#method-random_element_prob}{\code{BaseProvider$random_element_prob()}}
\item \href{#method-random_int}{\code{BaseProvider$random_int()}}
\item \href{#method-random_digit}{\code{BaseProvider$random_digit()}}
\item \href{#method-random_digit_not_zero}{\code{BaseProvider$random_digit_not_zero()}}
\item \href{#method-random_digit_or_empty}{\code{BaseProvider$random_digit_or_empty()}}
\item \href{#method-random_digit_not_zero_or_empty}{\code{BaseProvider$random_digit_not_zero_or_empty()}}
\item \href{#method-random_letter}{\code{BaseProvider$random_letter()}}
\item \href{#method-numerify}{\code{BaseProvider$numerify()}}
\item \href{#method-lexify}{\code{BaseProvider$lexify()}}
\item \href{#method-bothify}{\code{BaseProvider$bothify()}}
\item \href{#method-check_locale}{\code{BaseProvider$check_locale()}}
\item \href{#method-randomize_nb_elements}{\code{BaseProvider$randomize_nb_elements()}}
\item \href{#method-clone}{\code{BaseProvider$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-random_element"></a>}}
\if{latex}{\out{\hypertarget{method-random_element}{}}}
\subsection{Method \code{random_element()}}{
pick a random element from vector/list
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$random_element(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{vector or list}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a single element from x
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-random_element_prob"></a>}}
\if{latex}{\out{\hypertarget{method-random_element_prob}{}}}
\subsection{Method \code{random_element_prob()}}{
pick a random element with probability from vector/list
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$random_element_prob(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{vector or list}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-random_int"></a>}}
\if{latex}{\out{\hypertarget{method-random_int}{}}}
\subsection{Method \code{random_int()}}{
any number of random integers from a min, max
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$random_int(min = 0, max = 9999, size = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{min}}{the minimum value. default: 0}

\item{\code{max}}{the maximum value. default: 9999}

\item{\code{size}}{number of values to return. default: 1}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
random integer
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-random_digit"></a>}}
\if{latex}{\out{\hypertarget{method-random_digit}{}}}
\subsection{Method \code{random_digit()}}{
random integer between 0 and 9
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$random_digit()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-random_digit_not_zero"></a>}}
\if{latex}{\out{\hypertarget{method-random_digit_not_zero}{}}}
\subsection{Method \code{random_digit_not_zero()}}{
random integer between 1 and 9
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$random_digit_not_zero()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-random_digit_or_empty"></a>}}
\if{latex}{\out{\hypertarget{method-random_digit_or_empty}{}}}
\subsection{Method \code{random_digit_or_empty()}}{
random integer between 0 and 9 or empty character string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$random_digit_or_empty()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-random_digit_not_zero_or_empty"></a>}}
\if{latex}{\out{\hypertarget{method-random_digit_not_zero_or_empty}{}}}
\subsection{Method \code{random_digit_not_zero_or_empty()}}{
random integer between 1 and 9 or empty character string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$random_digit_not_zero_or_empty()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-random_letter"></a>}}
\if{latex}{\out{\hypertarget{method-random_letter}{}}}
\subsection{Method \code{random_letter()}}{
random letter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$random_letter()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-numerify"></a>}}
\if{latex}{\out{\hypertarget{method-numerify}{}}}
\subsection{Method \code{numerify()}}{
replace a template with numbers
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$numerify(text = "###")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{text}}{(character) a string}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-lexify"></a>}}
\if{latex}{\out{\hypertarget{method-lexify}{}}}
\subsection{Method \code{lexify()}}{
replace a template with letters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$lexify(text = "????")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{text}}{(character) a string}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-bothify"></a>}}
\if{latex}{\out{\hypertarget{method-bothify}{}}}
\subsection{Method \code{bothify()}}{
both numerify and lexify together
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$bothify(text = "## ??")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{text}}{(character) a string}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-check_locale"></a>}}
\if{latex}{\out{\hypertarget{method-check_locale}{}}}
\subsection{Method \code{check_locale()}}{
check a locale to see if it exists, if not, stop with
error message
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$check_locale(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{a locale name, e.g, 'bg_BG'}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
returns nothing if locale is supported; stops w/ message if not
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-randomize_nb_elements"></a>}}
\if{latex}{\out{\hypertarget{method-randomize_nb_elements}{}}}
\subsection{Method \code{randomize_nb_elements()}}{
Returns a random value near number
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$randomize_nb_elements(
  number = 10,
  le = FALSE,
  ge = FALSE,
  min = NULL,
  max = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{number}}{value to which the result must be near}

\item{\code{le}}{result must be lower or equal to number}

\item{\code{ge}}{result must be greater or equal to number}

\item{\code{min}}{the minimum value. default: \code{NULL}}

\item{\code{max}}{the maximum value. default: \code{NULL}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a random int near number
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BaseProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/credit_card.R
\name{ch_credit}
\alias{ch_credit}
\alias{ch_credit_card_provider}
\alias{ch_credit_card_number}
\alias{ch_credit_card_security_code}
\title{Create fake credit card data}
\usage{
ch_credit_card_provider(n = 1)

ch_credit_card_number(n = 1)

ch_credit_card_security_code(n = 1)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}
}
\description{
Create fake credit card data
}
\examples{
ch_credit_card_provider()
ch_credit_card_provider(n = 4)

ch_credit_card_number()
ch_credit_card_number(n = 10)
ch_credit_card_number(n = 500)

ch_credit_card_security_code()
ch_credit_card_security_code(n = 10)
ch_credit_card_security_code(n = 500)
}
\seealso{
\link{CreditCardProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datetime-provider.R
\name{DateTimeProvider}
\alias{DateTimeProvider}
\title{DateTimeProvider}
\description{
date and time methods
}
\examples{
z <- DateTimeProvider$new()
z$countries
z$centuries
z$century()
z$timezone()
z$unix_time()
z$date("\%Y-\%M-\%d")
z$date_time()
z$year()
z$iso8601("1932-02-12 05:32:12")
# z$iso8601("January 4, 1981")

# date time between a range of dates
(start_date <- Sys.time() - 604800L)
z$date_time_between(start_date = start_date)
# in the year 1900
z$date_time_between("1900-01-01 00:00:00 PST", "1900-12-31 00:00:00 PST")
z$date_time_between("1900-01-01", "1900-12-31")
}
\references{
https://en.wikipedia.org/wiki/Unix_time

https://en.wikipedia.org/wiki/Unix_time
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{DateTimeProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{centuries}}{(character) centuries in roman numerals}

\item{\code{countries}}{(list) countries list}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-unix_time}{\code{DateTimeProvider$unix_time()}}
\item \href{#method-date}{\code{DateTimeProvider$date()}}
\item \href{#method-date_time}{\code{DateTimeProvider$date_time()}}
\item \href{#method-date_time_fromtimestamp}{\code{DateTimeProvider$date_time_fromtimestamp()}}
\item \href{#method-iso8601}{\code{DateTimeProvider$iso8601()}}
\item \href{#method-year}{\code{DateTimeProvider$year()}}
\item \href{#method-century}{\code{DateTimeProvider$century()}}
\item \href{#method-timezone}{\code{DateTimeProvider$timezone()}}
\item \href{#method-date_time_between}{\code{DateTimeProvider$date_time_between()}}
\item \href{#method-clone}{\code{DateTimeProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-unix_time"></a>}}
\if{latex}{\out{\hypertarget{method-unix_time}{}}}
\subsection{Method \code{unix_time()}}{
Get a timestamp between January 1, 1970 and now, unless passed
explicit \code{start_date} or \code{end_date} values
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DateTimeProvider$unix_time(start_date = NULL, end_date = "now")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{start_date}}{start date, a valid date format}

\item{\code{end_date}}{start date, a valid date format, default: "now"}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-date"></a>}}
\if{latex}{\out{\hypertarget{method-date}{}}}
\subsection{Method \code{date()}}{
Generate a date between January 1, 1970 and now,
with given pattern
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DateTimeProvider$date(pattern = "\%Y-\%m-\%d")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{pattern}}{date pattern, default: \verb{\%Y-\%m-\%d}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-date_time"></a>}}
\if{latex}{\out{\hypertarget{method-date_time}{}}}
\subsection{Method \code{date_time()}}{
Generate a date time between January 1, 1970 and now
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DateTimeProvider$date_time(tzinfo = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{tzinfo}}{timezone, see \link{timezone}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-date_time_fromtimestamp"></a>}}
\if{latex}{\out{\hypertarget{method-date_time_fromtimestamp}{}}}
\subsection{Method \code{date_time_fromtimestamp()}}{
Generate a iso8601 format date
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DateTimeProvider$date_time_fromtimestamp(timestamp, tzinfo = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{timestamp}}{a timestamp}

\item{\code{tzinfo}}{timezone, see \link{timezone}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-iso8601"></a>}}
\if{latex}{\out{\hypertarget{method-iso8601}{}}}
\subsection{Method \code{iso8601()}}{
Generate a iso8601 format date
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DateTimeProvider$iso8601(date, tzinfo = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{date}}{a date, in a valid date format}

\item{\code{tzinfo}}{timezone, see \link{timezone}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-year"></a>}}
\if{latex}{\out{\hypertarget{method-year}{}}}
\subsection{Method \code{year()}}{
generate a year
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DateTimeProvider$year()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-century"></a>}}
\if{latex}{\out{\hypertarget{method-century}{}}}
\subsection{Method \code{century()}}{
generate a century
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DateTimeProvider$century()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-timezone"></a>}}
\if{latex}{\out{\hypertarget{method-timezone}{}}}
\subsection{Method \code{timezone()}}{
generate a timezone
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DateTimeProvider$timezone()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-date_time_between"></a>}}
\if{latex}{\out{\hypertarget{method-date_time_between}{}}}
\subsection{Method \code{date_time_between()}}{
Generate a date time based on a random date between
two given dates
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DateTimeProvider$date_time_between(start_date, end_date = "now", tzinfo = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{start_date}}{start date, a valid date format}

\item{\code{end_date}}{start date, a valid date format, default: "now"}

\item{\code{tzinfo}}{timezone, see \link{timezone}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DateTimeProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/job.R
\name{ch_job}
\alias{ch_job}
\title{Create fake jobs}
\usage{
ch_job(n = 1, locale = NULL)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}

\item{locale}{(character) the locale to use. Run
\code{JobProvider$new()$allowed_locales()} for locales supported
(default: en_US)}
}
\description{
Create fake jobs
}
\examples{
ch_job()
ch_job(10)
ch_job(500)

ch_job(locale = "da_DK", n = 10)
ch_job(locale = "fi_FI", n = 10)
ch_job(locale = "fr_FR", n = 10)
ch_job(locale = "fr_CH", n = 10)
ch_job(locale = "hr_HR", n = 10)
ch_job(locale = "fa_IR", n = 10)
ch_job(locale = "pl_PL", n = 10)
ch_job(locale = "ru_RU", n = 10)
ch_job(locale = "uk_UA", n = 10)
ch_job(locale = "zh_TW", n = 10)
}
\seealso{
\link{JobProvider}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doi-provider.R
\name{DOIProvider}
\alias{DOIProvider}
\title{DOIProvider}
\description{
DOIProvider

DOIProvider
}
\examples{
(z <- DOIProvider$new())
z$render()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{DOIProvider}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{funs}}{(list) list of functions to use to apply to DOI creation}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-render}{\code{DOIProvider$render()}}
\item \href{#method-clone}{\code{DOIProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-render"></a>}}
\if{latex}{\out{\hypertarget{method-render}{}}}
\subsection{Method \code{render()}}{
Make a random DOI
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DOIProvider$render()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DOIProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/element-provider.R
\name{ElementProvider}
\alias{ElementProvider}
\title{ElementProvider}
\description{
chemical elements methods
}
\details{
Data from Wikipedia at
\url{https://en.wikipedia.org/wiki/Chemical_element}
}
\examples{
z <- ElementProvider$new()
z$symbol()
z$element()
}
\keyword{internal}
\section{Super class}{
\code{\link[charlatan:BaseProvider]{charlatan::BaseProvider}} -> \code{ElementProvider}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-symbol}{\code{ElementProvider$symbol()}}
\item \href{#method-element}{\code{ElementProvider$element()}}
\item \href{#method-clone}{\code{ElementProvider$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="bothify">}\href{../../charlatan/html/BaseProvider.html#method-bothify}{\code{charlatan::BaseProvider$bothify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="check_locale">}\href{../../charlatan/html/BaseProvider.html#method-check_locale}{\code{charlatan::BaseProvider$check_locale()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="lexify">}\href{../../charlatan/html/BaseProvider.html#method-lexify}{\code{charlatan::BaseProvider$lexify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="numerify">}\href{../../charlatan/html/BaseProvider.html#method-numerify}{\code{charlatan::BaseProvider$numerify()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit">}\href{../../charlatan/html/BaseProvider.html#method-random_digit}{\code{charlatan::BaseProvider$random_digit()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero}{\code{charlatan::BaseProvider$random_digit_not_zero()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_not_zero_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_not_zero_or_empty}{\code{charlatan::BaseProvider$random_digit_not_zero_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_digit_or_empty">}\href{../../charlatan/html/BaseProvider.html#method-random_digit_or_empty}{\code{charlatan::BaseProvider$random_digit_or_empty()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element">}\href{../../charlatan/html/BaseProvider.html#method-random_element}{\code{charlatan::BaseProvider$random_element()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_element_prob">}\href{../../charlatan/html/BaseProvider.html#method-random_element_prob}{\code{charlatan::BaseProvider$random_element_prob()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_int">}\href{../../charlatan/html/BaseProvider.html#method-random_int}{\code{charlatan::BaseProvider$random_int()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="random_letter">}\href{../../charlatan/html/BaseProvider.html#method-random_letter}{\code{charlatan::BaseProvider$random_letter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="charlatan" data-topic="BaseProvider" data-id="randomize_nb_elements">}\href{../../charlatan/html/BaseProvider.html#method-randomize_nb_elements}{\code{charlatan::BaseProvider$randomize_nb_elements()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-symbol"></a>}}
\if{latex}{\out{\hypertarget{method-symbol}{}}}
\subsection{Method \code{symbol()}}{
Get a symbol
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ElementProvider$symbol()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-element"></a>}}
\if{latex}{\out{\hypertarget{method-element}{}}}
\subsection{Method \code{element()}}{
Get an element
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ElementProvider$element()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ElementProvider$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/available_locales.R
\docType{data}
\name{available_locales_df}
\alias{available_locales_df}
\title{Available locales}
\format{
A data frame with 45 rows and 4 variables:
\describe{
\item{Language}{language two letter code}
\item{Country}{country two letter code}
\item{Variant}{a variant code, if applicable}
\item{Name}{official locale two letter code}
}
}
\description{
A data.frame of locales available in \pkg{charlatan}
}
\seealso{
data.frame used in \code{\link[=charlatan_locales]{charlatan_locales()}}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinates.R
\name{coordinates}
\alias{coordinates}
\alias{ch_lon}
\alias{ch_lat}
\alias{ch_position}
\title{Create fake coordinates}
\usage{
ch_lon(n = 1)

ch_lat(n = 1)

ch_position(n = 1, bbox = NULL)
}
\arguments{
\item{n}{(integer) number of things to get, any non-negative integer}

\item{bbox}{a bounding box of the form \verb{[w,s,e,n]}}
}
\description{
Create fake coordinates
}
\examples{
ch_lon()
ch_lon(10)

ch_lat()
ch_lat(10)

ch_position()
ch_position(10)
ch_position(bbox = c(-120, 30, -110, 60))
}
\seealso{
\link{CoordinateProvider}
}
