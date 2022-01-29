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

*Wow, no problems at all. :)**Wow, no problems at all. :)*