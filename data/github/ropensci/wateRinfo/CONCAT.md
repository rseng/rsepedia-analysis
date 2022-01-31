<!-- README.md is generated from README.Rmd. Please edit that file and knit -->



# wateRinfo <img src="man/figures/logo.png" align="right" alt="" width="120">

<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/wateRinfo/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/wateRinfo/actions) [![codecov](https://codecov.io/gh/ropensci/wateRinfo/branch/master/graph/badge.svg?token=4MT1hZyw2l)](https://codecov.io/gh/ropensci/wateRinfo)
<!-- badges: end -->

wateRinfo facilitates access to [waterinfo.be](https://www.waterinfo.be/), a website managed by the [Flanders Environment Agency (VMM)](https://en.vmm.be/) and [Flanders Hydraulics Research](https://www.waterbouwkundiglaboratorium.be/). The website provides access to real-time water and weather related environmental variables for Flanders (Belgium), such as rainfall, air pressure, discharge, and water level. The package provides functions to search for stations and variables, and download time series.

To get started, see:

* [Function reference](https://docs.ropensci.org/wateRinfo/reference/index.html): an overview of all wateRinfo functions.
* [Articles](https://docs.ropensci.org/wateRinfo/articles/): tutorials on how to use the package.

## Installation

You can install wateRinfo from [GitHub](https://github.com/ropensci/wateRinfo) with:


```r
# install.packages("devtools")
devtools::install_github("ropensci/wateRinfo")
```

When successful, load it as usual:


```r
library(wateRinfo)
```

## Example

For a number of supported variables ([documented](https://www.waterinfo.be/download/9f5ee0c9-dafa-46de-958b-7cac46eb8c23?dl=0) by VMM), the stations providing time series data for a given variable can be listed with the command `get_stations()`.

If you want to know the supported variables, ask for the supported variables:


```r
supported_variables("en")
#>              variable_en
#> 1              discharge
#> 6        soil_saturation
#> 7          soil_moisture
#> 8  dew_point_temperature
#> 9     ground_temperature
#> 10           ground_heat
#> 11            irradiance
#> 12          air_pressure
#> 13 air_temperature_175cm
#> 14              rainfall
#> 20     relative_humidity
#> 21  evaporation_monteith
#> 25    evaporation_penman
#> 29        water_velocity
#> 34           water_level
#> 39     water_temperature
#> 40        wind_direction
#> 41            wind_speed
```

Listing the available air pressure stations:


```r
get_stations("air_pressure")
#>      ts_id station_latitude station_longitude station_id station_no            station_name
#> 1 78124042         51.20300          5.439589      12213   ME11_002             Overpelt_ME
#> 2 78005042         51.02263          2.970584      12206   ME01_003               Zarren_ME
#> 3 78039042         51.24379          4.266912      12208   ME04_001              Melsele_ME
#> 4 78073042         50.88663          4.094898      12210   ME07_006           Liedekerke_ME
#> 5 78107042         51.16224          4.845708      12212   ME10_011            Herentals_ME
#> 6 78022042         51.27226          3.728299      12207   ME03_017            Boekhoute_ME
#> 7 78056042         50.86149          3.411318      12209   ME05_019              Waregem_ME
#> 8 78090042         50.73795          5.141976      12211   ME09_012 Niel-bij-St.-Truiden_ME
#>   stationparameter_name parametertype_name ts_unitsymbol dataprovider
#> 1                    Pa                 Pa           hPa          VMM
#> 2                    Pa                 Pa           hPa          VMM
#> 3                    Pa                 Pa           hPa          VMM
#> 4                    Pa                 Pa           hPa          VMM
#> 5                    Pa                 Pa           hPa          VMM
#> 6                    Pa                 Pa           hPa          VMM
#> 7                    Pa                 Pa           hPa          VMM
#> 8                    Pa                 Pa           hPa          VMM
```

Each of the stations in the list for a given variable, are represented by a `ts_id`. These can be used to download the data of a given period with the command `get_timeseries_tsid()`, for example Overpelt (`ts_id = 78124042`):


```r
overpelt_pressure <- get_timeseries_tsid("78124042", 
                                         from = "2017-04-01", 
                                         to = "2017-04-02")
head(overpelt_pressure)
#>             Timestamp  Value Quality Code
#> 1 2017-04-01 00:00:00 1008.8          130
#> 2 2017-04-01 00:15:00 1008.7          130
#> 3 2017-04-01 00:30:00 1008.7          130
#> 4 2017-04-01 00:45:00 1008.6          130
#> 5 2017-04-01 01:00:00 1008.5          130
#> 6 2017-04-01 01:15:00 1008.4          130
```

Making a plot of the data with [`ggplot2`](https://ggplot2.tidyverse.org/):


```r
library(ggplot2)
ggplot(overpelt_pressure, aes(x = Timestamp, y = Value)) + 
    geom_line() + 
    xlab("") + ylab("hPa") + 
    scale_x_datetime(date_labels = "%H:%M\n%Y-%m-%d", date_breaks = "6 hours")
```

<img src="man/figures/README-plot_pressure-1.png" title="plot of chunk showplot1" alt="plot of chunk showplot1" width="80%" />

Another option is to check the available variables for a given station, with the function `get_variables()`. Let's consider again Overpelt (`ME11_002`) and check the first ten available variables at the Overpelt measurement station:


```r
vars_overpelt <- get_variables("ME11_002")
head(vars_overpelt, 10)
#>    station_name station_no    ts_id    ts_name parametertype_name stationparameter_name
#> 1   Overpelt_ME   ME11_002 78522042 HydJaarMax                 Ts                 SoilT
#> 2   Overpelt_ME   ME11_002 78523042 HydJaarMin                 Ts                 SoilT
#> 3   Overpelt_ME   ME11_002 78693042       P.15                 Ud                  WDir
#> 4   Overpelt_ME   ME11_002 94682042   MaandMin                 Ta                    Ta
#> 5   Overpelt_ME   ME11_002 78531042       P.10                 Ts                 SoilT
#> 6   Overpelt_ME   ME11_002 78518042     DagGem                 Ts                 SoilT
#> 7   Overpelt_ME   ME11_002 78521042 HydJaarGem                 Ts                 SoilT
#> 8   Overpelt_ME   ME11_002 78524042 KalJaarGem                 Ts                 SoilT
#> 9   Overpelt_ME   ME11_002 78533042       P.60                 Ts                 SoilT
#> 10  Overpelt_ME   ME11_002 78694042      Pv.15                 Ud                  WDir
```

Different pre-calculated variables are already available and a `ts_id` value is available for each of them to download the corresponding data. For example, `DagGem` (= daily mean values) of `RH` (= relative humidity), i.e. `ts_id = 78382042`:


```r
overpelt_rh_daily <- get_timeseries_tsid("78382042", 
                                         from = "2017-04-01", 
                                         to = "2017-04-30")
head(overpelt_rh_daily)
#>             Timestamp Value Quality Code
#> 1 2017-04-01 23:00:00 80.19          130
#> 2 2017-04-02 23:00:00 89.58          130
#> 3 2017-04-03 23:00:00 79.56          130
#> 4 2017-04-04 23:00:00 84.13          130
#> 5 2017-04-05 23:00:00 84.19          130
#> 6 2017-04-06 23:00:00 82.71          130
```


```r
ggplot(overpelt_rh_daily, aes(x = Timestamp, y = Value)) + 
    geom_line() + 
    xlab("") + ylab(" RH (%)") + 
    scale_x_datetime(date_labels = "%b-%d\n%Y", date_breaks = "5 days")
```

<img src="man/figures/README-plot_rh-1.png" title="plot of chunk showplot2" alt="plot of chunk showplot2" width="80%" />

Unfortunately, not all variables are documented, for which the check for the appropriate variable is not (yet) fully supported by the package.

More detailed tutorials are available in the package vignettes!

## Note on restrictions of the downloads

The amount of data downloaded from waterinfo.be is limited via a credit system. You do not need to get a token right away to download data. For limited and irregular downloads, a token will not be required.

When you require more extended data requests, please request a download token from the waterinfo.be site administrators via the e-mail address <hydrometrie@waterinfo.be> with a statement of which data and how frequently you would like to download data. You will then receive a client-credit code that can be used to obtain a token that is valid for 24 hours, after which the token can be refreshed with the same client-credit code.

Get token with client-credit code: (limited client-credit code for testing purposes)


```r
client <- paste0("MzJkY2VlY2UtODI2Yy00Yjk4LTljMmQtYjE2OTc4ZjBjYTZhOjRhZGE4",
                 "NzFhLTk1MjgtNGI0ZC1iZmQ1LWI1NzBjZThmNGQyZA==")
my_token <- get_token(client = client)
print(my_token)
#> Token:
#> eyJhbGciOiJIUzI1NiJ9.eyJqdGkiOiI5ODI5NDIzYS01NTNjLTQ3YTUtODUzNS1hZTBhN2FmMTFhN2MiLCJpYXQiOjE2MTk3Nzg0NDYsImlzcyI6Imh0dHA6Ly9sb2NhbGhvc3Q6ODA4MC9LaVdlYlBvcnRhbC9hdXRoIiwiYXVkIjoiMzJkY2VlY2UtODI2Yy00Yjk4LTljMmQtYjE2OTc4ZjBjYTZhIiwiZXhwIjoxNjE5ODY0ODQ2fQ.7pUqf8x0OxE-sA0PJUcKYGysl-DI5-KiodJ1ahfaMCA
#> 
#> Attributes:
#>  url: http://download.waterinfo.be/kiwis-auth/token
#>  type: Bearer
#>  expires: 2021-05-01 12:27:26 CEST
```

Receive information on the validity of the token:


```r
is.expired(my_token)
#> [1] FALSE
```

Check when the token expires:


```r
expires.in(my_token)
#> Time difference of 24 hours
```

Use token when retrieving data:


```r
get_stations(variable_name = "verdamping_monteith", token = my_token)
#>      ts_id station_latitude station_longitude station_id station_no            station_name
#> 1 94310042         51.02263          2.970584      12206   ME01_003               Zarren_ME
#> 2 94530042         51.16224          4.845708      12212   ME10_011            Herentals_ME
#> 3 94544042         51.20300          5.439589      12213   ME11_002             Overpelt_ME
#> 4 94516042         50.73795          5.141976      12211   ME09_012 Niel-bij-St.-Truiden_ME
#> 5 94488042         50.86149          3.411318      12209   ME05_019              Waregem_ME
#> 6 94502042         50.88663          4.094898      12210   ME07_006           Liedekerke_ME
#> 7 94474042         51.24379          4.266912      12208   ME04_001              Melsele_ME
#> 8 94460042         51.27226          3.728299      12207   ME03_017            Boekhoute_ME
#>   stationparameter_name parametertype_name ts_unitsymbol dataprovider
#> 1                   pET                PET            mm          VMM
#> 2                   pET                PET            mm          VMM
#> 3                   pET                PET            mm          VMM
#> 4                   pET                PET            mm          VMM
#> 5                   pET                PET            mm          VMM
#> 6                   pET                PET            mm          VMM
#> 7                   pET                PET            mm          VMM
#> 8                   pET                PET            mm          VMM
```

## Other clients

Besides this wateRinfo R client to gather data from [waterinfo.be](https://www.waterinfo.be/), there is also a Python client available. The [pywaterinfo](https://fluves.github.io/pywaterinfo/) package contains similar functionalities.

The [Flanders Hydraulics Research center](https://www.waterbouwkundiglaboratorium.be/en/) also distributes clients for R, Python and Matlab upon request to download the data they share on [waterinfo.be](https://www.waterinfo.be/). For more information, contact them directly via [hic@vlaanderen.be](mailto:hic@vlaanderen.be).

## Acknowledgements

This package is just a small wrapper around waterinfo.be to facilitate researchers and other stakeholders in downloading the data from [waterinfo.be](http://www.waterinfo.be). The availability of this data is made possible by *de Vlaamse Milieumaatschappij, Waterbouwkundig Laboratorium, Maritieme Dienstverlening & Kust, Waterwegen en Zeekanaal NV en De Scheepvaart NV*.

## Meta

* We welcome [contributions](.github/CONTRIBUTING.md) including bug reports.
* License: MIT
* Get citation information for `wateRinfo` in R doing `citation("wateRinfo")`.
* Please note that this project is released with a [Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
wateRinfo 0.3.0.9047 (2018-12-13)
=============================

### NEW FEATURES

* WateRinfo moved to ropensci, thanks to @ldecicco-USGS for the review and @karthik for editing the submission

### MINOR IMPROVEMENTS

* Change the URLs to `https` version, which are now supported by waterinfo.be
* Make sure data.frame content are characters and not factors
* Return supported frequencies as vector instead of single ling character

### DOCUMENTATION FIXES

* Update to usage of pkgdown 1.3.0
* Moved link references from INBO tp ropensci
* Add ropensci footer
* Mark Institute of Nature and Forest Research (INBO) as copyright holder
* Updated logo base color

Notice, we will add an in-development fourth component to the release tag with the development version. For more information, see [package guide](http://r-pkgs.had.co.nz/description.html#version).

wateRinfo 0.2.2 (2018-11-07)
=============================

### CRUCIAL QUICK FIX

* When no `datasource` is added to the waterinfo.be query, the server gives a html/text response containing a tomcat server error messag. This error is now captured properly. As the user normally only uses the wrapped API call functions, there should no difference on the user level.

### DOCUMENTATION FIXES

* Bring documentation up to date with the new `datasource` handling

wateRinfo 0.2.1 (2018-10-17)
=============================

The `datasource` of the non-VMM stations changed recently (our tests on CI broke). Moreover, the 
`ts_identifiers` for these stations changed as well. Both are adaptations on the data source (waterinfo.be) level. 

### CRUCIAL QUICK FIX

* Change the datasource to 4 for non-VMM stations
* Adapt the examples and vignettes to existing `ts_identifiers`

wateRinfo 0.2.0 (2018-10-01)
=============================

### NEW FEATURES

* Better support for changed station identifiers of the waterinfo.be dataproviders 

### MINOR IMPROVEMENTS

* Change print statements for message statements
* Add missing test to check token object type (remark: current missing tests are hard to introduce without having control on the original API behaviour)
* Test files have consistent names
* Update code style using styler

### DOCUMENTATION FIXES

* Describe the table output format explicitly in roxygen docs and add the quality codes of the waterinfo.be dataproviders (#27)
* Update contributing guidelines, thanks to @peterdesmet who did a great effort on providing welcoming, clear written and engaging contributing guidelines
* Provide code examples for all public functions
* Provide top-level documentation handle `?wateRinfo` (#23)
* Add keywords internal of internal functions not part of the documentation index (#24)


wateRinfo 0.1.2 (2018-05-03)
============================

* Added a `NEWS.md` file to track changes to the package (#25).
* Provide support for parsible non-existing dates (e.g. '2018-04-31')
* Add contributing guidelines


wateRinfo 0.1.1 (2017-11-29)
============================

* Add support for token handling as distributed by waterinfo.be



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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
# Contributing to wateRinfo

First of all, thanks for considering contributing to wateRinfo! 👍 It's people like you that make it rewarding for us - the project maintainers - to work on wateRinfo. 😊

wateRinfo is an open source project, maintained by people who care. We are not directly funded to do so.

[repo]: https://github.com/ropensci/wateRinfo
[issues]: https://github.com/ropensci/wateRinfo/issues
[new_issue]: https://github.com/ropensci/wateRinfo/issues/new
[website]: https://ropensci.github.io/wateRinfo
[citation]: https://ropensci.github.io/wateRinfo/authors.html
[email]: mailto:stijnvanhoey@gmail.com

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

### Share the love ❤️

Think wateRinfo is useful? Let others discover it, by telling them in person, via Twitter or a blog post.

Using wateRinfo for a paper you are writing? Consider [citing it][citation].

### Ask a question ⁉️

Using wateRinfo and got stuck? Browse the [documentation][website] to see if you can find a solution. Still stuck? Post your question as an [issue on GitHub][new_issue]. While we cannot offer user support, we'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [email][email].

### Propose an idea 💡

Have an idea for a new wateRinfo feature? Take a look at the [documentation][website] and [issue list][issues] to see if it isn't included or suggested yet. If not, suggest your idea as an [issue on GitHub][new_issue]. While we can't promise to implement your idea, it helps to:

* Explain in detail how it would work.
* Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Report a bug 🐛

Using wateRinfo and discovered a bug? That's annoying! Don't let others have the same experience and report it as an [issue on GitHub][new_issue] so we can fix it. A good bug report makes it easier for us to do so, so please include:

* Your operating system name and version (e.g. Mac OS 10.13.6).
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

### Improve the documentation 📖

Noticed a typo on the website? Think a function could use a better example? Good documentation makes all the difference, so your help to improve it is very welcome!

#### The website

[This website][website] is generated with [`pkgdown`](http://pkgdown.r-lib.org/). That means we don't have to write any html: content is pulled together from documentation in the code, vignettes, [Markdown](https://guides.github.com/features/mastering-markdown/) files, the package `DESCRIPTION` and `_pkgdown.yml` settings. If you know your way around `pkgdown`, you can [propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to improve documentation. If not, [report an issue][new_issue] and we can point you in the right direction.

#### Function documentation

Functions are described as comments near their code and translated to documentation using [`roxygen2`](https://klutometis.github.io/roxygen/). If you want to improve a function description:

1. Go to `R/` directory in the [code repository][repo].
2. Look for the file with the name of the function.
3. [Propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to update the function documentation in the roxygen comments (starting with `#'`).

### Contribute code 📝

Care to fix bugs or implement new functionality for wateRinfo? Awesome! 👏 Have a look at the [issue list][issues] and leave a comment on the things you want to work on. See also the development guidelines below.

## Development guidelines

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [this repo][repo] and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/).
2. If you have forked and cloned the project before and it has been a while since you worked on it, [pull changes from the original repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/) to your clone by using `git pull upstream master`.
3. Open the RStudio project file (`.Rproj`).
5. Make your changes:
    * Write your code.
    * Test your code (bonus points for adding unit tests).
    * Document your code (see function documentation above).
    * Do an `R CMD check` using `devtools::check()` and aim for 0 errors and warnings.
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request).
<!-- README.md is generated from README.Rmd. Please edit that file and knit -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-"
)
```

# wateRinfo <img src="man/figures/logo.png" align="right" alt="" width="120">

<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/wateRinfo/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/wateRinfo/actions) [![codecov](https://codecov.io/gh/ropensci/wateRinfo/branch/master/graph/badge.svg?token=4MT1hZyw2l)](https://codecov.io/gh/ropensci/wateRinfo)
<!-- badges: end -->

wateRinfo facilitates access to [waterinfo.be](https://www.waterinfo.be/), a website managed by the [Flanders Environment Agency (VMM)](https://en.vmm.be/) and [Flanders Hydraulics Research](https://www.waterbouwkundiglaboratorium.be/). The website provides access to real-time water and weather related environmental variables for Flanders (Belgium), such as rainfall, air pressure, discharge, and water level. The package provides functions to search for stations and variables, and download time series.

To get started, see:

* [Function reference](https://docs.ropensci.org/wateRinfo/reference/index.html): an overview of all wateRinfo functions.
* [Articles](https://docs.ropensci.org/wateRinfo/articles/): tutorials on how to use the package.

## Installation

You can install wateRinfo from [GitHub](https://github.com/ropensci/wateRinfo) with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/wateRinfo")
```

When successful, load it as usual:

```{r load_lib}
library(wateRinfo)
```

## Example

For a number of supported variables ([documented](https://www.waterinfo.be/download/9f5ee0c9-dafa-46de-958b-7cac46eb8c23?dl=0) by VMM), the stations providing time series data for a given variable can be listed with the command `get_stations()`.

If you want to know the supported variables, ask for the supported variables:

```{r var_listing}
supported_variables("en")
```

Listing the available air pressure stations:

```{r station_listing}
get_stations("air_pressure")
```

Each of the stations in the list for a given variable, are represented by a `ts_id`. These can be used to download the data of a given period with the command `get_timeseries_tsid()`, for example Overpelt (`ts_id = 78124042`):

```{r download_pressure}
overpelt_pressure <- get_timeseries_tsid("78124042", 
                                         from = "2017-04-01", 
                                         to = "2017-04-02")
head(overpelt_pressure)
```

Making a plot of the data with [`ggplot2`](https://ggplot2.tidyverse.org/):

```{r plot_pressure, echo = TRUE, eval = FALSE}
library(ggplot2)
ggplot(overpelt_pressure, aes(x = Timestamp, y = Value)) + 
    geom_line() + 
    xlab("") + ylab("hPa") + 
    scale_x_datetime(date_labels = "%H:%M\n%Y-%m-%d", date_breaks = "6 hours")
```

```{r showplot1, fig.width = 7, echo = FALSE, out.width = '80%'}
knitr::include_graphics("man/figures/README-plot_pressure-1.png")
```

Another option is to check the available variables for a given station, with the function `get_variables()`. Let's consider again Overpelt (`ME11_002`) and check the first ten available variables at the Overpelt measurement station:

```{r var_overpelt, message = FALSE}
vars_overpelt <- get_variables("ME11_002")
head(vars_overpelt, 10)
```

Different pre-calculated variables are already available and a `ts_id` value is available for each of them to download the corresponding data. For example, `DagGem` (= daily mean values) of `RH` (= relative humidity), i.e. `ts_id = 78382042`:

```{r download_rh}
overpelt_rh_daily <- get_timeseries_tsid("78382042", 
                                         from = "2017-04-01", 
                                         to = "2017-04-30")
head(overpelt_rh_daily)
```

```{r plot_rh, echo = TRUE, eval = FALSE}
ggplot(overpelt_rh_daily, aes(x = Timestamp, y = Value)) + 
    geom_line() + 
    xlab("") + ylab(" RH (%)") + 
    scale_x_datetime(date_labels = "%b-%d\n%Y", date_breaks = "5 days")
```

```{r showplot2, fig.width = 7, echo = FALSE, out.width = '80%'}
knitr::include_graphics("man/figures/README-plot_rh-1.png")
```

Unfortunately, not all variables are documented, for which the check for the appropriate variable is not (yet) fully supported by the package.

More detailed tutorials are available in the package vignettes!

## Note on restrictions of the downloads

The amount of data downloaded from waterinfo.be is limited via a credit system. You do not need to get a token right away to download data. For limited and irregular downloads, a token will not be required.

When you require more extended data requests, please request a download token from the waterinfo.be site administrators via the e-mail address <hydrometrie@waterinfo.be> with a statement of which data and how frequently you would like to download data. You will then receive a client-credit code that can be used to obtain a token that is valid for 24 hours, after which the token can be refreshed with the same client-credit code.

Get token with client-credit code: (limited client-credit code for testing purposes)

```{r token_receive}
client <- paste0("MzJkY2VlY2UtODI2Yy00Yjk4LTljMmQtYjE2OTc4ZjBjYTZhOjRhZGE4",
                 "NzFhLTk1MjgtNGI0ZC1iZmQ1LWI1NzBjZThmNGQyZA==")
my_token <- get_token(client = client)
print(my_token)
```

Receive information on the validity of the token:

```{r token_expired}
is.expired(my_token)
```

Check when the token expires:

```{r token_expires_when}
expires.in(my_token)
```

Use token when retrieving data:

```{r use_expires_when}
get_stations(variable_name = "verdamping_monteith", token = my_token)
```

## Other clients

Besides this wateRinfo R client to gather data from [waterinfo.be](https://www.waterinfo.be/), there is also a Python client available. The [pywaterinfo](https://fluves.github.io/pywaterinfo/) package contains similar functionalities.

The [Flanders Hydraulics Research center](https://www.waterbouwkundiglaboratorium.be/en/) also distributes clients for R, Python and Matlab upon request to download the data they share on [waterinfo.be](https://www.waterinfo.be/). For more information, contact them directly via [hic@vlaanderen.be](mailto:hic@vlaanderen.be).

## Acknowledgements

This package is just a small wrapper around waterinfo.be to facilitate researchers and other stakeholders in downloading the data from [waterinfo.be](http://www.waterinfo.be). The availability of this data is made possible by *de Vlaamse Milieumaatschappij, Waterbouwkundig Laboratorium, Maritieme Dienstverlening & Kust, Waterwegen en Zeekanaal NV en De Scheepvaart NV*.

## Meta

* We welcome [contributions](.github/CONTRIBUTING.md) including bug reports.
* License: MIT
* Get citation information for `wateRinfo` in R doing `citation("wateRinfo")`.
* Please note that this project is released with a [Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Define the date period to download"
author: "Stijn Van Hoey"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Define the date period to download}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
library(dplyr)
library(ggplot2)
```

## Introduction

When downloading time series with the `get_timeseries_tsid()` method, the `ts_id` argument provides the link with the variable, location and frequency of the time series, but not the extent/period to download.

The time period to download is defined by a combination of the arguments `from`, `to` and `period`. The usage is similar with the [VMM documentation](https://www.waterinfo.be/download/9f5ee0c9-dafa-46de-958b-7cac46eb8c23?dl=0) for the API itself. The main difference is that the `wateRinfo` package uses existing R functions to interpret the date strings given by the user before sending these to the API (as a formatted string according to `%Y-%m-%d %H:%M:%S`).

This vignette aims to briefly explain how to define the arguments. 

## Which combinations?

In order to define a period, a start and end date is required. Defining all three will result in an error, but any combination of `from/to`, `from/period` and `to/period` is allowed. Moreover, if only `period` or `from` are defined, the waterinfo.be API will automatically define `to` as the current time. Hence, defining *the last x days/months/years/...* can be achieved by only using the `period` option.

## How to define the from/to dates

The package will both except valid date strings as well as valid date objects (`POSIXct`, `POSIXt`) as input for the `from` and `to` arguments. When using a string value, it can be defined on different resolutions:

* "2017-01-01 11:00:00"
* "2017-01-01"
* "2017-01"
* "2017"

According to the [`lubridate`](https://lubridate.tidyverse.org/) package, these orders are accepted: `ymd_hms`, `ymd`, `ym`, `y`. As a result, also `"2017/01/01"`, `"2017 01 01"` or `"20170101"` are valid date string inputs. Make sure the order of year-month-day is respected. For example, `"01/01/2017"`, `"01-01-2017"` and `"01-2017"` are NOT valid. 

## How to define the period

The period string provides a flexible way to extract a time period starting (in combination with `from`) or ending (in combination with `to`) at a given moment. Moreover, by using only the `period` as argument, it will cover all cases where one is interested in *the last x days/months/years/...*.

Some examples are:

* `P3D` : period of three days
* `P2Y` : period of 2 years
* `PT6H` : period of 6 hours
* `P2DT6H` : period of 2 days and 6 hours
* ...

In general, the period string should be provided as `P#Y#M#DT#H#M#S`, where P defines `Period` (always required!) and each # is an integer value expressing *the number of...*. The codes define a specific time interval:

* `Y` - years
* `M` - months
* `D` - days
* `W` - weeks
* `H` - hours
* `M` - minutes
* `S` - seconds

`T` is required if codes about sub-day resolution (day, minutes, hours) is part of the period string. Furthermore, `D` and `W` are mutually exclusive.

More examples of valid period strings are:

* `P1DT12H` : period of 1 day and 12 hours
* `P2WT12H` : period of 2 weeks and 12 hours
* `P1Y6M3DT4H20M30S`: period of 1 year, six months, 3 days, 4 hours, 20 minutes and 30 seconds

## Examples

```{r loadlibrary, warning = FALSE}
library(wateRinfo)
```

When interested in irradiance (15min frequency) data, the following stations provide time series:

```{r irr_stats}
get_stations("irradiance")
```

Focusing on the data of Herentals, the `ts_id` to use is `78930042`. We have different options to define the period to get data from:

1. data about **the last day**, using `period` only:

```{r lastday, fig.width = 7}
irr_lastday <- get_timeseries_tsid("78930042", period = "P1D")
ggplot(irr_lastday, aes(Timestamp, Value)) +
    geom_line() + xlab("") + ylab("irradiance (W/m2)")
```

2. data about **the last 12 hours, 30 minutes**, using `period` only:

```{r lasthours, fig.width = 7}
irr_lasthours <- get_timeseries_tsid("78930042", period = "PT12H30M")
ggplot(irr_lasthours, aes(Timestamp, Value)) +
    geom_line() + xlab("") + ylab("irradiance (W/m2)")
```

3. historical data **from July till August 2014**, using `from` and `to` on month level

```{r historic, fig.width = 7}
irr_2014 <- get_timeseries_tsid("78930042", 
                                from = "2014-07-01", 
                                to = "2014-08-01")
ggplot(irr_2014, aes(Timestamp, Value)) +
    geom_line() + xlab("") + ylab("irradiance (W/m2)")
```

4. historical data for **one day from July 1st 2014**, using `from` and `period`

```{r day2014, fig.width = 7}
irr_2014day <- get_timeseries_tsid("78930042", 
                                from = "2014-07-01", 
                                period = "P1D")
ggplot(irr_2014day, aes(Timestamp, Value)) +
    geom_line() + xlab("") + ylab("irradiance (W/m2)")
```
---
title: "Download time series from multiple stations/variables"
author: "Stijn Van Hoey"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Download time series from multiple stations/variables}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

In many studies, the interest of the user is to download a batch of time series following on a selection criterion. Examples are:

* downloading air pressure data for the last day for all available measurement stations. 
* downloading all measured variables at a frequency of 15 minutes for a given measurement station.

In this vignette, this type of batch downloads is explained, using the available functions of the `wateRinfo` package in combination with already existing [tidyverse](https://www.tidyverse.org/) functionalities.

```{r load_libraries, message = FALSE, warning = FALSE}
library(dplyr)
library(ggplot2)
```

## Download all stations for a given variable

Consider the scenario: "downloading air pressure data for the last day for all available measurement stations". We can achieve this by downloading all the stations information providing air_pressure data (`get_stations()`) and for each of the `ts_id` values in the resulting data.frame, applying the `get_timeseries_tsid()` function:

```{r station_of_var, eval = FALSE}
# extract the available stations for a predefined variable
variable_of_interest <- "air_pressure"
stations <- get_stations(variable_of_interest)

# Download the data for a given period for each of the stations
air_pressure <- stations %>%
    group_by(ts_id) %>%
    do(get_timeseries_tsid(.$ts_id, period = "P1D", to = "2017-01-02")) %>%
    ungroup() %>%
    left_join(stations, by = "ts_id")
```

```{r load_saved_data, echo = FALSE}
air_pressure <- wateRinfo::air_pressure
```

As this results in a tidy data set, we can use the power of ggplot to plot the data of the individual measurement stations:

```{r plot_pressure, fig.width = 7, fig.height = 6}
# create a plot of the individual datasets
air_pressure %>% 
    ggplot(aes(x = Timestamp, y = Value)) + 
    geom_point() + xlab("1 Jan 2017") + 
    facet_wrap(c("station_name", "stationparameter_name")) + 
    scale_x_datetime(date_labels = "%H:%M",
                     date_breaks = "6 hours")
```

## Download set of variables from a station

Consider the scenario: "downloading all soil_moisture (in dutch: 'bodemvocht') variables at a frequency of 15 minutes for the measurement station Liedekerke". We can achieve this by downloading all the variables information of the Liedekerke station(`get_variables()`) using the station code of the waterinfo.be interface (`ME07_006`), filtering on the `P.15` time series and for each of the `ts_id` values, applying the `get_timeseries_tsid()` function:

```{r var_of_station, eval = FALSE}
liedekerke_stat <- "ME07_006"
variables <- get_variables(liedekerke_stat)

variables_to_download <- variables %>% 
    filter(parametertype_name == "Bodemvocht") %>%
    filter(ts_name == "P.15")

liedekerke <- variables_to_download %>%
    group_by(ts_id) %>%
    do(get_timeseries_tsid(.$ts_id, period = "P1M", from = "2017-01-01")) %>%
    ungroup() %>%
    left_join(variables, by = "ts_id")
```

```{r load_saved_data_liedekerke, echo = FALSE}
liedekerke <- wateRinfo::liedekerke
```

As this results in a tidy data set, we can use the power of ggplot to plot the data of the individual measurement stations:

```{r liedekerke_viz, fig.width = 7, fig.height = 6}
liedekerke %>% 
    ggplot(aes(x = Timestamp, y = Value)) + 
    geom_line() + xlab("") + ylab("bodemvocht") + 
    facet_wrap(c("ts_name", "stationparameter_name"), scales = "free") +
    scale_x_datetime(date_labels = "%d-%m\n%Y",
                     date_breaks = "10 days")    
```
---
title: "Introduction to downloading time series data from waterinfo.be"
author: "Stijn Van Hoey"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to downloading time series data from waterinfo.be}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
library(dplyr)
library(ggplot2)
```

## Introduction

The [waterinfo.be](https://www.waterinfo.be) API uses a system of identifiers, called `ts_id` to define individual time series. For example, the identifier `ts_id = 78073042` corresponds to the time series of air pressure data for the measurement station in Liedekerke, with a 15 min time resolution. Hence, the `ts_id` identifier defines a variable of interest from a measurement station of interest with a specific frequency (e.g. 15 min, hourly,...). The knowledge of the proper identifier is essential to be able to download the corresponding data.

```{r load_lib, echo = FALSE}
library(wateRinfo)
```

## Download with known `ts` identifier

In case you already know the `ts_id` identifier that defines your time serie, the package provides the function `get_timeseries_tsid()` to download a specific period of the time series. 

As an example, to download the air pressure time series data of Liedekerke with a 15 min resolution (`ts_id = 78073042`) for the first of January 2016:

```{r tsdownload}
my_data <- get_timeseries_tsid("78073042", from = "2016-01-01", to = "2016-01-02")
knitr::kable(head(my_data), align = "lcc")
```

For more information on defining the date period of the download, see [this vignette](https://ropensci.github.io/wateRinfo/articles/define_date_periods.html). Let's have a visual check of our data, using the `ggplot2` package:

```{r tsvisual, fig.width = 7}
ggplot(my_data, aes(Timestamp, Value)) + 
    geom_line()
```

As such, knowing the identifier is the most straightforward way of downloading a time series. In order to find these `ts_id` identifier, the package supports looking for identifiers based on a supported variable name (limited set of supported variables by VMM) or looking for identifiers by checking all variables for an individual station. These methods are explained in the next sections. 

## Search identifier based on variable name

For a number of variables, the [documentation of the waterinfo.be API](https://www.waterinfo.be/download/9f5ee0c9-dafa-46de-958b-7cac46eb8c23?dl=0) provides a direct overview option of all available VMM measurement stations, using the so-called `Timeseriesgroup_id`. For these variables, the package provides the function `get_stations()` to download an overview of available measurement stations and the related `ts_id` identifiers. The latter can be used to [download the time series](#download-with-known-ts-identifier). 

```{r statofvar}
get_stations("air_pressure")
```

By default, the expected frequency is the 15 min frequency of the time series. However, for some of the variables, multiple frequencies are supported by the API. The package provides a check on the supported variables and frequencies. An overview of the currently supported variables can be requested with the command `supported_variables()` (either in Dutch, `nl`, or in English, `en`). Actually, more variables are available with the API (see next section), but for each of these variables the `get_stations()` function is supported (i.e. the `Timeseriesgroup_id` is documented by VMM). 

```{r suppvar}
supported_variables("en") %>% 
    as.list()
```

To check which predefined frequencies are provided by the waterinfo.be API for a given variable, the `supported_frequencies()` function is available:

```{r suppfreq1}
supported_frequencies(variable_name = "air_pressure")
```

Hence, for air pressure data, only the 15 min resolution is supported. Compared to evaporation derived by the [Monteith equation](https://en.wikipedia.org/wiki/Penman%E2%80%93Monteith_equation):

```{r suppfreq2}
supported_frequencies(variable_name = "evaporation_monteith")
```

Multiple resolutions are available. Using the coarser time resolutions can be helpful when you want to download longer time series while keeping the number of records to download low (if the frequency would be sufficient for your analysis):

```{r statofvarwithfreq}
stations <- get_stations("evaporation_monteith", frequency = "year")
subset_of_columns <- stations %>% select(ts_id, station_no, station_name, 
                                         parametertype_name, ts_unitsymbol)
knitr::kable(subset_of_columns)
```

When interested in the data of `Herentals_ME`, we can use the corresponding `ts_id` to download the time series of PET with a yearly frequency and make a plot with `ggplot`:

```{r yearly_pet, fig.width = 7}
pet_yearly <- get_timeseries_tsid("94526042", period = "P10Y")
pet_yearly %>% 
    na.omit() %>%
    ggplot(aes(Timestamp, Value)) + 
    geom_bar(stat = "identity") +
    scale_x_datetime(date_labels = "%Y", date_breaks = "1 year") + 
    xlab("") + ylab("PET Herentals (mm)")
```

See [this vignette](https://ropensci.github.io/wateRinfo/articles/define_date_periods.html) to understand the `period = "P10Y"` format.

**Remark:** the `get_stations()` function only works for those measurement stations belonging to the VMM `meetnet` (network), related to the so-called `datasource = 1`. For [other networks](http://www.waterinfo.be/default.aspx?path=NL/HIC/Recent_toelichting), i.e. `datasource = 4`, the enlisting is not supported. Still, a search for data is provided starting from a given station name, as explained in the next section.

## Search identifier based on station name

In addition to the option to check the measurement stations that can provide data for a given variable, the package provides the function `get_variables()` to get an overview of the available variables for a given station, using the `station_no`. The advantage compared to the `ts_id` is that these `station_no` names are provided by the waterinfo.be website itself when exploring the data. When clicking on a measurement station on the map and checking the time series graph, the `station_no` is provided in the upper left corner in between brackets. 

![Waterinfo.be example printscreen of time series](./waterinfo_screen.png)

So, for the example in the figure, i.e. `station_no = zes42a-1066`, the available time series are retrieved by using the `get_variables()` command:

```{r varforstat}
available_variables <- get_variables("zes42a-1066")
available_variables %>% select(ts_id, station_name, ts_name, 
                               parametertype_name)
```

The available number of variables depends on the measurement station. The representation is not standardized and also depends on the type of [`meetnet`](http://www.waterinfo.be/default.aspx?path=NL/HIC/Recent_toelichting). Nevertheless, one can derive the required `ts_id` from the list when interpreting the field names. Remark that the datasource can be 4 instead of 1 for specific `meetnetten` (networks). The datasource to use is printed when asking the variables for a station. 

In order to download the 10 min time series water level data for the station in Sint-Amands tij/Zeeschelde, the `ts_id = 55419010` can be used in the `get_timeseries_tsid()` function, taking into account the `datasource = 4` (default is 1):

```{r downloadstamands, fig.width = 7}
tide_stamands <- get_timeseries_tsid("55419010", 
                                     from = "2017-06-01", to = "2017-06-05",
                                     datasource = 4)
ggplot(tide_stamands, aes(Timestamp, Value)) + 
    geom_line() + xlab("") + ylab("waterlevel")
```

For some measurement stations, the number of variables can be high (lots of precalculated derivative values) and extracting the required time series identifier is not always straightforward. For example, the dat Etikhove/Schuif/Nederaalbeek (`K06_221`), provides the following number of variables:

```{r check_bekken}
available_variables <- get_variables("K06_221")
nrow(available_variables)
```

As the measured variables at a small time resolution are of most interest, filtering on `P.15` (or `P.1`, `P.60`,...) will help to identify measured time series, for those stations belonging to the `meetnet` of VMM (`datasource = 1`):

```{r check_bekken_filter}
available_variables <- get_variables("K06_221")
available_variables %>% 
    filter(ts_name == "P.15")
```

Loading and visualizing the last day (period `P1D`) of available data for the water level downstream (`Hafw`):

```{r afw_etik_fig, fig.width = 7}
afw_etikhove <- get_timeseries_tsid("22302042", 
                                    period = "P1D",
                                    datasource = 1) # 1 is default

ggplot(afw_etikhove, aes(Timestamp, Value)) + 
    geom_line() + xlab("") + ylab("Volume")
```

We can do similar filtering to check for time series on other stations, for example the Molenbeek in Etikhove:

```{r check_etikhove}
available_variables <- get_variables("L06_347")
available_variables %>%
    filter(ts_name == "P.15")
```

And use the `ts_id` code representing discharge to create a plot of the discharge during a storm in 2010:

```{r disch_etik, fig.width = 7}
etikhove <- get_timeseries_tsid("276494042", 
                                from = "2010-11-09", to = "2010-11-16")

ggplot(etikhove, aes(Timestamp, Value)) + 
    geom_line() + xlab("") + ylab("Q (m3/s)")
```

## Check the URL used to request the data

As each of these data requests is actually a call to the [waterinfo.be](https://www.waterinfo.be/) API, the call can be tested in a browser as well. To retrieve the URL used to request certain data (using `get_variables()`, `get_stations()` or `get_timeseries_tsid()`), check the `comment()` attribute of the returned `data.frame`, for example:

```{r check_url_1}
air_stations <- get_stations("air_pressure")
comment(air_stations)
```

or 

```{r check_url_2}
etikhove <- get_timeseries_tsid("276494042", 
                                from = "2010-11-09", to = "2010-11-16")
comment(etikhove)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_variables.R
\name{get_variables}
\alias{get_variables}
\title{Get list of variables for a given station}
\format{
A data.frame with 6 variables:
\describe{
  \item{station_name}{Official name of the measurement station.}
  \item{station_no}{Station ID as provided on the waterinfo.be website.}
  \item{ts_id}{Unique timeseries identifier to access time series data
  corresponding to a combination of the station, measured variable and
  frequency.}
  \item{ts_name}{Timeseries identifier description name as provided by
  `waterinfo.be`.}
  \item{parametertype_name}{Measured variable description.}
  \item{stationparameter_name}{Station specific variable description.}
}
The URL of the specific request is provided as a comment attribute to the
returned data.frame. Use \code{comment(df)} to get the request URL.
}
\usage{
get_variables(station_no, token = NULL)
}
\arguments{
\item{station_no}{'stations-nummer' as it appears on the download page of
\href{https://www.waterinfo.be/default.aspx?path=NL/Rapporten/Downloaden}{
waterinfo.be}}

\item{token}{token to use with the call (optional, can be retrieved via
\code{\link{get_token}})}
}
\value{
data.frame with the station_name, station_no, ts_id, ts_name and
parametertype_name for each of the variables for this station.
}
\description{
Get list of variables for a given station
}
\examples{
variables_overpelt <- get_variables("ME11_002")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supported_variables.R
\name{supported_variables}
\alias{supported_variables}
\title{VMM supported timeseriesgroups variables}
\usage{
supported_variables(language = "nl")
}
\arguments{
\item{language}{char \code{nl} (dutch) or \code{en} (english) variable names}
}
\value{
data.frame containing the variable names in either english or dutch
}
\description{
Provide list of VMM supported variables in the timeseriesgroupID
in either dutch or english
}
\examples{
# Request supported variables in Dutch
supported_variables("nl")

# Request supported variables in English
supported_variables("en")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resolve_identifiers.R
\name{resolve_timeseriesgroupid}
\alias{resolve_timeseriesgroupid}
\title{Get timeseriesgroupID for a supported variable}
\usage{
resolve_timeseriesgroupid(variable_name, frequency = "15min")
}
\arguments{
\item{variable_name}{valid variable name, supported by VMM API}

\item{frequency}{valid frequency for the given variable}
}
\value{
list containing the \code{timeseriesgroup_id} of the variable
frequency combination
}
\description{
Translate the usage of available variables to the corresponding
timeseriesgroupID, based on the provided lookup table from VMM
}
\details{
Remark that this information is NOT based on a query, but on information
provided by the package itself to make variable names more readable

The lookup table is provided as external data of the package,
see inst/extdata
}
\examples{
resolve_timeseriesgroupid("rainfall", "15min")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_period.R
\name{parse_period}
\alias{parse_period}
\title{Check the from/to/period arguments}
\usage{
parse_period(from = NULL, to = NULL, period = NULL)
}
\arguments{
\item{from}{string representing date of datetime object}

\item{to}{string representing date of datetime object}

\item{period}{input string according to format required by waterinfo}
}
\value{
list with the relevant period/date information
}
\description{
Handle the information of provided date information on the period and provide
feedback to the user. Valid combinations of the arguments are:
from/to, from/period, to/period, period, from
}
\seealso{
check_period_format
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/call_waterinfo.R
\name{print.waterinfo_api}
\alias{print.waterinfo_api}
\title{Custom print function of the API request response}
\usage{
\method{print}{waterinfo_api}(x, ...)
}
\arguments{
\item{x}{waterinfo_api}

\item{...}{args further arguments passed to or from other methods.}
}
\description{
Custom print function of the API request response
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{air_pressure}
\alias{air_pressure}
\title{Air pressure data of January 1st, 2017}
\format{
A data frame with 710 rows and 13 variables:
\describe{
  \item{ts_id}{identifier of the downloaded time serie}
  \item{Timestamp}{datetime}
  \item{Value}{measured value of the variable}
  \item{Quality Code}{Quality code of the measurement}
  \item{station_latitude}{latitude coordinate}
  \item{station_longitude}{longitude coordinate}
  \item{station_id}{identifier of the measurement station}
  \item{station_no}{short code name of the measurement station}
  \item{station_name}{full name of the measurement station}
  \item{stationparameter_name}{parameter name on station level}
  \item{parametertype_name}{parameter type name}
  \item{ts_unitsymbol}{unit of the variable}
  \item{dataprovider}{provider of the time series value}
}
}
\source{
\url{https://www.waterinfo.be/}
}
\usage{
air_pressure
}
\description{
A dataset compiled by downloading 1 day of air pressure data for the
available stations of Waterinfo.be
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supported_variables.R
\name{is_supported_variable}
\alias{is_supported_variable}
\title{Check if variable is supported by VMM ts group id}
\usage{
is_supported_variable(variable_name)
}
\arguments{
\item{variable_name}{char}
}
\value{
Raise error when variable is not supported directly, otherwise NULL
}
\description{
Check if variable is supported by VMM ts group id
}
\examples{
is_supported_variable("wind_speed")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_period.R
\name{check_date_format}
\alias{check_date_format}
\title{Check if the date can be parsed to a datetime object in R}
\usage{
check_date_format(datetime)
}
\arguments{
\item{datetime}{string representation of the date}
}
\value{
POSIXct date-time object is date is valid representation
}
\description{
if the date is already a datetime object ("POSIXct" "POSIXt"), the object
itself is returned
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resolve_identifiers.R
\name{resolve_datasource}
\alias{resolve_datasource}
\title{Define the datasource using the station number}
\usage{
resolve_datasource(station_no)
}
\arguments{
\item{station_no}{'stations-nummer' as it appears on the download page of
\href{https://www.waterinfo.be/default.aspx?path=NL/Rapporten/Downloaden}{
waterinfo.be}}
}
\value{
integer 1 for VMM, 4 for other 'meetnetten' (HIC,...)
}
\description{
Using the 'stations-nummer' as provided on
\href{https://www.waterinfo.be/default.aspx?path=NL/Rapporten/Downloaden}{
waterinfo.be}, this function tries to identify the datasource to use for
the particular variable
}
\details{
Notice that VMM did not provide this in the official documentation, but this
has just been derived by checking the API response as such. A more automated
and less hard-coded approach would be beneficial, but this data is not
available at the moment.
}
\examples{
resolve_datasource('akl03e-1066')
resolve_datasource('K07_OM421')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_timeseries.R
\name{get_timeseries_tsid}
\alias{get_timeseries_tsid}
\title{Download timeseries data from waterinfo.be}
\format{
A data.frame with 3 variables:
\describe{
  \item{Timestamp}{Datetime of the measurement.}
  \item{Value}{Measured value.}
  \item{Quality Code}{Quality code of the measurement, dependent on the
  data source used:
     \itemize{
        \item{VMM Quality Code Interpretation (datasource 1)
           \itemize{
              \item{10/110 - Excellent}
              \item{30/100/130 - Good}
              \item{50/150 - Moderate}
              \item{70/170 - Poor}
              \item{80/180 - Estimated}
              \item{90/190 - Suspect}
              \item{220 - Default}
              \item{-1 - Missing}
           }
        }
        \item{HIC Quality Code Interpretation (datasource 2)
           \itemize{
              \item{40 - Good}
              \item{80 - Estimated}
              \item{120 - Suspect}
              \item{200 - Unchecked}
              \item{60 - Complete}
              \item{160 - Incomplete}
              \item{-1 - Missing}
           }
        }
        \item{Aggregated timeseries
           \itemize{
              \item{40 - Good}
              \item{100 - Estimated}
              \item{120 - Suspect}
              \item{200 - Unchecked}
              \item{-1 - Missing}
           }
        }
     }
   }
}
The URL of the specific request is provided as a comment attribute to the
returned data.frame. Use \code{comment(df)} to get the request URL.
}
\usage{
get_timeseries_tsid(
  ts_id,
  period = NULL,
  from = NULL,
  to = NULL,
  datasource = 1,
  token = NULL
)
}
\arguments{
\item{ts_id}{waterinfo.be database ts_id, defining a timeserie variable and
frequency it is defined.}

\item{period}{input string according to format required by waterinfo:
De period string is provided as P#Y#M#DT#H#M#S, with P defines `Period`,
each # is an integer value and the codes define the number of...
Y - years
M - months
D - days
T required if information about sub-day resolution is present
H - hours
D - days
M - minutes
S - seconds
Instead of D (days), the usage of W - weeks is possible as well
Examples of valid period strings: P3D, P1Y, P1DT12H, PT6H, P1Y6M3DT4H20M30S.}

\item{from}{date of datestring as start of the time series}

\item{to}{date of datestring as end of the time series}

\item{datasource}{int [0-4] defines the `meetnet` of which the measurement
station is part of. VMM based stations are net '1', MOW-HIC is net '2'}

\item{token}{token to use with the call (optional, can be retrieved via
\code{\link{get_token}})}
}
\value{
data.frame with the timestamps, values and quality code
}
\description{
Using the ts_id codes  and by providing a given date period, download the
corresponding time series from the waterinfo.be website
}
\examples{
get_timeseries_tsid("35055042", from = "2017-01-01", to = "2017-01-02")
get_timeseries_tsid("5156042", period = "P3D")
get_timeseries_tsid("55419010", from = "2017-06-01", to = "2017-06-03",
                    datasource = 4)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_stations.R
\name{get_stations}
\alias{get_stations}
\title{Get list of stations for a variable}
\format{
A data.frame with 10 variables:
\describe{
  \item{ts_id}{Unique timeseries identifier to access time series data
  corresponding to a combination of the station, measured variable and
  frequency.}
  \item{station_latitude}{Latitude coordinates of the station (WGS84)}
  \item{station_longitude}{Longitude coordinates of the station (WGS84)}
  \item{station_id}{Identifier of the station as used in the waterinfo
  backend}
  \item{station_no}{Station ID as provided on the waterinfo.be website.}
  \item{station_name}{Official name of the measurement station.}
  \item{stationparameter_name}{Station specific variable description.}
  \item{parametertype_name}{Measured variable description.}
  \item{ts_unitsymbol}{Unit of the variable.}
  \item{dataprovider}{Data provider of the time series data.}
}
The URL of the specific request is provided as a comment attribute to the
returned data.frame. Use \code{comment(df)} to get the request URL.
}
\usage{
get_stations(variable_name = NULL, frequency = "15min", token = NULL)
}
\arguments{
\item{variable_name}{char valid nam of available variable as timeseriesgroup}

\item{frequency}{char valid frequency for the given variable, for most
variables, the 15min frequency is available}

\item{token}{token to use with the call (optional, can be retrieved via
\code{\link{get_token}})}
}
\value{
data.frame with an overview of the available stations for the
requested variable
}
\description{
For a given timeseriesgroup (variable), provide a list of measurement
stations providing data. An overview of the variables is provided by the
function \code{\link{supported_variables}}.
}
\details{
For the moment, this only works for measurement stations of VMM (meetnet 1),
and stations from other measurement data sources are not included in the list
}
\examples{
get_stations('irradiance')
get_stations('soil_saturation')
}
\seealso{
supported_variables
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_period.R
\name{check_period_format}
\alias{check_period_format}
\title{Check period string format}
\usage{
check_period_format(period_string)
}
\arguments{
\item{period_string}{input string according to format required by waterinfo:
The period string is provided as P#Y#M#DT#H#M#S, with P defines `Period`,
each # is an integer value and the codes define the number of...
Y - years
M - months
D - days
T required if information about sub-day resolution is present
H - hours
D - days
M - minutes
S - seconds
Instead of D (days), the usage of W - weeks is possible as well
Examples of valid period strings: P3D, P1Y, P1DT12H, PT6H, P1Y6M3DT4H20M30S.}
}
\value{
str period string itself if valid
}
\description{
Check if the format of the period is conform the specifications of VMM
}
\examples{
check_period_format("P2DT6H") # period of 2 days and 6 hours
check_period_format("P3D") # period of 3 days
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/call_waterinfo.R
\name{call_waterinfo}
\alias{call_waterinfo}
\title{http call to waterinfo.be}
\usage{
call_waterinfo(query, base_url = "vmm", token = NULL)
}
\arguments{
\item{query}{list of query options to be used together with the base string}

\item{base_url}{str vmm | hic | pro, default download defined}

\item{token}{token to use with the call (optional, can be retrieved via
\code{\link{get_token}})}
}
\value{
waterinfo_api class object with content and info about call
}
\description{
General call used to request information and data from waterinfo.be,
providing error handling and json parsing
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_period.R
\name{isdatetime}
\alias{isdatetime}
\title{Check if the string input can be converted to a date, provides FALSE or date}
\usage{
isdatetime(datetime)
}
\arguments{
\item{datetime}{string representation of a date}
}
\value{
FALSE | "POSIXct" "POSIXt"
}
\description{
(acknowledgements to micstr/isdate.R)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/token.R
\name{get_token}
\alias{get_token}
\alias{token}
\alias{print.token}
\alias{show.token}
\alias{expires.in}
\alias{expires.in.token}
\alias{is.expired}
\alias{is.expired.token}
\title{Get waterinfo Token}
\usage{
get_token(
  client = NULL,
  client_id = NULL,
  client_secret = NULL,
  token_url = "http://download.waterinfo.be/kiwis-auth/token"
)

is.expired(token)

expires.in(token)
}
\arguments{
\item{client}{base64 encoded client containing the client id and client
secret, seperated by :}

\item{client_id}{client id string}

\item{client_secret}{client secret string}

\item{token_url}{url to get the token from}

\item{token}{a token object}
}
\value{
An object of class token containing the token string with the
token_url, token type and moment of expiration as attributes.
}
\description{
Retrieve a fresh waterinfo token. A token is not required to get started,
see Details section for more information.
}
\details{
Notice you do not need to get a token right away to download data.
For limited and irregular downloads, a token will not be required. The amount
of data downloaded from waterinfo.be is limited via a credit system.
When you require more extended data requests, request a download token.

Either client or client_id and client_secret need to be passed as
arguments. If provided, client is always used. Tokens remain valid for
24 hours, after which a fresh one must be acquired.
To limit load on the server, token objects should be reused as much as
possible until expiration in stead of creating fresh ones for each call.

The client_id and client_secret provided in the examples are for test
purposes, get your very own client via \email{hydrometrie@waterinfo.be}.
}
\examples{
# Get token via client_id and client_secret
client_id <- '32dceece-826c-4b98-9c2d-b16978f0ca6a'
client_secret <- '4ada871a-9528-4b4d-bfd5-b570ce8f4d2d'
my_token <- get_token(client_id = client_id,client_secret = client_secret)
print(my_token)

# get token via client
client <- paste0('MzJkY2VlY2UtODI2Yy00Yjk4LTljMmQtYjE2OTc4ZjBjYTZhOjRhZGE4',
                'NzFhLTk1MjgtNGI0ZC1iZmQ1LWI1NzBjZThmNGQyZA==')
my_token <- get_token(client = client)
print(my_token)
is.expired(my_token)
expires.in

# Use the token when requesting for data (i.e. get_* functions), e.g.
get_stations(variable_name = "verdamping_monteith", token = my_token)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supported_variables.R
\name{supported_frequencies}
\alias{supported_frequencies}
\title{VMM supported timeseriesgroups frequencies}
\usage{
supported_frequencies(variable_name)
}
\arguments{
\item{variable_name}{char name of a valid variable in either dutch or english}
}
\description{
Provide list of VMM supported frequencies for a given timeseriesgroupID
in either dutch or english
}
\examples{
supported_frequencies('rainfall')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{liedekerke}
\alias{liedekerke}
\title{Soil moisture data of Liedekerke, January 2017}
\format{
A data frame with 23,816 rows and 9 variables:
\describe{
  \item{ts_id}{identifier of the downloaded time serie}
  \item{Timestamp}{datetime}
  \item{Value}{measured value of the variable}
  \item{Quality Code}{Quality code of the measurement}
  \item{station_name}{full name of the measurement station}
  \item{station_no}{short code name of the measurement station}
  \item{ts_name}{type/frequency of the time serie}
  \item{parametertype_name}{parameter type name}
  \item{stationparameter_name}{parameter name on station level}
}
}
\source{
\url{https://www.waterinfo.be/}
}
\usage{
liedekerke
}
\description{
A dataset compiled by downloading 1 day of soil moisture data for the
Liedekerke  measurement station of Waterinfo.be
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wateRinfo.R
\docType{package}
\name{wateRinfo-package}
\alias{wateRinfo}
\alias{wateRinfo-package}
\title{wateRinfo: Download Time Series Data from Waterinfo.be}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

wateRinfo facilitates access to waterinfo.be (<https://www.waterinfo.be>), a website managed by the Flanders Environment Agency (VMM) and Flanders Hydraulics Research. The website provides access to real-time water and weather related environmental variables for Flanders (Belgium), such as rainfall, air pressure, discharge, and water level. The package provides functions to search for stations and variables, and download time series.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/wateRinfo}
  \item \url{https://docs.ropensci.org/wateRinfo}
  \item Report bugs at \url{https://github.com/ropensci/wateRinfo/issues}
}

}
\author{
\strong{Maintainer}: Stijn Van Hoey \email{stijnvanhoey@gmail.com} (\href{https://orcid.org/0000-0001-6413-3185}{ORCID})

Other contributors:
\itemize{
  \item Willem Maetens \email{w.maetens@vmm.be} [contributor]
  \item Peter Desmet \email{peter.desmet@inbo.be} (\href{https://orcid.org/0000-0002-8442-8025}{ORCID}) [contributor]
  \item Research Institute for Nature and Forest (INBO) \email{info@inbo.be} [copyright holder]
}

}
\keyword{internal}
