
# Subnational data for the COVID-19 outbreak <img src="man/figures/logo.png" align="right" alt="" width="120" />

[![R-CMD-check](https://github.com/epiforecasts/covidregionaldata/workflows/R-CMD-check/badge.svg)](https://github.com/epiforecasts/covidregionaldata/actions)
[![Codecov test
coverage](https://codecov.io/gh/epiforecasts/covidregionaldata/branch/master/graph/badge.svg)](https://codecov.io/gh/epiforecasts/covidregionaldata?branch=master)
[![Data
status](https://img.shields.io/badge/Data-status-lightblue.svg?style=flat)](https://epiforecasts.io/covidregionaldata/articles/supported-countries.html)
[![metacran
downloads](http://cranlogs.r-pkg.org/badges/grand-total/covidregionaldata?color=ff69b4)](https://cran.r-project.org/package=covidregionaldata)

[![MIT
license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/epiforecasts/covidregionaldata/blob/master/LICENSE.md/)
[![GitHub
contributors](https://img.shields.io/github/contributors/epiforecasts/covidregionaldata)](https://github.com/epiforecasts/covidregionaldata/graphs/contributors)
[![Discord](https://img.shields.io/discord/864828485981306890?logo=Discord)](https://discord.gg/9YPDDADVt3)
[![PRs
Welcome](https://img.shields.io/badge/PRs-welcome-yellow.svg)](https://makeapullrequest.com/)
[![GitHub
commits](https://img.shields.io/github/commits-since/epiforecasts/covidregionaldata/0.9.2.svg?color=orange)](https://GitHub.com/epiforecasts/covidregionaldata/commit/master/)

[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03290/status.svg)](https://doi.org/10.21105/joss.03290)
[![Zenodo](https://zenodo.org/badge/271601189.svg)](https://zenodo.org/badge/latestdoi/271601189)

Interface to subnational and national level COVID-19 data sourced from
both official sources, such as Public Health England in the UK, and from
other COVID-19 data collections, including the World Health Organisation
(WHO), European Centre for Disease Prevention and Control (ECDC), John
Hopkins University (JHU), Google Open Data and others. This package is
designed to streamline COVID-19 data extraction, cleaning, and
processing from a range of data sources in an open and transparent way.
This allows users to inspect and scrutinise the data, and tools used to
process it, at every step. For all countries supported, data includes a
daily time-series of cases and, wherever available, data on deaths,
hospitalisations, and tests. National level data is also supported using
a range of data sources as well as line list data and links to
intervention data sets.

## Installation

Install from CRAN:

``` r
install.packages("covidregionaldata")
```

Install the stable development version of the package with:

``` r
install.packages("covidregionaldata",
  repos = "https://epiforecasts.r-universe.dev"
)
```

Install the unstable development version of the package with:

``` r
remotes::install_github("epiforecasts/covidregionaldata")
```

## Quick start

[![Documentation](https://img.shields.io/badge/Documentation-lightgrey.svg?style=flat)](https://epiforecasts.io/covidregionaldata/)

Load `covidregionaldata`, `dplyr`, `scales`, and `ggplot2` (all used in
this quick start),

``` r
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)
```

### Setup data caching

This package can optionally use a data cache from `memoise` to locally
cache downloads. This can be enabled using the following (this will use
the temporary directory by default),

``` r
start_using_memoise()
#> Using a cache at: /tmp/RtmpPgZXiv
```

To stop using `memoise` use,

``` r
stop_using_memoise()
```

and to reset the cache (required to download new data),

``` r
reset_cache()
```

### National data

To get worldwide time-series data by country (sourced from the World
Health Organisation (WHO) by default but also optionally from the
European Centre for Disease Control (ECDC), John Hopkins University, or
the Google COVID-19 open data project), use:

``` r
nots <- get_national_data()
#> Downloading data from https://covid19.who.int/WHO-COVID-19-global-data.csv
#> Rows: 142911 Columns: 8
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr  (3): Country_code, Country, WHO_region
#> dbl  (4): New_cases, Cumulative_cases, New_deaths, Cumulative_deaths
#> date (1): Date_reported
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Cleaning data
#> Processing data
nots
#> # A tibble: 142,911 × 15
#>    date       un_region who_region country        iso_code cases_new cases_total
#>    <date>     <chr>     <chr>      <chr>          <chr>        <dbl>       <dbl>
#>  1 2020-01-03 Asia      EMRO       Afghanistan    AF               0           0
#>  2 2020-01-03 Europe    EURO       Albania        AL               0           0
#>  3 2020-01-03 Africa    AFRO       Algeria        DZ               0           0
#>  4 2020-01-03 Oceania   WPRO       American Samoa AS               0           0
#>  5 2020-01-03 Europe    EURO       Andorra        AD               0           0
#>  6 2020-01-03 Africa    AFRO       Angola         AO               0           0
#>  7 2020-01-03 Americas  AMRO       Anguilla       AI               0           0
#>  8 2020-01-03 Americas  AMRO       Antigua & Barbuda AG               0           0
#>  9 2020-01-03 Americas  AMRO       Argentina      AR               0           0
#> 10 2020-01-03 Asia      EURO       Armenia        AM               0           0
#> # … with 142,901 more rows, and 8 more variables: deaths_new <dbl>,
#> #   deaths_total <dbl>, recovered_new <dbl>, recovered_total <dbl>,
#> #   hosp_new <dbl>, hosp_total <dbl>, tested_new <dbl>, tested_total <dbl>
```

This can also be filtered for a country of interest,

``` r
g7 <- c(
  "United States", "United Kingdom", "France", "Germany",
  "Italy", "Canada", "Japan"
)
g7_nots <- get_national_data(countries = g7, verbose = FALSE)
```

Using this data we can compare case information between countries, for
example here is the number of deaths over time for each country in the
G7:

``` r
g7_nots %>%
  ggplot() +
  aes(x = date, y = deaths_new, col = country) +
  geom_line(alpha = 0.4) +
  labs(x = "Date", y = "Reported Covid-19 deaths") +
  scale_y_continuous(labels = comma) +
  theme_minimal() +
  theme(legend.position = "top") +
  guides(col = guide_legend(title = "Country"))
```

<img src="man/figures/README-g7_plot-1.png" width="100%" />

### Subnational data

To get time-series data for subnational regions of a specific country,
for example by level 1 region in the UK, use:

``` r
uk_nots <- get_regional_data(country = "UK", verbose = FALSE)
uk_nots
#> # A tibble: 7,501 × 26
#>    date       region   region_code cases_new cases_total deaths_new deaths_total
#>    <date>     <chr>    <chr>           <dbl>       <dbl>      <dbl>        <dbl>
#>  1 2020-01-30 East Mi… E12000004          NA          NA         NA           NA
#>  2 2020-01-30 East of… E12000006          NA          NA         NA           NA
#>  3 2020-01-30 England  E92000001           2           2         NA           NA
#>  4 2020-01-30 London   E12000007          NA          NA         NA           NA
#>  5 2020-01-30 North E… E12000001          NA          NA         NA           NA
#>  6 2020-01-30 North W… E12000002          NA          NA         NA           NA
#>  7 2020-01-30 Norther… N92000002          NA          NA         NA           NA
#>  8 2020-01-30 Scotland S92000003          NA          NA         NA           NA
#>  9 2020-01-30 South E… E12000008          NA          NA         NA           NA
#> 10 2020-01-30 South W… E12000009          NA          NA         NA           NA
#> # … with 7,491 more rows, and 19 more variables: recovered_new <dbl>,
#> #   recovered_total <dbl>, hosp_new <dbl>, hosp_total <dbl>, tested_new <dbl>,
#> #   tested_total <dbl>, areaType <chr>, cumCasesByPublishDate <dbl>,
#> #   cumCasesBySpecimenDate <dbl>, newCasesByPublishDate <dbl>,
#> #   newCasesBySpecimenDate <dbl>, cumDeaths28DaysByDeathDate <dbl>,
#> #   cumDeaths28DaysByPublishDate <dbl>, newDeaths28DaysByDeathDate <dbl>,
#> #   newDeaths28DaysByPublishDate <dbl>, …
```

Now we have the data we can create plots, for example the time-series of
the number of cases for each region:

``` r
uk_nots %>%
  filter(!(region %in% "England")) %>%
  ggplot() +
  aes(x = date, y = cases_new, col = region) +
  geom_line(alpha = 0.4) +
  labs(x = "Date", y = "Reported Covid-19 cases") +
  scale_y_continuous(labels = comma) +
  theme_minimal() +
  theme(legend.position = "top") +
  guides(col = guide_legend(title = "Region"))
```

<img src="man/figures/README-uk_plot-1.png" width="100%" />

See `get_available_datasets()` for supported regions and subregional
levels. To view what datasets we currently have subnationaldata for,
along with their current status, check the [supported
countries](https://epiforecasts.io/covidregionaldata/articles/supported-countries.html)
page or build the [supported countries
vignette](vignettes/supported-countries.Rmd).

For further examples see the [quick start
vignette](https://github.com/epiforecasts/covidregionaldata/blob/master/vignettes/quickstart.Rmd).
Additional subnational data are supported via the `JHU()` and `Google()`
classes. Use the `available_regions()` method once these data have been
downloaded and cleaned (see their examples) for subnational data they
internally support.

## Citation

If using `covidregionaldata` in your work please consider citing it
using the following,

    #> 
    #> To cite covidregionaldata in publications use:
    #> 
    #>   Joseph Palmer, Katharine Sherratt, Richard Martin-Nielsen, Jonnie
    #>   Bevan, Hamish Gibbs, Sebastian Funk and Sam Abbott (2021).
    #>   covidregionaldata: Subnational data for COVID-19 epidemiology, DOI:
    #>   10.21105/joss.03290
    #> 
    #> A BibTeX entry for LaTeX users is
    #> 
    #>   @Article{,
    #>     title = {covidregionaldata: Subnational data for COVID-19 epidemiology},
    #>     author = {Joseph Palmer and Katharine Sherratt and Richard Martin-Nielsen and Jonnie Bevan and Hamish Gibbs and Sebastian Funk and Sam Abbott},
    #>     journal = {Journal of Open Source Software},
    #>     year = {2021},
    #>     volume = {6},
    #>     number = {63},
    #>     pages = {3290},
    #>     doi = {10.21105/joss.03290},
    #>   }

## Development

[![Development](https://img.shields.io/badge/Wiki-lightblue.svg?style=flat)](https://github.com/epiforecasts/covidregionaldata/wiki/)

This package is the result of work from a number of contributors (see
contributors list
[here](https://epiforecasts.io/covidregionaldata/authors.html)). We
would like to thank the [CMMID COVID-19 working
group](https://cmmid.github.io/groups/ncov-group.html) for insightful
comments and feedback.

We welcome contributions and new contributors! We particularly
appreciate help adding new data sources for countries at sub-national
level, or work on priority problems in the
[issues](https://github.com/epiforecasts/covidregionaldata/issues).
Please check and add to the issues, and/or add a [pull
request](https://github.com/epiforecasts/covidregionaldata/pulls). For
more details, start with the [contributing
guide](https://github.com/epiforecasts/covidregionaldata/wiki/Contributing).
For details of the steps required to add support for a dataset see the
[adding data
guide](https://github.com/epiforecasts/covidregionaldata/wiki/Adding-Data).
# covidregionaldata 0.9.3

This release is currently under development

## New data sets

* Support for level 1 region data in Estonia (thanks to @RichardMN). See `?Estonia` for details.
* Support for level 1 region data in Vietnam (thanks to @biocyberman). See `?Vietnam` for details.

## Other changes

* Change the data source for Switzerland to draw data from the Swiss Federal Office of Public Health (FOPH)
* Updated the package logo to include the newly supported data sets.
* Reduced the number of package dependencies (@bisaloo and @RichardMN)
## Bug fixes

- Fixed a bug in the data sourced from Germany so that instead of treating it as a line list of individuals it is treated as a relatively finely resolved count data which needs to be summed up (by @sbfnk).

# covidregionaldata 0.9.2

This release adds support for the Covid19 Data Hub which includes Google and Apple mobility data amongst a large range of other data sets, data from the European Commission's Joint Research Centre which is at both the regional and national level, and individual sources for regional data from several countries. Package updates have been made in line with a software review at the [Journal of Open Source Software](https://github.com/openjournals/joss-reviews/issues/3290). Finally, this release exposes more of the testing infrastructure to users and adds a package hexsticker.

Thanks to @joseph-palmer, @RichardMN, and @kathsherratt for contributions towards this release.

## New features

* Support added for data sets from Covid19 Data Hub. This source aggregates a range of data at a national and subnational level and provides keys to link to mobility data provided by Apple and Google (by @joseph-palmer).
* Support added for data from the European Commission's Joint Research Centre (JRC). The source aggregates incidence data at the country and regional level for 34 UCPM Participating States plus Switzerland (by @joseph-palmer).
* Support added for data from the Netherlands provided by RVIM (English: National Institute for Public Health and the Environment). This source provides case, deaths and hospital admission data at the province and municipal levels (by @joseph-palmer).
* Support added for data from Switzerland and Liechtenstein collated by Canton Zurich (@OpenZH). This source provides case, deaths and hospital admission data at the canton level (by @RichardMN).
* Made package changes recomended in the JOSS review, including additional statements of need to the README, updates to the manuscript (paper.md) and fixes a bug of multiple sources for some countries. We are very grateful for the detailed feedback given by the JOSS reviewers and their help in improving this package.

## Changes to implemented data sources

* Increased the robustness of fetching UK NHS admissions by region. Rather than testing a single date for data we now look over the last 7 days and pick the most recent available data set (by @kathsherratt).

## Other changes

* Testing of classes updated to allow for at least one of `common_data_urls` or `level_data_urls` to be present. The previous default which forced the presence of `common_data_urls` meant that several classes had to define an empty field (by @joseph-palmer).
* Tests on data sets are now included as a method in `DataClass`. `test_regional-datasets` now calls the test function for all classes at each level. Data set specific tests (such as for NHS regions in the UK) are included as a `specific_tests` function within the country class, which is called by the parent (DataClass) `test` after performing standard checks. This allows all the code about a country to be defined in its own class. In addition, users can run tests interactively by calling the test method (e.g. `$test()`) (by @joseph-palmer)
* A function to create a template class and automatically add a github workflow file has been added. This makes adding a new data source for a country even easier as now you can call the function `make_new_data_source()` with the country / source name to add and it will set up the basic structure for you. There is also now a github check to make sure all new sources have a workflow set up (by @joseph-palmer).
* Adds `source_` fields to all data sets to help users properly attribute their data sources (by @RichardMN).

## Bug fixes

* An issue where the `Lithuania()` data set would ignore optional class specific arguments has been fixed (by @RichardMN).
* An issue where the `JHU()` source had multiple region codes for each country has been fixed, giving just one region code per country (by @joseph-palmer).

# covidregionaldata 0.9.1

This release adds support for data sets from John Hopkins University and the Google open data project. Both of these sources aggregate a range of data at national and subnational levels. It also contains a range of small fixes and improvements to documentation. Finally, this release adds optional data processing which will be extended in future releases (contributions warmly welcomed).

Thanks to @joseph-palmer, @RichardMN, and @kathsherratt for contributions towards this release.

## New features

* Support for data provided by John Hopkins University (by @joseph-palmer).
* Support for data provided by Google COVID-19 open data project (by @joseph-palmer).
* Added a `available_regions` method for all classes that shows level 1 regions with data available for the region of interest. This is of particular use when combined with the JHU or Google datasets where processing a large number of regions that are not required can take some time.
* Adds support for JHU or Google data to `get_national_data()`. This may also now be used to access lower level data from these sources  but it may be better to instead use the classes directly or via `initialise_dataclass()`.

## Other changes

* The optional downloading of NHS region data in the `UK()` has been improved to include both the dynamic data previously supported and the archive document now produced (by @kathsherratt).
* The examples for the `UK()` class have been expanded to better showcase the package functionality.
* The documentation and examples for `get_regional_data()`, `get_national_data()`, and `get_available_datasets()` has been expanded with a focus on increasing the visibility of the underlying package structure.
* The documentation and examples for `initialise_dataclass()`, `DataClass()`, and `CountryDataClass()` has been expanded and improved.
* Improvements to the linking of documentation for related functions and classes.
* Improvements to the documentation for contributors (by @RichardMN).
* Improvements to the `pkgdown` documentation to organise packages into appropriate subcategories.

# covidregionaldata 0.9.0

In this release `covidregionaldata` has been substantially retooled to be more robust, and to handle data in a more transparent way. Adding new data sets and functionality has also been made more streamlined. As this update is a substantial package refactor some breaking changes may be been inadvertently introduced. If requiring the old behaviour please install `covidregionaldata@0.8.3` from GitHub.

Thanks to @joseph-palmer, @RichardMN, and @kathsherratt for major contributions to this release. Thanks to @RichardMN for volunteering his time.

## New features

* Track data processing from raw to clean using the `step = TRUE` argument in `get_regional_data()`.
* Filter datasets for regions and countries of interest.
* Access the underlying methods for data sets and all steps in the data processing pipeline.

## Documentation

* All vignettes have been updated for readability.
* A quickstart has been added to the package README.

## Technical improvements

* `get_regional_data()` and `get_national_data()` now use R6 method dispatch. This is an internal change and so should have minimal user impact for users of the `get_` functions. However, all datasets are now available to be used as R6 methods (see `get_available_datasets`) which may allow for more modular use cases. These classes can also be initialised using `initialise_dataclass()` which is used internally by both `get_regional_data()` and `get_national_data()`.
* Unit testing has been separated from data downloading which is now tested individually by data set. This allows for contributors to more easily assess the impact of their additions and also allows us to publish data status updates for each data sets (see the README: <https://github.com/epiforecasts/covidregionaldata#readme>).

## Deprecated functions

* `get_available_datasets()` replaces `get_info_covidregionaldata()` to view available data. `get_info_covidregionaldata()` is deprecated.
* `get_interventions_data()` is deprecated. These data no longer update as of December 2020. Check for alternatives at <https://supertracker.spi.ox.ac.uk/policy-trackers/>
* `get_linelist` is deprecated. Linelist stopped updating June 2020. Up to date linelist data are now behind a login: access at <https://global.health/>. We are working on a solution for accessing with `covidregionaldata`.

## Data changes since 0.8.3

* Colombia now has capitalized region names.
* Germany level 2 region codes have been removed (previously was all NAs).
* India uses NA for unknown region codes, a change from IN-UN previously.
* Italy column region is now regioni.
* Mexico codes 'inegi_code' has been changed to 'inegi'.
* UK Level 1 'ons_region_code' is now 'region_code'.
* UK level 2 "ltla_code" is now "local_authority_code".
* `get_available_datasets()` now return an origin column rather than a country column and a type column rather than a get_data_function to better reflect the types of data supported.

# covidregionaldata 0.8.3

## New data sets

* Level 1 admin data for Cuba
* Level 1 admin data for South Africa

## Data set changes

* UK data - added option to get either lower tier or upper tier local authorities (level 2 regions).
* Updated Northern Ireland case data to be by specimen date by default rather than by date of report as was previously the case. This means that in the UK all data except for data streams from Wales and Scotland are by date of specimen.
* Switched to the WHO source as our default for national level data.

## New features

* Relevant up to date package information can be fetched using `get_info_covidregionaldata()`.
* Switched to using `vroom` for faster `csv` downloads.

## Other changes

* Replaced silently broken functions for converting cumulative data to daily and vice versa.
* Removed integration with Covid19R formatting.
* UK data - removed hospital admissions by age, and occupied mechanical ventilation beds. Currently, these don't fit into the existing data structure and are not available at lower level regions.
* Removed code that required ECDC variables to work.
* Update the ECDC source to pull data from the new weekly snapshots. Updated the variables. In a later update the downloading the now archived daily data will be made possible.

# covidregionaldata 0.8.2

* Updates the API backend used to extract UK data to V2. Adds a release date variable which can be used to return data releases from specified dates rather than the latest snapshot.
* Various fixes to maintain compatibility with data set sources.
* Adds a quickstart vignette with examples of exploratory data analysis.

# covidregionaldata 0.7.0

## Breaking changes

* `get_linelist`: argument `clean` changed to `clean_dates` to reflect slight change in use case.

## Changes

* Added new option to return UK data by NHS region. This will also return "first admissions" hospital data (excludes readmissions). Specify 'nhsregions = TRUE'. Default is FALSE, returning ONS regions as before.
* Fixed inconsistent reference dates for variables in UK data. cases_new and cases_total now by "Specimen date" (date of test), while deaths_new and deaths_total are by "Date of death", for all regions and nations.

* Additional delays added to `get_linelist` when `clean_dates = TRUE`.

# covidregionaldata 0.6.0

* Added whitespace trimming to all regional data functions.
* Fixed region codes for Colombia.
* Fixed region name cleaning for afghanistan.
* Updated UK data source and expanded available variables based on the newly implemented API.
* Enabled regional localisation to be optional.
* Minor quality of life changes.

# covidregionaldata 0.5.0

* Release candidate.

# covidregionaldata 0.4.0

* Added functions to extract regional cases for Spain, Japan, South Korea, the United Kingdom, and the United States.
* Improved data extraction from the ECDC.
* Added minimal unit tests.

# covidregionaldata 0.3.0

* Added a function to extract case counts by region for Italy.
* Function to extract ECDC cases.
* Added a function to extract case counts by region in Germany.
* Fixed cache reset.
* Added a `covidregionaldata` cache of the public linelist as a fall back option if the source is not available.

# covidregionaldata 0.2.0

* Added `memoise` functionality to automatically cache in the current directory all remote resources.
* Added a non-Hubei linelist function that can extract country and city linelists.
* Added a function to extract the latest WHO case counts.
* Added a function to summarise the current total number of cases.

# covidregionaldata 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Added `get_linelist` function
* Added functions to reparameterise Weibull and Gamma distributions with mean and standard deviations.
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at [INSERT CONTACT
METHOD]. All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0,
available at https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.
## Test environments
* local R installation, R 4.1.0
* ubuntu 16.04 (on travis-ci), R 4.1.0
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 note


---
title: "covidregionaldata: Subnational data for COVID-19 epidemiology"
authors:
  - name: Joseph Palmer^[co-first author]
    orcid: 0000-0002-5593-9352
    affiliation: 1
  - name: Katharine Sherratt^[co-first author]
    orcid: 0000-0003-2049-3423
    affiliation: 2
  - name: Richard Martin-Nielsen
    affiliation: 3
  - name: Jonnie Bevan
    affiliation: 4
  - name: Hamish Gibbs
    orcid: 0000-0003-4413-453X
    affiliation: 2
  - name: CMMID COVID-19 Working Group
    affiliation: 2
  - name: Sebastian Funk
    affiliation: 2
  - name: Sam Abbott^[corresponding author]
    affiliation: 2
    orcid: 0000-0001-8057-8037
affiliations:
 - name: Department of Biological Sciences, Royal Holloway University of London
   index: 1
 - name: Centre for Mathematical Modelling of Infectious Diseases, London School of Hygiene & Tropical Medicine
   index: 2
 - name: None
   index: 3
 - name: Tessella
   index: 4
date: "11 May 2021"
bibliography: paper.bib
tags:
  - R
  - COVID-19
  - Open data
  - rstats
  - Sars-Cov-2

output: articles::joss_article
journal: JOSS
link-citations: yes
---

# Summary

`covidregionaldata` is an R [@Rdev:2020] package that provides an interface to subnational and national level COVID-19 data. The package provides cleaned and verified COVID-19 test-positive case counts and, where available, counts of deaths, recoveries, and hospitalisations in a consistent and fully transparent framework. The package automates common processing steps while allowing researchers to easily and transparently trace the origin of the underlying data sources. It has been designed to allow users to easily extend the package's capabilities and contribute to shared data handling. All package code is archived on Zenodo and [GitHub](https://github.com/epiforecasts/covidregionaldata).

# Statement of need

The onset of the COVID-19 pandemic in late 2019 has placed pressure on public health and research communities to generate evidence that can help advise national and international policy in order to reduce transmission and mitigate harm. At the same time, there has been a renewed policy and public health emphasis on localised, subnational decision making and implementation [@Hale2021; @Liu2021]. This requires reliable sources of data disaggregated to a fine spatial scale, ideally with few and/or known sources of bias.

At a national level, epidemiological COVID-19 data is available to download from official sources such as the [World Health Organisation (WHO)](https://covid19.who.int/) [@WorldHealthOrganisation] or the [European Centre for Disease Prevention and Control (ECDC)](https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide%7D) [@EuropeanCentreforDiseasePreventionandControl]. Many government bodies provide a wider range of country specific data, such as [Public Health England in the United Kingdom](https://coronavirus.data.gov.uk/details/about-data) [@PublicHealthEngland], and this is often the only way to access data at a subnational scale, for example by state, district, or province.

Sometimes collated from a range of national and subnational sources, these data come in a variety of formats, requiring users to check and standardise data before it can be combined or processed for analysis.  This is a particularly time-consuming process for subnational data sets, which are often only available in the originating countries’ languages and require customised methods for downloading and processing. This generates potential for errors through programming mistakes, changes to a dependency package, or unexpected changes to a data source. This can lead to misrepresenting the data in ways which are difficult to identify. At best, an independent data processing workflow only slows down the pace of research and analysis, while at worst it can lead to misleading and erroneous results.

Because of these issues, it is important to develop robust tools that provide cleaned, checked and standardised data from multiple sources in a transparent manner. `covidregionaldata` provides easy access to clean data using a single-argument function, ready for analysing the epidemiology of COVID-19 from local to global scales, and in a framework that is easy to trace from raw data to the final standardised data set. Additional arguments to this function support users to, amongst other options, specify the spatial level of subnational data, return data with either standardised or country-specific variable names, or to access the full pipeline from raw to clean data. By default, cleaned and processed data is returned, however, the raw data from a source can also be returned. All data sources are checked daily via GitHub workflows and their status reported in the documentation section 'Data Status'. `covidregionaldata` largely depends on popular packages that many researchers are familiar with (such as the `tidyverse` suite [@Wickham2019]) and can therefore be easily adopted by researchers working in R. In addition to code coverage tests, we test and report the status of all data sets daily.

Currently, `covidregionaldata` provides subnational data collated by official government bodies or by credible non-governmental efforts for 15 countries, including the UK, India, USA, and Brazil. It also provides an interface to subnational data curated by Johns Hopkins University [@Dong2020], and the [Google COVID-19 open data project](https://github.com/GoogleCloudPlatform/covid-19-open-data) [@Wahltinez2020]. National-level data is provided from the World Health Organisation (WHO) [@WorldHealthOrganisation], European Centre for Disease Prevention and Control (ECDC) [@EuropeanCentreforDiseasePreventionandControl], Johns Hopkins University (JHU) [@Dong2020], and the Google COVID-19 open data project [@Wahltinez2020].

# State of the field

Multiple organisations  have built private COVID-19 data curation pipelines similar to that provided in `covidregionaldata`, including Johns Hopkins University (JHU) [@Dong2020], Google [@Wahltinez2020], and the COVID-19 Data Hub [@covid19datahub:2020]. However, most of these efforts aggregate the data they collate into a separate data stream, breaking the linkage with the raw data, and often do not fully surface their data processing pipeline for others to inspect. In contrast `covidregionaldata` provides a clear set of open and fully documented tools that directly operate on raw data where possible in order to make the full data cleaning process transparent to end users.

Other interfaces to COVID-19 data are available in R, though there are fewer that provide tools for downloading subnational data for multiple countries and none that are known to the authors provide a consistent cleaning pipeline of the data sources they support. COVID-19 Data Hub [@covid19datahub:2020] provides cleaning functions, a wrapper to a custom database hosted by COVID-19 Data Hub, and access to snapshots of data reported historically. `Covdata` [@covdata] provides weekly COVID-19 data updates as well as mobility and activity data from [Apple](https://covid19.apple.com/mobility) [@Apple] and [Google](https://www.google.com/covid19/mobility/data_documentation.html) [@Google]. `Sars2pack` [@sars2pack] provides interfaces to a large number of data sets curated by external organisations. To our knowledge, none of these packages provide an interface to individual country data sources or a consistent set of data handling tools for both raw and processed data.

`covidregionaldata` has been used by researchers to source standardised data for estimating the effective reproductive number of COVID-19 in real-time both nationally and subnationally [@Abbott2020]. It has also been used in analyses comparing effective reproduction numbers from different subnational data sources in the United Kingdom [@Sherratt2020], and estimating the increase in transmission related to the B.1.1.7 variant [@Davies2021]. As well as its use in research it has also been used to visualise and explore current trends in COVID-19 case, deaths, and hospitalisations.

# Acknowledgements

This package provides an interface to data sources which are often collected and maintained by individuals or small teams. Our work, both in this package and more generally, would not be possible without their efforts. Thanks to all contributors and package users who have otherwise provided feedback. Thanks to Tim Taylor for useful design discussions.

# Funding statement

This work was supported by a studentship to J.P. funded by the Biotechnology and Biological Sciences Research Council (BBSRC) grant nr. (BB/M011178/1). SEA, KS, and SF were funded by a Wellcome Trust Senior Research Fellowship to Sebastian Funk (210758/Z/18/Z).

# References
