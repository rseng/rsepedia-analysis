
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
---
output: github_document
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# Subnational data for the COVID-19 outbreak <img src="man/figures/logo.png" align="right" alt="" width="120" />

[![R-CMD-check](https://github.com/epiforecasts/covidregionaldata/workflows/R-CMD-check/badge.svg)](https://github.com/epiforecasts/covidregionaldata/actions) [![Codecov test coverage](https://codecov.io/gh/epiforecasts/covidregionaldata/branch/master/graph/badge.svg)](https://codecov.io/gh/epiforecasts/covidregionaldata?branch=master) [![Data status](https://img.shields.io/badge/Data-status-lightblue.svg?style=flat)](https://epiforecasts.io/covidregionaldata/articles/supported-countries.html) [![metacran downloads](http://cranlogs.r-pkg.org/badges/grand-total/covidregionaldata?color=ff69b4)](https://cran.r-project.org/package=covidregionaldata)

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/epiforecasts/covidregionaldata/blob/master/LICENSE.md/)   [![GitHub contributors](https://img.shields.io/github/contributors/epiforecasts/covidregionaldata)](https://github.com/epiforecasts/covidregionaldata/graphs/contributors)  [![Discord](https://img.shields.io/discord/864828485981306890?logo=Discord)](https://discord.gg/9YPDDADVt3)  [![PRs Welcome](https://img.shields.io/badge/PRs-welcome-yellow.svg)](https://makeapullrequest.com/) [![GitHub commits](https://img.shields.io/github/commits-since/epiforecasts/covidregionaldata/0.9.2.svg?color=orange)](https://GitHub.com/epiforecasts/covidregionaldata/commit/master/) 

[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03290/status.svg)](https://doi.org/10.21105/joss.03290) [![Zenodo](https://zenodo.org/badge/271601189.svg)](https://zenodo.org/badge/latestdoi/271601189)

Interface to subnational and national level COVID-19 data sourced from both official sources, such as Public Health England in the UK, and from other COVID-19 data collections, including the World Health Organisation (WHO), European Centre for Disease Prevention and Control (ECDC), John Hopkins University (JHU), Google Open Data and others. This package is designed to streamline COVID-19 data extraction, cleaning, and processing from a range of data sources in an open and transparent way. This allows users to inspect and scrutinise the data, and tools used to process it, at every step. For all countries supported, data includes a daily time-series of cases and, wherever available, data on deaths, hospitalisations, and tests. National level data is also supported using a range of data sources as well as line list data and links to intervention data sets.

## Installation

Install from CRAN:

```{r, eval = FALSE}
install.packages("covidregionaldata")
```

Install the stable development version of the package with:

```{r, eval = FALSE}
install.packages("covidregionaldata",
  repos = "https://epiforecasts.r-universe.dev"
)
```

Install the unstable development version of the package with:

```{r, eval = FALSE}
remotes::install_github("epiforecasts/covidregionaldata")
```

## Quick start

[![Documentation](https://img.shields.io/badge/Documentation-lightgrey.svg?style=flat)](https://epiforecasts.io/covidregionaldata/)


Load `covidregionaldata`, `dplyr`, `scales`, and `ggplot2` (all used in this quick start),

```{r, message = FALSE}
library(covidregionaldata)
library(dplyr)
library(ggplot2)
library(scales)
```

### Setup data caching

This package can optionally use a data cache from `memoise` to locally cache downloads. This can be enabled using the following (this will use the temporary directory by default),

```{r}
start_using_memoise()
```

To stop using `memoise` use,

```{r, eval = FALSE}
stop_using_memoise()
```

and to reset the cache (required to download new data),

```{r, eval = FALSE}
reset_cache()
```

### National data

To get worldwide time-series data by country (sourced from the World Health Organisation (WHO) by default but also optionally from the European Centre for Disease Control (ECDC), John Hopkins University, or the Google COVID-19 open data project), use:

```{r}
nots <- get_national_data()
nots
```

This can also be filtered for a country of interest,

```{r}
g7 <- c(
  "United States", "United Kingdom", "France", "Germany",
  "Italy", "Canada", "Japan"
)
g7_nots <- get_national_data(countries = g7, verbose = FALSE)
```

Using this data we can compare case information between countries, for example here is the number of deaths over time for each country in the G7:

```{r g7_plot, warning = FALSE, message = FALSE}
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

### Subnational data

To get time-series data for subnational regions of a specific country, for example by level 1 region in the UK, use:

```{r}
uk_nots <- get_regional_data(country = "UK", verbose = FALSE)
uk_nots
```

Now we have the data we can create plots, for example the time-series of the number of cases for each region:

```{r uk_plot, warning = FALSE, message = FALSE}
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

See `get_available_datasets()` for supported regions and subregional levels.
To view what datasets we currently have subnationaldata for, along with their current status, check the
[supported countries](https://epiforecasts.io/covidregionaldata/articles/supported-countries.html) page
or build the [supported countries vignette](vignettes/supported-countries.Rmd).

For further examples see the [quick start vignette](https://github.com/epiforecasts/covidregionaldata/blob/master/vignettes/quickstart.Rmd). Additional subnational data are supported via the `JHU()` and `Google()` classes. Use the `available_regions()` method once these data have been downloaded and cleaned (see their examples) for subnational data they internally support.

## Citation

If using `covidregionaldata` in your work please consider citing it using the following,

```{r, echo = FALSE}
citation("covidregionaldata")
```

## Development

[![Development](https://img.shields.io/badge/Wiki-lightblue.svg?style=flat)](https://github.com/epiforecasts/covidregionaldata/wiki/)

This package is the result of work from a number of contributors (see contributors list [here](https://epiforecasts.io/covidregionaldata/authors.html)). We would like to thank the [CMMID COVID-19 working group
](https://cmmid.github.io/groups/ncov-group.html) for insightful comments and feedback.

We welcome contributions and new contributors! We particularly appreciate help adding new data sources for countries at sub-national level, or work on priority problems in the [issues](https://github.com/epiforecasts/covidregionaldata/issues). Please check and add to the issues, and/or add a [pull request](https://github.com/epiforecasts/covidregionaldata/pulls). For more details, start with the [contributing guide](https://github.com/epiforecasts/covidregionaldata/wiki/Contributing). For details of the steps required to add support for a dataset see the [adding data guide](https://github.com/epiforecasts/covidregionaldata/wiki/Adding-Data).
---
title: "Package overview"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Package overview}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

## Installing and loading the package

The package can either be installed from CRAN, from our [`r-universe`](https://epiforecasts.r-universe.dev/ui) repository, or from GitHub. See the README for details. Once installed load the package using the following,

```{r, eval=TRUE}
library(covidregionaldata)
```

## Worldwide data

### Accessing national data

Both the World Health Organisation (WHO) and European Centre for Disease Control (ECDC) provide worldwide national data. Access national level data for any country using:

```{r, eval=FALSE}
get_national_data()
```

This returns daily new and cumulative (total) cases, and where available, deaths, hospitalisations, and tests. For a complete list of variables returned, see section 5, "Data glossary" below. See the documentation (`?get_national_data`) for details of optional arguments.

Data is returned with no gaps in the structure of the data by country over time, and NAs fill in where data are not available.

## Sub-national time-series data

### Accessing sub-national data

Access sub-national level data for a specific country over time using `get_regional_data()`. Use `get_available_datasets()` to explore the currently supported sub-national datasets and select the data set of interest using the `country` (selects the country of interest), and `level` (selects the spatial scale of the data) arguments of `get_regional_data`.

This function returns daily new and cumulative (total) cases, and where available, deaths, hospitalisations, and tests. For a complete list of variables returned, see section 5, "Data glossary" below. See the documentation (`?get_regional_data`) for details of optional arguments.

As for national level data any gaps in reported data are filled with NAs.

For example, data for France Level 1 regions over time can be accessed using:

```{r}
get_regional_data(country = "france")
```

This data then has the following format: 

```{r, echo=FALSE, eval=TRUE, message=FALSE}
start_using_memoise()
knitr::kable(
  tail(get_regional_data(country = "france"), n = 5)
)
```

Alternatively, the same data can be accessed using the underlying class as follows (the France object now contains data at each processing step and the methods used at each step), 

```{r, eval=FALSE, message=FALSE}
france <- France$new(get = TRUE)
france$return()
```
### Level 1 and Level 2 regions

All countries included in the package (see below,"Coverage") have data for regions at the admin-1 level, the largest administrative unit of the country (e.g. state in the USA). Some countries also have data for smaller areas at the admin-2 level (e.g. county in the USA).

Data for Level 2 units can be returned by using the `level = "2"` argument. The dataset will still show the corresponding Level 1 region.

An example of a country with Level 2 units is France, where Level 2 units are French departments:

```{r}
get_regional_data(country = "france", level = "2")
```

This data again has the following format: 

```{r, echo=FALSE, eval=TRUE, message=FALSE}
knitr::kable(
  tail(get_regional_data(country = "france", level = "2"), n = 5)
)
```

### Totals

For totalled data up to the most recent date available, use the `totals` argument.

```{r}
get_regional_data("france", totals = TRUE)
```

This data now has no date variable and reflects the latest total:

```{r, echo=FALSE, eval=TRUE, message=FALSE}
knitr::kable(
  tail(get_regional_data(country = "france", totals = TRUE), n = 5)
)
```



## Data glossary

#### Subnational data

The data columns that will be returned by `get_regional_data()` are listed below.

To standardise across countries and regions, the columns returned for each country will _always_ be the same. If the corresponding data was missing from the original source then that data field is filled with NA values (or 0 if accessing totals data).

Note that Date is not included if the `totals` argument is set to TRUE. Level 2 region/level 2 region code are not included if the `level = "1"`.

* `date`: the date that the counts were reported (YYYY-MM-DD).

* `level_1_region`: the level 1 region name. This column will be named differently for different countries (e.g. state, province).

* `level_1_region_code`: a standard code for the level 1 region. The column name reflects the specific administrative code used. Typically data returns the iso_3166_2 standard, although where not available the column will be named differently to reflect its source.

* `level_2_region`: the level 2 region name. This column will be named differently for different countries (e.g. city, county).

* `level_2_region_code`: a standard code for the level 2 region. The column will be named differently for different countries (e.g. `fips` in the USA).

* `cases_new`: new reported cases for that day.

* `cases_total`: total reported cases up to and including that day.

* `deaths_new`: new reported deaths for that day.

* `deaths_total`: total reported deaths up to and including that day.

* `recovered_new`: new reported recoveries for that day.

* `recovered_total`: total reported recoveries up to and including that day.

* `hosp_new`: new reported hospitalisations for that day.

* `hosp_total`: total reported hospitalisations up to and including that day (note this is cumulative total of new reported, _not_ total currently in hospital).

* `tested_new`: tests for that day.

* `tested_total`: total tests completed up to and including that day.

#### National data

In addition to the above, the following columns are included when using `get_national_data()`.

* `un_region`: country geographical region defined by the United Nations.

* `who_region`: only included when `source = "WHO"`. Country geographical region defined by WHO.

* `population_2019`: only included when `source = "ECDC"`. Total country population estimate in 2019.
---
title: "Supported countries and their support status"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Supported countries and their support status}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, include = FALSE}
library(covidregionaldata)
library(dplyr)
library(knitr)
library(dplyr)
library(ggplot2)
library(sf)
```

## Countries with subnational data

This map shows what countries have level 1 and level 2 subregion data directly from an official source within that country. Please note other countries may be provided through our interface to external data sources, such as `WHO()` and `JHU()`.

```{r, supported_region_plot, echo=FALSE, message=FALSE}

regional_countries <- get_available_datasets() %>%
  filter(type == "regional")

regional_countries_l2 <- regional_countries %>%
  filter(!(is.na(level_2_region)))

world <- map_data("world")

supported_countries <- mutate(
  world,
  fill = case_when(
    region %in% regional_countries_l2[["class"]] ~ "Level 2",
    region %in% regional_countries[["class"]] ~ "Level 1",
    TRUE ~ "Unsupported"
  )
)

ggplot(supported_countries, aes(long, lat, fill = fill, group = group)) +
  geom_polygon(color = "black", lwd = 0.1) +
  scale_fill_manual(
    name = "",
    values = c("#0072b2", "#cc79a7", "grey80")
  ) +
  theme_void() +
  theme(legend.position = "bottom") +
  coord_sf(ylim = c(-55, 80))
```

## Status

Dataset status is shown in the table below. Please see our [hosted page](https://epiforecasts.io/covidregionaldata/articles/supported-countries.html) for up to date information for the CRAN status of data sets. Please note that due to our release schedule datasets may remain non-functional if broken using the CRAN version for some time even if fixed on GitHub. Also note that transient issues may affect our testing of datasets and so our checks may occasionally show a spurious failure. 

```{r, echo = FALSE}
datasets <- get_available_datasets() %>%
  arrange(origin) %>%
  select(Origin = origin, Method = class) %>%
  mutate(
    `GitHub status` = paste0(
      "[![", Method,
      "](https://github.com/epiforecasts/covidregionaldata/workflows/",
      Method, "/badge.svg)]",
      "(https://github.com/epiforecasts/covidregionaldata/actions/workflows/", # nolint
      Method, ".yaml)"
    ),
    `CRAN status` = "*working*"
  )
kable(datasets)
```
---
title: "Testing"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Testing}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE
)
```


## Summary


This PR moves testing into the class methods, allowing tests to be ran interactively from within the country class object, and from within the existing testing infrastructure.

`DataClass` is the parent class of all country classes, such as `Italy`, `UK`, `USA`, etc. The generic testing function is the method `test` within `DataClass`.

Interactively, tests can be ran by doing the following:


```{r}
library(covidregionaldata)
ukdata <- UK$new(level = "1", verbose = FALSE)
# could use anything but I used the acutal test one here for simplicity
ukdata$test(snapshot_dir = "../tests/testthat/custom_data/")
```


Here, I break down the different components of this function to walk through how the tests are run.

```{r, echo = TRUE}
ukdata$test
```

For a given country, such as the UK, you would make this by calling something like `ukdata <- UK$new(level = "1")` and can run the steps one by one by calling the respective methods: `ukdata$download(); ukdata$clean(); ukdata$process()`. To run the tests you would call `ukdata$test(shapshot_dir = "place/2/save/snapshots")` By default, if `snapshot_dir` is not given it will use a generic `snapsots` directory. The `snapshot_dir` argument specifies a directory to save data snapshots to and if it doesn't exist it will make it for you. The file is constructed internally from the object you are testing and the data level, e.g. `Uk_level_1.rds`. Rather than running tests on your active class, the code first makes a clone of your class which is then used for tests: `self_copy <- self$clone()` The snapshot path is the path to a rds file where you either have stored some raw downloaded data to test, or where you want a downloaded snapshot of the data to end up. This is handled by the function `test_download()`. The cloned class is passed to this function, along with the download parameter, which dictates whether to overwrite the snapshot file provided.

```{r, echo = TRUE}
test_download
```

As shown by the code, if the data is to be downloaded (either through requesting this with the download parameter or by providing a path to a none existent file) then the download method is called on the class copy (`DataClass_obj`): `DataClass_obj$download()`. After this has downloaded the code tests the data is a data.frame, is not empty and has at least 2 columns. The code then takes the first 250 rows to work on, as using the whole data set could be very slow and for the purpose of testing is not needed. This sliced data is then saved to the snapshot file provided for use later on. If the data is not to be downloaded, the snapshot of data saved in the snapshot path is loaded as the raw data `DataClass_obj$data$raw`.

Once the data is downloaded, the cleaning methods are then tested `test_cleaning()`. Again this function takes the copied class as the argument and runs through tests on the clean data. These tests check the data is a data.frame, is not empty, has more than 2 columns and that the method `avaliable_regions` returns a list of characters. In addition, `expect_clean_cols` checks the date column is an s3 class and that region level column is a character in the cleaned data (`data$clean`).

```{r, echo = TRUE}
test_cleaning
```

Once cleaning has been tested `test_processing()` is called to process the data and run tests to check it all works. Again this function takes the clone of the class to work on, the same clone which has been called with the preceding functions. These tests check the data is a data.frame, is not empty, has more than 2 columns. In addition `expect_processed_cols` checks that processed data columns date, cases_new, cases_total, deaths_new, deaths_total and that region level have the correct types.

In processing there is an extra parameter called `all`. This parameter, if `TRUE` runs the processing step with both localised as `TRUE` and `FALSE`, making another copy of the object to check the localised data on so not to influence further tests. If `FALSE` processing is ran with localised set to whatever it is set at prior to `test()` being invoked.

```{r, echo = TRUE}
test_processing
```

After processing the return method is tested with `test_return` which check the data is a data.frame, is not empty, has more than 2 columns and that all columns contain data and are not just composed of NAs.

```{r, echo = TRUE}
test_return
```

These tests form the generic tests applied to all classes. However, country specific tests are then called by calling the method `specific_tests` if that country has specific tests defined. So for `Italy`, where there are no specific tests, no specific tests are called, but for `UK`, `WHO` and `ECDC` specific tests are ran though, which are defined in their own country class. These functions should take a clone of the class as an argument (`self_copy` ) and any additional arguments they may need, such as a path to NHS included data for `UK`.


## Integration with *testthat*

As well as interactive tests the `test()` method is also used by `testthat` when conducting package level tests but with the argument `all = TRUE`. This is done in the file `tests/testthat/custom_tests/test_regional_dataset.R`
---
title: "Quick start"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Quick start}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE
)
```

# Quick start to `covidregionaldata`

`covidregionaldata` is designed to extract national and subnational Covid-19 data from publicly available sources, both daily and aggregated over time.

For example, we could:

* Get the Covid-19 cases and deaths for all countries for each day since reporting began.
* Get the total Covid-19 cases and deaths for all countries up to today.
* Get the number of new cases for each state in the United States over time.
* Get the total number of deaths for all local authorities in the United Kingdom.

We will demonstrate some of these examples below, as well as some useful data manipulation and visualisation.

## Retrieving national cases and deaths over time

Let's say that we want to see the evolution of the pandemic over time in all countries. To do this, we need to get the number of cases (positive test results) and deaths for each country and date since the pandemic started. Firstly load `covidregionaldata`, then call `get_national_data()`:

```{r setup}
library(covidregionaldata)
all_countries <- get_national_data()
```

We can take a look at this data by printing it to the console:

```{r}
all_countries
```

We can see that this provides information about the number of cases and deaths for each country over time. As a tibble is returned, we can easily manipulate and plot this using `tidyverse` packages. For example, we can plot the total deaths in Italy over time.

```{r fig.width = 6, message = FALSE}
library(dplyr)
library(ggplot2)

all_countries %>%
  filter(country == "Italy") %>%
  ggplot() +
  aes(x = date, y = deaths_total) +
  geom_line() +
  labs(x = NULL, y = "All reported Covid-19 deaths, Italy") +
  theme_minimal()
```

We could have also filtered by country using `get_national_data(countries = "Italy")` to achieve the same result.

To plot the evolution of cases for several countries since the start of the pandemic, we could perform something like the following:

```{r fig.width = 6}
all_countries %>%
  filter(country %in% c(
    "Italy", "United Kingdom", "Spain",
    "United States"
  )) %>%
  ggplot() +
  aes(x = date, y = cases_total, colour = country) +
  geom_line() +
  labs(x = NULL, y = "All reported Covid-19 cases", colour = "Country") +
  theme_minimal() +
  theme(legend.position = "bottom")
```

## Retrieving total cases and deaths for all countries

Using `get_national_data(totals = TRUE)` we can obtain total cases and deaths for all nations up to the latest date reported. This is useful to get a snapshot of the current running total for each country. Note how the data is sorted by the total number of cases:

```{r}
all_countries_totals <- get_national_data(totals = TRUE, verbose = FALSE)
all_countries_totals
```

## Retrieving subnational region data

As well as national data, we can also retrieve subnational data, such as states and counties in the US, or regions and local authorities in the UK, and perform similar manipulations. See the [README](https://epiforecasts.io/covidregionaldata/) for a current list of countries with subnational data.

To do this, we use `get_regional_data()` (we could also use the underlying method as follows `us <- USA$new(get = TRUE); us$return()`). We'll get the state-level data for the United States over time:

```{r}
usa_states <- get_regional_data(country = "USA")
```

This retrieves the new and total cases and deaths for all states over time in the United States. Note that the `country` argument to the `get_regional_data()` function must be specified. Let's plot out how the number of cases has evolved for a selection of states:

```{r fig.width = 6, warning = FALSE}
usa_states %>%
  filter(state %in% c("New York", "Texas", "Florida")) %>%
  ggplot() +
  aes(x = date, y = cases_total, colour = state) +
  geom_line() +
  labs(x = NULL, y = "All reported Covid-19 cases", colour = "U.S. state") +
  theme_minimal() +
  theme(legend.position = "bottom")
```
The `totals` argument can also be used in `get_regional_data` to return the cumulative number of cases or deaths over time.

```{r,}
usa_states_totals <- get_regional_data(
  country = "USA", totals = TRUE,
  verbose = FALSE
)
```

Let's plot the total deaths for each state, ordered by total deaths:

```{r fig.width = 6}
usa_states_totals %>%
  ggplot() +
  aes(x = reorder(state, -deaths_total), y = deaths_total) +
  geom_bar(stat = "identity") +
  labs(x = "U.S. states", y = "All reported Covid-19 deaths") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.title.y = element_text(hjust = 1),
    legend.position = "bottom"
  )
```

Most countries have a hierarchy of political subdivisions. Where data are available, use the `level = "2"` argument to `get_regional_data`. We can dig down further into each US state by looking at data by county:

```{r, eval = FALSE}
usa_counties <- get_regional_data(country = "USA", level = "2", verbose = FALSE)
```

This returns data for each county on each day of reported data. As there are more than 3,000 counties in the US, with some reporting data since the 21st January, this makes for a very large download.

The same method will work for different countries where data are available, returning relevant sub-national units for each country. For example, in Germany data are available for Bundesland (federated state) and Kreise (district) areas.

Additional subnational data are supported via the `JHU()` and `Google()` classes. Use the `available_regions()` method once these data have been downloaded and cleaned (see their examples) for subnational data they internally support.

## Mapping data

We provide the relevant national and subnational ISO codes, or local alternatives, to enable mapping.  We can map Covid-19 data at national or subnational level, using any standard shapefile. In R, the `rnaturalearth` package already provides shapefiles for all country borders, and the largest (administrative level 1) subnational units for most countries.

For simplicity here we will use the `rworldmap` package, which uses data from `rnaturalearth` for national country borders and inbuilt mapping. We can join Covid-19 data to the map using the `iso_code`.

```{r fig.width = 6}
# Get latest worldwide WHO data
map_data <- get_national_data(totals = TRUE, verbose = FALSE) %>%
  rworldmap::joinCountryData2Map(
    joinCode = "ISO2",
    nameJoinColumn = "iso_code"
  )
# Produce map
rworldmap::mapCountryData(map_data,
  nameColumnToPlot = "deaths_total",
  catMethod = "fixedWidth",
  mapTitle = "Total Covid-19 deaths to date",
  addLegend = TRUE
)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{calculate_columns_from_existing_data}
\alias{calculate_columns_from_existing_data}
\title{Cumulative counts from daily counts or daily counts from cumulative,
dependent on which columns already exist}
\usage{
calculate_columns_from_existing_data(data)
}
\arguments{
\item{data}{A data frame}
}
\value{
A data frame with extra columns if required
}
\description{
Checks which columns are missing (cumulative/daily counts) and
if one is present and the other not then calculates the second from the
first.
}
\seealso{
Compulsory processing functions
\code{\link{add_extra_na_cols}()},
\code{\link{complete_cumulative_columns}()},
\code{\link{fill_empty_dates_with_na}()}
}
\concept{compulsory_processing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Covid19DataHub.R
\name{Covid19DataHub}
\alias{Covid19DataHub}
\title{R6 Class containing specific attributes and methods for Covid19 Data Hub}
\source{
\url{https://covid19datahub.io/articles/data.html}
}
\description{
Attributes and methods for COVID-19 data provided by the
Covid19 Data Hub
}
\details{
This dataset supports both national and subnational data sources
with national level data returned by default. National data is sourced from
John Hopkins University and so we recommend using the JHU class included in
this package. Subnational data is supported for a subset of countries which
can be found after cleaning using the \code{available_regions()} method,
see the examples for more details. These data sets are minimally cleaned
data files hosted by the team at COVID19 Data Hub so please see their
source repository for further details
(https://github.com/covid19datahub/COVID19/#data-sources)
If using for analysis checking the source for further details is
strongly advised.

If using this class please cite:
"Guidotti et al., (2020). COVID-19 Data Hub
Journal of Open Source Software, 5(51),
2376, https://doi.org/10.21105/joss.02376"
}
\examples{
# nolint start
\dontrun{
# set up a data cache
start_using_memoise()

# get all countries data
cv19dh <- Covid19DataHub$new(level = "1", get = TRUE)
cv19dh$return()

# show available regions with data at the second level of interest
cv19dh_level_2 <- Covid19DataHub$new(level = "2")
cv19dh_level_2$download()
cv19dh_level_2$clean()
cv19dh$available_regions()

# get all region data for the uk
cv19dh_level_2$filter("uk")
cv19dh_level_2$process()
cv19dh_level_2$return()

# get all regional data for the UK
uk <- Covid19DataHub$new(regions = "uk", level = "2", get = TRUE)
uk$return()

# get all subregional data for the UK
uk <- Covid19DataHub$new(regions = "uk", level = "3", get = TRUE)
uk$return()
}
# nolint end
}
\seealso{
Aggregated data sources
\code{\link{Google}},
\code{\link{JHU}}

National data sources
\code{\link{ECDC}},
\code{\link{Google}},
\code{\link{JHU}},
\code{\link{JRC}},
\code{\link{WHO}}

Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{aggregations}
\concept{dataset}
\concept{national}
\concept{subnational}
\section{Super classes}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{\link[covidregionaldata:CountryDataClass]{covidregionaldata::CountryDataClass}} -> \code{Covid19DataHub}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of country to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{level_data_urls}}{List of named links to raw data. The first, and
only entry, is be named main.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-clean_common}{\code{Covid19DataHub$clean_common()}}
\item \href{#method-clone}{\code{Covid19DataHub$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="set_region_codes">}\href{../../covidregionaldata/html/DataClass.html#method-set_region_codes}{\code{covidregionaldata::DataClass$set_region_codes()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="CountryDataClass" data-id="filter">}\href{../../covidregionaldata/html/CountryDataClass.html#method-filter}{\code{covidregionaldata::CountryDataClass$filter()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Covid19 Data Hub specific data cleaning.
This takes all the raw data, renames some columns and checks types.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Covid19DataHub$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Covid19DataHub$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/SouthAfrica.R
\name{SouthAfrica}
\alias{SouthAfrica}
\title{SouthAfrica Class for downloading, cleaning and processing notification data}
\source{
\url{https://github.com/dsfsi/covid19za/}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for South Africa.
}
\examples{
\dontrun{
region <- SouthAfrica$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{SouthAfrica}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{SouthAfrica$set_region_codes()}}
\item \href{#method-clean_common}{\code{SouthAfrica$clean_common()}}
\item \href{#method-clone}{\code{SouthAfrica$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SouthAfrica$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Province level data cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SouthAfrica$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SouthAfrica$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/shared-methods.R
\name{DataClass}
\alias{DataClass}
\title{R6 Class containing non-dataset specific methods}
\description{
A parent class containing non-dataset specific methods.
}
\details{
All data sets have shared methods for extracting geographic codes,
downloading, processing, and returning data. These functions are contained
within this parent class and so are accessible by all data sets which
inherit from here. Individual data sets can overwrite any functions or
fields providing they define a method with the same name, and can be
extended with additional functionality. See the individual method
documentaion for further details.
}
\seealso{
Data interface functions
\code{\link{CountryDataClass}},
\code{\link{get_available_datasets}()},
\code{\link{get_national_data}()},
\code{\link{get_regional_data}()},
\code{\link{initialise_dataclass}()}
}
\concept{interface}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{the origin of the data source. For regional data sources
this will usually be the name of the country.}

\item{\code{data}}{Once initialised, a list of named data frames: raw
(list of named raw data frames) clean (cleaned data) and processed
(processed data). Data is accessed using \verb{$data}.}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{region_name}}{string Name for the region column, e.g. 'region'.
This field is filled at initialisation with the region name for the
specified level (supported_region_names$level).}

\item{\code{code_name}}{string Name for the codes column, e.g. 'iso_3166_2'
Filled at initialisation with the code name associated with the
requested level (supported_region_codes$level).}

\item{\code{codes_lookup}}{string or tibble Region codes for the target origin
filled by origin specific codes in
\href{#method-set_region_codes}{\code{set_region_codes()}}}

\item{\code{data_urls}}{List of named common and shared url links to raw data.
Prefers shared if there is a name conflict.}

\item{\code{common_data_urls}}{List of named links to raw data that are common
across levels. The first entry should be named main.}

\item{\code{level_data_urls}}{List of named links to raw data that are level
specific. Any urls that share a name with a url from
\code{common_data_urls} will be selected preferentially. Each top level
list should be named after a supported level.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{level}}{target region level. This field is filled at initialisation
using user inputs or defaults in \verb{$new()}}

\item{\code{data_name}}{string. The country name followed by the level. E.g.
"Italy at level 1"}

\item{\code{totals}}{Boolean. If TRUE, returns totalled data per region
up to today's date. This field is filled at initialisation using user
inputs or defaults in \verb{$new()}}

\item{\code{localise}}{Boolean. Should region names be localised.
This field is filled at initialisation using user inputs or defaults
in \verb{$new()}}

\item{\code{verbose}}{Boolean. Display information at various stages.
This field is filled at initialisation. using user inputs or defaults
in \verb{$new()}}

\item{\code{steps}}{Boolean. Keep data from each processing step.
This field is filled at initialisation.using user inputs or defaults
in \verb{$new()}}

\item{\code{target_regions}}{A character vector of regions to filter for. Used
by the \verb{filter method}.}

\item{\code{process_fns}}{array, additional, user supplied functions to process
the data.}

\item{\code{filter_level}}{Character The level of the data to filter at.
Defaults to the target level.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{DataClass$set_region_codes()}}
\item \href{#method-new}{\code{DataClass$new()}}
\item \href{#method-download}{\code{DataClass$download()}}
\item \href{#method-download_JSON}{\code{DataClass$download_JSON()}}
\item \href{#method-clean}{\code{DataClass$clean()}}
\item \href{#method-clean_common}{\code{DataClass$clean_common()}}
\item \href{#method-available_regions}{\code{DataClass$available_regions()}}
\item \href{#method-filter}{\code{DataClass$filter()}}
\item \href{#method-process}{\code{DataClass$process()}}
\item \href{#method-get}{\code{DataClass$get()}}
\item \href{#method-return}{\code{DataClass$return()}}
\item \href{#method-summary}{\code{DataClass$summary()}}
\item \href{#method-test}{\code{DataClass$test()}}
\item \href{#method-clone}{\code{DataClass$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Place holder for custom country specific function to load
region codes.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize function used by all \code{DataClass} objects.
Set up the \code{DataClass} class with attributes set to input parameters.
Should only be called by a \code{DataClass} class object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$new(
  level = "1",
  filter_level,
  regions,
  totals = FALSE,
  localise = TRUE,
  verbose = TRUE,
  steps = FALSE,
  get = FALSE,
  process_fns
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{level}}{A character string indicating the target administrative
level of the data with the default being "1". Currently supported
options are level 1 ("1) and level 2 ("2").}

\item{\code{filter_level}}{A character string indicating the level to filter at.
Defaults to the level of the data if not specified and if not otherwise
defined in the class.
Use \code{get_available_datasets()} for supported options by dataset.}

\item{\code{regions}}{A character vector of target regions to be assigned to
the\code{target_regions} field if present.}

\item{\code{totals}}{Logical, defaults to FALSE. If TRUE, returns totalled
data per region up to today's date. If FALSE, returns the full dataset
stratified by date and region.}

\item{\code{localise}}{Logical, defaults to TRUE. Should region names be
localised.}

\item{\code{verbose}}{Logical, defaults to TRUE. Should verbose processing}

\item{\code{steps}}{Logical, defaults to FALSE. Should all processing and
cleaning steps be kept and output in a list.}

\item{\code{get}}{Logical, defaults to FALSE. Should the class \code{get} method be
called (this will download, clean, and process data at initialisation).}

\item{\code{process_fns}}{Array, additional functions to process the data.
Users can supply their own functions here which would act on clean data
and they will be called alongside our default processing functions.
The default optional function added is \code{set_negative_values_to_zero}.
if process_fns is not set (see \code{process_fns} field for all defaults).
If you want to keep this when supplying your own processing functions
remember to add it to your list also. If you feel you have created a
cool processing function that others could benefit from please submit a
Pull Request to our \href{https://github.com/epiforecasts/covidregionaldata}{github repository}
and we will consider adding it to the package.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-download"></a>}}
\if{latex}{\out{\hypertarget{method-download}{}}}
\subsection{Method \code{download()}}{
Download raw data from \code{data_urls}, stores a named list
of the \code{data_url} name and the corresponding raw data table in
\code{data$raw}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$download()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-download_JSON"></a>}}
\if{latex}{\out{\hypertarget{method-download_JSON}{}}}
\subsection{Method \code{download_JSON()}}{
Download raw data from \code{data_urls}, stores a named list
of the \code{data_url} name and the corresponding raw data table in
\code{data$raw}. Designed as a drop-in replacement for \code{download} so
it can be used in sub-classes.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$download_JSON()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean"></a>}}
\if{latex}{\out{\hypertarget{method-clean}{}}}
\subsection{Method \code{clean()}}{
Cleans raw data (corrects format, converts column types,
etc). Works on raw data and so should be called after
\href{#method-download}{\code{download()}}
Calls the specific class specific cleaning method (\code{clean_common})
followed by level specific cleaning methods.
\code{clean_level_[1/2]}. Cleaned data is stored in \code{data$clean}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$clean()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Cleaning methods that are common across a class.
By default this method is empty as if any code is required it should be
defined in a child class specific \code{clean_common} method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-available_regions"></a>}}
\if{latex}{\out{\hypertarget{method-available_regions}{}}}
\subsection{Method \code{available_regions()}}{
Show regions that are available to be used for
filtering operations. Can only be called once \code{clean()} has been
called. Filtering level is determined by checking the \code{filter_level}
field.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$available_regions(level)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{level}}{A character string indicating the level to filter at.
Defaults to using the \code{filter_level} field if not specified}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-filter"></a>}}
\if{latex}{\out{\hypertarget{method-filter}{}}}
\subsection{Method \code{filter()}}{
Filter cleaned data for a specific region  To be called
after \href{#method-clean}{\code{clean()}}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$filter(regions, level)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{regions}}{A character vector of target regions. Overrides the
current class setting for \code{target_regions}.}

\item{\code{level}}{Character The level of the data to filter at. Defaults to
the lowest level in the data.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-process"></a>}}
\if{latex}{\out{\hypertarget{method-process}{}}}
\subsection{Method \code{process()}}{
Processes data by adding and calculating absent columns.
Called on clean data (after \href{#method-clean}{\code{clean()}}).
Some countries may have data as new events (e.g. number of
new cases for that day) whilst others have a running total up to that
date. Processing calculates these based on what the data comes with
via the functions \code{region_dispatch()} and \code{process_internal()},
which does the following:
\itemize{
\item{Adds columns not present in the data \code{add_extra_na_cols()}}
\item{Ensures there are no negative values
\code{set_negative_values_to_zero()}}
\item{Removes NA dates \code{fill_empty_dates_with_na()}}
\item{Calculates cumulative data \code{complete_cumulative_columns()}}
\item{Calculates missing columns from existing ones
\code{calculate_columns_from_existing_data()}}
}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$process(process_fns)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{process_fns}}{Array, additional functions to process the data.
Users can supply their own functions here which would act on clean data
and they will be called alongside our default processing functions.
The default optional function added is \code{set_negative_values_to_zero}.
if process_fns is not set (see \code{process_fns} field for all defaults).}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get"></a>}}
\if{latex}{\out{\hypertarget{method-get}{}}}
\subsection{Method \code{get()}}{
Get data related to the data class. This runs each distinct
step in the workflow in order.
Internally calls \href{#method-download}{\code{download()}},
\href{#method-clean}{\code{clean()}},
\href{#method-filter}{\code{filter()}} and
\href{#method-process}{\code{process()}}
\code{download}, \code{clean}, \code{filter} and \code{process} methods.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$get()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-return"></a>}}
\if{latex}{\out{\hypertarget{method-return}{}}}
\subsection{Method \code{return()}}{
Return data. Designed to be called after
\href{#method-process}{\code{process()}}
this uses the steps argument to return either a
list of all the data preserved at each step or just the processed data.
For most datasets a custom method should not be needed.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$return()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-summary"></a>}}
\if{latex}{\out{\hypertarget{method-summary}{}}}
\subsection{Method \code{summary()}}{
Create a table of summary information for the data set
being processed.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$summary()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Returns a single row summary tibble containing the origin of the
data source, class, level 1 and 2 region names, the type of data,
the urls of the raw data and the columns present in the raw data.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-test"></a>}}
\if{latex}{\out{\hypertarget{method-test}{}}}
\subsection{Method \code{test()}}{
Run tests on a country class instance. Calling \code{test()} on a
class instance runs tests with the settings in use. For example, if you
set \code{level = "1"} and \code{localise = FALSE} the tests will be run on level 1
data which is not localised. Rather than downloading data for a test
users can provide a path to a snapshot file of data to test instead.
Tests are run on a clone of the class. This method calls generic tests
for all country class objects. It also calls country specific tests
which can be defined in an individual country class method called
\code{specific_tests()}. The snapshots contain the first 1000 rows of data.
For more details see the
\href{https://github.com/epiforecasts/covidregionaldata/blob/master/vignettes/testing.Rmd}{'testing' vignette}: \code{vignette(testing)}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataClass$test(
  download = FALSE,
  snapshot_dir = paste0(tempdir(), "/snapshots"),
  all = FALSE,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{download}}{logical. To download the data (TRUE) or use a snapshot
(FALSE). Defaults to FALSE.}

\item{\code{snapshot_dir}}{character_array the name of a directory to save the
downloaded data or read from. If not defined a directory called
'snapshots' will be created in the temp directory. Snapshots are saved as
rds files with the class name and level: e.g. \code{Italy_level_1.rds}.}

\item{\code{all}}{logical. Run tests with all settings (TRUE) or with those
defined in the current class instance (FALSE). Defaults to FALSE.}

\item{\code{...}}{Additional parameters to pass to \code{specific_tests}}
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
\if{html}{\out{<div class="r">}}\preformatted{DataClass$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/utils.R
\name{download_excel}
\alias{download_excel}
\title{Download Excel Documents}
\usage{
download_excel(url, archive, verbose = FALSE, transpose = TRUE, ...)
}
\arguments{
\item{url}{Character string containing the full URL to the Excel document.}

\item{archive}{Character string naming the file name to assign in the
temporary directory.}

\item{verbose}{Logical, defaults to \code{TRUE}. Should verbose processing
messages and warnings be returned.}

\item{transpose}{Logical, should the read in data be transposed}

\item{...}{Additional parameters to pass to \code{read_excel()}.}
}
\value{
A \code{data.frame}.
}
\description{
Download Excel Documents
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Germany.R
\name{Germany}
\alias{Germany}
\title{Germany Class for downloading, cleaning and processing notification data}
\source{
\url{https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region level 1 and 2 data for Germany.
}
\examples{
\dontrun{
region <- Germany$new(verbose = TRUE, steps = TRUE, level = "2", get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Germany}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data. The first, and
only entry, is be named main.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Germany$set_region_codes()}}
\item \href{#method-clean_common}{\code{Germany$clean_common()}}
\item \href{#method-clean_level_1}{\code{Germany$clean_level_1()}}
\item \href{#method-clean_level_2}{\code{Germany$clean_level_2()}}
\item \href{#method-clone}{\code{Germany$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Germany$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Common Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Germany$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
Bundesland Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Germany$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_2"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_2}{}}}
\subsection{Method \code{clean_level_2()}}{
Landkreis Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Germany$clean_level_2()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Germany$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/test-DataClass.R
\name{test_return}
\alias{test_return}
\title{Test return method works correctly}
\usage{
test_return(DataClass_obj)
}
\arguments{
\item{DataClass_obj}{The R6Class object to perform checks on.
Must be a \code{DataClass} or \code{DataClass} child object.}
}
\description{
Test data can be returned correctly using the return method.
return is invoked to generate returned data which is then checked to ensure
it is a data.frame, not empty and has at least 2 columns. Each column is then
checked to ensure it contains data and is not just composed of NAs.
}
\seealso{
Functions used for testing data is cleaned and processed correctly
\code{\link{expect_clean_cols}()},
\code{\link{expect_columns_contain_data}()},
\code{\link{expect_processed_cols}()},
\code{\link{test_cleaning}()},
\code{\link{test_download}()},
\code{\link{test_processing}()}
}
\concept{tests}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{region_dispatch}
\alias{region_dispatch}
\title{Control Grouping Variables used in process_internal}
\usage{
region_dispatch(level, all_levels, region_names, region_codes)
}
\arguments{
\item{level}{A character string indicating the current level.}

\item{all_levels}{A character vector indicating all the levels supported.}

\item{region_names}{A named list of region names named after the levels
supported.}

\item{region_codes}{A named list of region codes named after the levels
supported.}
}
\description{
Controls the grouping variables used in
\code{process_internal} based on the supported regions present in the
class.
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{lithuania_codes}
\alias{lithuania_codes}
\title{Region Codes for Lithuania Dataset.}
\format{
An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 61 rows and 6 columns.
}
\usage{
lithuania_codes
}
\value{
A tibble of region codes and related information,
including ISO 3166:2 codes for counties (apskritis)
and municipalities (savivaldybe), and noting which
municipalities are city municipalities or regional
municipalities.
}
\description{
The region codes for Lithuania
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{return_data}
\alias{return_data}
\title{Control data return}
\usage{
return_data(obj, class = FALSE)
}
\arguments{
\item{obj}{A Class based on a \code{DataClass}}

\item{class}{Logical, defaults to FALSE. If TRUE returns the
\code{DataClass} object rather than a tibble or a list of tibbles.
Overrides \code{steps}.}
}
\description{
Controls data return for \code{get_reigonal_data} and
\code{get_national_data}
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/USA.R
\name{USA}
\alias{USA}
\title{USA Class for downloading, cleaning and processing notification data}
\source{
\url{https://github.com/nytimes/covid-19-data/}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for USA.
}
\examples{
\dontrun{
region <- USA$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{USA}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{level_data_urls}}{List of named links to raw data that are level
specific.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{USA$set_region_codes()}}
\item \href{#method-clean_level_1}{\code{USA$clean_level_1()}}
\item \href{#method-clean_level_2}{\code{USA$clean_level_2()}}
\item \href{#method-clone}{\code{USA$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean_common">}\href{../../covidregionaldata/html/DataClass.html#method-clean_common}{\code{covidregionaldata::DataClass$clean_common()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{USA$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
State Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{USA$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_2"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_2}{}}}
\subsection{Method \code{clean_level_2()}}{
County Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{USA$clean_level_2()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{USA$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/Google.R
\name{Google}
\alias{Google}
\title{R6 Class containing specific attributes and methods for Google data}
\source{
\url{https://github.com/GoogleCloudPlatform/covid-19-open-data}
}
\description{
Google specific information for downloading, cleaning
and processing covid-19 region data for an example Country. The function
works the same as other national data sources, however, data from
Google supports three subregions (country, subregion and subregion2) which
can be accessed using the 'level' argument. There is also more data
available, such as hospitalisations data. The raw data comes as three
seperate data sets, "epidemiology" which is comprised of cases, tests and
deaths, "index", which holds information about countries linking the other
data sets, and "hospitalizations" which holds data about number of people
in hospital, ICU, etc.
}
\examples{
# nolint start
\dontrun{
# set up a data cache
start_using_memoise()

# get all countries
national <- Google$new(level = "1", get = TRUE)
national$return()

# show available regions with data at the second level of interest
google_level_2 <- Google$new(level = "2")
google_level_2$download()
google_level_2$clean()
google$available_regions()

# get all region data for the uk
google_level_2$filter("uk")
google_level_2$process()
google_level_2$return()

# get all regional data for the UK
uk <- Google$new(regions = "uk", level = "2", get = TRUE)
uk$return()

# get all subregional data for the UK
uk <- Google$new(regions = "uk", level = "3", get = TRUE)
uk$return()
}
# nolint end
}
\seealso{
Aggregated data sources
\code{\link{Covid19DataHub}},
\code{\link{JHU}}

National data sources
\code{\link{Covid19DataHub}},
\code{\link{ECDC}},
\code{\link{JHU}},
\code{\link{JRC}},
\code{\link{WHO}}

Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{aggregations}
\concept{dataset}
\concept{national}
\concept{subnational}
\section{Super classes}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{\link[covidregionaldata:CountryDataClass]{covidregionaldata::CountryDataClass}} -> \code{Google}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of country to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-clean_common}{\code{Google$clean_common()}}
\item \href{#method-clean_level_1}{\code{Google$clean_level_1()}}
\item \href{#method-clean_level_2}{\code{Google$clean_level_2()}}
\item \href{#method-new}{\code{Google$new()}}
\item \href{#method-clone}{\code{Google$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="set_region_codes">}\href{../../covidregionaldata/html/DataClass.html#method-set_region_codes}{\code{covidregionaldata::DataClass$set_region_codes()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="CountryDataClass" data-id="filter">}\href{../../covidregionaldata/html/CountryDataClass.html#method-filter}{\code{covidregionaldata::CountryDataClass$filter()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
GoogleData specific subregion2 level data cleaning. This
takes all the raw data, puts into a single data frame, renames some
columns and checks types.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Google$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
Google specific subregion level data cleaning. Takes the
data cleaned by \code{clean_common} and aggregates it to the country level
(level 1).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Google$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_2"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_2}{}}}
\subsection{Method \code{clean_level_2()}}{
Google specific subregion2 level data cleaning. Takes the
data cleaned by \code{clean_common} and aggregates it to the subregion level
(level 2).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Google$clean_level_2()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
custom initialize for Google
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Google$new(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{arguments to be passed to \code{DataClass} and initialize Google}
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
\if{html}{\out{<div class="r">}}\preformatted{Google$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/Italy.R
\name{Italy}
\alias{Italy}
\title{Italy Class for downloading, cleaning and processing notification data}
\source{
\url{https://github.com/pcm-dpc/COVID-19/}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for Italy.
}
\examples{
\dontrun{
region <- Italy$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Italy}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data. The first, and
only entry, is be named main.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Italy$set_region_codes()}}
\item \href{#method-clean_common}{\code{Italy$clean_common()}}
\item \href{#method-clone}{\code{Italy$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Italy$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
State level data cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Italy$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Italy$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/UK.R
\name{UK}
\alias{UK}
\title{United Kingdom Class for downloading, cleaning and processing notification
data.}
\source{
\url{https://coronavirus.data.gov.uk/details/download}

\url{https://coronavirus.data.gov.uk/details/download}
}
\description{
Extracts daily COVID-19 data for the UK, stratified by region
and nation. Additional options for this class are: to return subnational
English regions using NHS region boundaries instead of PHE boundaries
(nhsregions = TRUE), a release date to download from (release_date) and a
geographical resolution (resolution).
}
\examples{
\dontrun{
# setup a data cache
start_using_memoise()

# download, clean and process level 1 UK data with hospital admissions
region <- UK$new(level = "1", nhsregions = TRUE)
region$return()

# initialise level 2 data
utla <- UK$new(level = "2")

# download UTLA data
utla$download()

# clean UTLA data
utla$clean()

# inspect available level 1 regions
utla$available_regions(level = "1")

# filter data to the East of England
utla$filter("East of England")

# process UTLA data
utla$process()

# return processed and filtered data
utla$return()

# inspect all data steps
utla$data
}

## ------------------------------------------------
## Method `UK$new`
## ------------------------------------------------

\dontrun{
UK$new(
 level = 1, localise = TRUE,
 verbose = True, steps = FALSE,
 nhsregions = FALSE, release_date = NULL,
 resolution = "utla"
)
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{UK}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data. The first, and
only entry, is be named main.}

\item{\code{level_data_urls}}{List of named links to raw data that are level
specific.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}

\item{\code{query_filters}}{Set what filters to use to query the data}

\item{\code{nhsregions}}{Whether to include NHS regions in the data}

\item{\code{release_date}}{The release date for the data}

\item{\code{resolution}}{The resolution of the data to return}

\item{\code{authority_data}}{The raw data for creating authority lookup tables}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{UK$set_region_codes()}}
\item \href{#method-download}{\code{UK$download()}}
\item \href{#method-clean_level_1}{\code{UK$clean_level_1()}}
\item \href{#method-clean_level_2}{\code{UK$clean_level_2()}}
\item \href{#method-new}{\code{UK$new()}}
\item \href{#method-download_filter}{\code{UK$download_filter()}}
\item \href{#method-set_filters}{\code{UK$set_filters()}}
\item \href{#method-download_nhs_regions}{\code{UK$download_nhs_regions()}}
\item \href{#method-add_nhs_regions}{\code{UK$add_nhs_regions()}}
\item \href{#method-specific_tests}{\code{UK$specific_tests()}}
\item \href{#method-clone}{\code{UK$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean_common">}\href{../../covidregionaldata/html/DataClass.html#method-clean_common}{\code{covidregionaldata::DataClass$clean_common()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Specific function for getting region codes for UK .
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UK$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-download"></a>}}
\if{latex}{\out{\hypertarget{method-download}{}}}
\subsection{Method \code{download()}}{
UK specific \code{download()} function.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UK$download()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
Region Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UK$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_2"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_2}{}}}
\subsection{Method \code{clean_level_2()}}{
Level 2 Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UK$clean_level_2()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initalize the UK Class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UK$new(nhsregions = FALSE, release_date = NULL, resolution = "utla", ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nhsregions}}{Return subnational English regions using NHS region
boundaries instead of PHE boundaries.}

\item{\code{release_date}}{Date data was released. Default is to extract
latest release. Dates should be in the format "yyyy-mm-dd".}

\item{\code{resolution}}{"utla" (default) or "ltla", depending on which
geographical resolution is preferred}

\item{\code{...}}{Optional arguments passed to \code{\link[=DataClass]{DataClass()}} initalize.}
}
\if{html}{\out{</div>}}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
UK$new(
 level = 1, localise = TRUE,
 verbose = True, steps = FALSE,
 nhsregions = FALSE, release_date = NULL,
 resolution = "utla"
)
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-download_filter"></a>}}
\if{latex}{\out{\hypertarget{method-download_filter}{}}}
\subsection{Method \code{download_filter()}}{
Helper function for downloading data API
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UK$download_filter(filter)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{filter}}{region filters}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_filters"></a>}}
\if{latex}{\out{\hypertarget{method-set_filters}{}}}
\subsection{Method \code{set_filters()}}{
Set filters for UK data api query.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UK$set_filters()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-download_nhs_regions"></a>}}
\if{latex}{\out{\hypertarget{method-download_nhs_regions}{}}}
\subsection{Method \code{download_nhs_regions()}}{
Download NHS data for level 1 regions
Separate NHS data is available for "first" admissions, excluding
readmissions. This is available for England + English regions only.
Data are available separately for the periods 2020-08-01 to 2021-04-06,
and 2021-04-07 - present.
See: \url{https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/}
Section 2, "2. Estimated new hospital cases"
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UK$download_nhs_regions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nhs data.frame of nhs regions
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add_nhs_regions"></a>}}
\if{latex}{\out{\hypertarget{method-add_nhs_regions}{}}}
\subsection{Method \code{add_nhs_regions()}}{
Add NHS data for level 1 regions
Separate NHS data is available for "first" admissions, excluding
readmissions. This is available for England + English regions only.
See: \url{https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/}
Section 2, "2. Estimated new hospital cases"
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UK$add_nhs_regions(clean_data, nhs_data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{clean_data}}{Cleaned UK covid-19 data}

\item{\code{nhs_data}}{NHS region data}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specific_tests"></a>}}
\if{latex}{\out{\hypertarget{method-specific_tests}{}}}
\subsection{Method \code{specific_tests()}}{
Specific tests for UK data. In addition to generic tests ran
by \code{DataClass$test()} data for NHS regions are downloaded and ran through
the same generic checks (test_cleaning, test_processing, test_return). If
download = TRUE or a snapshot file is not found, the nhs data is
downloaded and saved to the snapshot location provided. If an existing
snapshot file is found then this data is used in the next tests.
Tests data can be downloaded, cleaned, processed and returned. Designed
to be ran from \code{test} and not ran directly.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UK$specific_tests(
  self_copy,
  download = FALSE,
  all = FALSE,
  snapshot_path = "",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{self_copy}}{R6class the object to test.}

\item{\code{download}}{logical. To download the data (TRUE) or use a snapshot
(FALSE). Defaults to FALSE.}

\item{\code{all}}{logical. Run tests with all settings (TRUE) or with those
defined in the current class instance (FALSE). Defaults to FALSE.}

\item{\code{snapshot_path}}{character_array the path to save the downloaded
snapshot to. Works on the snapshot path constructed by \code{test} but adds}

\item{\code{...}}{Additional parameters to pass to \code{specific_tests}}
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
\if{html}{\out{<div class="r">}}\preformatted{UK$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/datasets.R
\docType{data}
\name{mexico_codes}
\alias{mexico_codes}
\title{Region Codes for Mexico Dataset.}
\format{
An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 2489 rows and 4 columns.
}
\usage{
mexico_codes
}
\value{
A nested tibble of region codes and related information.
}
\description{
Details of the region codes used for the Mexico dataset.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{run_optional_processing_fns}
\alias{run_optional_processing_fns}
\title{Optional processing steps to run}
\usage{
run_optional_processing_fns(data, process_fns)
}
\arguments{
\item{data}{A data table}

\item{process_fns}{array, additional functions to be called after default
processing steps}
}
\description{
user supplied processing steps which are run after default steps
}
\seealso{
Functions used in the processing pipeline
\code{\link{process_internal}()},
\code{\link{run_default_processing_fns}()}
}
\concept{processing}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{csv_reader}
\alias{csv_reader}
\title{Custom CSV reading function}
\usage{
csv_reader(file, verbose = FALSE, guess_max = 1000, ...)
}
\arguments{
\item{file}{A URL or filepath to a CSV}

\item{verbose}{Logical, defaults to \code{TRUE}. Should verbose processing
messages and warnings be returned.}

\item{guess_max}{Maximum number of records to use for guessing column types.
Defaults to a 1000.}

\item{...}{extra parameters to be passed to vroom::vroom}
}
\value{
A data table
}
\description{
Checks for use of memoise and then uses vroom::vroom.
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{make_new_data_source}
\alias{make_new_data_source}
\title{Create new country class for a given source}
\usage{
make_new_data_source(
  source,
  type = "subnational",
  newfile_path = paste0("R/", source, ".R")
)
}
\arguments{
\item{source}{character_array The name of the class to create. Must start
with a capital letter (be upper camel case or an acronym in all caps such as
WHO).}

\item{type}{character_array the type of class to create, subnational or
National defaults to subnational. Regional classes are individual countries,
such as UK, Italy, India, etc. These inherit from \code{DataClass}, whilst
national classes are sources for multiple countries data, such as JRC, JHU,
Google, etc. These inherit from \code{CountryDataClass}.}

\item{newfile_path}{character_array the place to save the class file}
}
\description{
Makes a new regional or national country class with the name
provided as the source. This forms a basic template for the user to fill in
with the specific field values and cleaning functions required. This also
creates a github workflow file for the same country.
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_linelist.R
\name{get_linelist}
\alias{get_linelist}
\title{Get patient line list data}
\source{
\url{https://github.com/beoutbreakprepared/nCoV2019}
}
\usage{
get_linelist(clean = TRUE, report_delay_only = FALSE)
}
\arguments{
\item{clean}{Logical, defaults to \code{TRUE}.
Should the data returned be cleaned for use.}

\item{report_delay_only}{Logical, defaults to \code{FALSE}.
Should only certain variables (id, country, onset date, days' delay),
and observations (patients with a report delay) be returned}
}
\value{
A line list of reported cases of COVID-19
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}

Provides a public international patient line list from January 2020 to June
2020.

This version of the line list has stopped updating. The new version of the
line list is free but requires a login.

See: https://global.health/
}
\examples{
\dontrun{
# Get the complete linelist
get_linelist()

# Return the report delay only
get_linelist(report_delay_only = TRUE)
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Canada.R
\name{Canada}
\alias{Canada}
\title{Canada Class containing origin specific attributes and methods}
\source{
\url{https://health-infobase.canada.ca}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for Canada.
}
\examples{
\dontrun{
region <- Canada$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Canada}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data that are common
across levels.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Canada$set_region_codes()}}
\item \href{#method-clean_common}{\code{Canada$clean_common()}}
\item \href{#method-clone}{\code{Canada$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Canada$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Provincial Level Data
cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Canada$clean_common()}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{pass additional arguments}
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
\if{html}{\out{<div class="r">}}\preformatted{Canada$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/datasets.R
\docType{data}
\name{vietnam_codes}
\alias{vietnam_codes}
\title{Region Codes for Vietnam Dataset.}
\format{
An object of class \code{data.frame} with 63 rows and 2 columns.
}
\usage{
vietnam_codes
}
\value{
A tibble of region codes and related information.
}
\description{
The region codes for Viet Nam
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{complete_cumulative_columns}
\alias{complete_cumulative_columns}
\title{Completes cumulative columns if rows were added with NAs.}
\usage{
complete_cumulative_columns(data)
}
\arguments{
\item{data}{A data frame}
}
\value{
A data tibble with NAs filled in for cumulative data columns.
}
\description{
If a dataset had a row of NAs added to it (using
fill_empty_dates_with_na) then cumulative data columns will have NAs which
can cause issues later. This function fills these values with the previous
non-NA value.
}
\seealso{
Compulsory processing functions
\code{\link{add_extra_na_cols}()},
\code{\link{calculate_columns_from_existing_data}()},
\code{\link{fill_empty_dates_with_na}()}
}
\concept{compulsory_processing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-DataClass.R
\name{test_download}
\alias{test_download}
\title{Test download method works correctly}
\usage{
test_download(DataClass_obj, download, snapshot_path)
}
\arguments{
\item{DataClass_obj}{The R6Class object to perform checks on.
Must be a \code{DataClass} or \code{DataClass} child object.}

\item{download}{Logical check to download or use a snapshot of the data}

\item{snapshot_path}{character_array the path to save the downloaded
snapshot to.}
}
\description{
Test data can be downloaded if \code{download = TRUE}, or a requested
snapshot file is not found, and store a snap shot in the \code{snapshot_dir}. If
an existing snapshot file is found then load this data to use in future tests
}
\seealso{
Functions used for testing data is cleaned and processed correctly
\code{\link{expect_clean_cols}()},
\code{\link{expect_columns_contain_data}()},
\code{\link{expect_processed_cols}()},
\code{\link{test_cleaning}()},
\code{\link{test_processing}()},
\code{\link{test_return}()}
}
\concept{tests}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_regional_data.R
\name{get_regional_data}
\alias{get_regional_data}
\title{Get regional-level data}
\usage{
get_regional_data(
  country,
  level = "1",
  totals = FALSE,
  localise = TRUE,
  steps = FALSE,
  class = FALSE,
  verbose = TRUE,
  regions,
  include_level_2_regions = deprecated(),
  localise_regions = deprecated(),
  ...
)
}
\arguments{
\item{country}{A character string specifying the country to get data from.
Not case dependent. Name should be the English name. For a list of
options use \code{get_available_datasets()}.}

\item{level}{A character string indicating the target administrative level
of the data with the default being "1". Currently supported options are
level 1 ("1) and level 2 ("2"). Use \code{get_available_datasets()} for supported
options by dataset.}

\item{totals}{Logical, defaults to FALSE. If TRUE, returns totalled
data per region up to today's date. If FALSE, returns the full dataset
stratified by date and region.}

\item{localise}{Logical, defaults to TRUE. Should region names be localised.}

\item{steps}{Logical, defaults to FALSE. Should all processing and cleaning
steps be kept and output in a list.}

\item{class}{Logical, defaults to FALSE. If TRUE returns the
\code{DataClass} object rather than a tibble or a list of tibbles.
Overrides \code{steps}.}

\item{verbose}{Logical, defaults to \code{TRUE}. Should verbose processing
messages and warnings be returned.}

\item{regions}{A character vector of target regions to be assigned to the
\code{target_regions} field and used to filter the returned data.}

\item{include_level_2_regions}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Boolean. If TRUE, returns data stratified by
level 2 regions. If FALSE, stratified by Level 1. Note that Level 2 region
data is not always available. In these cases the user will get a warning
and the Level 1 data will be returned.}

\item{localise_regions}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Logical, defaults to TRUE. Should region names be localised.}

\item{...}{Additional arguments to pass to class specific functionality.}
}
\value{
A tibble with data related to cases, deaths, hospitalisations,
recoveries and testing stratified by regions within the given country.
}
\description{
Provides an interface to source specific classes which
support regional level data. For simple use cases this allows downloading
clean, standardised, regional-level COVID-19 data sets. Internally this uses
the \code{DataClass()} parent class which allows documented downloading, cleaning,
and processing. Optionally all steps of data processing can be returned
along with the functions used for processing but by default just the
finalised processed data is returned. See the examples for some potential
use cases and the links to lower level functions for more details and
options.
}
\examples{
\dontrun{
# set up a data cache
start_using_memoise()

# download data for Italy
get_regional_data("italy")

# return totals for Italy with no localisation
get_regional_data("italy", localise = FALSE, totals = TRUE)

# download data for the UK but return the class
uk <- get_regional_data("United Kingdom", class = TRUE)
uk

# return UK data from the class object]
uk$return()
}
}
\seealso{
\code{\link[=Italy]{Italy()}}, \code{\link[=UK]{UK()}}

Data interface functions
\code{\link{CountryDataClass}},
\code{\link{DataClass}},
\code{\link{get_available_datasets}()},
\code{\link{get_national_data}()},
\code{\link{initialise_dataclass}()}
}
\concept{interface}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/France.R
\name{France}
\alias{France}
\title{France Class containing origin specific attributes and methods}
\source{
\url{https://www.data.gouv.fr/fr/datasets/r/406c6a23-e283-4300-9484-54e78c8ae675}

\url{https://www.data.gouv.fr/fr/datasets/r/6fadff46-9efd-4c53-942a-54aca783c30c}

\url{https://www.data.gouv.fr/fr/datasets/r/001aca18-df6a-45c8-89e6-f82d689e6c01}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for France.
}
\examples{
\dontrun{
region <- France$new(level = "2", verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{France}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{level_data_urls}}{List of named links to raw data that are level
specific.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{France$set_region_codes()}}
\item \href{#method-clean_level_1}{\code{France$clean_level_1()}}
\item \href{#method-clean_level_2}{\code{France$clean_level_2()}}
\item \href{#method-clone}{\code{France$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean_common">}\href{../../covidregionaldata/html/DataClass.html#method-clean_common}{\code{covidregionaldata::DataClass$clean_common()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{France$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
Region Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{France$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_2"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_2}{}}}
\subsection{Method \code{clean_level_2()}}{
Department Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{France$clean_level_2()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{France$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/Lithuania.R
\name{Lithuania}
\alias{Lithuania}
\title{Lithuania Class for downloading, cleaning and processing notification data}
\source{
\url{https://hub.arcgis.com/datasets/d49a63c934be4f65a93b6273785a8449_0}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region level 1 and 2 data for Lithuania.
}
\section{OSP Data fields}{


The \href{https://osp.stat.gov.lt}{Official Statistics Portal} (OSP) provides
many data series in their table.

The full range of these vectors can be returned by setting
\code{all_osp_fields} to \code{TRUE}.

The following describes the data provided by the OSP.\tabular{ll}{
   field \tab description \cr
   \code{date} \tab the reporting day during which the events occurred or at the end of which the accounting was performed \cr
   \code{municipality_code} \code{*} \tab code of the municipality assigned to persons \cr
   \code{municipality_name} \code{+} \tab the name of the municipality assigned to the persons \cr
   \code{population} \tab population size according to the data of the beginning of 2021, according to the declared place of residence \cr
   \code{ab_pos_day} \tab Number of positive antibody test responses, days \cr
   \code{ab_neg_day} \tab Number of negative antibody test responses, days \cr
   \code{ab_tot_day} \tab Number of antibody tests, daily \cr
   \code{ab_prc_day} \tab Percentage of positive antibody test responses per day \cr
   \code{ag_pos_day} \tab Number of positive antigen test responses, daily \cr
   \code{ag_neg_day} \tab Number of negative antigen test responses, daily \cr
   \code{ag_tot_day} \tab Number of antigen tests, daily \cr
   \code{ag_prc_day} \tab Percentage of positive responses to antigen tests per day \cr
   \code{pcr_pos_day} \tab number of positive PCR test responses, daily \cr
   \code{pcr_neg_day} \tab Number of PCR test negative responses, daily \cr
   \code{pcr_tot_day} \tab number of PCR tests per day \cr
   \code{pcr_prc_day} \tab Percentage of positive PCR test responses per day \cr
   \code{dgn_pos_day} \tab Number of positive answers to diagnostic tests / tests, days \cr
   \code{dgn_neg_day} \tab Number of negative answers to diagnostic tests / tests, days \cr
   \code{dgn_prc_day} \tab Number of diagnostic examinations / tests, days \cr
   \code{dgn_tot_day} \tab Percentage of positive answers to diagnostic tests / tests per day \cr
   \code{dgn_tot_day_gmp} \tab Number of diagnostic examinations / tests of samples collected at mobile points, days \cr
   \code{daily_deaths_def1} \tab The number of new deaths per day according to the (narrowest) COVID death definition No. 1. \verb{#} \cr
   \code{daily_deaths_def2} \tab Number of new deaths per day according to COVID death definition No. 2. \verb{#} \cr
   \code{daily_deaths_def3} \tab Number of new deaths per day according to COVID death definition No. 3. \verb{#} \cr
   \code{daily_deaths_all} \tab Daily deaths in Lithuania (by date of death) \cr
   \code{incidence} + \tab Number of new COVID cases per day (laboratory or physician confirmed) \cr
   \code{cumulative_totals} + \tab Total number of COVID cases (laboratory or physician confirmed) \cr
   \code{active_de_jure} \tab Declared number of people with COVID \cr
   \code{active_sttstcl} \tab Statistical number of people with COVID \cr
   \code{dead_cases} \tab The number of dead persons who were ever diagnosed with COVID \cr
   \code{recovered_de_jure} \tab Declared number of recovered live persons \cr
   \code{recovered_sttstcl} \tab Statistical number of recovered live persons \cr
   \code{map_colors} \code{$} \tab The map colour-coding for the municipality, based on averages of test positivity and incidence per capita \cr
}


\code{*} The \code{municipality_code} is discarded since it does not correspond
to ISO-3166:2 codes used elsewhere in the package.

\code{+} These fields are renamed but returned unmodified.

\verb{#} Lithuania offers counts according to three
different definitions of whether a death is attributable to COVID-19.

\code{$} This field is not recalculated for counties and is deleted.
}

\section{Criteria for attributing deaths}{


Beginning in February 2021 the OSP publishes death counts according to
three different criteria, from most to least strictly attributed to
COVID-19.
\enumerate{
\item \emph{\code{of}} Number of deaths with COVID-19 (coronavirus infection) as
the leading cause of death. The indicator is calculated by summing
all registered records of medical form E106 (unique persons), in which
the main cause of death is IPC disease codes U07.1 or U07.2. Deaths
due to external causes are not included (ICD disease codes are V00-Y36,
or Y85-Y87, or Y89, or S00-T79, or T89-T98).
\item \emph{\code{with}} Number of deaths with COVID-19 (coronavirus infection) of
any cause of death.
The indicator is calculated by summing all registered records of the
medical form E106 (unique persons), in which the ICD disease codes
U07.1, U07.2, U07.3, U07.4, U07.5 are indicated as the main, direct,
intermediate cause of death or other important pathological condition,
or identified as related to COVID-19 disease (coronavirus infection).
Deaths due to external causes are not included (ICD disease codes
are V00-Y36, or Y85-Y87, or Y89, or S00-T79, or T89-T98).
\item \emph{\code{after}} Number of deaths from any cause of COVID-19 or COVID-19
deaths due to non-external causes within 28 days.
The indicator is calculated by summing all registered records of the
medical form E106 (unique persons), in which the ICD disease codes
U07.1, U07.2, U07.3, U07.4, U07 are indicated as the main, direct,
intermediate cause of death or other important pathological condition,
or identified as related to COVID-19 disease (coronavirus infection)
and all records of medical form E106 (unique individuals) where the
person died within the last 28 days after receiving a positive
diagnostic response to the SARS-CoV-2 test or had an entry in medical
form E025 with ICD disease code U07.2 or U07.1. Deaths due to external
causes are not included (ICD disease codes are V00-Y36, or Y85-Y87, or
Y89, or S00-T79, or T89-T98).
}

The number of deaths reported in the last day is preliminary and
increases by about 20-40\% in a few days. Such a "delay" in the data is
natural: for example, for those who died last night, a death certificate
is likely to be issued as soon as this report is published this morning.
}

\section{De jure and statistical counts}{


Beginning in February 2021 the OSP makes statistical estimates
of the number of recovered and active cases, since review of the data
showed that some cases individuals still considered as active cases
had recovered, but not documented or registered as such.

These are listed as by the OSP as \code{active_de_jure} and
\code{recovered_de_jure} (officially still considered sick),
and \code{active_sttstcl} and \code{recovered_sttstcl} (an estimate of how
many of these are still ill).
}

\examples{
\dontrun{
region <- Lithuania$new(verbose = TRUE, steps = TRUE, get = TRUE)
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Lithuania}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data that are common
across levels.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}

\item{\code{death_definition}}{which criteria of deaths attributed to
COVID to use}

\item{\code{recovered_definition}}{whether to use the official counts of
recovered cases or the statistical estimates provided by OSP}

\item{\code{all_osp_fields}}{whether to return all the data vectors provided
by OSP}

\item{\code{national_data}}{whether to return data rows for national results}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Lithuania$set_region_codes()}}
\item \href{#method-clean_common}{\code{Lithuania$clean_common()}}
\item \href{#method-clean_level_1}{\code{Lithuania$clean_level_1()}}
\item \href{#method-new}{\code{Lithuania$new()}}
\item \href{#method-clone}{\code{Lithuania$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lithuania$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Common data cleaning for both levels
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lithuania$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
Lithuania Specific County Level Data Cleaning

Aggregates data to the level 1 (county) regional level. Data is
provided by the source at the level 2 (municipality) regional level.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lithuania$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize the country
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lithuania$new(
  death_definition = "of",
  recovered_definition = "official",
  all_osp_fields = FALSE,
  national_data = FALSE,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{death_definition}}{A character string. Determines which criteria
for attributing deaths to COVID is used. Should be \code{"of"},
\code{"with"}, or \code{"after"}. Can also be \code{"daily_deaths_def1"},
\code{"daily_deaths_def2"}, or \code{"daily_deaths_def3"}. (Defaults
to \code{"of"}, the strictest definition.)}

\item{\code{recovered_definition}}{A character string. Determines whether
the count of officially-recovered (\emph{de jure}) cases is used, or
the statistical estimate provided by OSP. Should be \code{"official"}
or \code{"statistical"}. (Defaults to \code{"official"}.)}

\item{\code{all_osp_fields}}{A logical scalar. Should all the meaningful
data fields from the OSP source be returned? (Defaults \code{FALSE})}

\item{\code{national_data}}{A logical scalar. Should national values be
returned?  (Defaults \code{FALSE})}

\item{\code{...}}{Parameters passed to \code{\link[=DataClass]{DataClass()}} initalize}
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
\if{html}{\out{<div class="r">}}\preformatted{Lithuania$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/utils.R
\name{glue_level}
\alias{glue_level}
\title{Glue the spatial level into a variable name}
\usage{
glue_level(level)
}
\arguments{
\item{level}{A character string indicating the current level.}
}
\value{
A string in the form "level_1_region".
}
\description{
Glue the spatial level into a variable name
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-DataClass.R
\name{expect_columns_contain_data}
\alias{expect_columns_contain_data}
\title{Test that cleaned columns contain data/}
\usage{
expect_columns_contain_data(DataClass_obj)
}
\arguments{
\item{DataClass_obj}{The DataClass object (R6Class) to perform checks on.
Must be a \code{DataClass} or \code{DataClass} child object.}
}
\description{
Checks that cleaned columns cases, deaths, recovered and test
(new and total) are not entirely composed of NAs.
}
\seealso{
Functions used for testing data is cleaned and processed correctly
\code{\link{expect_clean_cols}()},
\code{\link{expect_processed_cols}()},
\code{\link{test_cleaning}()},
\code{\link{test_download}()},
\code{\link{test_processing}()},
\code{\link{test_return}()}
}
\concept{tests}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-DataClass.R
\name{expect_clean_cols}
\alias{expect_clean_cols}
\title{Test clean columns contain the correct data and types}
\usage{
expect_clean_cols(data, level)
}
\arguments{
\item{data}{The clean data to check}

\item{level}{character_array the level of the data to check}
}
\description{
Checks the date column is an s3 class and that region level
column is a character in the cleaned data (data$clean)
}
\seealso{
Functions used for testing data is cleaned and processed correctly
\code{\link{expect_columns_contain_data}()},
\code{\link{expect_processed_cols}()},
\code{\link{test_cleaning}()},
\code{\link{test_download}()},
\code{\link{test_processing}()},
\code{\link{test_return}()}
}
\concept{tests}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-DataClass.R
\name{test_processing}
\alias{test_processing}
\title{Test process method works correctly}
\usage{
test_processing(DataClass_obj, all = FALSE)
}
\arguments{
\item{DataClass_obj}{The R6Class object to perform checks on.
Must be a \code{DataClass} or \code{DataClass} child object.}

\item{all}{Logical. Run tests with all settings (TRUE) or with those
defined in the current class instance (FALSE). Defaults to FALSE.}
}
\description{
Test data can be processed correctly using the process method.
process is invoked to generate processed data which is then checked to ensure
it is a data.frame, which is not empty, has at least 2 columns and calls
\code{expect_processed_columns} to check each column types.
}
\seealso{
Functions used for testing data is cleaned and processed correctly
\code{\link{expect_clean_cols}()},
\code{\link{expect_columns_contain_data}()},
\code{\link{expect_processed_cols}()},
\code{\link{test_cleaning}()},
\code{\link{test_download}()},
\code{\link{test_return}()}
}
\concept{tests}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cuba.R
\name{Cuba}
\alias{Cuba}
\title{Cuba Class for downloading, cleaning and processing notification data}
\source{
\url{https://covid19cubadata.github.io/}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for Cuba
}
\examples{
\dontrun{
region <- Cuba$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Cuba}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Cuba$set_region_codes()}}
\item \href{#method-clean_common}{\code{Cuba$clean_common()}}
\item \href{#method-clone}{\code{Cuba$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cuba$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Cuba specific state level data cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cuba$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cuba$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/WHO.R
\name{WHO}
\alias{WHO}
\title{R6 Class containing specific attributes and methods for World Health
Organisation data}
\source{
\url{https://covid19.who.int/}
}
\description{
Information for downloading, cleaning and processing COVID-19
region data from the World Health Organisation
}
\examples{
\dontrun{
national <- WHO$new(verbose = TRUE, steps = TRUE, get = TRUE)
national$return()
}
}
\seealso{
National data sources
\code{\link{Covid19DataHub}},
\code{\link{ECDC}},
\code{\link{Google}},
\code{\link{JHU}},
\code{\link{JRC}}
}
\concept{dataset}
\concept{national}
\section{Super classes}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{\link[covidregionaldata:CountryDataClass]{covidregionaldata::CountryDataClass}} -> \code{WHO}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data. The first, and
only entry, is be named main.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-clean_common}{\code{WHO$clean_common()}}
\item \href{#method-return}{\code{WHO$return()}}
\item \href{#method-specific_tests}{\code{WHO$specific_tests()}}
\item \href{#method-clone}{\code{WHO$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="set_region_codes">}\href{../../covidregionaldata/html/DataClass.html#method-set_region_codes}{\code{covidregionaldata::DataClass$set_region_codes()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="CountryDataClass" data-id="filter">}\href{../../covidregionaldata/html/CountryDataClass.html#method-filter}{\code{covidregionaldata::CountryDataClass$filter()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
WHO specific data cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{WHO$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-return"></a>}}
\if{latex}{\out{\hypertarget{method-return}{}}}
\subsection{Method \code{return()}}{
Specific return settings for the WHO dataset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{WHO$return()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specific_tests"></a>}}
\if{latex}{\out{\hypertarget{method-specific_tests}{}}}
\subsection{Method \code{specific_tests()}}{
Run additional tests on WHO data. Tests that there is only
one row per country. Designed to be ran from \code{test} and not ran directly.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{WHO$specific_tests(self_copy, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{self_copy}}{R6class the object to test}

\item{\code{...}}{Extra params passed to specific download functions}
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
\if{html}{\out{<div class="r">}}\preformatted{WHO$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/datasets.R
\docType{data}
\name{colombia_codes}
\alias{colombia_codes}
\title{Region Codes for Colombia Dataset.}
\format{
An object of class \code{data.frame} with 33 rows and 2 columns.
}
\usage{
colombia_codes
}
\value{
A tibble of region codes and related information.
}
\description{
The region codes for Colombia
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ECDC.R
\name{ECDC}
\alias{ECDC}
\title{R6 Class containing specific attributes and methods for the European Centre
for Disease Prevention and Control dataset}
\source{
\url{https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide}
}
\description{
Information for downloading, cleaning
and processing the European Centre for
Disease Prevention and Control COVID-19 data.
}
\examples{
\dontrun{
national <- ECDC$new(verbose = TRUE, steps = TRUE, get = TRUE)
national$return()
}

}
\seealso{
National data sources
\code{\link{Covid19DataHub}},
\code{\link{Google}},
\code{\link{JHU}},
\code{\link{JRC}},
\code{\link{WHO}}
}
\concept{dataset}
\concept{national}
\section{Super classes}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{\link[covidregionaldata:CountryDataClass]{covidregionaldata::CountryDataClass}} -> \code{ECDC}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-clean_common}{\code{ECDC$clean_common()}}
\item \href{#method-return}{\code{ECDC$return()}}
\item \href{#method-specific_tests}{\code{ECDC$specific_tests()}}
\item \href{#method-clone}{\code{ECDC$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="set_region_codes">}\href{../../covidregionaldata/html/DataClass.html#method-set_region_codes}{\code{covidregionaldata::DataClass$set_region_codes()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="CountryDataClass" data-id="filter">}\href{../../covidregionaldata/html/CountryDataClass.html#method-filter}{\code{covidregionaldata::CountryDataClass$filter()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
ECDC specific state level data cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ECDC$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-return"></a>}}
\if{latex}{\out{\hypertarget{method-return}{}}}
\subsection{Method \code{return()}}{
Specific return settings for the ECDC dataset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ECDC$return()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specific_tests"></a>}}
\if{latex}{\out{\hypertarget{method-specific_tests}{}}}
\subsection{Method \code{specific_tests()}}{
Run additional tests on ECDC class. Tests ECDC has required
additional columns and that there is only one row per country. Designed
to be run from \code{test} and not run directly.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ECDC$specific_tests(self_copy, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{self_copy}}{R6class the object to test}

\item{\code{...}}{Extra params passed to specific download functions}
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
\if{html}{\out{<div class="r">}}\preformatted{ECDC$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/get_interventions_data.R
\name{get_interventions_data}
\alias{get_interventions_data}
\title{Get ACAPS Government Interventions dataset}
\source{
\url{https://www.acaps.org/covid-19-government-measures-dataset}
}
\usage{
get_interventions_data()
}
\value{
a dataframe of government interventions up to Dec 2020 from ACAPS
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}

Downloads the ACAPS Government Interventions dataset.
This function is deprecated: data are no longer updated as of December 2020.

Over 100 alternative datasets are available, covering government
interventions worldwide. Several include subnational level policy.
See: https://supertracker.spi.ox.ac.uk/policy-trackers/
}
\author{
Paul Campbell @paulcampbell91
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estonia.R
\name{Estonia}
\alias{Estonia}
\title{Estonia Class for downloading, cleaning and processing notification data}
\source{
\url{https://www.terviseamet.ee/et/koroonaviirus/avaandmed}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for Estonia
}
\examples{
\dontrun{
region <- Estonia$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Estonia}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Estonia$set_region_codes()}}
\item \href{#method-clean_common}{\code{Estonia$clean_common()}}
\item \href{#method-clone}{\code{Estonia$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Estonia$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Estonia specific state level data cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Estonia$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Estonia$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/utils.R
\name{reset_cache}
\alias{reset_cache}
\title{Reset Cache and Update all Local Data}
\usage{
reset_cache()
}
\value{
Null
}
\description{
Reset Cache and Update all Local Data
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{message_verbose}
\alias{message_verbose}
\title{Wrapper for message}
\usage{
message_verbose(verbose = TRUE, ...)
}
\arguments{
\item{verbose}{Logical, defaults to \code{TRUE}. Should verbose processing
messages and warnings be returned.}

\item{...}{Additional arguments passed to \code{message}.}
}
\description{
A wrapper for \code{message} that only prints output when
\code{verbose = TRUE}.
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Vietnam.R
\name{Vietnam}
\alias{Vietnam}
\title{Vietnam Class for downloading, cleaning and processing
notification data}
\source{
\url{https://covid19.ncsc.gov.vn}
}
\description{
Information for downloading, cleaning
and processing covid-19 region data for Vietnam.
}
\examples{
\dontrun{
region <- Vietnam$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Vietnam}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of country to fetch data for}

\item{\code{supported_levels}}{List of supported levels.}

\item{\code{supported_region_names}}{List of region names in order of level.}

\item{\code{supported_region_codes}}{List of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Vietnam$set_region_codes()}}
\item \href{#method-download}{\code{Vietnam$download()}}
\item \href{#method-clean_common}{\code{Vietnam$clean_common()}}
\item \href{#method-clone}{\code{Vietnam$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Vietnam$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-download"></a>}}
\if{latex}{\out{\hypertarget{method-download}{}}}
\subsection{Method \code{download()}}{
Download function to get raw data. Uses the
parent class JSON-specific method for downloads.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Vietnam$download()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Provincial Level Data
cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Vietnam$clean_common()}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{pass additional arguments}
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
\if{html}{\out{<div class="r">}}\preformatted{Vietnam$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/test-DataClass.R
\name{expect_processed_cols}
\alias{expect_processed_cols}
\title{Test that processed columns contain the correct data and types}
\usage{
expect_processed_cols(data, level = "1", localised = TRUE)
}
\arguments{
\item{data}{The data to check}

\item{level}{character_array the level of the data to check}

\item{localised}{logical to check localised data or not, defaults to
TRUE.}
}
\description{
Checks that processed data columns date, cases_new, cases_total,
deaths_new, deaths_total and that region level have the correct types.
}
\seealso{
Functions used for testing data is cleaned and processed correctly
\code{\link{expect_clean_cols}()},
\code{\link{expect_columns_contain_data}()},
\code{\link{test_cleaning}()},
\code{\link{test_download}()},
\code{\link{test_processing}()},
\code{\link{test_return}()}
}
\concept{tests}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_available_datasets.R
\name{get_available_datasets}
\alias{get_available_datasets}
\title{Get supported data sets}
\usage{
get_available_datasets(type, render = FALSE, namespace = "covidregionaldata")
}
\arguments{
\item{type}{A character vector indicating the types of data to
return. Current options include "national" (which are datasets at the
national level which inherit from \code{CountryDataClass}) and
"regional" (which are datasets at the regional level which inherit
directly from \code{DataClass()}).}

\item{render}{Logical If TRUE the supported data set table is built from the
available classes using \code{summary} methods. If FALSE the supported
data set table is taken from package data. Defaults to FALSE.}

\item{namespace}{Character string The name of the namespace to search for
class objects. Defaults to "covidregionaldata" as the package.}
}
\value{
A list of available data sets and the spatial aggregation data is
available for.
}
\description{
Returns data on what countries are available from
the data provided with this package either using a cached dataset or built
by searching the target namespace.
}
\examples{
# see all available datasets
get_available_datasets()

# see only national level datasets
get_available_datasets("national")

# see only regional level datasets
get_available_datasets("regional")

# render the data
get_available_datasets(render = TRUE)
}
\seealso{
Data interface functions
\code{\link{CountryDataClass}},
\code{\link{DataClass}},
\code{\link{get_national_data}()},
\code{\link{get_regional_data}()},
\code{\link{initialise_dataclass}()}
}
\concept{interface}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_level}
\alias{check_level}
\title{Checks a given level is supported}
\usage{
check_level(level, supported_levels)
}
\arguments{
\item{level}{A character string indicating the current level.}

\item{supported_levels}{A character vector of supported levels}
}
\description{
Checks a given level is supported
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Brazil.R
\name{Brazil}
\alias{Brazil}
\title{Brazil Class for downloading, cleaning and processing notification data}
\source{
\url{https://github.com/wcota/covid19br}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for Brazil.

Data available on Github, curated by Wesley Cota:
DOI 10.1590/SciELOPreprints.362
}
\examples{
\dontrun{
region <- Brazil$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Brazil}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data. Data is
available at the city level and is aggregated to provide state data.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Brazil$set_region_codes()}}
\item \href{#method-clean_common}{\code{Brazil$clean_common()}}
\item \href{#method-clean_level_1}{\code{Brazil$clean_level_1()}}
\item \href{#method-clean_level_2}{\code{Brazil$clean_level_2()}}
\item \href{#method-clone}{\code{Brazil$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Brazil$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Common data cleaning for both levels
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Brazil$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
State Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Brazil$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_2"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_2}{}}}
\subsection{Method \code{clean_level_2()}}{
City Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Brazil$clean_level_2()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Brazil$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/processing.R
\name{add_extra_na_cols}
\alias{add_extra_na_cols}
\title{Add extra columns filled with NA values to a dataset.}
\usage{
add_extra_na_cols(data)
}
\arguments{
\item{data}{A data frame}
}
\value{
A tibble with relevant NA columns added
}
\description{
Adds extra columns filled with NAs to a dataset.
This ensures that all datasets from the covidregionaldata package return
datasets of the same underlying structure (i.e. same columns).
}
\seealso{
Compulsory processing functions
\code{\link{calculate_columns_from_existing_data}()},
\code{\link{complete_cumulative_columns}()},
\code{\link{fill_empty_dates_with_na}()}
}
\concept{compulsory_processing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-DataClass.R
\name{test_cleaning}
\alias{test_cleaning}
\title{Test clean method works correctly}
\usage{
test_cleaning(DataClass_obj)
}
\arguments{
\item{DataClass_obj}{The R6Class object to perform checks on.
Must be a \code{DataClass} or \code{DataClass} child object.}
}
\description{
Test data can be cleaned properly. The clean method is invoked
to generate clean data. This data is checked to ensure it is a data.frame,
is not empty, has at least two columns and that columns are clean by calling
\code{expect_clean_cols}. Also tests that \code{avaliable_regions()} are not NA and
they are all characters.
}
\seealso{
Functions used for testing data is cleaned and processed correctly
\code{\link{expect_clean_cols}()},
\code{\link{expect_columns_contain_data}()},
\code{\link{expect_processed_cols}()},
\code{\link{test_download}()},
\code{\link{test_processing}()},
\code{\link{test_return}()}
}
\concept{tests}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{run_default_processing_fns}
\alias{run_default_processing_fns}
\title{Default processing steps to run}
\usage{
run_default_processing_fns(data)
}
\arguments{
\item{data}{A data table}
}
\description{
The default processing steps to which are always run. Runs on
clean data
}
\seealso{
Functions used in the processing pipeline
\code{\link{process_internal}()},
\code{\link{run_optional_processing_fns}()}
}
\concept{processing}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shared-methods.R
\name{CountryDataClass}
\alias{CountryDataClass}
\title{R6 Class containing national level methods}
\description{
Acts as parent class for national data classes, allowing them
to access general methods defined in \code{\link[=DataClass]{DataClass()}} but with additional
}
\details{
On top of the methods documented in \code{\link[=DataClass]{DataClass()}}, this class
implements a custom filter function that supports partial matching to
English country names using the \code{countrycode} package.
}
\seealso{
Data interface functions
\code{\link{DataClass}},
\code{\link{get_available_datasets}()},
\code{\link{get_national_data}()},
\code{\link{get_regional_data}()},
\code{\link{initialise_dataclass}()}
}
\concept{interface}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{CountryDataClass}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{filter_level}}{Character The level of the data to filter at.
Defaults to the country level of the data.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-filter}{\code{CountryDataClass$filter()}}
\item \href{#method-clone}{\code{CountryDataClass$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean_common">}\href{../../covidregionaldata/html/DataClass.html#method-clean_common}{\code{covidregionaldata::DataClass$clean_common()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="set_region_codes">}\href{../../covidregionaldata/html/DataClass.html#method-set_region_codes}{\code{covidregionaldata::DataClass$set_region_codes()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-filter"></a>}}
\if{latex}{\out{\hypertarget{method-filter}{}}}
\subsection{Method \code{filter()}}{
Filter method for country level data. Uses \code{countryname}
to match input countries with known names.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CountryDataClass$filter(countries, level)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{countries}}{A character vector of target countries. Overrides the
current class setting for \code{target_regions}. If the \code{filter_level} field
\code{level} argument is set to anything other than level 1 this is passed
directly to the parent \code{DataClass()} \code{filter()} method with no
alteration.}

\item{\code{level}}{Character The level of the data to filter at. Defaults to
the conuntry level if not specified.}
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
\if{html}{\out{<div class="r">}}\preformatted{CountryDataClass$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/get_info_covidregionaldata.R
\name{get_info_covidregionaldata}
\alias{get_info_covidregionaldata}
\title{Get available datasets}
\usage{
get_info_covidregionaldata()
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}

This function is deprecated. Please use \code{get_available_datasets()} instead.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{stop_using_memoise}
\alias{stop_using_memoise}
\title{Stop using useMemoise}
\usage{
stop_using_memoise()
}
\description{
Sets useMemoise in options to NULL, meaning memoise isn't used
when reading data in
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{fill_empty_dates_with_na}
\alias{fill_empty_dates_with_na}
\title{Add rows of NAs for dates where a region does not have any data}
\usage{
fill_empty_dates_with_na(data)
}
\arguments{
\item{data}{A data frame}
}
\value{
A tibble with rows of NAs added.
}
\description{
There are points, particularly early during data collection,
where data was not collected for all regions. This function finds dates
which have data for some regions, but not all, and adds rows of NAs for the
missing regions. This is mainly for reasons of completeness.
}
\seealso{
Compulsory processing functions
\code{\link{add_extra_na_cols}()},
\code{\link{calculate_columns_from_existing_data}()},
\code{\link{complete_cumulative_columns}()}
}
\concept{compulsory_processing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Netherlands.R
\name{Netherlands}
\alias{Netherlands}
\title{Netherlands Class for downloading, cleaning and processing notification data}
\source{
\url{https://data.rivm.nl/geonetwork/srv/dut/catalog.search#/metadata/5f6bc429-1596-490e-8618-1ed8fd768427?tab=relations}
}
\description{
Class for downloading, cleaning and processing COVID-19
sub-regional data for the Netherlands, provided by RVIM (English: National
Institute for Public Health and the Environment). This data contains number
of newly reported cases (that have tested positive), number of newly reported
hospital admissions and number of newly reported deaths going back to
27/02/2020. Data is provided at both the province and municipality level.
}
\examples{
\dontrun{
region <- Netherlands$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Netherlands}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data. The first, and
only entry, is be named main.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Netherlands$set_region_codes()}}
\item \href{#method-clean_common}{\code{Netherlands$clean_common()}}
\item \href{#method-clean_level_1}{\code{Netherlands$clean_level_1()}}
\item \href{#method-clone}{\code{Netherlands$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Netherlands$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Common cleaning steps to be applied to raw data, regardless
of level (province or municipality) for raw Netherlands data.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Netherlands$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
Netherlands specific province level data cleaning. Takes
the data cleaned by \code{clean_common} and aggregates it to the Province
level (level 1).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Netherlands$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Netherlands$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/JRC.R
\name{JRC}
\alias{JRC}
\title{R6 Class containing specific attributes and methods for European Commission's
Joint Research Centre data}
\source{
\url{https://github.com/ec-jrc/COVID-19}
}
\description{
Class for downloading, cleaning and processing COVID-19
region data from the European Commission's Joint Research Centre. Subnational
data (admin level 1) on numbers of contagious and fatalities by COVID-19,
collected directly from the National Authoritative sources (National
monitoring websites, when available). For more details see
https://github.com/ec-jrc/COVID-19
}
\examples{
\dontrun{
# get country level data
jrc_level_1 <- JRC$new(level = "1", verbose = TRUE, steps = TRUE, get = TRUE)
jrc_level_1$return()

# show available regions with data at the first level of interest (country)
jrc_level_1$available_regions()

# get region level data
jrc_level_2 <- JRC$new(level = "2", verbose = TRUE, steps = TRUE, get = TRUE)
jrc_level_2$return()

# show available regions with data at the second level of interest (region)
jrc_level_2$available_regions()
}
}
\seealso{
National data sources
\code{\link{Covid19DataHub}},
\code{\link{ECDC}},
\code{\link{Google}},
\code{\link{JHU}},
\code{\link{WHO}}
}
\concept{dataset}
\concept{national}
\section{Super classes}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{\link[covidregionaldata:CountryDataClass]{covidregionaldata::CountryDataClass}} -> \code{JRC}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{level_data_urls}}{List of named links to raw data.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-clean_common}{\code{JRC$clean_common()}}
\item \href{#method-clean_level_1}{\code{JRC$clean_level_1()}}
\item \href{#method-clean_level_2}{\code{JRC$clean_level_2()}}
\item \href{#method-clone}{\code{JRC$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="set_region_codes">}\href{../../covidregionaldata/html/DataClass.html#method-set_region_codes}{\code{covidregionaldata::DataClass$set_region_codes()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="CountryDataClass" data-id="filter">}\href{../../covidregionaldata/html/CountryDataClass.html#method-filter}{\code{covidregionaldata::CountryDataClass$filter()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
JRC specific data cleaning. The raw source data columns are
converted to the correct type and renamed appropriately to match the
standard for general processing.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JRC$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
JRC specific country level data cleaning. Selects country
level (level 1) columns from the data ready for further processing.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JRC$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_2"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_2}{}}}
\subsection{Method \code{clean_level_2()}}{
JRC specific region level data cleaning. Selects country
(level 1) and region (level 2) columns from the data ready for further
processing.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JRC$clean_level_2()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JRC$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/India.R
\name{India}
\alias{India}
\title{India Class for downloading, cleaning and processing notification data}
\source{
\url{https://www.covid19india.org}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for India.
}
\examples{
\dontrun{
region <- India$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{India}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{India$set_region_codes()}}
\item \href{#method-clean_common}{\code{India$clean_common()}}
\item \href{#method-get_desired_status}{\code{India$get_desired_status()}}
\item \href{#method-clone}{\code{India$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{India$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
India state level data cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{India$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_desired_status"></a>}}
\if{latex}{\out{\hypertarget{method-get_desired_status}{}}}
\subsection{Method \code{get_desired_status()}}{
Extract data from raw table
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{India$get_desired_status(status)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{status}}{The data to extract}
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
\if{html}{\out{<div class="r">}}\preformatted{India$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/shared-methods.R
\name{initialise_dataclass}
\alias{initialise_dataclass}
\title{Initialise a child class of DataClass if it exists}
\usage{
initialise_dataclass(
  class = character(),
  level = "1",
  totals = FALSE,
  localise = TRUE,
  regions,
  verbose = TRUE,
  steps = FALSE,
  get = FALSE,
  type = c("national", "regional"),
  ...
)
}
\arguments{
\item{class}{A character string specifying the \code{DataClass()} to initialise.
Not case dependent and matching is based on either the class name or the its
country definition. For a list of options use \code{get_available_datasets()}.}

\item{level}{A character string indicating the target administrative level
of the data with the default being "1". Currently supported options are
level 1 ("1) and level 2 ("2"). Use \code{get_available_datasets()} for supported
options by dataset.}

\item{totals}{Logical, defaults to FALSE. If TRUE, returns totalled
data per region up to today's date. If FALSE, returns the full dataset
stratified by date and region.}

\item{localise}{Logical, defaults to TRUE. Should region names be localised.}

\item{regions}{A character vector of target regions to be assigned to the
\code{target_regions} field and used to filter the returned data.}

\item{verbose}{Logical, defaults to \code{TRUE}. Should verbose processing
messages and warnings be returned.}

\item{steps}{Logical, defaults to FALSE. Should all processing and cleaning
steps be kept and output in a list.}

\item{get}{Logical, defaults to FALSE. Should the class \code{get} method be
called (this will download, clean, and process data at initialisation).}

\item{type}{A character vector indicating the types of data to
return. Current options include "national" (which are datasets at the
national level which inherit from \code{CountryDataClass}) and
"regional" (which are datasets at the regional level which inherit
directly from \code{DataClass()}).}

\item{...}{Additional arguments to pass to class specific functionality.}
}
\value{
An initialised version of the target class if available,
e.g. \code{Italy()}
}
\description{
This function initialises classes based on the \code{DataClass()}
which allows documented downloading, cleaning, and processing. See the
examples for some potential use cases and the \code{DataClass()} documentation
for more details.
}
\examples{
\dontrun{
# set up a cache to store data to avoid downloading repeatedly
start_using_memoise()

# check currently available datasets
get_available_datasets()

# initialise a data set in the United Kingdom
# at the UTLA level
utla <- UK$new(level = "2")

# download UTLA data
utla$download()

# clean UTLA data
utla$clean()

# inspect available level 1 regions
utla$available_regions(level = "1")

# filter data to the East of England
utla$filter("East of England")

# process UTLA data
utla$process()

# return processed and filtered data
utla$return()

# inspect all data steps
utla$data

# initialise Italian data, download, clean and process it
italy <- initialise_dataclass("Italy", get = TRUE)
italy$return()

# initialise ECDC data, fully process it, and return totals
ecdc <- initialise_dataclass("ecdc", get = TRUE, totals = TRUE)
ecdc$return()
}
}
\seealso{
Data interface functions
\code{\link{CountryDataClass}},
\code{\link{DataClass}},
\code{\link{get_available_datasets}()},
\code{\link{get_national_data}()},
\code{\link{get_regional_data}()}
}
\concept{interface}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{start_using_memoise}
\alias{start_using_memoise}
\title{Add useMemoise to options}
\usage{
start_using_memoise(path = tempdir(), verbose = TRUE)
}
\arguments{
\item{path}{Path to cache directory, defaults to a temporary directory.}

\item{verbose}{Logical, defaults to \code{TRUE}. Should verbose processing
messages and warnings be returned.}
}
\description{
Adds useMemoise to options meaning memoise is
used when reading data in.
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{dplyr}{\code{\link[dplyr:reexports]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{make_github_workflow}
\alias{make_github_workflow}
\title{Create github action for a given source}
\usage{
make_github_workflow(
  source,
  workflow_path = paste0(".github/workflows/", source, ".yaml"),
  cron = "36 12 * * *"
)
}
\arguments{
\item{source}{character_array The name of the class to create the workflow
for.}

\item{workflow_path}{character_array The path to where the workflow file
should be saved. Defaults to '.github/workflows/'}

\item{cron}{character_array the cron time to run the tests, defaults to
36 12 * * *, following the minute, hour, day(month), month and day(week)
format.}
}
\description{
Makes a github workflow yaml file for a given source to be used
as an action to check the data as a github action.
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mexico.R
\name{Mexico}
\alias{Mexico}
\title{Meixco Class for downloading, cleaning and processing notification data}
\source{
\url{https://datos.covid-19.conacyt.mx/}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for Mexico.

Notes on region codes:

Level 1 codes = ISO-3166-2,
source: https://en.wikipedia.org/wiki/ISO_3166-2:MX

Level 2 codes = INEGI Mexican official statistics geocoding,
source: raw data

Level 1 INEGI codes are the first 2 characters of Level 2 INEGI codes
}
\examples{
\dontrun{
region <- Mexico$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Mexico}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data.}

\item{\code{level_data_urls}}{List of named links to raw data that are level
specific.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Mexico$set_region_codes()}}
\item \href{#method-download}{\code{Mexico$download()}}
\item \href{#method-clean_common}{\code{Mexico$clean_common()}}
\item \href{#method-clean_level_1}{\code{Mexico$clean_level_1()}}
\item \href{#method-clean_level_2}{\code{Mexico$clean_level_2()}}
\item \href{#method-clone}{\code{Mexico$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Mexico$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-download"></a>}}
\if{latex}{\out{\hypertarget{method-download}{}}}
\subsection{Method \code{download()}}{
Data \code{download()} function for Mexico data. This replaces
the generic download function in \code{\link[=DataClass]{DataClass()}}. To get the latest data
use a PHP script from the website.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Mexico$download()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Common Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Mexico$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
Estados Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Mexico$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_2"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_2}{}}}
\subsection{Method \code{clean_level_2()}}{
Municipality Level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Mexico$clean_level_2()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Mexico$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/processing.R
\name{set_negative_values_to_zero}
\alias{set_negative_values_to_zero}
\title{Set negative data to 0}
\usage{
set_negative_values_to_zero(data)
}
\arguments{
\item{data}{A data frame}
}
\value{
A data frame with all relevant data > 0.
}
\description{
Set data values to 0 if they are negative in a dataset. Data in
the datasets should always be > 0.
}
\seealso{
Optional processing function
\code{\link{totalise_data}()}
}
\concept{optional_processing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{all_country_data}
\alias{all_country_data}
\title{Table of available datasets along with level and other information.
Rendered from the individual R6 class objects included in this package.}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 24 rows and 10 columns.
}
\usage{
all_country_data
}
\value{
A tibble of available datasets and related information.
}
\description{
Available datasets
}
\concept{utility}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Colombia.R
\name{Colombia}
\alias{Colombia}
\title{Colombia Class for downloading, cleaning and processing notification data}
\source{
\url{https://github.com/danielcs88/colombia_covid-19/}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for Colombia
}
\examples{
\dontrun{
region <- Colombia$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Colombia}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Colombia$set_region_codes()}}
\item \href{#method-clean_common}{\code{Colombia$clean_common()}}
\item \href{#method-clone}{\code{Colombia$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Colombia$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Colombia specific state level data cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Colombia$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Colombia$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/datasets.R
\docType{data}
\name{JHU_codes}
\alias{JHU_codes}
\title{Region Codes for JHU Dataset. Taken from the region codes provided as
part of the WHO dataset.}
\format{
An object of class \code{spec_tbl_df} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 4193 rows and 2 columns.
}
\usage{
JHU_codes
}
\value{
A tibble of region codes and related information.
}
\description{
The region codes for JHU
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{france_codes}
\alias{france_codes}
\title{Region Codes for France Dataset.}
\format{
An object of class \code{data.frame} with 104 rows and 5 columns.
}
\usage{
france_codes
}
\value{
A tibble of region codes and related information.
}
\description{
The region codes for France
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Switzerland.R
\name{Switzerland}
\alias{Switzerland}
\title{Switzerland Class for downloading, cleaning and processing notification data}
\description{
Information for downloading, cleaning
and processing COVID-19 region data for Switzerland
}
\examples{
\dontrun{
region <- Switzerland$new(verbose = TRUE, steps = TRUE, get = TRUE)
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Switzerland}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data. This url links
to a JSON file which provides the addresses for the most recently-updated
CSV files, which are then downloaded.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Switzerland$set_region_codes()}}
\item \href{#method-download}{\code{Switzerland$download()}}
\item \href{#method-clean_common}{\code{Switzerland$clean_common()}}
\item \href{#method-clone}{\code{Switzerland$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Switzerland$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-download"></a>}}
\if{latex}{\out{\hypertarget{method-download}{}}}
\subsection{Method \code{download()}}{
Download function to get raw data. Downloads
the updated list of CSV files using \code{download_JSON}, filters
that to identify the required CSV files, then uses the parent
method \code{download} to download the CSV files.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Switzerland$download()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
Switzerland specific state level data cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Switzerland$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Switzerland$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/JHU.R
\name{JHU}
\alias{JHU}
\title{R6 Class containing specific attributes and methods for John Hopkins
University data}
\source{
\url{https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data}
}
\description{
Attributes and methods for COVID-19 data used for the 2019
Novel Coronavirus Visual Dashboard operated by the Johns Hopkins University
Center for Systems Science and Engineering (JHU CSSE). Supported by ESRI
Living Atlas Team and the Johns Hopkins University Applied Physics Lab
(JHU APL)
}
\details{
This dataset support both national and subnational data sources
with national level data returned by default. Subnational data is supported
for a subset of countries which can be found after cleaning using the
\code{available_regions()} method, see the examples for more details. These data
sets are sourced, cleaned, standardised by the JHU team so please see the
source repository for further details. Note that unlike many other data sets
this means methods applied to this source are not being applied to raw
surveillance data but instead to already cleaned data. If using for
analysis checking the JHU source for further details is advisable.

If using this data please cite:
"Dong E, Du H, Gardner L. An interactive web-based dashboard to track
COVID-19 in real time.
Lancet Inf Dis. 20(5):533-534. doi: 10.1016/S1473-3099(20)30120-1"
}
\examples{
# nolint start
\dontrun{
# set up a data cache
start_using_memoise()

# get all countries data
jhu <- JHU$new(level = "1", get = TRUE)
jhu$return()

# show available regions with data at the second level of interest
jhu_level_2 <- JHU$new(level = "2")
jhu_level_2$download()
jhu_level_2$clean()
jhu$available_regions()

# get all region data for the uk
jhu_level_2$filter("uk")
jhu_level_2$process()
jhu_level_2$return()
}
# nolint end
}
\seealso{
Aggregated data sources
\code{\link{Covid19DataHub}},
\code{\link{Google}}

National data sources
\code{\link{Covid19DataHub}},
\code{\link{ECDC}},
\code{\link{Google}},
\code{\link{JRC}},
\code{\link{WHO}}

Subnational data sources
\code{\link{Belgium}},
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{aggregations}
\concept{dataset}
\concept{national}
\concept{subnational}
\section{Super classes}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{\link[covidregionaldata:CountryDataClass]{covidregionaldata::CountryDataClass}} -> \code{JHU}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of country to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.}

\item{\code{common_data_urls}}{List of named links to raw data. The first, and
only entry, is be named main.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{JHU$set_region_codes()}}
\item \href{#method-clean_common}{\code{JHU$clean_common()}}
\item \href{#method-clean_level_1}{\code{JHU$clean_level_1()}}
\item \href{#method-clone}{\code{JHU$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download">}\href{../../covidregionaldata/html/DataClass.html#method-download}{\code{covidregionaldata::DataClass$download()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="CountryDataClass" data-id="filter">}\href{../../covidregionaldata/html/CountryDataClass.html#method-filter}{\code{covidregionaldata::CountryDataClass$filter()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JHU$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_common"></a>}}
\if{latex}{\out{\hypertarget{method-clean_common}{}}}
\subsection{Method \code{clean_common()}}{
JHU specific data cleaning. Joins the raw data sets, checks
column types and renames where needed.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JHU$clean_common()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
JHU specific country level data cleaning. Aggregates the
data to the country (level 2) level.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JHU$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JHU$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/processing.R
\name{totalise_data}
\alias{totalise_data}
\title{Get totals data given the time series data.}
\usage{
totalise_data(data)
}
\arguments{
\item{data}{A data table}
}
\value{
A data table, totalled up
}
\description{
Get totals data given the time series data.
}
\seealso{
Optional processing function
\code{\link{set_negative_values_to_zero}()}
}
\concept{optional_processing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{process_internal}
\alias{process_internal}
\title{Internal Shared Regional Dataset Processing}
\usage{
process_internal(
  clean_data,
  level,
  group_vars,
  totals = FALSE,
  localise = TRUE,
  verbose = TRUE,
  process_fns
)
}
\arguments{
\item{clean_data}{The clean data for a class, e.g. \code{Italy$data$clean}}

\item{level}{The level of the data, e.g. 'level_1_region'}

\item{group_vars}{Grouping variables, used to
for grouping and to localise names. It is assumed that the first entry
indicates the main region variable and the second indicates the geocode for
this variable.}

\item{totals}{Logical, defaults to \code{FALSE}. If `TRUE``, returns totalled
data per region up to today's date. If FALSE, returns the full dataset
stratified by date and region.}

\item{localise}{Logical, defaults to \code{TRUE}. Should region names be
localised.}

\item{verbose}{Logical, defaults to \code{TRUE}. Should verbose processing
messages and warnings be returned.}

\item{process_fns}{array, additional functions to be called after default
processing steps}
}
\description{
Internal shared regional data cleaning designed to be called
by \code{process}.
}
\seealso{
Functions used in the processing pipeline
\code{\link{run_default_processing_fns}()},
\code{\link{run_optional_processing_fns}()}
}
\concept{processing}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{uk_codes}
\alias{uk_codes}
\title{Region Codes for UK Dataset.}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 429 rows and 4 columns.
}
\usage{
uk_codes
}
\value{
A tibble of region codes and related information.
}
\description{
The uk authority look table for providing region codes used for
level 2 UK data.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{json_reader}
\alias{json_reader}
\title{Custom JSON reading function}
\usage{
json_reader(file, verbose = FALSE, ...)
}
\arguments{
\item{file}{A URL or filepath to a JSON}

\item{verbose}{Logical, defaults to \code{TRUE}. Should verbose processing
messages and warnings be returned.}

\item{...}{extra parameters to be passed to jsonlite::fromJSON}
}
\value{
A data table
}
\description{
Checks for use of memoise and then uses vroom::vroom.
}
\concept{utility}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Belgium.R
\name{Belgium}
\alias{Belgium}
\title{Belgium Class for downloading, cleaning and processing notification data}
\source{
\url{https://epistat.sciensano.be/Data/COVID19BE_CASES_AGESEX.csv}
}
\description{
Information for downloading, cleaning
and processing COVID-19 region level 1 and 2 data for Belgium.
}
\examples{
\dontrun{
region <- Belgium$new(verbose = TRUE, steps = TRUE, get = TRUE, level = "2")
region$return()
}
}
\seealso{
Subnational data sources
\code{\link{Brazil}},
\code{\link{Canada}},
\code{\link{Colombia}},
\code{\link{Covid19DataHub}},
\code{\link{Cuba}},
\code{\link{Estonia}},
\code{\link{France}},
\code{\link{Germany}},
\code{\link{Google}},
\code{\link{India}},
\code{\link{Italy}},
\code{\link{JHU}},
\code{\link{Lithuania}},
\code{\link{Mexico}},
\code{\link{Netherlands}},
\code{\link{SouthAfrica}},
\code{\link{Switzerland}},
\code{\link{UK}},
\code{\link{USA}},
\code{\link{Vietnam}}
}
\concept{dataset}
\concept{subnational}
\section{Super class}{
\code{\link[covidregionaldata:DataClass]{covidregionaldata::DataClass}} -> \code{Belgium}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{origin}}{name of origin to fetch data for}

\item{\code{supported_levels}}{A list of supported levels.}

\item{\code{supported_region_names}}{A list of region names in order of level.}

\item{\code{supported_region_codes}}{A list of region codes in order of level.
ISO 3166-2 codes are used for both region and province levels in
Belgium, and for provinces these are marked as being
\code{iso_3166_2_province}}

\item{\code{common_data_urls}}{List of named links to raw data that are common
across levels.}

\item{\code{level_data_urls}}{List of named links to raw data specific to
each level of regions. For Belgium, there are only additional data for
level 1 regions.}

\item{\code{source_data_cols}}{existing columns within the raw data}

\item{\code{source_text}}{Plain text description of the source of the data}

\item{\code{source_url}}{Website address for explanation/introduction of the
data}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-set_region_codes}{\code{Belgium$set_region_codes()}}
\item \href{#method-download}{\code{Belgium$download()}}
\item \href{#method-clean_level_1}{\code{Belgium$clean_level_1()}}
\item \href{#method-clean_level_2}{\code{Belgium$clean_level_2()}}
\item \href{#method-clone}{\code{Belgium$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="available_regions">}\href{../../covidregionaldata/html/DataClass.html#method-available_regions}{\code{covidregionaldata::DataClass$available_regions()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean">}\href{../../covidregionaldata/html/DataClass.html#method-clean}{\code{covidregionaldata::DataClass$clean()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="clean_common">}\href{../../covidregionaldata/html/DataClass.html#method-clean_common}{\code{covidregionaldata::DataClass$clean_common()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="download_JSON">}\href{../../covidregionaldata/html/DataClass.html#method-download_JSON}{\code{covidregionaldata::DataClass$download_JSON()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="filter">}\href{../../covidregionaldata/html/DataClass.html#method-filter}{\code{covidregionaldata::DataClass$filter()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="get">}\href{../../covidregionaldata/html/DataClass.html#method-get}{\code{covidregionaldata::DataClass$get()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="initialize">}\href{../../covidregionaldata/html/DataClass.html#method-initialize}{\code{covidregionaldata::DataClass$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="process">}\href{../../covidregionaldata/html/DataClass.html#method-process}{\code{covidregionaldata::DataClass$process()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="return">}\href{../../covidregionaldata/html/DataClass.html#method-return}{\code{covidregionaldata::DataClass$return()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="summary">}\href{../../covidregionaldata/html/DataClass.html#method-summary}{\code{covidregionaldata::DataClass$summary()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="covidregionaldata" data-topic="DataClass" data-id="test">}\href{../../covidregionaldata/html/DataClass.html#method-test}{\code{covidregionaldata::DataClass$test()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_region_codes"></a>}}
\if{latex}{\out{\hypertarget{method-set_region_codes}{}}}
\subsection{Method \code{set_region_codes()}}{
Set up a table of region codes for clean data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Belgium$set_region_codes()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-download"></a>}}
\if{latex}{\out{\hypertarget{method-download}{}}}
\subsection{Method \code{download()}}{
Downloads data from source and (for Belgium)
applies an initial data patch.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Belgium$download()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_1"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_1}{}}}
\subsection{Method \code{clean_level_1()}}{
Region-level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Belgium$clean_level_1()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clean_level_2"></a>}}
\if{latex}{\out{\hypertarget{method-clean_level_2}{}}}
\subsection{Method \code{clean_level_2()}}{
Province-level Data Cleaning
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Belgium$clean_level_2()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Belgium$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/get_national_data.R
\name{get_national_data}
\alias{get_national_data}
\title{Get national-level data for countries globally from a range of sources}
\usage{
get_national_data(
  countries,
  source = "who",
  level = "1",
  totals = FALSE,
  steps = FALSE,
  class = FALSE,
  verbose = TRUE,
  country = deprecated(),
  ...
)
}
\arguments{
\item{countries}{A character vector specifying country names of interest.
Used to filter the data.}

\item{source}{A character string specifying the data source (not case
dependent). Defaults to WHO (the World Health Organisation). See
\code{get_available_datasets("national")} for all options.}

\item{level}{A character string indicating the target administrative level
of the data with the default being "1". Currently supported options are
level 1 ("1) and level 2 ("2"). Use \code{get_available_datasets()} for supported
options by dataset.}

\item{totals}{Logical, defaults to FALSE. If TRUE, returns totalled
data per region up to today's date. If FALSE, returns the full dataset
stratified by date and region.}

\item{steps}{Logical, defaults to FALSE. Should all processing and cleaning
steps be kept and output in a list.}

\item{class}{Logical, defaults to FALSE. If TRUE returns the
\code{DataClass} object rather than a tibble or a list of tibbles.
Overrides \code{steps}.}

\item{verbose}{Logical, defaults to \code{TRUE}. Should verbose processing
messages and warnings be returned.}

\item{country}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} A character string
specifying a country to filter for.}

\item{...}{Additional arguments to pass to class specific functionality.}
}
\value{
A tibble with data related to cases, deaths, hospitalisations,
recoveries and testing.
}
\description{
Provides an interface to source specific classes which
support national level data. For simple use cases this allows downloading
clean, standardised, national-level COVID-19 data sets. Internally this uses
the \code{CountryDataClass()} parent class which allows documented downloading,
cleaning, and processing. Optionally all steps of data processing can be
returned along with the functions used for processing but by default just
the finalised processed data is returned. See the examples for some
potential use cases and the links to lower level functions for more details
and options.
}
\examples{
\dontrun{
# set up a data cache
start_using_memoise()

# download all national data from the WHO
get_national_data(source = "who")

# download data for Canada keeping all processing steps
get_national_data(countries = "canada", source = "ecdc")

# download data for Canada from the JHU and return the full class
jhu <- get_national_data(countries = "canada", source = "jhu", class = TRUE)
jhu

# return the JHU data for canada
jhu$return()

# check which regions the JHU supports national data for
jhu$available_regions()

# filter instead for France (and then reprocess)
jhu$filter("France")
jhu$process()

# explore the structure of the stored JHU data
jhu$data
}
}
\seealso{
\code{\link[=WHO]{WHO()}}, \code{\link[=ECDC]{ECDC()}}, \code{\link[=JHU]{JHU()}}, \code{\link[=Google]{Google()}}

Data interface functions
\code{\link{CountryDataClass}},
\code{\link{DataClass}},
\code{\link{get_available_datasets}()},
\code{\link{get_regional_data}()},
\code{\link{initialise_dataclass}()}
}
\concept{interface}
