# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at mes335@cornell.edu. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
<!-- README.md is generated from README.Rmd. Please edit that file -->

# auk: eBird Data Extraction and Processing in R <img src="logo.png" align="right" width=140/>

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/auk)](https://cran.r-project.org/package=auk)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/auk?color=brightgreen)](http://www.r-pkg.org/pkg/auk)
[![R-CMD-check](https://github.com/CornellLabofOrnithology/auk/workflows/R-CMD-check/badge.svg)](https://github.com/CornellLabofOrnithology/auk/actions)
[![rOpenSci](https://badges.ropensci.org/136_status.svg)](https://github.com/ropensci/onboarding/issues/136)
<!-- badges: end -->

## Overview

[eBird](http://www.ebird.org) is an online tool for recording bird
observations. Since its inception, over 600 million records of bird
sightings (i.e. combinations of location, date, time, and bird species)
have been collected, making eBird one of the largest citizen science
projects in history and an extremely valuable resource for bird research
and conservation. The full eBird database is packaged as a text file and
available for download as the [eBird Basic Dataset
(EBD)](http://ebird.org/ebird/data/download). Due to the large size of
this dataset, it must be filtered to a smaller subset of desired
observations before reading into R. This filtering is most efficiently
done using AWK, a Unix utility and programming language for processing
column formatted text data. This package acts as a front end for AWK,
allowing users to filter eBird data before import into R.

For a comprehensive resource on using eBird data for modeling species
distributions, consult the free online book [Best Practices for Using
eBird
Data](https://cornelllabofornithology.github.io/ebird-best-practices/)
and the association paper *Analytical guidelines to increase the value
of community science data: An example using eBird data to estimate
species distributions* ([Johnston et
al. 2021](https://onlinelibrary.wiley.com/doi/10.1111/ddi.13271)).

## Installation

    # cran release
    install.packages("auk")

    # or install the development version from github
    # install.packages("remotes")
    remotes::install_github("CornellLabofOrnithology/auk")

`auk` requires the Unix utility AWK, which is available on most Linux
and Mac OS X machines. Windows users will first need to install
[Cygwin](https://www.cygwin.com) before using this package. Note that
**Cygwin must be installed in the default location**
(`C:/cygwin/bin/gawk.exe` or `C:/cygwin64/bin/gawk.exe`) in order for
`auk` to work.

## Vignette

Full details on using `auk` to produce both presence-only and
presence-absence data are outlined in the
[vignette](https://cornelllabofornithology.github.io/auk/articles/auk.html).

## Cheatsheet

An `auk` cheatsheet was developed by [Mickayla
Johnston](https://www.linkedin.com/in/mickayla-johnston/):

<a href="https://github.com/CornellLabofOrnithology/auk/blob/master/cheatsheet/auk-cheatsheet.pdf"><img src="cheatsheet/auk-cheatsheet.png" width=400/></a>

## `auk` and `rebird`

Those interested in eBird data may also want to consider
[`rebird`](https://github.com/ropensci/rebird), an R package that
provides an interface to the [eBird
APIs](https://confluence.cornell.edu/display/CLOISAPI/eBirdAPIs). The
functions in `rebird` are mostly limited to accessing recent
(i.e. within the last 30 days) observations, although `ebirdfreq()` does
provide historical frequency of observation data. In contrast, `auk`
gives access to the full set of ~ 500 million eBird observations. For
most ecological applications, users will require `auk`; however, for
some use cases, e.g. building tools for birders, `rebird` provides a
quick and easy way to access data.

## A note on versions

This package contains a current (as of the time of package release)
version of the [bird taxonomy used by
eBird](http://help.ebird.org/customer/portal/articles/1006825-the-ebird-taxonomy).
This taxonomy determines the species that can be reported in eBird and
therefore the species that users of `auk` can extract. eBird releases an
updated taxonomy once a year, typically in August, at which time `auk`
will be updated to include the current taxonomy. When using `auk`, users
should be careful to ensure that the version they’re using is in sync
with the eBird Basic Dataset they’re working with. This is most easily
accomplished by always using the must recent version of `auk` and the
most recent release of the dataset.

## Quick start

This package uses the command-line program AWK to extract subsets of the
eBird Basic Dataset for use in R. This is a multi-step process:

1.  Define a reference to the eBird data file.
2.  Define a set of spatial, temporal, or taxonomic filters. Each type
    of filter corresponds to a different function, e.g. `auk_species` to
    filter by species. At this stage the filters are only set up, no
    actual filtering is done until the next step.
3.  Filter the eBird data text file, producing a new text file with only
    the selected rows.
4.  Import this text file into R as a data frame.

Because the eBird dataset is so large, step 3 typically takes several
hours to run. Here’s a simple example that extract all Canada Jay
records from within Canada.

    library(auk)
    # path to the ebird data file, here a sample included in the package
    # get the path to the example data included in the package
    # in practice, provide path to ebd, e.g. f_in <- "data/ebd_relFeb-2018.txt
    f_in <- system.file("extdata/ebd-sample.txt", package = "auk")
    # output text file
    f_out <- "ebd_filtered_grja.txt"
    ebird_data <- f_in %>% 
      # 1. reference file
      auk_ebd() %>% 
      # 2. define filters
      auk_species(species = "Canada Jay") %>% 
      auk_country(country = "Canada") %>% 
      # 3. run filtering
      auk_filter(file = f_out) %>% 
      # 4. read text file into r data frame
      read_ebd()

For those not familiar with the pipe operator (`%>%`), the above code
could be rewritten:

    f_in <- system.file("extdata/ebd-sample.txt", package = "auk")
    f_out <- "ebd_filtered_grja.txt"
    ebd <- auk_ebd(f_in)
    ebd_filters <- auk_species(ebd, species = "Canada Jay")
    ebd_filters <- auk_country(ebd_filters, country = "Canada")
    ebd_filtered <- auk_filter(ebd_filters, file = f_out)
    ebd_df <- read_ebd(ebd_filtered)

## Usage

### Filtering

`auk` uses a [pipeline-based workflow](http://r4ds.had.co.nz/pipes.html)
for defining filters, which can then be compiled into an AWK script.
Users should start by defining a reference to the dataset file with
`auk_ebd()`. Then any of the following filters can be applied:

-   `auk_species()`: filter by species using common or scientific names.
-   `auk_country()`: filter by country using the standard English names
    or [ISO 2-letter country
    codes](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2).
-   `auk_state()`: filter by state using eBird state codes, see
    `?ebird_states`.
-   `auk_bcr()`: filter by [Bird Conservation Region
    (BCR)](http://nabci-us.org/resources/bird-conservation-regions/)
    using BCR codes, see `?bcr_codes`.
-   `auk_bbox()`: filter by spatial bounding box, i.e. a range of
    latitudes and longitudes in decimal degrees.
-   `auk_date()`: filter to checklists from a range of dates. To extract
    observations from a range of dates, regardless of year, use the
    wildcard “`*`” in place of the year,
    e.g. `date = c("*-05-01", "*-06-30")` for observations from May and
    June of any year.
-   `auk_last_edited()`: filter to checklists from a range of last
    edited dates, useful for extracting just new or recently edited
    data.
-   `auk_protocol()`: filter to checklists that following a specific
    search protocol, either stationary, traveling, or casual.
-   `auk_project()`: filter to checklists collected as part of a
    specific project (e.g. a breeding bird survey).
-   `auk_time()`: filter to checklists started during a range of
    times-of-day.
-   `auk_duration()`: filter to checklists with observation durations
    within a given range.
-   `auk_distance()`: filter to checklists with distances travelled
    within a given range.
-   `auk_breeding()`: only retain observations that have an associate
    breeding bird atlas code.
-   `auk_complete()`: only retain checklists in which the observer has
    specified that they recorded all species seen or heard. It is
    necessary to retain only complete records for the creation of
    presence-absence data, because the “absence”” information is
    inferred by the lack of reporting of a species on checklists.

Note that all of the functions listed above only modify the `auk_ebd`
object, in order to define the filters. Once the filters have been
defined, the filtering is actually conducted using `auk_filter()`.

    # sample data
    f <- system.file("extdata/ebd-sample.txt", package = "auk")
    # define an EBD reference and a set of filters
    ebd <- auk_ebd(f) %>% 
      # species: common and scientific names can be mixed
      auk_species(species = c("Canada Jay", "Cyanocitta cristata")) %>%
      # country: codes and names can be mixed; case insensitive
      auk_country(country = c("US", "Canada", "mexico")) %>%
      # bbox: long and lat in decimal degrees
      # formatted as `c(lng_min, lat_min, lng_max, lat_max)`
      auk_bbox(bbox = c(-100, 37, -80, 52)) %>%
      # date: use standard ISO date format `"YYYY-MM-DD"`
      auk_date(date = c("2012-01-01", "2012-12-31")) %>%
      # time: 24h format
      auk_time(start_time = c("06:00", "09:00")) %>%
      # duration: length in minutes of checklists
      auk_duration(duration = c(0, 60)) %>%
      # complete: all species seen or heard are recorded
      auk_complete()
    ebd
    #> Input 
    #>   EBD: /Users/mes335/projects/auk/inst/extdata/ebd-sample.txt 
    #> 
    #> Output 
    #>   Filters not executed
    #> 
    #> Filters 
    #>   Species: Cyanocitta cristata, Perisoreus canadensis
    #>   Countries: CA, MX, US
    #>   States: all
    #>   Counties: all
    #>   BCRs: all
    #>   Bounding box: Lon -100 - -80; Lat 37 - 52
    #>   Years: all
    #>   Date: 2012-01-01 - 2012-12-31
    #>   Start time: 06:00-09:00
    #>   Last edited date: all
    #>   Protocol: all
    #>   Project code: all
    #>   Duration: 0-60 minutes
    #>   Distance travelled: all
    #>   Records with breeding codes only: no
    #>   Complete checklists only: yes

In all cases, extensive checks are performed to ensure filters are
valid. For example, species are checked against the official [eBird
taxonomy](http://help.ebird.org/customer/portal/articles/1006825-the-ebird-taxonomy)
and countries are checked using the
[`countrycode`](https://github.com/vincentarelbundock/countrycode)
package.

Each of the functions described in the *Defining filters* section only
defines a filter. Once all of the required filters have been set,
`auk_filter()` should be used to compile them into an AWK script and
execute it to produce an output file. So, as an example of bringing all
of these steps together, the following commands will extract all Canada
Jay and Blue Jay records from Canada and save the results to a
tab-separated text file for subsequent use:

    output_file <- "ebd_filtered_blja-grja.txt"
    ebd_filtered <- system.file("extdata/ebd-sample.txt", package = "auk") %>% 
      auk_ebd() %>% 
      auk_species(species = c("Canada Jay", "Cyanocitta cristata")) %>% 
      auk_country(country = "Canada") %>% 
      auk_filter(file = output_file)

**Filtering the full dataset typically takes at least a couple hours**,
so set it running then go grab lunch!

### Reading

eBird Basic Dataset files can be read with `read_ebd()`:

    system.file("extdata/ebd-sample.txt", package = "auk") %>% 
      read_ebd() %>% 
      str()
    #> tibble [494 × 45] (S3: tbl_df/tbl/data.frame)
    #>  $ checklist_id             : chr [1:494] "S6852862" "S14432467" "S39033556" "S38303088" ...
    #>  $ global_unique_identifier : chr [1:494] "URN:CornellLabOfOrnithology:EBIRD:OBS97935965" "URN:CornellLabOfOrnithology:EBIRD:OBS201605886" "URN:CornellLabOfOrnithology:EBIRD:OBS530638734" "URN:CornellLabOfOrnithology:EBIRD:OBS520887169" ...
    #>  $ last_edited_date         : chr [1:494] "2016-02-22 14:59:49" "2013-06-16 17:34:19" "2017-09-06 13:13:34" "2017-07-24 15:17:16" ...
    #>  $ taxonomic_order          : num [1:494] 20145 20145 20145 20145 20145 ...
    #>  $ category                 : chr [1:494] "species" "species" "species" "species" ...
    #>  $ common_name              : chr [1:494] "Green Jay" "Green Jay" "Green Jay" "Green Jay" ...
    #>  $ scientific_name          : chr [1:494] "Cyanocorax yncas" "Cyanocorax yncas" "Cyanocorax yncas" "Cyanocorax yncas" ...
    #>  $ observation_count        : chr [1:494] "4" "2" "1" "1" ...
    #>  $ breeding_code            : chr [1:494] NA NA NA NA ...
    #>  $ breeding_category        : chr [1:494] NA NA NA NA ...
    #>  $ age_sex                  : chr [1:494] NA NA NA NA ...
    #>  $ country                  : chr [1:494] "Mexico" "Mexico" "Mexico" "Mexico" ...
    #>  $ country_code             : chr [1:494] "MX" "MX" "MX" "MX" ...
    #>  $ state                    : chr [1:494] "Yucatan" "Chiapas" "Chiapas" "Chiapas" ...
    #>  $ state_code               : chr [1:494] "MX-YUC" "MX-CHP" "MX-CHP" "MX-CHP" ...
    #>  $ county                   : chr [1:494] NA NA NA NA ...
    #>  $ county_code              : chr [1:494] NA NA NA NA ...
    #>  $ iba_code                 : chr [1:494] NA NA NA NA ...
    #>  $ bcr_code                 : int [1:494] 56 60 60 60 60 55 55 60 56 55 ...
    #>  $ usfws_code               : chr [1:494] NA NA NA NA ...
    #>  $ atlas_block              : chr [1:494] NA NA NA NA ...
    #>  $ locality                 : chr [1:494] "Yuc. Hacienda Chichen" "Berlin2_Punto_06" "07_020_LaConcordia_SanFrancsco_Magallanes_P01" "07_020_CerroBola_BuenaVista_tr3_trad_P03" ...
    #>  $ locality_id              : chr [1:494] "L989845" "L2224225" "L6247542" "L6120049" ...
    #>  $ locality_type            : chr [1:494] "P" "P" "P" "P" ...
    #>  $ latitude                 : num [1:494] 20.7 15.8 15.8 15.8 15.7 ...
    #>  $ longitude                : num [1:494] -88.6 -93 -93 -92.9 -92.9 ...
    #>  $ observation_date         : Date[1:494], format: "2010-09-05" "2011-08-18" "2012-02-02" ...
    #>  $ time_observations_started: chr [1:494] "06:30:00" "08:00:00" "09:13:00" "06:40:00" ...
    #>  $ observer_id              : chr [1:494] "obsr55719" "obsr313215" "obsr313215" "obsr313215" ...
    #>  $ sampling_event_identifier: chr [1:494] "S6852862" "S14432467" "S39033556" "S38303088" ...
    #>  $ protocol_type            : chr [1:494] "Traveling" "Traveling" "Stationary" "Stationary" ...
    #>  $ protocol_code            : chr [1:494] "P22" "P22" "P21" "P21" ...
    #>  $ project_code             : chr [1:494] "EBIRD" "EBIRD_MEX" "EBIRD_MEX" "EBIRD" ...
    #>  $ duration_minutes         : int [1:494] 90 10 10 10 10 120 30 10 80 30 ...
    #>  $ effort_distance_km       : num [1:494] 1 0.257 NA NA 0.257 ...
    #>  $ effort_area_ha           : num [1:494] NA NA NA NA NA NA NA NA NA NA ...
    #>  $ number_observers         : int [1:494] 3 1 1 1 1 2 2 1 13 2 ...
    #>  $ all_species_reported     : logi [1:494] TRUE TRUE TRUE TRUE TRUE TRUE ...
    #>  $ group_identifier         : chr [1:494] NA NA NA NA ...
    #>  $ has_media                : logi [1:494] FALSE FALSE FALSE FALSE FALSE FALSE ...
    #>  $ approved                 : logi [1:494] TRUE TRUE TRUE TRUE TRUE TRUE ...
    #>  $ reviewed                 : logi [1:494] FALSE FALSE FALSE FALSE FALSE FALSE ...
    #>  $ reason                   : chr [1:494] NA NA NA NA ...
    #>  $ trip_comments            : chr [1:494] NA "Alonso Gomez Hdz Monitoreo Comunitario, Transectos en Bosque de Pino Encino, La Concordia,1098 msnm" "Miguel Mndez Lpez" "Rogelio Lpez Encino" ...
    #>  $ species_comments         : chr [1:494] NA NA NA NA ...
    #>  - attr(*, "rollup")= logi TRUE

## Presence-absence data

For many applications, presence-only data are sufficient; however, for
modeling and analysis, presence-absence data are required. `auk`
includes functionality to produce presence-absence data from eBird
checklists. For full details, consult the vignette: `vignette("auk")`.

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](CONDUCT.md). By participating in this project you agree to
abide by its terms.

## Acknowledgements

This package is based on AWK scripts provided as part of the eBird Data
Workshop given by Wesley Hochachka, Daniel Fink, Tom Auer, and Frank La
Sorte at the 2016 NAOC on August 15, 2016.

`auk` benefited significantly from the [rOpenSci](https://ropensci.org/)
review process, including helpful suggestions from [Auriel
Fournier](https://github.com/aurielfournier) and [Edmund
Hart](https://github.com/emhart).

## References

    eBird Basic Dataset. Version: ebd_relFeb-2018. Cornell Lab of Ornithology, Ithaca, New York. May 2013.

[![](http://www.ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# auk 0.5.1

- drop `data.table` dependency, no longer needed with `readr` speed improvements
- fix bug arising from 'breeding bird atlas code' being renamed to 'breeding code' (issue #58)

# auk 0.5.0

- update to align with 2021 eBird taxonomy

# auk 0.4.4

- updates to align with readr 2.0

# auk 0.4.3

- `get_ebird_taxonomy()` now fails gracefully when eBird API is not accessible, fixing the CRAN check errors https://cran.r-project.org/web/checks/check_results_auk.html

# auk 0.4.2

- new `auk_county()` filter
- new `auk_year()` filter
- Drop taxonomy warnings since there was no taxonomy update this year

# auk 0.4.1

- Family common names now included in eBird taxonomy
- `auk_select()` now requires certain columns to be kept
- Better handling of file paths with `prefix` argument in `auk_split()`
- Fixed bug causing undescribed species to be dropped by `auk_rollup()`
- Add a `ll_digits` argument to `filter_repeat_visits()` to round lat/lng prior to identifying sites
- Change of default parameters to `filter_repeat_visits()`
- `auk_bbox()` now takes sf/raster spatial objects and grabs bbox from them

# auk 0.4.0

- Updated to 2019 eBird taxonomy
- `auk_observer()` filter added
- `tidyr::complete_()` deprecated, stopped using

# auk 0.3.3

- Dates can now wrap in `auk_date()`, e.g. use `date = c("*-12-01", "*-01-31")` for records from December or January
- Fixed bug preventing dropping of `age/sex` column
- Allow for a wider variety of protocols in `auk_protocol()`
- Addresing some deprecated functions from rlang
- Fixed bug causing `auk_set_awk_path()` to fail

# auk 0.3.2

- Work around for bug in system2() in some R versions: https://bugs.r-project.org/bugzilla/show_bug.cgi?id=17508
- Adding a filter for PROALAS checklists to `auk_protocol()`

# auk 0.3.1

- `rlang::UQ()` and `rlang::UQS()` deprecated, switching to `!!` and `!!!`
- `auk_unique()` now keeps track of all sampling event and observer IDs that comprise a group checklist

# auk 0.3.0

- Updated to 2018 taxonomy; new function `get_ebird_taxonomy()` to get taxonomy via the eBird API
- Better handling of taxonomy versions, many functions now take a `taxonomy_version` argument and use the eBird API to get the taxonomy
- `auk_getpath()` renamed `auk_get_awk_path()`, and added `auk_set_awk_path()`
- Added `auk_set_ebd_path()` and `auk_get_ebd_path()` to set and get the 
`EBD_PATH` environment variable. Now users only need to set this once and just 
refer to the file name, rather than specifying the full path every time.
- Functions to prepare data for occupancy modeling: `filter_repeat_visits()` and `format_unmarked_occu()`
- New `auk_bcr()` function to extract data from BCRs
- Added `bcr_codes` data frame to look up BCR names and codes
- "Area" protocol added to `auk_protocol()` filter.
- `auk_extent()` renamed `auk_bbox()`; `auk_extent()` deprecated and redirects to `auk_bbox()`
- `auk_zerofill()` now checks for complete checklists and gives option to not rollup
- `auk_rollup()` now gives the option of keeping higher taxa via `drop_higher` argument
- `auk_clean()` deprecated
- Fixed package load error when `EBD_PATH` is invalid
- Fixed bug when reading files with a blank column using `readr`

# auk 0.2.2

- Updated to work with EDB version 1.9
- Modified tests to be more general to all sample data
- `ebird_species()` now returns 6-letter species codes
- Fixed bug causing auk to fail on files downloaded via custom download form
- Fixed bug with `normalizePath()` use on Windows
- Fixed bug with `system2()` on Windows

# auk 0.2.1

- Patch release fixing a couple bugs
- Removed all non-ASCII characters from example files, closes [issue #14](https://github.com/CornellLabofOrnithology/auk/issues/14)
- Fixed issue with state filtering not working, closes [issue $16](https://github.com/CornellLabofOrnithology/auk/issues/16)

# auk 0.2.0

- New function, `auk_split()`, splits EBD up into multiple files by species
- New object, `auk_sampling`, and associated methods for working with the sampling data only
- New function, `auk_select()`, for selecting a subset of columns
- `auk_date()` now allows filtering date ranges across years using wildcards, e.g. `date = c("*-05-01", "*-06-30")` for observations from May and June of any year
- New function, `auk_state()` for filtering by state
- Now using AWK arrays to speed up country and species filtering; ~20% speed up when filtering on many species/countries
- Allow selection of a subset of columns when filtering
- Remove free text columns in `auk_clean()` to decrease file size
- Updated to work with Feb 2018 version of EBD
- Fixed broken dependency on `countrycode` package

# auk 0.1.0

- eBird taxonomy update to August 2017 version, users should download the most recent EBD to ensure the taxonomy is in sync with the new package
- Manually set AWK path with environment variable `AWK_PATH` in `.Renviron` file 
- `auk_distance`, `auk_breeding`, `auk_protocol`, and `auk_project` filters added
- Users can now specify a subset of columns to return when calling auk_filter using the keep and drop arguments
- Many changes suggested by rOpenSci package peer review process, see https://github.com/ropensci/onboarding/issues/136 for details
- New vignette added to aid those wanting to contribute to package development

# auk 0.0.2

- Patch release converting ebird_taxonomy to ASCII to pass CRAN checks

# auk 0.0.1

- First CRAN release# CONTRIBUTING

## Please contribute!

We love collaboration.

## Bugs?

- Submit an issue on the Issues page [here](https://github.com/CornellLabofOrnithology/auk/issues)

## Code contributions

- Fork this repo to your Github account
- Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/auk.git`
- Make sure to track progress upstream (i.e., on our version of `auk` at `CornellLabofOrnithology/auk`) by doing `git remote add upstream https://github.com/CornellLabofOrnithology/auk.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
- Make your changes (bonus points for making changes on a new branch)
- If you alter package functionality at all (e.g., the code itself, not just documentation) please do write some tests to cove the new functionality.
- Push up to your account
- Submit a pull request to home base at `CornellLabofOrnithology/auk`

### Thanks for contributing!# auk 0.5.1

- drop `data.table` dependency, no longer needed with `readr` speed improvements
- fix bug arising from 'breeding bird atlas code' being renamed to 'breeding code' (issue #58)

# Test environments

- local OS X install, R 4.1
- OS X (github actions), R 4.1
- Windows (github actions), R 4.1
- ubuntu 14.04 (github actions), R 4.1
- win-builder (devel and release)

# R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTEs
---
output: md_document
editor_options: 
  chunk_output_type: console
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# auk: eBird Data Extraction and Processing in R  <img src="logo.png" align="right" width=140/>

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/auk)](https://cran.r-project.org/package=auk)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/auk?color=brightgreen)](http://www.r-pkg.org/pkg/auk)
[![R-CMD-check](https://github.com/CornellLabofOrnithology/auk/workflows/R-CMD-check/badge.svg)](https://github.com/CornellLabofOrnithology/auk/actions)
[![rOpenSci](https://badges.ropensci.org/136_status.svg)](https://github.com/ropensci/onboarding/issues/136)
<!-- badges: end -->

## Overview

[eBird](http://www.ebird.org) is an online tool for recording bird observations. Since its inception, over 600 million records of bird sightings (i.e. combinations of location, date, time, and bird species) have been collected, making eBird one of the largest citizen science projects in history and an extremely valuable resource for bird research and conservation. The full eBird database is packaged as a text file and available for download as the [eBird Basic Dataset (EBD)](http://ebird.org/ebird/data/download). Due to the large size of this dataset, it must be filtered to a smaller subset of desired observations before reading into R. This filtering is most efficiently done using AWK, a Unix utility and programming language for processing column formatted text data. This package acts as a front end for AWK, allowing users to filter eBird data before import into R.

For a comprehensive resource on using eBird data for modeling species distributions, consult the free online book [Best Practices for Using eBird Data](https://cornelllabofornithology.github.io/ebird-best-practices/) and the association paper _Analytical guidelines to increase the value of community science data: An example using eBird data to estimate species distributions_ ([Johnston et al. 2021](https://onlinelibrary.wiley.com/doi/10.1111/ddi.13271)).

## Installation

```{r gh-install, eval=FALSE}
# cran release
install.packages("auk")

# or install the development version from github
# install.packages("remotes")
remotes::install_github("CornellLabofOrnithology/auk")
```

`auk` requires the Unix utility AWK, which is available on most Linux and Mac OS X machines. Windows users will first need to install [Cygwin](https://www.cygwin.com) before using this package. Note that **Cygwin must be installed in the default location** (`C:/cygwin/bin/gawk.exe` or `C:/cygwin64/bin/gawk.exe`) in order for `auk` to work.

## Vignette

Full details on using `auk` to produce both presence-only and presence-absence data are outlined in the [vignette](https://cornelllabofornithology.github.io/auk/articles/auk.html).

## Cheatsheet

An `auk` cheatsheet was developed by [Mickayla Johnston](https://www.linkedin.com/in/mickayla-johnston/):

<a href="https://github.com/CornellLabofOrnithology/auk/blob/master/cheatsheet/auk-cheatsheet.pdf"><img src="cheatsheet/auk-cheatsheet.png" width=400/></a>



## `auk` and `rebird`

Those interested in eBird data may also want to consider [`rebird`](https://github.com/ropensci/rebird), an R package that provides an interface to the [eBird APIs](https://confluence.cornell.edu/display/CLOISAPI/eBirdAPIs). The functions in `rebird` are mostly limited to accessing recent (i.e. within the last 30 days) observations, although `ebirdfreq()` does provide historical frequency of observation data. In contrast, `auk` gives access to the full set of ~ 500 million eBird observations. For most ecological applications, users will require `auk`; however, for some use cases, e.g. building tools for birders, `rebird` provides a quick and easy way to access data.

## A note on versions

This package contains a current (as of the time of package release) version of the [bird taxonomy used by eBird](http://help.ebird.org/customer/portal/articles/1006825-the-ebird-taxonomy). This taxonomy determines the species that can be reported in eBird and therefore the species that users of `auk` can extract. eBird releases an updated taxonomy once a year, typically in August, at which time `auk` will be updated to include the current taxonomy. When using `auk`, users should be careful to ensure that the version they're using is in sync with the eBird Basic Dataset they're working with. This is most easily accomplished by always using the must recent version of `auk` and the most recent release of the dataset.

## Quick start

This package uses the command-line program AWK to extract subsets of the eBird Basic Dataset for use in R. This is a multi-step process:

1. Define a reference to the eBird data file.
2. Define a set of spatial, temporal, or taxonomic filters. Each type of filter corresponds to a different function, e.g. `auk_species` to filter by species. At this stage the filters are only set up, no actual filtering is done until the next step.
3. Filter the eBird data text file, producing a new text file with only the selected rows.
4. Import this text file into R as a data frame.

Because the eBird dataset is so large, step 3 typically takes several hours to run. Here's a simple example that extract all Canada Jay records from within Canada.

```{r packages, include=FALSE}
library(auk)
```


```{r quickstart, eval = FALSE}
library(auk)
# path to the ebird data file, here a sample included in the package
# get the path to the example data included in the package
# in practice, provide path to ebd, e.g. f_in <- "data/ebd_relFeb-2018.txt
f_in <- system.file("extdata/ebd-sample.txt", package = "auk")
# output text file
f_out <- "ebd_filtered_grja.txt"
ebird_data <- f_in %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_species(species = "Canada Jay") %>% 
  auk_country(country = "Canada") %>% 
  # 3. run filtering
  auk_filter(file = f_out) %>% 
  # 4. read text file into r data frame
  read_ebd()
```

For those not familiar with the pipe operator (`%>%`), the above code could be rewritten:

```{r quickstart-nopipes, eval = FALSE}
f_in <- system.file("extdata/ebd-sample.txt", package = "auk")
f_out <- "ebd_filtered_grja.txt"
ebd <- auk_ebd(f_in)
ebd_filters <- auk_species(ebd, species = "Canada Jay")
ebd_filters <- auk_country(ebd_filters, country = "Canada")
ebd_filtered <- auk_filter(ebd_filters, file = f_out)
ebd_df <- read_ebd(ebd_filtered)
```

## Usage

### Filtering

`auk` uses a [pipeline-based workflow](http://r4ds.had.co.nz/pipes.html) for defining filters, which can then be compiled into an AWK script. Users should start by defining a reference to the dataset file with `auk_ebd()`. Then any of the following filters can be applied:

- `auk_species()`: filter by species using common or scientific names.
- `auk_country()`: filter by country using the standard English names or [ISO 2-letter country codes](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2).
- `auk_state()`: filter by state using eBird state codes, see `?ebird_states`.
- `auk_bcr()`: filter by [Bird Conservation Region (BCR)](http://nabci-us.org/resources/bird-conservation-regions/) using BCR codes, see `?bcr_codes`.
- `auk_bbox()`: filter by spatial bounding box, i.e. a range of latitudes and longitudes in decimal degrees.
- `auk_date()`: filter to checklists from a range of dates. To extract observations from a range of dates, regardless of year, use the wildcard "`*`" in place of the year, e.g. `date = c("*-05-01", "*-06-30")` for observations from May and June of any year.
- `auk_last_edited()`: filter to checklists from a range of last edited dates, useful for extracting just new or recently edited data.
- `auk_protocol()`: filter to checklists that following a specific search protocol, either stationary, traveling, or casual.
- `auk_project()`: filter to checklists collected as part of a specific project (e.g. a breeding bird survey).
- `auk_time()`: filter to checklists started during a range of times-of-day.
- `auk_duration()`: filter to checklists with observation durations within a given range.
- `auk_distance()`: filter to checklists with distances travelled within a given range.
- `auk_breeding()`: only retain observations that have an associate breeding bird atlas code.
- `auk_complete()`: only retain checklists in which the observer has specified that they recorded all species seen or heard. It is necessary to retain only complete records for the creation of presence-absence data, because the "absence"" information is inferred by the lack of reporting of a species on checklists. 

Note that all of the functions listed above only modify the `auk_ebd` object, in order to define the filters. Once the filters have been defined, the filtering is actually conducted using `auk_filter()`.

```{r auk-filter}
# sample data
f <- system.file("extdata/ebd-sample.txt", package = "auk")
# define an EBD reference and a set of filters
ebd <- auk_ebd(f) %>% 
  # species: common and scientific names can be mixed
  auk_species(species = c("Canada Jay", "Cyanocitta cristata")) %>%
  # country: codes and names can be mixed; case insensitive
  auk_country(country = c("US", "Canada", "mexico")) %>%
  # bbox: long and lat in decimal degrees
  # formatted as `c(lng_min, lat_min, lng_max, lat_max)`
  auk_bbox(bbox = c(-100, 37, -80, 52)) %>%
  # date: use standard ISO date format `"YYYY-MM-DD"`
  auk_date(date = c("2012-01-01", "2012-12-31")) %>%
  # time: 24h format
  auk_time(start_time = c("06:00", "09:00")) %>%
  # duration: length in minutes of checklists
  auk_duration(duration = c(0, 60)) %>%
  # complete: all species seen or heard are recorded
  auk_complete()
ebd
```

In all cases, extensive checks are performed to ensure filters are valid. For example, species are checked against the official [eBird taxonomy](http://help.ebird.org/customer/portal/articles/1006825-the-ebird-taxonomy) and countries are checked using the [`countrycode`](https://github.com/vincentarelbundock/countrycode) package.

Each of the functions described in the *Defining filters* section only defines a filter. Once all of the required filters have been set, `auk_filter()` should be used to compile them into an AWK script and execute it to produce an output file. So, as an example of bringing all of these steps together, the following commands will extract all Canada Jay and Blue Jay records from Canada and save the results to a tab-separated text file for subsequent use:

```{r auk-complete, eval = FALSE}
output_file <- "ebd_filtered_blja-grja.txt"
ebd_filtered <- system.file("extdata/ebd-sample.txt", package = "auk") %>% 
  auk_ebd() %>% 
  auk_species(species = c("Canada Jay", "Cyanocitta cristata")) %>% 
  auk_country(country = "Canada") %>% 
  auk_filter(file = output_file)
```

**Filtering the full dataset typically takes at least a couple hours**, so set it running then go grab lunch!

### Reading

eBird Basic Dataset files can be read with `read_ebd()`:

```{r read}
system.file("extdata/ebd-sample.txt", package = "auk") %>% 
  read_ebd() %>% 
  str()
```

## Presence-absence data

For many applications, presence-only data are sufficient; however, for modeling and analysis, presence-absence data are required. `auk` includes functionality to produce presence-absence data from eBird checklists. For full details, consult the vignette: `vignette("auk")`.

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

## Acknowledgements

This package is based on AWK scripts provided as part of the eBird Data Workshop given by Wesley Hochachka, Daniel Fink, Tom Auer, and Frank La Sorte at the 2016 NAOC on August 15, 2016.

`auk` benefited significantly from the [rOpenSci](https://ropensci.org/) review process, including helpful suggestions from 
[Auriel Fournier](https://github.com/aurielfournier) and [Edmund Hart](https://github.com/emhart).

## References

```
eBird Basic Dataset. Version: ebd_relFeb-2018. Cornell Lab of Ornithology, Ithaca, New York. May 2013.
```

[![](http://www.ropensci.org/public_images/github_footer.png)](http://ropensci.org)---
title: "Introduction to auk"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to auk}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE, error = FALSE, message = FALSE
)
suppressPackageStartupMessages(library(auk))
suppressPackageStartupMessages(library(dplyr))
```

[eBird](http://www.ebird.org) is an online tool for recording bird observations. Since its inception, nearly 500 million records of bird sightings (i.e. combinations of location, date, time, and bird species) have been collected, making eBird one of the largest citizen science projects in history and an extremely valuable resource for bird research and conservation. The full eBird database is packaged as a text file and available for download as the [eBird Basic Dataset (EBD)](http://ebird.org/ebird/data/download). Due to the large size of this dataset, it must be filtered to a smaller subset of desired observations before reading into R. This filtering is most efficiently done using AWK, a Unix utility and programming language for processing column formatted text data. This package acts as a front end for AWK, allowing users to filter eBird data before import into R.

This vignette is divided into three sections. The first section provides background on the eBird data and motivation for the development of this package. The second section outlines the use of `auk` for filtering text file to produce a presence-only dataset. The final section demonstrates how `auk` can be used to produce zero-filled, presence-absence (or more correctly presence–non-detection) data, a necessity for many modeling and analysis applications.

## Quick start

This package uses the command-line program AWK to extract subsets of the eBird Basic Dataset for use in R. This is a multi-step process:

1. Define a reference to the eBird data file.
2. Define a set of spatial, temporal, or taxonomic filters. Each type of filter corresponds to a different function, e.g. `auk_species` to filter by species. At this stage the filters are only set up, no actual filtering is done until the next step.
3. Filter the eBird data text file, producing a new text file with only the selected rows.
4. Import this text file into R as a data frame.

Because the eBird dataset is so large, step 3 typically takes several hours to run. Here's a simple example that extract all Canada Jay records from within Canada.

```{r quickstart, eval = FALSE}
library(auk)
# path to the ebird data file, here a sample included in the package
# in practice, provide path to ebd, e.g. input_file <- "data/ebd_relFeb-2018.txt"
input_file <- system.file("extdata/ebd-sample.txt", package = "auk")
# output text file
output_file <- "ebd_filtered_grja.txt"
ebird_data <- input_file %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_species(species = "Canada Jay") %>% 
  auk_country(country = "Canada") %>% 
  # 3. run filtering
  auk_filter(file = output_file) %>% 
  # 4. read text file into r data frame
  read_ebd()
```

For those not familiar with the pipe operator (`%>%`), the above code could be rewritten:

```{r quickstart-nopipes, eval = FALSE}
input_file <- system.file("extdata/ebd-sample.txt", package = "auk")
output_file <- "ebd_filtered_grja.txt"
ebd <- auk_ebd(input_file)
ebd_filters <- auk_species(ebd, species = "Canada Jay")
ebd_filters <- auk_country(ebd_filters, country = "Canada")
ebd_filtered <- auk_filter(ebd_filters, file = output_file)
ebd_df <- read_ebd(ebd_filtered)
```

## Background

### The eBird Basic Dataset

The eBird database currently contains nearly 500 million bird observations, and this rate of increase is accelerating as new users join eBird. These data are an extremely valuable tool both for basic science and conservation; however, given the sheer amount of data, accessing eBird data poses a unique challenge. Currently, access to the complete set of eBird observations is provided via the eBird Basic Dataset (EBD). This is a tab-separated text file, released quarterly, containing all validated bird sightings in the eBird database at the time of release. Each row corresponds to the sighting of a single species within a checklist and, in addition to the species and number of individuals reported, information is provided at the checklist level (location, time, date, search effort, etc.).

In addition, eBird provides a Sampling Event Data file that contains the checklist-level data for every valid checklist submitted to eBird, including checklists for which no species of birds were reported. In this file, each row corresponds to a checklist and only the checklist-level variables are included, not the associated bird data. While the eBird Basic Dataset provides presence-only data, it can be combined with the Sampling Event Data file to produce presence-absence data. This process is described below.

For full metadata on the both datasets, consult the documentation provided when the [files are downloaded](http://ebird.org/ebird/data/download).

## `auk` vs. `rebird`

Those interested in eBird data may also want to consider [`rebird`](https://docs.ropensci.org/rebird/), an R package that provides an interface to the [eBird APIs](https://confluence.cornell.edu/display/CLOISAPI/eBirdAPIs). The functions in `rebird` are mostly limited to accessing recent (i.e. within the last 30 days) observations, although `ebirdfreq()` does provide historical frequency of observation data. In contrast, `auk` gives access to the full set of ~ 500 million eBird observations. For most ecological applications, users will require `auk`; however, for some use cases, e.g. building tools for birders, `rebird` provides a quick and easy way to access data.

### Data access

To access eBird data, begin by [creating an eBird account and signing in](https://secure.birds.cornell.edu/cassso/login). Then visit the [Download Data](http://ebird.org/ebird/data/download) page. eBird data access is free; however, you will need to [request access](http://ebird.org/ebird/data/request) in order to obtain access to the EBD. Filling out the access request form allows eBird to keep track of the number of people using the data and obtain information on the applications for which the data are used

Once you have access to the data, proceed to the [download page](http://ebird.org/ebird/data/download/ebd). There are two download options: prepackage download and custom download. Downloading the prepackaged option gives you access to the full global dataset. If you choose this route, you'll likely want to download both the EBD (~ 25 GB) and corresponding Sampling Event Data (~ 2.5 GB). If you know you're likely to only need data for a single species, or a small region, you can request a custom download be prepared consisting of only a subset of the data. This will result in significantly smaller files; however, note that custom requests that would result in huge numbers of checklists (e.g. all records from the US) won't work. In either case, download and decompress the files.

### Example data

This package comes with two example datasets. The first is suitable for practicing filtering the EBD and producing presence-only data. It's a sample of 500 records from the EBD. It contains data from North and Central America from 2010-2012 on 4 jay species: Canada Jay, Blue Jay, Steller's Jay, and Green Jay. It can be accessed with:

```{r example-data-1, eval = FALSE}
library(auk)
library(dplyr)
system.file("extdata/ebd-sample.txt", package = "auk")
```

The second is suitable for producing zero-filled, presence-absence data. It contains every sighting from Singapore in 2012 of Collared Kingfisher, White-throated Kingfisher, and Blue-eared Kingfisher. The full Sampling Event Data file is also included, and contains all checklists from Singapore in 2012. These files can be accessed with:

```{r example-data-2, eval = FALSE}
# ebd
system.file("extdata/zerofill-ex_ebd.txt", package = "auk")
# sampling event data
system.file("extdata/zerofill-ex_sampling.txt", package = "auk")
```

**Important note:** in this vignette, `system.file()` is used to return the path to the example data included in this package. When using `auk` in practice, provide the path to the location of the EBD on your computer, this could be a relative path, e.g. `"data/ebd_relFeb-2018.txt"`, or an absolute path, e.g. `"~/ebird/ebd_relFeb-2018/ebd_relFeb-2018.txt"`.

### AWK

R typically works with objects in memory and, as a result, there is a hard limit on the size of objects that can be brought into R. Because eBird contains nearly 500 million sightings, the eBird Basic Dataset is an inherently large file (~150 GB uncompressed) and therefore impossible to manipulate directly in R. Thus it is generally necessary to create a subset of the file outside of R, then import this smaller subset for analysis.

AWK is a Unix utility and programming language for processing column formatted text data. It is highly flexible and extremely fast, making it a valuable tool for pre-processing the eBird data in order to create the smaller subset of data that is required. Users of the data can use AWK to produce a smaller file, subsetting the full text file taxonomically, spatially, or temporally, in order to produce a smaller file that can then be loaded in to R for visualization, analysis, and modelling. 

Although AWK is a powerful tool, it has three disadvantages: it requires learning the syntax of a new language, it is only accessible via the command line, and it results in a portion of your workflow existing outside of R. This package is a wrapper for AWK specifically designed for filtering eBird data The goal is to ease the use of the this data by removing the hurdle of learning and using AWK.

Linux and Mac users should already have AWK installed on their machines, however, Windows uses will need to install [Cygwin](https://www.cygwin.com) to gain access to AWK. Note that **Cygwin should be installed in the default location** (`C:/cygwin/bin/gawk.exe` or `C:/cygwin64/bin/gawk.exe`) in order for `auk` to work. To check that AWK is installed and can be found run `auk_getpath()`.

If AWK is installed in a non-standard location, or can't be found by `auk`, you can manually set the path to AWK. To do so, set the `AWK_PATH` environment in your `.Renviron` file. For example, Mac and Linux users might add the following line:

```
AWK_PATH=/usr/bin/awk
```

while Windows users might add:

```
AWK_PATH=C:/cygwin64/bin/gawk.exe
```

### A note on versions

This package contains a current (as of the time of package release) version of the [bird taxonomy used by eBird](https://ebird.org/science/use-ebird-data/the-ebird-taxonomy). This taxonomy determines the species that can be reported in eBird and therefore the species that users of `auk` can extract from the EBD. eBird releases an updated taxonomy once a year, typically in August, at which time `auk` will be updated to include the current taxonomy. When using `auk`, users should be careful to ensure that the version they're using is in sync with the EBD file they're working with. This is most easily accomplished by always using the most recent version of `auk` and the most recent release of the eBird Basic Dataset

## Presence data

The most common use of the eBird data is to produce a set of bird sightings, i.e. where and when was a given species seen. For example, this type of data could be used to produce a map of sighting locations, or to determine if a given bird has been seen in an area of interest. For more analytic work, such as species distribution modeling, presence and absence data are likely preferred (see Guillera-Arroita et al. 2015). Producing presence-absence data will be covered in the next section.

### The `auk_ebd` object

This package uses an `auk_ebd` object to keep track of the input data file, any filters defined, and the output file that is produced after filtering has been executed. By keeping everything wrapped up in one object, the user can keep track of exactly what set of input data and filters produced any given output data. To set up the initial `auk_ebd` object, use `auk_ebd()`:

```{r auk-ebd}
ebd <- system.file("extdata/ebd-sample.txt", package = "auk") %>% 
  auk_ebd()
ebd
```

### Defining filters

`auk` uses a [pipeline-based workflow](https://r4ds.had.co.nz/pipes.html) for defining filters, which can then be compiled into an AWK script. Any of the following filters can be applied:

- `auk_species()`: filter by species using common or scientific names.
- `auk_country()`: filter by country using the standard English names or [ISO 2-letter country codes](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2).
- `auk_state()`: filter by state using the eBird state codes, see `?ebird_states`.
- `auk_bcr()`: filter by [Bird Conservation Region (BCR)](https://nabci-us.org/resources/bird-conservation-regions/) using BCR codes, see `?bcr_codes`.
- `auk_bbox()`: filter by spatial bounding box, i.e. a range of latitudes and longitudes in decimal degrees.
- `auk_date()`: filter to checklists from a range of dates. To extract observations from a range of dates, regardless of year, use the wildcard "`*`" in place of the year, e.g. `date = c("*-05-01", "*-06-30")` for observations from May and June of any year.
- `auk_last_edited()`: filter to checklists from a range of last edited dates, useful for extracting just new or recently edited data.
- `auk_protocol()`: filter to checklists that following a specific search protocol, either stationary, traveling, or casual.
- `auk_project()`: filter to checklists collected as part of a specific project (e.g. a breeding bird survey).
- `auk_time()`: filter to checklists started during a range of times-of-day.
- `auk_duration()`: filter to checklists with observation durations within a given range.
- `auk_distance()`: filter to checklists with distances travelled within a given range.
- `auk_breeding()`: only retain observations that have an associate breeding bird atlas code.
- `auk_complete()`: only retain checklists in which the observer has specified that they recorded all species seen or heard. It is necessary to retain only complete records for the creation of presence-absence data, because the "absence" information is inferred by the lack of reporting of a species on checklists. 

Note that all of the functions listed above only modify the `auk_ebd` object, in order to define the filters. Once the filters have been defined, the filtering is actually conducted using `auk_filter()`.

```{r auk-filter}
ebd_filters <- ebd %>% 
  # species: common and scientific names can be mixed
  auk_species(species = c("Canada Jay", "Cyanocitta cristata")) %>%
  # country: codes and names can be mixed; case insensitive
  auk_country(country = c("US", "Canada", "mexico")) %>%
  # bbox: long and lat in decimal degrees
  # formatted as `c(lng_min, lat_min, lng_max, lat_max)`
  auk_bbox(bbox = c(-100, 37, -80, 52)) %>%
  # date: use standard ISO date format `"YYYY-MM-DD"`
  auk_date(date = c("2012-01-01", "2012-12-31")) %>%
  # time: 24h format
  auk_time(start_time = c("06:00", "09:00")) %>%
  # duration: length in minutes of checklists
  auk_duration(duration = c(0, 60)) %>%
  # complete: all species seen or heard are recorded
  auk_complete()
ebd_filters
```

In all cases, extensive checks are performed to ensure filters are valid. For example, species are checked against the official [eBird taxonomy](https://ebird.org/science/use-ebird-data/the-ebird-taxonomy) and countries are checked using the [`countrycode`](https://github.com/vincentarelbundock/countrycode) package. This is particularly important because filtering is a time consuming process, so catching errors in advance can avoid several hours of wasted time.

### Executing filters

Each of the functions described in the *Defining filters* section only defines a filter. Once all of the required filters have been set, `auk_filter()` should be used to compile them into an AWK script and execute it to produce an output file. So, as an example of bringing all of these steps together, the following commands will extract all Canada Jay and Blue Jay records from Canada and save the results to a tab-separated text file for subsequent use:

```{r auk-complete, eval = FALSE}
output_file <- "ebd_filtered_blja-grja.txt"
ebd_jays <- system.file("extdata/ebd-sample.txt", package = "auk") %>% 
  auk_ebd() %>% 
  auk_species(species = c("Canada Jay", "Cyanocitta cristata")) %>% 
  auk_country(country = "Canada") %>% 
  auk_filter(file = output_file)
```

**Filtering the full EBD typically takes at least a couple hours**, so set it running then go grab lunch!

### Reading

eBird Basic Dataset files can be read with `read_ebd()`. This is a wrapper around `readr::read_delim()`. `read_ebd()` uses `stringsAsFactors = FALSE`, `quote = ""`, sets column classes, and converts variable names to `snake_case`.

```{r read}
system.file("extdata/ebd-sample.txt", package = "auk") %>% 
  read_ebd() %>% 
  glimpse()
```

`auk_filter()` returns an `auk_ebd` object with the output file paths stored in it. This `auk_ebd` object can then be passed directly to `auk_read()`, allowing for a complete pipeline. For example, we can create an `auk_ebd` object, define filters, filter with AWK, and read in the results all at once.

```{r read-auk-ebd, eval = FALSE}
output_file <- "ebd_filtered_blja-grja.txt"
ebd_df <- system.file("extdata/ebd-sample.txt", package = "auk") %>% 
  auk_ebd() %>% 
  auk_species(species = c("Canada Jay", "Cyanocitta cristata")) %>% 
  auk_country(country = "Canada") %>% 
  auk_filter(file = output_file) %>% 
  read_ebd()
```

### Saving the AWK command

The AWK script can be saved for future reference by providing an output filename to `awk_file`. In addition, by setting `execute = FALSE` the AWK script will be generated but not run.

```{r awk-script}
awk_script <- system.file("extdata/ebd-sample.txt", package = "auk") %>% 
  auk_ebd() %>% 
  auk_species(species = c("Canada Jay", "Cyanocitta cristata")) %>% 
  auk_country(country = "Canada") %>% 
  auk_filter(awk_file = "awk-script.txt", execute = FALSE)
# read back in and prepare for printing
awk_file <- readLines(awk_script)
unlink("awk-script.txt")
awk_file[!grepl("^[[:space:]]*$", awk_file)] %>% 
  paste0(collapse = "\n") %>% 
  cat()
```

### Group checklists

eBird allows observers birding together to share checklists. This process creates a new copy of the original checklist for each observer with whom the original checklist was shared; these copies can then be tweaked to add or remove some species that weren’t seen by the entire group, or altering the sampling-event data. For most applications, it's best to remove these duplicate (or near-duplicate) checklists. `auk_unique()` removes duplicates resulting from group checklists by selecting the observation with the lowest `sampling_event_identifier` (a unique ID for each checklist); this is the original checklists from which shared copies were generated. In addition to removing duplicates, a `checklist_id` field is added, which is equal to the `sampling_event_identifier` for non-group checklists and the `group_identifier` for grouped checklists. After running `auk_unique()`, every species will have a single entry for each `checklist_id`.

`read_ebd()` automatically runs `auk_unique()`, however, we can use `unique = FALSE` then manually run `auk_unique()`.

```{r auk-unique}
# read in an ebd file and don't automatically remove duplicates
ebd_dupes <- system.file("extdata/ebd-sample.txt", package = "auk") %>%
  read_ebd(unique = FALSE)
# remove duplicates
ebd_unique <- auk_unique(ebd_dupes)
# compare number of rows
nrow(ebd_dupes)
nrow(ebd_unique)
```

### Taxonomic rollup

The eBird Basic Dataset includes both true species and other taxa, including domestics, hybrids, subspecies, "spuhs", and recognizable forms. In some cases, a checklist may contain multiple records for the same species, for example, both Audubon's and Myrtle Yellow-rumped Warblers, as well as some records that are not resolvable to species, for example, "warbler sp.". For most use cases, a single record for each species on each checklist is desired. The function `ebd_rollup()` addresses these cases by removing taxa not identifiable to species and rolling up taxa identified below species level to a single record for each species in each checklist. 

```{r auk-rollup}
# read in sample data without rolling up
ebd <- system.file("extdata/ebd-rollup-ex.txt", package = "auk") %>%
  read_ebd(rollup = FALSE)
# apply roll up
ebd_ru <- auk_rollup(ebd)

# all taxa not identifiable to species are dropped
# taxa below species have been rolled up to species
unique(ebd$category)
unique(ebd_ru$category)

# yellow-rump warbler subspecies rollup
# without rollup, there are three observations
ebd %>%
  filter(common_name == "Yellow-rumped Warbler") %>%
  select(checklist_id, category, common_name, subspecies_common_name,
         observation_count)
# with rollup, they have been combined
ebd_ru %>%
  filter(common_name == "Yellow-rumped Warbler") %>%
  select(checklist_id, category, common_name, observation_count)
```

By default, `read_ebd()` calls `ebd_rollup()` when importing an eBird Basic Dataset file. To avoid this, and retain subspecies, use `read_ebd(rollup = FALSE)`.

## Zero-filled, presence-absence data

For many applications, presence-only data are sufficient; however, for modeling and analysis, presence-absence data are required. eBird observers only explicitly collect presence data, but they have the option of flagging their checklist as "complete" meaning that they are reporting all the species they saw or heard, and identified. Therefore, given a list of positive sightings (the basic dataset) and a list of all checklists (the sampling event data) it is possible to infer absences by filling zeros for all species not explicitly reported. This section of the vignette describes functions for producing zero-filled, presence-absence data.

### Filtering

When preparing to create zero-filled data, the eBird Basic Dataset and sampling event data must be filtered to the same set of checklists to ensure consistency. To ensure these two datasets are synced, provide *both* to `auk_ebd`, then filter as described in the previous section. This will ensure that all the filters applied to the ebd (except species) will be applied to the sampling event data so that we'll be working with the same set of checklists. It is critical that `auk_compete()` is called, since complete checklists are a requirement for zero-filling.

For example, the following filters to only include sightings of Collared Kingfisher between 6 and 10am:

```{r ebd-zf}
# to produce zero-filled data, provide an EBD and sampling event data file
f_ebd <- system.file("extdata/zerofill-ex_ebd.txt", package = "auk")
f_smp <- system.file("extdata/zerofill-ex_sampling.txt", package = "auk")
filters <- auk_ebd(f_ebd, file_sampling = f_smp) %>% 
  auk_species("Collared Kingfisher") %>% 
  auk_time(c("06:00", "10:00")) %>% 
  auk_complete()
filters
```

As with presence-only data, call `auk_filter()` to actually run AWK. Output files must be provided for both the EBD and sampling event data.

```{r zf-filter-fake, echo = FALSE}
# needed to allow building vignette on machines without awk
ebd_sed_filtered <- filters
ebd_sed_filtered$output <- "ebd-filtered.txt"
ebd_sed_filtered$output_sampling <- "sampling-filtered.txt"
```

```{r zf-filter, eval = -1}
ebd_sed_filtered <- auk_filter(filters, 
                               file = "ebd-filtered.txt",
                               file_sampling = "sampling-filtered.txt")
ebd_sed_filtered
```

### Reading and zero-filling

The filtered datasets can now be combined into a zero-filled, presence-absence dataset using `auk_zerofill()`.

```{r auk-zf-fake, echo = FALSE}
# needed to allow building vignette on machines without awk
fake_ebd <- read_ebd(f_ebd)
fake_smp <- read_sampling(f_smp)
# filter in R to fake AWK call
fake_ebd <- subset(
  fake_ebd, 
  all_species_reported & 
    scientific_name %in% filters$filters$species & 
    time_observations_started >= filters$filters$time[1] & 
    time_observations_started <= filters$filters$time[2])
fake_smp <- subset(
  fake_smp, 
  all_species_reported & 
    time_observations_started >= filters$filters$time[1] & 
    time_observations_started <= filters$filters$time[2])
ebd_zf <- auk_zerofill(fake_ebd, fake_smp)
```

```{r auk-zf, eval = -1}
ebd_zf <- auk_zerofill(ebd_sed_filtered)
ebd_zf
```

Filenames or data frames of the basic dataset and sampling event data can also be passed to `auk_zerofill()`; see the documentation for these cases. By default, `auk_zerofill()` returns an `auk_zerofill` object consisting of two data frames that can be linked by a common `checklist_id` field: 

- `ebd_zf$sampling_events` contains the checklist information
- `ebd_zf$observations` contains the species counts and a binary presence-absence variable

```{r zf-components}
head(ebd_zf$observations)
glimpse(ebd_zf$sampling_events)
```

This format is efficient for storage because the checklist information isn't duplicated, however, a single flat data frame is often required for analysis. To collapse the two data frames together use `collapse_zerofill()`, or call `auk_zerofill()` with `collapse = TRUE`.

```{r zf-collapse, eval = -1}
ebd_zf_df <- auk_zerofill(ebd_filtered, collapse = TRUE)
ebd_zf_df <- collapse_zerofill(ebd_zf)
class(ebd_zf_df)
ebd_zf_df
```

## Acknowledgements

This package is based on the AWK scripts provided in a presentation given by Wesley Hochachka, Daniel Fink, Tom Auer, and Frank La Sorte at the 2016 NAOC eBird Data Workshop on August 15, 2016.

`auk` benefited significantly from the [rOpenSci](https://ropensci.org/) review process, including helpful suggestions from 
[Auriel Fournier](https://github.com/aurielfournier) and [Edmund Hart](https://github.com/emhart).

## References

```
eBird Basic Dataset. Version: ebd_relFeb-2018. Cornell Lab of Ornithology, Ithaca, New York. May 2013.

Guillera-Arroita, G., J.J. Lahoz-Monfort, J. Elith, A. Gordon, H. Kujala, P.E. Lentini, M.A. McCarthy, R. Tingley, and B.A. Wintle. 2015. Is my species distribution model fit for purpose? Matching data and models to applications. Global Ecology and Biogeography 24:276-292.
```
---
title: "auk development"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{auk development}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

This vignette describes the process of updating and extending `auk`. Three topics are covered: updating `auk` when a new eBird taxonomy is released, extending `auk` to include new filters, and CRAN submission.

## Updating the eBird taxonomy

The species, and other taxa, available for entry into the eBird database is dependent on the [eBird taxonomy](https://ebird.org/science/use-ebird-data/the-ebird-taxonomy). Every August, the eBird team updates this taxonomy to reflect name changes splits, merges, new species, or any other changes. Historical eBird records are then updated accordingly and subsequent EBD files reflect this updated taxonomy. The `auk` package stores a copy of this taxonomy as the data frame `ebird_taxonomy`, and uses it both for filtering by species (`auk_species()`) and for taxonomic roll-up (`auk_rollup()`). Therefore, `auk` must be updated when a new eBird taxonomy is released. This section described how this is done. It is best to do this after the new taxonomy **and** the new EBD have both been released, otherwise the taxonomy and EBD will be out of sync.

When the eBird taxonomy is updated, the new version can be downloaded from the [eBird website](https://ebird.org/science/use-ebird-data/the-ebird-taxonomy). The taxonomy can be downloaded in csv or Excel format, **be sure to download the Excel file** because the csv file has character encoding issues. Copy this file to `data-raw/`. At this point, you should check that this new taxonomy has the same format as the previous file, which will also be in this directory. Ensure that the same columns are present and that they're named the same.

The file `data-raw/ebird-taxonomy.r` prepares the taxonomy as a data frame to be stored in the package. Open this file and edit the `read_xlsx()` call to point to the new file you just downloaded. Run the code, then open the `ebird_taxonomy` data frame to inspect it and make sure there's no glaring issues. One potential error that should be investigated is non-ASCII characters. Some common names have accented characters (e.g. Rüppell's Griffon, Gyps rueppelli), which can cause problems. `ebird-taxonomy.r` converts these characters to their unaccented equivalents (e.g. Ruppell's Griffon). Check that this record, or others with accented characters, has been properly converted.

Next, update `auk_version_date()` (`R/auk-version-date.r`) to reflect the date of the new taxonomy and the new EBD.

Finally, build the package (`devtools::build()`) and run `R CMD check` (`devtools::check()`). If everything looks good, commit to git and push to GitHub.

## Adding new filters

The primary functionality of `auk` is to apply filters to the EBD to extract a subset of records that can be imported into R and further analyzed. Individual filters are defined by a particular function (e.g. `auk_date()` or `auk_country()`) and correspond to subsetting on a particular column (e.g. "OBSERVATION DATE" and "COUNTRY CODE", respectively). Defining a new filter is a fairly complicated process, involving carefully updating many components of the package, and should only be attempted by experienced R programmers. To add a filter called `color`, the following steps are required:

1. Update `auk_ebd()` (in file `R/auk-ebd.r`) to define the column number for the new filter, create a placeholder in the `auk_ebd` object to store the filtering criteria, and update the `auk_ebd` print method for the new filter.
2. Create a new function `auk_color()` (in file `R/auk-color.r`) that defines the new filter. As a starting point, use one of the other filtering functions. For example to filter on a range of numeric values, start with `auk_duration()`, to filter on a logical (true/false) variable use `auk_complete()`, or to filter on a discrete, categorical variable use `auk_country()`. Be sure to apply extensive checking on the validity of inputs and update the documentation, including examples.
3. Update `auk_filter()` (in file `R/auk-filter.r`) to incorporate the filtering criteria into the AWK script. Again, use an existing filter as a template.
4. Create unit tests for the new filter by creating a new `test_that()` block in `tests/testthat/test_filters.r`. Again, use an existing filter as a template.
5. Update `README.md` and `vignettes/auk.Rmd` to add the new filter to the list of potential filters.
6. Build, test, check, and push to GitHub

### 1. Update `auk_ebd()`

Near the top of the code for `auk_ebd()`, a data frame named `filter_cols` is defined which specifies which columns have associated filters. Add a new row to this data frame and set `name` as the name of the column in the file header that will be filtered on and `id` as the name of the filter. For example, if you're creating a filter called `auk_color()` that filters on the column "FEATHER COLOR", then set `id = "color"` and `name = "feather color"`. Ideally, similar filters should be grouped together in this data frame, so insert the new row accordingly.

For filters that don't apply to the sampling event data file, i.e. filters at the species level rather than the checklist level, add the id to the character vector `not_in_sampling`. For example, modify the code to read: `not_in_sampling <- c("species", "breeding", "color")`.

Next, at the end of the code for `auk_ebd()`, the `auk_ebd` object is created and returned with the statement beginning with `structure(...`. This object should have placeholders for every filter. So, add a new element to the list, naming the variable after the `id` in the above data frame, putting it in the same order as in the above data frame, and choosing a sensible data type. For example, if `color` is a categorical variable, add a new list element `color = character()`, and if it's a numeric variable, add `color = numeric()`.

Finally, within `auk-ebd.r` a `print.auk_ebd()` method is defined, which you'll need to update to print the filter in a sensible way. Here you're best to find another filter with a similar format and use that as a template. Again, be sure to put the print code for the filter in the right order. For example, for a categorical filter allow multiple potential values, you may way something like:

```{r print-filter, eval=FALSE}
# color filter
cat("  Feather color: ")
if (length(x$filters$color) == 0) {
  cat("all")
} else {
  cat(paste(x$filters$color, collapse = ", "))
}
cat("\n")
```

### 2. Create filter function

Create a new function that will allow users to define a filter. Be sure to following the naming conventions used, for our color example, the function should be named `auk_color()` and it should be in a file called `auk-color.r`. It's easiest to use an existing function as a template here. In general, the function should take two argument, the `auk_ebd` object to modify, and an argument with the filter criteria, e.g. `auk_color(x, color)`. Note how the name of the function matches the name of the second argument. The function should be edited to include the following:

1. Extensive checks on the incoming arguments. Remember that filtering with AWK takes multiple hours, so it's best to catch any errors early, prior to running `auk_filter()`. At the very least, check data types and, where possible, check that values are valid (e.g. `color` should be in `c("red", "green", "blue", ...)`). Provide informative error or warning messages where appropriate.
2. Setting the filter criteria in the `auk_ebd` object. This is generally as simple as `x$filters$color = color`.
3. Thorough documentation. Document all the arguments and provide examples with and without the pipe operator (`%>%`).

### 3. Update `auk_filter()`

The actual work of filtering is done by `auk_filter()`, which generates an AWK script, then calls AWK. This function must be updated to parse the filters defined using the function you created in step 2 into AWK code. In the code for `auk_filter()`, there are two calls to the internal function `awk_translate()`, which is an internal function defined in the same file. It's `awk_translate()` that you'll need to edit. This function has a series of code blocks each of which prepares the AWK code for a different filter. Find an existing filter that is similar to the new one you're creating and copy it over to the correct spot (remember to preserve the ordering of the filters). For the `auk_color()` example, the code chunk would look like:

```{r awk-code, eval=FALSE}
  # color filter
  if (length(filters$color) == 0) {
    filter_strings$color <- ""
  } else {
    idx <- col_idx$index[col_idx$id == "color"]
    condition <- paste0("$", idx, " == \"", filters$color, "\"",
                        collapse = " || ")
    filter_strings$color <- str_interp(awk_if, list(condition = condition))
  }
```

When given a sampling event data file in addition to a EBD file, `auk_filter()` will filter both files. By default `auk_filter()` will apply all filters to both files, however, some filters (e.g. species) are only appropriate for the EBD. To address this, prior to calling `auk_translate()` for the sampling data, reset the species-specific filters. In the case of color, which is a species specific variable, modify the code as follows:

```{r species-specific, eval=FALSE}
s_filters <- x$filters
s_filters$species <- character()
## ADD THIS LINE
s_filters$color <- character()
##
awk_script_sampling <- awk_translate(filters = s_filters,
                                     col_idx = x$col_idx_sampling,
                                     sep = sep,
                                     select = select_cols)
```

Finally, at the end of the `auk-filter.r` file, there's a string named `awk_filter`, which defines the template for the AWK script. Each filter has a line in this string (e.g. `${species}`) where the AWK code will be inserted. You'll need to add a line in this file for your new filter: `${color}`.

### 4. Unit tests

Now that you've successfully created the filter, play around with it a bit to make sure it works as expected. Once you feel the filter is working, it's time to formalize this testing process by defining unit tests. Open the file `tests/testthat/test_filters.r` and you'll notice a series of calls like `test_that("auk_species", ...`, each of which contains tests for a specific filter.

Using an existing test block as an example, define a new block (again, put it in the correct order relative to the other filters). Consult the [Testing chapter](https://r-pkgs.org/tests.html) of Hadley Wickham's [R packages book](https://r-pkgs.org/) for details on defining good unit tests. At the very least, define tests to make sure that typical use works as expected, that errors are caught when input is invalid, and that edge cases are correctly handled.

### 5. Update vignette and README

Both the vignette (`vignettes/auk.Rmd`) and README (`README.Rmd`) have sections giving a short description of each filter. Add the new filter you've created here.

### 6. Build, test, check, and push to GitHub

Carry out the following final steps:

1. Run `devtools::document()` to generate package documentation
2. Run `devtools::build()` to build and install the package
3. Run `devtools::check()` to run the units tests and variety of other checks via `R CMD check`
5. Build the vignettes with `devtools::build_vignettes()`
6. Build the package website with `pkgdown::build_site()`
7. Commit to git, then push to GitHub

## CRAN submission

Minor updates to `auk` can be pushed to GitHub giving users the option of installing the development version from there. However, at least once a year, when a new eBird taxonomy is released, a new version of `auk` should be released on CRAN. For full details on this process, consult Hadley Wickham's [R Packages book](https://r-pkgs.org/release.html), however, I'll provide a quick guide here. Once The package has been updated following the instructions from the above sections:

1. Check the package. Run `devtools::check()` to run `R CMD check` locally. Check that a Windows binary can be built by running `devtools::build_win()`. The results will be emailed to you within about 30 minutes. Also, this package uses continuous integration to automatically check the package on Linux, Mac, and Windows whenever it's pushed to GitHub. Check the badges at the top of the GitHub repo to ensure the builds are passing. Any NOTEs, ERRORs, or WARNINGs returned by R CMD check must be fixed before submission to CRAN.
2. Increment the version number in the `DESCRIPTION` file.
3. Update `NEWS.md` to note any new features or changes.
4. Build the package with `devtools::build()`, the vignettes with `devtools::build_vignettes()`, and the website with `pkgdown::build_site()`.
5. Commit to git and push to GitHub.
6. Submit to CRAN with `devtools::release()`

At this point, you'll need to wait for binaries of your package to build, which could take a couple days. It's possible that problems will arise during this process and your package will be rejected, in which case, you'll need to fix any problems and resubmit.

Once the package is on CRAN, create a new release on GitHub and tag it with the version number.% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-date.R
\name{auk_date}
\alias{auk_date}
\title{Filter the eBird data by date}
\usage{
auk_date(x, date)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{date}{character or date; date range to filter by, provided either as a
character vector in the format \code{"2015-12-31"} or a vector of Date objects.
To filter on a range of dates, regardless of year, use \code{"*"} in place of
the year.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on a range of dates.
This function only defines the filter and, once all filters have been
defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and perform the
filtering.
}
\details{
To select observations from a range of dates, regardless of year,
the  wildcard \code{"*"} can be used in place of the year. For example, using
\code{date = c("*-05-01", "*-06-30")} will return observations from May and June
of \emph{any year}. When using wildcards, dates can wrap around the year end.

This function can also work with on an \code{auk_sampling} object if the user only
wishes to filter the sampling event data.
}
\examples{
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_date(date = c("2010-01-01", "2010-12-31"))
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_date(ebd, date = c("2010-01-01", "2010-12-31"))

# the * wildcard can be used in place of year to select dates from all years
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  # may-june records from all years
  auk_date(date = c("*-05-01", "*-06-30"))
  
# dates can also wrap around the end of the year
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  # dec-jan records from all years
  auk_date(date = c("*-12-01", "*-01-31"))
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-ebird-taxonomy.R
\name{get_ebird_taxonomy}
\alias{get_ebird_taxonomy}
\title{Get eBird taxonomy via the eBird API}
\usage{
get_ebird_taxonomy(version, locale)
}
\arguments{
\item{version}{integer; the version (i.e. year) of the taxonomy. The eBird
taxonomy is updated once a year in August. Leave this parameter blank to
get the current taxonomy.}

\item{locale}{character; the \href{https://support.ebird.org/support/solutions/articles/48000804865-bird-names-in-ebird}{locale for the common names},
defaults to English.}
}
\value{
A data frame of all species in the eBird taxonomy, consisting of the
following columns:
\itemize{
\item \code{scientific_name}: scientific name.
\item \code{common_name}: common name, defaults to English, but different languages
can be selected using the \code{locale} parameter.
\item \code{species_code}: a unique alphanumeric code identifying each species.
\item \code{category}: whether the entry is for a species or another
field-identifiable taxon, such as \code{spuh}, \code{slash}, \code{hybrid}, etc.
\item \code{taxon_order}: numeric value used to sort rows in taxonomic order.
\item \code{order}: the scientific name of the order that the species belongs to.
\item \code{family}: the scientific name of the family that the species belongs to.
\item \code{report_as}: for taxa that can be resolved to true species (i.e. species,
subspecies, and recognizable forms), this field links to the corresponding
species code. For taxa that can't be resolved, this field is \code{NA}.
}
}
\description{
Get the taxonomy used in eBird via the eBird API.
}
\examples{
\dontrun{
get_ebird_taxonomy()
}
}
\seealso{
Other helpers: 
\code{\link{auk_ebd_version}()},
\code{\link{auk_version}()},
\code{\link{ebird_species}()}
}
\concept{helpers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-select.R
\name{auk_select}
\alias{auk_select}
\title{Select a subset of columns}
\usage{
auk_select(x, select, file, sep = "\\t", overwrite = FALSE)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{select}{character; a character vector specifying the names of the
columns to select. Columns should be as they appear in the header of the
EBD; however, names are not case sensitive and spaces may be replaced by
underscores, e.g. \code{"COMMON NAME"}, \code{"common name"}, and \code{"common_NAME"} are
all valid.}

\item{file}{character; output file.}

\item{sep}{character; the input field separator, the eBird file is tab
separated by default. Must only be a single character and space delimited
is not allowed since spaces appear in many of the fields.}

\item{overwrite}{logical; overwrite output file if it already exists}
}
\value{
Invisibly returns the filename of the output file.
}
\description{
Select a subset of columns from the eBird Basic Dataset (EBD) or the sampling
events file. Subsetting the columns can significantly decrease file size.
}
\examples{
\dontrun{
# select a minimal set of columns
out_file <- tempfile()
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
cols <- c("latitude", "longitude",
          "group identifier", "sampling event identifier", 
          "scientific name", "observation count",
          "observer_id")
selected <- auk_select(ebd, select = cols, file = out_file)
str(read_ebd(selected))
}
}
\seealso{
Other text: 
\code{\link{auk_clean}()},
\code{\link{auk_split}()}
}
\concept{text}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-observer.R
\name{auk_observer}
\alias{auk_observer}
\title{Filter the eBird data by observer}
\usage{
auk_observer(x, observer_id)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{observer_id}{character or integer; observers to filter by. Observer IDs
can be provided either as integer (e.g. 12345) or character with the "obsr"
prefix as they appear in the EBD (e.g. "obsr12345").}
}
\value{
An \code{auk_ebd} or `auk_sampling`` object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on a set of
observer IDs This function only defines the filter and, once all filters have
been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and perform the
filtering.
}
\examples{
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_observer("obsr313215")
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_observer(ebd, observer = 313215)
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter-repeat-visits.R
\name{filter_repeat_visits}
\alias{filter_repeat_visits}
\title{Filter observations to repeat visits for hierarchical modeling}
\usage{
filter_repeat_visits(
  x,
  min_obs = 2L,
  max_obs = 10L,
  annual_closure = TRUE,
  n_days = NULL,
  date_var = "observation_date",
  site_vars = c("locality_id", "observer_id"),
  ll_digits = 6L
)
}
\arguments{
\item{x}{\code{data.frame}; observation data, e.g. data from the eBird Basic
Dataset (EBD) zero-filled with \code{\link[=auk_zerofill]{auk_zerofill()}}. This function will also
work with an \code{auk_zerofill} object, in which case it will be converted to
a data frame with \code{\link[=collapse_zerofill]{collapse_zerofill()}}.
\strong{Note that these data must for a single species}.}

\item{min_obs}{integer; minimum number of observations required for each
site.}

\item{max_obs}{integer; maximum number of observations allowed for each site.}

\item{annual_closure}{logical; whether the entire year should be treated as
the period of closure (the default). This can be useful, for example, if
the data have been subset to a period of closure prior to calling
\code{\link[=filter_repeat_visits]{filter_repeat_visits()}}.}

\item{n_days}{integer; number of days defining the temporal length of
closure. If \code{annual_closure = TRUE} closure periods will be split at year
boundaries. If \code{annual_closure = FALSE} the closure periods will ignore
year boundaries.}

\item{date_var}{character; column name of the variable in \code{x} containing the
date. This column should either be in \code{Date} format or convertible to
\code{Date} format with \code{\link[=as.Date]{as.Date()}}.}

\item{site_vars}{character; names of one of more columns in \code{x} that define a
site, typically the location (e.g. latitude/longitude) and observer ID.}

\item{ll_digits}{integer; the number of digits to round latitude and longitude
to. If latitude and/or longitude are used as \code{site_vars}, it's usually best
to round them prior to identifying sites, otherwise locations that are only
slightly offset (e.g. a few centimeters) will be treated as different. This
argument can also be used to group sites together that are close but not
identical. Note that 1 degree of latitude is approximately 100 km, so the
default value of 6 for \code{ll_digits} is equivalent to about 10 cm.}
}
\value{
A \code{data.frame} filtered to only retain observations from sites with
the allowed number of observations within the period of closure. The
results will be sorted such that sites are together and in chronological
order. The following variables are added to the data frame:
\itemize{
\item \code{site}: a unique identifier for each "site" corresponding to all the
variables in \code{site_vars} and \code{closure_id} concatenated together with
underscore separators.
\item \code{closure_id}: a unique ID for each closure period. If \code{annual_closure =   TRUE} this ID will include the year. If \code{n_days} is used an index given the
number of blocks of \code{n_days} days since the earliest observation will be
included. Note that in this case, there may be gaps in the IDs.
\item \code{n_observations}: number of observations at each site after all
filtering.
}
}
\description{
Hierarchical modeling of abundance and occurrence requires repeat visits to
sites to estimate detectability. These visits should be all be within a
period of closure, i.e. when the population can be assumed to be closed.
eBird data, and many other data sources, do not explicitly follow this
protocol; however, subsets of the data can be extracted to produce data
suitable for hierarchical modeling. This function extracts a subset of
observation data that have a desired number of repeat visits within a period
of closure.
}
\details{
In addition to specifying the minimum and maximum number of
observations per site, users must specify the variables in the dataset that
define a "site". This is typically a combination of IDs defining the
geographic site and the unique observer (repeat visits are meant to be
conducted by the same observer). Finally, the closure period must be
defined, which is a period within which the population of the focal species
can reasonably be assumed to be closed. This can be done using a
combination of the \code{n_days} and \code{annual_closure} arguments.
}
\examples{
# read and zero-fill the ebd data
f_ebd <- system.file("extdata/zerofill-ex_ebd.txt", package = "auk")
f_smpl <- system.file("extdata/zerofill-ex_sampling.txt", package = "auk")
# data must be for a single species
ebd_zf <- auk_zerofill(x = f_ebd, sampling_events = f_smpl,
                       species = "Collared Kingfisher",
                       collapse = TRUE)
filter_repeat_visits(ebd_zf, n_days = 30)
}
\seealso{
Other modeling: 
\code{\link{format_unmarked_occu}()}
}
\concept{modeling}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-bcr.R
\name{auk_bcr}
\alias{auk_bcr}
\title{Filter the eBird data by Bird Conservation Region}
\usage{
auk_bcr(x, bcr, replace = FALSE)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{bcr}{integer; BCRs to filter by. BCRs are identified by an integer,
from 1 to 66, that can be looked up in the \link{bcr_codes} table.}

\item{replace}{logical; multiple calls to \code{auk_state()} are additive,
unless \code{replace = FALSE}, in which case the previous list of states to
filter by will be removed and replaced by that in the current call.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) to extract data for a set
of \href{https://nabci-us.org/resources/bird-conservation-regions/}{Bird Conservation Regions} (BCRs).
BCRs are ecologically distinct regions in North America with similar bird
communities, habitats, and resource management issues. This function only
defines the filter and, once all filters have been defined, \code{\link[=auk_filter]{auk_filter()}}
should be used to call AWK and perform the filtering.
}
\details{
This function can also work with on an \code{auk_sampling} object if the
user only wishes to filter the sampling event data.
}
\examples{
# bcr codes can be looked up in bcr_codes
dplyr::filter(bcr_codes, bcr_name == "Central Hardwoods")
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_bcr(bcr = 24)
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_bcr(ebd, bcr = 24)
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-bbox.R
\name{auk_extent}
\alias{auk_extent}
\title{Filter the eBird data by spatial extent}
\usage{
auk_extent(x, extent)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{extent}{numeric; spatial extent expressed as the range of latitudes and
longitudes in decimal degrees: \code{c(lng_min, lat_min, lng_max, lat_max)}.
Note that longitudes in the Western Hemisphere and latitudes sound of the
equator should be given as negative numbers.}
}
\value{
An \code{auk_ebd} object.
}
\description{
\strong{Deprecated}, use \code{\link[=auk_bbox]{auk_bbox()}} instead.
}
\examples{
# fliter to locations roughly in the Pacific Northwest
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_bbox(bbox = c(-125, 37, -120, 52))
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_bbox(ebd, bbox = c(-125, 37, -120, 52))
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format-unmarked-occu.R
\name{format_unmarked_occu}
\alias{format_unmarked_occu}
\title{Format EBD data for occupancy modeling with \code{unmarked}}
\usage{
format_unmarked_occu(
  x,
  site_id = "site",
  response = "species_observed",
  site_covs,
  obs_covs
)
}
\arguments{
\item{x}{\code{data.frame}; observation data, e.g. from the eBird Basic Dataset
(EBD), for \strong{a single species}, that has been filtered to those with
repeat visits by \code{\link[=filter_repeat_visits]{filter_repeat_visits()}}.}

\item{site_id}{character; a unique idenitifer for each "site", typically
identifying observations from a unique location by the same observer
within a period of temporal closure. Data output from
\code{\link[=filter_repeat_visits]{filter_repeat_visits()}} will have a \code{.site_id} variable that meets these
requirements.}

\item{response}{character; the variable that will act as the response in
modeling efforts, typically a binary variable indicating presence or
absence or a count of individuals seen.}

\item{site_covs}{character; the variables that will act as site-level
covariates, i.e. covariates that vary at the site level, for example,
latitude/longitude or habitat predictors. If this parameter is missing, it
will be assumed that any variable that is not an observation-level
covariate (\code{obs_covs}) or the \code{site_id}, is a site-level covariate.}

\item{obs_covs}{character; the variables that will act as observation-level
covariates, i.e. covariates that vary within sites, at the level of
observations, for example, time or length of observation.}
}
\value{
A data frame that can be processed by \code{\link[unmarked:formatWideLong]{unmarked::formatWide()}}.
Each row will correspond to a unqiue site and, assuming there are a maximum
of \code{N} observations per site, columns will be as follows:
\enumerate{
\item The unique site identifier, named "site".
\item \code{N} response columns, one for each observation, named "y.1", ..., "y.N".
\item Columns for each of the site-level covariates.
\item Groups of \code{N} columns of observation-level covariates, one column per
covariate per observation, names "covariate_name.1", ...,
"covariate_name.N".
}
}
\description{
Prepare a data frame of species observations for ingestion into the package
\code{unmarked} for hierarchical modeling of abundance and occurrence. The
function \code{\link[unmarked:formatWideLong]{unmarked::formatWide()}} takes a data frame and converts it to one
of several \code{unmarked} objects, which can then be used for modeling. This
function converts data from a format in which each row is an observation
(e.g. as in the eBird Basic Dataset) to the esoteric format required by
\code{\link[unmarked:formatWideLong]{unmarked::formatWide()}} in which each row is a site.
}
\details{
Hierarchical modeling requires repeat observations at each "site" to
estimate detectability. A "site" is typically defined as a geographic
location visited by the same observer within a period of temporal closure.
To define these sites and filter out observations that do not correspond to
repeat visits, users should use \code{\link[=filter_repeat_visits]{filter_repeat_visits()}}, then pass the
output to this function.

\code{\link[=format_unmarked_occu]{format_unmarked_occu()}} is designed to prepare data to be converted into
an \code{unmarkedFrameOccu} object for occupancy modeling with
\code{\link[unmarked:occu]{unmarked::occu()}}; however, it can also be used to prepare data for
conversion to an \code{unmarkedFramePCount} object for abundance modeling with
\code{\link[unmarked:pcount]{unmarked::pcount()}}.
}
\examples{
# read and zero-fill the ebd data
f_ebd <- system.file("extdata/zerofill-ex_ebd.txt", package = "auk")
f_smpl <- system.file("extdata/zerofill-ex_sampling.txt", package = "auk")
# data must be for a single species
ebd_zf <- auk_zerofill(x = f_ebd, sampling_events = f_smpl,
                       species = "Collared Kingfisher",
                       collapse = TRUE)
occ <- filter_repeat_visits(ebd_zf, n_days = 30)
# format for unmarked
# typically one would join in habitat covariates prior to this step
occ_wide <- format_unmarked_occu(occ,
                                 response = "species_observed",
                                 site_covs = c("latitude", "longitude"),
                                 obs_covs = c("effort_distance_km", 
                                              "duration_minutes"))
# create an unmarked object
if (requireNamespace("unmarked", quietly = TRUE)) {
  occ_um <- unmarked::formatWide(occ_wide, type = "unmarkedFrameOccu")
  unmarked::summary(occ_um)
}

# this function can also be used for abundance modeling
abd <- ebd_zf \%>\% 
  # convert count to integer, drop records with no count
  dplyr::mutate(observation_count = as.integer(observation_count)) \%>\% 
  dplyr::filter(!is.na(observation_count)) \%>\% 
  # filter to repeated visits
  filter_repeat_visits(n_days = 30)
# prepare for conversion to unmarkedFramePCount object
abd_wide <- format_unmarked_occu(abd,
                                 response = "observation_count",
                                 site_covs = c("latitude", "longitude"),
                                 obs_covs = c("effort_distance_km", 
                                              "duration_minutes"))
# create an unmarked object
if (requireNamespace("unmarked", quietly = TRUE)) {
  abd_um <- unmarked::formatWide(abd_wide, type = "unmarkedFrameOccu")
  unmarked::summary(abd_um)
}
}
\seealso{
Other modeling: 
\code{\link{filter_repeat_visits}()}
}
\concept{modeling}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-state.R
\name{auk_state}
\alias{auk_state}
\title{Filter the eBird data by state}
\usage{
auk_state(x, state, replace = FALSE)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{state}{character; states to filter by. eBird uses 4 to 6 character
state codes consisting of two parts, the 2-letter ISO country code and a
1-3 character state code, separated by a dash. For example, \code{"US-NY"}
corresponds to New York State in the United States. Refer to the data frame
\link{ebird_states} for look up state codes.}

\item{replace}{logical; multiple calls to \code{auk_state()} are additive,
unless \code{replace = FALSE}, in which case the previous list of states to
filter by will be removed and replaced by that in the current call.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on a set of
states. This function only defines the filter and, once all filters have
been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and perform the
filtering.
}
\details{
It is not possible to filter by both country and state, so calling
\code{auk_state()} will reset the country filter to all countries, and vice versa.

This function can also work with on an \code{auk_sampling} object if the user only
wishes to filter the sampling event data.
}
\examples{
# state codes for a given country can be looked up in ebird_states
dplyr::filter(ebird_states, country == "Costa Rica")
# choose texas, united states and puntarenas, cost rica
states <- c("US-TX", "CR-P")
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_state(states)
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_state(ebd, states)
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-rollup.R
\name{auk_rollup}
\alias{auk_rollup}
\title{Roll up eBird taxonomy to species}
\usage{
auk_rollup(x, taxonomy_version, drop_higher = TRUE)
}
\arguments{
\item{x}{data.frame; data frame of eBird data, typically as imported by
\code{\link[=read_ebd]{read_ebd()}}}

\item{taxonomy_version}{integer; the version (i.e. year) of the taxonomy. In
most cases, this should be left empty to use the version of the taxonomy
included in the package. See \code{\link[=get_ebird_taxonomy]{get_ebird_taxonomy()}}.}

\item{drop_higher}{logical; whether to remove taxa above species during the
rollup process, e.g. "spuhs" like "duck sp.".}
}
\value{
A data frame of the eBird data with taxonomic rollup applied.
}
\description{
The eBird Basic Dataset (EBD) includes both true species and every other
field-identifiable taxon that could be relevant for birders to report. This
includes taxa not identifiable to a species (e.g. hybrids) and taxa reported
below the species level (e.g. subspecies). This function produces a list of
observations of true species, by removing the former and rolling the latter
up to the species level. In the resulting EBD data.frame,
\code{category} will be \code{"species"} for all records and the subspecies fields will
be dropped. By default, \code{\link[=read_ebd]{read_ebd()}} calls \code{ebd_rollup()} when importing an
eBird data file.
}
\details{
When rolling observations up to species level the observed counts
are summed across any taxa that resolve to the same species. However, if
any of these taxa have a count of "X" (i.e. the observer did not enter a
count), then the rolled up record will get an "X" as well. For example, if
an observer saw 3 Myrtle and 2 Audubon's Warblers, this will roll up to 5
Yellow-rumped Warblers. However, if an "X" was entered for Myrtle, this
would roll up to "X" for Yellow-rumped Warbler.

The eBird taxonomy groups taxa into eight different categories. These
categories, and the way they are treated by \code{\link[=auk_rollup]{auk_rollup()}} are as follows:
\itemize{
\item \strong{Species:} e.g., Mallard. Combined with lower level taxa if present on
the same checklist.
\item \strong{ISSF or Identifiable Sub-specific Group:} Identifiable subspecies or
group of subspecies, e.g., Mallard (Mexican). Rolled-up to species level.
\item \strong{Intergrade:} Hybrid between two ISSF (subspecies or subspecies
groups), e.g., Mallard (Mexican intergrade. Rolled-up to species level.
\item \strong{Form:} Miscellaneous other taxa, including recently-described species
yet to be accepted or distinctive forms that are not universally accepted
(Red-tailed Hawk (Northern), Upland Goose (Bar-breasted)). If the checklist
contains multiple taxa corresponding to the same species, the lower level
taxa are rolled up, otherwise these records are left as is.
\item \strong{Spuh:}  Genus or identification at broad level -- e.g., duck sp.,
dabbling duck sp.. Dropped by \code{auk_rollup()}.
\item \strong{Slash:} Identification to Species-pair e.g., American Black
Duck/Mallard). Dropped by \code{auk_rollup()}.
\item \strong{Hybrid:} Hybrid between two species, e.g., American Black Duck x
Mallard (hybrid). Dropped by \code{auk_rollup()}.
\item \strong{Domestic:} Distinctly-plumaged domesticated varieties that may be
free-flying (these do not count on personal lists) e.g., Mallard (Domestic
type). Dropped by \code{auk_rollup()}.
}

The rollup process is based on the eBird taxonomy, which is updated once a
year in August. The \code{auk} package includes a copy of the eBird taxonomy,
current at the time of release; however, if the EBD and \code{auk} versions are
not aligned, you may need to explicitly specify which version of the
taxonomy to use, in which case the eBird API will be queried to get the
correct version of the taxonomy.
}
\examples{
# get the path to the example data included in the package
# in practice, provide path to ebd, e.g. f <- "data/ebd_relFeb-2018.txt
f <- system.file("extdata/ebd-rollup-ex.txt", package = "auk")
# read in data without rolling up
ebd <- read_ebd(f, rollup = FALSE)
# rollup
ebd_ru <- auk_rollup(ebd)
# keep higher taxa
ebd_higher <- auk_rollup(ebd, drop_higher = FALSE)

# all taxa not identifiable to species are dropped
unique(ebd$category)
unique(ebd_ru$category)
unique(ebd_higher$category)

# yellow-rump warbler subspecies rollup
library(dplyr)
# without rollup, there are three observations
ebd \%>\%
  filter(common_name == "Yellow-rumped Warbler") \%>\% 
  select(checklist_id, category, common_name, subspecies_common_name, 
         observation_count)
# with rollup, they have been combined
ebd_ru \%>\%
  filter(common_name == "Yellow-rumped Warbler") \%>\% 
  select(checklist_id, category, common_name, observation_count)
}
\references{
Consult the \href{https://ebird.org/science/use-ebird-data/the-ebird-taxonomy}{eBird taxonomy}
page for further details.
}
\seealso{
Other pre: 
\code{\link{auk_unique}()}
}
\concept{pre}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{valid_protocols}
\alias{valid_protocols}
\title{Valid Protocols}
\format{
A vector with 42 elements.
}
\usage{
valid_protocols
}
\description{
A vector of valid protocol names, e.g. "Traveling", "Stationary", etc.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-protocol.R
\name{auk_protocol}
\alias{auk_protocol}
\title{Filter the eBird data by protocol}
\usage{
auk_protocol(x, protocol)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{protocol}{character. Many protocols exist in the database, however, the
most commonly used are:
\itemize{
\item Stationary
\item Traveling
\item Area
\item Incidental
}

A complete list of valid protocols is contained within the vector
\code{valid_protocols} within this package. Multiple protocols are allowed at
the same time.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Filter to just data collected following a specific search protocol:
stationary, traveling, or casual. This function only defines the filter and,
once all filters have been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK
and perform the filtering.
}
\details{
This function can also work with on an \code{auk_sampling} object if the
user only wishes to filter the sampling event data.
}
\examples{
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_protocol("Stationary")
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_protocol(ebd, "Stationary")
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ebird_taxonomy}
\alias{ebird_taxonomy}
\title{eBird Taxonomy}
\format{
A data frame with eight variables and 16,248 rows:
\itemize{
\item \code{scientific_name}: scientific name.
\item \code{common_name}: common name, defaults to English, but different languages
can be selected using the \code{locale} parameter.
\item \code{species_code}: a unique alphanumeric code identifying each species.
\item \code{category}: whether the entry is for a species or another
field-identifiable taxon, such as \code{spuh}, \code{slash}, \code{hybrid}, etc.
\item \code{taxon_order}: numeric value used to sort rows in taxonomic order.
\item \code{order}: the scientific name of the order that the species belongs to.
\item \code{family}: the scientific name of the family that the species belongs to.
\item \code{report_as}: for taxa that can be resolved to true species (i.e. species,
subspecies, and recognizable forms), this field links to the corresponding
species code. For taxa that can't be resolved, this field is \code{NA}.
}

For further details, see \url{https://support.ebird.org/support/solutions/articles/48000837816-the-ebird-taxonomy}
}
\usage{
ebird_taxonomy
}
\description{
A simplified version of the taxonomy used by eBird. Includes proper species
as well as various other categories such as \code{spuh} (e.g. \emph{duck sp.}) and
\emph{slash} (e.g. \emph{American Black Duck/Mallard}). This taxonomy is based on the
Clements Checklist, which is updated annually, typically in the late summer.
Non-ASCII characters (e.g. those with accents) have been converted to ASCII
equivalents in this data frame.
}
\seealso{
Other data: 
\code{\link{bcr_codes}},
\code{\link{ebird_states}}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-ebd-version.R
\name{auk_ebd_version}
\alias{auk_ebd_version}
\title{Get the EBD version and associated taxonomy version}
\usage{
auk_ebd_version(x, check_exists = TRUE)
}
\arguments{
\item{x}{filename of EBD of sampling event data file, \code{auk_ebd} object, or
\code{auk_sampling} object.}

\item{check_exists}{logical; should the file be checked for existence before
processing. If \code{check_exists = TRUE} and the file does not exists, the
function will raise an error.}
}
\value{
A list with two elements:
\itemize{
\item \code{ebd_version}: a date object specifying the release date of the EBD.
\item \code{taxonomy_version}: the year of the taxonomy used in this EBD.
}

Both elements will be NA if an EBD version cannot be extracted from the
filename.
}
\description{
Based on the filename of eBird Basic Dataset (EBD) or sampling event data,
determine the version (i.e. release date) of this EBD. Also determine the
corresponding taxonomy version. The eBird taxonomy is updated annually in
August.
}
\examples{
auk_ebd_version("ebd_relAug-2018.txt", check_exists = FALSE)
}
\seealso{
Other helpers: 
\code{\link{auk_version}()},
\code{\link{ebird_species}()},
\code{\link{get_ebird_taxonomy}()}
}
\concept{helpers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-version.R
\name{auk_version}
\alias{auk_version}
\title{Versions of auk, the EBD, and the eBird taxonomy}
\usage{
auk_version()
}
\value{
A list with three elements:
\itemize{
\item \code{auk_version}: the version of \code{auk}, e.g. \code{"auk 0.4.1"}.
\item \code{ebd_version}: a date object specifying the release date of the EBD
version that this \code{auk} version is designed to work with.
\item \code{taxonomy_version}: the year of the taxonomy built in to this version of
\code{auk}, i.e. the one stored in \link{ebird_taxonomy}.
}
}
\description{
This package depends on the version of the EBD and on the eBird taxonomy. Use
this function to determine the currently installed version of \code{auk}, the
version of the EBD that this \code{auk} version works with, and the version of the
eBird taxonomy included in the packages. The EBD is update quarterly, in
March, June, September, and December, while the taxonomy is updated annually
in August or September. To ensure proper functioning, always use the latest
version of the auk package and the EBD.
}
\examples{
auk_version()
}
\seealso{
Other helpers: 
\code{\link{auk_ebd_version}()},
\code{\link{ebird_species}()},
\code{\link{get_ebird_taxonomy}()}
}
\concept{helpers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-project.R
\name{auk_project}
\alias{auk_project}
\title{Filter the eBird data by project code}
\usage{
auk_project(x, project)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{project}{character; project code to filter by (e.g. \code{"EBIRD_MEX"}).
Multiple codes are accepted.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Some eBird records are collected as part of a particular project (e.g. the
Virginia Breeding Bird Survey) and have an associated project code in the
eBird dataset (e.g. EBIRD_ATL_VA). This function only defines the filter and,
once all filters have been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK
and perform the filtering.
}
\details{
This function can also work with on an \code{auk_sampling} object if the
user only wishes to filter the sampling event data.
}
\examples{
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_project("EBIRD_MEX")
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_project(ebd, "EBIRD_MEX")
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{bcr_codes}
\alias{bcr_codes}
\title{BCR Codes}
\format{
A data frame with two variables and 66 rows:
\itemize{
\item \code{bcr_code}: integer code from 1 to 66.
\item \code{bcr_name}: name of BCR.
}
}
\usage{
bcr_codes
}
\description{
A data frame of Bird Conservation Region (BCR) codes. BCRs are ecologically
distinct regions in North America with similar bird communities, habitats,
and resource management issues. These codes are required to filter by BCR
using \code{\link[=auk_bcr]{auk_bcr()}}.
}
\seealso{
Other data: 
\code{\link{ebird_states}},
\code{\link{ebird_taxonomy}}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-set-awk-path.R
\name{auk_set_awk_path}
\alias{auk_set_awk_path}
\title{Set a custom path to AWK executable}
\usage{
auk_set_awk_path(path, overwrite = FALSE)
}
\arguments{
\item{path}{character; path to the AWK executable on your system, e.g.
\code{"C:/cygwin64/bin/gawk.exe"} or \code{"/usr/bin/awk"}.}

\item{overwrite}{logical; should the existing \code{AWK_PATH} be overwritten if it
has already been set in .Renviron.}
}
\value{
Edits .Renviron, then returns the AWK path invisibly.
}
\description{
If AWK has been installed in a non-standard location, the environment
variable \code{AWK_PATH} must be set to specify the location of the executable.
Use this function to set \code{AWK_PATH} in your .Renviron file. \strong{Most users
should NOT set \code{AWK_PATH}, only do so if you have installed AWK in
non-standard location and \code{auk} cannot find it.}
}
\examples{
\dontrun{
auk_set_awk_path("/usr/bin/awk")
}
}
\seealso{
Other paths: 
\code{\link{auk_get_awk_path}()},
\code{\link{auk_get_ebd_path}()},
\code{\link{auk_set_ebd_path}()}
}
\concept{paths}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-clean.R
\name{auk_clean}
\alias{auk_clean}
\title{Clean an eBird data file (Deprecated)}
\usage{
auk_clean(f_in, f_out, sep = "\\t", remove_text = FALSE, overwrite = FALSE)
}
\arguments{
\item{f_in}{character; input file. If file is not found as specified, it will
be looked for in the directory specified by the \code{EBD_PATH} environment
variable.}

\item{f_out}{character; output file.}

\item{sep}{character; the input field separator, the basic dataset is tab
separated by default. Must only be a single character and space delimited
is not allowed since spaces appear in many of the fields.}

\item{remove_text}{logical; whether all free text entry columns should be
removed. These columns include comments, location names, and observer
names. These columns cause import errors due to special characters and
increase the file size, yet are rarely valuable for analytical
applications, so may be removed. Setting this argument to \code{TRUE} can lead
to a significant reduction in file size.}

\item{overwrite}{logical; overwrite output file if it already exists.}
}
\value{
If AWK ran without errors, the output filename is returned, however,
if an error was encountered the exit code is returned.
}
\description{
This function is no longer required by current versions of the eBird Basic
Dataset (EBD).
}
\examples{
\dontrun{
# get the path to the example data included in the package
f <- system.file("extdata/ebd-sample.txt", package = "auk")
# output to a temp file for example
# in practice, provide path to output file
# e.g. f_out <- "output/ebd_clean.txt"
f_out <- tempfile()

# clean file to remove problem rows
# note: this function is deprecated and no longer does anything
auk_clean(f, f_out)
}
}
\seealso{
Other text: 
\code{\link{auk_select}()},
\code{\link{auk_split}()}
}
\concept{text}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-zerofill.R
\name{auk_zerofill}
\alias{auk_zerofill}
\alias{auk_zerofill.data.frame}
\alias{auk_zerofill.character}
\alias{auk_zerofill.auk_ebd}
\alias{collapse_zerofill}
\title{Read and zero-fill an eBird data file}
\usage{
auk_zerofill(x, ...)

\method{auk_zerofill}{data.frame}(
  x,
  sampling_events,
  species,
  taxonomy_version,
  collapse = FALSE,
  unique = TRUE,
  rollup = TRUE,
  drop_higher = TRUE,
  complete = TRUE,
  ...
)

\method{auk_zerofill}{character}(
  x,
  sampling_events,
  species,
  taxonomy_version,
  collapse = FALSE,
  unique = TRUE,
  rollup = TRUE,
  drop_higher = TRUE,
  complete = TRUE,
  sep = "\\t",
  ...
)

\method{auk_zerofill}{auk_ebd}(
  x,
  species,
  taxonomy_version,
  collapse = FALSE,
  unique = TRUE,
  rollup = TRUE,
  drop_higher = TRUE,
  complete = TRUE,
  sep = "\\t",
  ...
)

collapse_zerofill(x)
}
\arguments{
\item{x}{filename, \code{data.frame} of eBird observations, or \code{auk_ebd} object
with associated output files as created by \code{\link[=auk_filter]{auk_filter()}}. If a filename is
provided, it must point to the EBD and the \code{sampling_events} argument must
point to the sampling event data file. If a \code{data.frame} is provided it
should have been imported with \code{\link[=read_ebd]{read_ebd()}}, to ensure the variables names
have been set correctly, and it must have been passed through
\code{\link[=auk_unique]{auk_unique()}} to ensure duplicate group checklists have been removed.}

\item{...}{additional arguments passed to methods.}

\item{sampling_events}{character or \code{data.frame}; filename for the sampling
event data or a \code{data.frame} of the same data. If a \code{data.frame} is
provided it should have been imported with \code{\link[=read_sampling]{read_sampling()}}, to ensure the
variables names have been set correctly, and it must have been passed
through \code{\link[=auk_unique]{auk_unique()}} to ensure duplicate group checklists have been
removed.}

\item{species}{character; species to include in zero-filled dataset, provided
as scientific or English common names, or a mixture of both. These names
must match the official eBird Taxomony (\link{ebird_taxonomy}). To include all
species, leave this argument blank.}

\item{taxonomy_version}{integer; the version (i.e. year) of the taxonomy. In
most cases, this should be left empty to use the version of the taxonomy
included in the package. See \code{\link[=get_ebird_taxonomy]{get_ebird_taxonomy()}}.}

\item{collapse}{logical; whether to call \code{collapse_zerofill()} to return a
data frame rather than an \code{auk_zerofill} object.}

\item{unique}{logical; should \code{\link[=auk_unique]{auk_unique()}} be run on the input data if it
hasn't already.}

\item{rollup}{logical; should \code{\link[=auk_rollup]{auk_rollup()}} be run on the input data if it
hasn't already.}

\item{drop_higher}{logical; whether to remove taxa above species during the
rollup process, e.g. "spuhs" like "duck sp.". See \code{\link[=auk_rollup]{auk_rollup()}}.}

\item{complete}{logical; if \code{TRUE} (the default) all checklists are required
to be complete prior to zero-filling.}

\item{sep}{character; single character used to separate fields within a row.}
}
\value{
By default, an \code{auk_zerofill} object, or a data frame if \code{collapse = TRUE}.
}
\description{
Read an eBird Basic Dataset (EBD) file, and associated sampling event data
file, to produce a zero-filled, presence-absence dataset. The EBD contains
bird sightings and the sampling event data is a set of all checklists, they
can be combined to infer absence data by assuming any species not reported on
a checklist was had a count of zero.
}
\details{
\code{auk_zerofill()} generates an \code{auk_zerofill} object consisting of a list with
elements \code{observations} and \code{sampling_events}. \code{observations} is a data frame
giving counts and binary presence/absence data for each species.
\code{sampling_events} is a data frame with checklist level information. The two
data frames can be connected via the \code{checklist_id} field. This format is
efficient for storage since the checklist columns are not duplicated for each
species, however, working with the data often requires joining the two data
frames together.

To return a data frame, set \code{collapse = TRUE}. Alternatively,
\code{zerofill_collapse()} generates a data frame from an \code{auk_zerofill} object,
by joining the two data frames together to produce a single data frame in
which each row provides both checklist and species information for a
sighting.

The list of species is checked against the eBird taxonomy for validity. This
taxonomy is updated once a year in August. The \code{auk} package includes a copy
of the eBird taxonomy, current at the time of release; however, if the EBD
and \code{auk} versions are not aligned, you may need to explicitly specify which
version of the taxonomy to use, in which case the eBird API will be queried
to get the correct version of the taxonomy.
}
\section{Methods (by class)}{
\itemize{
\item \code{data.frame}: EBD data frame.

\item \code{character}: Filename of EBD.

\item \code{auk_ebd}: \code{auk_ebd} object output from \code{\link[=auk_filter]{auk_filter()}}. Must
have had a sampling event data file set in the original call to
\code{\link[=auk_ebd]{auk_ebd()}}.
}}

\examples{
# read and zero-fill the ebd data
f_ebd <- system.file("extdata/zerofill-ex_ebd.txt", package = "auk")
f_smpl <- system.file("extdata/zerofill-ex_sampling.txt", package = "auk")
auk_zerofill(x = f_ebd, sampling_events = f_smpl)

# use the species argument to only include a subset of species
auk_zerofill(x = f_ebd, sampling_events = f_smpl,
             species = "Collared Kingfisher")

# to return a data frame use collapse = TRUE
ebd_df <- auk_zerofill(x = f_ebd, sampling_events = f_smpl, collapse = TRUE)
}
\seealso{
Other import: 
\code{\link{read_ebd}()}
}
\concept{import}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-get-awk-path.R
\name{auk_get_awk_path}
\alias{auk_get_awk_path}
\title{OS specific path to AWK executable}
\usage{
auk_get_awk_path()
}
\value{
Path to AWK or \code{NA} if AWK wasn't found.
}
\description{
Return the OS specific path to AWK (e.g. \code{"C:/cygwin64/bin/gawk.exe"} or
\code{"/usr/bin/awk"}), or highlights if it's not installed. To manually set the
path to AWK, set the \code{AWK_PATH} environment variable in your \code{.Renviron}
file, which can be accomplished with the helper function
\code{auk_set_awk_path(path)}.
}
\examples{
auk_get_awk_path()
}
\seealso{
Other paths: 
\code{\link{auk_get_ebd_path}()},
\code{\link{auk_set_awk_path}()},
\code{\link{auk_set_ebd_path}()}
}
\concept{paths}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-complete.R
\name{auk_complete}
\alias{auk_complete}
\title{Filter out incomplete checklists from the eBird data}
\usage{
auk_complete(x)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) to only keep complete
checklists, i.e. those for which all birds seen or heard were recorded. These
checklists are the most valuable for scientific uses since they provide
presence and absence data.This function only defines the filter and, once all
filters have been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and
perform the filtering.
}
\details{
This function can also work with on an \code{auk_sampling} object if the
user only wishes to filter the sampling event data.
}
\examples{
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_complete()
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-breeding.R
\name{auk_breeding}
\alias{auk_breeding}
\title{Filter to only include observations with breeding codes}
\usage{
auk_breeding(x)
}
\arguments{
\item{x}{\code{auk_ebd} object; reference to basic dataset file created by
\code{\link[=auk_ebd]{auk_ebd()}}.}
}
\value{
An \code{auk_ebd} object.
}
\description{
eBird users have the option of specifying breeding bird atlas codes for their
observations, for example, if nesting building behaviour is observed. Use
this filter to select only those observations with an associated breeding
code. This function only defines the filter and, once all filters have been
defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and perform the filtering.
}
\examples{
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_breeding()
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-split.R
\name{auk_split}
\alias{auk_split}
\title{Split an eBird data file by species}
\usage{
auk_split(
  file,
  species,
  prefix,
  taxonomy_version,
  sep = "\\t",
  overwrite = FALSE
)
}
\arguments{
\item{file}{character; input file.}

\item{species}{character; species to filter and split by, provided as
scientific or English common names, or a mixture of both. These names must
match the official eBird Taxomony (\link{ebird_taxonomy}).}

\item{prefix}{character; a file and directory prefix. For example, if
splitting by species "A" and "B" and \code{prefix = "data/ebd_"}, the resulting
files will be "data/ebd_A.txt" and "data/ebd_B.txt".}

\item{taxonomy_version}{integer; the version (i.e. year) of the taxonomy. In
most cases, this should be left empty to use the version of the taxonomy
included in the package. See \code{\link[=get_ebird_taxonomy]{get_ebird_taxonomy()}}.}

\item{sep}{character; the input field separator, the eBird file is tab
separated by default. Must only be a single character and space delimited
is not allowed since spaces appear in many of the fields.}

\item{overwrite}{logical; overwrite output files if they already exists.}
}
\value{
A vector of output filenames, one for each species.
}
\description{
Given an eBird Basic Dataset (EBD) and a list of species, split the file into
multiple text files, one for each species. This function is typically used
after \code{\link[=auk_filter]{auk_filter()}} has been applied if the resulting file is too large to
be read in all at once.
}
\details{
The list of species is checked against the eBird taxonomy for
validity. This taxonomy is updated once a year in August. The \code{auk} package
includes a copy of the eBird taxonomy, current at the time of release;
however, if the EBD and \code{auk} versions are not aligned, you may need to
explicitly specify which version of the taxonomy to use, in which case
the eBird API will be queried to get the correct version of the taxonomy.
}
\examples{
\dontrun{
species <- c("Canada Jay", "Cyanocitta stelleri")
# get the path to the example data included in the package
# in practice, provide path to a filtered ebd file
# e.g. f <- "data/ebd_filtered.txt
f <- system.file("extdata/ebd-sample.txt", package = "auk")
# output to a temporary directory for example
# in practice, provide the path to the output location
# e.g. prefix <- "output/ebd_"
prefix <- file.path(tempdir(), "ebd_")
species_files <- auk_split(f, species = species, prefix = prefix)
}
}
\seealso{
Other text: 
\code{\link{auk_clean}()},
\code{\link{auk_select}()}
}
\concept{text}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-last-edited.R
\name{auk_last_edited}
\alias{auk_last_edited}
\title{Filter the eBird data by last edited date}
\usage{
auk_last_edited(x, date)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{date}{character or date; date range to filter by, provided either as a
character vector in the format \code{"2015-12-31"} or a vector of Date objects.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on a range of last
edited dates. Last edited date is typically used to extract just new or
recently edited data. This function only defines the filter and, once all
filters have been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and
perform the filtering.
}
\details{
This function can also work with on an \code{auk_sampling} object if the
user only wishes to filter the sampling event data.
}
\examples{
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_last_edited(date = c("2010-01-01", "2010-12-31"))
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-bbox.R
\name{auk_bbox}
\alias{auk_bbox}
\title{Filter the eBird data by spatial bounding box}
\usage{
auk_bbox(x, bbox)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{bbox}{numeric or \code{sf} or \verb{Raster*} object; spatial bounding box
expressed as the range of latitudes and longitudes in decimal degrees:
\code{c(lng_min, lat_min, lng_max, lat_max)}. Note that longitudes in the
Western Hemisphere and latitudes sound of the equator should be given as
negative numbers. Alternatively, a spatial object from either the \code{sf} or
\code{raster} packages can be provided and the bounding box will be extracted
from this object.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on spatial bounding
box. This function only defines the filter and, once all filters have been
defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and perform the filtering.
}
\details{
This function can also work with on an \code{auk_sampling} object if the
user only wishes to filter the sampling event data.
}
\examples{
# fliter to locations roughly in the Pacific Northwest
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_bbox(bbox = c(-125, 37, -120, 52))
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_bbox(ebd, bbox = c(-125, 37, -120, 52))
}
\seealso{
Other filter: 
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-country.R
\name{auk_country}
\alias{auk_country}
\title{Filter the eBird data by country}
\usage{
auk_country(x, country, replace = FALSE)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{country}{character; countries to filter by. Countries can either be
expressed as English names or
\href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}{ISO 2-letter country codes}.
English names are matched via regular expressions using
\link[=countrycode]{countrycode}, so there is some flexibility in names.}

\item{replace}{logical; multiple calls to \code{auk_country()} are additive,
unless \code{replace = FALSE}, in which case the previous list of countries to
filter by will be removed and replaced by that in the current call.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on a set of
countries. This function only defines the filter and, once all filters have
been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and perform the
filtering.
}
\details{
This function can also work with on an \code{auk_sampling} object if the
user only wishes to filter the sampling event data.
}
\examples{
# country names and ISO2 codes can be mixed
# not case sensitive
country <- c("CA", "United States", "mexico")
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_country(country)
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_country(ebd, country)
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-distance.R
\name{auk_distance}
\alias{auk_distance}
\title{Filter eBird data by distance travelled}
\usage{
auk_distance(x, distance, distance_units)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{distance}{integer; 2 element vector specifying the range of distances
to filter by. The default is to accept distances in kilometers, use
\code{distance_units = "miles"} for miles.}

\item{distance_units}{character; whether distances are provided in kilometers
(the default) or miles.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on the distance
travelled on the checklist. This function only defines the filter and, once
all filters have been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and
perform the filtering. Note that stationary checklists (i.e. point counts)
have no distance associated with them, however, since these checklists can
be assumed to have 0 distance they will be kept if 0 is in the range defined
by \code{distance}.
}
\details{
This function can also work with on an \code{auk_sampling} object if the
user only wishes to filter the sampling event data.
}
\examples{
# only keep checklists that are less than 10 km long
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_distance(distance = c(0, 10))
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_distance(ebd, distance = c(0, 10))
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-species.R
\name{auk_species}
\alias{auk_species}
\title{Filter the eBird data by species}
\usage{
auk_species(x, species, taxonomy_version, replace = FALSE)
}
\arguments{
\item{x}{\code{auk_ebd} object; reference to object created by \code{\link[=auk_ebd]{auk_ebd()}}.}

\item{species}{character; species to filter by, provided as scientific or
English common names, or a mixture of both. These names must match the
official eBird Taxomony (\link{ebird_taxonomy}).}

\item{taxonomy_version}{integer; the version (i.e. year) of the taxonomy. In
most cases, this should be left empty to use the version of the taxonomy
included in the package. See \code{\link[=get_ebird_taxonomy]{get_ebird_taxonomy()}}.}

\item{replace}{logical; multiple calls to \code{auk_species()} are additive,
unless \code{replace = FALSE}, in which case the previous list of species to
filter by will be removed and replaced by that in the current call.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on species. This
function only defines the filter and, once all filters have been defined,
\code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and perform the filtering.
}
\details{
The list of species is checked against the eBird taxonomy for
validity. This taxonomy is updated once a year in August. The \code{auk} package
includes a copy of the eBird taxonomy, current at the time of release;
however, if the EBD and \code{auk} versions are not aligned, you may need to
explicitly specify which version of the taxonomy to use, in which case
the eBird API will be queried to get the correct version of the taxonomy.
}
\examples{
# common and scientific names can be mixed
species <- c("Canada Jay", "Pluvialis squatarola")
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_species(species)
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_species(ebd, species)
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-ebd.R
\name{auk_ebd}
\alias{auk_ebd}
\title{Reference to eBird data file}
\usage{
auk_ebd(file, file_sampling, sep = "\\t")
}
\arguments{
\item{file}{character; input file. If file is not found as specified, it will
be looked for in the directory specified by the \code{EBD_PATH} environment
variable.}

\item{file_sampling}{character; optional input sampling event data (i.e.
checklists) file, required if you intend to zero-fill the data to produce a
presence-absence data set. This file consists of just effort information
for every eBird checklist. Any species not appearing in the EBD for a given
checklist is implicitly considered to have a count of 0. This file should
be downloaded at the same time as the basic dataset to ensure they are in
sync. If file is not found as specified, it will be looked for in the
directory specified by the \code{EBD_PATH} environment variable.}

\item{sep}{character; the input field separator, the eBird data are tab
separated so this should generally not be modified. Must only be a single
character and space delimited is not allowed since spaces appear in many of
the fields.}
}
\value{
An \code{auk_ebd} object storing the file reference and the desired
filters once created with other package functions.
}
\description{
Create a reference to an eBird Basic Dataset (EBD) file in preparation for
filtering using AWK.
}
\details{
eBird data can be downloaded as a tab-separated text file from the
\href{http://ebird.org/ebird/data/download}{eBird website} after submitting a
request for access. As of February 2017, this file is nearly 150 GB making it
challenging to work with. If you're only interested in a single species or a
small region it is possible to submit a custom download request. This
approach is suggested to speed up processing time.

There are two potential pathways for preparing eBird data. Users wishing to
produce presence only data, should download the
\href{http://ebird.org/ebird/data/download/}{eBird Basic Dataset} and reference
this file when calling \code{auk_ebd()}. Users wishing to produce zero-filled,
presence absence data should additionally download the sampling event data
file associated with the basic dataset This file contains only checklist
information and can be used to infer absences. The sampling event data file
should be provided to \code{auk_ebd()} via the \code{file_sampling} argument. For
further details consult the vignettes.
}
\examples{
# get the path to the example data included in the package
# in practice, provide path to ebd, e.g. f <- "data/ebd_relFeb-2018.txt
f <- system.file("extdata/ebd-sample.txt", package = "auk")
auk_ebd(f)
# to produce zero-filled data, provide a checklist file
f_ebd <- system.file("extdata/zerofill-ex_ebd.txt", package = "auk")
f_cl <- system.file("extdata/zerofill-ex_sampling.txt", package = "auk")
auk_ebd(f_ebd, file_sampling = f_cl)
}
\seealso{
Other objects: 
\code{\link{auk_sampling}()}
}
\concept{objects}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipe.R
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
  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-time.R
\name{auk_time}
\alias{auk_time}
\title{Filter the eBird data by checklist start time}
\usage{
auk_time(x, start_time)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{start_time}{character; 2 element character vector giving the range of
times in 24 hour format, e.g. \code{"06:30"} or \code{"16:22"}.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on a range of start
times for the checklist. This function only defines the filter and, once all
filters have been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and
perform the filtering.
}
\details{
This function can also work with on an \code{auk_sampling} object if the
user only wishes to filter the sampling event data.
}
\examples{
# only keep checklists started between 6 and 8 in the morning
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_time(start_time = c("06:00", "08:00"))
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_time(ebd, start_time = c("06:00", "08:00"))
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-sampling.R
\name{auk_sampling}
\alias{auk_sampling}
\title{Reference to eBird sampling event file}
\usage{
auk_sampling(file, sep = "\\t")
}
\arguments{
\item{file}{character; input sampling event data file, which contains
checklist data from eBird.}

\item{sep}{character; the input field separator, the eBird data are tab
separated so this should generally not be modified. Must only be a single
character and space delimited is not allowed since spaces appear in many of
the fields.}
}
\value{
An \code{auk_sampling} object storing the file reference and the desired
filters once created with other package functions.
}
\description{
Create a reference to an eBird sampling event file in preparation for
filtering using AWK. For working with the sightings data use \code{auk_ebd()},
only use \code{auk_sampling()} if you intend to only work with checklist-level
data.
}
\details{
eBird data can be downloaded as a tab-separated text file from the
\href{http://ebird.org/ebird/data/download}{eBird website} after submitting a
request for access. In the eBird Basic Dataset (EBD) each row corresponds
to a observation of a single bird species on a single checklist, while the
sampling event data file contains a single row for every checklist. This
function creates an R object to reference only the sampling data.
}
\examples{
# get the path to the example data included in the package
# in practice, provide path to the sampling event data
# e.g. f <- "data/ebd_sampling_relFeb-2018.txt"
f <- system.file("extdata/zerofill-ex_sampling.txt", package = "auk")
auk_sampling(f)
}
\seealso{
Other objects: 
\code{\link{auk_ebd}()}
}
\concept{objects}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-duration.R
\name{auk_duration}
\alias{auk_duration}
\title{Filter the eBird data by duration}
\usage{
auk_duration(x, duration)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{duration}{integer; 2 element vector specifying the range of durations
in minutes to filter by.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on the duration of
the checklist. This function only defines the filter and, once all filters
have been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and perform the
filtering. Note that checklists with no effort, such as incidental
observations, will be excluded if this filter is used since they have no
associated duration information.
}
\details{
This function can also work with on an \code{auk_sampling} object if the
user only wishes to filter the sampling event data.
}
\examples{
# only keep checklists that are less than an hour long
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_duration(duration = c(0, 60))
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_duration(ebd, duration = c(0, 60))
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-filter.R
\name{auk_filter}
\alias{auk_filter}
\alias{auk_filter.auk_ebd}
\alias{auk_filter.auk_sampling}
\title{Filter the eBird file using AWK}
\usage{
auk_filter(x, file, ...)

\method{auk_filter}{auk_ebd}(
  x,
  file,
  file_sampling,
  keep,
  drop,
  awk_file,
  sep = "\\t",
  filter_sampling = TRUE,
  execute = TRUE,
  overwrite = FALSE,
  ...
)

\method{auk_filter}{auk_sampling}(
  x,
  file,
  keep,
  drop,
  awk_file,
  sep = "\\t",
  execute = TRUE,
  overwrite = FALSE,
  ...
)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{file}{character; output file.}

\item{...}{arguments passed on to methods.}

\item{file_sampling}{character; optional output file for sampling data.}

\item{keep}{character; a character vector specifying the names of the columns
to keep in the output file. Columns should be as they appear in the header
of the EBD; however, names are not case sensitive and spaces may be
replaced by underscores, e.g. \code{"COMMON NAME"}, \code{"common name"}, and
\code{"common_NAME"} are all valid.}

\item{drop}{character; a character vector of columns to drop in the same
format as \code{keep}. Ignored if \code{keep} is supplied.}

\item{awk_file}{character; output file to optionally save the awk script to.}

\item{sep}{character; the input field separator, the eBird file is tab
separated by default. Must only be a single character and space delimited
is not allowed since spaces appear in many of the fields.}

\item{filter_sampling}{logical; whether the sampling event data should also
be filtered.}

\item{execute}{logical; whether to execute the awk script, or output it to a
file for manual execution. If this flag is \code{FALSE}, \code{awk_file} must be
provided.}

\item{overwrite}{logical; overwrite output file if it already exists}
}
\value{
An \code{auk_ebd} object with the output files set. If \code{execute = FALSE},
then the path to the AWK script is returned instead.
}
\description{
Convert the filters defined in an \code{auk_ebd} object into an AWK script and run
this script to produce a filtered eBird Reference Dataset (ERD). The initial
creation of the \code{auk_ebd} object should be done with \code{\link[=auk_ebd]{auk_ebd()}} and filters
can be defined using the various other functions in this package, e.g.
\code{\link[=auk_species]{auk_species()}} or \code{\link[=auk_country]{auk_country()}}. \strong{Note that this function typically takes
at least a couple hours to run on the full dataset}
}
\details{
If a sampling file is provided in the \link[=auk_ebd]{auk_ebd} object, this
function will filter both the eBird Basic Dataset and the sampling data using
the same set of filters. This ensures that the files are in sync, i.e. that
they contain data on the same set of checklists.

The AWK script can be saved for future reference by providing an output
filename to \code{awk_file}. The default behavior of this function is to generate
and run the AWK script, however, by setting \code{execute = FALSE} the AWK script
will be generated but not run. In this case, \code{file} is ignored and \code{awk_file}
must be specified.

Calling this function requires that the command line utility AWK is
installed. Linux and Mac machines should have AWK by default, Windows users
will likely need to install \href{https://www.cygwin.com}{Cygwin}.
}
\section{Methods (by class)}{
\itemize{
\item \code{auk_ebd}: \code{auk_ebd} object

\item \code{auk_sampling}: \code{auk_sampling} object
}}

\examples{
# get the path to the example data included in the package
# in practice, provide path to ebd, e.g. f <- "data/ebd_relFeb-2018.txt"
f <- system.file("extdata/ebd-sample.txt", package = "auk")
# define filters
filters <- auk_ebd(f) \%>\%
  auk_species(species = c("Canada Jay", "Blue Jay")) \%>\%
  auk_country(country = c("US", "Canada")) \%>\%
  auk_bbox(bbox = c(-100, 37, -80, 52)) \%>\%
  auk_date(date = c("2012-01-01", "2012-12-31")) \%>\%
  auk_time(start_time = c("06:00", "09:00")) \%>\%
  auk_duration(duration = c(0, 60)) \%>\%
  auk_complete()
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
filters <- auk_species(ebd, species = c("Canada Jay", "Blue Jay"))
filters <- auk_country(filters, country = c("US", "Canada"))
filters <- auk_bbox(filters, bbox = c(-100, 37, -80, 52))
filters <- auk_date(filters, date = c("2012-01-01", "2012-12-31"))
filters <- auk_time(filters, start_time = c("06:00", "09:00"))
filters <- auk_duration(filters, duration = c(0, 60))
filters <- auk_complete(filters)

# apply filters
\dontrun{
# output to a temp file for example
# in practice, provide path to output file
# e.g. f_out <- "output/ebd_filtered.txt"
f_out <- tempfile()
filtered <- auk_filter(filters, file = f_out)
str(read_ebd(filtered))
}
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-unique.R
\name{auk_unique}
\alias{auk_unique}
\title{Remove duplicate group checklists}
\usage{
auk_unique(
  x,
  group_id = "group_identifier",
  checklist_id = "sampling_event_identifier",
  species_id = "scientific_name",
  observer_id = "observer_id",
  checklists_only = FALSE
)
}
\arguments{
\item{x}{data.frame; the EBD data frame, typically as imported by
\code{\link[=read_ebd]{read_ebd()}}.}

\item{group_id}{character; the name of the group ID column.}

\item{checklist_id}{character; the name of the checklist ID column, each
checklist within a group will get a unique value for this field. The record
with the lowest \code{checklist_id} will be picked as the unique record within
each group. In the output dataset, this field will be updated to have a
full list of the checklist IDs that went into this group checklist.}

\item{species_id}{character; the name of the column identifying species
uniquely. This is required to ensure that removing duplicates is done
independently for each species. Note that this will not treat sub-species
independently and, if that behavior is desired, the user will have to
generate a column uniquely identifying species and subspecies and pass that
column's name to this argument.}

\item{observer_id}{character; the name of the column identifying the owner
of this instance of the group checklist. In the output dataset, the full
list of observer IDs will be stored (comma separated) in the new
\code{observer_id} field. The order of these IDs will match the order of the
comma separated checklist IDs.}

\item{checklists_only}{logical; whether the dataset provided only contains
checklist information as with the sampling event data file. If this
argument is \code{TRUE}, then the \code{species_id} argument is ignored and removing
of duplicated records is done at the checklist level not the species level.}
}
\value{
A data frame with unique observations, and an additional field,
\code{checklist_id}, which is a combination of the sampling event and group IDs.
}
\description{
eBird checklists can be shared among a group of multiple observers, in which
case observations will be duplicated in the database. This functions removes
these duplicates from the eBird Basic Dataset (EBD) or the EBD sampling event
data (with \code{checklists_only = TRUE}), creating a set of unique bird
observations. This function is called automatically by \code{\link[=read_ebd]{read_ebd()}} and
\code{\link[=read_sampling]{read_sampling()}}.
}
\details{
This function chooses the checklist within in each that has the
lowest value for the field specified by \code{checklist_id}. A new column is
also created, \code{checklist_id}, whose value is the taken from the field
specified in the \code{checklist_id} parameter for non-group checklists and from
the field specified by the \code{group_id} parameter for grouped checklists.

All the checklist and observer IDs for the checklists that comprise a given
group checklist will be retained as a comma separated string ordered by
checklist ID.
}
\examples{
# read in an ebd file and don't automatically remove duplicates
f <- system.file("extdata/ebd-sample.txt", package = "auk")
ebd <- read_ebd(f, unique = FALSE)
# remove duplicates
ebd_unique <- auk_unique(ebd)
nrow(ebd)
nrow(ebd_unique)
}
\seealso{
Other pre: 
\code{\link{auk_rollup}()}
}
\concept{pre}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-get-ebd-path.R
\name{auk_get_ebd_path}
\alias{auk_get_ebd_path}
\title{Return EBD data path}
\usage{
auk_get_ebd_path()
}
\value{
The path stored in the \code{EBD_PATH} environment variable.
}
\description{
Returns the environment variable \code{EBD_PATH}, which users are encouraged to
set to the directory that stores the eBird Basic Dataset (EBD) text files.
}
\examples{
auk_get_ebd_path()
}
\seealso{
Other paths: 
\code{\link{auk_get_awk_path}()},
\code{\link{auk_set_awk_path}()},
\code{\link{auk_set_ebd_path}()}
}
\concept{paths}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebird-species.R
\name{ebird_species}
\alias{ebird_species}
\title{Lookup species in eBird taxonomy}
\usage{
ebird_species(
  x,
  type = c("scientific", "common", "code", "all"),
  taxonomy_version
)
}
\arguments{
\item{x}{character; species to look up, provided as scientific or
English common names, or a mixture of both. Case insensitive.}

\item{type}{character; whether to return scientific names (\code{scientific}),
English common names (\code{common}), or 6-letter eBird species codes (\code{code}).
Alternatively, use \code{all} to return a data frame with the all the taxonomy
information.}

\item{taxonomy_version}{integer; the version (i.e. year) of the taxonomy.
Leave empty to use the version of the taxonomy included in the package.
See \code{\link[=get_ebird_taxonomy]{get_ebird_taxonomy()}}.}
}
\value{
Character vector of species identified by scientific name, common
name, or species code. If \code{type = "all"} a data frame of the taxonomy of
the requested species is returned.
}
\description{
Given a list of common or scientific names, check that they appear in the
official eBird taxonomy and convert them all to scientific names, common
names, or species codes. Un-matched species are returned as \code{NA}.
}
\examples{
# mix common and scientific names, case-insensitive
species <- c("Blackburnian Warbler", "Poecile atricapillus",
             "american dipper", "Caribou")
# note that species not in the ebird taxonomy return NA
ebird_species(species)

# use taxonomy_version to query older taxonomy versions
\dontrun{
ebird_species("Cordillera Azul Antbird")
ebird_species("Cordillera Azul Antbird", taxonomy_version = 2017)
}
}
\seealso{
Other helpers: 
\code{\link{auk_ebd_version}()},
\code{\link{auk_version}()},
\code{\link{get_ebird_taxonomy}()}
}
\concept{helpers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ebird_states}
\alias{ebird_states}
\title{eBird States}
\format{
A data frame with four variables and 3,145 rows:
\itemize{
\item \code{country}: short form of English country name.
\item \code{country_code}: 2-letter ISO country code.
\item \code{state}: state name.
\item \code{state_code}: 4 to 6 character state code.
}
}
\usage{
ebird_states
}
\description{
A data frame of state codes used by eBird. These codes are 4 to 6 characters,
consisting of two parts, the 2-letter ISO country code and a 1-3 character
state code, separated by a dash. For example, \code{"US-NY"} corresponds to New
York State in the United States. These state codes are required to filter by
state using \code{\link[=auk_state]{auk_state()}}.
}
\details{
Note that some countries are not broken into states in eBird and therefore do
not appear in this data frame.
}
\seealso{
Other data: 
\code{\link{bcr_codes}},
\code{\link{ebird_taxonomy}}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-year.R
\name{auk_year}
\alias{auk_year}
\title{Filter the eBird data to a set of years}
\usage{
auk_year(x, year, replace = FALSE)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{year}{integer; years to filter to.}

\item{replace}{logical; multiple calls to \code{auk_year()} are additive,
unless \code{replace = FALSE}, in which case the previous list of years to
filter by will be removed and replaced by that in the current call.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on a set of
years. This function only defines the filter and, once all filters have
been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and perform the
filtering.
}
\details{
For filtering to a range of dates use \code{auk_date()}; however,
sometimes the goal is to extract data for a given year or set of years, in
which case \code{auk_year()} is simpler. In addition, \code{auk_year()} can be used
to get data from discontiguous sets of years (e.g. 2010 and 2012, but not
2011), which is not possible with \code{auk_date()}. Finally, \code{auk_year()} can
be used in conjunction with \code{auk_date()} to extract data from a given range
of dates within a set of years (see example below).

This function can also work with on an \code{auk_sampling} object if the user
only wishes to filter the sampling event data.
}
\examples{
# years to filter to
years <- c(2010, 2012)
# set up filter
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_year(year = years)
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_year(ebd, years)

# filter to may and june of 2010 and 2012
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_year(year = c(2010, 2012)) \%>\% 
  auk_date(date = c("*-05-01", "*-06-30"))
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_county}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-county.R
\name{auk_county}
\alias{auk_county}
\title{Filter the eBird data by county}
\usage{
auk_county(x, county, replace = FALSE)
}
\arguments{
\item{x}{\code{auk_ebd} or \code{auk_sampling} object; reference to file created by
\code{\link[=auk_ebd]{auk_ebd()}} or \code{\link[=auk_sampling]{auk_sampling()}}.}

\item{county}{character; counties to filter by. eBird uses county codes
consisting of three parts, the 2-letter ISO country code, a 1-3 character
state code, and a county code, all separated by a dash. For example,
\code{"US-NY-109"} corresponds to Tompkins, NY, US. The easiest way to find a
county code is to find the corresponding \href{https://ebird.org/explore}{explore region} page and look at the URL.}

\item{replace}{logical; multiple calls to \code{auk_county()} are additive,
unless \code{replace = FALSE}, in which case the previous list of states to
filter by will be removed and replaced by that in the current call.}
}
\value{
An \code{auk_ebd} object.
}
\description{
Define a filter for the eBird Basic Dataset (EBD) based on a set of
counties This function only defines the filter and, once all filters have
been defined, \code{\link[=auk_filter]{auk_filter()}} should be used to call AWK and perform the
filtering.
}
\details{
It is not possible to filter by both county as well as country or
state, so calling \code{auk_county()} will reset these filters to all countries
and states, and vice versa.

This function can also work with on an \code{auk_sampling} object if the user only
wishes to filter the sampling event data.
}
\examples{
# choose tompkins county, ny, united states
system.file("extdata/ebd-sample.txt", package = "auk") \%>\%
  auk_ebd() \%>\%
  auk_county("US-NY-109")
  
# alternatively, without pipes
ebd <- auk_ebd(system.file("extdata/ebd-sample.txt", package = "auk"))
auk_county(ebd, "US-NY-109")
}
\seealso{
Other filter: 
\code{\link{auk_bbox}()},
\code{\link{auk_bcr}()},
\code{\link{auk_breeding}()},
\code{\link{auk_complete}()},
\code{\link{auk_country}()},
\code{\link{auk_date}()},
\code{\link{auk_distance}()},
\code{\link{auk_duration}()},
\code{\link{auk_extent}()},
\code{\link{auk_filter}()},
\code{\link{auk_last_edited}()},
\code{\link{auk_observer}()},
\code{\link{auk_project}()},
\code{\link{auk_protocol}()},
\code{\link{auk_species}()},
\code{\link{auk_state}()},
\code{\link{auk_time}()},
\code{\link{auk_year}()}
}
\concept{filter}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_ebd}
\alias{read_ebd}
\alias{read_ebd.character}
\alias{read_ebd.auk_ebd}
\alias{read_sampling}
\alias{read_sampling.character}
\alias{read_sampling.auk_ebd}
\alias{read_sampling.auk_sampling}
\title{Read an EBD file}
\usage{
read_ebd(x, sep = "\\t", unique = TRUE, rollup = TRUE)

\method{read_ebd}{character}(x, sep = "\\t", unique = TRUE, rollup = TRUE)

\method{read_ebd}{auk_ebd}(x, sep = "\\t", unique = TRUE, rollup = TRUE)

read_sampling(x, sep = "\\t", unique = TRUE)

\method{read_sampling}{character}(x, sep = "\\t", unique = TRUE)

\method{read_sampling}{auk_ebd}(x, sep = "\\t", unique = TRUE)

\method{read_sampling}{auk_sampling}(x, sep = "\\t", unique = TRUE)
}
\arguments{
\item{x}{filename or \code{auk_ebd} object with associated output
files as created by \code{\link[=auk_filter]{auk_filter()}}.}

\item{sep}{character; single character used to separate fields within a row.}

\item{unique}{logical; should duplicate grouped checklists be removed. If
\code{unique = TRUE}, \code{\link[=auk_unique]{auk_unique()}} is called on the EBD before returning.}

\item{rollup}{logical; should taxonomic rollup to species level be applied.
If \code{rollup = TRUE}, \code{\link[=auk_rollup]{auk_rollup()}} is called on the EBD before returning.
Note that this process can be time consuming for large files, try turning
rollup off if reading is taking too long.}
}
\value{
A data frame of EBD observations. An additional column,
\code{checklist_id}, is added to output files if \code{unique = TRUE}, that uniquely
identifies the checklist from which the observation came. This field is
equal to \code{sampling_event_identifier} for non-group checklists, and
\code{group_identifier} for group checklists.
}
\description{
Read an eBird Basic Dataset file using \code{\link[readr:read_delim]{readr::read_delim()}}. \code{read_ebd()}
reads the EBD itself, while read_sampling()` reads a sampling event data
file.
}
\details{
This functions performs the following processing steps:
\itemize{
\item Data types for columns are manually set based on column names used in the
February 2017 EBD. If variables are added or names are changed in later
releases, any new variables will have data types inferred by the import
function used.
\item Variables names are converted to \code{snake_case}.
\item Duplicate observations resulting from group checklists are removed using
\code{\link[=auk_unique]{auk_unique()}}, unless \code{unique = FALSE}.
}
}
\section{Methods (by class)}{
\itemize{
\item \code{character}: Filename of EBD.

\item \code{auk_ebd}: \code{auk_ebd} object output from \code{\link[=auk_filter]{auk_filter()}}

\item \code{character}: Filename of sampling event data file

\item \code{auk_ebd}: \code{auk_ebd} object output from \code{\link[=auk_filter]{auk_filter()}}. Must have
had a sampling event data file set in the original call to \code{\link[=auk_ebd]{auk_ebd()}}.

\item \code{auk_sampling}: \code{auk_sampling} object output from \code{\link[=auk_filter]{auk_filter()}}.
}}

\examples{
f <- system.file("extdata/ebd-sample.txt", package = "auk")
read_ebd(f)
# read a sampling event data file
x <- system.file("extdata/zerofill-ex_sampling.txt", package = "auk") \%>\%
  read_sampling()
}
\seealso{
Other import: 
\code{\link{auk_zerofill}()}
}
\concept{import}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-package.R
\docType{package}
\name{auk}
\alias{auk}
\title{\code{auk}: eBird Data Extraction and Processing in R}
\description{
Tools for extracting and processing eBird data from the eBird Basic Dataset
(EBD).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auk-set-ebd-path.R
\name{auk_set_ebd_path}
\alias{auk_set_ebd_path}
\title{Set the path to EBD text files}
\usage{
auk_set_ebd_path(path, overwrite = FALSE)
}
\arguments{
\item{path}{character; directory where the EBD text files are stored, e.g.
\code{"/home/matt/ebd"}.}

\item{overwrite}{logical; should the existing \code{EBD_PATH} be overwritten if it
has already been set in .Renviron.}
}
\value{
Edits .Renviron, then returns the AWK path invisibly.
}
\description{
Users of \code{auk} are encouraged to set the path to the directory containing the
eBird Basic Dataset (EBD) text files in the \code{EBD_PATH} environment variable.
All functions referencing the EBD or sampling event data files will check in
this directory to find the files, thus avoiding the need to specify the full
path every time. This will increase the portability of your code. Use this
function to set \code{EBD_PATH} in your .Renviron file; it is
also possible to manually edit the file.
}
\examples{
\dontrun{
auk_set_ebd_path("/home/matt/ebd")
}
}
\seealso{
Other paths: 
\code{\link{auk_get_awk_path}()},
\code{\link{auk_get_ebd_path}()},
\code{\link{auk_set_awk_path}()}
}
\concept{paths}
