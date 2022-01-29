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
